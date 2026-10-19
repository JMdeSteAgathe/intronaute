# Intronaute (French_approach_updated) — Pipeline Guide

**Tool:** [intronaute](https://github.com/JMdeSteAgathe/intronaute) v1.1.0 (vendored) — minor (U12) intron retention from RNA-seq BAM cohorts  
**Document version:** 2026-10-15  

---

## 1. Purpose and scientific model

This pipeline detects **genome-wide retention of minor (U12-type) introns** from **cohort RNA-seq BAMs**.

There is **no formal case/control design** in the statistics. Each sample is scored against the **cohort itself**:

- For every intron, each sample’s adjusted coverage (or junction) measure is compared to the **robust median and MAD** of all samples (or all other samples in small cohorts).
- Affected patients with minor-spliceosome defects are expected to be **rare** and to show **sparse but coordinated** elevation across many different U12 introns—not the same introns in every patient.

The **primary sample-level test** is the **Fisher tail test** (`tail_p` / `tail_q`): within each sample, are hits on **minor** introns enriched relative to hits on **matched major (U2) control** introns from the same genes? That contrast is **self-calibrating** for library quality and generic intronic noise.

---

## 2. Control introns and reference exons

Control introns and exons are two different parts of the same measurement. Every target, minor or control, has both an intron body and the exons beside it. "Control" names a second class of intron. The exons are the local expression reference for whichever intron they flank.

**Exons** are the reference sequence next to an intron: by default the last 300 bp of the upstream exon and the first 300 bp of the downstream exon, trimmed 5 bp back from the splice site so the count is not the boundary itself. Their depth is how much that gene is expressed there.

**Control introns** are one or two major (U2) introns from the same gene as a minor intron. They are chosen to have a similar length, and they are skipped if they overlap a known minor intron. Each one is its own target (`kind = control`, id like `CT00040_1_ARPC5`), with its own intron core, its own flanking exons, and its own splice junctions. They estimate ordinary intronic background in that gene: degradation, DNA contamination, nuclear pre-mRNA. A minor-spliceosome defect should raise minor introns above that background.

The coverage score is the ratio of those two depths, for every target:

\[
y = \log_2 \frac{\text{intron-core depth} + 0.5}{\text{flanking-exon depth} + 0.5}
\]

High exon depth with little intron depth means the gene is expressed and the intron was spliced. High intron depth relative to those exons means retention. If the flanking exons are shallower than 8, \(y\) is left missing.

### Prepare

Prepare builds both from the U12 intron list and the GTF.

Exons from the GTF are used twice. Overlapping exon sequence is cut out of the intron core, including antisense exons by default, so an overlapping exon is not counted as retention. Separately, the flanking windows above are stored as `exon_blocks`, `left_exon_blocks`, and `right_exon_blocks`.

Control introns are selected here and run through that same geometry. A target is usable only if the intron core is at least 40 bp and the exon reference is at least 50 bp. Each control stays linked to its minor intron by `pair_id`.

### Discover

Discover does not use either. It only lists BAMs.

### Quantify

Quantify counts both the same way, for minor and control targets.

Reads on `exon_blocks` become `exon_depth_mean` and the left and right exon depths. Reads on the intron core become `intron_depth_mean`. Reads at the two splice boundaries, which sit between exon and intron, become the junction counts: exact exon–exon splices (`n_EE`), unspliced reads crossing the boundary (`n_EI_left`, `n_EI_right`), and alternative donors and acceptors.

A few well-covered reference exons are also what the pipeline uses to guess whether the library is stranded.

### Analyse

Analyse is where their jobs split.

Exon depth is the denominator of \(y\), and the gate that drops a target when the flanks are too shallow. Sample QC uses it directly: median exonic depth, and the fraction of targets deep enough to score.

Control-intron \(y\) values are the within-sample background. Each sample's median control \(y\) is subtracted from every target, minor and control, so a sample with globally messy introns does not look retained. The junction ratio gets its own subtraction of the same kind: the median \(y_J\) of that sample's measurable control introns, computed from the junction ratios, is subtracted from every target. The correlation used to combine coverage z and junction z is estimated on the control introns.

The hit rule is then applied to both classes. Minor hits are `n_hits`. The same calls on control introns are `n_hits_control`. The Fisher test, enrichment, and MSRI all ask whether the minor introns are high relative to that sample's own control introns. The heatmap colors each minor intron by its z minus that sample's median control-intron z. The report then draws a PCA of the minor-intron z.


## 3. Repository layout

```
French_approach_updated/
  intronaute/           Python package (prepare, quantify, analyse, report)
  config/
    config.yaml         Paths plus optional flags (--config)
    all_bams_updated.tsv  Cohort BAM list + clinical/metadata columns
  config.example.json   Full parameter template (YAML or JSON via --config)
  submit_intronaute.sh  Submits intronaute.sh
  intronaute.sh         Slurm settings, then calls the intronaute tool
  requirements.txt      pip install notes (slurm wheels)
  pyproject.toml        Package metadata (intronaute 1.1.0)
  docs/
    PIPELINE_GUIDE.md   This document
  output/with_metadata/
    annotation/         Target table from prepare
    results/            One analysis. This is the outdir.
      counts/           Cached per-sample counts
      sample_mapped_PCA_minor_introns/     Sample map from minor-intron z
      samples_ranked_minor_introns/        Sample scores, hit table, and rank plots
      most_retained_minor_introns/         Heatmaps of retained minor introns
      shared_minor_intron_genes/           UpSet plots and shared-gene lists
      quality_control/                     Depth, measurability, and intronic background
      report.html
      samplesheet.tsv
```

`submit_intronaute.sh` only submits the job. `intronaute.sh` holds the Slurm settings and calls the `intronaute` tool. Slurm logs for a run are written to `<outdir>/Logs/`.

Upstream generic documentation also lives in `README.md` (install, simulation notes, hit modes).

---

## 4. End-to-end workflow

Typical production:

```bash
cd /path/to/your/intronaute/pipeline
./submit_intronaute.sh --config config/config.yaml
./submit_intronaute.sh --config config/config.yaml analyse
./submit_intronaute.sh --local --config config/config.yaml all
```

**Stages:**

| Stage | Command | Heavy? | Cached? |
|-------|---------|--------|---------|
| Prepare | `intronaute prepare` | ~minutes | `targets.tsv.gz` |
| Discover | `intronaute discover` | light | Run only when `config.yaml` has `bam_dir` or `bam_dirs` and no `samplesheet` |
| Quantify | `intronaute quantify` | **hours** (BAM I/O) | Per-sample `counts/<sample>.tsv.gz` |
| Analyse | `intronaute analyse` | ~minutes | Recomputable anytime |

**Configured inputs** (from `config/config.yaml`):

- **Introns:** `French_approach/metadata/GRCh38_ensembl110_u12-type_introns.ods`
- **GTF:** `French_approach/metadata/gencode.v47.chr_patch_hapl_scaff.basic.annotation.gtf.gz`
- **Samplesheet or BAM directories:** `config/all_bams_updated.tsv` (~170 STAR BAMs, hg38), or `bam_dir` / `bam_dirs` when no samplesheet is set
- **Results outdir:** e.g. `output/with_metadata/results/`
- **Targets:** `output/with_metadata/annotation/targets.tsv.gz` (derived from parent of `outdir`)

Default **hit mode** in the submit script: **`combined`** (coverage + junction, with veto).

---

## 5. Call graph: what calls what

### Submit scripts

- **`submit_intronaute.sh`** requires **`--config`** and submits **`intronaute.sh`**. That job script holds the Slurm CPUs, memory, time, and **`WORKERS`**, then calls the **`intronaute`** tool. It sources **`load_deps.sh`** from the script directory when the project virtualenv is not already active. Slurm logs go in **`<outdir>/Logs`**.

### Python entry

- **`python -m intronaute`** → **`intronaute/cli.py:main()`** → subcommand `func` (`cmd_prepare`, `cmd_discover`, `cmd_quantify`, `cmd_analyse`, `cmd_run`, `cmd_selftest`).

### Module dependencies (simplified)

```
cli.py
  prepare  → regions.build_targets
  discover → bamio (index check), pandas
  quantify → quantify.load_targets, quantify.quantify_bam
              → bamio, fastbam (or pysam)
  analyse  → stats.load_counts, stats.* scoring, stats.pca
              → report.fig_*, report.write_html
```

Supporting modules: **`intervals.py`** (interval algebra), **`fastbam.py`** (pure-Python BAM reader), **`selftest.py`**, **`_minibam.py`**.

---

## 6. Steps

### 6.1 Step 1: prepare — build targets

**Purpose:** Turn a **U12 intron list** + **GTF** into a table of genomic **measurement targets**.

**For each listed minor intron:**

1. **Intron core** — intron minus regions exonic in other transcripts/genes (avoid measuring exons or antisense overlap as “retention”).
2. **Reference exons** — default **300 bp** at the 3′ end of the upstream exon and 5′ end of the downstream exon (flanking normalization).
3. **Control introns** — **1–2 major (U2)** introns from the **same gene**, length-matched, for per-sample background subtraction in analyse.

**Output:** `targets.tsv.gz` with columns including `target_id`, `kind` (`minor` / `control`), coordinates, block strings, `usable`, gene metadata.

**Coordinate conventions:** BED is 0-based half-open; TSV/ODS default 1-based inclusive unless `--coord-base 0`.

---

### 6.2 Step 2: discover — BAM inventory (optional)

**Purpose:** Build a samplesheet from BAM files when `config.yaml` has no `samplesheet` and does have `bam_dir` or `bam_dirs`.

Discover walks those directories, finds `.bam` files, records whether a `.bai` index exists, and writes `samplesheet.tsv`. Columns are `sample`, `bam`, `run`, `has_index`, and `size_GB`. If a samplesheet is already set, discover is not run, and quantify and analyse use that file. A prepare-only run does not need either key.

**`bam_dir` is one directory. `bam_dirs` is one directory or several.** Use `bam_dirs` when the BAMs are split across folders. Two forms are accepted.

A block list, one path per line:

```yaml
bam_dirs:
  - "/path/to/first/bam/folder"
  - "/path/to/second/bam/folder"
```

A flow list on one line:

```yaml
bam_dirs: ["/path/to/first/bam/folder", "/path/to/second/bam/folder"]
```

A single path is also valid, either as `bam_dir: "/path/to/bam/folder"` or as `bam_dirs` with one entry. If both keys are set, the directories are combined and a repeated path is kept once. A path that does not start with `/` is resolved from the directory that contains `submit_intronaute.sh`.

Discover walks each directory, including folders nested inside it. The **run** for a new sample is the directory name (the last path component that is not a generic label such as `bam` or `star`). The same sample name in two directories is kept distinct by appending that run. While `samplesheet` is set, discover is skipped, so `bam_dirs` is not read until that line is removed.

**`run` is the label supplied on the sheet.** If the sheet has no `run` column, a column named `batch` or `batches` is read as `run`. Capitalization does not matter (`batch`, `Batch`, `BATCH`, `batches`, `Batches`, `BATCHES`). If discover is run again and a sample already has a non-empty `run`, that value is kept. The BAM directory name is written only when `run` is empty. Extra columns already on the sheet are kept and lined up by `sample`. QC and the run table group on `run`.

---

### 6.3 Step 3: quantify — per-BAM counting

**Purpose:** For each sample, read only **indexed regions** needed for all targets (merged chunks); never scan whole BAMs.

**Per target, stored in `counts/<sample>.tsv.gz`.** One row per minor intron and per control intron. `intron_start` and `intron_end` are the prepare coordinates: 0-based and half-open, so `intron_len = intron_end - intron_start`.

| Column | What it is |
|--------|------------|
| `target_id` | Target id from prepare. Minor introns look like `MI00767_SMAP2`. Controls look like `CT00767_1_SMAP2`. |
| `kind` | `minor` (U12 intron) or `control` (matched U2 intron from the same gene). |
| `pair_id` | Id shared by a minor intron and the control introns chosen for it, for example `MI00767`. |
| `gene_id` | Host gene id from the GTF. |
| `gene_name` | Host gene name. |
| `chrom` | Chromosome. |
| `intron_start` | Start of the full annotated intron, 0-based. |
| `intron_end` | End of the full annotated intron, half-open: the first base after the intron. |
| `strand` | Host-gene strand, `+` or `-`. |
| `intron_len` | Full intron length in bp (`intron_end - intron_start`). This is the annotated intron, not the trimmed core that is counted. |
| `core_len` | Bases of intron core included in the depth arrays. The core is the intron after end-trimming and after subtracting overlapping exons. |
| `exon_len` | Bases of flanking reference exons included in the depth arrays. |
| `intron_depth_mean` | Mean read depth across the intron core. |
| `intron_depth_med` | Median read depth across the intron core. |
| `intron_frac_cov` | Fraction of core bases with depth greater than 0. |
| `exon_depth_mean` | Mean read depth across both flanking exon windows, pooled. This is the denominator of the coverage ratio in analyse. |
| `exon_depth_med` | Median read depth across those same exon bases. |
| `exon_depth_left` | Mean read depth on the genomic-left exon window only. |
| `exon_depth_right` | Mean read depth on the genomic-right exon window only. |
| `n_intron_reads` | Reads whose alignment overlaps the intron core. |
| `n_EE` | Reads spliced on this exact intron: the gap runs from `intron_start` to `intron_end`. These are the exon–exon junctions that removed this intron. |
| `n_EI_left` | Reads that cover the genomic-left boundary (`intron_start`) and extend at least `min_overhang` bp (default 6) into both the exon and the intron. Unspliced evidence at that boundary. |
| `n_EI_right` | The same count at the genomic-right boundary (`intron_end`). |
| `n_alt_donor` | Spliced reads that use the genomic-left boundary as one end of a junction and some other coordinate as the other end. This is alternative splicing at that site, not retention. The name means genomic left, which is the donor only on the plus strand. |
| `n_alt_acceptor` | Spliced reads that use the genomic-right boundary as one end of a junction and some other coordinate as the other end. Genomic right is the acceptor only on the plus strand. |
| `truncated` | `1` if counting stopped early because a fetched region exceeded `max_reads_per_region` (default 300,000 reads). Depths from that region are incomplete. `0` otherwise. |

**Left and right are genomic, not donor and acceptor.** `n_EI_left` and `n_alt_donor` are the boundary at `intron_start`. `n_EI_right` and `n_alt_acceptor` are the boundary at `intron_end`. On a plus-strand gene the left end is the donor and the right end is the acceptor. On a minus-strand gene those ends swap: genomic left is the acceptor and genomic right is the donor. Analyse swaps the count-file columns into `n_EI_donor`, `n_EI_acceptor`, `n_alt_5p`, and `n_alt_3p` only when it writes `hits.tsv.gz`. Those four names are not columns in the count file.

**Sidecar** `counts/<sample>.json`, one file per sample:

| Field | What it is |
|-------|------------|
| `bam` | Absolute path of the BAM that was counted. |
| `bam_size` | BAM size in bytes at quantification time. The cache misses if this changes. |
| `bam_mtime` | BAM modification time, seconds since the epoch. Recorded only. |
| `n_targets` | Rows written to the count file. |
| `n_chunks` | Merged genomic regions read from the BAM. |
| `n_truncated_chunks` | Regions that hit the read cap. |
| `library` | `forward`, `reverse`, or `unstranded`. With `library=auto`, this is detected from reference exons. |
| `min_mapq` | Minimum mapping quality kept (default 10). |
| `min_overhang` | Bases required on each side of a boundary before an unspliced crossing counts (default 6). |
| `drop_dup` | Whether duplicate reads were removed. |
| `backend` | BAM reader used (`fastbam`, `pysam`, or `bamnostic`). |
| `seconds` | Wall time for this sample. |
| `median_exon_depth` | Median of `exon_depth_mean` across the targets in this file. |

**Parallelism:** `--workers`. `intronaute.sh` sets `WORKERS` and passes that number.

**Cache:** Skip sample if output exists and BAM + quant params unchanged (`--force` to recompute).

**Filters:** MAPQ (default 10), duplicate removal, stranded read filtering when library is stranded (`library=auto` detects per sample).

---

### 6.4 Step 4: analyse — scores, PCA, report

**Pipeline inside `cmd_analyse`:**

1. Load counts → matrices (samples × targets).
2. Load target metadata. After the scores are computed, copy samplesheet columns onto the score table (`run`, plus any extra columns such as `diagnosis`). Those columns are labels for the table and the PCA colour menu. They are not used to compute the scores.
3. Build **coverage** `y`, `y_adj`, `z_cov`, `d_cov` (defined below).
4. If junction columns present: **junction** `y_J`, `theta`, `z_j`, `d_j`, `z_comb`, correlation **ρ** on controls.
5. Define **measurable** introns (default ≥50% samples with finite `y`; defined below).
6. **Call hits** per `--hit-mode`.
7. **Sample scores** (Fisher tail, enrichment, MSRI, Stouffer, and the other per-sample summaries defined below).
8. **PCA** on measurable minor introns. It does not change the hit calls. The map uses the minor-intron **z** (`z_comb` in combined mode). Each sample is one row. Each measurable minor intron is one column. Control introns are not columns. The sample-wide control subtraction inside `y_adj` is already inside this z.
9. Write **`report.html`** and the result directories inside **`outdir`**. This is where the plots and tables are created:

   - **`sample_mapped_PCA_minor_introns/`** — sample map from the minor-intron z-scores, one picture per colour, plus the PC tables.
   - **`samples_ranked_minor_introns/`** — one row of scores per sample, the minor-intron hit table, and the rank plots.
   - **`most_retained_minor_introns/`** — heatmaps of minor-intron retention.
   - **`shared_minor_intron_genes/`** — UpSet plots and the lists of genes shared across samples.
   - **`quality_control/`** — depth, how much of each sample could be measured, and intronic background.

   The count cache in **`counts/`** is written earlier, by quantify.

#### 6.4.1 Coverage quantities

`y`, `y_adj`, `z_cov`, and `d_cov` are four successive versions of the same coverage measurement, one number per sample and intron. The paragraphs under each equation in §6.4.5 are the full account: what the controls are, which median is across samples, why the pseudocount is 0.5, why the junction floor is 10 reads, and what the combined z is doing.

**`y`** is the raw coverage ratio. For sample *i* and intron *j* it is

\[
y_{i,j} = \log_2 \frac{\text{intron depth}_{i,j} + 0.5}{\text{flanking exon depth}_{i,j} + 0.5}
\]

Intron depth is the mean depth on the intron core from Step 3. Exon depth is the mean depth on the flanking exon windows. A retained intron has reads inside the intron relative to its own exons, so `y` is higher. If the flanking exon depth is below 8, `y` is missing.

**`y_adj`** is that ratio after removing the sample's own intronic background. The pipeline takes the median `y` of that sample's measurable U2 control introns and subtracts it from every intron in the sample:

\[
y_{\text{adj},i,j} = y_{i,j} - \operatorname{median}_{k \in \text{controls}} y_{i,k}
\]

A sample with degradation, genomic DNA, or a high nuclear fraction has high `y` on ordinary introns too. Subtracting that median leaves the excess over the sample's own background.

**`z_cov`** is the robust z-score of `y_adj`, one intron at a time, across samples. It asks how unusual this sample is for this intron compared with everyone else. That median is across the cohort for this intron. It is not a median inside one sample, and it is not split by run or diagnosis. See the robust z formula in §6.4.5.

\[
z_{\text{cov},i,j} = 0.6745 \times \frac{y_{\text{adj},i,j} - \operatorname{median}_j}{\operatorname{MAD}_j}
\]

The 0.6745 factor puts the score on a standard-deviation scale when the values are roughly normal. The MAD is floored so a tiny spread cannot create a huge z. For this cohort of about 170 samples the median and MAD include the sample itself. `z_cov` is this coverage z. The junction track has its own z, `z_j`.

**`d_cov`** is the same comparison without dividing by the spread. It is how far `y_adj` sits above that intron's cohort median, in log2 units:

\[
d_{\text{cov},i,j} = y_{\text{adj},i,j} - \operatorname{median}_j
\]

A coverage hit needs both a high z and a real shift: `z_cov` at least 3 and `d_cov` at least 0.5, which is about 1.41 times the cohort median.

#### 6.4.2 Measurable introns

An intron is **measurable** when enough samples actually have a coverage ratio for it. The default is half the cohort.

`y` is a real number only when that sample's flanking exon depth is at least 8. Below that, `y` is missing, so that sample cannot be used for that intron. Step 4 then counts, for each intron, the fraction of samples with a real `y`. The intron is measurable when that fraction is at least **`min_frac_measurable`**.

**Default ≥50%** means `min_frac_measurable` is **0.5** unless you pass `--min-frac-measurable`. In a cohort of about 170 samples, at least about 85 samples must have exon depth ≥ 8 on that intron. Introns below that cutoff are dropped: their `y` is cleared for every sample, and they are left out of the z-scores, the hit calls, and the PCA. The cohort median and MAD are not computed from a handful of samples.

#### 6.4.3 Junction quantities

`theta`, `y_J`, `z_j`, `d_j`, `z_comb`, and **ρ** are the junction-side scores. They are computed per intron from the splice-boundary reads, rather than from intron-body depth. In the Step 4 list above, the junction z-score is `z_j` and the correlation is **ρ**.

**`theta`** is the fraction of reads at the two splice sites that are unspliced. EI is the unspliced reads that cross either boundary (`n_EI_left` + `n_EI_right`). S is the spliced reads at those same sites: the exact exon–exon junction counted at both ends (`2 × n_EE`), plus alternative partners (`n_alt_donor` + `n_alt_acceptor`).

\[
\theta = \frac{\text{EI}}{\text{EI} + S}
\]

A spliced intron has theta near 0. A retained intron has theta near 1. If EI + S is below 10 reads, theta is left missing.

**`y_J`** is that same comparison on the log2 scale used for the coverage ratio:

\[
y_J = \log_2 \frac{\text{EI} + 0.5}{S + 0.5}
\]

High `y_J` means unspliced boundary reads outweigh spliced ones. Any intron with EI + S below 10 is set to missing. What remains is scored in the same two steps as the coverage ratio. The full formulas are in the junction measure in §6.4.5.

**Control subtraction.** For each sample, the median `y_J` of that sample's measurable U2 control introns is subtracted from every intron's `y_J`. That median is computed from the junction ratios, so a sample with messy splice junctions on ordinary introns does not look retained. If fewer than 20 control introns are measurable, nothing is subtracted. The result is the adjusted junction ratio, `y_{J,adj}`.

**Robust z and effect size.** For one intron, line up `y_{J,adj}` across the cohort. The middle of that list is the cohort median, and the MAD is the spread around it.

\[
z_j = 0.6745 \times \frac{y_{J,\text{adj}} - \text{Median}}{\text{MAD}}, \qquad
d_j = y_{J,\text{adj}} - \text{Median}
\]

**`z_j`** is how unusual this sample's splice-site retention is for this intron. **`d_j`** is the same gap left in log2 units, so it stays a size. A value of 0 matches the cohort middle. A value of 0.5 means the adjusted unspliced-to-spliced ratio is about 1.41 times that middle. A junction hit needs both: `z_j` at least 3, and `d_j` at least **`min_delta_junction`** (default 0.5).

**`z_comb`** turns the two z-scores into one. `z_cov` says the intron body has more coverage than usual. `z_j` says the splice sites look unspliced. Coverage and junctions both rise when the same intron is retained, so adding them raw would double-count that shared evidence. The full indexed formula is in the combined z in §6.4.5.

\[
z_{\text{comb}} = \frac{z_{\text{cov}} + z_j}{\sqrt{2 + 2\rho}}
\]

**ρ** is how strongly the two scores move together. It is measured on the U2 control introns, where that shared movement is background rather than a U12 signal, and then clipped to the range 0–0.95. Each score is clipped at ±6 before the correlation, so a few extreme introns do not set ρ. If fewer than 200 control introns have both scores, ρ is treated as 0.

The square-root term is larger when ρ is larger. When the two scores normally agree, their sum is shrunk so agreeing does not look like extra evidence. When only one of the two z-scores is measurable, `z_comb` is that one score. In the combined hit mode used for the heatmap, `z_comb` is the z that is plotted and tested. It is one number, on about the same scale as a single z, for both the intron body and the splice sites looking unusual. A large positive value means this sample is high on both, after removing the agreement those two measures show on ordinary introns.

#### 6.4.4 Sample scores

The quantities above are one number per sample and per intron. The sample scores collapse that matrix to one row per sample, written to `sample_scores.tsv`. In combined mode, the z that enters MSRI, Stouffer, and the medians is `z_comb`. A hit is whatever the chosen hit mode calls. The same call is counted on minor introns and on U2 control introns.

Read a row by asking three questions, in this order:

1. Are minor introns called more often than this sample's own control introns? That is `n_hits` against `n_hits_control`, then enrichment and the Fisher tail.
2. Is the shift broad, or only a few introns? That is `median_z_minor`.
3. Do the minor introns sit above the control introns, not merely above zero? That is MSRI. Stouffer does not make this comparison.

Three rows from the without-metadata run show the difference. ID350075 is the strongest minor-intron outlier. ID029422 is a real but weaker one. G467 has a lot of minor hits and even more control hits.

| Sample | Minor hits | Control hits | Enrichment | `tail_q` | Median z, minor | MSRI | Stouffer |
|--------|------------|--------------|------------|----------|-----------------|------|----------|
| ID350075 | 623 | 50 | 26 | 0 | 5.31 | 5.13 | 143 |
| ID029422 | 68 | 13 | 10 | \(8 \times 10^{-21}\) | 1.23 | 1.36 | 31 |
| G467 | 107 | 149 | 1.5 | 0.020 | 0.46 | 0.08 | 15 |

ID350075 and ID029422 pass all three questions. G467 fails the first and the third: its minor hits are not in excess of its control hits, and its minor z-scores are not higher than its control z-scores. A raw hit count would have ranked G467 above ID029422. The other columns do not.

**`n_hits`** is how many minor introns are called in that sample. **`n_hits_control`** is the same call on the control introns. A minor-spliceosome defect should raise the minor count above the control count. A noisy library, degradation, or DNA contamination raises both, because those affect U2 introns too. G467 is that second pattern: 107 minor hits, 149 control hits.

**`n_minor_measurable`** and **`n_control_measurable`** are how many introns of each class have a finite z, the denominator of the hit rate. ID350075 has 718 measurable minor introns, so 623 hits is 623/718, not 623 out of an unknown total. The Fisher test and the enrichment are computed only when both denominators are at least 20. Below that, the rate is too unstable and `tail_p` and enrichment are left missing.

**Enrichment** is the minor hit rate divided by the control hit rate in the same sample:

\[
\text{enrichment} = \frac{(n_{\text{hits}} + 0.5)\,/\,(n_{\text{minor}} + 1)}{(n_{\text{hits, control}} + 0.5)\,/\,(n_{\text{control}} + 1)}
\]

The added 0.5 and 1 keep a zero count from making the ratio zero or infinite. Enrichment near 1 means minor introns are called about as often as control introns. ID350075 is about 26. G467 is about 1.5, close to no enrichment. The five high-hit samples sit between about 14 and 33.

**Fisher tail** (`tail_p`) is the one-sided Fisher exact test of that same contrast. For one sample the table is:

|  | called | not called |
|--|--------|------------|
| minor introns | `n_hits` | measurable minors minus `n_hits` |
| control introns | `n_hits_control` | measurable controls minus `n_hits_control` |

For ID350075 the cells are 623 and 95 on the minor row, and 50 and 1,460 on the control row. The test asks whether a split that lopsided would turn up if minor and control introns shared one hit rate. The alternative is "greater," so only an excess of minor hits counts. A small `tail_p` means the minor calls are too numerous to be that sample's background rate.

**`tail_q`** is the Benjamini–Hochberg false-discovery rate of `tail_p` across every sample in the cohort. With about 170 samples, some will have a small `tail_p` by chance. `tail_q` is the column to quote. ID350075 and ID029422 have `tail_q` far below any usual cutoff. G467's `tail_q` is 0.02: a hint, not the same evidence. This is the main sample-level test. Prefer it, and the comparison of `n_hits` with `n_hits_control`, over the raw hit count.

**`median_z_minor`** and **`median_z_control`** are the median z of each class. They describe the typical intron, not the strongest one. A median near 0 is the cohort background: half the introns are above the cohort and half are below. ID350075 has a minor median of 5.3, so the typical minor intron is already a strong outlier, and a control median near 0, so its U2 introns are ordinary. Retention there is broad. G467's two medians are both near 0.5, so its typical intron, minor or control, looks like the rest of the cohort. Its 107 hits are a tail, and the control introns have a tail of the same kind.

**MSRI** (minor spliceosome retention index) answers whether the minor introns, as a group, sit above the control introns. It is the mean z of the minor introns minus the mean z of the control introns, after clipping each z to the range −4 to 8. Clipping keeps one extreme intron from dominating the mean. A trim would instead throw away the high scores, which is the signal when only some introns are retained. A positive MSRI means this sample's minor introns sit higher than its own control introns. ID350075 is 5.1. G467 is 0.08, effectively no minor-versus-control shift, which matches its enrichment near 1.

**Stouffer** adds the minor-intron z-scores and divides by the square root of how many were measurable:

\[
\text{Stouffer} = \frac{\sum z_{\text{minor}}}{\sqrt{n_{\text{minor measurable}}}}
\]

Many introns shifted up by a little, or fewer introns shifted up by a lot, both raise it. Under a pure null the values scatter around 0. It does not subtract the control introns, so a sample whose every intron is a bit high, minor and control alike, still gets a positive Stouffer. G467's Stouffer is 15, which looks elevated on its own, while its MSRI is about 0 and its enrichment is about 1.5. That is why Stouffer is the tie-break on the ranking tables, and not the test. The minor-versus-control contrast is MSRI and the Fisher tail.

**`n_hits_coverage_only_rule`** and **`n_hits_junction_only_rule`** are written when junction counts exist. Each is the minor-intron hit count you would get from that rule alone, at the same thresholds. The word "only" means "this rule by itself," not "introns called by one rule and not the other." ID350075 has 555 coverage-rule hits and 544 junction-rule hits, and 623 combined hits. Combined can exceed either single rule because an intron with only one of the two measures available can still be called. ID029422 has 69 coverage-rule hits and 24 junction-rule hits, close to its 68 combined hits, so its calls are mostly the coverage rule with the junction veto still satisfied.

---

#### 6.4.5 Formulas

Each block below is one equation and then what that equation does. The symbols, the direction of each median, the floors, and why the constants are there are in the paragraphs under the formula.

##### Coverage ratio

For sample *i*, intron *j*:

\[
y_{i,j} = \log_2 \frac{\text{intron depth}_{i,j} + a}{\text{flanking exon depth}_{i,j} + a}
\]

with pseudo-count **a = 0.5**. If flanking exon depth &lt; **`min_exon_depth`** (default 8), \(y_{i,j} = \text{NaN}\) (not measurable).

The same **0.5** is added in the junction ratio, \(\log_2((\text{EI}+0.5)/(S+0.5))\). It keeps the log2 ratio finite when a depth or a junction count is zero. Without it, an intron depth of 0 is \(\log_2(0)\), which is undefined, and a spliced count \(S = 0\) divides by zero. Adding the same small number on the top and the bottom turns those cases into a large but finite score: no intronic reads becomes a very negative \(y\), and reads that are all unspliced becomes a very positive \(y_J\).

The value 0.5 is half a read, the usual continuity correction for a log ratio of counts. Next to the floors already required for a real measurement (exon depth at least 8, and at least 10 reads at the splice sites), half a read barely moves a well-supported ratio. Samples below those floors are set to missing. The pseudocount does not turn a failed measurement into a real one.

##### Control-intron normalization (within sample)

\[
y_{\text{adj},i,j} = y_{i,j} - \operatorname{median}_{k \in \text{controls}} y_{i,k}
\]

This subtracts one number per sample: the typical coverage ratio of that sample's ordinary introns. A U12 intron then has to sit above that sample's own background.

**What the controls are.** They are major (U2) introns, the ordinary spliceosome introns, not other U12 introns. Prepare picks them from the GTF before any BAM is read. For each listed U12 intron it looks in the same gene and keeps up to two other introns that:

- are at least 60 bp long
- do not overlap a known U12 intron
- are closest in length to that U12 intron

Each control is measured the same way as a U12 intron: depth in the intron core divided by depth in the flanking exons, which produces its own `y`.

**How they are used.** The formula does not subtract only the one or two controls matched to that gene. For each sample it takes every measurable control intron in the cohort, computes the median of their `y` values, and subtracts that single median from every intron in that sample, U12 and control alike. If a sample has fewer than 20 measurable controls, nothing is subtracted.

**What that does.** Degradation, leftover genomic DNA, or a high fraction of nuclear RNA raises intron coverage across the board. The control median is that sample-wide rise. After subtraction, `y_adj` is how far an intron sits above the sample's ordinary introns. A minor-spliceosome defect should lift U12 introns above that line. The ordinary introns should stay near it.

##### Robust z (one intron, across the cohort)

The median here is **per intron, across samples**. It is not the median inside one sample, and it is not split by run, diagnosis, or any other group.

For intron *j*, line up `y_adj` from every sample that has a real value for that intron. **Median**_j is the middle of that list. The MAD is the middle of the distances from that middle. Sample *i* is then scored by how far its own `y_adj` sits from that cohort middle:

\[
z_{i,j} = 0.6745 \times \frac{y_{\text{adj},i,j} - \text{Median}_j}{\text{MAD}_j}
\]

The 0.6745 factor puts the score on a standard-deviation scale when the values are roughly normal.

The control-intron normalization was the other median: one number from that sample's control introns, subtracted from every intron in the sample. The robust z looks the other way: one intron, compared across the cohort.

For this cohort of about 170 samples, the sample being scored is included in **Median**_j. Leave-one-out, which drops that sample from the median and the MAD, applies only when the cohort has **`loo_max`** or fewer samples (default 60).

**MAD floor:** If an intron's values barely vary, the MAD is raised to at least 0.05 (or to a low quantile of the MADs across introns, whichever is larger). A tiny spread cannot turn a small difference into a huge z.

##### Effect size (coverage)

**`d_cov`** is how much higher this sample's adjusted coverage is than the typical sample for that same intron. It uses the same **Median**_j from the robust z: the middle of `y_adj` for intron *j* across the cohort.

\[
d_{\text{cov},i,j} = y_{\text{adj},i,j} - \text{Median}_j
\]

`z_cov` takes that same gap and divides it by how spread out the cohort is. **`d_cov`** leaves the gap in log2 units, so it stays a size, not a rarity. A value of 0 means this sample matches the cohort middle. A value of 0.5 means its adjusted intron-to-exon ratio is about 1.41 times that middle.

A coverage hit needs both. `z_cov` at least 3 means the sample is unusual. **`d_cov`** at least **`min_delta`** (default 0.5) means the shift is large enough to count. A tiny difference can still look rare if the intron barely varies across samples. This cutoff blocks those.

##### Junction measure (FRASER-like)

This is the splice-site version of the coverage chain above. Coverage asks whether the intron body is deep relative to the flanking exons. The junction measure asks whether the reads at the two splice sites are unspliced. The sentence "then the same control subtraction and robust z → `z_j`, `d_j`" means: once `y_J` exists, run it through those same two steps. The outputs are **`z_j`** and **`d_j`**.

**What goes into the ratio.** Both splice sites are pooled into two counts.

- **EI** is the unspliced reads that cross either boundary: `n_EI_left` + `n_EI_right`.
- **S** is the spliced reads at those same sites: the exact exon–exon junction counted at both ends (`2 × n_EE`), plus alternative partners (`n_alt_donor` + `n_alt_acceptor`).

\[
\theta_{i,j} = \frac{\text{EI}}{\text{EI} + S}, \quad
y_{J,i,j} = \log_2 \frac{\text{EI} + 0.5}{S + 0.5}
\]

The 0.5 is the same pseudocount as in the coverage ratio. **`theta`** is the fraction of those reads that are unspliced. A spliced intron has theta near 0. A retained intron has theta near 1. **`y_J`** is that same comparison on the log2 scale. High `y_J` means unspliced boundary reads outweigh spliced ones.

**The 10-read floor.** `y_J` and `theta` are ratios of two counts. With only a few reads, one read moving from spliced to unspliced swings the ratio a lot. A ratio built from 3 reads is treated as noise. The default requires **EI + S** at least **`min_junction_reads`** (10) before the junction score is kept. Those 10 reads are the total at both sites together. Below 10, `y_J` and `theta` are set to missing, so that sample does not get a junction z for that intron. The cutoff is `--min-junction-reads`.

What remains is scored like this.

**Control subtraction**, the junction version of the control-intron normalization above. For each sample, take the median `y_J` of that sample's measurable U2 control introns and subtract it from every intron's `y_J`, U12 and control alike. That median is computed from the junction ratios. It is a new number for the junction track. A sample with messy splice junctions on ordinary introns does not look retained, because that sample-wide rise is removed. If fewer than 20 control introns are measurable on the junction track, nothing is subtracted. The result is the adjusted junction ratio:

\[
y_{J,\text{adj},i,j} = y_{J,i,j} - \operatorname{median}_{k \in \text{controls}} y_{J,i,k}
\]

**Robust z and effect size**, the junction versions of the robust z and the effect size above. For one intron, line up `y_{J,adj}` from every sample that has a real value. The middle of that list is the cohort median. The MAD is the middle of the distances from that middle. This median is per intron, across samples. It is the other direction from the control subtraction, which was one number inside a single sample.

\[
z_{J,i,j} = 0.6745 \times \frac{y_{J,\text{adj},i,j} - \text{Median}_j}{\text{MAD}_j}, \qquad
d_{J,i,j} = y_{J,\text{adj},i,j} - \text{Median}_j
\]

The 0.6745 factor puts `z_j` on a standard-deviation scale when the values are roughly normal. For this cohort of about 170 samples, the sample being scored is included in the median and the MAD. Leave-one-out, which drops that sample, applies only when the cohort has **`loo_max`** or fewer samples (default 60). If an intron's junction ratios barely vary, the MAD is raised to at least 0.05 (or to a low quantile of the MADs across introns, whichever is larger). A tiny spread cannot turn a small difference into a huge z.

**`z_j`** is how unusual this sample is for this intron: the gap divided by the cohort spread. **`d_j`** is that same gap left in log2 units, so it stays a size. A value of 0 means this sample matches the cohort middle for that intron. A value of 0.5 means the adjusted unspliced-to-spliced ratio is about 1.41 times that middle (\(2^{0.5} \approx 1.41\)).

A junction hit needs both. `z_j` at least 3 means the sample is unusual. **`d_j`** at least **`min_delta_junction`** (default 0.5) means the shift is large enough to count. A tiny difference can still look rare if the intron barely varies across samples. This cutoff blocks those.

##### Combined z

The control introns here are the same ordinary U2 introns from the control-intron normalization above. They are the background set: introns that a minor-spliceosome defect is not expected to retain.

Each intron already has two scores. `z_cov` says the intron body has more coverage than usual. `z_j` says the splice sites look unspliced. `z_comb` turns those two into one score:

\[
z_{\text{comb}} = \frac{z_{\text{cov}} + z_J}{\sqrt{2 + 2\rho}}
\]

Adding them raw would double-count. Both scores often move together, because they are two views of the same transcripts. **ρ** is how strongly they move together. It is measured on the control introns, where that shared movement is background rather than a U12 signal. Each score is clipped at ±6 before the correlation, so a few extreme introns do not set ρ. ρ itself is then kept between 0 and 0.95, so the two scores are never treated as the same piece of evidence. If fewer than 200 control introns have both scores, ρ is treated as 0.

The square-root term is larger when ρ is larger. When the two scores normally agree, their sum is shrunk so agreeing does not look like extra evidence. When only one score exists for that sample and intron, `z_comb` is just that one score.

**`z_comb`** is one number, on about the same scale as a single z, for "both the intron body and the splice sites look unusual." A large positive value means this sample is high on both, after removing the agreement those two measures show on ordinary introns.

---

## 7. Hit modes and per-intron calls

| `hit_mode` | Intron is a hit when… | Explanation |
|---|---|---|
| `coverage` | \(z_{\text{cov}} \geq z_{\text{thr}}\) and \(d_{\text{cov}} \geq \text{min\_delta}\) | The intron body has to look both rare and large. `z_cov` is how unusual this sample's adjusted intron-to-exon coverage is for this intron, compared with the cohort (§6.4.5). `z_thr` is that rarity cutoff, default 3. `d_cov` is the size of the same shift in log2 units (§6.4.5). `min_delta` is that size cutoff, default 0.5, about 1.41 times the cohort-middle ratio (\(2^{0.5} \approx 1.41\)). Both must clear their cutoffs. If either score is missing, the intron is not a hit. This mode uses the intron-body coverage only. |
| `junction` | \(z_J \geq z_{\text{thr}}\) and \(d_J \geq \text{min\_delta\_junction}\) | The splice sites have to look both rare and large. `z_J` is how unusual this sample's adjusted unspliced-to-spliced ratio is at the two splice sites (§6.4.5). `z_thr` is the same rarity cutoff, default 3. `d_J` is the size of that junction shift in log2 units. `min_delta_junction` is that size cutoff, default 0.5, about 1.41 times the cohort-middle junction ratio (\(2^{0.5} \approx 1.41\)). Both must clear their cutoffs. If the junction score is missing, the intron is not a hit. This mode uses the splice-site reads only. |
| `combined` | \(z_{\text{comb}} \geq z_{\text{thr}}\), at least one effect guard passes, and the intron is **not vetoed** | Three checks all have to pass. `z_comb` is the single score from §6.4.5 that merges intron-body coverage and splice-site evidence, on about the same scale as one z. It must be a real number and at least `z_thr` (default 3). An effect guard is a size check: `d_cov` at least `min_delta` (default 0.5), or `d_J` at least `min_delta_junction` (default 0.5). At least one of those sizes must pass. **Not vetoed** means the junction veto below did not fire. If `z_J` exists and is below `junction_veto` (default 1.0), the intron is not a hit, whatever the coverage. |
| `both` | Coverage rule **and** junction rule both pass | The coverage row and the junction row must both be true for the same sample and intron. Each uses `z_thr` (default 3) and its own size cutoff (`min_delta` for coverage, `min_delta_junction` for junctions). This is the strictest mode. It is meant to be paired with a lower `--z-threshold`, because asking both measures to clear 3 on their own removes most events. |

**Veto (`junction_veto`, default 1.0).** This applies in **combined** mode. High coverage inside an intron is not called retention when the splice sites themselves look normal.

A combined-mode hit has to pass three checks at once:

1. `z_comb` is a real number and at least `z_thr` (default 3).
2. At least one effect-size guard passes. The coverage shift `d_cov` is at least `min_delta` (default 0.5), or the junction shift `d_J` is at least `min_delta_junction` (default 0.5).
3. The intron is not vetoed.

**Junctions are measurable** when this sample has a real `z_J` for this intron. That requires enough reads at the two splice sites to build the junction ratio in §6.4.5. Unspliced crossings plus spliced reads, EI + S, must be at least `min_junction_reads` (default 10). The intron also has to be measurable on the junction side, so a cohort median and MAD exist and `z_J` comes out finite. If those reads are too few, `z_J` is missing. A missing score makes `z_J < veto` false, so this veto does not apply. A coverage-driven `z_comb` can still be a hit when the splice sites were too sparse to score.

**`z_J` below the veto** means those splice-site reads do not look unusually unspliced. The default veto is **1.0**. That number is a z cutoff, the same kind of score as `z_thr`, not a read count. A junction z under 1 sits near the cohort middle, or below it. The reads at the splice sites are being spliced about as often as in the rest of the cohort. The hit threshold for calling an intron unusual is 3, so the veto sits much lower: it blocks the call when the splice sites fail to show even a mild unspliced shift. A `z_J` between 1 and 3 does not trigger the veto. Combined mode can still call that intron when `z_comb` is at least 3 and a size guard passes. That is the case where coverage is high and the splice sites are only mildly elevated.

**Treated as background** means the intron is not called a hit. It is left out of `n_hits`. Under the same rule it is left out of `n_hits_control`. This holds even when the intron body is full of reads and `z_comb` is high. Reads piled in the middle of an intron, without unspliced reads crossing the splice sites, can come from degradation or leftover genomic DNA. Those reads raise intron-body coverage without raising the unspliced fraction at the exon–intron boundaries. The veto treats that coverage as generic intronic noise. Real retention should show up at the splice sites as well: unspliced reads that cross from the exon into the intron.

The same hit definition, including this veto in combined mode, applies to minor introns and to control introns. Minor hits are counted in **`n_hits`**. Control hits are counted in **`n_hits_control`**.

### 7.1 Where these calls are saved

Analyse writes the hit calls into **`samples_ranked_minor_introns/`** inside the directory named by **`outdir`**.

| File | What it holds |
|------|----------------|
| `hits.tsv.gz` | One row per called **minor** intron, for one sample. Gene, coordinates, the calling score `z`, the adjusted ratio `y_adj`, depths, and the splice-site counts. |
| `sample_scores.tsv` | One row per sample. `n_hits` is how many minor introns passed the hit rule. `n_hits_control` is how many U2 control introns passed the same rule. |

**`hits.tsv.gz` is the per-intron list.** Each row is one sample and one minor intron that passed the rule in the table above. It carries `gene_name`, `gene_id`, `chrom`, `intron_start`, `intron_end`, `strand`, and `intron_len`. `z` is the score the chosen hit mode used (`z_cov`, `z_J`, or `z_comb`). `y_adj` is the adjusted coverage ratio, or the adjusted junction ratio when the mode is `junction`. The row also has the raw depths: `intron_depth_mean`, `exon_depth_mean`, `exon_depth_left`, `exon_depth_right`, `n_intron_reads`, and `intron_frac_cov`. When junctions were counted, the splice-site reads are relabeled by transcript direction: `n_EI_donor` and `n_EI_acceptor` are the unspliced crossings, `n_EE` is the exon–exon spliced reads, and `n_alt_5p` and `n_alt_3p` are the alternate splice sites. Those rows also include `z_coverage`, `z_junction`, and `theta`. In `combined` or `both` mode the row includes `z_combined` as well.

Rows are sorted by sample, then by `z` from high to low. Each sample keeps only its top **`--max-events-per-sample`** minor hits (default **200**). **`n_hits` is the full count, before that cap.** A sample with 250 minor hits has `n_hits` of 250 and 200 rows in `hits.tsv.gz`. A sample with no minor hits has `n_hits` of 0 and no rows in that file. If the whole cohort has no minor hits, `hits.tsv.gz` is not written.

Control-intron hits use the same rule, including the veto in combined mode. They are **not** given rows in `hits.tsv.gz`. Their total is only the column `n_hits_control`.

**`sample_scores.tsv` has every sample**, including samples with zero hits. `n_hits` counts minor introns that passed. `n_hits_control` counts control introns that passed. Enrichment, `tail_p`, `tail_q`, `MSRI`, and `stouffer` are built from those two counts and from the z scores. Those columns are described in §8. The pictures of these scores (`samples_ranked_minor_introns/score_distribution.png`, `most_retained_minor_introns/heatmap_n_hits.png`, `report.html`) are described in §9.1. The intron-level calls themselves are the two files above.

---

## 8. Per-sample scores (output columns)

| Column | Meaning |
|--------|---------|
| `sample` | Sample ID |
| `run` | Run label from the samplesheet. A column named `batch` or `batches` (any capitalization) is read as `run` when `run` is absent. This is the group on the QC figure and in the run table. Discover fills it from the BAM directory only when the cell is empty. |
| *(other samplesheet columns)* | Copied after scoring, for example `diagnosis`, `minor_spliceopathy`, `tissue`, `size_GB`. Not used in the median, MAD, or Fisher test. `bam` and `has_index` are left off. |
| `n_hits` | Full count of minor introns called hits, before the 200-row cap in `hits.tsv.gz` |
| `n_hits_control` | Same rule on U2 control introns. A count only; those introns are not rows in `hits.tsv.gz` |
| `n_minor_measurable` | Minor introns with finite z in analysed set |
| `n_control_measurable` | Control introns measurable |
| `enrichment` | (hits/minor meas.) / (control hits/control meas.) |
| `tail_p` | Fisher exact, one-sided: minor hits vs control hits **same sample** |
| `tail_q` | BH-FDR on `tail_p` |
| `MSRI` | Winsorised mean z(minor) − winsorised mean z(control) |
| `stouffer` | Σ z(minor) / √n_measurable |
| `median_z_minor`, `median_z_control` | Location summaries |
| `n_hits_coverage_only_rule`, `n_hits_junction_only_rule` | If junctions present: side-by-side rule counts |

**Ranking:** `sample_scores.tsv` is sorted by `n_hits`, then `stouffer`. The report switch in §9.1 lists the top samples by largest `n_hits`, largest `n_hits_control`, smallest `tail_q`, or largest `enrichment`. Missing `tail_q` and `enrichment` stay at the bottom. Every one of those tables breaks ties on Stouffer, high first. **Interpretation:** Prefer **`tail_q`** and **`n_hits` vs `n_hits_control`** over raw hit count alone.

---

## 9. Outputs and file reference

A run writes two directories under `output/with_metadata/`. **`annotation/`** is the target table from prepare. **`results/`** is one analysis, and it is the `outdir` in `config/config.yaml`. What each hit file contains is in §7.

The same pictures are in `results/report.html` and as PNG files in the subdirectories below. The HTML report is the one you can search, colour, zoom, and switch. A PNG in a subdirectory is that same picture with none of those controls.

| Directory or file | What it holds |
|---|---|
| `annotation/` | Target table built by prepare. |
| `results/counts/` | Cached per-sample counts. Keep this directory unless you are re-counting the BAMs. |
| `results/sample_mapped_PCA_minor_introns/` | Sample map from minor-intron z-scores: one picture per colour, plus the PC tables. |
| `results/samples_ranked_minor_introns/` | One row of scores per sample, the minor-intron hit table, and the rank plots. |
| `results/most_retained_minor_introns/` | Heatmaps of minor-intron retention. |
| `results/shared_minor_intron_genes/` | UpSet plots and the lists of genes shared across samples. |
| `results/quality_control/` | Depth, how much of each sample could be measured, and intronic background. |
| `results/report.html` | The report, with the same figures as the PNGs in the subdirectories. Search, colour, zoom, and the ranking switch work on this page. |
| `results/samplesheet.tsv` | The sample list used for that run. |

---

### 9.1 Report and figures — what they show

`report.html` opens with a contents bar. Each link jumps to that section, and the bar stays at the top of the page while you scroll. A section is listed only when that part of the report was written.

#### Section 1 — PCA map (`report.html`)

The report draws one map. A point is one sample. The colour menu is `n_hits`, −log10(`tail_q`), log2(`enrichment`), `MSRI`, `stouffer`, and each samplesheet column. Numeric columns use a continuous ramp. A category such as `diagnosis` uses a discrete palette. Typing in the box highlights a sample. Scrolling zooms, dragging pans, and clicking a point pins its label. Pin top N labels the samples with the most minor-intron hits. Samplesheet columns colour a point. They do not decide where it sits. Hit calls are computed before the map and are not changed by it.

The axes are PC1 and PC2 of the minor-intron z matrix. In combined mode that z is `z_comb`. Each sample is one row. Each measurable minor intron is one column. Control introns are not in this matrix. Where a sample has no z for an intron, that cell is filled with the cohort median for that intron, and extreme values are clipped so one intron cannot take over an axis. The sign is set so further right means a higher z across minor introns. PC2 is the next pattern in that same profile, often a sequencing batch or a different subset of introns. The files are in `sample_mapped_PCA_minor_introns/`.

The map already includes the sample-wide control subtraction inside `y_adj`, because that subtraction happens before z is computed.

#### `score_distribution.png`

This is report section **2. Highest-ranking samples**. Above the plots, a switch lists the top samples by largest `n_hits`, largest `n_hits_control`, smallest `tail_q`, or largest `enrichment`. The paragraph above the switch defines each of those scores. **Enrichment** is the minor-intron hit rate divided by the control-intron hit rate in the same sample. Each rate is (hits + 0.5) / (measurable introns of that class + 1). A value near 1 means the two classes are called about equally often. A larger value means the hits are concentrated on minor introns. Samples with fewer than 20 measurable introns of either class have no enrichment, and those rows stay at the bottom of that table. The same rule puts a missing `tail_q` at the bottom of the `tail_q` table. Ties on every ranking break on Stouffer, high first. The report draws the two panels below that switch on a canvas. Search highlights a sample on both plots. Click a dot to pin its name. The PNG file is `samples_ranked_minor_introns/score_distribution.png`.

- **Left:** Samples ranked by **`n_hits`** (log rank).
- **Right:** Log–log scatter **`n_hits_control`** vs **`n_hits`**; dashed line = no enrichment. Points **above** the line = more minor than control hits.
- **Controls:** colour by a score or a samplesheet column, search for a sample, and pin the top N samples by `n_hits`.

#### Heatmaps

Three files share one layout. The report shows one at a time and lets you switch.

- **Color:** robust **z** of that minor intron minus the median robust z of the same sample's U2 control introns, drawn from −6 to +6. Red means the minor intron is retained above the cohort and above that sample's own ordinary introns. A sample whose non-minor introns are retained moves toward white. Blue is below that contrast. The scale is a legend under the left side of each heatmap, and the same legend sits under section 3 of the report. Hit calls still use the unadjusted z.
- **Columns:** every measurable minor intron, labeled by gene. A gene with several minor introns has one column per intron. Columns are ordered by the largest adjusted z in that heatmap's rows, strongest on the left. The report keeps the page width and scrolls sideways. Switching the ranking can change that order.
- **`heatmap_n_hits.png`:** top 30 samples by `n_hits`.
- **`heatmap_n_hits_control.png`:** top 30 samples by `n_hits_control` (generic intronic hits).
- **`heatmap_tail_q.png`:** top 30 samples by smallest `tail_q` (strongest minor-versus-control enrichment).

Those three files are in `most_retained_minor_introns/`.

#### Shared genes (UpSet)

`shared_minor_intron_genes/upset_*.png`, with the same switch in `report.html`. `upset_shared_genes.tsv` lists the gene names.

Each person is the set of genes with minor-specific retention. A gene counts when at least one of its minor introns is a hit and every control intron that shares that minor intron's `pair_id` is not a hit. A hit on the matched U2 intron is widespread retention, and that pair does not count. Another minor intron in the same gene can still qualify the gene when its own paired controls are not hits. A gene counts once even if several of its minor introns qualify. The full hit matrix is used, not the 200-row cap in `hits.tsv.gz`.

A column of dots is one combination of people. The bar above it is how many genes have minor-specific retention in exactly those people. Blue means two or more people. Grey means one person. The bars on the left are each person's gene count. The table under the plot names the genes in the blue combinations. Every group uses the same bar and dot size. A group with many combinations is wider than the page, and the report scrolls that plot sideways.

The switch offers a group only when that group has at least one person:

- **minor spliceopathy:** `minor_spliceopathy` is any text other than a placeholder or a negative. Placeholders, compared without regard to case, are blank, `0`, `control`, `N/A`, `NA`, null, and `none`. `NULL`, `null`, and a missing cell are the same placeholder. Negatives are `false`, `no`, and `non-MS`. `MS` and `ms` are the same included value, and so are `Yes` and `yes`, `TRUE` and `true`, and `1`. Any other label is included the same way.
- **Diagnosis:** one group for each value in `diagnosis` other than those same placeholders. `RNU4ATAC`, `RNU12`, and `RNU6ATAC` are included, and so is any other diagnosis. Values that differ only by case are one group and one button, so `Yes` and `yes` are the same diagnosis. The button uses a single spelling.

A group with one person does not get an UpSet plot. That same part of the report lists the genes with a minor-intron hit whose paired control introns are not hits. The heading is the sample id plus the diagnosis when the sample has one, and **minor spliceopathy** otherwise. If that one person is both the only spliceopathy sample and the only person with a diagnosis, the list is shown once, under the diagnosis. The gene names are also written to `singleton_minor_hit_genes.tsv`.

The next section of `report.html`, **Build an UpSet plot**, is the same comparison with samples you choose. The menu labels each sample with diagnosis and spliceopathy status. Buttons for **minor spliceopathy** and for each diagnosis fill the selection when those labels are present. A placeholder diagnosis does not get a button. Two or more names draw the plot and list the shared genes. One name lists that person's genes with minor-specific retention.

#### Top events table

An **event** is one minor intron in one sample that passed the hit rule in §7. Control-intron hits are not rows here. The **`z`** column is the score that rule uses: **`z_cov`** in coverage mode, **`z_j`** in junction mode, and **`z_comb`** in combined mode. In **`both`** mode the call still requires the coverage rule and the junction rule, and **`z`** is **`z_comb`**, the merged score the rows are ranked by.

The table is the 60 events with the highest **`z`** across the cohort, taken from **`hits.tsv.gz`**. One sample can fill several rows. Columns are gene, coordinates, z, depths, **`n_EI_donor`**, **`n_EI_acceptor`**. High intronic signal without those boundary-crossing reads is more likely an overlapping transcript or a mapping artifact than retention.

#### `qc.png`

Three panels by **`run`**: median exonic depth; fraction measurable; control-intron background offset. The report draws the same three panels on a canvas and embeds the PNG. Colour and search work as on the score plots. Pin top N marks the samples with the most minor-intron hits on those axes. Flags low depth, poor measurability, or high generic intronic noise.

#### Run table (if multiple runs)

Median/max hits per run — batch structure.

#### Parameters

Section 8 lists every key `config/config.yaml` can set: the input paths, the prepare flags, the quantify flags, and the analyse flags. Analyse rows are the values used for this report. `min_mapq`, `min_overhang`, `keep_duplicates`, `library`, and `engine` are read from the count metadata when those files agree, because that is what the cached counts used. The other prepare and quantify rows come from the config file when the key is set, and from the default when it is left commented out.

---

## 10. Configuration (YAML, JSON, Slurm)

### 10.1 `config/config.yaml`

`submit_intronaute.sh` requires `--config /path/to/config.yaml`. It reads the paths from that file and passes them to `intronaute.sh`, which calls every `intronaute` step. Required keys are `outdir`, `introns`, `gtf`, and either `samplesheet` or `bam_dir` / `bam_dirs` (see [§6.2](#62-step-2-discover--bam-inventory-optional)). Add `analyse` to run only the stats and report. Slurm starts in the directory that contains the scripts, where `.venv` and `load_deps.sh` live. Job logs are written to `<outdir>/Logs`.

`bam_dirs` is how to give discover more than one BAM folder. A block list or a flow list both work, and a one-entry list is allowed:

```yaml
bam_dirs:
  - "/path/to/first/bam/folder"
  - "/path/to/second/bam/folder"
```

```yaml
bam_dirs: ["/path/to/first/bam/folder", "/path/to/second/bam/folder"]
```

`bam_dir` remains a single directory. Setting both keys combines them and drops a repeated path. Relative paths are resolved from the script directory. A `samplesheet` line skips discover, so these keys are unused until that line is removed.

Every other flag may be set in the same file. The names match the command-line flags with underscores (`min_mapq`, `z_threshold`, `hit_mode`, `loo_max`, and the rest). A key that is present replaces that flag's default. A key that is absent keeps the default. A flag on the `intronaute` command line wins over the file. Commented lines in `config/config.yaml` list the flags and their defaults.

Uncomment a threshold or `hit_mode` in the YAML to replace that default. A commented line keeps the default. For `hit_mode`, that default is `combined` when you use `submit_intronaute.sh`. `--hit-mode` on the submit command wins over the file. `WORKERS` in `intronaute.sh` is passed as `--workers`, so it wins over a `workers` key in the YAML.

Changing `min_mapq`, `min_overhang`, `keep_duplicates`, or `library` applies on the next quantify. Pass `--force` (or `--force-quantify` on the submit script) so cached counts are rebuilt.

### 10.2 `config.example.json`

The same keys, written out with values, for `intronaute <step> --config config.example.json`. YAML and JSON both work. Command-line flags win.

### 10.3 Slurm

`intronaute.sh` holds the job settings: `#SBATCH` CPUs, memory, time, partition, and account, plus `WORKERS`. `submit_intronaute.sh` submits that file. `--force-prepare`, `--force-quantify`, and `--hit-mode` are arguments of the submit command.

Install (slurm): see **`requirements.txt`** — prefer binary wheels, then `pip install -e . --no-deps`, then `intronaute selftest`.

---

## 11. Performance, caching, and re-runs

- **Quantify:** Dominates runtime (NFS latency); use multiple workers; local BAM copy optional.
- **Re-analyse only:** Change `--hit-mode` or thresholds without re-reading BAMs (needs junction columns in cache for non-coverage modes).
- **Re-quantify `--force`:** Required after changing GTF/intron list, quant params, or after **1.1.0 strandedness fix** on stranded libraries (1.0.0 caches).
- **Prepare `--force-prepare`:** New target table when annotation or intron list changes.

Reference scale (upstream README): ~218 samples × ~3121 targets → analyse ~41 s, ~0.9 GB RAM.

---

## 12. Caveats and interpretation

1. **Cohort size:** Unstable below ~20 samples; LOO below ~60. ~170 samples is adequate.
2. **No case/control in model:** Labels validate; they are not used in Fisher/MAD reference (avoid leakage).
3. **Compare `n_hits` to `n_hits_control`:** Similar values → generic noise, not U12-specific signal.
4. **Per-intron hits:** Weaker than aggregate **`tail_q`**; confirm in IGV / RT-PCR.
5. **Mixed strandedness** in one cohort is suboptimal; check reported library type per sample in count JSON.
6. **Index required:** No `.bai` → sample skipped or unusable.
7. **Metadata is a label, not a covariate:** samplesheet columns such as `minor_spliceopathy` are on `sample_scores.tsv` and in the PCA colour menu. They do not enter the Fisher test or the cohort median.

---

## Appendix A: Symbols and notation

| Symbol | Meaning |
|--------|---------|
| \(i\) | Sample index |
| \(j\) | Intron / target index |
| \(a\) | Pseudo-count (0.5) added to depths and junction counts |
| \(y_{i,j}\) | Log₂ intron/exon coverage ratio (raw) |
| \(y_{\text{adj},i,j}\) | \(y\) after subtracting sample’s control-intron median |
| \(y_{J,i,j}\) | Log₂ junction odds ratio (EI vs S) |
| \(\theta_{i,j}\) | EI / (EI + S), junction “retention fraction” at splice sites |
| \(\text{Median}_j\) | Cohort median of \(y_{\text{adj}}\) for intron \(j\) |
| \(\text{MAD}_j\) | Median absolute deviation from \(\text{Median}_j\) |
| \(z_{i,j}\) | Robust z-score (coverage track unless stated) |
| \(z_{\text{cov}}, z_J, z_{\text{comb}}\) | Coverage, junction, combined z |
| \(d_{\text{cov}}, d_J\) | Effect size vs cohort median (\(y_{\text{adj}} - \text{Median}\)) |
| \(\rho\) | Correlation between \(z_{\text{cov}}\) and \(z_J\) on control introns |
| \(z_{\text{thr}}\) | Hit z threshold (default 3) |
| EI | Reads crossing splice sites without splicing (donor + acceptor) |
| S | Spliced reads at donor + acceptor |
| PC1, PC2, … | Principal components of the minor-intron z matrix |
| NaN | Not measurable (low depth, low junction support, etc.) |

---

## Appendix B: Abbreviations and acronyms

| Abbr. | Expansion |
|-------|-----------|
| **BAM** | Binary Alignment Map (aligned reads) |
| **BAI** | BAM index (required for region queries) |
| **BH** | Benjamini–Hochberg (FDR on `tail_p`) |
| **BGZF** | Block gzip (BAM compression) |
| **CIGAR** | Read alignment operation string |
| **EI** | Exon–intron (unspliced boundary-crossing reads) |
| **EE** | Exon–exon (correctly spliced junction reads) |
| **FDR** | False discovery rate |
| **GTF** | Gene transfer format (annotation) |
| **LOO** | Leave-one-out (robust reference for small n) |
| **MAD** | Median absolute deviation |
| **MAPQ** | Mapping quality |
| **MS** | Minor spliceopathy (metadata label in samplesheet) |
| **MSRI** | Minor spliceosome retention index (winsorised z contrast) |
| **ODS** | OpenDocument spreadsheet (intron list input) |
| **PCA** | Principal component analysis |
| **RT-PCR** | Reverse transcription PCR (validation) |
| **Slurm** | Cluster job scheduler |
| **STAR** | Spliced read aligner (BAM source) |
| **TSV** | Tab-separated values |
| **U2** | Major spliceosome (canonical introns) |
| **U12** | Minor spliceosome (minor-type introns) |
| **FRASER** | Finding RAre Splicing Events in RNA-seq |
| **GENCODE** | Reference gene annotation |
| **IGV** | Integrative Genomics Viewer |
| **NYGC** | New York Genome Center (batch label) |
| **RNU12** | Small nuclear RNA gene (minor spliceosome; clinical label) |

**Column name shorthand in hits:**

| Column | Meaning |
|--------|---------|
| `n_EE` | Exact exon–exon junction reads |
| `n_EI_left`, `n_EI_right` | Genomic left/right boundary crossings (raw in counts) |
| `n_EI_donor`, `n_EI_acceptor` | Biologically oriented 5′/3′ splice site crossings |
| `n_alt_5p`, `n_alt_3p` | Alternative splicing at one boundary (not retention) |

**Target `kind`:**

| Value | Meaning |
|-------|---------|
| `minor` | U12 intron target |
| `control` | Matched U2 intron in same gene |

---

## Appendix C: CLI command quick reference

```bash
intronaute prepare   --introns LIST --gtf GTF [--targets PATH]
intronaute discover  --bam-dir DIR [DIR ...] [--samplesheet PATH]
intronaute quantify  --targets PATH --samplesheet PATH --outdir DIR [--workers N] [--force]
intronaute analyse   --targets PATH --samplesheet PATH --outdir DIR [--hit-mode MODE] [--config YAML_OR_JSON]
intronaute run       # prepare + discover + quantify + analyse
intronaute selftest  # BAM engine check
```

**Aliases:** `analyze` = `analyse`.

---

## Appendix D: Default thresholds

| Parameter | Default | Role |
|-----------|---------|------|
| `min_exon_depth` | 8 | Minimum flanking exon depth to measure coverage |
| `min_frac_measurable` | 0.5 | Fraction of cohort needed for intron to enter analyse set |
| `z_threshold` | 3.0 | Hit z cutoff |
| `min_delta` | 0.5 | Min log₂ effect above cohort median (coverage) |
| `min_junction_reads` | 10 | Min EI+S for junction measure |
| `min_delta_junction` | 0.5 | Junction effect-size guard |
| `junction_veto` | 1.0 | Combined mode: veto if \(z_J\) below this |
| `loo_max` | 60 | Exact LOO below this sample count |
| `min_cohort` | 20 | Warning threshold |
| `exon_window` | 300 bp | Flanking exon window (prepare) |
| `min_mapq` | 10 | Read filter (quantify) |

---

*End of document.*
