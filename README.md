# Intronaute

Detect **minor-intron (U12) retention** across a cohort of RNA-seq BAM files.

Built for locked-down machines: **no compiler, no admin rights, no `pysam`**.
The BAM reader is bundled and written in pure Python, so a plain
`pip install` is all you need — which matters because there is no `pysam` wheel
for Windows.

There is **no case/control design**. The cohort is its own background: every
sample is scored against the robust median and MAD of all the others,
intron by intron. This works because a minor-spliceosome defect is rare and
hits a set of introns that differs from patient to patient — so an affected
sample stands out as an outlier against its own cohort. Point it at a few
hundred to a few thousand BAMs and the reference only gets better.

---

## Install

```bash
git clone https://github.com/JMdeSteAgathe/intronaute.git
cd intronaute
python -m venv .venv
# Linux/macOS:  source .venv/bin/activate
# Windows:      .venv\Scripts\activate
pip install -e .
intronaute selftest
```

`selftest` unpacks a tiny embedded BAM and checks the counts against known
values. Run it first — two seconds, and it validates the whole BAM path on your
machine. It should print `BAM engine : fastbam` (or `pysam` if you have it).

Optional: `pip install pysam` on Linux/macOS for a slightly faster engine.
`bamnostic` is only a fallback for `.csi` indexes and is **not** recommended —
it was measured here to miss spliced reads whose junction reaches into the
queried region, which under-counts coverage.

---

## Quick start

```bash
# 1. build the target table from your minor-intron list + a GTF  (once, ~1 min)
intronaute prepare --introns examples/example_introns.bed \
                   --gtf gencode.v50.basic.annotation.gtf.gz

# 2. find the BAMs -> samplesheet.tsv  (one or many folders)
intronaute discover --bam-dir /data/run1/BAM /data/run2/BAM

# 3. count (parallel, cached, resumable)
intronaute quantify --workers 8

# 4. statistics, PCA, interactive HTML report
intronaute analyse
```

Or `intronaute run --bam-dir ...` for all four. Every option can also live in a
JSON config: `intronaute analyse --config myconfig.json` (command-line flags
take precedence).

### Input formats for the intron list

| format | coordinates | notes |
|---|---|---|
| `.bed`, `.bed.gz` | 0-based half-open | columns 4/6 used as gene/strand when present |
| `.tsv`, `.csv`, `.txt` | 1-based inclusive by default | header matched loosely: `chrom`/`chr`/`chromosome`, `start`/`chromStart`, `end`/`stop`… |
| `.ods`, `.xlsx` | 1-based inclusive by default | `--sheet` to pick a sheet |
| any of the above | — | a single `pos`/`region` column such as `chr1:1000-2000` also works |

Gene and strand columns are **optional**: when missing they are taken from the
host gene found in the GTF. Use `--coord-base 0` if your TSV is BED-like.

A quick sanity check when adopting a new list: the counts of usable introns
should match across formats. On a 1095-intron U12 list, BED / BED6 / TSV all
gave 1095 intervals and 1006 usable targets.

---

## What it measures

For every minor intron:

1. **coverage ratio** — intronic coverage normalised by the exons immediately
   flanking it (the last 300 bp of the upstream exon and the first 300 bp of the
   downstream one, not the whole gene, so it is insensitive to 3′ bias):
   `y = log2((intron depth + 0.5) / (flanking exon depth + 0.5))`
2. **junction counts** — correctly spliced reads (exact junction), reads
   crossing each splice boundary (direct evidence of retention), and reads using
   an alternative donor or acceptor.
3. **within-sample normalisation by major (U2) control introns.** Two
   size-matched major introns are picked from the *same gene* as each minor
   intron. Subtracting their median absorbs whatever makes intronic background
   vary between samples — RNA degradation, genomic-DNA carry-over, nuclear
   pre-mRNA, depth. Without it, a slightly degraded sample looks like a global
   positive.
4. **robust z per intron** against the cohort:
   `z = 0.6745 × (y_adj − median) / MAD`
5. **aggregation across the whole minor-intron set** — see below.

### A "hit"

An intron counts as a hit in a sample when **all three** hold:

- it is **measurable**: flanking-exon depth ≥ 8 (`--min-exon-depth`), else `NaN`
  and it counts in neither numerator nor denominator;
- `z ≥ 3` (`--z-threshold`);
- `y_adj − cohort median ≥ 0.5` (`--min-delta`), i.e. the intron/exon ratio is
  at least **1.41× (2^0.5)** above the cohort. This guards against introns so
  stable that a tiny MAD inflates `z` for a biologically meaningless difference.

`n_hits_control` applies **exactly the same rule** to the major introns, which
makes the ratio of the two directly interpretable.

### Per-sample scores

| column | meaning |
|---|---|
| `n_hits` | minor introns called — the default ranking column, and your gene list |
| `tail_p` / `tail_q` | **the main test.** Fisher, one-sided: hits on minor introns vs hits on control introns **of the same sample**. Self-calibrated, so it absorbs that sample's own quality and batch |
| `enrichment` | ratio of those two hit rates — readable effect size |
| `stouffer` | `Σz / √n`, sensitive to a small but coordinated shift |
| `MSRI` | winsorised mean of `z(minor)` minus that of `z(control)` |

The tail test is the right one here because the expected biology is **sparse**:
a few dozen strongly affected introns among a thousand. A location test such as
Mann-Whitney has almost no power in that regime — a few percent of shifted
values does not move the median.

**Always compare `n_hits` with `n_hits_control`.** Similar values mean generic
intronic noise, not a minor-spliceosome defect. One control sample in testing
showed 11 minor-intron hits but 18 control-intron hits: the tail test correctly
ranked it at `q = 1.0` where a raw hit count would have flagged it.

---

## Output

The report opens on an **interactive PCA** of the cohort: every sample is a
point, positioned by its minor-intron retention profile.

- type in the search box to highlight samples or whole runs (space-separated
  terms, substring match)
- scroll to zoom, drag to pan, click a point to pin its label
- recolour by `n_hits`, `-log10(tail_q)`, `log2(enrichment)`, `MSRI` or
  `stouffer`
- "pin top 10" labels the strongest samples

It is plain JavaScript on a canvas with no CDN, so the single `report.html`
works offline, handles thousands of points, and can be emailed as-is.

| file | contents |
|---|---|
| `results/report.html` | self-contained report, interactive PCA first |
| `results/sample_scores.tsv` | one row per sample, ranked |
| `results/hits.tsv.gz` | detailed rows for called introns only |
| `results/pca_coordinates.tsv` | PC1…PC6 + scores per sample (for your own plots) |
| `results/pca_loadings.tsv` | per-intron weights: which genes drive each axis |
| `results/pca_variance.tsv` | variance explained per component |
| `results/counts/*.tsv.gz` | raw per-sample counts — **this is the cache, do not delete** |

### Reading the results

1. **Quality control first.** A sample with low exonic depth, or an intronic
   background far above the rest, deserves caution — that is a library problem,
   not necessarily a splicing one.
2. **Check the run panel.** If runs separate, you have a batch effect. `tail_q`
   is largely immune to it since it calibrates within each sample.
3. **Then the PCA and the scores.** An affected sample sits to the right on PC1,
   with `tail_q < 0.05`, an enrichment of several tens, and `n_hits` far above
   its own `n_hits_control`.
4. **Then the events**, to know which genes to look at. Check that
   `n_EI_donor` / `n_EI_acceptor` (reads crossing each boundary) are non-zero:
   intronic coverage without boundary-crossing reads points to an overlapping
   transcript or a mapping artefact rather than true retention.

Per-intron calls are **much less reliable than the aggregate**. Treat them as
candidates to confirm (IGV, then RT-PCR).

`n_alt_5p` / `n_alt_3p` count spliced reads sharing only **one** of the intron's
two boundaries — they measure *alternative splicing*, not retention. A high
value is the classic signature of a nearby cryptic U2 site taking over when a
minor intron fails to splice; a high value with no intronic coverage usually
means your GTF disagrees with the expressed isoform. These four columns are
oriented **biologically** (5′ = donor), including on the minus strand.

---

## Validation

The statistics were developed against simulations with known truth. Two results
worth knowing about, because both were failures before they were fixes:

**Sparse signal, large cohort** — 198 samples over 14 runs with a deliberate
batch effect, 5 affected samples retaining 30 introns each out of 379 measurable
(8% of introns). Final: the 5 affected ranked 1–5, `tail_q < 0.05` for exactly
those 5, **0 false positives out of 193**, AUC 1.000.

Getting there exposed three design errors invisible at small scale:

1. a **symmetric trimmed mean** in `MSRI` discarded the top 10% of values — that
   is, exactly the signal, when a patient retains only 8% of introns. Replaced
   by a winsorised mean: AUC 0.85 → 1.00.
2. **Mann-Whitney** had no power on a sparse signal (no `q < 0.05` at all while
   the signal was overwhelming). Replaced by the Fisher tail test.
3. **PC1's sign is arbitrary** in an SVD and came out inverted. Now oriented so
   that further right always means more retention.

**PCA input matters** — running it on `y_adj` with per-intron rescaling gives
AUC 0.72; on the z-score matrix without rescaling, AUC 1.00. z already puts
every intron on a comparable scale, so re-standardising gives the noise of
stable introns the same weight as the signal. Values are winsorised so one wild
intron cannot capture PC1 on its own.

**The bundled BAM reader** produces counts identical to `pysam`, bit for bit,
across 3121 targets (max difference 0.0 on every numeric column), at comparable
speed.

---

## Performance

Across ~3200 targets the tool reads only about **6 Mb of genome** (versus 16 Mb
if it read from the upstream exon to the downstream one), in ~1900
coordinate-sorted regions:

- targeted access via the `.bai` index; the BAM is never scanned end to end;
- neighbouring targets are merged into one read region, so each BGZF block is
  decompressed once;
- for huge introns (500 kb+), only 2 kb windows either side of the two splice
  sites are measured — the middle is never read;
- one `int32` coverage array per region, so memory stays flat;
- **per-sample caching**: re-running recomputes nothing. If a network share
  drops mid-run, just relaunch the same command.

Analysis of 218 samples × 3121 targets: **41 s and 0.9 GB of RAM**. Matrices are
held as float32 numpy arrays rather than a long-format DataFrame, so thousands
of samples stay in a few hundred MB.

Reading over SMB/NFS, latency dominates — use `--workers 8` or copy BAMs
locally.

---

## Notes and caveats

- **An index is mandatory.** No `.bai` next to the `.bam` means no region
  queries; `discover` reports which files are affected.
- **Cohort size.** Below ~20 samples the median/MAD reference is unstable and
  scores are unreliable; the tool warns. Below 60 samples it switches to an
  exact leave-one-out so a sample cannot mask itself (`--loo-max`).
- **Strandedness** is detected per sample. Mixing stranded and unstranded runs
  in one cohort is not ideal; check the reported value during `quantify`.
- **Antisense coverage** (a 3′UTR of an overlapping gene) is removed twice,
  independently: exonic segments are masked out of the intron core regardless of
  strand, and antisense reads are filtered at quantification. `--mask-strand
  sense` masks only the same strand and recovers ~1.5% of measurable length, but
  is valid **only if every library is stranded** — `quantify` warns if it is not.
  Note that an overlapping gene on the *same* strand is not separable by strand
  at all; that is exactly what the annotation-level masking is for, and it stays
  active either way.
- **Overlapping paired-end mates** slightly inflate depth, but equally in intron
  and exons, so the ratio is essentially unaffected.
- Introns are dropped when too short after trimming, entirely exonic in another
  transcript, or without a usable reference exon (`usable` column of the target
  table).

---

## Layout

```
intronaute/
  fastbam.py    minimal pure-Python BAM/BAI reader (BGZF + index + CIGAR)
  bamio.py      engine abstraction, CIGAR, filters, strandedness
  regions.py    intron list + GTF -> targets (core, reference exons, controls)
  quantify.py   per-sample counting, merged regions, caching
  stats.py      normalisation, robust z, aggregation, PCA
  report.py     figures + interactive HTML report
  cli.py        command line
  selftest.py   self-check against an embedded mini-BAM
examples/       small BED and TSV intron lists showing accepted formats
```

Want to try other aggregation statistics? Start from `hits.tsv.gz` and
`pca_coordinates.tsv` rather than touching the counting step.

## License

MIT.
