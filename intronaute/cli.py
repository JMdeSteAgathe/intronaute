"""Command line interface.

  intronaute prepare   # build the target table from an intron list + a GTF
  intronaute discover  # scan BAM folders -> samplesheet.tsv
  intronaute quantify  # count in every BAM (parallel, cached)
  intronaute analyse   # statistics, PCA, interactive HTML report
  intronaute run       # all of the above
  intronaute selftest  # check that BAM reading works on this machine
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from typing import List, Optional

DEFAULT_OUT = "results"
SUFFIX_RE = re.compile(
    r"([-_.]?Aligned\.sortedByCoord\.out|\.sortedByCoord|\.sorted|\.dedup|"
    r"\.markdup|\.recal|\.coordSorted)+$", re.I)
RUN_SKIP = {"bam", "bams", "03_analyse", "analyse", "analysis", "alignment",
            "align", "star", "results", "out", "output", "mapping", "02_mapping"}


def _p(msg: str = "") -> None:
    print(msg, flush=True)


_HITS_RDS_R = r"""
hits <- read.delim(gzfile(Sys.getenv("HITS_TSV")), sep = "\t", check.names = FALSE,
                   stringsAsFactors = FALSE, quote = "", comment.char = "")
int_cols <- c("n_intron_reads", "n_EE", "n_EI_donor", "n_EI_acceptor",
              "n_alt_5p", "n_alt_3p", "intron_start", "intron_end", "intron_len")
for (col in intersect(int_cols, names(hits))) {
  storage.mode(hits[[col]]) <- "integer"
}
saveRDS(hits, Sys.getenv("HITS_RDS"), version = 3)
"""


def _rscript() -> str:
    import shutil
    found = shutil.which("Rscript")
    if found:
        return found
    home = os.path.expanduser("~")
    for cand in (os.path.join(home, ".pixi", "envs", "r-base", "bin", "Rscript"),
                 os.path.join(home, ".pixi", "bin", "Rscript")):
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    raise SystemExit("Rscript not found, so hits.rds was not written. "
                     "hits.tsv.gz is already in the output directory.")


def _write_hits_rds(tsv_gz: str, rds_path: str) -> None:
    """Same rows as hits.tsv.gz, stored as an R data frame."""
    import subprocess
    env = os.environ.copy()
    env["HITS_TSV"] = tsv_gz
    env["HITS_RDS"] = rds_path
    proc = subprocess.run(
        [_rscript(), "-e", _HITS_RDS_R],
        capture_output=True, text=True, env=env)
    if proc.returncode != 0:
        detail = (proc.stderr or proc.stdout or "").strip()
        raise SystemExit(f"Rscript failed while writing {rds_path}:\n{detail}")


def sample_name(bam: str) -> str:
    base = re.sub(r"\.bam$", "", os.path.basename(bam), flags=re.I)
    return SUFFIX_RE.sub("", base) or base


def run_name(bam_dir: str) -> str:
    parts = [x for x in re.split(r"[\\/]+", os.path.normpath(bam_dir)) if x]
    for x in reversed(parts):
        if x.lower() not in RUN_SKIP:
            return x
    return parts[-1] if parts else "run"


def _run_is_set(value) -> bool:
    """True when a samplesheet cell already carries a run label."""
    if value is None:
        return False
    try:
        import pandas as pd
        if pd.isna(value):
            return False
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    return bool(text) and text.lower() != "nan"


_BATCH_AS_RUN = ("batch", "batches")


def _rename_batch_column(df):
    """Use a batch column as run when the samplesheet has no run column.

    Accepted headers, ignoring case: batch, Batch, BATCH, batches, Batches,
    BATCHES. A sheet that already has run is left unchanged.
    """
    if "run" in df.columns:
        return df
    chosen = None
    for name in _BATCH_AS_RUN:
        for col in df.columns:
            if str(col).strip().lower() == name:
                chosen = col
                break
        if chosen is not None:
            break
    if chosen is None:
        return df
    _p(f"      samplesheet column '{chosen}' read as run")
    return df.rename(columns={chosen: "run"})


def _merge_existing_samplesheet(df, path: str):
    """Keep a user-defined run, and any extra columns, from an existing sheet.

    The directory name is the run only for samples whose run cell is missing.
    A batch column is treated as run when the sheet has no run column.
    """
    import pandas as pd

    if not path or not os.path.exists(path):
        return df, 0
    old = _rename_batch_column(pd.read_csv(path, sep="\t", dtype=str))
    if "sample" not in old.columns:
        return df, 0
    old = old.drop_duplicates("sample").set_index("sample")
    n_kept = 0
    if "run" in old.columns:
        kept = []
        for sample, directory_run in zip(df["sample"], df["run"]):
            if sample in old.index and _run_is_set(old.at[sample, "run"]):
                kept.append(str(old.at[sample, "run"]).strip())
                n_kept += 1
            else:
                kept.append(directory_run)
        df = df.copy()
        df["run"] = kept
    extra = [c for c in old.columns if c not in df.columns]
    if extra:
        df = df.merge(old[extra], left_on="sample", right_index=True, how="left")
    ordered = [c for c in ("sample", "bam", "run", "has_index", "size_GB")
               if c in df.columns]
    ordered += [c for c in df.columns if c not in ordered]
    return df[ordered], n_kept


# --------------------------------------------------------------------------- #
def cmd_prepare(a) -> None:
    from .regions import build_targets

    build_targets(introns_path=a.introns, gtf=a.gtf, out_tsv=a.targets,
                  exon_window=a.exon_window, intron_pad=a.intron_pad,
                  exon_pad=a.exon_pad, min_core_len=a.min_core_len,
                  min_exon_ref=a.min_exon_ref, max_core_window=a.max_core_window,
                  n_controls=a.n_controls, mask_strand=a.mask_strand,
                  coord_base=a.coord_base, sheet=a.sheet)


def cmd_discover(a) -> None:
    import pandas as pd

    from . import bamio

    dirs = getattr(a, "bam_dir", None)
    if not dirs:
        raise SystemExit("--bam-dir is required (on the command line or in --config)")
    rows = []
    for d in dirs:
        run = run_name(d)
        if not os.path.isdir(bamio.longpath(d)):
            _p(f"[warn] directory missing or unreadable: {d}")
            continue
        found = 0
        for root, _sub, files in os.walk(bamio.longpath(d)):
            for f in sorted(files):
                if not f.lower().endswith(".bam"):
                    continue
                p = os.path.join(root, f)
                rows.append(dict(sample=sample_name(p), bam=p, run=run,
                                 has_index=int(bamio.has_index(p)),
                                 size_GB=round(os.path.getsize(p) / 1e9, 2)))
                found += 1
        _p(f"  {found:>5} BAM   run={run}")

    if not rows:
        raise SystemExit("No BAM found. Check the paths and your access rights.")

    df = pd.DataFrame(rows)
    df, n_kept = _merge_existing_samplesheet(df, a.samplesheet)
    if n_kept:
        _p(f"  kept a user-defined run for {n_kept} sample(s); "
           f"the directory name fills only an empty run")
    dup = df["sample"].duplicated(keep=False)
    if dup.any():
        n0 = df.loc[dup, "sample"].nunique()
        df.loc[dup, "sample"] = df.loc[dup, "sample"] + "__" + df.loc[dup, "run"]
        still = df["sample"].duplicated(keep=False)
        if still.any():
            df.loc[still, "sample"] = (
                df.loc[still, "sample"] + "_"
                + (df.loc[still].groupby("sample").cumcount() + 1).astype(str))
        _p(f"\n[info] {n0} duplicated sample name(s) across runs: "
           f"disambiguated with the run name")
    assert not df["sample"].duplicated().any()

    noidx = df[df["has_index"] == 0]
    if len(noidx):
        _p(f"\n[WARNING] {len(noidx)} BAM without a .bai index -> unusable:")
        for s in noidx["sample"].head(8):
            _p(f"    {s}")
        if len(noidx) > 8:
            _p(f"    ... and {len(noidx) - 8} more")

    os.makedirs(os.path.dirname(os.path.abspath(a.samplesheet)) or ".", exist_ok=True)
    df.to_csv(a.samplesheet, sep="\t", index=False)
    _p(f"\n-> {a.samplesheet}  ({len(df)} samples, {df['run'].nunique()} runs, "
       f"{int(df['has_index'].sum())} usable, {df['size_GB'].sum():.0f} GB total)")


def _read_samplesheet(path: str):
    import pandas as pd

    if not os.path.exists(path):
        raise SystemExit(f"Samplesheet not found: {path}\n"
                         f"Run `intronaute discover` first.")
    df = _rename_batch_column(pd.read_csv(path, sep="\t"))
    for c in ("sample", "bam"):
        if c not in df.columns:
            raise SystemExit(f"Missing column '{c}' in {path}")
    if "run" not in df.columns:
        df["run"] = "run"
    if "has_index" in df.columns:
        df = df[df["has_index"].fillna(1).astype(int) == 1]
    return df.reset_index(drop=True)


# Paths and the index filter are inputs, not labels. Everything else on the
# samplesheet is copied onto the score table after scoring.
_SHEET_KEEP_OUT = {"sample", "bam", "has_index"}


def _score_annotations(ss, samples):
    """Samplesheet columns aligned to `samples`, for the score table only.

    These columns are not read by the median, MAD, or Fisher calculations.
    """
    import pandas as pd

    cols = [c for c in ss.columns if c not in _SHEET_KEEP_OUT]
    ann = ss.drop_duplicates("sample").set_index("sample").reindex(list(samples))
    if "run" not in ann.columns:
        ann["run"] = "run"
        cols = ["run", *cols]
    else:
        ann["run"] = ann["run"].fillna("run")
    return ann.loc[:, cols]


def _annotation_kind(series) -> str:
    """'num' for a continuous colour ramp, 'cat' for a discrete palette."""
    import pandas as pd

    if not pd.api.types.is_numeric_dtype(series):
        return "cat"
    n = int(series.notna().sum())
    nunq = int(series.nunique(dropna=True))
    # Few distinct numbers are groups (a 0/1 flag, a batch id), not a scale.
    if n and nunq <= 20 and nunq < 0.5 * n:
        return "cat"
    return "num"


def _rank_tables(scores, columns, n: int):
    """Top-n sample tables for the report rankings."""
    def take(column, ascending):
        ordered = scores.sort_values(
            [column, "stouffer"], ascending=[ascending, False],
            na_position="last")
        return ordered.head(n)[columns]

    return [
        ("n_hits", "n_hits", take("n_hits", False)),
        ("n_hits_control", "n_hits_control", take("n_hits_control", False)),
        ("tail_q", "tail_q", take("tail_q", True)),
        ("enrichment", "enrichment", take("enrichment", False)),
    ]


def _colour_menu(ann) -> list:
    menu = [
        {"id": "n_hits", "label": "colour: n_hits", "kind": "score"},
        {"id": "tail_q", "label": "colour: -log10(tail_q)", "kind": "score"},
        {"id": "enrichment", "label": "colour: log2(enrichment)", "kind": "score"},
        {"id": "MSRI", "label": "colour: MSRI", "kind": "score"},
        {"id": "stouffer", "label": "colour: stouffer", "kind": "score"},
    ]
    for c in ann.columns:
        menu.append({"id": str(c), "label": f"colour: {c}",
                     "kind": _annotation_kind(ann[c])})
    return menu


def _json_cell(v):
    import math

    import numpy as np
    import pandas as pd

    if v is None:
        return None
    try:
        if pd.isna(v):
            return None
    except (TypeError, ValueError):
        pass
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, (int, np.integer)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        x = float(v)
        return x if math.isfinite(x) else None
    return str(v)


def cmd_quantify(a) -> None:
    from concurrent.futures import ProcessPoolExecutor, as_completed

    from . import quantify as Q

    ss = _read_samplesheet(a.samplesheet)
    counts_dir = os.path.join(a.outdir, "counts")
    os.makedirs(counts_dir, exist_ok=True)
    params = dict(min_mapq=a.min_mapq, min_overhang=a.min_overhang,
                  drop_dup=not a.keep_duplicates, library=a.library,
                  backend=a.engine, max_reads_per_chunk=a.max_reads_per_region)

    jobs = []
    for r in ss.itertuples(index=False):
        out = Q.os.path.join(counts_dir, f"{r.sample}.tsv.gz")
        if not a.force and Q.is_cached(out, r.bam, params):
            continue
        if not os.path.exists(r.bam):
            _p(f"  [missing] {r.sample}: {r.bam}")
            continue
        jobs.append((r.sample, r.bam, a.targets, out,
                     dict(params, min_methods=0, kinds=None)))

    cached = len(ss) - len(jobs)
    if cached:
        _p(f"  {cached} sample(s) already in cache")
    if not jobs:
        _p("Nothing to compute.")
        return

    nw = a.workers if a.workers > 0 else max(1, min(8, (os.cpu_count() or 2) - 1))
    nw = min(nw, len(jobs))
    _p(f"Quantifying {len(jobs)} sample(s) with {nw} process(es)...")
    t0, done, failed = time.time(), 0, 0
    if nw == 1:
        for j in jobs:
            s, meta, err = Q._worker(j)
            done += 1
            failed += _log(s, meta, err, done, len(jobs), t0)
    else:
        with ProcessPoolExecutor(max_workers=nw) as ex:
            futs = {ex.submit(Q._worker, j): j[0] for j in jobs}
            for f in as_completed(futs):
                s, meta, err = f.result()
                done += 1
                failed += _log(s, meta, err, done, len(jobs), t0)
    _p(f"Done in {time.time()-t0:.0f}s ({failed} failed) -> {counts_dir}")


def _log(s, meta, err, done, total, t0) -> int:
    if err:
        _p(f"  [{done}/{total}] FAILED {s}\n{err}")
        return 1
    if done % max(1, total // 50) == 0 or done == total:
        eta = (time.time() - t0) / done * (total - done)
        _p(f"  [{done}/{total}] {s}: depth={meta['median_exon_depth']:.0f} "
           f"lib={meta['library']} {meta['seconds']:.0f}s  ETA {eta/60:.0f} min")
    return 0


# Every key config.yaml can set, in the same order as that file.
_CONFIG_GROUPS = (
    ("input", ("samplesheet", "bam_dir", "outdir", "introns", "gtf", "targets")),
    ("prepare", ("coord_base", "sheet", "exon_window", "intron_pad", "exon_pad",
                 "min_core_len", "min_exon_ref", "max_core_window", "n_controls",
                 "mask_strand", "force_prepare")),
    ("quantify", ("min_mapq", "min_overhang", "keep_duplicates", "library",
                  "engine", "workers", "max_reads_per_region", "force")),
    ("analyse", ("min_exon_depth", "min_frac_measurable", "z_threshold",
                 "min_delta", "hit_mode", "min_junction_reads",
                 "min_delta_junction", "junction_veto", "loo_max", "min_cohort",
                 "pca_clip", "pca_scale", "score_column", "top_samples",
                 "heatmap_samples", "heatmap_introns", "max_events_per_sample")),
)

# Count sidecars store the filters that produced the numbers in the report.
_COUNT_PARAM = {
    "min_mapq": "min_mapq",
    "min_overhang": "min_overhang",
    "library": "library",
    "engine": "backend",
}


def _param_text(value) -> str:
    if isinstance(value, (list, tuple)):
        return ", ".join(_param_text(v) for v in value)
    if value is None:
        return "None"
    return str(value)


def _flag_defaults() -> dict:
    defaults = {}
    parser = build_parser()
    for act in parser._actions:
        if not isinstance(act, argparse._SubParsersAction):
            continue
        for sp in act.choices.values():
            for aa in sp._actions:
                if aa.dest in ("help", "config") or aa.dest in defaults:
                    continue
                defaults[aa.dest] = aa.default
    return defaults


def _config_values(path) -> dict:
    if not path:
        return {}
    out = {}
    for key, value in _load_config_file(path).items():
        if not isinstance(key, str) or key.startswith("_"):
            continue
        out[key.replace("-", "_")] = value
    return out


def _recorded_text(counts) -> Optional[str]:
    """One recorded value, or each value with the number of samples that used it."""
    if not counts:
        return None
    if len(counts) == 1:
        return _param_text(next(iter(counts)))
    return ", ".join(f"{_param_text(v)} ({n})" for v, n in counts.most_common())


def _count_param_values(outdir: str) -> dict:
    """Quantify settings recorded beside the cached counts.

    ``keep_duplicates`` is the inverse of the recorded ``drop_dup``.
    """
    from collections import Counter

    counts = os.path.join(outdir, "counts")
    if not os.path.isdir(counts):
        return {}
    found = {k: Counter() for k in ("min_mapq", "min_overhang", "keep_duplicates",
                                    "library", "engine")}
    for name in os.listdir(counts):
        if not name.endswith(".json"):
            continue
        try:
            with open(os.path.join(counts, name), encoding="utf-8") as fh:
                meta = json.load(fh)
        except (OSError, json.JSONDecodeError):
            continue
        for dest, src in _COUNT_PARAM.items():
            if src in meta:
                found[dest][meta[src]] += 1
        if "drop_dup" in meta:
            found["keep_duplicates"][not bool(meta["drop_dup"])] += 1
    return {k: text for k, vals in found.items()
            if (text := _recorded_text(vals)) is not None}


def _parameter_table(a):
    """Every configuration key, with the value that applies to this report."""
    import pandas as pd

    defaults = _flag_defaults()
    cfg = _config_values(getattr(a, "config", None))
    recorded = _count_param_values(a.outdir)
    rows = [("config", _param_text(getattr(a, "config", None)), "command"),
            ("cmd", _param_text(getattr(a, "cmd", None)), "command")]
    for stage, keys in _CONFIG_GROUPS:
        for key in keys:
            if key in recorded:
                value = recorded[key]
            elif hasattr(a, key):
                value = getattr(a, key)
            elif key in cfg:
                value = cfg[key]
            else:
                value = defaults.get(key)
            rows.append((key, _param_text(value), stage))
    return pd.DataFrame(rows, columns=["parameter", "value", "stage"])


def cmd_analyse(a) -> None:
    import numpy as np
    import pandas as pd

    from . import layout as L
    from . import report as R
    from . import stats as S

    ss = _read_samplesheet(a.samplesheet)
    counts_dir = os.path.join(a.outdir, "counts")

    _p("[1/5] reading counts...")
    _p(f"      samplesheet: {os.path.abspath(a.samplesheet)} ({len(ss)} samples)")
    import glob

    on_disk = {os.path.basename(f)[:-7]
               for f in glob.glob(os.path.join(counts_dir, "*.tsv.gz"))}
    extra = on_disk - set(ss["sample"])
    if extra:
        _p(f"      [WARNING] {len(on_disk)} samples quantified in {counts_dir} "
           f"but only {len(ss)} listed in the samplesheet; {len(extra)} will be "
           f"IGNORED (e.g. {sorted(extra)[:3]}).")
        _p("      -> stale samplesheet: re-run `discover`. Do NOT delete "
           "counts/, it is the cache.")

    mode = a.hit_mode
    cnt = S.load_counts(counts_dir, list(ss["sample"]))
    intron, exon, target_ids, samples = cnt.intron, cnt.exon, cnt.target_ids, cnt.samples
    has_junc = cnt.ei is not None
    if mode != "coverage" and not has_junc:
        raise SystemExit(f"--hit-mode {mode} needs junction counts, which this "
                         f"counts cache lacks. Re-run `quantify --force`.")
    meta = S.load_meta(a.targets).reindex(target_ids)
    # Labels stay out of the reference. They are attached to `scores` below,
    # after hits and the Fisher test have already been computed.
    ann = _score_annotations(ss, samples)
    _p(f"      {len(samples)} samples x {len(target_ids)} targets")
    if len(samples) < a.min_cohort:
        _p(f"      [WARNING] only {len(samples)} samples. This tool scores each "
           f"sample against the cohort itself; below ~{a.min_cohort} samples the "
           f"median/MAD reference is unstable and the scores are unreliable.")

    kind = meta["kind"].to_numpy()
    is_minor = kind == "minor"
    is_control = kind == "control"

    _p(f"[2/5] normalising and scoring (hit mode: {mode})...")
    y = S.build_y(intron, exon, a.min_exon_depth)
    measurable = np.isfinite(y).mean(axis=0) >= a.min_frac_measurable
    y[:, ~measurable] = np.nan
    y_adj, offsets = S.normalise_by_controls(y, is_control & measurable)
    z_cov = S.robust_z(y_adj, loo_max=a.loo_max)
    d_cov = S.delta_from_cohort(y_adj)

    z_j = d_j = z_comb = theta = None
    rho = float("nan")
    measurable_j = np.zeros_like(measurable)
    offsets_j = np.zeros(len(samples), dtype=np.float32)
    if has_junc:
        yj, theta = S.build_y_junction(cnt.ei, cnt.split, a.min_junction_reads)
        measurable_j = np.isfinite(yj).mean(axis=0) >= a.min_frac_measurable
        yj[:, ~measurable_j] = np.nan
        yj_adj, offsets_j = S.normalise_by_controls(yj, is_control & measurable_j)
        z_j = S.robust_z(yj_adj, loo_max=a.loo_max,
                         mad_floor_q=S.JUNCTION_MAD_FLOOR_Q)
        d_j = S.delta_from_cohort(yj_adj)
        z_comb, rho = S.combine_z(z_cov, z_j, is_control & measurable & measurable_j)

    # the analysed set: an intron is in if the chosen mode can measure it
    if mode == "coverage":
        usable = measurable
    elif mode == "junction":
        usable = measurable_j
    elif mode == "both":
        usable = measurable & measurable_j
    else:
        usable = measurable | measurable_j
    z = {"coverage": z_cov, "junction": z_j}.get(mode, z_comb)
    y_rep = y_adj if mode != "junction" else yj_adj   # for hits.tsv 'y_adj'

    hit = S.call_hits(mode, a.z_threshold, z_cov, d_cov, a.min_delta,
                      z_j=z_j, d_j=d_j, min_delta_j=a.min_delta_junction,
                      z_comb=z_comb, junction_veto=a.junction_veto)
    scores, hit = S.sample_scores(z, y_rep, is_minor & usable,
                                  is_control & usable, samples,
                                  z_thr=a.z_threshold, hit=hit)
    if has_junc:
        # single-measure counts side by side, same thresholds, for comparison
        for name, m in (("coverage", "coverage"), ("junction", "junction")):
            h = S.call_hits(m, a.z_threshold, z_cov, d_cov, a.min_delta,
                            z_j=z_j, d_j=d_j, min_delta_j=a.min_delta_junction)
            scores[f"n_hits_{name}_only_rule"] = h[:, is_minor & usable].sum(axis=1)
        scores["n_minor_junction_measurable"] = np.isfinite(
            z_j[:, is_minor & measurable_j]).sum(axis=1)
    overlap = [c for c in ann.columns if c in scores.columns]
    if overlap:
        _p(f"      [WARNING] samplesheet column(s) {overlap} match score "
           f"columns; the score values were kept")
        ann = ann.drop(columns=overlap)
    runs = ann["run"] if "run" in ann.columns else pd.Series(
        ["run"] * len(samples), index=list(samples))
    for c in ann.columns:
        scores[c] = ann[c].to_numpy()
    scores = scores.sort_values(["n_hits", "stouffer"], ascending=False
                                ).reset_index(drop=True)
    measurable = usable
    n_minor = int((is_minor & measurable).sum())
    n_ctrl = int((is_control & measurable).sum())
    _p(f"      {n_minor} minor / {n_ctrl} control introns measurable")
    if has_junc:
        _p(f"      junctions measurable: {int((is_minor & measurable_j).sum())} "
           f"minor / {int((is_control & measurable_j).sum())} control "
           f"(>= {a.min_junction_reads:g} reads at the splice sites); "
           f"corr(z_cov, z_junction) on controls = {rho:.2f}")

    _p("[3/5] PCA...")
    sc, load, expl, shape = S.pca(z, samples, is_minor & measurable,
                                  n_comp=min(6, len(samples)), clip=a.pca_clip,
                                  scale=a.pca_scale)
    _p(f"      {shape[1]} introns; PC1={100*expl[0]:.1f}% PC2={100*expl[1]:.1f}%")

    _p("[4/5] writing tables...")
    dirs = L.prepare(a.outdir)
    sr = scores.set_index("sample")
    pca_cols = ["run", "n_hits", "n_hits_control", "enrichment",
                "tail_q", "MSRI", "stouffer"]
    pca_cols += [c for c in ann.columns if c not in pca_cols and c in sr.columns]
    acp = sc.join(sr[pca_cols])
    acp.index.name = "sample"
    acp.to_csv(os.path.join(dirs["pca"], "pca_coordinates.tsv"), sep="\t")
    pd.DataFrame({"component": list(sc.columns),
                  "variance_explained": [float(v) for v in expl]}).to_csv(
        os.path.join(dirs["pca"], "pca_variance.tsv"), sep="\t", index=False)
    kept = np.where(is_minor & measurable)[0][:load.shape[1]]
    pd.DataFrame(load.T, columns=list(sc.columns),
                 index=[target_ids[i] for i in kept]).assign(
        gene_name=[meta.iloc[i]["gene_name"] for i in kept]).to_csv(
        os.path.join(dirs["pca"], "pca_loadings.tsv"), sep="\t")
    scores.to_csv(os.path.join(dirs["ranked"], "sample_scores.tsv"),
                  sep="\t", index=False)

    extra = {}
    if has_junc:
        extra = dict(z_coverage=z_cov, z_junction=z_j, theta=theta)
        if mode in ("combined", "both"):
            extra["z_combined"] = z_comb
    events = S.hit_details(counts_dir, samples, target_ids, hit, z, y_rep, meta,
                           is_minor, max_per_sample=a.max_events_per_sample,
                           extra=extra)
    if len(events):
        hits_gz = os.path.join(dirs["ranked"], "hits.tsv.gz")
        events.to_csv(hits_gz, sep="\t", index=False, compression="gzip")
        _write_hits_rds(hits_gz, os.path.join(dirs["ranked"], "hits.rds"))

    _p("[5/5] figures and report...")
    with np.errstate(all="ignore"):
        qc = pd.DataFrame(dict(
            sample=samples, run=runs.to_numpy(),
            median_exon_depth=np.nanmedian(exon, axis=1),
            frac_measurable=(exon >= a.min_exon_depth).mean(axis=1),
            control_background=offsets))
    if has_junc:
        qc["control_junction_background"] = offsets_j
    f_scores = R.fig_score_distribution(scores, dirs["ranked"], col=a.score_column)
    hm_kw = dict(n_samples=a.heatmap_samples, n_introns=a.heatmap_introns)
    # Heatmap color is minor-intron z minus that sample's own control-intron z,
    # so a sample whose ordinary introns are retained is less red.
    ctrl_z = (scores.set_index("sample")["median_z_control"]
              .reindex(list(samples)).to_numpy(dtype=float))
    ctrl_z = np.where(np.isfinite(ctrl_z), ctrl_z, 0.0)
    z_hm = z - ctrl_z[:, None]
    f_hm = R.fig_heatmap(
        z_hm, samples, target_ids, is_minor & measurable, meta, scores,
        dirs["heat"], rank_col="n_hits", filename="heatmap_n_hits.png",
        rank_label="n_hits", **hm_kw)
    f_hm_control = R.fig_heatmap(
        z_hm, samples, target_ids, is_minor & measurable, meta, scores,
        dirs["heat"], rank_col="n_hits_control",
        filename="heatmap_n_hits_control.png", rank_label="n_hits_control",
        **hm_kw)
    f_hm_q = R.fig_heatmap(
        z_hm, samples, target_ids, is_minor & measurable, meta, scores,
        dirs["heat"], rank_col="tail_q", ascending=True,
        filename="heatmap_tail_q.png", rank_label="smallest tail_q", **hm_kw)
    f_qc = R.fig_qc(qc, dirs["qc"])
    gene_sets = R.minor_hit_genes(hit, samples, target_ids, meta, is_minor)
    upset_items = []
    upset_tables = []
    singleton_rows = []
    groups = R.upset_groups(scores)
    diagnosis_singletons = {ids[0] for gid, _lab, ids in groups
                            if len(ids) == 1 and gid != "spliceopathy"}

    def drop_upset_png(gid):
        stale = os.path.join(dirs["shared"], L.upset_filename(gid))
        if os.path.exists(stale):
            os.remove(stale)

    diag_labels = (
        R._label_spellings(scores["diagnosis"])
        if "diagnosis" in scores.columns else {})

    def diagnosis_of(sample):
        if "diagnosis" not in scores.columns:
            return ""
        hit_rows = scores.loc[scores["sample"].astype(str) == str(sample),
                              "diagnosis"]
        if hit_rows.empty:
            return ""
        return diag_labels.get(R._cell_text(hit_rows.iloc[0]).casefold(), "")

    for gid, label, ids in groups:
        if len(ids) < 2:
            drop_upset_png(gid)
            sample = ids[0]
            if gid == "spliceopathy" and sample in diagnosis_singletons:
                continue
            kind = R.singleton_kind(gid, label, diagnosis_of(sample))
            genes = sorted(gene_sets.get(sample, ()))
            upset_items.append((
                gid, f"{sample} ({kind})", "", "",
                R.singleton_gene_html(sample, kind, genes)))
            singleton_rows.extend(
                {"sample": sample, "label": kind, "gene": g} for g in genes)
            continue
        rows = R.upset_intersections(gene_sets, ids)
        upset_name = L.upset_filename(gid)
        path = R.fig_upset(
            gene_sets, ids, dirs["shared"], upset_name,
            f"Genes with minor-specific retention shared by {label} "
            f"(n={len(ids)})")
        table = R.upset_shared_table(rows)
        extra = (R._table(table) if len(table)
                 else "<div class='note'>No gene is shared by two or more "
                      "of these individuals.</div>")
        note = ("Each row of dots is one person. A column is one combination "
                "of people. The bar above it is the number of genes with a "
                "minor-intron hit, and no hit on the paired control, in "
                "exactly those people. Blue combinations are shared by two "
                "or more people. Grey combinations are private to one person. "
                "Intersection bars run from the most genes on the left to the "
                "fewest on the right. The bars on the left are how many genes "
                "each person has.")
        upset_items.append((gid, label, path, note, extra))
        if len(table):
            block = table.copy()
            block.insert(0, "group", label)
            upset_tables.append(block)
    if upset_tables:
        pd.concat(upset_tables, ignore_index=True).to_csv(
            os.path.join(dirs["shared"], "upset_shared_genes.tsv"),
            sep="\t", index=False)
    else:
        stale_tsv = os.path.join(dirs["shared"], "upset_shared_genes.tsv")
        if os.path.exists(stale_tsv):
            os.remove(stale_tsv)
    if singleton_rows:
        pd.DataFrame(singleton_rows).to_csv(
            os.path.join(dirs["shared"], "singleton_minor_hit_genes.tsv"),
            sep="\t", index=False)
    else:
        stale_single = os.path.join(dirs["shared"], "singleton_minor_hit_genes.tsv")
        if os.path.exists(stale_single):
            os.remove(stale_single)

    def col(name):
        return sr[name].reindex(samples).to_numpy()

    def clean(v):
        return None if (v is None or not np.isfinite(v)) else round(float(v), 4)

    ann_cols = [c for c in ann.columns if c in scores.columns]
    colours = _colour_menu(ann[ann_cols] if ann_cols else ann.iloc[:, 0:0])
    run_names = list(dict.fromkeys(str(r) for r in runs.tolist()))
    run_pos = {name: i for i, name in enumerate(run_names)}

    def qc_val(column, i, nd):
        v = qc.iloc[i][column]
        try:
            v = float(v)
        except (TypeError, ValueError):
            return None
        if not np.isfinite(v):
            return None
        return round(v, nd)

    points = []
    for i, s in enumerate(samples):
        points.append(dict(
            n=s, sh=(s if len(s) <= 22 else s[:21] + "\u2026"),
            x=round(float(sc.iloc[i, 0]), 3),
            y=round(float(sc.iloc[i, 1]), 3),
            h=int(col("n_hits")[i]), hc=int(col("n_hits_control")[i]),
            e=clean(col("enrichment")[i]), q=clean(col("tail_q")[i]),
            m=clean(col("MSRI")[i]), st=clean(col("stouffer")[i]),
            r=str(runs.iloc[i]),
            ri=run_pos[str(runs.iloc[i])],
            depth=qc_val("median_exon_depth", i, 2),
            frac=qc_val("frac_measurable", i, 4),
            bg=qc_val("control_background", i, 3),
            meta={c: _json_cell(ann.iloc[i][c]) for c in ann_cols}))

    show = ["sample", "run", "n_hits", "n_hits_control", "enrichment", "tail_q",
            "MSRI", "stouffer", "n_minor_measurable"]
    if has_junc:
        show += ["n_hits_coverage_only_rule", "n_hits_junction_only_rule"]
    show += [c for c in ann_cols if c not in show]
    run_tbl = pd.DataFrame()
    if runs.nunique() > 1:
        run_tbl = (scores.groupby("run")
                   .agg(n=("sample", "size"),
                        n_hits_median=("n_hits", "median"),
                        n_hits_max=("n_hits", "max"),
                        control_hits_median=("n_hits_control", "median"))
                   .sort_values("n_hits_median", ascending=False).reset_index())

    ev_cols = [c for c in ["sample", "gene_name", "chrom", "intron_start",
                           "intron_end", "strand", "z", "z_coverage",
                           "z_junction", "theta", "y_adj",
                           "intron_depth_mean", "exon_depth_mean", "n_EE",
                           "n_EI_donor", "n_EI_acceptor", "n_alt_5p", "n_alt_3p"]
               if not len(events) or c in events.columns]
    top_ev = (events.head(0) if not len(events) else
              events.sort_values("z", ascending=False, na_position="last")
              .head(60)[ev_cols])

    xlab = f"PC1 ({100 * expl[0]:.1f} %)"
    ylab = f"PC2 ({100 * expl[1]:.1f} %)"
    L.write_colour_plots(scores, sc, qc, colours, xlab, ylab, dirs)
    L.remove_legacy(a.outdir)

    R.write_html(os.path.join(a.outdir, "report.html"), dict(
        n_samples=len(samples), n_minor=n_minor, n_control=n_ctrl, hit_mode=mode,
        points=points, colours=colours, run_names=run_names,
        xlab=xlab, ylab=ylab,
        top_samples=scores.head(a.top_samples)[show],
        rankings=_rank_tables(scores, show, a.top_samples),
        fig_scores=f_scores, fig_heatmap=f_hm, fig_qc=f_qc,
        fig_heatmap_control=f_hm_control, fig_heatmap_q=f_hm_q,
        upset=upset_items,
        upset_catalog=R.upset_catalog(scores, gene_sets),
        top_events=top_ev, run_table=run_tbl,
        params=_parameter_table(a)))

    _p("")
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        _p(scores.head(15)[show].to_string(index=False))
    _p("")
    _p(f"Report : {os.path.join(a.outdir, 'report.html')}")
    _p(f"Tables : {dirs['ranked']}/sample_scores.tsv, "
       f"{dirs['pca']}/pca_coordinates.tsv, "
       f"{dirs['ranked']}/hits.tsv.gz, {dirs['ranked']}/hits.rds")


def cmd_run(a) -> None:
    if a.force_prepare or not os.path.exists(a.targets):
        _p("### prepare ###")
        cmd_prepare(a)
    if getattr(a, "bam_dir", None) and (a.force or not os.path.exists(a.samplesheet)):
        _p("\n### discover ###")
        cmd_discover(a)
    _p("\n### quantify ###")
    cmd_quantify(a)
    _p("\n### analyse ###")
    cmd_analyse(a)


def cmd_selftest(a) -> None:
    from .selftest import run_selftest

    run_selftest()


# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="intronaute",
        description="Minor-intron (U12) retention detector for RNA-seq BAM cohorts.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--config", help="YAML or JSON file of parameters. Every key "
                                   "is optional and replaces the default for "
                                   "that flag. Command-line flags win.")
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp):
        sp.add_argument("--config", help="YAML or JSON file of parameters. Every "
                                       "key is optional and replaces the "
                                       "default for that flag. Command-line "
                                       "flags win.")
        sp.add_argument("--outdir", default=DEFAULT_OUT)
        sp.add_argument("--targets", default=os.path.join("annotation",
                                                          "targets.tsv.gz"))
        sp.add_argument("--samplesheet", default="samplesheet.tsv")

    def prep_args(sp):
        sp.add_argument("--introns", help="minor-intron list: BED, TSV, CSV, "
                                          "ODS or XLSX")
        sp.add_argument("--gtf", help="GENCODE/Ensembl GTF (may be .gz)")
        sp.add_argument("--coord-base", default="auto", choices=["auto", "0", "1"],
                        help="coordinate base of the intron list; BED is always "
                             "0-based, other formats default to 1-based")
        sp.add_argument("--sheet", default=None, help="sheet name for ODS/XLSX")
        sp.add_argument("--exon-window", type=int, default=300)
        sp.add_argument("--intron-pad", type=int, default=12)
        sp.add_argument("--exon-pad", type=int, default=5)
        sp.add_argument("--min-core-len", type=int, default=40)
        sp.add_argument("--min-exon-ref", type=int, default=50)
        sp.add_argument("--max-core-window", type=int, default=2000)
        sp.add_argument("--n-controls", type=int, default=2)
        sp.add_argument("--mask-strand", default="both", choices=["both", "sense"])
        sp.add_argument("--force-prepare", action="store_true")

    def quant_args(sp):
        sp.add_argument("--min-mapq", type=int, default=10)
        sp.add_argument("--min-overhang", type=int, default=6)
        sp.add_argument("--keep-duplicates", action="store_true")
        sp.add_argument("--library", default="auto",
                        choices=["auto", "reverse", "forward", "unstranded"])
        sp.add_argument("--engine", default="auto",
                        choices=["auto", "fastbam", "pysam", "bamnostic"])
        sp.add_argument("--workers", type=int, default=0, help="0 = auto")
        sp.add_argument("--max-reads-per-region", type=int, default=300000)
        sp.add_argument("--force", action="store_true", help="ignore the cache")

    def ana_args(sp):
        sp.add_argument("--min-exon-depth", type=float, default=8.0)
        sp.add_argument("--min-frac-measurable", type=float, default=0.5)
        sp.add_argument("--z-threshold", type=float, default=3.0)
        sp.add_argument("--min-delta", type=float, default=0.5)
        sp.add_argument("--hit-mode", default="coverage",
                        choices=["coverage", "junction", "combined", "both"],
                        help="what defines a hit: intronic coverage ratio, "
                             "splice-site junction ratio (FRASER-like theta), "
                             "their correlation-corrected combination, or both "
                             "rules at once")
        sp.add_argument("--min-junction-reads", type=float, default=10,
                        help="junction ratio is NaN below this many reads "
                             "(crossing + spliced) at the two splice sites")
        sp.add_argument("--min-delta-junction", type=float, default=0.5,
                        help="effect-size guard for the junction ratio, log2 "
                             "odds above the cohort median")
        sp.add_argument("--junction-veto", type=float, default=1.0,
                        help="combined mode: never call an intron whose "
                             "junctions are measurable with z_junction below "
                             "this. Use -inf to disable")
        sp.add_argument("--loo-max", type=int, default=60,
                        help="exact leave-one-out below this cohort size")
        sp.add_argument("--min-cohort", type=int, default=20)
        sp.add_argument("--pca-clip", type=float, default=8.0)
        sp.add_argument("--pca-scale", action="store_true")
        sp.add_argument("--score-column", default="n_hits")
        sp.add_argument("--top-samples", type=int, default=25)
        sp.add_argument("--heatmap-samples", type=int, default=30)
        sp.add_argument("--heatmap-introns", type=int, default=0,
                        help="minor introns as heatmap columns; 0 means every "
                             "measurable minor intron")
        sp.add_argument("--max-events-per-sample", type=int, default=200)

    sp = sub.add_parser("prepare", help="build the target table")
    common(sp), prep_args(sp)
    sp.set_defaults(func=cmd_prepare)

    sp = sub.add_parser("discover", help="scan BAM folders")
    common(sp)
    sp.add_argument("--bam-dir", nargs="+")
    sp.set_defaults(func=cmd_discover)

    sp = sub.add_parser("quantify", help="count in the BAMs")
    common(sp), quant_args(sp)
    sp.set_defaults(func=cmd_quantify)

    sp = sub.add_parser("analyse", aliases=["analyze"],
                        help="statistics, PCA, report")
    common(sp), ana_args(sp)
    sp.set_defaults(func=cmd_analyse)

    sp = sub.add_parser("run", help="prepare + discover + quantify + analyse")
    common(sp), prep_args(sp), quant_args(sp), ana_args(sp)
    sp.add_argument("--bam-dir", nargs="+")
    sp.set_defaults(func=cmd_run)

    sp = sub.add_parser("selftest", help="check BAM reading on this machine")
    sp.set_defaults(func=cmd_selftest)
    return p


def _all_dests(parser) -> set:
    out = set()
    for act in parser._actions:
        if isinstance(act, argparse._SubParsersAction):
            for s in act.choices.values():
                out |= {x.dest for x in s._actions}
        else:
            out.add(act.dest)
    return out


def _dests_given(parser, argv) -> set:
    opt = {}
    for act in parser._actions:
        if isinstance(act, argparse._SubParsersAction):
            for s in act.choices.values():
                for aa in s._actions:
                    for o in aa.option_strings:
                        opt[o] = aa.dest
        else:
            for o in act.option_strings:
                opt[o] = act.dest
    return {opt[t.split("=")[0]] for t in argv if t.split("=")[0] in opt}


class _BareKey:
    """A YAML key that was written with no value and no following list."""


def _strip_yaml_comment(line: str) -> str:
    out = []
    quote = None
    for ch in line:
        if quote:
            out.append(ch)
            if ch == quote:
                quote = None
        elif ch in ("'", '"'):
            quote = ch
            out.append(ch)
        elif ch == "#":
            break
        else:
            out.append(ch)
    return "".join(out).rstrip()


def _parse_yaml_scalar(raw: str):
    s = raw.strip()
    if len(s) >= 2 and s[0] == s[-1] and s[0] in ("'", '"'):
        return s[1:-1]
    low = s.lower()
    if low in ("true", "yes", "y", "on"):
        return True
    if low in ("false", "no", "n", "off"):
        return False
    if low in ("null", "~", "none"):
        return None
    if low in ("inf", "+inf", ".inf"):
        return float("inf")
    if low in ("-inf", "-.inf"):
        return float("-inf")
    if low == "nan":
        return float("nan")
    try:
        if any(c in s for c in ".eE"):
            return float(s)
        return int(s)
    except ValueError:
        return s


def _as_path_list(value) -> List[str]:
    """One path, a block list, or a flow list ``[a, b]``."""
    if value is None or isinstance(value, _BareKey):
        return []
    if isinstance(value, (list, tuple)):
        out: List[str] = []
        for item in value:
            out.extend(_as_path_list(item))
        return out
    text = str(value).strip()
    if text.startswith("[") and text.endswith("]"):
        return _split_flow_list(text[1:-1])
    return [text] if text else []


def _split_flow_list(inner: str) -> List[str]:
    parts: List[str] = []
    buf: List[str] = []
    quote = None
    for ch in inner:
        if quote:
            if ch == quote:
                quote = None
            else:
                buf.append(ch)
        elif ch in ("'", '"'):
            quote = ch
        elif ch == ",":
            part = "".join(buf).strip()
            if part:
                parts.append(part)
            buf = []
        else:
            buf.append(ch)
    part = "".join(buf).strip()
    if part:
        parts.append(part)
    return parts


def _fold_bam_dirs(cfg: dict) -> dict:
    """Accept bam_dir and bam_dirs. Both become the bam_dir list."""
    if "bam_dir" not in cfg and "bam_dirs" not in cfg:
        return cfg
    seen: List[str] = []
    have = set()
    for key in ("bam_dir", "bam_dirs"):
        if key not in cfg:
            continue
        for path in _as_path_list(cfg[key]):
            if path not in have:
                have.add(path)
                seen.append(path)
    cfg = dict(cfg)
    cfg.pop("bam_dir", None)
    cfg.pop("bam_dirs", None)
    if seen:
        cfg["bam_dir"] = seen
    return cfg


def _load_yaml_mapping(text: str) -> dict:
    """Mapping subset used by config.yaml: scalars, comments, and lists.

    No nested mappings. A bare ``key:`` must be followed by ``- item`` lines.
    """
    data = {}
    current = None
    for lineno, raw in enumerate(text.splitlines(), 1):
        line = _strip_yaml_comment(raw)
        if not line.strip():
            continue
        indent = len(line) - len(line.lstrip())
        body = line.strip()
        if body.startswith("-"):
            if current is None:
                raise SystemExit(f"{lineno}: list item is not under a key: {raw}")
            item = body[1:].strip()
            if item == "":
                raise SystemExit(f"{lineno}: empty list item")
            bucket = data.get(current)
            if not isinstance(bucket, list):
                data[current] = []
            data[current].append(_parse_yaml_scalar(item))
            continue
        if ":" not in body or indent:
            raise SystemExit(f"{lineno}: expected 'key: value' in config: {raw}")
        key, _, val = body.partition(":")
        key = key.strip()
        if not key:
            raise SystemExit(f"{lineno}: empty key in config")
        val = val.strip()
        if val == "":
            data[key] = _BareKey()
            current = key
        else:
            data[key] = _parse_yaml_scalar(val)
            current = None
    bare = [k for k, v in data.items() if isinstance(v, _BareKey)]
    if bare:
        raise SystemExit("Config key has no value: " + ", ".join(bare))
    return data


def _load_config_file(path: str) -> dict:
    try:
        with open(path, encoding="utf-8") as fh:
            text = fh.read()
    except OSError as exc:
        raise SystemExit(f"Could not read config {path}: {exc}") from exc
    ext = os.path.splitext(path)[1].lower()
    if ext in (".yaml", ".yml"):
        cfg = _load_yaml_mapping(text)
    elif ext == ".json":
        try:
            cfg = json.loads(text)
        except json.JSONDecodeError as exc:
            raise SystemExit(f"Could not parse {path} as JSON: {exc}") from exc
    else:
        stripped = text.lstrip()
        if stripped.startswith("{"):
            try:
                cfg = json.loads(text)
            except json.JSONDecodeError as exc:
                raise SystemExit(f"Could not parse {path} as JSON: {exc}") from exc
        else:
            cfg = _load_yaml_mapping(text)
    if not isinstance(cfg, dict):
        raise SystemExit(f"Config {path} must be a mapping of flag names to values")
    return cfg


def _subcommand_actions(parser, cmd: str) -> dict:
    for act in parser._actions:
        if isinstance(act, argparse._SubParsersAction) and cmd in act.choices:
            return {a.dest: a for a in act.choices[cmd]._actions}
    return {}


def _coerce_config_value(dest: str, action, value):
    if isinstance(action, argparse._StoreTrueAction):
        if isinstance(value, bool):
            return value
        if isinstance(value, (int, float)) and value in (0, 1):
            return bool(value)
        raise SystemExit(
            f"Config key {dest} must be true or false, got {value!r}")
    if action.nargs in ("+", "*"):
        if value is None:
            return None
        items = [value] if isinstance(value, str) else list(value)
        if action.type:
            items = [action.type(v) for v in items]
        return items
    if value is not None and action.type is not None:
        try:
            value = action.type(value)
        except (TypeError, ValueError) as exc:
            raise SystemExit(
                f"Config key {dest} is not a valid {action.type.__name__}: "
                f"{value!r}") from exc
    if action.choices is not None and value not in action.choices:
        if str(value) in action.choices:
            value = str(value)
        else:
            allowed = ", ".join(str(c) for c in action.choices)
            raise SystemExit(
                f"Config key {dest}={value!r} is not one of: {allowed}")
    return value


def _apply_config(parser, args, argv) -> None:
    """Replace defaults with keys present in --config. Absent keys stay default.

    A flag also given on the command line is left as the command line set it.
    Keys that belong to another subcommand are ignored. Keys starting with
    ``_`` are comments.
    """
    path = getattr(args, "config", None)
    if not path:
        return
    cfg = _fold_bam_dirs(_load_config_file(path))
    actions = _subcommand_actions(parser, args.cmd)
    known = _all_dests(parser)
    given = _dests_given(parser, argv)
    for key, value in cfg.items():
        if not isinstance(key, str) or key.startswith("_"):
            continue
        dest = key.replace("-", "_")
        if dest in ("config", "cmd", "func") or dest in given:
            continue
        action = actions.get(dest)
        if action is None:
            if dest not in known:
                _p(f"[warn] unknown config parameter ignored: {dest}")
            continue
        setattr(args, dest, _coerce_config_value(dest, action, value))


def main(argv: Optional[List[str]] = None) -> None:
    p = build_parser()
    argv = list(sys.argv[1:] if argv is None else argv)
    a = p.parse_args(argv)
    _apply_config(p, a, argv)
    a.func(a)
