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


def sample_name(bam: str) -> str:
    base = re.sub(r"\.bam$", "", os.path.basename(bam), flags=re.I)
    return SUFFIX_RE.sub("", base) or base


def run_name(bam_dir: str) -> str:
    parts = [x for x in re.split(r"[\\/]+", os.path.normpath(bam_dir)) if x]
    for x in reversed(parts):
        if x.lower() not in RUN_SKIP:
            return x
    return parts[-1] if parts else "run"


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
    df = pd.read_csv(path, sep="\t")
    for c in ("sample", "bam"):
        if c not in df.columns:
            raise SystemExit(f"Missing column '{c}' in {path}")
    if "run" not in df.columns:
        df["run"] = "run"
    if "has_index" in df.columns:
        df = df[df["has_index"].fillna(1).astype(int) == 1]
    return df.reset_index(drop=True)


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


def cmd_analyse(a) -> None:
    import numpy as np
    import pandas as pd

    from . import report as R
    from . import stats as S

    ss = _read_samplesheet(a.samplesheet)
    counts_dir = os.path.join(a.outdir, "counts")
    figdir = os.path.join(a.outdir, "figures")
    os.makedirs(figdir, exist_ok=True)

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

    intron, exon, target_ids, samples = S.load_matrix(counts_dir, list(ss["sample"]))
    meta = S.load_meta(a.targets).reindex(target_ids)
    runs = ss.set_index("sample")["run"].reindex(samples).fillna("run")
    _p(f"      {len(samples)} samples x {len(target_ids)} targets")
    if len(samples) < a.min_cohort:
        _p(f"      [WARNING] only {len(samples)} samples. This tool scores each "
           f"sample against the cohort itself; below ~{a.min_cohort} samples the "
           f"median/MAD reference is unstable and the scores are unreliable.")

    kind = meta["kind"].to_numpy()
    is_minor = kind == "minor"
    is_control = kind == "control"

    _p("[2/5] normalising and scoring...")
    y = S.build_y(intron, exon, a.min_exon_depth)
    measurable = np.isfinite(y).mean(axis=0) >= a.min_frac_measurable
    y[:, ~measurable] = np.nan
    y_adj, offsets = S.normalise_by_controls(y, is_control & measurable)
    z = S.robust_z(y_adj, loo_max=a.loo_max)
    scores, hit = S.sample_scores(z, y_adj, is_minor & measurable,
                                  is_control & measurable, samples,
                                  z_thr=a.z_threshold, min_delta=a.min_delta)
    scores["run"] = runs.to_numpy()
    scores = scores.sort_values(["n_hits", "stouffer"], ascending=False
                                ).reset_index(drop=True)
    n_minor = int((is_minor & measurable).sum())
    n_ctrl = int((is_control & measurable).sum())
    _p(f"      {n_minor} minor / {n_ctrl} control introns measurable")

    _p("[3/5] PCA...")
    sc, load, expl, shape = S.pca(z, samples, is_minor & measurable,
                                  n_comp=min(6, len(samples)), clip=a.pca_clip,
                                  scale=a.pca_scale)
    _p(f"      {shape[1]} introns; PC1={100*expl[0]:.1f}% PC2={100*expl[1]:.1f}%")

    _p("[4/5] writing tables...")
    sr = scores.set_index("sample")
    acp = sc.join(sr[[c for c in ("run", "n_hits", "n_hits_control", "enrichment",
                                  "tail_q", "MSRI", "stouffer")]])
    acp.index.name = "sample"
    acp.to_csv(os.path.join(a.outdir, "pca_coordinates.tsv"), sep="\t")
    pd.DataFrame({"component": list(sc.columns),
                  "variance_explained": [float(v) for v in expl]}).to_csv(
        os.path.join(a.outdir, "pca_variance.tsv"), sep="\t", index=False)
    kept = np.where(is_minor & measurable)[0][:load.shape[1]]
    pd.DataFrame(load.T, columns=list(sc.columns),
                 index=[target_ids[i] for i in kept]).assign(
        gene_name=[meta.iloc[i]["gene_name"] for i in kept]).to_csv(
        os.path.join(a.outdir, "pca_loadings.tsv"), sep="\t")
    scores.to_csv(os.path.join(a.outdir, "sample_scores.tsv"), sep="\t", index=False)

    events = S.hit_details(counts_dir, samples, target_ids, hit, z, y_adj, meta,
                           is_minor, max_per_sample=a.max_events_per_sample)
    if len(events):
        events.to_csv(os.path.join(a.outdir, "hits.tsv.gz"), sep="\t",
                      index=False, compression="gzip")

    _p("[5/5] figures and report...")
    with np.errstate(all="ignore"):
        qc = pd.DataFrame(dict(
            sample=samples, run=runs.to_numpy(),
            median_exon_depth=np.nanmedian(exon, axis=1),
            frac_measurable=(exon >= a.min_exon_depth).mean(axis=1),
            control_background=offsets))
    f_scores = R.fig_score_distribution(scores, figdir, col=a.score_column)
    f_hm = R.fig_heatmap(z, samples, target_ids, is_minor & measurable, meta,
                         scores, figdir, n_samples=a.heatmap_samples,
                         n_introns=a.heatmap_introns)
    f_qc = R.fig_qc(qc, figdir)

    def col(name):
        return sr[name].reindex(samples).to_numpy()

    def clean(v):
        return None if (v is None or not np.isfinite(v)) else round(float(v), 4)

    points = [dict(n=s, sh=(s if len(s) <= 22 else s[:21] + "\u2026"),
                   x=round(float(sc.iloc[i, 0]), 3),
                   y=round(float(sc.iloc[i, 1]), 3),
                   h=int(col("n_hits")[i]), hc=int(col("n_hits_control")[i]),
                   e=clean(col("enrichment")[i]), q=clean(col("tail_q")[i]),
                   m=clean(col("MSRI")[i]), st=clean(col("stouffer")[i]),
                   r=str(runs.iloc[i]))
              for i, s in enumerate(samples)]

    show = ["sample", "run", "n_hits", "n_hits_control", "enrichment", "tail_q",
            "MSRI", "stouffer", "n_minor_measurable"]
    run_tbl = pd.DataFrame()
    if runs.nunique() > 1:
        run_tbl = (scores.groupby("run")
                   .agg(n=("sample", "size"),
                        n_hits_median=("n_hits", "median"),
                        n_hits_max=("n_hits", "max"),
                        control_hits_median=("n_hits_control", "median"))
                   .sort_values("n_hits_median", ascending=False).reset_index())

    ev_cols = [c for c in ["sample", "gene_name", "chrom", "intron_start",
                           "intron_end", "strand", "z", "y_adj",
                           "intron_depth_mean", "exon_depth_mean", "n_EE",
                           "n_EI_donor", "n_EI_acceptor", "n_alt_5p", "n_alt_3p"]
               if not len(events) or c in events.columns]
    top_ev = (events.head(0) if not len(events) else
              events.groupby("sample").head(3).head(60)[ev_cols])

    R.write_html(os.path.join(a.outdir, "report.html"), dict(
        n_samples=len(samples), n_minor=n_minor, n_control=n_ctrl,
        points=points, xlab=f"PC1 ({100*expl[0]:.1f} %)",
        ylab=f"PC2 ({100*expl[1]:.1f} %)",
        top_samples=scores.head(a.top_samples)[show],
        fig_scores=f_scores, fig_heatmap=f_hm, fig_qc=f_qc,
        top_events=top_ev, run_table=run_tbl,
        params=pd.DataFrame([(k, str(v)) for k, v in vars(a).items()
                             if k != "func"], columns=["parameter", "value"])))

    _p("")
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        _p(scores.head(15)[show].to_string(index=False))
    _p("")
    _p(f"Report : {os.path.join(a.outdir, 'report.html')}")
    _p(f"Tables : {a.outdir}/sample_scores.tsv, pca_coordinates.tsv, hits.tsv.gz")


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
    p.add_argument("--config", help="JSON file of parameters")
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp):
        sp.add_argument("--config", help="JSON file of parameters")
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
        sp.add_argument("--loo-max", type=int, default=60,
                        help="exact leave-one-out below this cohort size")
        sp.add_argument("--min-cohort", type=int, default=20)
        sp.add_argument("--pca-clip", type=float, default=8.0)
        sp.add_argument("--pca-scale", action="store_true")
        sp.add_argument("--score-column", default="n_hits")
        sp.add_argument("--top-samples", type=int, default=25)
        sp.add_argument("--heatmap-samples", type=int, default=30)
        sp.add_argument("--heatmap-introns", type=int, default=45)
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


def main(argv: Optional[List[str]] = None) -> None:
    p = build_parser()
    argv = list(sys.argv[1:] if argv is None else argv)
    a = p.parse_args(argv)
    if getattr(a, "config", None):
        with open(a.config, encoding="utf-8") as fh:
            cfg = json.load(fh)
        known, given = _all_dests(p), _dests_given(p, argv)
        for k, v in cfg.items():
            if k.startswith("_"):
                continue
            k = k.replace("-", "_")
            if k in given:          # the command line wins over the config file
                continue
            if hasattr(a, k):
                setattr(a, k, v)
            elif k not in known:
                _p(f"[warn] unknown config parameter ignored: {k}")
    a.func(a)
