"""Normalisation, robust z-scores, per-sample aggregation and PCA.

There is no case/control design here. The cohort IS the background: every
sample is scored against the robust median and MAD of all the others, intron by
intron. That works because a minor-spliceosome defect is rare and affects a set
of introns that differs from patient to patient, so an affected sample is an
outlier against its own cohort.

Measurement chain
-----------------
1. per intron and per sample
       y = log2( (intron depth + a) / (flanking-exon depth + a) )
   i.e. intronic coverage normalised by the exons immediately around it.

2. within-sample normalisation by the MAJOR (U2) control introns
       y_adj = y - median( y over measurable control introns )
   This absorbs whatever makes intronic background vary between samples: RNA
   degradation, genomic-DNA carry-over, nuclear pre-mRNA fraction, depth.
   Without it a slightly degraded sample looks like a global positive.

3. robust z per intron against the cohort
       z = 0.6745 * (y_adj - median) / MAD

4. per-sample aggregation -- the part that gives sensitivity. Patients do not
   share the same affected introns, so a per-intron test has no power; we
   aggregate over the whole minor-intron set:
     - n_hits      : minor introns passing the z and effect-size thresholds
     - tail_p      : Fisher test, n_hits on minor introns vs n_hits on the
                     control introns OF THE SAME SAMPLE. Self-calibrated, and
                     the right test for a sparse signal (a location test such
                     as Mann-Whitney has almost no power when only a few
                     percent of values shift).
     - enrichment  : ratio of those two hit rates
     - stouffer    : sum(z)/sqrt(n), sensitive to a small but coordinated shift
     - MSRI        : winsorised mean of z(minor) minus that of z(control)

Memory: matrices are held as float32 numpy arrays, never as a long-format
DataFrame, so a cohort of several thousand BAMs stays around a few hundred MB.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

PSEUDO = 0.5
NUM_COLS = ["target_id", "intron_depth_mean", "exon_depth_mean"]
DETAIL_COLS = ["target_id", "intron_depth_mean", "exon_depth_mean",
               "exon_depth_left", "exon_depth_right", "n_intron_reads",
               "n_EE", "n_EI_left", "n_EI_right", "n_alt_donor",
               "n_alt_acceptor", "intron_frac_cov"]


def counts_path(counts_dir: str, sample: str) -> str:
    return os.path.join(counts_dir, f"{sample}.tsv.gz")


def load_meta(targets_tsv: str) -> pd.DataFrame:
    m = pd.read_csv(targets_tsv, sep="\t", dtype={"chrom": str})
    m = m[m["usable"] == 1] if "usable" in m.columns else m
    return m.set_index("target_id")


def load_matrix(counts_dir: str, samples: Sequence[str],
                progress: int = 250) -> Tuple[np.ndarray, np.ndarray, List[str],
                                              List[str]]:
    """Read every sample's counts into two float32 matrices.

    Returns (intron_depth, exon_depth, target_ids, samples_found).
    """
    ref_ids: Optional[List[str]] = None
    rows_i: List[np.ndarray] = []
    rows_e: List[np.ndarray] = []
    found: List[str] = []
    for k, s in enumerate(samples):
        p = counts_path(counts_dir, s)
        if not os.path.exists(p):
            continue
        d = pd.read_csv(p, sep="\t", usecols=NUM_COLS)
        if ref_ids is None:
            ref_ids = list(d["target_id"])
            index = {t: i for i, t in enumerate(ref_ids)}
        if list(d["target_id"]) != ref_ids:
            # different target order/set: realign rather than fail
            d = d.set_index("target_id").reindex(ref_ids).reset_index()
        rows_i.append(d["intron_depth_mean"].to_numpy(np.float32))
        rows_e.append(d["exon_depth_mean"].to_numpy(np.float32))
        found.append(s)
        if progress and k and k % progress == 0:
            print(f"      {k}/{len(samples)} samples read", flush=True)
    if not found:
        raise SystemExit("No count file found. Run `quantify` first.")
    return np.vstack(rows_i), np.vstack(rows_e), ref_ids, found


def build_y(intron: np.ndarray, exon: np.ndarray, min_exon_depth: float
            ) -> np.ndarray:
    y = np.log2((intron + PSEUDO) / (exon + PSEUDO))
    y[exon < min_exon_depth] = np.nan
    return y


def normalise_by_controls(y: np.ndarray, is_control: np.ndarray) -> Tuple[
        np.ndarray, np.ndarray]:
    """Subtract each sample's own control-intron median."""
    if is_control.sum() < 20:
        return y, np.zeros(y.shape[0], dtype=np.float32)
    with np.errstate(all="ignore"):
        off = np.nanmedian(y[:, is_control], axis=1).astype(np.float32)
    off = np.nan_to_num(off)
    return y - off[:, None], off


def robust_z(y: np.ndarray, loo_max: int = 60, mad_floor_q: float = 0.2
             ) -> np.ndarray:
    """Robust z per column against the rest of the cohort.

    For a small cohort we do an exact leave-one-out, so a sample cannot mask
    itself. Above `loo_max` samples we use the plain median/MAD: removing one
    sample out of hundreds shifts a median by far less than the MAD floor, and
    exact LOO would cost one full pass per sample.
    """
    n = y.shape[0]
    with np.errstate(all="ignore"):
        med = np.nanmedian(y, axis=0)
        mad = np.nanmedian(np.abs(y - med), axis=0)
    floor = float(np.nanquantile(mad, mad_floor_q))
    floor = max(floor if np.isfinite(floor) else 0.05, 0.05)
    mad = np.maximum(np.nan_to_num(mad, nan=floor), floor)

    if n > loo_max:
        return (0.6745 * (y - med) / mad).astype(np.float32)

    z = np.empty_like(y)
    for i in range(n):
        keep = np.arange(n) != i
        with np.errstate(all="ignore"):
            m = np.nanmedian(y[keep], axis=0)
            d = np.nanmedian(np.abs(y[keep] - m), axis=0)
        d = np.maximum(np.nan_to_num(d, nan=floor), floor)
        z[i] = 0.6745 * (y[i] - m) / d
    return z.astype(np.float32)


def winsor_mean(x: np.ndarray, lo: float = -4.0, hi: float = 8.0) -> float:
    """Mean after clipping -- NOT a symmetric trimmed mean.

    A 10% symmetric trim would discard the top 10% of values, i.e. exactly the
    signal when a patient retains only a few percent of minor introns. Measured
    on simulation: AUC 0.85 with trimming vs 1.00 with clipping.
    """
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan")
    return float(np.clip(x, lo, hi).mean())


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    out = np.full(p.shape, np.nan)
    q = p[ok]
    n = q.size
    if n == 0:
        return out
    o = np.argsort(q)
    r = np.empty(n)
    r[o] = np.minimum.accumulate((q[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out[ok] = np.clip(r, 0, 1)
    return out


def sample_scores(z: np.ndarray, y_adj: np.ndarray, is_minor: np.ndarray,
                  is_control: np.ndarray, samples: Sequence[str],
                  z_thr: float = 3.0, min_delta: float = 0.5) -> Tuple[
                      pd.DataFrame, np.ndarray]:
    """One row per sample. Also returns the boolean hit matrix."""
    from scipy import stats

    with np.errstate(all="ignore"):
        med_ref = np.nanmedian(y_adj, axis=0)
    delta = y_adj - med_ref
    hit = np.isfinite(z) & (z >= z_thr) & np.isfinite(delta) & (delta >= min_delta)

    zm, zc = z[:, is_minor], z[:, is_control]
    hm = hit[:, is_minor].sum(axis=1)
    hc = hit[:, is_control].sum(axis=1)
    nm = np.isfinite(zm).sum(axis=1)
    nc = np.isfinite(zc).sum(axis=1)

    tail_p = np.full(len(samples), np.nan)
    enrich = np.full(len(samples), np.nan)
    for i in range(len(samples)):
        if nm[i] >= 20 and nc[i] >= 20:
            _odds, tail_p[i] = stats.fisher_exact(
                [[int(hm[i]), int(nm[i] - hm[i])],
                 [int(hc[i]), int(nc[i] - hc[i])]], alternative="greater")
            enrich[i] = ((hm[i] + 0.5) / (nm[i] + 1)) / ((hc[i] + 0.5) / (nc[i] + 1))

    with np.errstate(all="ignore"):
        df = pd.DataFrame(dict(
            sample=list(samples),
            n_hits=hm.astype(int), n_hits_control=hc.astype(int),
            enrichment=enrich, tail_p=tail_p,
            n_minor_measurable=nm.astype(int),
            n_control_measurable=nc.astype(int),
            MSRI=[winsor_mean(zm[i]) - winsor_mean(zc[i]) for i in range(len(samples))],
            stouffer=np.nansum(zm, axis=1) / np.sqrt(np.maximum(nm, 1)),
            median_z_minor=np.nanmedian(zm, axis=1),
            median_z_control=np.nanmedian(zc, axis=1),
        ))
    df["tail_q"] = bh(df["tail_p"].to_numpy())
    return df, hit


def pca(z: np.ndarray, samples: Sequence[str], keep: np.ndarray,
        n_comp: int = 6, clip: float = 8.0, scale: bool = False):
    """PCA by SVD (numpy only, no scikit-learn).

    Run this on the z-score matrix, NOT on y_adj: z already puts every intron on
    a comparable scale, so re-standardising column by column would give the
    noise of stable introns the same weight as the signal and the PCA would stop
    separating anything (measured: AUC 0.72 with scaling vs 1.00 without).

    Values are winsorised: otherwise one wildly outlying intron captures PC1 on
    its own.
    """
    X = z[:, keep].astype(np.float64)
    col_ok = np.isfinite(X).sum(axis=0) >= max(3, int(0.6 * X.shape[0]))
    X = X[:, col_ok]
    with np.errstate(all="ignore"):
        colmed = np.nanmedian(X, axis=0)
    bad = ~np.isfinite(X)
    X[bad] = np.take(np.nan_to_num(colmed), np.where(bad)[1])
    X = X[:, X.std(axis=0) > 1e-9]
    if X.shape[1] < 3 or X.shape[0] < 3:
        raise SystemExit("Not enough usable data for the PCA.")
    if clip:
        X = np.clip(X, -clip / 2.0, clip)
    X = X - X.mean(axis=0)
    if scale:
        sd = X.std(axis=0, ddof=1)
        sd[sd == 0] = 1.0
        X = X / sd
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    n_comp = min(n_comp, S.size)
    scores = U[:, :n_comp] * S[:n_comp]
    var = (S ** 2) / max(X.shape[0] - 1, 1)
    expl = var / var.sum()

    # SVD component signs are arbitrary: orient so that "further right" always
    # means "more retention", otherwise the figure flips between runs.
    agg = X.mean(axis=1)
    for j in range(n_comp):
        if np.corrcoef(scores[:, j], agg)[0, 1] < 0:
            scores[:, j] = -scores[:, j]
            Vt[j] = -Vt[j]

    sc = pd.DataFrame(scores, index=list(samples),
                      columns=[f"PC{i+1}" for i in range(n_comp)])
    return sc, Vt[:n_comp], expl[:n_comp], X.shape


def hit_details(counts_dir: str, samples: Sequence[str], target_ids: Sequence[str],
                hit: np.ndarray, z: np.ndarray, y_adj: np.ndarray,
                meta: pd.DataFrame, is_minor: np.ndarray,
                max_per_sample: int = 200) -> pd.DataFrame:
    """Second pass: detailed rows for hits only.

    Writing every intron x every sample would be tens of millions of rows for a
    large cohort, so only the calls are detailed.
    """
    tid = np.asarray(target_ids)
    strand = meta.reindex(tid)["strand"].to_numpy()
    out = []
    for i, s in enumerate(samples):
        idx = np.where(hit[i] & is_minor)[0]
        if idx.size == 0:
            continue
        idx = idx[np.argsort(-z[i, idx])][:max_per_sample]
        p = counts_path(counts_dir, s)
        cols = [c for c in DETAIL_COLS
                if c in pd.read_csv(p, sep="\t", nrows=0).columns]
        d = pd.read_csv(p, sep="\t", usecols=cols).set_index("target_id")
        sub = d.reindex(tid[idx]).reset_index()
        sub.insert(0, "sample", s)
        sub["z"] = z[i, idx]
        sub["y_adj"] = y_adj[i, idx]
        minus = strand[idx] == "-"
        # Raw columns are labelled by GENOMIC left/right. For a minus-strand
        # gene the genomic-left boundary is the acceptor, not the donor, so we
        # relabel here (post hoc, to keep the counts cache valid).
        sub["n_EI_donor"] = np.where(minus, sub["n_EI_right"], sub["n_EI_left"])
        sub["n_EI_acceptor"] = np.where(minus, sub["n_EI_left"], sub["n_EI_right"])
        sub["n_alt_5p"] = np.where(minus, sub["n_alt_acceptor"], sub["n_alt_donor"])
        sub["n_alt_3p"] = np.where(minus, sub["n_alt_donor"], sub["n_alt_acceptor"])
        sub = sub.drop(columns=["n_EI_left", "n_EI_right", "n_alt_donor",
                                "n_alt_acceptor"], errors="ignore")
        out.append(sub)
    if not out:
        return pd.DataFrame()
    res = pd.concat(out, ignore_index=True)
    keep = ["gene_name", "gene_id", "chrom", "intron_start", "intron_end",
            "strand", "intron_len"]
    res = res.merge(meta.reset_index()[["target_id"] + keep], on="target_id",
                    how="left")
    return res.sort_values(["sample", "z"], ascending=[True, False])
