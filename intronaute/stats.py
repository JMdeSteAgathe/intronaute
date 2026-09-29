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

Junction measure (optional, `--hit-mode`)
-----------------------------------------
Coverage alone has two weaknesses: in exon-capture libraries the probes barely
reach the intron body, and in noisy samples intronic coverage can rise with no
retention at all (reads piling up mid-intron without touching either splice
site). The junction measure only looks at the two splice sites:

    EI = reads crossing the donor + reads crossing the acceptor (unspliced,
         >= min_overhang bp on each side of the boundary)
    S  = spliced reads using the donor + spliced reads using the acceptor
         (the exact junction counts at both sites, plus alternative partners)
    theta = EI / (EI + S)                 -- FRASER's theta, pooled over both sites
    y_J   = log2((EI + a) / (S + a))      -- same scale as y, measurable when
                                             EI + S >= min_junction_reads

y_J then goes through exactly the same chain as y: control-intron offset,
robust z, effect-size guard. The two z are combined per intron as
    z_comb = (z_cov + z_J) / sqrt(2 + 2 rho)
where rho is their correlation measured on the CONTROL introns of the cohort,
so z_comb stays on the same scale as a single z under the null. When only one
measure is available for a cell, z_comb falls back to that one.

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


JUNC_COLS = ["n_EE", "n_EI_left", "n_EI_right", "n_alt_donor", "n_alt_acceptor"]
HIT_MODES = ("coverage", "junction", "combined", "both")
# MAD floor quantile for the junction z. The coverage floor (0.2) is fine for
# coverage, whose MAD mostly reflects biology. Junction MADs are dominated by
# counting noise, which shrinks with depth: a floor taken at the 20th
# percentile comes from shallow introns and caps the z of well-covered ones.
# Measured on simulation: sensitivity 0.04 at 0.2, 0.19 at 0.05, no gain below.
JUNCTION_MAD_FLOOR_Q = 0.05


def counts_path(counts_dir: str, sample: str) -> str:
    return os.path.join(counts_dir, f"{sample}.tsv.gz")


def load_meta(targets_tsv: str) -> pd.DataFrame:
    m = pd.read_csv(targets_tsv, sep="\t", dtype={"chrom": str})
    m = m[m["usable"] == 1] if "usable" in m.columns else m
    return m.set_index("target_id")


class Counts:
    """Per-sample x per-target float32 matrices read from the counts cache."""

    def __init__(self, intron, exon, ei, split, target_ids, samples):
        self.intron, self.exon = intron, exon
        self.ei, self.split = ei, split          # None if the cache lacks them
        self.target_ids, self.samples = target_ids, samples


def load_counts(counts_dir: str, samples: Sequence[str], junctions: bool = True,
                progress: int = 250) -> Counts:
    """Read every sample's counts.

    With `junctions`, also build two junction matrices:
      ei    = n_EI_left + n_EI_right                       (unspliced, crossing)
      split = 2 * n_EE + n_alt_donor + n_alt_acceptor      (spliced, per site)
    i.e. per splice site, the reads that do not splice vs those that do, summed
    over the two sites of the intron.
    """
    ref_ids: Optional[List[str]] = None
    rows = {k: [] for k in ("i", "e", "ei", "sp")}
    found: List[str] = []
    want = list(NUM_COLS) + (JUNC_COLS if junctions else [])
    have_junc = junctions
    for k, s in enumerate(samples):
        p = counts_path(counts_dir, s)
        if not os.path.exists(p):
            continue
        header = pd.read_csv(p, sep="\t", nrows=0).columns
        if have_junc and not set(JUNC_COLS) <= set(header):
            print(f"      [warn] {s}: junction columns missing from the counts "
                  f"cache; junction measure disabled", flush=True)
            have_junc = False
        d = pd.read_csv(p, sep="\t", usecols=[c for c in want if c in header])
        if ref_ids is None:
            ref_ids = list(d["target_id"])
        if list(d["target_id"]) != ref_ids:
            # different target order/set: realign rather than fail
            d = d.set_index("target_id").reindex(ref_ids).reset_index()
        rows["i"].append(d["intron_depth_mean"].to_numpy(np.float32))
        rows["e"].append(d["exon_depth_mean"].to_numpy(np.float32))
        if have_junc:
            rows["ei"].append((d["n_EI_left"] + d["n_EI_right"]).to_numpy(np.float32))
            rows["sp"].append((2 * d["n_EE"] + d["n_alt_donor"]
                               + d["n_alt_acceptor"]).to_numpy(np.float32))
        found.append(s)
        if progress and k and k % progress == 0:
            print(f"      {k}/{len(samples)} samples read", flush=True)
    if not found:
        raise SystemExit("No count file found. Run `quantify` first.")
    ok = have_junc and len(rows["ei"]) == len(found)
    return Counts(np.vstack(rows["i"]), np.vstack(rows["e"]),
                  np.vstack(rows["ei"]) if ok else None,
                  np.vstack(rows["sp"]) if ok else None, ref_ids, found)


def load_matrix(counts_dir: str, samples: Sequence[str],
                progress: int = 250) -> Tuple[np.ndarray, np.ndarray, List[str],
                                              List[str]]:
    """Backward-compatible wrapper: (intron_depth, exon_depth, target_ids, samples)."""
    c = load_counts(counts_dir, samples, junctions=False, progress=progress)
    return c.intron, c.exon, c.target_ids, c.samples


def build_y(intron: np.ndarray, exon: np.ndarray, min_exon_depth: float
            ) -> np.ndarray:
    y = np.log2((intron + PSEUDO) / (exon + PSEUDO))
    y[exon < min_exon_depth] = np.nan
    return y


def build_y_junction(ei: np.ndarray, split: np.ndarray, min_reads: float
                     ) -> Tuple[np.ndarray, np.ndarray]:
    """-> (y_J, theta). NaN where fewer than `min_reads` reads touch the two
    splice sites, since a ratio from 3 reads is noise."""
    with np.errstate(all="ignore"):
        yj = np.log2((ei + PSEUDO) / (split + PSEUDO)).astype(np.float32)
        theta = (ei / (ei + split)).astype(np.float32)
    low = (ei + split) < min_reads
    yj[low] = np.nan
    theta[low] = np.nan
    return yj, theta


def delta_from_cohort(y_adj: np.ndarray) -> np.ndarray:
    """Effect size: distance to the cohort median of the same intron."""
    with np.errstate(all="ignore"):
        return y_adj - np.nanmedian(y_adj, axis=0)


def combine_z(z_cov: np.ndarray, z_j: np.ndarray, is_control: np.ndarray,
              clip: float = 6.0, min_pairs: int = 200) -> Tuple[np.ndarray, float]:
    """Stouffer-like combination, corrected for the correlation of the two z.

    Both measures see the same retained transcripts, so they are correlated;
    dividing by sqrt(2) would inflate the combined score under the null.
    rho is estimated on control introns (the null set) with values clipped so
    a few true events cannot drive it.
    """
    both = np.isfinite(z_cov) & np.isfinite(z_j)
    m = both & is_control[None, :]
    rho = 0.0
    if m.sum() >= min_pairs:
        a = np.clip(z_cov[m], -clip, clip)
        b = np.clip(z_j[m], -clip, clip)
        if a.std() > 0 and b.std() > 0:
            rho = float(np.clip(np.corrcoef(a, b)[0, 1], 0.0, 0.95))
    out = np.where(np.isfinite(z_cov), z_cov, z_j).astype(np.float32)
    out[both] = ((z_cov[both] + z_j[both]) / np.sqrt(2.0 + 2.0 * rho)
                 ).astype(np.float32)
    return out, rho


def call_hits(mode: str, z_thr: float, z_cov: np.ndarray, d_cov: np.ndarray,
              min_delta: float, z_j: Optional[np.ndarray] = None,
              d_j: Optional[np.ndarray] = None, min_delta_j: float = 0.5,
              z_comb: Optional[np.ndarray] = None,
              junction_veto: float = 1.0) -> np.ndarray:
    """Boolean hit matrix. NaN comparisons are False, so an unmeasurable cell
    is never a hit.

    coverage : z_cov >= thr and d_cov >= min_delta             (original rule)
    junction : z_J   >= thr and d_J   >= min_delta_j
    combined : z_comb >= thr, at least one measure passing its effect-size
               guard, and NOT vetoed: if the junctions are measurable and look
               normal (z_J < junction_veto) the intron is not called whatever
               its coverage -- intronic coverage that never touches a splice
               site is background, not retention.
    both     : coverage rule AND junction rule, each at z_thr (strictest; meant
               to be used with a lower --z-threshold)
    """
    def rule(z, d, md):
        with np.errstate(invalid="ignore"):
            return np.isfinite(z) & (z >= z_thr) & np.isfinite(d) & (d >= md)

    if mode == "coverage":
        return rule(z_cov, d_cov, min_delta)
    if z_j is None:
        raise ValueError(f"hit mode '{mode}' needs junction counts")
    if mode == "junction":
        return rule(z_j, d_j, min_delta_j)
    if mode == "both":
        return rule(z_cov, d_cov, min_delta) & rule(z_j, d_j, min_delta_j)
    if mode == "combined":
        with np.errstate(invalid="ignore"):
            effect = ((np.isfinite(d_cov) & (d_cov >= min_delta))
                      | (np.isfinite(d_j) & (d_j >= min_delta_j)))
            vetoed = np.isfinite(z_j) & (z_j < junction_veto)
            return np.isfinite(z_comb) & (z_comb >= z_thr) & effect & ~vetoed
    raise ValueError(f"unknown hit mode: {mode} (choose from {HIT_MODES})")


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
                  z_thr: float = 3.0, min_delta: float = 0.5,
                  hit: Optional[np.ndarray] = None) -> Tuple[
                      pd.DataFrame, np.ndarray]:
    """One row per sample. Also returns the boolean hit matrix.

    `z` is the score used for stouffer / MSRI. If `hit` is None the original
    coverage rule is applied to (z, y_adj); otherwise `hit` (from call_hits)
    is used as is, so the tail test counts whatever the chosen mode calls, on
    minor and control introns alike.
    """
    from scipy import stats

    if hit is None:
        hit = call_hits("coverage", z_thr, z, delta_from_cohort(y_adj), min_delta)

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
                max_per_sample: int = 200,
                extra: Optional[Dict[str, np.ndarray]] = None) -> pd.DataFrame:
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
        for name, mat in (extra or {}).items():
            sub[name] = mat[i, idx]
        minus = strand[idx] == "-"
        # Raw columns are labelled by GENOMIC left/right. For a minus-strand
        # gene the genomic-left boundary is the acceptor, not the donor, so we
        # relabel here (post hoc, to keep the counts cache valid).
        for lo, hi, five, three in (("n_EI_left", "n_EI_right", "n_EI_donor",
                                     "n_EI_acceptor"),
                                    ("n_alt_donor", "n_alt_acceptor", "n_alt_5p",
                                     "n_alt_3p")):
            if lo in sub.columns and hi in sub.columns:
                sub[five] = np.where(minus, sub[hi], sub[lo])
                sub[three] = np.where(minus, sub[lo], sub[hi])
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
