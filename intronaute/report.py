"""Figures and a self-contained HTML report.

The report opens with a contents bar, then an interactive PCA: a canvas scatter with a
search box, zoom/pan, hover tooltips and click-to-pin labels. It is plain
JavaScript with no CDN, so the single .html file works offline and can be
emailed as-is. Canvas (not SVG) because a cohort can hold several thousand
points.
"""
from __future__ import annotations

import base64
import datetime
import html
import json
import os
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def _save(fig, path: str, dpi: int = 140) -> str:
    fig.tight_layout()
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return path


def _b64(path: str) -> str:
    with open(path, "rb") as fh:
        return base64.b64encode(fh.read()).decode()


# --------------------------------------------------------------------------- #
# static figures
# --------------------------------------------------------------------------- #
def fig_score_distribution(scores: pd.DataFrame, outdir: str,
                           col: str = "n_hits", n_label: int = 12) -> str:
    d = scores.sort_values(col, ascending=False).reset_index(drop=True)
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4))

    v = d[col].to_numpy(float)
    axes[0].plot(np.arange(1, len(v) + 1), v, lw=1.1, color="#4a5568")
    axes[0].scatter(np.arange(1, min(n_label, len(v)) + 1), v[:n_label],
                    s=28, color="#c1272d", zorder=3)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("rank (log scale)")
    axes[0].set_ylabel(col)
    axes[0].set_title(f"{col} across the cohort (n={len(d)})", fontsize=10)
    axes[0].spines[["top", "right"]].set_visible(False)

    m = d["n_hits"].to_numpy(float)
    c = d["n_hits_control"].to_numpy(float)
    axes[1].scatter(c + 0.5, m + 0.5, s=16, alpha=.5, color="#2b6cb0",
                    edgecolor="none")
    lim = max(m.max(), c.max()) + 2
    axes[1].plot([0.5, lim], [0.5, lim], lw=.9, color="#a0aec0", ls="--")
    for _, r in d.head(n_label).iterrows():
        axes[1].scatter(r["n_hits_control"] + .5, r["n_hits"] + .5, s=34,
                        facecolor="none", edgecolor="#c1272d", lw=1.4, zorder=3)
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("hits on CONTROL introns (per sample)")
    axes[1].set_ylabel("hits on MINOR introns")
    axes[1].set_title("Minor vs control hits\n(dashed line = no enrichment)",
                      fontsize=10)
    axes[1].spines[["top", "right"]].set_visible(False)
    return _save(fig, os.path.join(outdir, "score_distribution.png"))


def fig_heatmap(z: np.ndarray, samples: Sequence[str], target_ids: Sequence[str],
                is_minor: np.ndarray, meta: pd.DataFrame, scores: pd.DataFrame,
                outdir: str, n_samples: int = 30, n_introns: int = 45,
                vmax: float = 6.0, rank_col: str = "n_hits",
                ascending: bool = False, filename: str = "heatmap.png",
                rank_label: str = "n_hits",
                colorbar_label: str = "robust z above this sample's control introns"
                ) -> str:
    """Top samples by `rank_col`, columns = minor introns with the highest z in those rows."""
    order = scores.sort_values(rank_col, ascending=ascending, na_position="last")[
        "sample"].head(n_samples)
    pos = {s: i for i, s in enumerate(samples)}
    rows = [pos[s] for s in order if s in pos]
    if not rows:
        return ""
    cols = np.where(is_minor)[0]
    sub = z[np.ix_(rows, cols)]
    rank = np.nanmax(sub, axis=0)
    show_all = n_introns is None or n_introns <= 0
    n_keep = cols.size if show_all else n_introns
    top = cols[np.argsort(-np.nan_to_num(rank, nan=-np.inf))[:n_keep]]
    M = z[np.ix_(rows, top)]
    labels = [str(meta.loc[target_ids[j], "gene_name"]) for j in top]
    wide = len(top) > 80
    fig, ax = plt.subplots(figsize=(max(7.0, (0.18 if wide else 0.24) * len(top)),
                                    max(3.4, 0.28 * len(rows) + (1.4 if wide else 0))))
    im = ax.imshow(M, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(np.arange(len(top)))
    ax.set_xticklabels(labels, rotation=90, fontsize=6.5)
    ax.set_yticks(np.arange(len(rows)))
    ax.set_yticklabels([samples[i] for i in rows], fontsize=7)
    intron_txt = (f"all {len(top)} measurable minor introns" if show_all
                  else f"top {len(top)} minor introns")
    ax.set_title(f"Top {len(rows)} samples by {rank_label} x {intron_txt} "
                 f"(z minus this sample's control z)", fontsize=10)
    # The heatmap is often wider than the page. A scale on the right sits at
    # the last columns. A short bar under the left edge stays with the samples.
    fig_w = float(fig.get_size_inches()[0])
    shrink = min(0.55, 4.0 / max(fig_w, 4.0))
    cbar = fig.colorbar(
        im, ax=ax, orientation="horizontal", location="bottom",
        shrink=shrink, pad=0.22, anchor=(0.0, 1.0), aspect=28)
    cbar.set_label(colorbar_label, fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    return _save(fig, os.path.join(outdir, filename), dpi=90 if wide else 140)


def fig_qc(qc: pd.DataFrame, outdir: str) -> str:
    """Depth, measurable fraction and intronic background, grouped by run."""
    runs = sorted(qc["run"].astype(str).unique())
    many = len(qc) > 40
    fig, axes = plt.subplots(1, 3, figsize=(13.5, max(3.4, .30 * len(runs) + 1.4)))
    rng = np.random.default_rng(0)
    cols = ["median_exon_depth", "frac_measurable", "control_background"]
    labs = ["median exonic depth", "fraction of targets measurable",
            "intronic background (log2, control introns)"]
    for ax, c, lab in zip(axes, cols, labs):
        if many:
            for yi, r in enumerate(runs):
                v = pd.to_numeric(qc.loc[qc["run"].astype(str) == r, c],
                                  errors="coerce").to_numpy()
                ax.scatter(v, yi + rng.uniform(-.16, .16, v.size), s=13,
                           alpha=.7, color="#2b6cb0", edgecolor="none")
                if np.isfinite(v).any():
                    ax.plot([np.nanmedian(v)] * 2, [yi - .3, yi + .3],
                            color="#111", lw=1.4)
            ax.set_yticks(range(len(runs)))
            ax.set_yticklabels(runs, fontsize=7)
        else:
            d = qc.sort_values(c)
            ax.barh(np.arange(len(d)), pd.to_numeric(d[c], errors="coerce"),
                    color="#2b6cb0")
            ax.set_yticks(np.arange(len(d)))
            ax.set_yticklabels(d["sample"], fontsize=6.5)
        ax.set_xlabel(lab, fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_title("Quality control" + (" by run (bar = run median)" if many
                                           else ""), fontsize=10)
    return _save(fig, os.path.join(outdir, "qc.png"))


# Cells that are not a real label. Matching ignores case and surrounding spaces.
_BLANK_LABELS = frozenset({
    "", "0", "control", "n/a", "na", "null", "none", "nan",
})
# minor_spliceopathy values that say this person is not a spliceopathy case.
_SPLICE_NEGATIVE = frozenset({"false", "no", "non-ms", "nonms"})


def _cell_text(value) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def _is_blank_label(value) -> bool:
    """True for a missing label: blank, 0, control, N/A, NA, null, or none."""
    text = _cell_text(value).casefold()
    if text in _BLANK_LABELS:
        return True
    try:
        return float(text) == 0.0
    except ValueError:
        return False


def _choose_spelling(spellings: Sequence[str]) -> str:
    """One spelling for labels that differ only by case.

    The most common spelling wins. A tie keeps the spelling with more
    capitals, so MS is kept over ms when both appear equally often.
    """
    counts: Dict[str, int] = {}
    for text in spellings:
        counts[text] = counts.get(text, 0) + 1
    best = spellings[0]
    best_key = None
    for i, text in enumerate(spellings):
        key = (counts[text], sum(ch.isupper() for ch in text), -i)
        if best_key is None or key > best_key:
            best_key = key
            best = text
    return best


def _label_spellings(values) -> Dict[str, str]:
    """Map a case-folded label to the spelling shown in the report.

    Placeholders, including NULL, null, and a missing cell, are left out, so
    they share no label. Yes and yes share one spelling, and so do MS and ms.
    """
    buckets: Dict[str, List[str]] = {}
    for value in values:
        if _is_blank_label(value):
            continue
        text = _cell_text(value)
        buckets.setdefault(text.casefold(), []).append(text)
    return {key: _choose_spelling(seen) for key, seen in buckets.items()}


def _spliceopathy_mask(series: pd.Series) -> pd.Series:
    """True when minor_spliceopathy names a spliceopathy case.

    A person is left out when the cell is blank, 0, control, N/A, NA, null,
    or none, and when it is false, no, or non-MS. Yes and yes are the same,
    and so are MS and ms. NULL, null, and a missing cell are the same.
    The comparison ignores case and surrounding spaces. A boolean True counts.
    Any other text counts, so a label other than MS is still included.
    """
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    text = series.map(lambda v: _cell_text(v).casefold())
    return ~text.isin(_BLANK_LABELS | _SPLICE_NEGATIVE) & text.ne("")


def minor_hit_genes(hit: np.ndarray, samples: Sequence[str],
                    target_ids: Sequence[str], meta: pd.DataFrame,
                    is_minor: np.ndarray) -> Dict[str, set]:
    """Genes with minor-intron retention that the paired control does not share.

    A gene is kept when at least one minor intron is a hit and every control
    intron with that minor intron's pair_id is not a hit. A hit on the matched
    U2 intron is widespread retention and does not count. A second minor intron
    in the same gene whose control is also a hit does not remove the gene.
    """
    indexed = meta.reindex(list(target_ids))
    genes = indexed["gene_name"].to_numpy()
    pairs = indexed["pair_id"].astype(str).to_numpy()
    minor = np.asarray(is_minor, dtype=bool)
    called = np.asarray(hit, dtype=bool)
    control = ~minor & ~np.isin(pairs, ["", "nan"])
    out: Dict[str, set] = {}
    for i, sample in enumerate(samples):
        row = called[i]
        blocked = set(pairs[row & control])
        keep = row & minor
        if blocked:
            keep = keep & ~np.isin(pairs, list(blocked))
        names = genes[keep]
        out[str(sample)] = {
            str(g) for g in names
            if pd.notna(g) and str(g).strip() and str(g).lower() != "nan"
        }
    return out


def _diagnosis_groups(scores: pd.DataFrame) -> List[tuple]:
    """One group per diagnosis label that is not a placeholder.

    RNU12, RNU4ATAC, and RNU6ATAC are included, and so is any other text.
    Blank, 0, control, N/A, NA, null, and none are left out. Spelling that
    differs only by case is one group.
    """
    buckets: Dict[str, dict] = {}
    for sample, value in zip(scores["sample"].astype(str), scores["diagnosis"]):
        if _is_blank_label(value):
            continue
        label = _cell_text(value)
        key = label.casefold()
        bucket = buckets.setdefault(key, {"spellings": [], "ids": []})
        bucket["spellings"].append(label)
        bucket["ids"].append(sample)
    groups = []
    for key in sorted(buckets):
        bucket = buckets[key]
        label = _choose_spelling(bucket["spellings"])
        groups.append((label, label, bucket["ids"]))
    return groups


def upset_groups(scores: pd.DataFrame) -> List[tuple]:
    """Groups that have at least one individual. Empty diagnoses are omitted."""
    if "sample" not in scores.columns:
        return []
    groups = []
    if "minor_spliceopathy" in scores.columns:
        ids = scores.loc[_spliceopathy_mask(scores["minor_spliceopathy"]),
                         "sample"].astype(str).tolist()
        if ids:
            groups.append(("spliceopathy", "minor spliceopathy", ids))
    if "diagnosis" in scores.columns:
        groups.extend(_diagnosis_groups(scores))
    return groups


def upset_intersections(gene_sets: Dict[str, set],
                        sample_ids: Sequence[str]) -> List[dict]:
    """One row per distinct combination of individuals that share a gene."""
    buckets: Dict[str, set] = {}
    for sample in sample_ids:
        for gene in gene_sets.get(sample, ()):
            buckets.setdefault(gene, set()).add(sample)
    grouped: Dict[tuple, List[str]] = {}
    for gene, members in buckets.items():
        key = tuple(s for s in sample_ids if s in members)
        grouped.setdefault(key, []).append(gene)
    rows = []
    for members, genes in grouped.items():
        rows.append({
            "samples": members,
            "n_individuals": len(members),
            "genes": tuple(sorted(genes)),
            "n_genes": len(genes),
        })
    rows.sort(key=lambda r: (-r["n_genes"], -r["n_individuals"], r["samples"]))
    return rows


def fig_upset(gene_sets: Dict[str, set], sample_ids: Sequence[str],
              outdir: str, filename: str, title: str,
              max_intersections: int = 40) -> str:
    """UpSet of genes with minor-specific retention. Each individual is one set."""
    rows = upset_intersections(gene_sets, list(sample_ids))
    shown = rows[:max_intersections]
    n_sets = max(len(sample_ids), 1)
    n_int = max(len(shown), 1)
    # Fixed inches per column and per person, so a wide group has the same
    # bar and dot size as a small group. The report scrolls the wide figure.
    left_in, col_in = 2.4, 0.72
    top_in, row_in = 2.3, 0.48
    fig = plt.figure(
        figsize=(left_in + col_in * n_int, top_in + row_in * n_sets),
        layout="constrained")
    gs = fig.add_gridspec(
        2, 2, width_ratios=[left_in, col_in * n_int],
        height_ratios=[top_in, row_in * n_sets])
    ax_size = fig.add_subplot(gs[1, 0])
    ax_dots = fig.add_subplot(gs[1, 1])
    ax_bars = fig.add_subplot(gs[0, 1], sharex=ax_dots)
    ax_size.sharey(ax_dots)

    ys = np.arange(n_sets)
    sizes = [len(gene_sets.get(s, ())) for s in sample_ids]
    ax_size.barh(ys, sizes, color="#4a5568", height=0.62)
    ax_size.invert_xaxis()
    ax_size.tick_params(labelleft=False, left=False)
    ax_size.set_xlabel("genes with minor-specific retention", fontsize=8)
    ax_size.spines[["top", "right"]].set_visible(False)
    ax_dots.set_yticks(ys)
    ax_dots.set_yticklabels(list(sample_ids), fontsize=8)
    ax_dots.set_xlim(-0.6, n_int - 0.4)
    ax_dots.set_ylim(-0.5, n_sets - 0.5)
    ax_dots.invert_yaxis()

    if shown:
        for x, row in enumerate(shown):
            member_y = [sample_ids.index(s) for s in row["samples"]]
            shared = row["n_individuals"] >= 2
            color = "#2b6cb0" if shared else "#a0aec0"
            if len(member_y) >= 2:
                ax_dots.plot([x, x], [min(member_y), max(member_y)],
                             color=color, lw=1.6, zorder=1)
            ax_dots.scatter([x] * len(member_y), member_y, s=36, color=color,
                            zorder=2)
        counts = [row["n_genes"] for row in shown]
        colors = ["#2b6cb0" if row["n_individuals"] >= 2 else "#a0aec0"
                  for row in shown]
        ax_bars.bar(np.arange(len(shown)), counts, color=colors, width=0.72)
        for x, n in enumerate(counts):
            ax_bars.text(x, n, str(n), ha="center", va="bottom", fontsize=7,
                         color="#2d3748")
        ax_bars.set_ylim(0, max(counts) * 1.18 + 1)
        ax_bars.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    else:
        ax_bars.text(0.5, 0.5, "no minor-specific genes", ha="center", va="center",
                     transform=ax_bars.transAxes, color="#718096", fontsize=10)

    ax_bars.tick_params(labelbottom=False, bottom=False)
    ax_bars.set_ylabel("genes", fontsize=8)
    ax_bars.spines[["top", "right"]].set_visible(False)
    ax_dots.set_xticks([])
    ax_dots.tick_params(bottom=False)
    ax_dots.spines[["top", "right", "bottom"]].set_visible(False)
    note = ""
    if len(rows) > len(shown):
        note = f" (tallest {len(shown)} of {len(rows)} combinations)"
    ax_bars.set_title(title + note, fontsize=11)
    path = os.path.join(outdir, filename)
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return path


def singleton_kind(group_id: str, group_label: str, diagnosis: str) -> str:
    """Label a one-person group by diagnosis when that diagnosis is a real label."""
    if group_id != "spliceopathy":
        return group_id or group_label
    if not _is_blank_label(diagnosis):
        return _cell_text(diagnosis)
    return group_label


def upset_catalog(scores: pd.DataFrame, gene_sets: Dict[str, set]) -> List[dict]:
    """One entry per sample, labeled by diagnosis and spliceopathy status."""
    if "sample" not in scores.columns:
        return []
    splice = ( _spliceopathy_mask(scores["minor_spliceopathy"])
               if "minor_spliceopathy" in scores.columns
               else pd.Series(False, index=scores.index))
    diag_labels = (
        _label_spellings(scores["diagnosis"]) if "diagnosis" in scores.columns
        else {})
    splice_labels = (
        _label_spellings(scores["minor_spliceopathy"])
        if "minor_spliceopathy" in scores.columns else {})
    rows = []
    for i, record in scores.iterrows():
        sample = str(record["sample"])
        diagnosis = ""
        if "diagnosis" in scores.columns:
            diagnosis = diag_labels.get(
                _cell_text(record["diagnosis"]).casefold(), "")
        is_ms = bool(splice.loc[i])
        splice_label = "minor spliceopathy" if is_ms else ""
        if not splice_label and "minor_spliceopathy" in scores.columns:
            splice_label = splice_labels.get(
                _cell_text(record["minor_spliceopathy"]).casefold(), "")
        parts = [p for p in (diagnosis, splice_label) if p]
        label = sample if not parts else f"{sample} — {', '.join(parts)}"
        rows.append({
            "id": sample,
            "label": label,
            "diagnosis": diagnosis,
            "spliceopathy": is_ms,
            "genes": sorted(gene_sets.get(sample, ())),
        })
    rows.sort(key=lambda d: (not d["spliceopathy"], d["diagnosis"], d["id"]))
    return rows


def singleton_gene_html(sample: str, kind: str, genes: Sequence[str]) -> str:
    """Gene list used in place of an UpSet plot for a one-person group."""
    names = list(genes)
    items = "".join(f"<li>{html.escape(g)}</li>" for g in names)
    listing = (f"<ul class='genelist'>{items}</ul>" if names
               else "<div class='note'>No genes with a minor-intron hit "
                    "whose paired control is not a hit.</div>")
    return (f"<h3>{html.escape(sample)} ({html.escape(kind)})</h3>"
            f"<div class='sub'>{len(names)} genes with a minor-intron hit "
            "and no hit on the paired control intron.</div>"
            + listing)


def upset_shared_table(rows: Sequence[dict]) -> pd.DataFrame:
    """Genes found in two or more individuals of one group."""
    shared = [r for r in rows if r["n_individuals"] >= 2]
    if not shared:
        return pd.DataFrame(columns=["individuals", "n_individuals", "n_genes",
                                     "genes"])
    return pd.DataFrame({
        "individuals": [", ".join(r["samples"]) for r in shared],
        "n_individuals": [r["n_individuals"] for r in shared],
        "n_genes": [r["n_genes"] for r in shared],
        "genes": [", ".join(r["genes"]) for r in shared],
    })


# --------------------------------------------------------------------------- #
# interactive PCA
# --------------------------------------------------------------------------- #
_PCA_JS = r"""
(function(){
  const D = __DATA__, VAR = __VAR__, SPEC = __SPEC__;
  const cv = document.getElementById('pca'), tip = document.getElementById('tip');
  const ctx = cv.getContext('2d');
  let view = null, hi = new Set(), pinned = new Set(), metric = 'n_hits';
  let hover = -1, drag = null, catOrder = [];

  const PAL = ['#c1272d','#2b6cb0','#2f855a','#b7791f','#6b46c1','#c05621',
               '#2c7a7b','#97266d','#2d3748','#dd6b20','#3182ce','#38a169',
               '#d53f8c','#805ad5','#4a5568','#9b2c2c','#2a4365','#276749',
               '#744210','#44337a'];
  const finite = a => a.filter(v => v !== null && isFinite(v));
  function esc(s){
    return String(s).replace(/[&<>"']/g, ch => ({
      '&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[ch]));
  }
  function ramp(t){
    t = Math.max(0, Math.min(1, t));
    const stops = [[237,242,247],[144,180,214],[ 99,133,186],[214,120, 86],[193, 39, 45]];
    const x = t*(stops.length-1), i = Math.min(Math.floor(x), stops.length-2), f = x-i;
    const a = stops[i], b = stops[i+1];
    return `rgb(${Math.round(a[0]+f*(b[0]-a[0]))},${Math.round(a[1]+f*(b[1]-a[1]))},${Math.round(a[2]+f*(b[2]-a[2]))})`;
  }
  function metaOf(d){ return d.meta || {}; }
  function catKey(d){
    const v = metaOf(d)[metric];
    if(v === null || v === undefined || v === '') return 'NA';
    return String(v);
  }
  function catFill(key){
    if(key === 'NA') return '#cbd5e0';
    const i = catOrder.indexOf(key);
    if(i < 0) return '#cbd5e0';
    if(i < PAL.length) return PAL[i];
    return `hsl(${Math.round((i * 137.508) % 360)},62%,42%)`;
  }
  function metricVal(d){
    if(SPEC[metric] === 'num'){
      const v = metaOf(d)[metric];
      return (v === null || v === undefined || v === '') ? null : +v;
    }
    if(SPEC[metric] === 'cat') return null;
    if(metric === 'tail_q') return (d.q === null || d.q <= 0) ? 320 : -Math.log10(d.q);
    if(metric === 'enrichment') return d.e === null ? null : Math.log2(Math.max(d.e, 0.05));
    return d[{n_hits:'h', MSRI:'m', stouffer:'st'}[metric]];
  }
  let lo = 0, hiV = 1;
  function rescale(){
    if(SPEC[metric] === 'cat'){
      const seen = new Set();
      catOrder = [];
      D.forEach(d => { const k = catKey(d); if(!seen.has(k)){ seen.add(k); catOrder.push(k); } });
      catOrder.sort((a,b) => {
        if(a === 'NA') return 1;
        if(b === 'NA') return -1;
        return a.localeCompare(b, undefined, {numeric:true, sensitivity:'base'});
      });
      return;
    }
    catOrder = [];
    const v = finite(D.map(metricVal));
    v.sort((a,b)=>a-b);
    lo = v.length ? v[Math.floor(0.02*(v.length-1))] : 0;
    hiV = v.length ? v[Math.floor(0.98*(v.length-1))] : 1;
    // A handful of extreme samples among hundreds makes the 98th percentile
    // collapse onto the minimum; fall back to the true maximum so the ramp
    // still carries information.
    if(hiV <= lo && v.length) hiV = v[v.length-1];
    if(hiV <= lo) hiV = lo + 1;
  }
  function resetView(){
    const xs = D.map(d=>d.x), ys = D.map(d=>d.y);
    const x0=Math.min(...xs), x1=Math.max(...xs), y0=Math.min(...ys), y1=Math.max(...ys);
    const px=(x1-x0)*0.08+1e-9, py=(y1-y0)*0.08+1e-9;
    view = {x0:x0-px, x1:x1+px, y0:y0-py, y1:y1+py};
    draw();
  }
  function size(){
    const r = cv.getBoundingClientRect(), dpr = window.devicePixelRatio||1;
    cv.width = r.width*dpr; cv.height = r.height*dpr;
    ctx.setTransform(dpr,0,0,dpr,0,0);
    return r;
  }
  function sx(x,r){ return (x-view.x0)/(view.x1-view.x0)*(r.width-70)+55; }
  function sy(y,r){ return r.height-38-(y-view.y0)/(view.y1-view.y0)*(r.height-68); }

  function draw(){
    const r = size();
    ctx.clearRect(0,0,r.width,r.height);
    // axes
    ctx.strokeStyle='#edf2f7'; ctx.lineWidth=1;
    if(view.x0<0&&view.x1>0){ctx.beginPath();ctx.moveTo(sx(0,r),12);ctx.lineTo(sx(0,r),r.height-38);ctx.stroke();}
    if(view.y0<0&&view.y1>0){ctx.beginPath();ctx.moveTo(55,sy(0,r));ctx.lineTo(r.width-15,sy(0,r));ctx.stroke();}
    ctx.fillStyle='#4a5568'; ctx.font='12px system-ui,sans-serif';
    ctx.textAlign='center'; ctx.fillText(VAR.xlab, r.width/2, r.height-10);
    ctx.save(); ctx.translate(14,r.height/2); ctx.rotate(-Math.PI/2);
    ctx.fillText(VAR.ylab,0,0); ctx.restore();

    const labelled = [];
    for(let i=0;i<D.length;i++){
      const d=D[i], X=sx(d.x,r), Y=sy(d.y,r);
      if(X<40||X>r.width||Y<0||Y>r.height) continue;
      const on = hi.has(i)||pinned.has(i), v = metricVal(d);
      ctx.beginPath();
      ctx.arc(X,Y, on?6.5:(hover===i?5.5:3.4), 0, 6.2832);
      ctx.fillStyle = SPEC[metric] === 'cat' ? catFill(catKey(d))
        : ((v===null||!isFinite(v)) ? '#cbd5e0' : ramp((v-lo)/(hiV-lo)));
      ctx.globalAlpha = (hi.size&&!on) ? 0.28 : 0.9;
      ctx.fill();
      if(on){ ctx.globalAlpha=1; ctx.lineWidth=1.8; ctx.strokeStyle='#111'; ctx.stroke();
              labelled.push([i,X,Y]); }
      ctx.globalAlpha=1;
    }
    ctx.font='bold 12px system-ui,sans-serif'; ctx.textAlign='left';
    for(const [i,X,Y] of labelled.slice(0,40)){
      const t = D[i].sh+'  ('+D[i].h+')', w = ctx.measureText(t).width;
      let lx = X+11, ly = Y-9;
      if(lx+w+10>r.width) lx = X-w-13;
      ctx.fillStyle='rgba(255,255,255,.93)';
      ctx.strokeStyle='#111'; ctx.lineWidth=.8;
      ctx.beginPath(); ctx.roundRect(lx-5,ly-13,w+10,19,5); ctx.fill(); ctx.stroke();
      ctx.fillStyle='#111'; ctx.fillText(t,lx,ly+1);
    }
    const legend = document.getElementById('legend');
    if(SPEC[metric] === 'cat'){
      legend.innerHTML = catOrder.map(k =>
        `<span class="sw" style="background:${catFill(k)}"></span>${esc(k)}`
      ).join('') + `&nbsp;&nbsp;<span style="color:#718096">${esc(metric)}</span>`;
    } else {
      legend.innerHTML =
        `<span class="sw" style="background:${ramp(0)}"></span>${VAR.fmt(lo)}` +
        `<span class="sw" style="background:${ramp(.5)}"></span>` +
        `<span class="sw" style="background:${ramp(1)}"></span>${VAR.fmt(hiV)}` +
        `&nbsp;&nbsp;<span style="color:#718096">${esc(metric)}</span>`;
    }
  }

  function nearest(mx,my){
    const r = cv.getBoundingClientRect();
    let best=-1, bd=225;
    for(let i=0;i<D.length;i++){
      const dx=sx(D[i].x,r)-mx, dy=sy(D[i].y,r)-my, d2=dx*dx+dy*dy;
      if(d2<bd){ bd=d2; best=i; }
    }
    return best;
  }
  cv.addEventListener('mousemove', ev=>{
    const r=cv.getBoundingClientRect(), mx=ev.clientX-r.left, my=ev.clientY-r.top;
    if(drag){
      const kx=(view.x1-view.x0)/(r.width-70), ky=(view.y1-view.y0)/(r.height-68);
      view.x0-=(mx-drag.x)*kx; view.x1-=(mx-drag.x)*kx;
      view.y0+=(my-drag.y)*ky; view.y1+=(my-drag.y)*ky;
      drag={x:mx,y:my}; draw(); return;
    }
    const i=nearest(mx,my);
    if(i!==hover){ hover=i; draw(); }
    if(i>=0){
      const d=D[i];
      const metaBits = Object.keys(metaOf(d)).filter(k => k !== 'run').map(k => {
        const v = metaOf(d)[k];
        const shown = (v===null||v===undefined||v==='') ? 'NA' : v;
        return `${esc(k)}: ${esc(shown)}`;
      });
      tip.innerHTML = `<b>${esc(d.n)}</b><br>run: ${esc(d.r)}<br>minor hits: <b>${d.h}</b>`+
        ` &nbsp;control hits: ${d.hc}<br>enrichment: ${d.e===null?'NA':d.e.toFixed(1)}`+
        `&times; &nbsp;q: ${d.q===null?'NA':d.q.toExponential(1)}<br>`+
        `MSRI ${d.m===null?'NA':d.m.toFixed(2)} &nbsp;PC1 ${d.x.toFixed(1)} &nbsp;PC2 ${d.y.toFixed(1)}`+
        (metaBits.length ? '<br>'+metaBits.join('<br>') : '');
      tip.style.display='block';
      tip.style.left=Math.min(mx+16, r.width-250)+'px';
      tip.style.top=(my+14)+'px';
    } else tip.style.display='none';
  });
  cv.addEventListener('mouseleave', ()=>{ tip.style.display='none'; hover=-1; drag=null; draw(); });
  cv.addEventListener('mousedown', ev=>{
    const r=cv.getBoundingClientRect(), mx=ev.clientX-r.left, my=ev.clientY-r.top;
    const i=nearest(mx,my);
    if(i>=0){ pinned.has(i)?pinned.delete(i):pinned.add(i); draw(); }
    else drag={x:mx,y:my};
  });
  window.addEventListener('mouseup', ()=>{ drag=null; });
  cv.addEventListener('wheel', ev=>{
    ev.preventDefault();
    const r=cv.getBoundingClientRect(), mx=ev.clientX-r.left, my=ev.clientY-r.top;
    const f = ev.deltaY>0 ? 1.15 : 1/1.15;
    const cx = view.x0+(mx-55)/(r.width-70)*(view.x1-view.x0);
    const cy = view.y0+(r.height-38-my)/(r.height-68)*(view.y1-view.y0);
    view = {x0:cx+(view.x0-cx)*f, x1:cx+(view.x1-cx)*f,
            y0:cy+(view.y0-cy)*f, y1:cy+(view.y1-cy)*f};
    draw();
  }, {passive:false});

  const box = document.getElementById('search');
  function runSearch(){
    const q = box.value.trim().toLowerCase();
    hi = new Set();
    if(q){
      const terms = q.split(/[\s,;]+/).filter(Boolean);
      D.forEach((d,i)=>{
        const extra = Object.values(metaOf(d)).join(' ');
        const n=(d.n+' '+d.r+' '+extra).toLowerCase();
        if(terms.some(t=>n.includes(t))) hi.add(i);
      });
    }
    document.getElementById('nfound').textContent =
      q ? `${hi.size} match${hi.size===1?'':'es'}` : '';
    draw();
  }
  box.addEventListener('input', runSearch);
  document.getElementById('metric').addEventListener('change', ev=>{
    metric = ev.target.value; rescale(); draw();
  });
  document.getElementById('topn').addEventListener('click', ()=>{
    const idx = D.map((d,i)=>[d.h,i]).sort((a,b)=>b[0]-a[0]).slice(0,10).map(p=>p[1]);
    pinned = new Set(idx); draw();
  });
  document.getElementById('clear').addEventListener('click', ()=>{
    pinned = new Set(); hi = new Set(); box.value=''; 
    document.getElementById('nfound').textContent=''; resetView();
  });
  window.addEventListener('resize', draw);
  rescale(); resetView();
})();
"""

CSS = """
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
 margin:0 auto;max-width:1180px;padding:26px;color:#1a202c;line-height:1.5}
h1{font-size:22px;margin-bottom:2px}
h2{font-size:16px;margin-top:34px;border-bottom:1px solid #e2e8f0;padding-bottom:5px;
 scroll-margin-top:72px}
.toc{position:sticky;top:0;z-index:20;background:#fff;border:1px solid #e2e8f0;
 border-radius:6px;padding:8px 12px;margin:14px 0 4px;display:flex;flex-wrap:wrap;
 gap:4px 14px;box-shadow:0 1px 4px rgba(0,0,0,.06)}
.toc a{font-size:12.5px;color:#2b6cb0;text-decoration:none}
.toc a:hover{text-decoration:underline}
.sub{color:#718096;font-size:12.5px;margin-bottom:16px}
table{border-collapse:collapse;font-size:11.5px;margin:10px 0;width:100%}
th,td{border:1px solid #e2e8f0;padding:4px 7px;text-align:right}
th{background:#f7fafc;font-weight:600;text-align:center}
td:first-child,th:first-child{text-align:left}
img{max-width:100%;margin:10px 0;border:1px solid #edf2f7;border-radius:4px}
.hscroll{overflow-x:auto;overflow-y:hidden;max-width:100%;margin:10px 0;
 border:1px solid #edf2f7;border-radius:4px}
.hscroll img{max-width:none;width:auto;height:auto;margin:0;border:none;
 border-radius:0}
.note{background:#f7fafc;border-left:3px solid #4a5568;padding:9px 13px;
 font-size:12.5px;margin:14px 0}
.warn{background:#fffaf0;border-left:3px solid #b7791f}
code{background:#f1f5f9;padding:1px 4px;border-radius:3px;font-size:11px}
.wrap{position:relative;border:1px solid #e2e8f0;border-radius:6px;padding:6px}
.wrap canvas{width:100%;height:400px;display:block;cursor:crosshair}
.wrap canvas.tall{height:520px}
.tip{position:absolute;display:none;background:#fff;border:1px solid #cbd5e0;
 border-radius:5px;padding:7px 10px;font-size:11.5px;pointer-events:none;
 box-shadow:0 2px 8px rgba(0,0,0,.13);z-index:5;max-width:360px}
.bar{display:flex;gap:9px;align-items:center;flex-wrap:wrap;margin:6px 0 10px}
.bar input[type=text]{flex:1;min-width:220px;padding:6px 10px;font-size:13px;
 border:1px solid #cbd5e0;border-radius:5px}
.bar select,.bar button{padding:6px 10px;font-size:12.5px;border:1px solid #cbd5e0;
 border-radius:5px;background:#fff;cursor:pointer}
.bar button.on{background:#2b6cb0;color:#fff;border-color:#2b6cb0}
.pinlab{font-size:12.5px;color:#4a5568}
.pinlab input{width:3.4em;margin-left:4px;padding:4px 6px;font-size:12.5px;
 border:1px solid #cbd5e0;border-radius:5px}
.sw{display:inline-block;width:20px;height:11px;margin:0 4px;vertical-align:-1px;
 border:1px solid #cbd5e0}
.legend{font-size:11.5px;color:#4a5568;margin-top:4px;line-height:1.8}
.hm-legend{max-width:440px;margin:4px 0 16px}
.hm-bar{height:14px;border:1px solid #cbd5e0;border-radius:2px;
 background:linear-gradient(90deg,#053061,#2065ab,#4393c3,#90c4dd,#d1e5f0,#f7f6f6,#fddbc7,#f3a481,#d6604d,#b1182b,#67001f)}
.hm-ticks{display:flex;justify-content:space-between;font-size:11.5px;color:#4a5568;
 margin-top:3px}
.hm-cap{font-size:12px;color:#4a5568;margin-top:2px}
.genelist{columns:4;column-gap:28px;font-size:12.5px;margin:8px 0 16px;padding-left:1.2em}
.genelist li{break-inside:avoid;padding:1px 0}
.upset-tools{display:flex;gap:8px;flex-wrap:wrap;margin:8px 0}
.upset-tools button.on{background:#2b6cb0;color:#fff;border-color:#2b6cb0}
select.upset-pick{width:100%;max-width:720px;font-size:13px;padding:4px;
 border:1px solid #cbd5e0;border-radius:5px;background:#fff}
.upset-scroll{overflow-x:auto;margin:8px 0}
h3{font-size:14px;margin:8px 0 4px}
.nfound{font-size:12px;color:#718096;min-width:74px}
.cap{font-size:12px;color:#4a5568;margin:0 0 4px}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:12px}
.grid3{display:grid;grid-template-columns:1fr 1fr 1fr;gap:10px}
@media (max-width:980px){.grid2,.grid3{grid-template-columns:1fr}}
"""


_N_HITS = (
    "<b>n_hits</b> is how many minor introns are called hits in that sample."
)
_N_HITS_CONTROL = (
    "<b>n_hits_control</b> is that same hit rule counted on the sample's own "
    "U2 control introns. A minor-spliceosome signal raises <b>n_hits</b> above "
    "<b>n_hits_control</b>. A noisy library raises both."
)
_TAIL_Q = (
    "<b>tail_q</b> is the Benjamini–Hochberg false-discovery rate, across the "
    "samples, of the one-sided Fisher test that this sample's minor-intron hit "
    "rate is higher than its own control-intron hit rate. A smaller "
    "<b>tail_q</b> is stronger evidence. Samples with fewer than 20 measurable "
    "introns of either class have no <b>tail_q</b>."
)
_ENRICHMENT = (
    "<b>enrichment</b> is the minor-intron hit rate divided by the "
    "control-intron hit rate in that same sample. Each rate is "
    "(hits + 0.5) / (measurable introns of that class + 1). Those added "
    "counts keep a zero from making the ratio zero or infinite. A value near "
    "1 means minor introns are called about as often as that sample's own "
    "control introns. A larger value means the hits are concentrated on minor "
    "introns. Samples with fewer than 20 measurable introns of either class "
    "have no <b>enrichment</b>."
)


def _section(toc: List[tuple], anchor: str, title: str) -> str:
    toc.append((anchor, title))
    return f"<h2 id='{html.escape(anchor, quote=True)}'>{html.escape(title)}</h2>"


def _contents(toc: List[tuple]) -> str:
    links = "".join(
        f"<a href='#{html.escape(anchor, quote=True)}'>{html.escape(title)}</a>"
        for anchor, title in toc)
    return f"<nav class='toc' aria-label='Contents'>{links}</nav>"


def _defs(*keys: str) -> str:
    text = {"n_hits": _N_HITS, "n_hits_control": _N_HITS_CONTROL,
            "tail_q": _TAIL_Q, "enrichment": _ENRICHMENT}
    return " ".join(text[k] for k in keys)


def _param_value(ctx: dict, key: str, default: str) -> str:
    params = ctx.get("params")
    if params is None or "parameter" not in getattr(params, "columns", []):
        return default
    rows = params.loc[params["parameter"] == key, "value"]
    if len(rows) == 0 or pd.isna(rows.iloc[0]):
        return default
    return str(rows.iloc[0])


def _as_float(value) -> Optional[float]:
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


def _fmt_cut(value) -> str:
    number = _as_float(value)
    if number is None or number != number or abs(number) == float("inf"):
        return html.escape(str(value))
    if number == int(number) and abs(number) < 1e12:
        return str(int(number))
    return f"{number:.4g}"


def _fold_clause(value, middle: str) -> str:
    """log2 cutoff as a fold change over the cohort middle, when it is finite."""
    number = _as_float(value)
    if number is None or number != number or abs(number) == float("inf"):
        return ""
    shown = f"{2 ** number:.2f}".rstrip("0").rstrip(".")
    return f", about {shown} times the cohort-middle {middle}"


def _hit_definition(ctx: dict) -> str:
    """What an intron call requires in the hit mode used for this report."""
    mode = str(ctx.get("hit_mode") or _param_value(ctx, "hit_mode", "coverage"))
    z = _fmt_cut(_param_value(ctx, "z_threshold", "3"))
    d_cov = _param_value(ctx, "min_delta", "0.5")
    d_j = _param_value(ctx, "min_delta_junction", "0.5")
    veto = _param_value(ctx, "junction_veto", "1")
    d_cov_s = _fmt_cut(d_cov)
    d_j_s = _fmt_cut(d_j)
    cov_fold = _fold_clause(d_cov, "intron-to-exon ratio")
    j_fold = _fold_clause(d_j, "unspliced-to-spliced ratio")
    same = (" The same rule is applied to minor introns and to U2 control "
            "introns.")
    if mode == "coverage":
        return (
            "A <b>hit</b> is one intron in one sample. This report uses hit "
            f"mode <b>coverage</b>, and both checks have to pass. <b>z_cov</b>, "
            "how unusual this sample's adjusted intron-to-exon coverage is "
            f"for this intron, is at least <b>{z}</b>. The size of that shift, "
            f"<b>d_cov</b>, is at least <b>{d_cov_s}</b>{cov_fold}. A missing "
            "score on either check keeps the intron out of the hit set."
            + same)
    if mode == "junction":
        return (
            "A <b>hit</b> is one intron in one sample. This report uses hit "
            f"mode <b>junction</b>, and both checks have to pass. <b>z_j</b>, "
            "how unusual this sample's adjusted unspliced-to-spliced ratio is "
            f"at the two splice sites, is at least <b>{z}</b>. The size of "
            f"that shift, <b>d_j</b>, is at least <b>{d_j_s}</b>{j_fold}. A "
            "missing score on either check keeps the intron out of the hit set."
            + same)
    if mode == "both":
        return (
            "A <b>hit</b> is one intron in one sample. This report uses hit "
            f"mode <b>both</b>, so the coverage rule and the junction rule "
            f"both have to pass. Coverage requires <b>z_cov</b> at least "
            f"<b>{z}</b> and <b>d_cov</b> at least <b>{d_cov_s}</b>{cov_fold}. "
            f"Junctions require <b>z_j</b> at least <b>{z}</b> and <b>d_j</b> "
            f"at least <b>{d_j_s}</b>{j_fold}. A missing score on either rule "
            "keeps the intron out of the hit set." + same)
    if mode == "combined":
        veto_n = _as_float(veto)
        if veto_n is not None and veto_n == float("-inf"):
            veto_text = (" The junction veto is off for this report, so a low "
                         "splice-site score does not block the call.")
        else:
            veto_text = (
                " When splice-site reads are measurable, <b>z_j</b> is at "
                f"least <b>{_fmt_cut(veto)}</b>. A lower splice-site score "
                "keeps the intron out of the hit set, because reads in the "
                "intron body with ordinary splice sites count as background. "
                "When those reads are too few to score, that check is unused "
                "and a coverage-driven <b>z_comb</b> can still be a hit.")
        return (
            "A <b>hit</b> is one intron in one sample. This report uses hit "
            f"mode <b>combined</b>, and all three checks have to pass. "
            f"<b>z_comb</b> is at least <b>{z}</b>. That score merges how "
            "unusual the intron-body coverage is with how unusual the splice "
            "sites are, and it is whichever of those two scores is measurable "
            "when only one exists. The shift is large on at least one measure: "
            f"<b>d_cov</b> is at least <b>{d_cov_s}</b>{cov_fold}, or "
            f"<b>d_j</b> is at least <b>{d_j_s}</b>{j_fold}."
            + veto_text + same)
    return (
        "A <b>hit</b> is one intron in one sample, called with hit mode "
        f"<b>{html.escape(mode)}</b>. The cutoffs for that call are in "
        "Parameters." + same)


def _pca_score_name(mode: str) -> str:
    return {"coverage": "z_cov", "junction": "z_j"}.get(mode, "z_comb")


def _top_events_scoring(ctx: dict) -> str:
    """How a row gets into the top-events table, and what its z is."""
    mode = str(ctx.get("hit_mode") or _param_value(ctx, "hit_mode", "coverage"))
    score = _pca_score_name(mode)
    if mode == "both":
        z_text = (
            f"The <b>z</b> column is <b>{score}</b>, the merged coverage and "
            "junction score. The call itself still requires both rules above, "
            "and the rows are ordered by this merged score.")
    elif mode == "junction":
        z_text = (
            f"The <b>z</b> column is <b>{score}</b>, the splice-site score the "
            "hit rule tests. The rows are ordered by that score.")
    elif mode == "coverage":
        z_text = (
            f"The <b>z</b> column is <b>{score}</b>, the intron-body score the "
            "hit rule tests. The rows are ordered by that score.")
    elif mode == "combined":
        z_text = (
            f"The <b>z</b> column is <b>{score}</b>, the merged score the hit "
            "rule tests. The rows are ordered by that score.")
    else:
        z_text = (
            f"The <b>z</b> column is <b>{html.escape(score)}</b>. The rows are "
            "ordered by that score.")
    return (
        "An <b>event</b> is one minor intron in one sample that passed the "
        "hit rule above. Control-intron hits are not listed. " + z_text +
        " This table keeps the 60 events with the largest <b>z</b> in the "
        "cohort, so one sample can occupy several rows. Per-intron calls are "
        "much less reliable than the aggregate score; treat them as candidates "
        "to confirm (IGV, then RT-PCR). <code>n_EI_donor</code> and "
        "<code>n_EI_acceptor</code> are reads crossing each splice boundary. "
        "Retention without them is more likely an overlapping transcript or a "
        "mapping artefact.")


def _pca_position(ctx: dict) -> str:
    """How a point is placed, as distinct from the hit calls in section 4."""
    mode = str(ctx.get("hit_mode") or _param_value(ctx, "hit_mode", "coverage"))
    score = _pca_score_name(mode)
    z = _fmt_cut(_param_value(ctx, "z_threshold", "3"))
    clip = _as_float(_param_value(ctx, "pca_clip", "8"))
    if mode == "combined":
        still = (" Introns below the size cutoff in section 4, and introns "
                 "kept out by the splice-site veto, still move the point.")
    elif mode == "coverage":
        still = (" Introns below the <b>d_cov</b> size cutoff in section 4 "
                 "still move the point.")
    elif mode == "junction":
        still = (" Introns below the <b>d_j</b> size cutoff in section 4 "
                 "still move the point.")
    elif mode == "both":
        still = (" An intron moves the point from its <b>z_comb</b> even when "
                 "it fails the coverage rule or the junction rule in section 4.")
    else:
        still = (" Introns that fail the hit rule in section 4 still move "
                 "the point.")
    if clip is not None and clip == clip and abs(clip) != float("inf") and clip > 0:
        lo = _fmt_cut(-clip / 2.0).replace("-", "−")
        clip_text = (
            f" Extreme values are clipped from {lo} to {_fmt_cut(clip)} "
            "before the axes are fit.")
    else:
        clip_text = ""
    return (
        "Each point is placed from the continuous score of every measurable "
        f"minor intron. This report uses <b>{score}</b>, the same score that "
        f"has to reach <b>{z}</b> before an intron can be a hit. The map keeps "
        "that full value." + still + " Control introns are left out of this "
        "map. A missing score is filled with that intron's cohort median."
        + clip_text + " Colour by <b>n_hits</b>, and Pin top N, count only the "
        "minor introns that pass the hit rule in section 4.")


def _pca_caption(ctx: dict) -> str:
    mode = str(ctx.get("hit_mode") or _param_value(ctx, "hit_mode", "coverage"))
    score = _pca_score_name(mode)
    return (
        "Each point is one sample, placed from the continuous "
        f"<b>{score}</b> of every measurable minor intron. Introns that fail "
        "the hit rule in section 4 still move the point. Further right on PC1 "
        "means more coordinated retention across minor introns. Further up "
        "on PC2 is a second pattern in that profile, often a batch or a "
        "different subset of introns.")


def _table(df: pd.DataFrame, fmt: str = "{:.3g}") -> str:
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda v: "" if pd.isna(v) else fmt.format(v))
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in d.columns)
    body = "".join("<tr>" + "".join(f"<td>{html.escape(str(v))}</td>" for v in row)
                   + "</tr>" for row in d.itertuples(index=False))
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


_DEFAULT_COLOURS = [
    {"id": "n_hits", "label": "colour: n_hits", "kind": "score"},
    {"id": "tail_q", "label": "colour: -log10(tail_q)", "kind": "score"},
    {"id": "enrichment", "label": "colour: log2(enrichment)", "kind": "score"},
    {"id": "MSRI", "label": "colour: MSRI", "kind": "score"},
    {"id": "stouffer", "label": "colour: stouffer", "kind": "score"},
]


def _js_json(obj) -> str:
    # Keep a value like "</script>" from closing the inline script.
    return json.dumps(obj, separators=(",", ":")).replace("<", "\\u003c")


_PLOTS_JS = r"""
(function(){
  const D = __DATA__, SPEC = __SPEC__, CFG = __CFG__;
  const PAL = ['#c1272d','#2b6cb0','#2f855a','#b7791f','#6b46c1','#c05621',
               '#2c7a7b','#97266d','#2d3748','#dd6b20','#3182ce','#38a169',
               '#d53f8c','#805ad5','#4a5568','#9b2c2c','#2a4365','#276749',
               '#744210','#44337a'];
  const rank = new Map();
  D.map((d,i)=>[d.h,i]).sort((a,b)=>b[0]-a[0] || a[1]-b[1])
    .forEach((p,k)=>rank.set(p[1], k+1));
  function jitter(name){
    let h = 2166136261;
    for(let i=0;i<name.length;i++){ h ^= name.charCodeAt(i); h = Math.imul(h, 16777619); }
    return ((h >>> 0) % 1000) / 1000 - 0.5;
  }
  function esc(s){
    return String(s).replace(/[&<>"']/g, ch => ({
      '&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[ch]));
  }
  function metaOf(d){ return d.meta || {}; }
  function fmtRaw(v){
    if(v===null || v===undefined || v==='') return 'NA';
    if(typeof v === 'number'){
      if(!isFinite(v)) return 'NA';
      if(Number.isInteger(v)) return String(v);
      return (Math.abs(v)>=100 || (Math.abs(v)<0.01 && v!==0)) ? v.toExponential(1) : v.toFixed(2);
    }
    return esc(v);
  }
  function tipHtml(d, plot){
    const bits = [`<b>${esc(d.n)}</b>`, `run: ${esc(d.r)}`,
      `minor hits: <b>${d.h}</b> &nbsp;control hits: ${d.hc}`];
    if(plot.xName) bits.push(`${esc(plot.xName)}: ${fmtRaw(plot.xRaw(d))}`);
    if(plot.yName && plot.yName !== 'minor hits') bits.push(`${esc(plot.yName)}: ${fmtRaw(plot.yRaw(d))}`);
    Object.keys(metaOf(d)).filter(k => k !== 'run').forEach(k => {
      const v = metaOf(d)[k];
      bits.push(`${esc(k)}: ${fmtRaw(v)}`);
    });
    return bits.join('<br>');
  }

  function mountGroup(group){
    const plots = group.plots.map(p => {
      const cv = document.getElementById(p.id);
      return Object.assign({cv, tip: document.getElementById(p.id+'-tip'),
        legend: document.getElementById(p.id+'-leg'), ctx: cv.getContext('2d'),
        view: null, hover: -1, drag: null, margin: p.margin || 56}, p);
    });
    let hi = new Set(), pinned = new Set(), metric = 'n_hits', catOrder = [];
    let lo = 0, hiV = 1;
    const finite = a => a.filter(v => v !== null && isFinite(v));
    function ramp(t){
      t = Math.max(0, Math.min(1, t));
      const stops = [[237,242,247],[144,180,214],[99,133,186],[214,120,86],[193,39,45]];
      const x = t*(stops.length-1), i = Math.min(Math.floor(x), stops.length-2), f = x-i;
      const a = stops[i], b = stops[i+1];
      return `rgb(${Math.round(a[0]+f*(b[0]-a[0]))},${Math.round(a[1]+f*(b[1]-a[1]))},${Math.round(a[2]+f*(b[2]-a[2]))})`;
    }
    function catKey(d){
      const v = metaOf(d)[metric];
      if(v === null || v === undefined || v === '') return 'NA';
      return String(v);
    }
    function catFill(key){
      if(key === 'NA') return '#cbd5e0';
      const i = catOrder.indexOf(key);
      if(i < 0) return '#cbd5e0';
      if(i < PAL.length) return PAL[i];
      return `hsl(${Math.round((i * 137.508) % 360)},62%,42%)`;
    }
    function metricVal(d){
      if(SPEC[metric] === 'num'){
        const v = metaOf(d)[metric];
        return (v === null || v === undefined || v === '') ? null : +v;
      }
      if(SPEC[metric] === 'cat') return null;
      if(metric === 'tail_q') return (d.q === null || d.q <= 0) ? 320 : -Math.log10(d.q);
      if(metric === 'enrichment') return d.e === null ? null : Math.log2(Math.max(d.e, 0.05));
      return d[{n_hits:'h', MSRI:'m', stouffer:'st'}[metric]];
    }
    function tx(v, log){
      if(v === null || v === undefined || !isFinite(v)) return null;
      if(!log) return v;
      if(v <= 0) return null;
      return Math.log10(v);
    }
    function rescale(){
      if(SPEC[metric] === 'cat'){
        const seen = new Set();
        catOrder = [];
        D.forEach(d => { const k = catKey(d); if(!seen.has(k)){ seen.add(k); catOrder.push(k); } });
        catOrder.sort((a,b) => {
          if(a === 'NA') return 1;
          if(b === 'NA') return -1;
          return a.localeCompare(b, undefined, {numeric:true, sensitivity:'base'});
        });
        return;
      }
      catOrder = [];
      const vals = finite(D.map(metricVal)).sort((a,b)=>a-b);
      lo = vals.length ? vals[Math.floor(0.02*(vals.length-1))] : 0;
      hiV = vals.length ? vals[Math.floor(0.98*(vals.length-1))] : 1;
      if(hiV <= lo && vals.length) hiV = vals[vals.length-1];
      if(hiV <= lo) hiV = lo + 1;
    }
    function bounds(plot){
      const xs = [], ys = [];
      D.forEach((d,i) => {
        const x = tx(plot.x(d,i), plot.xlog), y = tx(plot.y(d,i), plot.ylog);
        if(x===null || y===null) return;
        xs.push(x); ys.push(y);
      });
      if(!xs.length) return {x0:-1,x1:1,y0:-1,y1:1};
      let x0=Math.min(...xs), x1=Math.max(...xs), y0=Math.min(...ys), y1=Math.max(...ys);
      const px=(x1-x0)*0.08+1e-9, py=(y1-y0)*0.08+1e-9;
      return {x0:x0-px, x1:x1+px, y0:y0-py, y1:y1+py};
    }
    function size(plot){
      const r = plot.cv.getBoundingClientRect(), dpr = window.devicePixelRatio||1;
      plot.cv.width = r.width*dpr; plot.cv.height = r.height*dpr;
      plot.ctx.setTransform(dpr,0,0,dpr,0,0);
      return r;
    }
    function frame(plot, r){
      const left = plot.margin, right = 16;
      const bottom = plot.axes ? 58 : 36;
      const top = plot.axes ? 16 : 20;
      return {left:left, right:right, top:top, bottom:bottom,
              w: Math.max(r.width-left-right, 1),
              h: Math.max(r.height-top-bottom, 1)};
    }
    function sx(plot, x, r){
      const f = frame(plot, r);
      return f.left + (x-plot.view.x0)/(plot.view.x1-plot.view.x0)*f.w;
    }
    function sy(plot, y, r){
      const f = frame(plot, r);
      return r.height-f.bottom - (y-plot.view.y0)/(plot.view.y1-plot.view.y0)*f.h;
    }
    function niceStep(span, target){
      if(!(span > 0)) return 1;
      const raw = span / target;
      const mag = Math.pow(10, Math.floor(Math.log10(raw)));
      const err = raw / mag;
      if(err >= 7.5) return 10*mag;
      if(err >= 3.5) return 5*mag;
      if(err >= 1.5) return 2*mag;
      return mag;
    }
    function niceTicks(lo, hi, target){
      const step = niceStep(hi-lo, target);
      const start = Math.ceil((lo - step*1e-9) / step) * step;
      const out = [];
      for(let v = start; v <= hi + step*1e-6 && out.length < 20; v += step){
        const rounded = Math.round(v / step) * step;
        out.push(Math.abs(rounded) < step*1e-8 ? 0 : rounded);
      }
      return out;
    }
    function fmtTick(v){
      const a = Math.abs(v);
      if(a === 0) return '0';
      if(a >= 1000 || a < 0.01) return v.toExponential(0);
      if(Math.abs(v - Math.round(v)) < 1e-6) return String(Math.round(v));
      return (Math.round(v*10)/10).toFixed(1);
    }
    function fmtCount(v){
      const a = Math.abs(v);
      if(!(a >= 0)) return '';
      if(Math.abs(v - Math.round(v)) < 1e-6 * Math.max(1, a)) return String(Math.round(v));
      if(a >= 100) return String(Math.round(v));
      if(a < 0.01) return v.toExponential(0);
      return (Math.round(v * 10) / 10).toFixed(1);
    }
    function axisTicks(lo, hi, target, log){
      if(!log) return niceTicks(lo, hi, target).map(v => ({pos:v, text:fmtTick(v), major:true}));
      const marks = (hi - lo) > 3 ? [1] : [1, 2, 5];
      const out = [];
      const e0 = Math.floor(lo) - 1, e1 = Math.ceil(hi) + 1;
      for(let e = e0; e <= e1 && out.length < 48; e++){
        marks.forEach(m => {
          const data = m * Math.pow(10, e);
          const pos = Math.log10(data);
          if(pos >= lo && pos <= hi) out.push({pos:pos, text:fmtCount(data), major:m === 1});
        });
      }
      if(out.filter(t => t.major).length >= 2) return out;
      return niceTicks(lo, hi, target).map(v => (
        {pos:v, text:fmtCount(Math.pow(10, v)), major:true}));
    }
    function drawTickRun(ticks, pixel, inside, drawMark, drawLabel, collapse){
      const visible = [];
      ticks.forEach(t => {
        const p = pixel(t.pos);
        if(!inside(p)) return;
        drawMark(p, !!t.major);
        visible.push({p:p, t:t});
      });
      if(!collapse){
        visible.forEach(v => drawLabel(v.p, v.t.text));
        return;
      }
      const labeled = [];
      const room = p => labeled.every(q => Math.abs(q - p) > 28);
      visible.filter(v => v.t.major).forEach(v => {
        if(room(v.p)){ drawLabel(v.p, v.t.text); labeled.push(v.p); }
      });
      visible.filter(v => !v.t.major).forEach(v => {
        if(room(v.p)){ drawLabel(v.p, v.t.text); labeled.push(v.p); }
      });
    }
    function drawAxes(plot, r, ctx){
      if(!plot.axes) return;
      const f = frame(plot, r);
      const yAxis = r.height - f.bottom;
      ctx.save();
      ctx.strokeStyle = '#718096';
      ctx.fillStyle = '#4a5568';
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(f.left, f.top);
      ctx.lineTo(f.left, yAxis);
      ctx.lineTo(f.left + f.w, yAxis);
      ctx.stroke();
      ctx.font = '11px system-ui,sans-serif';
      const xTicks = axisTicks(plot.view.x0, plot.view.x1, 6, !!plot.xlog);
      ctx.textAlign = 'center';
      ctx.textBaseline = 'top';
      drawTickRun(xTicks,
        v => sx(plot, v, r),
        X => X >= f.left + 8 && X <= f.left + f.w - 8,
        (X, major) => { ctx.beginPath(); ctx.moveTo(X, yAxis); ctx.lineTo(X, yAxis + (major ? 5 : 3)); ctx.stroke(); },
        (X, text) => { ctx.fillText(text, X, yAxis + 7); },
        !!plot.xlog);
      if(!plot.yTicks){
        const yTicksNum = axisTicks(plot.view.y0, plot.view.y1, 5, !!plot.ylog);
        ctx.textAlign = 'right';
        ctx.textBaseline = 'middle';
        drawTickRun(yTicksNum,
          v => sy(plot, v, r),
          Y => Y >= f.top + 8 && Y <= yAxis - 8,
          (Y, major) => { ctx.beginPath(); ctx.moveTo(f.left, Y); ctx.lineTo(f.left - (major ? 5 : 3), Y); ctx.stroke(); },
          (Y, text) => { ctx.fillText(text, f.left - 8, Y); },
          !!plot.ylog);
      }
      ctx.restore();
    }
    function legendHtml(){
      if(SPEC[metric] === 'cat'){
        return catOrder.map(k =>
          `<span class="sw" style="background:${catFill(k)}"></span>${esc(k)}`
        ).join('') + `&nbsp;&nbsp;<span style="color:#718096">${esc(metric)}</span>`;
      }
      const fmt = v => (Math.abs(v)>=100 || (Math.abs(v)<0.01 && v!==0)) ? v.toExponential(1) : v.toFixed(1);
      return `<span class="sw" style="background:${ramp(0)}"></span>${fmt(lo)}` +
        `<span class="sw" style="background:${ramp(.5)}"></span>` +
        `<span class="sw" style="background:${ramp(1)}"></span>${fmt(hiV)}` +
        `&nbsp;&nbsp;<span style="color:#718096">${esc(metric)}</span>`;
    }
    function fillOf(d){
      if(SPEC[metric] === 'cat') return catFill(catKey(d));
      const v = metricVal(d);
      return (v===null || !isFinite(v)) ? '#cbd5e0' : ramp((v-lo)/(hiV-lo));
    }
    function drawPlot(plot){
      if(!plot.view) plot.view = bounds(plot);
      const r = size(plot), ctx = plot.ctx;
      const f = frame(plot, r);
      ctx.clearRect(0,0,r.width,r.height);
      drawAxes(plot, r, ctx);
      ctx.fillStyle='#4a5568'; ctx.font='12px system-ui,sans-serif';
      ctx.textAlign='center'; ctx.textBaseline='alphabetic';
      ctx.fillText(plot.xlab||'', f.left + f.w/2, r.height-8);
      if(plot.ylab){
        ctx.save(); ctx.translate(14, (f.top + r.height - f.bottom)/2); ctx.rotate(-Math.PI/2);
        ctx.textAlign='center'; ctx.textBaseline='middle';
        ctx.fillText(plot.ylab, 0, 0); ctx.restore();
      }
      if(plot.yTicks){
        ctx.textAlign='right'; ctx.font='11px system-ui,sans-serif';
        plot.yTicks.forEach(t => {
          const label = t.text.length>24 ? t.text.slice(0,23)+'\u2026' : t.text;
          ctx.fillText(label, plot.margin-8, sy(plot, t.y, r)+4);
        });
      }
      if(plot.xyline){
        const lo = Math.max(plot.view.x0, plot.view.y0);
        const hi = Math.min(plot.view.x1, plot.view.y1);
        if(hi > lo){
          ctx.beginPath(); ctx.moveTo(sx(plot, lo, r), sy(plot, lo, r));
          ctx.lineTo(sx(plot, hi, r), sy(plot, hi, r));
          ctx.strokeStyle='#a0aec0'; ctx.lineWidth=1; ctx.setLineDash([4,3]); ctx.stroke();
          ctx.setLineDash([]);
        }
      }
      const labelled = [];
      for(let i=0;i<D.length;i++){
        const d = D[i];
        const Xv = tx(plot.x(d,i), plot.xlog), Yv = tx(plot.y(d,i), plot.ylog);
        if(Xv===null || Yv===null) continue;
        const X = sx(plot, Xv, r), Y = sy(plot, Yv, r);
        if(plot.axes){
          if(X<f.left || X>f.left+f.w || Y<f.top || Y>r.height-f.bottom) continue;
        } else if(X<plot.margin-8 || X>r.width || Y<0 || Y>r.height-20) continue;
        const on = hi.has(i) || pinned.has(i);
        ctx.beginPath();
        ctx.arc(X, Y, on?6.5:(plot.hover===i?5.5:3.4), 0, 6.2832);
        ctx.fillStyle = fillOf(d);
        ctx.globalAlpha = (hi.size && !on) ? 0.28 : 0.9;
        ctx.fill();
        if(on){ ctx.globalAlpha=1; ctx.lineWidth=1.8; ctx.strokeStyle='#111'; ctx.stroke();
                labelled.push([i,X,Y]); }
        ctx.globalAlpha=1;
      }
      ctx.font='bold 12px system-ui,sans-serif'; ctx.textAlign='left';
      labelled.slice(0,40).forEach(([i,X,Y]) => {
        const t = D[i].sh+'  ('+D[i].h+')', w = ctx.measureText(t).width;
        let lx = X+11, ly = Y-9;
        if(lx+w+10>r.width) lx = X-w-13;
        ctx.fillStyle='rgba(255,255,255,.93)';
        ctx.strokeStyle='#111'; ctx.lineWidth=.8;
        ctx.beginPath(); ctx.roundRect(lx-5, ly-13, w+10, 19, 5); ctx.fill(); ctx.stroke();
        ctx.fillStyle='#111'; ctx.fillText(t, lx, ly+1);
      });
      plot.legend.innerHTML = legendHtml();
    }
    function drawAll(){ plots.forEach(drawPlot); }
    function nearest(plot, mx, my){
      const r = plot.cv.getBoundingClientRect();
      let best=-1, bd=225;
      for(let i=0;i<D.length;i++){
        const Xv = tx(plot.x(D[i],i), plot.xlog), Yv = tx(plot.y(D[i],i), plot.ylog);
        if(Xv===null || Yv===null) continue;
        const dx=sx(plot, Xv, r)-mx, dy=sy(plot, Yv, r)-my, d2=dx*dx+dy*dy;
        if(d2<bd){ bd=d2; best=i; }
      }
      return best;
    }
    plots.forEach(plot => {
      plot.cv.addEventListener('mousemove', ev => {
        const r=plot.cv.getBoundingClientRect(), mx=ev.clientX-r.left, my=ev.clientY-r.top;
        if(plot.drag){
          const f = frame(plot, r);
          const kx=(plot.view.x1-plot.view.x0)/f.w;
          const ky=(plot.view.y1-plot.view.y0)/f.h;
          plot.view.x0 -= (mx-plot.drag.x)*kx; plot.view.x1 -= (mx-plot.drag.x)*kx;
          plot.view.y0 += (my-plot.drag.y)*ky; plot.view.y1 += (my-plot.drag.y)*ky;
          plot.drag={x:mx,y:my}; drawPlot(plot); return;
        }
        const i=nearest(plot, mx, my);
        if(i!==plot.hover){ plot.hover=i; drawPlot(plot); }
        if(i>=0){
          plot.tip.innerHTML = tipHtml(D[i], plot);
          plot.tip.style.display='block';
          plot.tip.style.left=Math.min(mx+16, r.width-250)+'px';
          plot.tip.style.top=(my+14)+'px';
        } else plot.tip.style.display='none';
      });
      plot.cv.addEventListener('mouseleave', () => {
        plot.tip.style.display='none'; plot.hover=-1; plot.drag=null; drawPlot(plot);
      });
      plot.cv.addEventListener('mousedown', ev => {
        const r=plot.cv.getBoundingClientRect(), mx=ev.clientX-r.left, my=ev.clientY-r.top;
        const i=nearest(plot, mx, my);
        if(i>=0){ pinned.has(i)?pinned.delete(i):pinned.add(i); drawAll(); }
        else plot.drag={x:mx,y:my};
      });
      window.addEventListener('mouseup', () => { plot.drag=null; });
      plot.cv.addEventListener('wheel', ev => {
        ev.preventDefault();
        const r=plot.cv.getBoundingClientRect(), mx=ev.clientX-r.left, my=ev.clientY-r.top;
        const fscale = ev.deltaY>0 ? 1.15 : 1/1.15;
        const box = frame(plot, r);
        const cx = plot.view.x0+(mx-box.left)/box.w*(plot.view.x1-plot.view.x0);
        const cy = plot.view.y0+(r.height-box.bottom-my)/box.h*(plot.view.y1-plot.view.y0);
        plot.view = {x0:cx+(plot.view.x0-cx)*fscale, x1:cx+(plot.view.x1-cx)*fscale,
                     y0:cy+(plot.view.y0-cy)*fscale, y1:cy+(plot.view.y1-cy)*fscale};
        drawPlot(plot);
      }, {passive:false});
    });
    const box = document.getElementById(group.prefix+'-search');
    function runSearch(){
      const q = box.value.trim().toLowerCase();
      hi = new Set();
      if(q){
        const terms = q.split(/[\s,;]+/).filter(Boolean);
        D.forEach((d,i)=>{
          const extra = Object.values(metaOf(d)).join(' ');
          const n=(d.n+' '+d.r+' '+extra).toLowerCase();
          if(terms.some(t=>n.includes(t))) hi.add(i);
        });
      }
      document.getElementById(group.prefix+'-nfound').textContent =
        q ? `${hi.size} match${hi.size===1?'':'es'}` : '';
      drawAll();
    }
    box.addEventListener('input', runSearch);
    document.getElementById(group.prefix+'-metric').addEventListener('change', ev => {
      metric = ev.target.value; rescale(); drawAll();
    });
    const pinBtn = document.getElementById(group.prefix+'-pin');
    if(pinBtn){
      pinBtn.addEventListener('click', () => {
        const raw = parseInt(document.getElementById(group.prefix+'-topx').value, 10);
        const n = Math.max(1, Math.min(40, isFinite(raw) ? raw : 10));
        pinned = new Set(D.map((d,i)=>[d.h,i]).sort((a,b)=>b[0]-a[0]).slice(0,n).map(p=>p[1]));
        drawAll();
      });
    }
    document.getElementById(group.prefix+'-clear').addEventListener('click', () => {
      pinned = new Set(); hi = new Set(); box.value='';
      document.getElementById(group.prefix+'-nfound').textContent='';
      plots.forEach(p => { p.view = bounds(p); });
      drawAll();
    });
    window.addEventListener('resize', drawAll);
    rescale();
    plots.forEach(p => { p.view = bounds(p); });
    drawAll();
  }

  const manyRuns = (CFG.runs||[]).length > 1;
  const yBand = manyRuns ? 0.55 : 0.8;
  function qcY(d){ return (d.ri||0) + jitter(d.n)*yBand; }
  const yTicks = manyRuns ? CFG.runs.map((name,i)=>({y:i, text:name})) : null;
  const qcMargin = manyRuns ? 150 : 48;
  mountGroup({
    prefix:'pca',
    plots:[{id:'pca', xlab:CFG.xlab, ylab:CFG.ylab, xlog:false, ylog:false,
            axes:true, margin:72,
            x:(d)=>d.x, y:(d)=>d.y, xRaw:(d)=>d.x, yRaw:(d)=>d.y,
            xName:'PC1', yName:'PC2'}]
  });
  mountGroup({
    prefix:'score',
    plots:[
      {id:'rank', xlab:'rank (log scale)', ylab:'minor-intron hits',
       xlog:true, ylog:false, axes:true, margin:84,
       x:(d,i)=>rank.get(i), y:(d)=>d.h, xRaw:(d,i)=>rank.get(i), yRaw:(d)=>d.h,
       xName:'rank', yName:'minor hits'},
      {id:'enrich', xlab:'hits on control introns', ylab:'hits on minor introns',
       xlog:true, ylog:true, xyline:true, axes:true, margin:84,
       x:(d)=>d.hc+0.5, y:(d)=>d.h+0.5, xRaw:(d)=>d.hc, yRaw:(d)=>d.h,
       xName:'control hits', yName:'minor hits'}
    ]
  });
  mountGroup({
    prefix:'qc',
    plots:[
      {id:'qcdepth', xlab:'median exonic depth', ylab: manyRuns?'run':'',
       axes:true,
       x:(d)=>d.depth, y:qcY, xRaw:(d)=>d.depth, yRaw:(d)=>CFG.runs[d.ri]||d.r,
       xName:'median exonic depth', yName:'run', yTicks:yTicks, margin:qcMargin},
      {id:'qcfrac', xlab:'fraction of targets measurable', ylab: manyRuns?'run':'',
       axes:true,
       x:(d)=>d.frac, y:qcY, xRaw:(d)=>d.frac, yRaw:(d)=>CFG.runs[d.ri]||d.r,
       xName:'fraction measurable', yName:'run', yTicks:yTicks, margin:qcMargin},
      {id:'qcbg', xlab:'intronic background (log2)', ylab: manyRuns?'run':'',
       axes:true,
       x:(d)=>d.bg, y:qcY, xRaw:(d)=>d.bg, yRaw:(d)=>CFG.runs[d.ri]||d.r,
       xName:'intronic background', yName:'run', yTicks:yTicks, margin:qcMargin}
    ]
  });
})();
"""


def _bar(prefix: str, options: str, placeholder: str, pin: bool = True) -> str:
    pin_html = ""
    if pin:
        pin_html = (
            f"<label class='pinlab'>pin top "
            f"<input id='{prefix}-topx' type='number' min='1' max='40' value='10' "
            f"aria-label='how many top samples to pin'></label>"
            f"<button type='button' id='{prefix}-pin'>pin</button>"
        )
    return (
        "<div class='bar'>"
        f"<input type='text' id='{prefix}-search' placeholder='{html.escape(placeholder)}'>"
        f"<span class='nfound' id='{prefix}-nfound'></span>"
        f"<select id='{prefix}-metric'>{options}</select>"
        f"{pin_html}"
        f"<button type='button' id='{prefix}-clear'>reset</button>"
        "</div>"
    )


def _panels(specs, tall: bool = False) -> str:
    cols = {1: "grid1", 2: "grid2", 3: "grid3"}.get(len(specs), "grid1")
    klass = "tall" if tall else ""
    bits = []
    for cid, cap, meaning in specs:
        bits.append(
            "<div class='panel'>"
            f"<div class='cap'>{cap}</div>"
            "<div class='wrap'>"
            f"<canvas id='{cid}' class='{klass}'></canvas>"
            f"<div class='tip' id='{cid}-tip'></div>"
            "</div>"
            f"<div class='legend' id='{cid}-leg'></div>"
            f"<div class='note'>{meaning}</div>"
            "</div>"
        )
    return f"<div class='{cols}'>{''.join(bits)}</div>"


def _img(path: str, meaning: str = "", scroll: bool = False) -> str:
    if not path:
        return ""
    tag = f"<img src='data:image/png;base64,{_b64(path)}'>"
    if scroll:
        tag = f"<div class='hscroll'>{tag}</div>"
    note = f"<div class='note'>{meaning}</div>" if meaning else ""
    return tag + note


def _html_toggle(items: list, prefix: str) -> str:
    """Show one HTML block at a time. Items are (id, button label, body html)."""
    if not items:
        return ""
    root = f"toggle-{prefix}"
    buttons = "".join(
        f"<button type='button' class='{'on' if n == 0 else ''}' "
        f"data-pane='{html.escape(i)}'>{html.escape(lab)}</button>"
        for n, (i, lab, _body) in enumerate(items))
    panes = []
    for n, (i, _lab, body) in enumerate(items):
        hidden = "" if n == 0 else " hidden"
        panes.append(
            f"<div class='tog-pane' data-pane='{html.escape(i)}'{hidden}>"
            f"{body}</div>")
    script = (
        "<script>(function(){"
        f"const root=document.getElementById('{root}');"
        "if(!root) return;"
        "const btns=root.querySelectorAll('button');"
        "btns.forEach(btn=>btn.addEventListener('click',()=>{"
        "btns.forEach(b=>b.classList.toggle('on', b===btn));"
        "root.querySelectorAll('.tog-pane').forEach(p=>{"
        "p.hidden=p.dataset.pane!==btn.dataset.pane;});}));})();</script>"
    )
    return (f"<div id='{root}'><div class='bar'>{buttons}</div>"
            + "".join(panes) + script + "</div>")


def _figure_toggle(items: list, prefix: str, scroll: bool = False) -> str:
    """Show one figure at a time. Items are (id, label, path, note[, extra html])."""
    shown = []
    for item in items:
        extra = item[4] if len(item) > 4 else ""
        if not item[2] and not extra:
            continue
        shown.append((item[0], item[1], item[2], item[3], extra))
    if not shown:
        return ""
    root = f"toggle-{prefix}"
    buttons = "".join(
        f"<button type='button' class='{'on' if n == 0 else ''}' "
        f"data-pane='{html.escape(i)}'>{html.escape(lab)}</button>"
        for n, (i, lab, _path, _note, _extra) in enumerate(shown))
    panes = []
    for n, (i, _lab, path, note, extra) in enumerate(shown):
        hidden = "" if n == 0 else " hidden"
        panes.append(
            f"<div class='tog-pane' data-pane='{html.escape(i)}'{hidden}>"
            f"{_img(path, note, scroll=scroll)}{extra}</div>")
    script = (
        "<script>(function(){"
        f"const root=document.getElementById('{root}');"
        "if(!root) return;"
        "const btns=root.querySelectorAll('button');"
        "btns.forEach(btn=>btn.addEventListener('click',()=>{"
        "btns.forEach(b=>b.classList.toggle('on', b===btn));"
        "root.querySelectorAll('.tog-pane').forEach(p=>{"
        "p.hidden=p.dataset.pane!==btn.dataset.pane;});}));})();</script>"
    )
    bar = f"<div class='bar'>{buttons}</div>" if len(shown) > 1 else ""
    return f"<div id='{root}'>{bar}" + "".join(panes) + script + "</div>"


def _heatmap_legend(label: str = "robust z above this sample's control introns"
                    ) -> str:
    """Color key for the retention heatmaps. Same scale on every ranking."""
    return (
        "<div class='hm-legend'>"
        "<div class='hm-bar' role='img' "
        f"aria-label='{html.escape(label)}, from −6 blue to +6 red'></div>"
        "<div class='hm-ticks'><span>−6 below</span><span>0</span>"
        "<span>+6 retained</span></div>"
        f"<div class='hm-cap'>{html.escape(label)}. Blue means this minor "
        "intron is below the cohort and below that sample's ordinary introns. "
        "Red means it is retained above both. The same scale is on each figure."
        "</div></div>")


def _heatmap_toggle(items: list) -> str:
    """Show one heatmap at a time. `items` are (id, button label, path, note)."""
    return _figure_toggle(items, "hm", scroll=True) + _heatmap_legend()


def interactive_upset_html(catalog: Sequence[dict]) -> str:
    """Picker for a custom UpSet plot. Samples are labeled in the menu."""
    if not catalog:
        return ""
    options = "".join(
        f"<option value='{html.escape(row['id'], quote=True)}'>"
        f"{html.escape(row['label'])}</option>"
        for row in catalog)
    diagnoses = []
    for row in catalog:
        diag = row.get("diagnosis") or ""
        if diag and diag not in diagnoses:
            diagnoses.append(diag)
    buttons = []
    if any(row.get("spliceopathy") for row in catalog):
        buttons.append(
            "<button type='button' data-pick='spliceopathy'>"
            "minor spliceopathy</button>")
    seen_diag = set()
    for diag in diagnoses:
        if _is_blank_label(diag):
            continue
        key = _cell_text(diag).casefold()
        if key in seen_diag:
            continue
        seen_diag.add(key)
        buttons.append(
            f"<button type='button' data-pick='diag:{html.escape(diag, quote=True)}'>"
            f"{html.escape(diag)}</button>")
    data = _js_json(catalog)
    return f"""
<div class='upset-tools'>
  <span class='pinlab'>Select a group, or pick names below.</span>
  {''.join(buttons)}
</div>
<input type='text' id='upset-filter' placeholder='filter samples' style='max-width:720px;width:100%;margin:4px 0;padding:6px 10px;font-size:13px;border:1px solid #cbd5e0;border-radius:5px'>
<select id='upset-pick' class='upset-pick' multiple size='10'>
{options}
</select>
<div class='upset-tools'>
  <button type='button' id='upset-draw'>Draw</button>
  <button type='button' id='upset-clear'>Clear</button>
</div>
<div id='upset-out' class='note'>Select two or more samples, then Draw. One sample lists that person's genes. Hold Ctrl (or Command) to select several names.</div>
<script type='application/json' id='upset-catalog'>{data}</script>
<script>
(function(){{
  const catalog = JSON.parse(document.getElementById('upset-catalog').textContent);
  const byId = new Map(catalog.map(d => [d.id, d]));
  const pick = document.getElementById('upset-pick');
  const out = document.getElementById('upset-out');
  const filter = document.getElementById('upset-filter');
  function selected(){{
    return Array.from(pick.selectedOptions).map(o => byId.get(o.value)).filter(Boolean);
  }}
  function selectIds(ids){{
    const want = new Set(ids);
    Array.from(pick.options).forEach(o => {{ o.selected = want.has(o.value); }});
    draw();
  }}
  document.querySelectorAll('[data-pick]').forEach(btn => {{
    btn.addEventListener('click', () => {{
      const key = btn.dataset.pick;
      const ids = catalog.filter(d => key === 'spliceopathy' ? d.spliceopathy
        : (d.diagnosis || '').trim().toLowerCase() === key.slice(5).trim().toLowerCase()).map(d => d.id);
      selectIds(ids);
    }});
  }});
  filter.addEventListener('input', () => {{
    const q = filter.value.trim().toLowerCase();
    Array.from(pick.options).forEach(o => {{
      o.hidden = q && !o.textContent.toLowerCase().includes(q);
    }});
  }});
  document.getElementById('upset-clear').addEventListener('click', () => {{
    Array.from(pick.options).forEach(o => {{ o.selected = false; }});
    out.className = 'note';
    out.innerHTML = 'Select two or more samples, then Draw.';
  }});
  document.getElementById('upset-draw').addEventListener('click', draw);
  pick.addEventListener('change', draw);

  function intersections(people){{
    const order = people.map(p => p.id);
    const buckets = new Map();
    people.forEach(p => p.genes.forEach(g => {{
      if(!buckets.has(g)) buckets.set(g, []);
      buckets.get(g).push(p.id);
    }}));
    const grouped = new Map();
    buckets.forEach((members, gene) => {{
      const key = order.filter(id => members.includes(id)).join('\\t');
      if(!grouped.has(key)) grouped.set(key, {{samples: order.filter(id => members.includes(id)), genes: []}});
      grouped.get(key).genes.push(gene);
    }});
    const rows = Array.from(grouped.values());
    rows.forEach(r => r.genes.sort());
    rows.sort((a,b) => b.genes.length - a.genes.length || b.samples.length - a.samples.length || a.samples.join().localeCompare(b.samples.join()));
    return rows;
  }}
  function esc(s){{ return String(s).replace(/[&<>]/g, c => ({{'&':'&amp;','<':'&lt;','>':'&gt;'}}[c])); }}
  function geneList(person){{
    const items = person.genes.map(g => '<li>'+esc(g)+'</li>').join('');
    const listing = items ? '<ul class="genelist">'+items+'</ul>' : '<div class="note">No genes with a minor-intron hit whose paired control is not a hit.</div>';
    return '<h3>'+esc(person.label)+'</h3><div class="sub">'+person.genes.length+' genes with a minor-intron hit and no hit on the paired control intron.</div>'+listing;
  }}
  function sharedTable(rows){{
    const shared = rows.filter(r => r.samples.length >= 2);
    if(!shared.length) return '<div class="note">No gene is shared by two or more of these samples.</div>';
    const shown = shared.slice(0, 200);
    const more = shared.length > shown.length ? '<div class="sub">Tallest '+shown.length+' of '+shared.length+' shared combinations.</div>' : '';
    const body = shown.map(r => '<tr><td>'+esc(r.samples.join(', '))+'</td><td>'+r.samples.length+'</td><td>'+r.genes.length+'</td><td>'+esc(r.genes.join(', '))+'</td></tr>').join('');
    return more + '<table><thead><tr><th>individuals</th><th>n_individuals</th><th>n_genes</th><th>genes</th></tr></thead><tbody>'+body+'</tbody></table>';
  }}
  function upsetSvg(people, rows){{
    const shown = rows.slice(0, 40);
    const nS = people.length, nI = Math.max(shown.length, 1);
    const rowH = 22, colW = 36, barW = 120, nameW = 108, top = 150, pad = 16;
    const left = barW + nameW;
    const w = left + nI * colW + pad, h = top + nS * rowH + pad;
    const sizes = people.map(p => p.genes.length);
    const maxSize = Math.max(1, ...sizes);
    const counts = shown.map(r => r.genes.length);
    const maxCount = Math.max(1, ...counts, 1);
    let svg = '<svg xmlns="http://www.w3.org/2000/svg" width="'+w+'" height="'+h+'" viewBox="0 0 '+w+' '+h+'">';
    svg += '<text x="'+(left + (nI*colW)/2)+'" y="16" text-anchor="middle" font-size="13" font-family="sans-serif">Genes with minor-specific retention ('+nS+' samples)</text>';
    people.forEach((p, i) => {{
      const y = top + i * rowH + rowH/2;
      const bw = (sizes[i] / maxSize) * (barW - 16);
      svg += '<rect x="'+(barW - 8 - bw)+'" y="'+(y-6)+'" width="'+bw+'" height="12" fill="#4a5568"/>';
      svg += '<text x="'+(left-8)+'" y="'+(y+4)+'" text-anchor="end" font-size="11" font-family="sans-serif">'+esc(p.id)+'</text>';
    }});
    svg += '<text x="8" y="'+(h-6)+'" font-size="10" font-family="sans-serif" fill="#4a5568">genes per person</text>';
    shown.forEach((r, x) => {{
      const cx = left + x * colW + colW/2;
      const shared = r.samples.length >= 2;
      const color = shared ? '#2b6cb0' : '#a0aec0';
      const bh = (r.genes.length / maxCount) * (top - 46);
      svg += '<rect x="'+(cx-10)+'" y="'+(top - 24 - bh)+'" width="20" height="'+bh+'" fill="'+color+'"/>';
      svg += '<text x="'+cx+'" y="'+(top - 28 - bh)+'" text-anchor="middle" font-size="10" font-family="sans-serif">'+r.genes.length+'</text>';
      const ys = r.samples.map(id => people.findIndex(p => p.id === id)).filter(i => i >= 0);
      if(ys.length >= 2){{
        const y1 = top + Math.min(...ys) * rowH + rowH/2;
        const y2 = top + Math.max(...ys) * rowH + rowH/2;
        svg += '<line x1="'+cx+'" y1="'+y1+'" x2="'+cx+'" y2="'+y2+'" stroke="'+color+'" stroke-width="2"/>';
      }}
      ys.forEach(i => {{
        const cy = top + i * rowH + rowH/2;
        svg += '<circle cx="'+cx+'" cy="'+cy+'" r="5" fill="'+color+'"/>';
      }});
    }});
    svg += '</svg>';
    const note = rows.length > shown.length ? '<div class="sub">Tallest '+shown.length+' of '+rows.length+' combinations. The table lists every shared combination.</div>' : '';
    return '<div class="upset-scroll">'+svg+'</div>'+note;
  }}
  function draw(){{
    const people = selected();
    out.className = '';
    if(!people.length){{
      out.className = 'note';
      out.textContent = 'Select two or more samples, then Draw.';
      return;
    }}
    if(people.length === 1){{
      out.innerHTML = geneList(people[0]);
      return;
    }}
    const rows = intersections(people);
    out.innerHTML = upsetSvg(people, rows) + sharedTable(rows);
  }}
}})();
</script>
"""


def write_html(path: str, ctx: dict) -> str:
    colours = ctx.get("colours") or _DEFAULT_COLOURS
    pts = _js_json(ctx["points"])
    spec = _js_json({item["id"]: item["kind"] for item in colours})
    cfg = _js_json({"xlab": ctx["xlab"], "ylab": ctx["ylab"],
                    "runs": ctx.get("run_names") or []})
    js = (_PLOTS_JS.replace("__DATA__", pts).replace("__SPEC__", spec)
          .replace("__CFG__", cfg))
    options = "".join(
        f"<option value='{html.escape(item['id'], quote=True)}'>"
        f"{html.escape(item['label'])}</option>"
        for item in colours)

    p: List[str] = [
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
        "<title>Intronaute report</title>",
        f"<style>{CSS}</style></head><body>",
        "<h1>Intronaute &mdash; minor-intron retention</h1>",
        f"<div class='sub'>{ctx['n_samples']} samples &middot; "
        f"{ctx['n_minor']} minor introns &middot; {ctx['n_control']} control "
        f"introns &middot; hit mode <b>{ctx.get('hit_mode', 'coverage')}</b> "
        f"&middot; {datetime.datetime.now():%Y-%m-%d %H:%M}</div>",
        "__TOC__",
    ]
    toc: List[tuple] = []
    p += [
        _section(toc, "sec-pca", "1. Sample map (PCA of minor introns)"),
        f"<div class='note'>{_pca_position(ctx)}</div>",
        "<div class='sub'><b>Type in the box</b> to highlight a sample, "
        "scroll to zoom, drag to pan, click a point to pin its label. "
        "Colour follows the selected score or a samplesheet column. "
        f"{_defs('n_hits', 'tail_q')} Colour by <b>n_hits</b> and Pin top N "
        "use the hit count from section 4. The place of each point stays the "
        "continuous score above. "
        "Colour by −log10(<b>tail_q</b>) stretches a small false-discovery rate "
        "so the strongest samples are the high end of the ramp. "
        "Samplesheet columns are labels only: they are not used in the scores.</div>",
        _bar("pca", options, "search sample, run, or annotation"),
        _panels([("pca", "PC1 vs PC2", _pca_caption(ctx))], tall=True),
    ]

    rankings = ctx.get("rankings")
    if not rankings and ctx.get("top_samples") is not None:
        rankings = [("n_hits", "n_hits", _table(ctx["top_samples"]))]
    rank_notes = {
        "n_hits": (
            f"{_defs('n_hits')} This table is the samples with the largest "
            "<b>n_hits</b>. Ties break on Stouffer, high first."),
        "n_hits_control": (
            f"{_defs('n_hits_control')} This table is the samples with the "
            "largest <b>n_hits_control</b>. Ties break on Stouffer, high first."),
        "tail_q": (
            f"{_defs('tail_q')} This table is the samples with the smallest "
            "<b>tail_q</b>. Missing values stay at the bottom. Ties break on "
            "Stouffer, high first."),
        "enrichment": (
            f"{_defs('enrichment')} This table is the samples with the largest "
            "<b>enrichment</b>. Missing values stay at the bottom. Ties break "
            "on Stouffer, high first."),
    }
    rank_panes = []
    for key, label, body in rankings or []:
        note = rank_notes.get(key, "")
        table = body if isinstance(body, str) else _table(body)
        rank_panes.append((key, label, f"<div class='note'>{note}</div>{table}"))
    p += [_section(toc, "sec-rank", "2. Highest-ranking samples"),
          "<div class='sub'>Switch the table among four rankings. "
          f"{_defs('n_hits', 'n_hits_control', 'tail_q', 'enrichment')} "
          "The two plots are the same panels as "
          "<code>score_distribution.png</code>. Search highlights a sample on "
          "both. Click a dot to pin its name, and click it again to unpin. "
          "Colour and pin top N apply to both. Pin top N uses <b>n_hits</b>.</div>",
          _html_toggle(rank_panes, "rank"),
          _bar("score", options, "search a sample"),
          _panels([
              ("rank", "Rank by n_hits",
               _defs("n_hits") + " Each point is one sample. Rank 1, the "
               "leftmost point, has the largest <b>n_hits</b>. The curve "
               "falls as that count drops. Search highlights a sample. "
               "Click a dot to pin its name."),
              ("enrich", "n_hits vs n_hits_control",
               _defs("n_hits", "n_hits_control") + " Each point is one sample. "
               "The dashed line is equal minor and control hits. Points above "
               "the line have <b>n_hits</b> above <b>n_hits_control</b>. "
               "Search highlights a sample. Click a dot to pin its name."),
          ])]

    heatmaps = [
        ("hits", "rank by n_hits", ctx.get("fig_heatmap"),
         _defs("n_hits") + " Each cell is one minor intron in one sample. "
         "Red means that intron is retained above the rest of the cohort and "
         "above that sample's own control introns. Blue means it is below. "
         "Rows are the samples with the largest <b>n_hits</b>. Every "
         "measurable minor intron is a column, labeled by its gene. A gene "
         "with several minor introns has one column for each. Columns are "
         "ordered by the largest adjusted z in these rows, strongest on the "
         "left. Slide the figure sideways to see the rest."),
        ("control", "rank by n_hits_control", ctx.get("fig_heatmap_control"),
         _defs("n_hits_control") + " Rows are the samples with the largest "
         "<b>n_hits_control</b>. The columns are again every measurable minor "
         "intron, reordered by the largest adjusted z in these rows."),
        ("q", "rank by tail_q", ctx.get("fig_heatmap_q"),
         _defs("tail_q") + " Rows are the samples with the smallest "
         "<b>tail_q</b>. A sample can appear here with a smaller "
         "<b>n_hits</b> than on the first heatmap when those hits are "
         "specific to minor introns. The columns are again every measurable "
         "minor intron, reordered by the largest adjusted z in these rows."),
    ]
    if any(item[2] for item in heatmaps):
        p += [_section(toc, "sec-heatmap", "3. Most retained introns"),
              "<div class='sub'>Each cell is one minor intron in one sample. "
              "The color is that intron's robust z minus the median robust z "
              "of the same sample's U2 control introns. When the non-minor "
              "introns in that sample are retained too, the cell moves toward "
              "white by that amount. Red means this minor intron is retained "
              "above the cohort and above that sample's own ordinary introns. "
              "Blue means it is below them. The scale is −6 to +6. Every "
              "measurable minor intron is shown, labeled by its gene. Slide "
              "the figure sideways. The strongest value in these rows is on "
              "the left. Switching the ranking changes the samples and can "
              "change that order. It does not change the color rule. Hit "
              "calls are still the unadjusted z.</div>",
              _heatmap_toggle(heatmaps)]

    if ctx.get("upset"):
        has_plot = any(item[2] for item in ctx["upset"])
        if has_plot:
            shared_title = "4. Shared minor-intron genes"
            shared_sub = (
                "Each person is the set of genes in which a minor intron is a "
                "hit and the control introns paired to that minor intron are "
                "not. A gene is left out when retention is also called on the "
                "matched U2 intron. A blue bar is a combination of people who "
                "share those genes. A grey bar is genes found in only one "
                "person. The number on a bar is how many genes are in that "
                "combination. Bars run from the most genes on the left to the "
                "fewest on the right. The table under the plot names the shared "
                "genes, in that same order. "
                "A group with one person is a list of that person's genes "
                "instead of an UpSet plot. The list is labeled with the "
                "diagnosis when the sample has one, and with minor spliceopathy "
                "otherwise. A diagnosis of control, N/A, NA, null, none, 0, "
                "or a blank cell does not get its own group.")
        else:
            shared_title = "4. Minor-intron genes"
            shared_sub = (
                "Each group here has one person, so there is no UpSet plot. "
                "The list is every gene in which a minor intron is a hit and "
                "the control introns paired to that minor intron are not. A "
                "gene is left out when retention is also called on the matched "
                "U2 intron. The label is the diagnosis when that sample has "
                "one, and minor spliceopathy otherwise. A diagnosis of control, "
                "N/A, NA, null, none, 0, or a blank cell does not get its "
                "own group.")
        p += [_section(toc, "sec-shared", shared_title),
              f"<div class='note'>{_hit_definition(ctx)}</div>",
              f"<div class='sub'>{shared_sub}</div>",
              _figure_toggle(ctx["upset"], "up")]

    if ctx.get("upset_catalog"):
        p += [_section(toc, "sec-upset-custom",
                       "5. Build an UpSet plot"),
              "<div class='sub'>Choose the samples to compare. Each name is "
              "labeled with diagnosis and spliceopathy status. "
              "<b>minor spliceopathy</b> and each diagnosis button fill the "
              "selection. A diagnosis button is shown for every label other "
              "than control, N/A, NA, null, none, 0, or blank. Two or more "
              "samples draw an UpSet plot and list "
              "every gene combination those people share. Bars run from the "
              "most genes on the left to the fewest on the right. One sample lists "
              "that person's genes with a minor-intron hit and no hit on the "
              "paired control intron. A gene is left out when retention is "
              "also called on the matched U2 intron.</div>",
              interactive_upset_html(ctx["upset_catalog"])]

    if len(ctx.get("top_events", [])):
        p += [_section(toc, "sec-events", "6. Top events"),
              f"<div class='note'>{_hit_definition(ctx)}</div>",
              f"<div class='sub'>{_top_events_scoring(ctx)}</div>",
              _table(ctx["top_events"])]

    p += [_section(toc, "sec-qc", "7. Quality control"),
          "<div class='sub'>A sample with low exonic depth, or with an intronic "
          "background far above the rest, should be read with caution. If runs "
          "separate here, that is a batch effect, and "
          f"{_defs('tail_q')} These are the three panels of <code>qc.png</code>. "
          "Colour and search work as on the other plots. "
          f"{_defs('n_hits')} Pin top N marks the samples with the largest "
          "<b>n_hits</b>, so a high-ranking sample can be found on the "
          "depth axis.</div>",
          _bar("qc", options, "search a sample"),
          _panels([
              ("qcdepth", "Median exonic depth",
               "Each point is one sample. This is the typical read depth "
               "on the exons used as the denominator of the retention "
               "score. A low value means the library may be too shallow "
               "for the scores to be trusted."),
              ("qcfrac", "Fraction of targets measurable",
               "Each point is one sample. This is the share of targets "
               "whose flanking exons were deep enough to score. A low "
               "value means much of that sample could not be measured."),
              ("qcbg", "Intronic background (log2, control introns)",
               "Each point is one sample. This is the control-intron "
               "offset subtracted before scoring. A high value means "
               "generic intronic background, such as degradation or DNA "
               "contamination, rather than U12-specific retention."),
          ]),
          _img(ctx.get("fig_qc", ""),
               "Left to right: median exonic depth, the fraction of targets "
               "deep enough to score, and the control-intron background. "
               "Low depth or a low measurable fraction means the sample is "
               "hard to score. A high background means generic intronic "
               "noise. When samples are split by run, the bar is that "
               "run's median.")]
    if len(ctx.get("run_table", [])):
        p.append(
            "<div class='note'>"
            f"{_defs('n_hits', 'n_hits_control')} "
            "One row is one <b>run</b>. <b>n</b> is how many samples are in "
            "that run. <b>n_hits_median</b> and <b>n_hits_max</b> are the "
            "median and the maximum <b>n_hits</b> in the run. "
            "<b>control_hits_median</b> is the median <b>n_hits_control</b>."
            "</div>" + _table(ctx["run_table"]))

    p += [_section(toc, "sec-params", "8. Parameters"),
          "<div class='sub'>Every key the configuration file can set. "
          "Analyse rows are the values used for this report. "
          "<b>min_mapq</b>, <b>min_overhang</b>, <b>keep_duplicates</b>, "
          "<b>library</b>, and <b>engine</b> are the values recorded with "
          "the cached counts. When those files disagree, the row lists each "
          "value and how many samples used it. The other prepare and "
          "quantify rows are taken from the config file when that key is "
          "set, and from the default when it is not.</div>",
          _table(ctx["params"]),
          f"<script>{js}</script>", "</body></html>"]
    p[p.index("__TOC__")] = _contents(toc)

    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(p))
    return path
