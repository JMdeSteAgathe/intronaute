"""On-disk layout written next to report.html.

The HTML report embeds its figures, so these paths are for the files on disk.
Colour toggles match the report menus: each score and each samplesheet column
is one PNG.
"""
from __future__ import annotations

import colorsys
import os
import re
import shutil
from typing import Dict, Sequence

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import ListedColormap, Normalize  # noqa: E402

PCA_DIR = "sample_mapped_PCA_minor_introns"
PCA_GENE_DIR = "sample_mapped_PCA_same_gene_control"
RANKED_DIR = "samples_ranked_minor_introns"
HEAT_DIR = "most_retained_minor_introns"
SHARED_DIR = "shared_minor_intron_genes"
QC_DIR = "quality_control"

# diagnosis is the short label used in the file names (dx_intronic_background).
_SLUG = {"diagnosis": "dx"}

_PAL = [
    "#c1272d", "#2b6cb0", "#2f855a", "#b7791f", "#6b46c1", "#c05621",
    "#2c7a7b", "#97266d", "#2d3748", "#dd6b20", "#3182ce", "#38a169",
    "#d53f8c", "#805ad5", "#4a5568", "#9b2c2c", "#2a4365", "#276749",
    "#744210", "#44337a",
]
_NA = "#cbd5e0"
_STOPS = np.array(
    [[237, 242, 247], [144, 180, 214], [99, 133, 186], [214, 120, 86],
     [193, 39, 45]], dtype=float) / 255.0

_LEGACY_ROOT = (
    "pca_coordinates.tsv", "pca_variance.tsv", "pca_loadings.tsv",
    "sample_scores.tsv", "hits.tsv.gz", "hits.rds",
    "upset_shared_genes.tsv", "singleton_minor_hit_genes.tsv",
)

_UPSET_FILE = {
    "spliceopathy": "upset_minor_spliceopathy.png",
    "RNU4ATAC": "upset_RNU4ATAC.png",
    "RNU12": "upset_RNU12.png",
    "RNU6ATAC": "upset_RNU6ATAC.png",
}


def _clear_pngs(directory: str) -> None:
    for name in os.listdir(directory):
        if name.lower().endswith(".png"):
            os.remove(os.path.join(directory, name))


def colour_slug(metric_id: str) -> str:
    if metric_id in _SLUG:
        return _SLUG[metric_id]
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", str(metric_id)).strip("_")
    return slug or "colour"


def upset_filename(group_id: str) -> str:
    return _UPSET_FILE.get(group_id, f"upset_{colour_slug(group_id)}.png")


def prepare(outdir: str) -> Dict[str, str]:
    """Create the result folders and drop PNGs left by an older run."""
    ranked = os.path.join(outdir, RANKED_DIR)
    qc = os.path.join(outdir, QC_DIR)
    dirs = {
        "pca": os.path.join(outdir, PCA_DIR),
        "ranked": ranked,
        "ranked_hits": os.path.join(ranked, "n_hits"),
        "ranked_vs": os.path.join(ranked, "n_hits_vs_n_hits_control"),
        "heat": os.path.join(outdir, HEAT_DIR),
        "shared": os.path.join(outdir, SHARED_DIR),
        "qc": qc,
        "qc_depth": os.path.join(qc, "median_exonic_depth"),
        "qc_frac": os.path.join(qc, "frac_targets_meas"),
        "qc_bg": os.path.join(qc, "intronic_background"),
    }
    for path in dirs.values():
        os.makedirs(path, exist_ok=True)
        _clear_pngs(path)
    _remove_figures_dir(outdir)
    return dirs


def copy_png(src: str, dest: str) -> str:
    if not src:
        return ""
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    shutil.copy2(src, dest)
    return dest


def _remove_figures_dir(outdir: str) -> None:
    path = os.path.join(outdir, "figures")
    if os.path.isdir(path):
        shutil.rmtree(path)


def remove_legacy(outdir: str) -> None:
    """Drop tables and the old figures folder that used to sit beside report.html."""
    for name in _LEGACY_ROOT:
        path = os.path.join(outdir, name)
        if os.path.isfile(path):
            os.remove(path)
    _remove_figures_dir(outdir)
    gene_dir = os.path.join(outdir, PCA_GENE_DIR)
    if os.path.isdir(gene_dir):
        shutil.rmtree(gene_dir)


def write_colour_plots(scores: pd.DataFrame, pca: pd.DataFrame, qc: pd.DataFrame,
                       colours: Sequence[dict], xlab: str, ylab: str,
                       dirs: Dict[str, str]) -> None:
    """One PNG per report colour toggle, for the PCA, the two score plots, and QC."""
    view = _view_table(scores, pca, qc)
    runs = list(dict.fromkeys(view["run"].astype(str)))
    run_index = {name: i for i, name in enumerate(runs)}
    jitter = _jitter(view["sample"].astype(str).tolist())
    y_run = np.array([run_index[r] for r in view["run"].astype(str)], dtype=float)
    y_run = y_run + jitter * (0.55 if len(runs) > 1 else 0.8)

    for item in colours:
        slug = colour_slug(item["id"])
        label = str(item.get("label", item["id"]))
        label = label.split(": ", 1)[-1]
        colors, legend = _colours_for(view, item)

        _scatter(
            view["PC1"].to_numpy(float), view["PC2"].to_numpy(float), colors,
            os.path.join(dirs["pca"], f"pca_{slug}.png"),
            title=f"PCA of minor introns, coloured by {label}",
            xlabel=xlab, ylabel=ylab, legend=legend)
        rank = np.arange(1, len(view) + 1)
        _scatter(
            rank.astype(float), view["n_hits"].to_numpy(float), colors,
            os.path.join(dirs["ranked_hits"], f"rank_{slug}.png"),
            title=f"Minor-intron hits by rank, coloured by {label}",
            xlabel="rank (log scale)", ylabel="n_hits",
            xlog=True, legend=legend)
        _scatter(
            view["n_hits_control"].to_numpy(float) + 0.5,
            view["n_hits"].to_numpy(float) + 0.5, colors,
            os.path.join(dirs["ranked_vs"], f"scatter_{slug}.png"),
            title=f"Minor vs control hits, coloured by {label}",
            xlabel="hits on control introns",
            ylabel="hits on minor introns",
            xlog=True, ylog=True, diagonal=True, legend=legend)

        qc_specs = (
            ("qc_depth", "median_exon_depth", "median_exonic_depth",
             "Median exonic depth"),
            ("qc_frac", "frac_measurable", "frac_targets_meas",
             "Fraction of targets measurable"),
            ("qc_bg", "control_background", "intronic_background",
             "Intronic background"),
        )
        for key, column, stem, nice in qc_specs:
            _scatter(
                pd.to_numeric(view[column], errors="coerce").to_numpy(float),
                y_run, colors,
                os.path.join(dirs[key], f"qc_{slug}_{stem}.png"),
                title=f"{nice}, coloured by {label}",
                xlabel=nice, ylabel="run" if len(runs) > 1 else "",
                yticks=(list(range(len(runs))), runs) if len(runs) > 1 else None,
                legend=legend)


def _view_table(scores: pd.DataFrame, pca: pd.DataFrame,
                qc: pd.DataFrame) -> pd.DataFrame:
    base = scores.copy()
    if "sample" not in base.columns:
        raise ValueError("scores have no sample column")
    base = base.drop_duplicates("sample").set_index("sample")
    base.index = base.index.astype(str)
    pcs = pca[["PC1", "PC2"]].copy()
    pcs.index = pcs.index.astype(str)
    q = qc.drop_duplicates("sample").set_index("sample")
    q.index = q.index.astype(str)
    view = base.join(pcs, how="left").join(
        q[["run", "median_exon_depth", "frac_measurable", "control_background"]],
        how="left", rsuffix="_qc")
    if "run_qc" in view.columns:
        view["run"] = view["run"].where(view["run"].notna(), view["run_qc"])
        view = view.drop(columns=["run_qc"])
    view = view.reset_index()
    # Rank plots follow n_hits, high first, matching the report.
    view = view.sort_values(["n_hits", "stouffer"], ascending=False)
    view["run"] = view["run"].fillna("run").astype(str)
    return view.reset_index(drop=True)


def _jitter(names: Sequence[str]) -> np.ndarray:
    """Stable vertical jitter in [-0.5, 0.5], one value per sample name."""
    out = np.empty(len(names), dtype=float)
    for i, name in enumerate(names):
        h = 2166136261
        for ch in name:
            h ^= ord(ch)
            h = (h * 16777619) & 0xFFFFFFFF
        out[i] = (h % 1000) / 1000.0 - 0.5
    return out


def _metric_values(view: pd.DataFrame, item: dict) -> np.ndarray:
    metric = item["id"]
    kind = item.get("kind", "num")
    if kind == "score":
        if metric == "tail_q":
            q = pd.to_numeric(view["tail_q"], errors="coerce").to_numpy(float)
            out = np.full(q.shape, np.nan)
            positive = np.isfinite(q) & (q > 0)
            out[positive] = -np.log10(q[positive])
            out[np.isfinite(q) & (q <= 0)] = 320.0
            return out
        if metric == "enrichment":
            e = pd.to_numeric(view["enrichment"], errors="coerce").to_numpy(float)
            return np.where(np.isfinite(e), np.log2(np.maximum(e, 0.05)), np.nan)
        column = {"n_hits": "n_hits", "MSRI": "MSRI", "stouffer": "stouffer"}.get(
            metric, metric)
        return pd.to_numeric(view[column], errors="coerce").to_numpy(float)
    if kind == "cat":
        text = view[metric].map(lambda v: "" if pd.isna(v) else str(v).strip())
        return text.to_numpy()
    return pd.to_numeric(view[metric], errors="coerce").to_numpy(float)


def _colours_for(view: pd.DataFrame, item: dict):
    values = _metric_values(view, item)
    if item.get("kind") == "cat":
        return _categorical_colors(values)
    return _continuous_colors(values, str(item.get("label", item["id"])).split(": ", 1)[-1])


def _ramp(t: np.ndarray) -> np.ndarray:
    t = np.clip(np.asarray(t, float), 0.0, 1.0)
    x = t * (len(_STOPS) - 1)
    i = np.minimum(np.floor(x).astype(int), len(_STOPS) - 2)
    f = (x - i)[:, None]
    return _STOPS[i] * (1.0 - f) + _STOPS[i + 1] * f


def _continuous_colors(values: np.ndarray, label: str):
    v = np.asarray(values, float)
    finite = v[np.isfinite(v)]
    if finite.size == 0:
        lo, hi = 0.0, 1.0
    else:
        lo = float(np.quantile(finite, 0.02))
        hi = float(np.quantile(finite, 0.98))
        if hi <= lo:
            hi = float(np.max(finite))
        if hi <= lo:
            hi = lo + 1.0
    t = (v - lo) / (hi - lo)
    colors = _ramp(np.where(np.isfinite(v), t, 0.0))
    colors[~np.isfinite(v)] = np.array([203, 213, 224]) / 255.0
    return colors, ("continuous", label, lo, hi)


_NULLISH = frozenset({"", "nan", "null", "none", "na", "n/a", "<na>"})


def _choose_category_spelling(spellings: list) -> str:
    """Most common spelling, then the one with more capitals."""
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


def _category_key(value):
    """Case-folded category, with null-like cells collapsed to NA."""
    if value is None:
        return "NA"
    try:
        if pd.isna(value):
            return "NA"
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    folded = text.casefold()
    if folded in _NULLISH:
        return "NA"
    return folded


def _categorical_colors(values: np.ndarray):
    spellings: Dict[str, list] = {}
    keys = []
    for value in values:
        key = _category_key(value)
        keys.append(key)
        if key == "NA":
            continue
        spellings.setdefault(key, []).append(str(value).strip())
    display = {key: _choose_category_spelling(seen)
               for key, seen in spellings.items()}
    shown = ["NA" if key == "NA" else display[key] for key in keys]
    keys = shown
    order = sorted({k for k in keys if k != "NA"})
    if any(k == "NA" for k in keys):
        order.append("NA")
    palette = {}
    for i, key in enumerate(order):
        if key == "NA":
            palette[key] = _hex_rgb(_NA)
        elif i < len(_PAL):
            palette[key] = _hex_rgb(_PAL[i])
        else:
            hue = (i * 137.508) % 360
            palette[key] = colorsys.hls_to_rgb(hue / 360.0, 0.42, 0.62)
    colors = np.array([palette[k] for k in keys], dtype=float)
    return colors, ("categorical", order, palette)


def _hex_rgb(color: str):
    text = color.lstrip("#")
    return tuple(int(text[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


def _scatter(x, y, colors, path: str, title: str, xlabel: str, ylabel: str,
             legend, xlog: bool = False, ylog: bool = False,
             diagonal: bool = False, yticks=None) -> str:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    ax.scatter(x[ok], y[ok], c=np.asarray(colors)[ok], s=22, alpha=0.85,
               edgecolor="none", zorder=3)
    if diagonal and ok.any():
        lo = float(np.nanmin([x[ok].min(), y[ok].min()]))
        hi = float(np.nanmax([x[ok].max(), y[ok].max()]))
        ax.plot([lo, hi], [lo, hi], lw=0.9, color="#a0aec0", ls="--", zorder=1)
    if xlog:
        ax.set_xscale("log")
    if ylog:
        ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if yticks is not None:
        ax.set_yticks(yticks[0])
        ax.set_yticklabels(yticks[1], fontsize=8)
    ax.set_title(title, fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    _add_legend(fig, ax, legend)
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    return path


def _add_legend(fig, ax, legend) -> None:
    kind = legend[0]
    if kind == "continuous":
        label, lo, hi = legend[1], legend[2], legend[3]
        gradient = _ramp(np.linspace(0, 1, 256))
        cmap = ListedColormap(gradient)
        scale = plt.cm.ScalarMappable(
            norm=Normalize(vmin=lo, vmax=hi), cmap=cmap)
        scale.set_array([])
        bar = fig.colorbar(scale, ax=ax, fraction=0.046, pad=0.04)
        bar.set_label(label, fontsize=8)
        return
    _kind, order, palette = legend
    handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=palette[k],
                   markeredgecolor="none", markersize=7, label=k)
        for k in order
    ]
    if not handles:
        return
    ax.legend(handles=handles, frameon=False, fontsize=7, loc="best",
              title=None)
