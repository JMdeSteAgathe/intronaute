"""Figures and a self-contained HTML report.

The first thing in the report is an interactive PCA: a canvas scatter with a
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


def _save(fig, path: str) -> str:
    fig.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
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
                vmax: float = 6.0) -> str:
    order = scores.sort_values("n_hits", ascending=False)["sample"].head(n_samples)
    pos = {s: i for i, s in enumerate(samples)}
    rows = [pos[s] for s in order if s in pos]
    if not rows:
        return ""
    cols = np.where(is_minor)[0]
    sub = z[np.ix_(rows, cols)]
    rank = np.nanmax(sub, axis=0)
    top = cols[np.argsort(-np.nan_to_num(rank, nan=-9))[:n_introns]]
    M = z[np.ix_(rows, top)]
    labels = [str(meta.loc[target_ids[j], "gene_name"]) for j in top]

    fig, ax = plt.subplots(figsize=(max(7.0, .24 * len(top)),
                                    max(3.4, .28 * len(rows))))
    im = ax.imshow(M, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(np.arange(len(top)))
    ax.set_xticklabels(labels, rotation=90, fontsize=6.5)
    ax.set_yticks(np.arange(len(rows)))
    ax.set_yticklabels([samples[i] for i in rows], fontsize=7)
    ax.set_title(f"Top {len(rows)} samples x top {len(top)} minor introns "
                 f"(robust z)", fontsize=10)
    fig.colorbar(im, ax=ax, shrink=.6, label="robust z")
    return _save(fig, os.path.join(outdir, "heatmap.png"))


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


# --------------------------------------------------------------------------- #
# interactive PCA
# --------------------------------------------------------------------------- #
_PCA_JS = r"""
(function(){
  const D = __DATA__, VAR = __VAR__;
  const cv = document.getElementById('pca'), tip = document.getElementById('tip');
  const ctx = cv.getContext('2d');
  let view = null, hi = new Set(), pinned = new Set(), metric = 'n_hits';
  let hover = -1, drag = null;

  const finite = a => a.filter(v => v !== null && isFinite(v));
  function ramp(t){
    t = Math.max(0, Math.min(1, t));
    const stops = [[237,242,247],[144,180,214],[ 99,133,186],[214,120, 86],[193, 39, 45]];
    const x = t*(stops.length-1), i = Math.min(Math.floor(x), stops.length-2), f = x-i;
    const a = stops[i], b = stops[i+1];
    return `rgb(${Math.round(a[0]+f*(b[0]-a[0]))},${Math.round(a[1]+f*(b[1]-a[1]))},${Math.round(a[2]+f*(b[2]-a[2]))})`;
  }
  function metricVal(d){
    if(metric === 'tail_q') return (d.q === null || d.q <= 0) ? 320 : -Math.log10(d.q);
    if(metric === 'enrichment') return d.e === null ? null : Math.log2(Math.max(d.e, 0.05));
    return d[{n_hits:'h', MSRI:'m', stouffer:'st'}[metric]];
  }
  let lo = 0, hiV = 1;
  function rescale(){
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
      ctx.fillStyle = (v===null||!isFinite(v)) ? '#cbd5e0' : ramp((v-lo)/(hiV-lo));
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
    document.getElementById('legend').innerHTML =
      `<span class="sw" style="background:${ramp(0)}"></span>${VAR.fmt(lo)}` +
      `<span class="sw" style="background:${ramp(.5)}"></span>` +
      `<span class="sw" style="background:${ramp(1)}"></span>${VAR.fmt(hiV)}` +
      `&nbsp;&nbsp;<span style="color:#718096">${metric}</span>`;
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
      tip.innerHTML = `<b>${d.n}</b><br>run: ${d.r}<br>minor hits: <b>${d.h}</b>`+
        ` &nbsp;control hits: ${d.hc}<br>enrichment: ${d.e===null?'NA':d.e.toFixed(1)}`+
        `&times; &nbsp;q: ${d.q===null?'NA':d.q.toExponential(1)}<br>`+
        `MSRI ${d.m===null?'NA':d.m.toFixed(2)} &nbsp;PC1 ${d.x.toFixed(1)} &nbsp;PC2 ${d.y.toFixed(1)}`;
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
      D.forEach((d,i)=>{ const n=(d.n+' '+d.r).toLowerCase();
                         if(terms.some(t=>n.includes(t))) hi.add(i); });
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
h2{font-size:16px;margin-top:34px;border-bottom:1px solid #e2e8f0;padding-bottom:5px}
.sub{color:#718096;font-size:12.5px;margin-bottom:16px}
table{border-collapse:collapse;font-size:11.5px;margin:10px 0;width:100%}
th,td{border:1px solid #e2e8f0;padding:4px 7px;text-align:right}
th{background:#f7fafc;font-weight:600;text-align:center}
td:first-child,th:first-child{text-align:left}
img{max-width:100%;margin:10px 0;border:1px solid #edf2f7;border-radius:4px}
.note{background:#f7fafc;border-left:3px solid #4a5568;padding:9px 13px;
 font-size:12.5px;margin:14px 0}
.warn{background:#fffaf0;border-left:3px solid #b7791f}
code{background:#f1f5f9;padding:1px 4px;border-radius:3px;font-size:11px}
#wrap{position:relative;border:1px solid #e2e8f0;border-radius:6px;padding:6px}
#pca{width:100%;height:520px;display:block;cursor:crosshair}
#tip{position:absolute;display:none;background:#fff;border:1px solid #cbd5e0;
 border-radius:5px;padding:7px 10px;font-size:11.5px;pointer-events:none;
 box-shadow:0 2px 8px rgba(0,0,0,.13);z-index:5;max-width:260px}
.bar{display:flex;gap:9px;align-items:center;flex-wrap:wrap;margin:6px 0 10px}
.bar input[type=text]{flex:1;min-width:220px;padding:6px 10px;font-size:13px;
 border:1px solid #cbd5e0;border-radius:5px}
.bar select,.bar button{padding:6px 10px;font-size:12.5px;border:1px solid #cbd5e0;
 border-radius:5px;background:#fff;cursor:pointer}
.sw{display:inline-block;width:20px;height:11px;margin:0 4px;vertical-align:-1px;
 border:1px solid #cbd5e0}
#legend{font-size:11.5px;color:#4a5568;margin-top:4px}
#nfound{font-size:12px;color:#718096;min-width:74px}
"""


def _table(df: pd.DataFrame, fmt: str = "{:.3g}") -> str:
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda v: "" if pd.isna(v) else fmt.format(v))
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in d.columns)
    body = "".join("<tr>" + "".join(f"<td>{html.escape(str(v))}</td>" for v in row)
                   + "</tr>" for row in d.itertuples(index=False))
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def write_html(path: str, ctx: dict) -> str:
    pts = json.dumps(ctx["points"], separators=(",", ":"))
    var = ("{xlab:%s,ylab:%s,fmt:function(v){return (Math.abs(v)>=100||"
           "(Math.abs(v)<0.01&&v!==0))?v.toExponential(1):v.toFixed(1);}}"
           % (json.dumps(ctx["xlab"]), json.dumps(ctx["ylab"])))
    js = _PCA_JS.replace("__DATA__", pts).replace("__VAR__", var)

    p: List[str] = [
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
        "<title>Intronaute report</title>",
        f"<style>{CSS}</style></head><body>",
        "<h1>Intronaute &mdash; minor-intron retention</h1>",
        f"<div class='sub'>{ctx['n_samples']} samples &middot; "
        f"{ctx['n_minor']} minor introns &middot; {ctx['n_control']} control "
        f"introns &middot; {datetime.datetime.now():%Y-%m-%d %H:%M}</div>",

        "<h2>1. Sample map (PCA of minor introns)</h2>",
        "<div class='sub'>Every point is a sample, positioned by its minor-intron "
        "retention profile. <b>Type in the box</b> to highlight samples, "
        "scroll to zoom, drag to pan, click a point to pin its label. "
        "Colour follows the selected metric.</div>",
        "<div class='bar'>",
        "<input type='text' id='search' placeholder="
        "'search sample or run (space-separated terms)'>",
        "<span id='nfound'></span>",
        "<select id='metric'>"
        "<option value='n_hits'>colour: n_hits</option>"
        "<option value='tail_q'>colour: -log10(tail_q)</option>"
        "<option value='enrichment'>colour: log2(enrichment)</option>"
        "<option value='MSRI'>colour: MSRI</option>"
        "<option value='stouffer'>colour: stouffer</option></select>",
        "<button id='topn'>pin top 10</button>",
        "<button id='clear'>reset</button>",
        "</div>",
        "<div id='wrap'><canvas id='pca'></canvas><div id='tip'></div></div>",
        "<div id='legend'></div>",
        f"<script>{js}</script>",
    ]

    p += ["<h2>2. Highest-ranking samples</h2>",
          "<div class='sub'>Compare <code>n_hits</code> with "
          "<code>n_hits_control</code>: similar values mean generic intronic "
          "noise, not a minor-spliceosome defect. <code>tail_q</code> tests "
          "exactly that contrast within each sample.</div>",
          _table(ctx["top_samples"]),
          f"<img src='data:image/png;base64,{_b64(ctx['fig_scores'])}'>"]

    if ctx.get("fig_heatmap"):
        p += ["<h2>3. Most retained introns</h2>",
              f"<img src='data:image/png;base64,{_b64(ctx['fig_heatmap'])}'>"]

    if len(ctx.get("top_events", [])):
        p += ["<h2>4. Top events</h2>",
              "<div class='sub'>Per-intron calls are much less reliable than the "
              "aggregate score; treat them as candidates to confirm (IGV, then "
              "RT-PCR). <code>n_EI_donor</code>/<code>n_EI_acceptor</code> are "
              "reads crossing each splice boundary: retention without them is "
              "more likely an overlapping transcript or a mapping artefact.</div>",
              _table(ctx["top_events"])]

    p += ["<h2>5. Quality control</h2>",
          "<div class='sub'>A sample with low exonic depth, or with an intronic "
          "background far above the rest, should be read with caution. If runs "
          "separate here, that is a batch effect, and the within-sample "
          "<code>tail_q</code> is the statistic to trust.</div>",
          f"<img src='data:image/png;base64,{_b64(ctx['fig_qc'])}'>"]
    if len(ctx.get("run_table", [])):
        p.append(_table(ctx["run_table"]))

    p += ["<h2>6. Parameters</h2>", _table(ctx["params"]), "</body></html>"]

    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(p))
    return path
