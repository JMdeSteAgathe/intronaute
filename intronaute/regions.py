"""Turn a minor-intron list + a GTF into the table of measurement targets.

Accepted input formats for the intron list
------------------------------------------
BED / BED.gz      chrom, start (0-based), end, [name], [score], [strand]
TSV / CSV / TXT   a header with chrom/start/end columns under any of the usual
                  spellings, plus optional strand / gene columns; OR a single
                  column holding "chr1:1000-2000"
ODS / XLSX        same column logic, from a spreadsheet

Coordinate base: BED is always 0-based half-open. For the other formats the
default is 1-based inclusive (what genome browsers and most published intron
tables use); override with --coord-base 0 if yours is BED-like.

Strand and gene are optional: when missing they are taken from the host gene
found in the GTF.

For every intron we derive
  - "intron core" blocks: the intron, trimmed at both ends, MINUS anything that
    is exonic in another transcript or another gene (otherwise we would be
    measuring exon, not retention);
  - "reference exon" blocks: a window of the last N bp of the upstream exon and
    the first N bp of the downstream exon, cleaned the same way;
  - the coordinates of the expected exon-exon junction.

We also pick, for each minor intron, 1-3 major (U2) introns from the SAME gene
with a comparable length. These internal controls are what separates a genuine
minor-spliceosome defect from ordinary intronic background (RNA degradation,
DNA contamination, nuclear pre-mRNA).

Internal coordinates are 0-based half-open throughout.
"""
from __future__ import annotations

import gzip
import math
import os
import re
import sys
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from . import intervals as iv

BIN = 4096  # bin size for the GTF pre-filter

_POS_RE = re.compile(r"^(?P<chrom>[\w.]+)\s*:\s*(?P<start>[\d,_]+)\s*[-.]+\s*"
                     r"(?P<end>[\d,_]+)")

# Column aliases, compared after stripping everything but [a-z0-9]: so
# "Gene symbol", "gene_symbol" and "GENE-SYMBOL" all collapse to "genesymbol".
_ALIASES = {
    "chrom": ("chrom", "chr", "chromosome", "seqname", "seqnames", "contig",
              "reference", "chrname", "chromname"),
    "start": ("start", "chromstart", "begin", "intronstart", "istart", "from",
              "startpos", "left"),
    "end": ("end", "chromend", "stop", "intronend", "iend", "to", "endpos",
            "right"),
    "strand": ("strand", "sense", "brin"),
    "gene": ("gene", "genename", "genesymbol", "symbol", "geneid", "name",
             "hgnc", "hgncsymbol"),
    "pos": ("pos", "position", "coord", "coords", "coordinates", "region",
            "locus", "interval", "intron"),
}


def _canon(col: str) -> Optional[str]:
    c = re.sub(r"[^a-z0-9]", "", str(col).strip().lower())
    for key, names in _ALIASES.items():
        if c in names:
            return key
    return None


def _open_maybe_gz(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


# --------------------------------------------------------------------------- #
# Reading the intron list
# --------------------------------------------------------------------------- #
def load_introns(path: str, coord_base: str = "auto",
                 sheet: Optional[str] = None) -> List[dict]:
    """Read a minor-intron list from BED / TSV / CSV / ODS / XLSX."""
    ext = "".join(os.path.splitext(path.replace(".gz", ""))[1:]).lower()
    is_bed = ext in (".bed", ".bed6", ".bed12")

    if ext in (".ods", ".xlsx", ".xlsm", ".xls"):
        import pandas as pd

        engine = "odf" if ext == ".ods" else None
        df = pd.read_excel(path, engine=engine, sheet_name=sheet or 0)
        rows = _from_frame(df, base=(0 if coord_base == "0" else 1))
    elif is_bed:
        rows = _from_bed(path)
    else:
        import pandas as pd

        df = pd.read_csv(path, sep=None, engine="python", comment="#"
                         if _looks_commented(path) else None)
        rows = _from_frame(df, base=(0 if coord_base == "0" else 1))

    if coord_base == "0" and not is_bed:
        pass  # already applied
    if not rows:
        raise SystemExit(f"No usable interval found in {path}")
    return rows


def _looks_commented(path: str) -> bool:
    with _open_maybe_gz(path) as fh:
        return fh.readline().startswith("#")


def _from_bed(path: str) -> List[dict]:
    rows = []
    with _open_maybe_gz(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                f = line.split()
            if len(f) < 3:
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                continue  # header line
            rows.append(dict(
                chrom=f[0], start=start, end=end,                  # BED: 0-based
                strand=(f[5].strip() if len(f) > 5 and f[5].strip() in "+-" else ""),
                gene=(f[3].strip() if len(f) > 3 else ""),
                extra=""))
    return rows


def _from_frame(df, base: int) -> List[dict]:
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    found: Dict[str, str] = {}
    for c in df.columns:
        k = _canon(c)
        if k and k not in found:
            found[k] = c

    rows: List[dict] = []
    if {"chrom", "start", "end"} <= set(found):
        for r in df.itertuples(index=False):
            d = dict(zip(df.columns, r))
            try:
                s, e = int(d[found["start"]]), int(d[found["end"]])
            except (TypeError, ValueError):
                continue
            rows.append(dict(
                chrom=str(d[found["chrom"]]).strip(),
                start=s - (1 if base == 1 else 0), end=e,
                strand=str(d.get(found.get("strand", ""), "") or "").strip(),
                gene=str(d.get(found.get("gene", ""), "") or "").strip(),
                extra=""))
    elif "pos" in found:
        for r in df.itertuples(index=False):
            d = dict(zip(df.columns, r))
            m = _POS_RE.match(str(d[found["pos"]]).strip())
            if not m:
                continue
            s = int(m.group("start").replace(",", "").replace("_", ""))
            e = int(m.group("end").replace(",", "").replace("_", ""))
            rows.append(dict(
                chrom=m.group("chrom"),
                start=s - (1 if base == 1 else 0), end=e,
                strand=str(d.get(found.get("strand", ""), "") or "").strip(),
                gene=str(d.get(found.get("gene", ""), "") or "").strip(),
                extra=""))
    else:
        raise SystemExit(
            "Could not find coordinate columns.\n"
            f"  columns present : {list(df.columns)}\n"
            "  expected either chrom/start/end (any usual spelling) or a single\n"
            "  column such as 'pos' holding 'chr1:1000-2000'.")
    return rows


# --------------------------------------------------------------------------- #
# GTF (two passes, spatially filtered to stay light on RAM)
# --------------------------------------------------------------------------- #
def _attr(attrs: str, key: str) -> Optional[str]:
    k = key + ' "'
    i = attrs.find(k)
    if i < 0:
        return None
    i += len(k)
    return attrs[i:attrs.find('"', i)]


def read_gtf_genes(gtf: str) -> List[dict]:
    """Pass 1: only `gene` lines (~80k) -- very light."""
    out = []
    with _open_maybe_gz(gtf) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            gid = _attr(f[8], "gene_id") or ""
            out.append(dict(
                chrom=f[0], start=int(f[3]) - 1, end=int(f[4]), strand=f[6],
                gene_id=gid.split(".")[0],
                gene_name=_attr(f[8], "gene_name") or gid.split(".")[0],
                gene_type=_attr(f[8], "gene_type")
                or _attr(f[8], "gene_biotype") or ""))
    return out


def read_gtf_exons_in_roi(gtf: str, roi_bins: Dict[str, set]) -> List[tuple]:
    """Pass 2: exons intersecting the regions of interest.

    The bin pre-filter avoids holding 3 million exons in memory.
    """
    out = []
    with _open_maybe_gz(gtf) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            bins = roi_bins.get(f[0])
            if not bins:
                continue
            s, e = int(f[3]) - 1, int(f[4])
            if not any(b in bins for b in range(s // BIN, (e - 1) // BIN + 1)):
                continue
            a = f[8]
            out.append((f[0], s, e, (_attr(a, "gene_id") or "").split(".")[0],
                        _attr(a, "transcript_id") or "", f[6]))
    return out


# --------------------------------------------------------------------------- #
def _pick_flanks(tx_introns, tx_exons, start, end):
    """Find the exons flanking intron [start,end).

    Transcripts containing EXACTLY this intron win; otherwise fall back to the
    nearest exon on each side.
    """
    left, right, exact = [], [], False
    for tid, introns in tx_introns.items():
        if (start, end) in introns:
            exact = True
            for a, b in tx_exons[tid]:
                if b == start:
                    left.append((a, b))
                if a == end:
                    right.append((a, b))
    if exact and left and right:
        return iv.merge(left), iv.merge(right), True

    best_l = best_r = None
    for _tid, ex in tx_exons.items():
        for a, b in ex:
            if b <= start and (best_l is None or b > best_l[1]):
                best_l = (a, b)
            if a >= end and (best_r is None or a < best_r[0]):
                best_r = (a, b)
    return (iv.merge([best_l]) if best_l else [],
            iv.merge([best_r]) if best_r else [], False)


def build_targets(introns_path: str, gtf: str, out_tsv: str,
                  exon_window: int = 300, intron_pad: int = 12,
                  exon_pad: int = 5, min_core_len: int = 40,
                  min_exon_ref: int = 50, max_core_window: int = 2000,
                  n_controls: int = 2, min_ctrl_len: int = 60,
                  mask_strand: str = "both", coord_base: str = "auto",
                  sheet: Optional[str] = None) -> str:
    import pandas as pd

    print("[1/5] reading the minor-intron list...", flush=True)
    minor = load_introns(introns_path, coord_base=coord_base, sheet=sheet)
    print(f"      {len(minor)} intervals", flush=True)

    print("[2/5] reading genes from the GTF...", flush=True)
    genes = read_gtf_genes(gtf)
    by_id = {g["gene_id"]: g for g in genes}
    by_name: Dict[str, List[dict]] = defaultdict(list)
    for g in genes:
        by_name[g["gene_name"]].append(g)
    by_chrom: Dict[str, List[dict]] = defaultdict(list)
    for g in genes:
        by_chrom[g["chrom"]].append(g)
    print(f"      {len(genes)} genes", flush=True)

    # --- host gene + regions of interest ---------------------------------- #
    roi_bins: Dict[str, set] = defaultdict(set)
    host: Dict[int, dict] = {}
    n_nohost = n_gene_from_overlap = n_strand_from_gene = 0
    for k, r in enumerate(minor):
        g = by_id.get(r["gene"].split(".")[0]) if r["gene"] else None
        if g is None and r["gene"]:
            g = next((x for x in by_name.get(r["gene"], [])
                      if x["chrom"] == r["chrom"]), None)
        if g is None:
            # no usable gene column: take the host gene by overlap
            cands = [x for x in by_chrom.get(r["chrom"], ())
                     if x["start"] <= r["start"] and r["end"] <= x["end"]
                     and (not r["strand"] or x["strand"] == r["strand"])]
            if cands:
                g = min(cands, key=lambda x: x["end"] - x["start"])
                n_gene_from_overlap += 1
        if g is None:
            g = dict(chrom=r["chrom"], start=r["start"], end=r["end"],
                     strand=r["strand"] or "+", gene_id=r["gene"] or f"NA_{k}",
                     gene_name=r["gene"] or f"intron_{k}", gene_type="")
            n_nohost += 1
        if not r["strand"]:
            r["strand"] = g["strand"]
            n_strand_from_gene += 1
        host[k] = g
        s, e = min(g["start"], r["start"]) - 5000, max(g["end"], r["end"]) + 5000
        for b in range(max(0, s) // BIN, (e - 1) // BIN + 1):
            roi_bins[r["chrom"]].add(b)

    if n_gene_from_overlap:
        print(f"      {n_gene_from_overlap} host genes inferred by overlap", flush=True)
    if n_strand_from_gene:
        print(f"      {n_strand_from_gene} strands taken from the host gene", flush=True)
    if n_nohost:
        print(f"      [warn] {n_nohost} introns with no host gene in the GTF",
              flush=True)

    print("[3/5] extracting exons in the regions of interest "
          "(streaming the GTF)...", flush=True)
    exons = read_gtf_exons_in_roi(gtf, roi_bins)
    print(f"      {len(exons)} exons kept", flush=True)

    exon_bins: Dict[Tuple[str, int], List[tuple]] = defaultdict(list)
    tx_exons_by_gene: Dict[str, Dict[str, List[Tuple[int, int]]]] = defaultdict(
        lambda: defaultdict(list))
    for chrom, s, e, gid, tid, strd in exons:
        for b in range(s // BIN, (e - 1) // BIN + 1):
            exon_bins[(chrom, b)].append((s, e, gid, strd))
        if tid:
            tx_exons_by_gene[gid][tid].append((s, e))

    def exons_over(chrom, s, e, exclude_gene=None, strand=None):
        out = []
        for b in range(s // BIN, (e - 1) // BIN + 1):
            for xs, xe, gid, strd in exon_bins.get((chrom, b), ()):
                if xe <= s or xs >= e:
                    continue
                if exclude_gene is not None and gid == exclude_gene:
                    continue
                if strand is not None and strd != strand:
                    continue
                out.append((xs, xe))
        return iv.merge(out)

    minor_set: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for r in minor:
        minor_set[r["chrom"]].append((r["start"], r["end"]))
    for c in minor_set:
        minor_set[c] = iv.merge(minor_set[c])

    print("[4/5] building targets (intron core, reference exons, "
          "U2 control introns)...", flush=True)
    rows: List[dict] = []
    n_exact = n_drop = 0

    for k, r in enumerate(minor):
        g = host[k]
        gid, chrom = g["gene_id"], r["chrom"]
        tx_exons = {t: iv.merge(v) for t, v in tx_exons_by_gene.get(gid, {}).items()}
        tx_introns = {t: [(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)]
                      for t, ex in tx_exons.items()}

        left, right, exact = _pick_flanks(tx_introns, tx_exons, r["start"], r["end"])
        n_exact += int(exact)

        rec = _make_target(
            kind="minor", target_id=f"MI{k:05d}_{g['gene_name']}",
            pair_id=f"MI{k:05d}", gene=g, chrom=chrom, start=r["start"],
            end=r["end"], strand=r["strand"], left=left, right=right,
            exons_over=exons_over, exon_window=exon_window,
            intron_pad=intron_pad, exon_pad=exon_pad, min_core_len=min_core_len,
            max_core_window=max_core_window, min_exon_ref=min_exon_ref,
            mask_strand=mask_strand)
        rec["flanks_exact"] = int(exact)
        n_drop += int(not rec["usable"])
        rows.append(rec)

        if n_controls <= 0:
            continue
        cand = {}
        for _t, introns in tx_introns.items():
            for a, b in introns:
                if (a, b) == (r["start"], r["end"]) or b - a < min_ctrl_len:
                    continue
                if any(a < ms[1] and ms[0] < b for ms in minor_set.get(chrom, ())):
                    continue  # overlaps a known minor intron
                cand[(a, b)] = 1
        if not cand:
            continue
        tgt = math.log(max(r["end"] - r["start"], 1))
        taken = 0
        for (a, b) in sorted(cand, key=lambda ab: abs(math.log(ab[1] - ab[0]) - tgt)):
            if taken >= n_controls:
                break
            lft, rgt, ex2 = _pick_flanks(tx_introns, tx_exons, a, b)
            crec = _make_target(
                kind="control", target_id=f"CT{k:05d}_{taken+1}_{g['gene_name']}",
                pair_id=f"MI{k:05d}", gene=g, chrom=chrom, start=a, end=b,
                strand=r["strand"], left=lft, right=rgt, exons_over=exons_over,
                exon_window=exon_window, intron_pad=intron_pad,
                exon_pad=exon_pad, min_core_len=min_core_len,
                max_core_window=max_core_window, min_exon_ref=min_exon_ref,
                mask_strand=mask_strand)
            if not crec["usable"]:
                continue
            crec["flanks_exact"] = int(ex2)
            rows.append(crec)
            taken += 1

    cols = ["target_id", "kind", "pair_id", "gene_id", "gene_name", "gene_type",
            "chrom", "intron_start", "intron_end", "strand", "intron_len",
            "core_blocks", "core_len", "exon_blocks", "exon_len",
            "left_exon_blocks", "right_exon_blocks", "fetch_blocks",
            "fetch_start", "fetch_end", "usable", "flanks_exact",
            "exon_ref_masked", "mask_strand"]
    df = pd.DataFrame(rows)[cols].sort_values(
        ["chrom", "intron_start", "target_id"]).reset_index(drop=True)
    os.makedirs(os.path.dirname(os.path.abspath(out_tsv)) or ".", exist_ok=True)
    df.to_csv(out_tsv, sep="\t", index=False)

    nm = int((df["kind"] == "minor").sum())
    print("[5/5] done.")
    print(f"      minor introns      : {nm} "
          f"({int(((df['kind'] == 'minor') & (df['usable'] == 1)).sum())} usable)")
    print(f"      exact flanking exons: {n_exact}/{len(minor)}")
    print(f"      U2 control introns : {int((df['kind'] == 'control').sum())}")
    print(f"      total targets      : {len(df)}")
    print(f"      -> {out_tsv}")
    return out_tsv


def _make_target(kind, target_id, pair_id, gene, chrom, start, end, strand,
                 left, right, exons_over, exon_window, intron_pad, exon_pad,
                 min_core_len, max_core_window, min_exon_ref, mask_strand):
    gid = gene["gene_id"]
    ilen = end - start

    # --- intron core ------------------------------------------------------- #
    # Adaptive trimming: on a 70 bp intron, removing 12 bp from each end would
    # leave almost nothing to measure.
    pad = min(intron_pad, max(4, ilen // 6))
    core = [(start + pad, end - pad)] if ilen > 2 * pad + 10 else []
    # Huge introns: measure only the windows next to the two splice sites.
    # That is where retention shows, and it avoids reading 500 kb per target.
    if core and max_core_window and ilen > 2 * max_core_window:
        a, b = core[0]
        core = [(a, a + max_core_window), (b - max_core_window, b)]
    if core:
        # Remove anything exonic in another transcript / another gene, else we
        # would measure exon rather than retention.
        #   mask_strand='both'  (default) also masks ANTISENSE exons (typically
        #     a 3'UTR of an overlapping gene). Safe even for unstranded data.
        #   mask_strand='sense' masks only the same strand and relies on the
        #     strandedness filter at quantification time. Recovers measurable
        #     length but is ONLY valid if every library is stranded.
        core = iv.subtract(core, exons_over(
            chrom, core[0][0], core[-1][1],
            strand=(strand if mask_strand == "sense" else None)))
    core_len = iv.total_len(core)

    # --- reference exons --------------------------------------------------- #
    lb = iv.window_from_edge(left, exon_window, "right") if left else []
    rb = iv.window_from_edge(right, exon_window, "left") if right else []
    lb = iv.trim_inner(lb, exon_pad, "right")
    rb = iv.trim_inner(rb, exon_pad, "left")
    masked = 1
    lb2 = iv.subtract(lb, exons_over(chrom, *iv.span(lb), exclude_gene=gid)) if lb else []
    rb2 = iv.subtract(rb, exons_over(chrom, *iv.span(rb), exclude_gene=gid)) if rb else []
    if iv.total_len(lb2) + iv.total_len(rb2) >= min_exon_ref:
        lb, rb = lb2, rb2
    else:
        masked = 0  # masking would destroy the reference: keep it unmasked
    exon_blocks = iv.merge(list(lb) + list(rb))
    exon_len = iv.total_len(exon_blocks)

    # --- what to actually read from the BAM -------------------------------- #
    fetch = iv.merge(list(core) + list(exon_blocks)
                     + [(start - pad - 2, start + pad + 2),
                        (end - pad - 2, end + pad + 2)])

    return dict(
        target_id=target_id, kind=kind, pair_id=pair_id, gene_id=gid,
        gene_name=gene["gene_name"], gene_type=gene.get("gene_type", ""),
        chrom=chrom, intron_start=start, intron_end=end, strand=strand,
        intron_len=ilen, core_blocks=iv.blocks_to_str(core), core_len=core_len,
        exon_blocks=iv.blocks_to_str(exon_blocks), exon_len=exon_len,
        left_exon_blocks=iv.blocks_to_str(lb),
        right_exon_blocks=iv.blocks_to_str(rb),
        fetch_blocks=iv.blocks_to_str(fetch),
        fetch_start=fetch[0][0], fetch_end=fetch[-1][1],
        exon_ref_masked=masked, mask_strand=mask_strand, flanks_exact=0,
        usable=1 if (core_len >= min_core_len and exon_len >= min_exon_ref) else 0)
