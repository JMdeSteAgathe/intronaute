"""Per-sample quantification: intronic coverage vs adjacent exons, plus
splice-junction counts.

Designed to be cheap on a modest machine reading over a network share:
  - the BAM is NEVER read end to end: targeted access through the .bai index;
  - nearby targets (introns of one gene plus their controls) are merged into a
    single coordinate-sorted read interval, so each BGZF block is decompressed
    once;
  - one int32 coverage array per merged region (a few hundred kB), so memory
    stays flat regardless of cohort size;
  - results are cached per sample: re-running recomputes nothing.
"""
from __future__ import annotations

import bisect
import gzip
import json
import os
import sys
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from . import bamio
from . import intervals as iv


# --------------------------------------------------------------------------- #
@dataclass
class Target:
    target_id: str
    kind: str
    pair_id: str
    gene_id: str
    gene_name: str
    chrom: str
    intron_start: int
    intron_end: int
    strand: str
    core_blocks: List[Tuple[int, int]]
    exon_blocks: List[Tuple[int, int]]
    left_exon_blocks: List[Tuple[int, int]] = field(default_factory=list)
    right_exon_blocks: List[Tuple[int, int]] = field(default_factory=list)
    fetch_blocks: List[Tuple[int, int]] = field(default_factory=list)
    mask_strand: str = "both"
    fetch_start: int = 0
    fetch_end: int = 0

    def __hash__(self):
        return hash(self.target_id)

    def __eq__(self, other):
        return isinstance(other, Target) and other.target_id == self.target_id


def load_targets(tsv: str, kinds: Optional[Sequence[str]] = None,
                 usable_only: bool = True, min_methods: int = 0) -> List[Target]:
    import pandas as pd

    df = pd.read_csv(tsv, sep="\t", dtype={"chrom": str})
    if usable_only and "usable" in df.columns:
        df = df[df["usable"] == 1]
    if kinds:
        df = df[df["kind"].isin(kinds)]
    if min_methods > 0:
        keep = (df["kind"] != "u12") | (df["n_methods"].fillna(0) >= min_methods)
        pairs = set(df.loc[keep & (df["kind"] == "u12"), "pair_id"])
        df = df[keep & df["pair_id"].isin(pairs)]
    out = []
    for r in df.itertuples(index=False):
        out.append(
            Target(
                target_id=r.target_id, kind=r.kind, pair_id=r.pair_id,
                gene_id=r.gene_id, gene_name=r.gene_name, chrom=str(r.chrom),
                intron_start=int(r.intron_start), intron_end=int(r.intron_end),
                strand=str(r.strand),
                core_blocks=iv.str_to_blocks(r.core_blocks),
                exon_blocks=iv.str_to_blocks(r.exon_blocks),
                left_exon_blocks=iv.str_to_blocks(r.left_exon_blocks),
                right_exon_blocks=iv.str_to_blocks(r.right_exon_blocks),
                fetch_blocks=iv.str_to_blocks(getattr(r, "fetch_blocks", "")) or
                iv.merge(iv.str_to_blocks(r.core_blocks)
                         + iv.str_to_blocks(r.exon_blocks)),
                mask_strand=str(getattr(r, "mask_strand", "both")),
                fetch_start=int(r.fetch_start), fetch_end=int(r.fetch_end),
            )
        )
    return out


def make_chunks(targets: Sequence[Target], max_gap: int = 3000,
                max_span: int = 400_000) -> List[Tuple[str, int, int, List[Target]]]:
    """Build the regions to read.

    We start from the blocks we actually MEASURE (intron core + reference exons
    + splice boundaries) rather than the intron-to-exon span, so the middle of a
    500 kb intron is never read. Blocks closer than `max_gap` are merged so each
    BGZF block is decompressed once.

    `max_gap` must stay well above the read length: that guarantees a read
    cannot fall into two regions at once (no double counting).
    """
    # The key includes the STRAND: a region mixing two overlapping genes on
    # opposite strands would make the strandedness filter wrong.
    by_chrom: Dict[Tuple[str, str], List[Tuple[int, int]]] = {}
    for t in targets:
        by_chrom.setdefault((t.chrom, t.strand), []).extend(t.fetch_blocks)

    regions: Dict[Tuple[str, str], List[Tuple[int, int]]] = {}
    for chrom, blocks in by_chrom.items():
        merged = iv.merge([(max(0, a), b) for a, b in blocks])
        out: List[List[int]] = []
        for a, b in merged:
            if out and a - out[-1][1] <= max_gap and (b - out[-1][0]) <= max_span:
                out[-1][1] = b
            else:
                out.append([a, b])
        regions[chrom] = [(a, b) for a, b in out]

    chunks: List[Tuple[str, int, int, List[Target]]] = []
    index: Dict[Tuple[Tuple[str, str], int], List[Target]] = {}
    for key, regs in regions.items():
        chrom, strand = key
        starts = [r[0] for r in regs]
        for i, (a, b) in enumerate(regs):
            index[(key, i)] = []
        for t in targets:
            if t.chrom != chrom or t.strand != strand:
                continue
            for ba, bb in t.fetch_blocks:
                i = max(0, bisect.bisect_right(starts, ba) - 1)
                while i < len(regs) and regs[i][0] < bb:
                    if regs[i][1] > ba and t not in index[(key, i)]:
                        index[(key, i)].append(t)
                    i += 1
        for i, (a, b) in enumerate(regs):
            if index[(key, i)]:
                chunks.append((chrom, a, b, index[(key, i)]))
    chunks.sort(key=lambda c: (c[0], c[1]))
    return chunks


# --------------------------------------------------------------------------- #
COUNT_COLS = [
    "target_id", "kind", "pair_id", "gene_id", "gene_name", "chrom",
    "intron_start", "intron_end", "strand", "intron_len",
    "core_len", "exon_len",
    "intron_depth_mean", "intron_depth_med", "intron_frac_cov",
    "exon_depth_mean", "exon_depth_med",
    "exon_depth_left", "exon_depth_right",
    "n_intron_reads", "n_EE", "n_EI_left", "n_EI_right",
    "n_alt_donor", "n_alt_acceptor", "truncated",
]


def _cat(segs) -> np.ndarray:
    return np.concatenate(segs) if segs else np.zeros(0, dtype=np.int32)


def _finalize(t: Target, d: dict) -> list:
    dI, dE = _cat(d["dI"]), _cat(d["dE"])
    dL, dR = _cat(d["dL"]), _cat(d["dR"])
    return [
        t.target_id, t.kind, t.pair_id, t.gene_id, t.gene_name, t.chrom,
        t.intron_start, t.intron_end, t.strand, t.intron_end - t.intron_start,
        int(dI.size), int(dE.size),
        round(float(dI.mean()), 4) if dI.size else 0.0,
        float(np.median(dI)) if dI.size else 0.0,
        round(float((dI > 0).mean()), 4) if dI.size else 0.0,
        round(float(dE.mean()), 4) if dE.size else 0.0,
        float(np.median(dE)) if dE.size else 0.0,
        round(float(dL.mean()), 4) if dL.size else 0.0,
        round(float(dR.mean()), 4) if dR.size else 0.0,
        d["n_core"], d["ee"], d["eil"], d["eir"], d["alt_d"], d["alt_a"],
        d["trunc"],
    ]


def quantify_bam(
    bam: str,
    targets: Sequence[Target],
    out_tsv: str,
    min_mapq: int = 10,
    min_overhang: int = 6,
    drop_dup: bool = True,
    library: str = "auto",
    backend: str = "auto",
    max_reads_per_chunk: int = 300_000,
    progress_every: int = 250,
    quiet: bool = False,
) -> dict:
    t0 = time.time()

    if library == "auto":
        library, frac = bamio.detect_library_type(
            bam, [t for t in targets if t.kind == "u12"][:400],
            min_mapq=min_mapq, prefer=backend)
        if not quiet:
            print(f"    detected strandedness: {library} (sense fraction={frac:.2f})", flush=True)

    # Guard: an annotation built with --mask-strand sense assumes the
    # strandedness filter removes antisense reads. With an unstranded library an
    # overlapping antisense 3'UTR would be counted as retention.
    if library == "unstranded" and any(t.mask_strand == "sense" for t in targets):
        print("    [WARNING] annotation built with --mask-strand sense but this "
              "library is NOT stranded: antisense coverage (a 3'UTR of an "
              "overlapping gene) will be counted as retention. Rebuild the "
              "annotation with --mask-strand both.", flush=True)

    chunks = make_chunks(targets)
    rows: List[list] = []
    n_trunc = 0

    # Regions remaining per target -> finalise (and free) as soon as a target's
    # last region has been processed, so memory stays bounded.
    remaining: Dict[str, int] = {}
    for _c, _a, _b, ts in chunks:
        for t in ts:
            remaining[t.target_id] = remaining.get(t.target_id, 0) + 1
    acc: Dict[str, dict] = {}

    with bamio.BamReader(bam, prefer=backend) as fh:
        missing = [c for c in {t.chrom for t in targets} if fh.resolve_chrom(c) is None]
        if missing and not quiet:
            print(f"    [warn] {len(missing)} chromosomes not in the BAM"
                  f" (e.g. {missing[:3]})", flush=True)

        for ci, (chrom, cs, ce, ts) in enumerate(chunks):
            if not quiet and progress_every and ci % progress_every == 0 and ci:
                print(f"    {ci}/{len(chunks)} regions ({time.time()-t0:.0f}s)",
                      flush=True)

            n = ce - cs
            cov = np.zeros(n, dtype=np.int32)        # aligned coverage (M/=/X)
            cross = np.zeros(n + 2, dtype=np.int32)  # blocks spanning a position
            junc: Dict[Tuple[int, int], int] = {}
            cores = [(t.core_blocks, iv.span(t.core_blocks) if t.core_blocks else None)
                     for t in ts]
            nreads_core = [0] * len(ts)
            strand0 = ts[0].strand
            stranded = library != "unstranded"
            nr = 0
            trunc = 0

            for flag, pos, mapq, cig in fh.fetch(chrom, cs, ce):
                if not bamio.passes_filters(flag, mapq, min_mapq, drop_dup):
                    continue
                nr += 1
                if nr > max_reads_per_chunk:
                    trunc = 1
                    break
                if stranded and not bamio.read_is_sense(flag, strand0, library):
                    continue
                blocks, gaps = bamio.read_blocks(pos, cig)
                if not blocks:
                    continue
                for bs, be in blocks:
                    a, b = max(bs, cs) - cs, min(be, ce) - cs
                    if b > a:
                        cov[a:b] += 1
                    ca = max(bs + min_overhang - cs, 0)
                    cb = min(be - min_overhang - cs, n + 1)
                    if cb > ca:
                        cross[ca:cb] += 1
                for g in gaps:
                    # A junction belongs to the region containing its start, so
                    # a spliced read touching two distant regions is not counted
                    # twice.
                    if cs <= g[0] < ce:
                        junc[g] = junc.get(g, 0) + 1
                r0, r1 = blocks[0][0], blocks[-1][1]
                for i, (cb, sp) in enumerate(cores):
                    if sp is None or r1 <= sp[0] or r0 >= sp[1]:
                        continue
                    if any(min(be, y) > max(bs, x)
                           for bs, be in blocks for x, y in cb):
                        nreads_core[i] += 1
            n_trunc += trunc

            def depth(blocks) -> np.ndarray:
                segs = [cov[max(a - cs, 0):max(min(b, ce) - cs, 0)]
                        for a, b in blocks if b > cs and a < ce]
                segs = [s for s in segs if s.size]
                return np.concatenate(segs) if segs else np.zeros(0, dtype=np.int32)

            for i, t in enumerate(ts):
                d = acc.setdefault(t.target_id, dict(
                    dI=[], dE=[], dL=[], dR=[], n_core=0, ee=0,
                    eil=0, eir=0, alt_d=0, alt_a=0, trunc=0))
                for key, blk in (("dI", t.core_blocks), ("dE", t.exon_blocks),
                                 ("dL", t.left_exon_blocks),
                                 ("dR", t.right_exon_blocks)):
                    v = depth(blk)
                    if v.size:
                        d[key].append(v)
                istart, iend = t.intron_start, t.intron_end
                d["n_core"] += nreads_core[i]
                d["trunc"] = max(d["trunc"], trunc)
                if cs <= istart < ce:
                    d["ee"] += junc.get((istart, iend), 0)
                    d["alt_d"] += sum(v for (a, b), v in junc.items()
                                      if a == istart and b != iend)
                    d["eil"] = int(cross[istart - cs])
                if cs <= iend < ce:
                    d["alt_a"] += sum(v for (a, b), v in junc.items()
                                      if b == iend and a != istart)
                    d["eir"] = int(cross[iend - cs])
                elif iend == ce:
                    d["eir"] = int(cross[iend - cs])

                remaining[t.target_id] -= 1
                if remaining[t.target_id] == 0:
                    rows.append(_finalize(t, acc.pop(t.target_id)))

            del cov, cross

    for tid, d in list(acc.items()):  # safety net: targets not yet finalised
        t = next(x for x in targets if x.target_id == tid)
        rows.append(_finalize(t, d))

    import pandas as pd

    df = pd.DataFrame(rows, columns=COUNT_COLS)
    os.makedirs(os.path.dirname(os.path.abspath(out_tsv)) or ".", exist_ok=True)
    df.to_csv(out_tsv, sep="\t", index=False, compression="gzip")

    meta = dict(
        bam=os.path.abspath(bam),
        bam_size=os.path.getsize(bam),
        bam_mtime=os.path.getmtime(bam),
        n_targets=len(df),
        n_chunks=len(chunks),
        n_truncated_chunks=n_trunc,
        library=library,
        min_mapq=min_mapq,
        min_overhang=min_overhang,
        drop_dup=drop_dup,
        backend=bamio.BACKEND,
        seconds=round(time.time() - t0, 1),
        median_exon_depth=float(np.median(df["exon_depth_mean"])) if len(df) else 0.0,
    )
    with open(out_tsv.replace(".tsv.gz", ".json"), "w") as fo:
        json.dump(meta, fo, indent=1)
    return meta


# --------------------------------------------------------------------------- #
# Multiprocessing worker (must be picklable at module level)
# --------------------------------------------------------------------------- #
def _worker(args) -> Tuple[str, Optional[dict], Optional[str]]:
    sample, bam, targets_tsv, out_tsv, params = args
    try:
        targets = load_targets(
            targets_tsv, usable_only=True,
            min_methods=params.pop("min_methods", 0),
            kinds=params.pop("kinds", None),
        )
        meta = quantify_bam(bam, targets, out_tsv, quiet=True, progress_every=0, **params)
        return sample, meta, None
    except Exception:  # one broken BAM must not kill the whole run
        import traceback

        return sample, None, traceback.format_exc(limit=4)


def is_cached(out_tsv: str, bam: str, params: dict) -> bool:
    js = out_tsv.replace(".tsv.gz", ".json")
    if not (os.path.exists(out_tsv) and os.path.exists(js)):
        return False
    try:
        with open(js) as fh:
            m = json.load(fh)
        if int(m.get("bam_size", -1)) != os.path.getsize(bam):
            return False
        for k in ("min_mapq", "min_overhang", "drop_dup"):
            if k in params and m.get(k) != params[k]:
                return False
        return True
    except Exception:
        return False
