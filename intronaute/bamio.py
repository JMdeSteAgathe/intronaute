"""BAM access without a compiler or admin rights.

Three read engines, all exposing the SAME interface: `fetch()` yields tuples of
`(flag, position, mapq, cigartuples)`. Nothing else is decoded, which is all we
need to measure coverage and is much faster.

  fastbam   : the minimal BAM/BAI reader bundled here, pure Python (default)
  pysam     : used when importable (Linux/macOS/WSL), slightly faster
  bamnostic : opt-in fallback only -- measured to under-count coverage for
              reads whose splice junction spans into the queried region

On Windows there is no pysam wheel, so `fastbam` does the work. It has been
verified to give counts identical to pysam, bit for bit.
"""
from __future__ import annotations

import os
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

BACKEND: Optional[str] = None

# SAM FLAG bits
F_PAIRED, F_UNMAP, F_REVERSE, F_READ2 = 0x1, 0x4, 0x10, 0x80
F_SECONDARY, F_QCFAIL, F_DUP, F_SUPPL = 0x100, 0x200, 0x400, 0x800
BAD_BASE = F_UNMAP | F_SECONDARY | F_SUPPL | F_QCFAIL


def longpath(p: str) -> str:
    """Work around the legacy 260-character path limit on Windows.

    Useful for deep UNC paths such as \\\\server\\share\\...; a no-op elsewhere.
    """
    if os.name != "nt" or len(p) < 240 or p.startswith("\\\\?\\"):
        return p
    if p.startswith("\\\\"):
        return "\\\\?\\UNC\\" + p[2:]
    return "\\\\?\\" + os.path.abspath(p)


def index_paths(bam: str) -> List[str]:
    return [bam + ".bai", os.path.splitext(bam)[0] + ".bai",
            bam + ".csi", os.path.splitext(bam)[0] + ".csi"]


def has_index(bam: str) -> bool:
    return any(os.path.exists(longpath(p)) for p in index_paths(bam))


def find_bai(bam: str) -> Optional[str]:
    for p in (bam + ".bai", os.path.splitext(bam)[0] + ".bai"):
        if os.path.exists(longpath(p)):
            return longpath(p)
    return None


class BamReader:
    """Uniform wrapper. `fetch` yields (flag, pos, mapq, cigartuples)."""

    def __init__(self, path: str, prefer: str = "auto"):
        self.path = path
        p = longpath(path)
        if not os.path.exists(p):
            raise FileNotFoundError(path)
        if not has_index(path):
            raise FileNotFoundError(
                f"No BAM index found for {path}\n"
                f"  expected: {path}.bai\n"
                f"  (region queries are impossible without an index; "
                f"run `samtools index`)")

        order = {"auto": ["pysam", "fastbam"]}.get(prefer, [prefer])
        errs = []
        for name in order:
            try:
                getattr(self, f"_open_{name}")(p)
                self.backend = name
                break
            except Exception as e:
                errs.append(f"{name}: {e}")
        else:
            raise SystemExit("No usable BAM read engine:\n  " + "\n  ".join(errs))
        global BACKEND
        BACKEND = self.backend
        self._chrom_map: Dict[str, Optional[str]] = {}

    # ---- engines ---------------------------------------------------------- #
    def _open_fastbam(self, p: str) -> None:
        from .fastbam import FastBam

        bai = find_bai(self.path)
        if bai is None:
            raise ValueError("no .bai index (.csi is not supported by fastbam)")
        self._fh = FastBam(p, index=bai)
        self.references = list(self._fh.references)
        self._fetch = self._fh.fetch

    def _open_pysam(self, p: str) -> None:
        import pysam

        self._fh = pysam.AlignmentFile(p, "rb")
        self.references = list(self._fh.references)

        def f(c, s, e):
            for r in self._fh.fetch(c, s, e):
                yield (r.flag, r.reference_start, r.mapping_quality,
                       list(r.cigartuples or ()))

        self._fetch = f

    def _open_bamnostic(self, p: str) -> None:
        # Opt-in fallback only (--engine bamnostic), e.g. for a .csi index.
        # Measured here: bamnostic misses spliced reads whose junction reaches
        # into the queried region from far away, which under-counts coverage.
        import bamnostic

        self._fh = bamnostic.AlignmentFile(p, "rb")
        # bamnostic coerces "2" to an int in header['SQ'] but keeps strings in
        # .references
        self.references = [str(x) for x in self._fh.references]
        print("    [warn] bamnostic engine: coverage is under-estimated for "
              "reads with long splice junctions; use only as a fallback.",
              flush=True)

        def f(c, s, e):
            for r in self._fh.fetch(contig=c, start=s, stop=e):
                yield (r.flag, r.reference_start, r.mapping_quality,
                       list(r.cigartuples or ()))

        self._fetch = f

    # ---- chromosome naming (chr1 / 1 / MT) -------------------------------- #
    def resolve_chrom(self, chrom: str) -> Optional[str]:
        if chrom in self._chrom_map:
            return self._chrom_map[chrom]
        refs = set(self.references)
        cands = [chrom]
        if chrom.startswith("chr"):
            bare = chrom[3:]
            cands += [bare, "MT" if bare == "M" else ""]
        else:
            cands += ["chr" + chrom, "chrM" if chrom == "MT" else ""]
        hit = next((c for c in cands if c and c in refs), None)
        self._chrom_map[chrom] = hit
        return hit

    def fetch(self, chrom: str, start: int, end: int) -> Iterator[
            Tuple[int, int, int, list]]:
        ref = self.resolve_chrom(chrom)
        if ref is None:
            return iter(())
        try:
            return self._fetch(ref, max(0, start), end)
        except (IndexError, ValueError, KeyError, StopIteration):
            # bamnostic raises IndexError on an empty region; not an error
            return iter(())

    def close(self) -> None:
        try:
            self._fh.close()
        except Exception:
            pass

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()


# --------------------------------------------------------------------------- #
# CIGAR
# --------------------------------------------------------------------------- #
_ALIGNED = (0, 7, 8)  # M = X : aligned to the reference


def read_blocks(pos: int, cigartuples: Sequence[Tuple[int, int]]):
    """-> (aligned blocks [(s,e)...], splice gaps [(s,e)...])

    A deletion (D) does not count as covered; an N gap is a splice junction.
    """
    blocks: List[Tuple[int, int]] = []
    gaps: List[Tuple[int, int]] = []
    for op, ln in cigartuples:
        if op in _ALIGNED:
            if blocks and blocks[-1][1] == pos:
                blocks[-1] = (blocks[-1][0], pos + ln)
            else:
                blocks.append((pos, pos + ln))
            pos += ln
        elif op == 3:
            gaps.append((pos, pos + ln))
            pos += ln
        elif op == 2:
            pos += ln
        # I / S / H / P do not consume the reference
    return blocks, gaps


def passes_filters(flag: int, mapq: int, min_mapq: int, drop_dup: bool) -> bool:
    bad = BAD_BASE | (F_DUP if drop_dup else 0)
    return not (flag & bad) and mapq >= min_mapq


def read_is_sense(flag: int, gene_strand: str, library: str) -> bool:
    """Does this read come from the gene's strand?

    library: 'reverse' (dUTP / Illumina stranded mRNA, the common case),
             'forward', or 'unstranded' (everything passes).
    """
    if library == "unstranded":
        return True
    rev = bool(flag & F_REVERSE)
    if (flag & F_PAIRED) and (flag & F_READ2):
        rev = not rev
    if library == "reverse":
        rev = not rev
    return ("-" if rev else "+") == gene_strand


def detect_library_type(bam: str, targets, n_probe: int = 120,
                        min_mapq: int = 10, prefer: str = "auto") -> Tuple[str, float]:
    """Estimate strandedness from a few well-covered reference exons."""
    fwd = rev = 0
    with BamReader(bam, prefer=prefer) as fh:
        used = 0
        for t in targets:
            if used >= n_probe:
                break
            if t.strand not in ("+", "-") or not t.exon_blocks:
                continue
            s, e = t.exon_blocks[0]
            n = 0
            for flag, _pos, mapq, _cig in fh.fetch(t.chrom, s, e):
                if not passes_filters(flag, mapq, min_mapq, True):
                    continue
                n += 1
                if n > 200:
                    break
                r = bool(flag & F_REVERSE)
                if (flag & F_PAIRED) and (flag & F_READ2):
                    r = not r
                if ("-" if r else "+") == t.strand:
                    fwd += 1
                else:
                    rev += 1
            if n:
                used += 1
    tot = fwd + rev
    if tot < 200:
        return "unstranded", 0.5
    frac = fwd / tot
    if frac > 0.8:
        return "forward", frac
    if frac < 0.2:
        return "reverse", frac
    return "unstranded", frac
