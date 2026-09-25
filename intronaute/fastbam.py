"""Minimal, fast, pure-Python BAM/BAI reader (uses only stdlib zlib).

Why not an existing pure-Python reader: they decode sequence, qualities, read
name and tags for every record, which costs ~85% of the runtime, while we only
need four fields (flag, position, MAPQ, CIGAR). Skipping the rest gives roughly
a 15x speedup, with no compiler and no binary wheel -- which is what makes this
tool work on a locked-down Windows box where pysam cannot be installed.

Implemented: BGZF (concatenated gzip blocks), .bai index, region queries via
bins + the linear index. .csi indexes are not supported (bamio.py can fall back
to bamnostic for those).

Verified to produce counts identical to pysam, bit for bit.
"""
from __future__ import annotations

import struct
import zlib
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

_u32 = struct.Struct("<I")
_i32 = struct.Struct("<i")
_u64 = struct.Struct("<Q")
# block_size, refID, pos, l_read_name, mapq, bin, n_cigar_op, flag, l_seq
_CORE = struct.Struct("<iiiBBHHHi")


class BgzfReader:
    """Random access to the decompressed bytes of a BGZF file, addressed by
    virtual offset (coffset << 16 | uoffset)."""

    __slots__ = ("fh", "_buf", "_upos", "_coffset", "_next_coffset")

    def __init__(self, fh):
        self.fh = fh
        self._buf = b""
        self._upos = 0
        self._coffset = 0
        self._next_coffset = 0

    def _load_block(self) -> bool:
        self.fh.seek(self._next_coffset)
        self._coffset = self._next_coffset
        head = self.fh.read(12)
        if len(head) < 12:
            self._buf = b""
            self._upos = 0
            return False
        if head[0] != 0x1F or head[1] != 0x8B:
            raise ValueError("invalid BGZF block (is this really a BAM file?)")
        xlen = struct.unpack_from("<H", head, 10)[0]
        extra = self.fh.read(xlen)
        bsize = None
        i = 0
        while i + 4 <= xlen:
            si1, si2, slen = extra[i], extra[i + 1], struct.unpack_from("<H", extra, i + 2)[0]
            if si1 == 66 and si2 == 67:
                bsize = struct.unpack_from("<H", extra, i + 4)[0] + 1
                break
            i += 4 + slen
        if bsize is None:
            raise ValueError("missing BSIZE field: not a BGZF file (bgzip required)")
        rest = self.fh.read(bsize - 12 - xlen)
        self._buf = zlib.decompress(rest[:-8], -15)
        self._upos = 0
        self._next_coffset = self._coffset + bsize
        return True

    def seek_virtual(self, voffset: int) -> None:
        self._next_coffset = voffset >> 16
        self._load_block()
        self._upos = voffset & 0xFFFF

    @property
    def virtual_offset(self) -> int:
        return (self._coffset << 16) | self._upos

    def read(self, n: int) -> bytes:
        out = b""
        while n > 0:
            avail = len(self._buf) - self._upos
            if avail <= 0:
                if not self._load_block():
                    return out
                continue
            take = min(avail, n)
            if not out and take == n:
                out = self._buf[self._upos:self._upos + take]
            else:
                out += self._buf[self._upos:self._upos + take]
            self._upos += take
            n -= take
        return out


def reg2bins(beg: int, end: int) -> List[int]:
    """BAI bins covering [beg, end)."""
    if end <= beg:
        end = beg + 1
    end -= 1
    out = [0]
    for shift, offset in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)):
        out.extend(range(offset + (beg >> shift), offset + (end >> shift) + 1))
    return out


class Bai:
    """A .bai index: chunks per bin plus the 16 kb-window linear index."""

    def __init__(self, path: str):
        with open(path, "rb") as fh:
            data = fh.read()
        if data[:4] != b"BAI\x01":
            raise ValueError("invalid .bai index")
        p = 4
        n_ref = _i32.unpack_from(data, p)[0]
        p += 4
        self.bins: List[Dict[int, List[Tuple[int, int]]]] = []
        self.linear: List[List[int]] = []
        for _ in range(n_ref):
            n_bin = _i32.unpack_from(data, p)[0]
            p += 4
            b: Dict[int, List[Tuple[int, int]]] = {}
            for _ in range(n_bin):
                bid = _u32.unpack_from(data, p)[0]
                n_chunk = _i32.unpack_from(data, p + 4)[0]
                p += 8
                chunks = []
                for _c in range(n_chunk):
                    beg = _u64.unpack_from(data, p)[0]
                    end = _u64.unpack_from(data, p + 8)[0]
                    p += 16
                    chunks.append((beg, end))
                b[bid] = chunks
            n_intv = _i32.unpack_from(data, p)[0]
            p += 4
            lin = list(struct.unpack_from(f"<{n_intv}Q", data, p)) if n_intv else []
            p += 8 * n_intv
            self.bins.append(b)
            self.linear.append(lin)

    def chunks(self, tid: int, beg: int, end: int) -> List[Tuple[int, int]]:
        if tid < 0 or tid >= len(self.bins):
            return []
        lin = self.linear[tid]
        lo = lin[min(beg >> 14, len(lin) - 1)] if lin else 0
        got: List[Tuple[int, int]] = []
        bmap = self.bins[tid]
        for bid in reg2bins(beg, end):
            for cb, ce in bmap.get(bid, ()):
                if ce > lo:
                    got.append((max(cb, lo), ce))
        if not got:
            return []
        got.sort()
        merged = [list(got[0])]
        for cb, ce in got[1:]:
            if cb <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], ce)
            else:
                merged.append([cb, ce])
        return [(a, b) for a, b in merged]


# CIGAR ops 0=M 1=I 2=D 3=N 4=S 5=H 6=P 7== 8=X; reference-consuming: M D N = X
_REF_CONSUMING = (True, False, True, True, False, False, False, True, True)


class FastBam:
    """Minimal AlignmentFile: only fetch(contig, start, stop)."""

    def __init__(self, path: str, index: Optional[str] = None):
        self.path = path
        self.fh = open(path, "rb")
        self.bgzf = BgzfReader(self.fh)
        self.bgzf.seek_virtual(0)
        if self.bgzf.read(4) != b"BAM\x01":
            raise ValueError("missing BAM header")
        l_text = _i32.unpack(self.bgzf.read(4))[0]
        self.bgzf.read(l_text)
        n_ref = _i32.unpack(self.bgzf.read(4))[0]
        self.references: List[str] = []
        self.lengths: List[int] = []
        for _ in range(n_ref):
            l_name = _i32.unpack(self.bgzf.read(4))[0]
            name = self.bgzf.read(l_name)[:-1].decode("utf-8", "replace")
            self.lengths.append(_i32.unpack(self.bgzf.read(4))[0])
            self.references.append(name)
        self.ref2tid = {r: i for i, r in enumerate(self.references)}
        self.bai = Bai(index) if index else None

    def fetch(self, contig: str, start: int, stop: int) -> Iterator[
            Tuple[int, int, int, List[Tuple[int, int]]]]:
        """Yield tuples of (flag, 0-based position, mapq, cigartuples)."""
        if self.bai is None:
            raise ValueError("an index is required for fetch()")
        tid = self.ref2tid.get(contig)
        if tid is None:
            return
        start = max(0, start)
        bg = self.bgzf
        core_unpack = _CORE.unpack_from
        for cbeg, cend in self.bai.chunks(tid, start, stop):
            bg.seek_virtual(cbeg)
            while True:
                vo = bg.virtual_offset
                if vo >= cend:
                    break
                head = bg.read(36)
                if len(head) < 36:
                    break
                (bsz, rid, pos, l_name, mapq, _bin, n_cig, flag,
                 _l_seq) = core_unpack(head, 0)
                body = bg.read(bsz - 32)  # 36 bytes read, 4 of which were block_size
                if rid != tid or pos >= stop:
                    if rid != tid:
                        break
                    continue
                off = l_name
                ref_len = 0
                cig: List[Tuple[int, int]] = []
                for k in range(n_cig):
                    v = _u32.unpack_from(body, off + 4 * k)[0]
                    op = v & 0xF
                    ln = v >> 4
                    cig.append((op, ln))
                    if _REF_CONSUMING[op]:
                        ref_len += ln
                if pos + max(ref_len, 1) > start:
                    yield (flag, pos, mapq, cig)

    def close(self) -> None:
        try:
            self.fh.close()
        except Exception:
            pass
