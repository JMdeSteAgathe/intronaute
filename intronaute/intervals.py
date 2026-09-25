"""Small interval helpers. Coordinates are 0-based, half-open: [start, end).

Pure Python / stdlib only: no dependencies, no compiler.
"""
from __future__ import annotations

from typing import Iterable, List, Sequence, Tuple

Interval = Tuple[int, int]


def merge(iv: Iterable[Interval]) -> List[Interval]:
    """Merge overlapping or adjacent intervals."""
    s = sorted((a, b) for a, b in iv if b > a)
    if not s:
        return []
    out = [list(s[0])]
    for a, b in s[1:]:
        if a <= out[-1][1]:
            if b > out[-1][1]:
                out[-1][1] = b
        else:
            out.append([a, b])
    return [(a, b) for a, b in out]


def subtract(base: Sequence[Interval], holes: Sequence[Interval]) -> List[Interval]:
    """Subtract `holes` from `base`."""
    holes = merge(holes)
    out: List[Interval] = []
    for a, b in merge(base):
        cur = a
        for ha, hb in holes:
            if hb <= cur:
                continue
            if ha >= b:
                break
            if ha > cur:
                out.append((cur, min(ha, b)))
            cur = max(cur, hb)
            if cur >= b:
                break
        if cur < b:
            out.append((cur, b))
    return out


def intersect(a: Sequence[Interval], b: Sequence[Interval]) -> List[Interval]:
    a, b = merge(a), merge(b)
    out: List[Interval] = []
    i = j = 0
    while i < len(a) and j < len(b):
        s = max(a[i][0], b[j][0])
        e = min(a[i][1], b[j][1])
        if e > s:
            out.append((s, e))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out


def total_len(iv: Sequence[Interval]) -> int:
    return sum(b - a for a, b in merge(iv))


def overlaps(a: Interval, b: Interval) -> bool:
    return a[0] < b[1] and b[0] < a[1]


def trim_inner(iv: Sequence[Interval], pad: int, side: str) -> List[Interval]:
    """Trim `pad` bp from the inner edge.

    side='right': drop the last `pad` bp (right edge of the rightmost block)
    side='left' : drop the first `pad` bp
    """
    iv = merge(iv)
    if not iv or pad <= 0:
        return list(iv)
    if side == "right":
        limit = iv[-1][1] - pad
        return [(a, min(b, limit)) for a, b in iv if min(b, limit) > a]
    limit = iv[0][0] + pad
    return [(max(a, limit), b) for a, b in iv if b > max(a, limit)]


def window_from_edge(iv: Sequence[Interval], width: int, side: str) -> List[Interval]:
    """Keep at most `width` bp measured from the inner edge.

    side='right': the last `width` bp of the union (upstream exon of an intron)
    side='left' : the first `width` bp
    """
    iv = merge(iv)
    if not iv:
        return []
    out: List[Interval] = []
    remaining = width
    if side == "right":
        for a, b in reversed(iv):
            if remaining <= 0:
                break
            take = min(b - a, remaining)
            out.append((b - take, b))
            remaining -= take
    else:
        for a, b in iv:
            if remaining <= 0:
                break
            take = min(b - a, remaining)
            out.append((a, a + take))
            remaining -= take
    return merge(out)


def blocks_to_str(iv: Sequence[Interval]) -> str:
    return ";".join(f"{a}-{b}" for a, b in merge(iv))


def str_to_blocks(s: str) -> List[Interval]:
    if not s or (isinstance(s, float)):
        return []
    s = str(s).strip()
    if not s or s == "nan":
        return []
    out = []
    for part in s.split(";"):
        a, b = part.split("-")
        out.append((int(a), int(b)))
    return out


def span(iv: Sequence[Interval]) -> Interval:
    iv = merge(iv)
    return (iv[0][0], iv[-1][1])
