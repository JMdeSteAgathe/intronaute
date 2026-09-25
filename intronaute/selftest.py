"""Self-test: validates BAM reading and the counting logic against a tiny
embedded BAM. Run this first on any new machine; it takes about two seconds.

    intronaute selftest
"""
from __future__ import annotations

import base64
import os
import sys
import tempfile

EXPECTED = {
    # target_id: (n_EE, n_EI_left, n_EI_right, n_intron_reads, core_len)
    "MI_test_chr1": (12, 4, 3, 12, 476),
    "MI_test_chr2": (6, 0, 0, 3, 376),
}

TARGETS = """target_id\tkind\tpair_id\tgene_id\tgene_name\tgene_type\tchrom\tintron_start\tintron_end\tstrand\tintron_len\tcore_blocks\tcore_len\texon_blocks\texon_len\tleft_exon_blocks\tright_exon_blocks\tfetch_start\tfetch_end\tusable\tflanks_exact\tdonor\tacceptor\tmethods\tn_methods
MI_test_chr1\tminor\tP1\tENSGTEST1\tTESTGENE1\tprotein_coding\tchr1\t1200\t1700\t+\t500\t1212-1688\t476\t900-1195;1705-2000\t590\t900-1195\t1705-2000\t900\t2000\t1\t1\tAT\tAC\ttest\t1
MI_test_chr2\tminor\tP2\tENSGTEST2\tTESTGENE2\tprotein_coding\tchr2\t500\t900\t+\t400\t512-888\t376\t250-495;905-1150\t491\t250-495\t905-1150\t250\t1150\t1\t1\tGT\tAG\ttest\t1
CT_test_chr1\tcontrol\tP1\tENSGTEST1\tTESTGENE1\tprotein_coding\tchr1\t1200\t1700\t+\t500\t1212-1688\t476\t900-1195;1705-2000\t590\t900-1195\t1705-2000\t900\t2000\t1\t1\t\t\t\t0
"""


def run_selftest() -> None:
    from . import bamio
    from . import quantify as Q
    from ._minibam import MINI_BAI_B64, MINI_BAM_B64

    print(f"Intronaute self-test")
    print(f"  python      : {sys.version.split()[0]} ({sys.platform})")
    print(f"  executable  : {sys.executable}")
    ok = True
    for mod in ("numpy", "pandas", "scipy", "matplotlib"):
        try:
            m = __import__(mod)
            print(f"  {mod:<12}: {getattr(m, '__version__', '?')}")
        except Exception as e:
            print(f"  {mod:<12}: MISSING ({e})")
            ok = False
    # bamnostic is OPTIONAL: the default engine (fastbam) ships with the
    # package. bamnostic is only an explicit fallback.
    try:
        import bamnostic as _bn

        print(f"  bamnostic   : {getattr(_bn, '__version__', '?')} (optional)")
    except Exception:
        print("  bamnostic   : absent (optional, not required)")
    if not ok:
        raise SystemExit("\nMissing dependencies: pip install -r requirements.txt")

    tmp = tempfile.mkdtemp(prefix="intronaute_selftest_")
    bam = os.path.join(tmp, "mini.bam")
    with open(bam, "wb") as fh:
        fh.write(base64.b64decode(MINI_BAM_B64))
    with open(bam + ".bai", "wb") as fh:
        fh.write(base64.b64decode(MINI_BAI_B64))
    tsv = os.path.join(tmp, "targets.tsv")
    with open(tsv, "w") as fh:
        fh.write(TARGETS)

    reader = bamio.BamReader(bam)
    print(f"  BAM engine  : {reader.backend}")
    print(f"  contigs     : {reader.references}")
    assert reader.resolve_chrom("chr2") == "2", "chromosome name mapping"
    reader.close()

    targets = Q.load_targets(tsv)
    out = os.path.join(tmp, "mini.tsv.gz")
    meta = Q.quantify_bam(bam, targets, out, library="unstranded", quiet=True,
                          progress_every=0)

    import pandas as pd

    df = pd.read_csv(out, sep="\t").set_index("target_id")
    fails = []
    for tid, (ee, eil, eir, nin, core) in EXPECTED.items():
        r = df.loc[tid]
        got = (int(r.n_EE), int(r.n_EI_left), int(r.n_EI_right),
               int(r.n_intron_reads), int(r.core_len))
        status = "OK  " if got == (ee, eil, eir, nin, core) else "FAIL"
        if status != "OK  ":
            fails.append((tid, got, (ee, eil, eir, nin, core)))
        print(f"  [{status}] {tid}: EE={got[0]} EIg={got[1]} EId={got[2]} "
              f"intron_reads={got[3]} core={got[4]}bp")

    r = df.loc["MI_test_chr1"]
    ratio = r.intron_depth_mean / r.exon_depth_mean
    print(f"  intron/exon coverage (chr1) = {ratio:.3f} "
          f"(intron={r.intron_depth_mean:.1f}, exons={r.exon_depth_mean:.1f})")
    if not (0.05 < ratio < 1.5):
        fails.append(("ratio", ratio, "0.05-1.5"))

    if fails:
        print("\nSELF-TEST FAILED:")
        for f in fails:
            print("   ", f)
        raise SystemExit(1)
    print(f"\nAll good. BAM reading works here ({meta['seconds']}s).")
    print(f"(temporary files: {tmp})")
