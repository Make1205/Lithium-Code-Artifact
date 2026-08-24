#!/usr/bin/env python3
"""check_rans.py - reproduce the test_rans.c deterministic dump byte-for-byte
using the rans.py reference, proving the C encoder matches the Python golden
model.  Parses the tables out of ref/irs_rans.h.

Run:  python3 check_rans.py        # prints expected len+hash per mode
"""
import os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rans
import autogen

NZ, NH = 64, 48


def xorshift_stream(seed):
    st = seed & ((1 << 64) - 1)
    M = (1 << 64) - 1
    while True:
        st ^= (st << 13) & M
        st ^= st >> 7
        st ^= (st << 17) & M
        yield (st >> 27) & 0xffffffff


def parse_table(blk, name):
    fr = [int(x) for x in re.search(r"%s_FREQ\[\d+\] = \{([^}]*)\}" % name, blk).group(1).split(",")]
    lo = int(re.search(r"%s_LO\s+\((-?\d+)\)" % name, blk).group(1))
    return rans.make_table(fr, lo)


def parse_mode(txt, sec):
    blk = re.search(r"#if LITHIUM_MODE == %s\b(.*?)#endif" % sec, txt, re.S).group(1)
    return parse_table(blk, "RANS_ZH"), parse_table(blk, "RANS_HINT")


def fnv1a(b):
    h = 1469598103934665603
    M = (1 << 64) - 1
    for byte in b:
        h = ((h ^ byte) * 1099511628211) & M
    return h


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    txt = autogen.read_region(os.path.normpath(os.path.join(here, "..", "ref", "irs_rans.h")),
                              "rans_tables")
    for sec in ("120", "180", "260"):
        zt, ht = parse_mode(txt, sec)
        g = xorshift_stream(0xCAFEBABEDEADBEEF)
        zh = [zt["sym_lo"] + (next(g) % zt["n"]) for _ in range(NZ)]
        hi = [ht["sym_lo"] + (next(g) % ht["n"]) for _ in range(NH)]
        enc = rans.encode([(zt, zh), (ht, hi)])
        print(f"mode={sec} DUMP len={len(enc)} hash={fnv1a(enc)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
