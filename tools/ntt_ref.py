#!/usr/bin/env python3
"""ntt_ref.py - bit-exact reference for IRS ntt.c (forward / inverse / basemul),
uint16 Montgomery, R = 2^16.

Mirrors ntt.c EXACTLY (same butterfly order, same montgomery_reduce16, same
addm16/subm16, same zetas table) so that C intermediate arrays can be compared
element-by-element.  Also checks the algebraic identities against an independent
schoolbook negacyclic multiply.

Reads the generated zetas straight out of ref/ntt.h so the reference
and the C build share one source of truth.

Run:  python3 ntt_ref.py
"""
import os, re, random

R = 1 << 16
MASK16 = R - 1

MODES = {
    "120": dict(q=57601, n=256, base_deg=2),
    "180": dict(q=30977, n=256, base_deg=2),
    "260": dict(q=64513, n=512, base_deg=1),
}


def parse_constants(path):
    """Pull per-mode {MONT,QINV,INVNTT_F,ZETAS_SIZE,zetas[]} from ntt.h."""
    txt = open(path).read()
    out = {}
    for sec in ("120", "180", "260"):
        m = re.search(r"#if LITHIUM_MODE == %s\b(.*?)#endif /\* LITHIUM_MODE == %s" % (sec, sec),
                      txt, re.S)
        if not m:
            continue
        blk = m.group(1)
        def gi(name):
            mm = re.search(
                r"#\s*define\s+LITHIUM_%s\s+\(?(-?\d+)\s*U?\)?" % name, blk
            )
            if not mm:
                raise ValueError(f"missing LITHIUM_{name} in mode {sec}")
            return int(mm.group(1))
        mont = gi("MONT"); qinv = gi("QINV"); invf = gi("INVNTT_F"); zs = gi("NTT_ZETAS_SIZE")
        zm = re.search(r"zetas\[NTT_ZETAS_SIZE\] = \{(.*?)\};", blk, re.S)
        zetas = [int(x) for x in re.findall(r"-?\d+", zm.group(1))]
        out[sec] = dict(MONT=mont, QINV=qinv, INVNTT_F=invf, ZETAS_SIZE=zs, zetas=zetas)
    return out


def make_mreduce(q, qinv):
    def mreduce(a):                          # 0 <= a < q*R
        m = (a * qinv) & MASK16
        t = (a + m * q) >> 16
        return t - q if t >= q else t
    return mreduce


def addm(a, b, q):
    s = a + b
    return s - q if s >= q else s


def subm(a, b, q):
    return a - b if a >= b else a + q - b


def ntt(a, q, zetas, base_deg, mreduce):
    a = a[:]
    n = len(a); k = 0; length = n // 2
    while length >= base_deg:
        start = 0
        while start < n:
            k += 1; zeta = zetas[k]
            for j in range(start, start + length):
                t = mreduce(zeta * a[j + length])
                a[j + length] = subm(a[j], t, q)
                a[j] = addm(a[j], t, q)
            start = j + length + 1
        length >>= 1
    return a


def invntt_tomont(a, q, zetas, base_deg, zsize, invf, mreduce):
    a = a[:]
    n = len(a); k = zsize; length = base_deg
    while length < n:
        start = 0
        while start < n:
            k -= 1; zeta = (q - zetas[k]) % q          # -zeta mod q
            for j in range(start, start + length):
                t = a[j]
                a[j] = addm(t, a[j + length], q)
                a[j + length] = mreduce(zeta * subm(t, a[j + length], q))
            start = j + length + 1
        length <<= 1
    return [mreduce(invf * x) for x in a]


def basemul_full(ah, bh, q, mreduce):
    return [mreduce(ah[i] * bh[i]) for i in range(len(ah))]


def basemul_incomplete(ah, bh, q, zetas, zsize, mreduce):
    n = len(ah); r = [0] * n; half = zsize // 2
    def bm(o, a, b, zeta):
        t = mreduce(a[1] * b[1]); t = mreduce(t * zeta); t = addm(t, mreduce(a[0] * b[0]), q); r[o] = t
        t = mreduce(a[0] * b[1]); t = addm(t, mreduce(a[1] * b[0]), q); r[o + 1] = t
    for i in range(n // 4):
        zeta = zetas[half + i]
        bm(4 * i,     ah[4*i:4*i+2],   bh[4*i:4*i+2],   zeta)
        bm(4 * i + 2, ah[4*i+2:4*i+4], bh[4*i+2:4*i+4], (q - zeta) % q)
    return r


def schoolbook(a, b, q):
    n = len(a); t = [0] * (2 * n)
    for i in range(n):
        for j in range(n):
            t[i + j] += a[i] * b[j]
    return [(t[i] - t[i + n]) % q for i in range(n)]


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    consts = parse_constants(os.path.normpath(os.path.join(here, "..", "ref", "ntt.h")))
    print("ntt_ref.py - IRS uint16 NTT reference (R=2^16) + identity check\n")
    allok = True
    for sec, p in MODES.items():
        q, n, bd = p["q"], p["n"], p["base_deg"]
        c = consts[sec]
        mont, qinv, zsize, invf, zetas = c["MONT"], c["QINV"], c["ZETAS_SIZE"], c["INVNTT_F"], c["zetas"]
        mr = make_mreduce(q, qinv)
        rng = random.Random(0x1234 ^ q)
        bad_rt = bad_mul = 0
        for _ in range(200):
            a = [rng.randrange(q) for _ in range(n)]          # [0,q)
            b = [rng.randrange(q) for _ in range(n)]
            # roundtrip: invntt_tomont(ntt(a)) == a*R ; strip one R via mreduce
            ah = ntt(a, q, zetas, bd, mr)
            rt = invntt_tomont(ah, q, zetas, bd, zsize, invf, mr)
            if any(mr(x) != a[i] for i, x in enumerate(rt)):
                bad_rt += 1
            # multiply
            bh = ntt(b, q, zetas, bd, mr)
            ch = basemul_full(ah, bh, q, mr) if bd == 1 else basemul_incomplete(ah, bh, q, zetas, zsize, mr)
            prod = invntt_tomont(ch, q, zetas, bd, zsize, invf, mr)
            ref = schoolbook(a, b, q)
            if any(prod[i] != ref[i] for i in range(n)):
                bad_mul += 1
        ok = bad_rt == 0 and bad_mul == 0
        allok &= ok
        print(f"--- IRS-{sec} (q={q}, n={n}, base_deg={bd}, zsize={zsize}) ---")
        print(f"  roundtrip bad={bad_rt}  multiply-vs-schoolbook bad={bad_mul}  "
              f"-> {'OK' if ok else 'FAIL'}")
    print("\n" + ("ALL OK" if allok else "FAILURES"))
    return 0 if allok else 1


if __name__ == "__main__":
    raise SystemExit(main())
