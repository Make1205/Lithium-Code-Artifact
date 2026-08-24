#!/usr/bin/env python3
"""rans.py - byte-exact reference rANS engine for IRS (two-table single stream).

GOLDEN reference the C irs_rans.c must match byte-for-byte.

Engine: 32-bit state x, L = 2^23 renorm lower bound, 8-bit renormalization,
prob_bits = 10 (each frequency table sums to 2^10 = 1024).  The merged stream
encodes the z-head sub-stream (table g1) then the hint sub-stream (table g2);
internally symbols are pushed in reverse (hint first, then z-head) so decode
pulls z-head then hint.

A "table" is dict(freqs, cdf, sym_lo, n) over a contiguous signed alphabet
[sym_lo, sym_lo+n-1]; cdf has n+1 prefix sums (cdf[-1] = 1024).

Run:  python3 rans.py     # self-test (2-table single-stream roundtrip)
"""
RANS_L = 1 << 23
PROB_BITS = 10
PROB_SCALE = 1 << PROB_BITS
RANS_N = 2  # interleaved streams (mirror: irs_rans.c RANS_N)


def make_table(freqs, sym_lo):
    cdf = [0]
    for f in freqs:
        cdf.append(cdf[-1] + f)
    assert cdf[-1] == PROB_SCALE, f"freqs sum {cdf[-1]} != {PROB_SCALE}"
    return dict(freqs=freqs, cdf=cdf, sym_lo=sym_lo, n=len(freqs))


def encode(streams):
    """streams: list of (table, [symbols]); decode pulls them in the SAME order.
    Flat sequence S[t] = concatenation of the streams' symbols; symbol t is
    coded on rANS state t % RANS_N, pushed in reverse t (decode mirrors it)."""
    flat = [(table, s) for table, syms in streams for s in syms]
    x = [RANS_L] * RANS_N
    out = bytearray()
    for t in reversed(range(len(flat))):
        s = t % RANS_N
        table, sym = flat[t]
        slot = sym - table["sym_lo"]
        freq = table["freqs"][slot]; start = table["cdf"][slot]
        x_max = ((RANS_L >> PROB_BITS) << 8) * freq
        while x[s] >= x_max:
            out.append(x[s] & 0xff); x[s] >>= 8
        x[s] = ((x[s] // freq) << PROB_BITS) + (x[s] % freq) + start
    for s in range(RANS_N):  # flush state 0 first (-> last in the stream)
        for _ in range(4):
            out.append(x[s] & 0xff); x[s] >>= 8
    out.reverse()
    return bytes(out)


def _slot(table, val):
    cdf = table["cdf"]
    for s in range(table["n"]):
        if cdf[s] <= val < cdf[s + 1]:
            return s
    raise ValueError("bad rANS slot")


def decode(data, plan):
    """plan: list of (table, count) in stream order. Returns list of sym-lists."""
    flat_tab = [table for table, count in plan for _ in range(count)]
    T = len(flat_tab)
    if len(data) < 4 * RANS_N:
        raise ValueError("truncated rANS stream")
    x = [0] * RANS_N
    bp = 0
    for s in range(RANS_N - 1, -1, -1):  # init: stream front is state N-1
        v = 0
        for _ in range(4):
            v = (v << 8) | data[bp]; bp += 1
        x[s] = v
    syms = [0] * T
    for t in range(T):
        s = t % RANS_N
        table = flat_tab[t]
        val = x[s] & (PROB_SCALE - 1)
        slot = _slot(table, val)
        syms[t] = table["sym_lo"] + slot
        x[s] = table["freqs"][slot] * (x[s] >> PROB_BITS) + val - table["cdf"][slot]
        while x[s] < RANS_L:
            if bp >= len(data):
                raise ValueError("truncated rANS stream")
            x[s] = (x[s] << 8) | data[bp]; bp += 1
    if bp != len(data):
        raise ValueError("trailing bytes in rANS stream")
    if any(v != RANS_L for v in x):
        raise ValueError("non-canonical final rANS state")
    out, off = [], 0
    for table, count in plan:
        out.append(syms[off:off + count]); off += count
    return out


def _mk(n, lo):
    import math
    raw = [max(1, int(800 * math.exp(-((i - n // 2) ** 2) / (2 * (n / 6) ** 2)))) for i in range(n)]
    s = sum(raw)
    raw = [max(1, round(r * PROB_SCALE / s)) for r in raw]
    raw[raw.index(max(raw))] += PROB_SCALE - sum(raw)
    return make_table(raw, lo)


def _selftest():
    import random
    rng = random.Random(1234)
    t1 = _mk(33, -16); t2 = _mk(13, -6)
    for trial in range(300):
        a = [rng.randint(-16, 16) for _ in range(rng.randint(1, 60))]
        b = [rng.randint(-6, 6) for _ in range(rng.randint(1, 40))]
        enc = encode([(t1, a), (t2, b)])
        da, db = decode(enc, [(t1, len(a)), (t2, len(b))])
        assert da == a and db == b, (trial, a, da)
        for bad in (enc + b"\x42", enc[:-1]):
            try:
                decode(bad, [(t1, len(a)), (t2, len(b))])
            except ValueError:
                pass
            else:
                raise AssertionError("non-canonical stream accepted")
    print("rans.py self-test: roundtrip OK (300 trials, 2-table single stream)")
    return 0


if __name__ == "__main__":
    raise SystemExit(_selftest())
