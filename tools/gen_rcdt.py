#!/usr/bin/env python3
"""gen_rcdt.py - regenerate and verify the IRS RCDT tables.

The script has three jobs:
  1. regenerate the exact 93-bit PDT/CDT/RCDT byte tables from the design
     parameters;
  2. verify the checked-in limb tables are byte-for-byte equal to those tables;
  3. verify the induced finite-precision distributions satisfy the documented
     Renyi-divergence bounds.

The generated C tables are patched into the @@AUTOGEN:rcdt@@ regions of
ref/sampler.c, avx2/sampler.c, and avx512/sampler.c.
"""

import math
import os
from decimal import Decimal, getcontext

import autogen

getcontext().prec = 300

THETA = 93
LIMB_BITS = 31
LIMBS = 3
P = 1 << THETA

# (sigma, Renyi order, minimum -log2(R_a - 1), rows)
Z_TABLES = {
    "120": (Decimal("1.421875"), 249, 80.0, [
        (0x7A81A80C, 0x1C543995, 0x47E904C9),
        (0x526DB044, 0x0EB15A73, 0x1C1C254B),
        (0x66D499C6, 0x3EEB9B72, 0x0740B175),
        (0x6A76244D, 0x28DEE8C3, 0x01323963),
        (0x722B35B5, 0x5CA79891, 0x001FAD77),
        (0x6862CF2D, 0x0122DAD4, 0x00020847),
        (0x0DA12494, 0x1F1A637A, 0x00001490),
        (0x4AC0D4AA, 0x54FA2079, 0x0000007F),
        (0x20D97E2F, 0x72907420, 0x00000001),
        (0x6BFAF40C, 0x02335EAE, 0x00000000),
        (0x61F1476A, 0x00031F01, 0x00000000),
        (0x4E2B723F, 0x000002B3, 0x00000000),
        (0x369F3BB7, 0x00000001, 0x00000000),
        (0x003AC8E8, 0x00000000, 0x00000000),
    ]),
    "180": (Decimal("1.59375"), 369, 80.0, [
        (0x3E0EF881, 0x33B014DC, 0x4CBF84DA),
        (0x1154B70D, 0x67DBF313, 0x22A77549),
        (0x0F30D703, 0x5237B27A, 0x0B553B20),
        (0x4A623BA8, 0x7C42BC28, 0x029E0304),
        (0x08D8DB0F, 0x29AEC850, 0x006B85CB),
        (0x4D5B61AA, 0x269E0E6D, 0x000BDE2C),
        (0x503F9525, 0x1FA37291, 0x0000E51F),
        (0x40FB6AC4, 0x2C8F7338, 0x00000BC1),
        (0x53D2D77E, 0x5EF14273, 0x00000068),
        (0x69A49D58, 0x3BF8203A, 0x00000002),
        (0x70E8EE5A, 0x05093725, 0x00000000),
        (0x7B4E0978, 0x000DE26C, 0x00000000),
        (0x0D770E6D, 0x000019DA, 0x00000000),
        (0x3F5868BA, 0x00000020, 0x00000000),
        (0x0DC88430, 0x00000000, 0x00000000),
        (0x0007E0E1, 0x00000000, 0x00000000),
    ]),
    "260": (Decimal("2.0078125"), 529, 80.0, [
        (0x135CE1E7, 0x0CE9F6A3, 0x5590CA61),
        (0x51107A58, 0x63FA476D, 0x3014BAD6),
        (0x78F93DEA, 0x75915F0D, 0x163E3905),
        (0x66F3987B, 0x1F70F6E1, 0x085880F0),
        (0x29576A36, 0x240B0CBD, 0x02834FEF),
        (0x3F9AD1B7, 0x0AB033B3, 0x009A4A22),
        (0x76ECDA55, 0x6751D5D7, 0x001D5157),
        (0x265706D2, 0x7A9C4D13, 0x00046592),
        (0x1BB09EFC, 0x1EDF58E5, 0x000084DA),
        (0x64C7828F, 0x6964CD5D, 0x00000C50),
        (0x0FBF9931, 0x199CB932, 0x000000E5),
        (0x62A8475D, 0x068F7A91, 0x0000000D),
        (0x7237064D, 0x4A72297E, 0x00000000),
        (0x4E9A9404, 0x029837FC, 0x00000000),
        (0x562C8A76, 0x00121871, 0x00000000),
        (0x456780E5, 0x0000629C, 0x00000000),
        (0x628E5D0A, 0x000001A3, 0x00000000),
        (0x39BBD11D, 0x00000005, 0x00000000),
        (0x071304D0, 0x00000000, 0x00000000),
        (0x000E3F77, 0x00000000, 0x00000000),
    ]),
}

NOISE_TABLE = (Decimal("1.0"), 529, 88.0, [
    (0x057FB7AB, 0x2DF4C635, 0x36FED273),
    (0x055455E2, 0x2E17C41B, 0x0AB7433D),
    (0x22A8A96C, 0x45896D14, 0x00D5F51C),
    (0x44D9083F, 0x399BDCFB, 0x000656E6),
    (0x36EE7761, 0x11C429AF, 0x000011E7),
    (0x61F254F6, 0x574A372E, 0x00000012),
    (0x7B1A38BC, 0x03977C96, 0x00000000),
    (0x76A6C503, 0x00008224, 0x00000000),
    (0x6397C3AB, 0x00000006, 0x00000000),
    (0x00109FA6, 0x00000000, 0x00000000),
])


def val(limbs):
    return limbs[0] + (limbs[1] << LIMB_BITS) + (limbs[2] << (2 * LIMB_BITS))


def limbs(v):
    mask = (1 << LIMB_BITS) - 1
    return tuple((v >> (i * LIMB_BITS)) & mask for i in range(LIMBS))


def rho(k, sigma):
    x = Decimal(k)
    return Decimal.exp(-(x * x) / (Decimal(2) * sigma * sigma))


def ideal_half_distribution(sigma, max_val):
    weights = [rho(k, sigma) for k in range(max_val)]
    total = sum(weights)
    return [w / total for w in weights]


def build_tables(sigma, w):
    """Build PDT/CDT/RCDT exactly as documented in tools/BaseSampler.ipynb.

    The implementation table drops the terminal zero RCDT entry.  Therefore a
    table with n rows corresponds to support {0,...,n} and w == n.
    """
    weights = [rho(k, sigma) for k in range(w + 1)]
    total = sum(weights)
    dist = [x / total for x in weights]

    pdt = [0] * (w + 1)
    tail_sum = 0
    for z in range(1, w + 1):
        pdt[z] = int(dist[z] * P)
        tail_sum += pdt[z]
    pdt[0] = P - tail_sum

    cdt = []
    acc = 0
    for x in pdt:
        acc += x
        cdt.append(acc)
    assert cdt[-1] == P

    rcdt = [P - x for x in cdt]
    assert rcdt[-1] == 0
    return pdt, cdt, rcdt


def pdt_from_rcdt_rows(rows):
    rcdt = [val(r) for r in rows] + [0]
    cdt = [P - x for x in rcdt]
    pdt = [cdt[0]] + [cdt[i] - cdt[i - 1] for i in range(1, len(cdt))]
    return pdt


def renyi_divergence(pdt, sigma, alpha):
    q = ideal_half_distribution(sigma, 100 * max(math.ceil(float(sigma)), 2))
    a = Decimal(alpha)
    s = Decimal(0)
    for x, px_int in enumerate(pdt):
        if px_int == 0:
            continue
        px = Decimal(px_int) / P
        s += (px ** a) / (q[x] ** (a - 1))
    return s ** (Decimal(1) / (a - 1))


def renyi_inf(pdt, sigma):
    q = ideal_half_distribution(sigma, 100 * max(math.ceil(float(sigma)), 2))
    worst = Decimal(0)
    for x, px_int in enumerate(pdt):
        if px_int == 0:
            continue
        worst = max(worst, (Decimal(px_int) / P) / q[x])
    return worst


def bits_over_one(r):
    d = r - Decimal(1)
    if d <= 0:
        return float("inf")
    return -float(d.ln() / Decimal(2).ln())


def check_table(name, sigma, alpha, min_bits, rows):
    vals = [val(r) for r in rows]
    limb_ok = all(all(0 <= x < (1 << LIMB_BITS) for x in r) for r in rows)
    mono = all(vals[i] > vals[i + 1] for i in range(len(vals) - 1))

    pdt_expected, _, rcdt_expected = build_tables(sigma, len(rows))
    expected_rows = [limbs(x) for x in rcdt_expected[:-1]]
    byte_exact = rows == expected_rows

    pdt_actual = pdt_from_rcdt_rows(rows)
    ra = renyi_divergence(pdt_actual, sigma, alpha)
    rinf = renyi_inf(pdt_actual, sigma)
    ra_bits = bits_over_one(ra)
    rinf_bits = bits_over_one(rinf)
    renyi_ok = ra_bits >= min_bits

    # Keep the old sanity check too: exact tables should approximate ideal
    # tails to within tiny absolute error; this catches support/sigma mixups.
    q = ideal_half_distribution(sigma, 100 * max(math.ceil(float(sigma)), 2))
    tails = [sum(q[i + 1:]) for i in range(len(rows))]
    worst_tail = max(abs((Decimal(vals[i]) / P) - tails[i]) for i in range(len(rows)))

    ok = limb_ok and mono and byte_exact and renyi_ok
    return {
        "ok": ok,
        "limb_ok": limb_ok,
        "mono": mono,
        "byte_exact": byte_exact,
        "renyi_ok": renyi_ok,
        "ra_bits": ra_bits,
        "rinf_bits": rinf_bits,
        "min_bits": min_bits,
        "worst_tail": worst_tail,
        "first_mismatch": next(
            (i for i, (a, b) in enumerate(zip(rows, expected_rows)) if a != b), None
        ),
    }


BODY = '''\
/* 93-bit reverse-CDT tables, 3 x 31-bit limbs (low, mid, high); entry i encodes
 * Pr[X >= i+1] * 2^93 (half-discrete-gaussian).  See tools/gen_rcdt.py.
 * The generator checks byte-exact table regeneration and Renyi bounds.
 *   RCDT_Z[RCDT_ENTRIES][3] : BLISS base sampler, sigma2 = r/128
 *   RCDT_NOISE[10][3]           : key-noise sampler, sigma = 1 (shared) */
{zblocks}

const uint32_t RCDT_NOISE[NOISE_RCDT_ENTRIES][3] = {{
{noise}
}};'''


def fmt_rows(rows, indent="  "):
    return "\n".join(indent + "{0x%08XU, 0x%08XU, 0x%08XU}," % r for r in rows)


def report(name, result):
    status = "OK" if result["ok"] else "FAIL"
    mismatch = ""
    if result["first_mismatch"] is not None:
        mismatch = f" first_mismatch={result['first_mismatch']}"
    print(
        f"  {name}: limbs={result['limb_ok']} monotone={result['mono']} "
        f"byte_exact={result['byte_exact']} "
        f"R_a=1+2^-{result['ra_bits']:.2f} "
        f"R_inf=1+2^-{result['rinf_bits']:.2f} "
        f"R_a_bound>=2^-{result['min_bits']:.0f} "
        f"max|tail-table|={float(result['worst_tail']):.2e}{mismatch} -> {status}"
    )


def main():
    print("gen_rcdt.py - transcribe + verify IRS RCDT tables")
    print("  checks: monotone limbs, byte-exact regeneration, Renyi bounds\n")
    allok = True
    zblocks = []
    for sec, (sigma, alpha, min_bits, rows) in Z_TABLES.items():
        result = check_table(f"z IRS-{sec}", sigma, alpha, min_bits, rows)
        allok &= result["ok"]
        report(f"z IRS-{sec} (sigma={sigma}, alpha={alpha}, rows={len(rows)})", result)
        zblocks.append(
            f"#if LITHIUM_MODE == {sec}\n"
            f"const uint32_t RCDT_Z[RCDT_ENTRIES][3] = {{\n"
            f"{fmt_rows(rows)}\n}};\n#endif")

    sigma, alpha, min_bits, rows = NOISE_TABLE
    result = check_table("noise", sigma, alpha, min_bits, rows)
    allok &= result["ok"]
    report(f"noise (sigma={sigma}, alpha={alpha}, rows={len(rows)})", result)

    body = BODY.format(zblocks="\n\n".join(zblocks), noise=fmt_rows(NOISE_TABLE[3]))
    here = os.path.dirname(os.path.abspath(__file__))
    for rel in (("ref", "sampler.c"), ("avx2", "sampler.c"), ("avx512", "sampler.c")):
        path = os.path.normpath(os.path.join(here, "..", *rel))
        if os.path.exists(path):
            autogen.patch_region(path, "rcdt", body)
            print(f"\nPatched @@AUTOGEN:rcdt@@ region in {path}")
    print("ALL OK" if allok else "FAILURES")
    return 0 if allok else 1


if __name__ == "__main__":
    raise SystemExit(main())
