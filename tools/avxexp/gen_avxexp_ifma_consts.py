#!/usr/bin/env python3
"""
gen_avxexp_ifma_consts.py - reproducible generator + double-check for the
AVX-512 IFMA integer exp(-a) constants used by
avx512/approx_exp_avx512ifma.c.

Patches the @@AUTOGEN:approx_exp_avx512ifma@@ region in
avx512/approx_exp_avx512ifma.c with:
  - POLY_C_Q51[6] : degree-5 minimax of exp on [31/32 ln2, ln2), coefficients
                    in Q51 (true value = POLY_C_Q51[i] * 2^-51).  Provenance:
                    the same Sollya fpminimax fit as the AVX2 float coeffs
                    (run_avxexp.sh / run_polyapprox.sh, logged under
                    tools/avxexp/LogAVXApprox/), rounded to Q51 (<= 52 bits so
                    they are valid vpmadd52 operands).  Embedded below as the
                    source of truth and validated to beat 2^-43 here.
  - TAB_Q51[32]   : round(2^{-j/32} * 2^51), closed form, exactly reproducible.
  - the range-reduction scalars (ln2, 32/ln2, ln2/32, 2^51, 2^-60).

Run with no args to (re)generate the header.  Run with --check to ALSO assert
the max relative error over [0, AMAX] beats 2^-43 (the HARD IRS accuracy bound, see
docs/design-notes/AVXExp.md) and exit non-zero otherwise.  The avx512 `make
check-consts AVX512IFMAEXP=1` gate regenerates + requires git-clean + the gate.

The reconstruction below mirrors the C kernel EXACTLY (Q51 integer Horner via
the hi-only mulsh52 = 2*floor(p*X/2^52), the 32-entry 2^{-j/32} table, and the
2^{-n} scale) so the gate measures the kernel's true accuracy, not an idealized
polynomial.
"""
import sys, os, math

# --- source-of-truth Q51 coefficients (Sollya fpminimax, degree 5,
# [31/32 ln2, ln2)); true value = CQ[i] * 2^-51.  C[1] ~ 1.0022*2^51 saturates
# 52 bits (the widest vpmadd52 operand). ---
CQ = [
    2251231781277982,
    2256720101499690,
    1108248280973626,
    408721085595457,
    58973132430121,
    37125483605256,
]

B = 32
Q = 51
LN2 = math.log(2.0)
AMAX = 7.369  # sampler exponent upper bound (a in [0, AMAX])


def TQ(j):
    return int(round(math.ldexp(2.0 ** (-j / 32.0), Q)))  # round(2^{-j/32} * 2^51)


TAB = [TQ(j) for j in range(32)]


def mulsh(a, b, s):
    # Mirror the hi-only C mulsh52: 2*floor(a*b / 2^(s+1)), i.e. vpmadd52hi
    # << 1, dropping the bit-s carry of the low product (<= 1 ULP(Qs) per
    # multiply).  Keeps the gate measuring the kernel's TRUE accuracy.
    return ((a * b) >> (s + 1)) << 1


def reconstruct(a):
    """Replicate the Q51 IFMA kernel exactly: float range reduction, Q51 integer
    Horner + table, 2^{-n} scale (mirrors ifma_ref.c approx_exp_q51)."""
    m = math.ceil(B * a / LN2) + (B - 1)
    n = m >> 5
    j = m & 31
    rl = (m / B) * LN2 - a
    X = int(round(math.ldexp(rl, Q)))  # r_l in Q51
    p = CQ[5]
    for i in range(4, -1, -1):
        p = mulsh(p, X, Q) + CQ[i]
    ps = mulsh(p, TAB[j], Q)
    return math.ldexp(float(ps), -Q - n)


def max_rel_err():
    worst = 0.0
    steps = int(AMAX / 0.0002) + 1
    for s in range(steps + 1):
        a = min(s * 0.0002, AMAX)
        ref = math.exp(-a)
        got = reconstruct(a)
        if ref != 0.0:
            r = abs((got - ref) / ref)
            if r > worst:
                worst = r
    return worst


def max_n():
    """The Q63 output uses exp_q63 = ps << (12 - n), exact only for n <= 12.
    Assert the whole domain stays within that headroom."""
    mx = 0
    steps = int(AMAX / 0.0002) + 1
    for s in range(steps + 1):
        a = min(s * 0.0002, AMAX)
        m = math.ceil(B * a / LN2) + (B - 1)
        mx = max(mx, m >> 5)
    return mx


def emit(path):
    """Patch the @@AUTOGEN:approx_exp_avx512ifma@@ region in
    approx_exp_avx512ifma.c with the Q51 constants used by the IFMA kernel."""
    err = max_rel_err()
    nmx = max_n()
    L = []
    L.append("/* Q51 integer exp(-a) constants for the AVX-512 IFMA sampler exp.")
    L.append(" * POLY_C_Q51[i]: degree-5 minimax of exp on [31/32 ln2, ln2), Q51")
    L.append(" *   (true = POLY_C_Q51[i]*2^-51; Sollya fpminimax, same fit as the")
    L.append(" *   AVX2 float coeffs, logged under tools/avxexp/LogAVXApprox/).")
    L.append(" * TAB_Q51[j]: round(2^{-j/32} * 2^51), j=0..31 (closed form).")
    L.append(" * Validated: max rel error over [0,%.3f] = 2^%.2f (HARD gate 2^-43)."
             % (AMAX, math.log2(err)))
    L.append(" * Domain headroom: max n = %d (<= 12, so exp_q63 = ps<<(12-n) is exact). */"
             % nmx)
    L.append("/* range reduction (B=32): m=ceil(a*INV)+31, r_l=m*LB-a; X=r_l*2^Q. */")
    L.append("#define APPROX_EXP_IFMA_Q 51")
    L.append("#define APPROX_EXP_IFMA_LN2  %s   /* hex %s */" % (repr(LN2), LN2.hex()))
    L.append("#define APPROX_EXP_IFMA_INV  %s   /* hex %s  (32/ln2) */"
             % (repr(32.0 / LN2), (32.0 / LN2).hex()))
    L.append("#define APPROX_EXP_IFMA_LB   %s   /* hex %s  (ln2/32) */"
             % (repr(LN2 / 32.0), (LN2 / 32.0).hex()))
    L.append("#define APPROX_EXP_IFMA_SC   9007199254740992.0   /* 2^53 unused */")
    L.append("#define APPROX_EXP_IFMA_SCQ  2251799813685248.0   /* 2^51 */")
    L.append("#define APPROX_EXP_IFMA_INV60 8.673617379884035e-19 /* 2^-60 */")
    L.append("")
    L.append("static const uint64_t POLY_C_Q51[6] = {")
    for i in range(0, 6, 2):
        L.append("    %16dULL, %16dULL," % (CQ[i], CQ[i + 1]))
    L.append("};")
    L.append("")
    L.append("static const uint64_t TAB_Q51[32] = { /* round(2^{-j/32} * 2^51) */")
    for j in range(0, 32, 2):
        L.append("    %16dULL, %16dULL," % (TAB[j], TAB[j + 1]))
    L.append("};")
    body = "\n".join(L) + "\n"
    irs_tools = os.path.normpath(os.path.join(os.path.dirname(path), "..", "tools"))
    sys.path.insert(0, irs_tools)
    import autogen
    autogen.patch_region(path, "approx_exp_avx512ifma", body)
    return err, nmx


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.normpath(os.path.join(
        here, "..", "..", "avx512", "approx_exp_avx512ifma.c"))
    err, nmx = emit(out)
    bits = math.log2(err)
    print("gen_avxexp_ifma_consts: wrote %s" % out)
    print("  max rel error over [0,%.3f] = 2^%.2f (gate 2^-43); max n = %d (<=12)"
          % (AMAX, bits, nmx))
    if "--check" in sys.argv:
        if bits > -43.0:
            print("  FAIL: accuracy 2^%.2f does not beat 2^-43" % bits)
            sys.exit(1)
        if nmx > 12:
            print("  FAIL: max n = %d > 12 (Q63 left-shift headroom exceeded)" % nmx)
            sys.exit(1)
        print("  PASS: 2^-43 accuracy gate + n<=12 headroom met")
    return 0


if __name__ == "__main__":
    sys.exit(main())
