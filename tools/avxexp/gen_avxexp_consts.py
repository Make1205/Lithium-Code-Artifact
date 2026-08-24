#!/usr/bin/env python3
"""
gen_avxexp_consts.py - reproducible generator + double-check for the AVX2 float
exp(-a) constants used by avx2/approx_exp_avx2.c.

Patches the @@AUTOGEN:approx_exp_avx2@@ region in avx2/approx_exp_avx2.h
with:
  - APPROX_EXP_AVX2_C[7] : degree-6 binary64 minimax of exp on [7/8 ln2, ln2)
                            (provenance: Sollya fpminimax via run_avxexp.sh, logged
                             under tools/avxexp/LogAVXApprox/; embedded below as the
                             source of truth and validated to beat 2^-43 here).
  - APPROX_EXP_AVX2_T[8]: TD[j] = 2^{-j/8}, closed form, exactly reproducible.
  - the derived range-reduction scalars (LN2, 8/LN2, LN2/8, 2^63), each emitted
    as the exact double of its closed form (no magic numbers).

Run with no args to (re)generate the header.  Run with --check to ALSO assert the
max relative error over [0, AMAX] beats 2^-43 and exit non-zero otherwise.  The
`make check-consts` gate regenerates and requires the header to be git-clean
(reproducible) and the accuracy gate to pass.

Constants are emitted via float.hex() in a comment plus the shortest round-trip
decimal (repr), so the C literal and this mirror are bit-for-bit identical.
"""
import sys, os, math

# --- source-of-truth coefficients (Sollya fpminimax, degree DEG, [(B-1)/B ln2,
# ln2)) ---  C[i] is the coefficient of r^i in P(r) ~ exp(r).  B=8 / deg-6:
# the AVX2 table lookup is a constant-time scan (no cheap 32-double permute on
# AVX2, unlike AVX-512's vpermt2q), so a smaller B (shorter scan) + one more
# Horner term is the throughput optimum -- same B=8/deg-6 choice as the scalar
# x64 exp.  Regenerate the coeffs with:
#   LO='7/8*log(2)' HI='log(2)' DMIN=6 DMAX=6 TARGET=43 ./tools/avxexp/run_avxexp.sh
CD = [
    1.000017029688622,    # hex 0x3ff00011db5f4a84
    0.9998181124150078,   # hex 0x3feffe828dd8ab4c
    0.5008292820587139,   # hex 0x3fe006cb216a4862
    0.16457890339420803,  # hex 0x3fc510ebe7d84701
    0.044787736908069656, # hex 0x3fa6ee6b129063ca
    0.00559066059537778,  # hex 0x3f76e63b86b91249
    0.0026598531260834946,# hex 0x3f65ca1dc609c3c0
]

B = 8
LOG2B = 3                           # log2(B)
DEG = len(CD) - 1                   # polynomial degree (6)
LN2 = math.log(2.0)                 # nearest double to ln 2
AMAX = 7.369                        # sampler exponent upper bound (a in [0, AMAX])

def TD(j):
    return 2.0 ** (-j / B)          # exact closed form, rounded to double

def reconstruct(a):
    """Replicate the kernel in double: exp(-a) via B range reduction + Horner."""
    m = math.ceil(B * a / LN2) + (B - 1)
    n = m >> LOG2B
    j = m & (B - 1)
    rl = (m / B) * LN2 - a
    p = CD[DEG]
    for i in range(DEG - 1, -1, -1):
        p = p * rl + CD[i]
    return p * TD(j) * (2.0 ** (-n))

def max_rel_err():
    worst = 0.0
    a = 0.0
    # dense sweep incl. a=0 and a=AMAX
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

def cfloat(x):
    """Round-trip C double literal."""
    return repr(float(x))

def emit(path):
    """Patch the @@AUTOGEN:approx_exp_avx2@@ region in approx_exp_avx2.h with the
    float-exp constants (the header also hosts the hand-written kernel
    prototypes, so we patch a region rather than own the whole file)."""
    err = max_rel_err()
    lines = []
    lines.append("/* Float exp(-a) constants (binary64) for the AVX2 sampler exp.")
    lines.append(" * C[i]: degree-%d minimax of exp on [%d/%d ln2, ln2) (Sollya fpminimax,"
                 % (DEG, B - 1, B))
    lines.append(" *       logged under tools/avxexp/LogAVXApprox/).")
    lines.append(" * T[j]: 2^{-j/%d}, j=0..%d (closed form)." % (B, B - 1))
    lines.append(" * Validated: max rel error over [0,%.3f] = 2^%.2f (gate 2^-43). */" % (AMAX, math.log2(err)))
    lines.append("/* range reduction: B-entry table, degree-DEG Horner.  m=ceil(a*INV)+(B-1),")
    lines.append(" * n=m>>LOG2B, j=m&(B-1), r_l=m*LB-a; exp(-a)=2^-n * 2^{-j/B} * P(r_l). */")
    lines.append("#define APPROX_EXP_AVX2_B       %d" % B)
    lines.append("#define APPROX_EXP_AVX2_LOG2B   %d" % LOG2B)
    lines.append("#define APPROX_EXP_AVX2_DEG     %d" % DEG)
    lines.append("#define APPROX_EXP_AVX2_LN2     %s   /* hex %s */" % (cfloat(LN2), float(LN2).hex()))
    lines.append("#define APPROX_EXP_AVX2_INV     %s   /* hex %s  (B/ln2) */" % (cfloat(B/LN2), (B/LN2).hex()))
    lines.append("#define APPROX_EXP_AVX2_LB      %s   /* hex %s  (ln2/B) */" % (cfloat(LN2/B), (LN2/B).hex()))
    lines.append("#define APPROX_EXP_AVX2_SCALE63 9223372036854775808.0       /* 2^63 */")
    lines.append("")
    lines.append("static const double APPROX_EXP_AVX2_C[%d] = {" % (DEG + 1))
    for c in CD:
        lines.append("  %-22s /* hex %s */," % (cfloat(c), float(c).hex()))
    lines.append("};")
    lines.append("")
    lines.append("static const double APPROX_EXP_AVX2_T[%d] = {" % B)
    for j in range(B):
        lines.append("  %-22s /* 2^-%d/%d  hex %s */," % (cfloat(TD(j)), j, B, TD(j).hex()))
    lines.append("};")
    body = "\n".join(lines) + "\n"
    # Patch the region in approx_exp_avx2.h + clang-format (team autogen helper).
    irs_tools = os.path.normpath(os.path.join(os.path.dirname(path), "..", "tools"))
    sys.path.insert(0, irs_tools)
    import autogen
    autogen.patch_region(path, "approx_exp_avx2", body)
    return err

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    out = os.path.normpath(os.path.join(here, "..", "..", "avx2", "approx_exp_avx2.h"))
    err = emit(out)
    bits = math.log2(err)
    print("gen_avxexp_consts: wrote %s" % out)
    print("  max rel error over [0,%.3f] = 2^%.2f (gate 2^-43)" % (AMAX, bits))
    if "--check" in sys.argv:
        if bits > -43.0:
            print("  FAIL: accuracy 2^%.2f does not beat 2^-43" % bits)
            sys.exit(1)
        print("  PASS: accuracy gate met")
    return 0

if __name__ == "__main__":
    sys.exit(main())
