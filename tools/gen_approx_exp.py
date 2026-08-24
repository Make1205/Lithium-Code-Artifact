#!/usr/bin/env python3
"""gen_approx_exp.py - generate + validate the IRS fixed-point exp constants and
patch ref/approx_exp.h, PRINTING every derivation (Python-reproduction
principle: every constant here is either recomputed in pure Python from
ln2/factorials, or parsed bit-exactly from the GALACTICS tool log).

TWO exp routines share the generated constants:

  approx_exp(a_q60)       -- SAMPLER hot path.  exp(-a)*2^63, a = a_q60/2^60 in
                             [0, ~8).  Scheme = the lookup-table ApproxExp scheme, table size
                             B = 2^BETA (here BETA=3 -> B=8):
                                 exp(-a) = 2^{-n} * 2^{-j/B} * exp(r_l),
                                 m = ceil(B*a/ln2)+(B-1), n=m>>BETA, j=m&(B-1),
                                 r_l = (m/B)*ln2 - a  in [(B-1)/B ln2, ln2),
                             with exp(r_l) by a Sobolev-optimal degree-D Horner in
                             Q63 and 2^{-j/B} from a B-entry constant-time table.
                             The B-entry scan dominates cost, so the (B,deg) x
                             k-way sweep (docs/design-notes/LITHIUM_x64_ApproxExp.md) picks
                             the smallest table meeting 2^-43: B=8, degree-6.

  approx_exp_wide(a_q52)  -- irs.c rejection recurrence.  exp(-a)*2^63,
                             a in [0,40].  Uses the same B=8 table range
                             reduction as approx_exp, with Q52 input and a
                             degree-7 IRS polynomial.

  exp_small_signed(d_q52) -- Q62 table-reduced exp(d), |d| < ln2/8.  IRS uses
                             0 <= d <= max(BK_SQ/r^2) ~= 0.064.

Constant provenance (so every number is reproducible, never hand-typed):
  * APPROX_EXP_RCP/LB, APPROX_EXP_TABLE_T, per-mode I80
    -> recomputed below in pure Python (decimal).  Re-running this script
       reproduces the generated block byte-for-byte.
  * EXP_POLY[] and EXP_IRS_POLY[] -> PARSED from GALACTICS logs
        tools/LogPolyApprox/expx_(B-1)divBmulln2_ln2_43_64.txt   (sampler Q63)
        tools/LogPolyApprox/expx_(B-1)divBmulln2_ln2_53_64.txt   (IRS Q63)

Run:  python3 gen_approx_exp.py
"""
import os, re, math, io
import autogen
from decimal import Decimal, getcontext

getcontext().prec = 300

# ============================================================
# Design parameters (sampler scheme = the lookup-table ApproxExp scheme, table size B = 2^BETA).
#
# BETA was chosen by the (B, degree) x k-way sweep in docs/design-notes/LITHIUM_x64_ApproxExp.md:
# the cost is dominated by the B-entry constant-time table scan, so a SMALLER
# table is faster.  At the 2^-43 target the minimal degree per B is
#   B=32->deg5, B=16->deg5, B=8->deg6, B=4->deg7.
# B=8 (BETA=3) is the sweep-optimal: degree-6, one-cache-line 64 B table,
# ~1.5x the 4-way exp throughput of the B=32 with ~6 bits of margin.
# ============================================================
BETA = 3
B = 1 << BETA                       # = 8
RCP_SHIFT = 57                      # RCP = floor((B/ln2) * 2^57)
LB_SHIFT = 64 + BETA               # LB = round((ln2/B) * 2^(64+BETA))
I_SHIFT = 80                       # per-mode I = floor(2^80 / (2 r^2))
POLY_Q = 63                        # Q63 coefficients for the minimax poly
# POLY_DEG is derived from the parsed GALACTICS log (degree-agnostic).

# ============================================================
# IRS wide/small table-reduced parameters
# ============================================================
Q52_CEIL_SHIFT = 45                # Q52 input: mulh result has 45 frac bits
Q52_ARG_SHIFT = 64 + BETA - 52     # lift Q52 input into the LB/R scale
MODES = {"120": 182, "180": 204, "260": 257}   # r per security level
BK_SQ = {"120": 2116, "180": 2632, "260": 4193}

# High-precision ln2 as a literal (same digits used by the checked reference generator
# so APPROX_EXP_RCP/LB/TABLE_T come out bit-identical across regenerations).
LN2_DEC = Decimal(
    "0.69314718055994530941723212145817656807550013436025525412068000949339"
    "362196969471560586332699641868754200148102057068573368552023575813055"
    "703267075163507596193072757082837143519030703862389167347112335011")
STEP_Q52 = int((LN2_DEC / Decimal(B) * Decimal(1 << 52)).quantize(Decimal(1)))

# Float ln2 for the numerical validation gates.
LN2 = math.log(2)

# Path to the GALACTICS log carrying the Q63 coefficients.  The interval is
# [(B-1)/B*ln2, ln2], so the log name is derived from B (run_polyapprox sanitizes
# '(B-1)/B*log(2)' to '(B-1)divBmulln2').
_HERE = os.path.dirname(os.path.abspath(__file__))
_LOG_NAME = "expx_%ddiv%dmulln2_ln2_43_64.txt" % (B - 1, B)
_IRS_LOG_NAME = "expx_%ddiv%dmulln2_ln2_53_64.txt" % (B - 1, B)
POLY_LOG = os.environ.get(
    "EXP_POLY_LOG",
    os.path.normpath(os.path.join(
        _HERE, "LogPolyApprox", _LOG_NAME)))
IRS_POLY_LOG = os.environ.get(
    "EXP_IRS_POLY_LOG",
    os.path.normpath(os.path.join(
        _HERE, "LogPolyApprox", _IRS_LOG_NAME)))


# ============================================================
# Pure-Python reproduction of every recomputable constant
# ============================================================
def fact(k):
    f = 1
    for i in range(2, k + 1):
        f *= i
    return f


def make_rcp():
    """RCP = floor((B / ln2) * 2^RCP_SHIFT). int() truncates == floor for >0."""
    return int(Decimal(B) / LN2_DEC * Decimal(1 << RCP_SHIFT))


def make_lb():
    """LB = round((ln2 / B) * 2^LB_SHIFT)."""
    return int((LN2_DEC / Decimal(B) * Decimal(1 << LB_SHIFT)).quantize(Decimal(1)))


def make_table():
    """T[j] = round(2^{-j/B} * 2^63) = round(exp(-j*ln2/B) * 2^63), j=0..B-1."""
    two63 = Decimal(1 << 63)
    tab = []
    for j in range(B):
        val = two63 * (Decimal(-j) * LN2_DEC / Decimal(B)).exp()
        tab.append(int(val.quantize(Decimal(1))))
    assert tab[0] == (1 << 63), f"T[0] expected 2^63, got {tab[0]}"
    return tab


def make_recips():
    """Per-mode I80 = floor(2^80 / (2 r^2)), split HI/LO."""
    recips = {}
    for sec, r in MODES.items():
        two_r2 = 2 * r * r
        I80 = (1 << I_SHIFT) // two_r2
        recips[sec] = (two_r2, I80, I80 >> 64, I80 & ((1 << 64) - 1))
    return recips


# ============================================================
# Parse the degree-5 Q63 minimax coefficients from the GALACTICS log
# ============================================================
def parse_poly_q63(log_path):
    """Read C0..C{POLY_DEG} (Q63) from the 'Running unified Q63' -> 'latex exact
    values poly' block of the GALACTICS log.  Bit-exact, no hand entry."""
    if not os.path.exists(log_path):
        raise SystemExit(
            f"[gen_approx_exp] missing poly log: {log_path}\n"
            "  regenerate it with:\n"
            f"    ./tools/run_polyapprox.sh --f 'exp(x)' --a '{B-1}/{B}*log(2)' "
            "--b 'log(2)' --target 43 --nbits 64 --qbits 63,62,61,60,59,58,57,56 "
            "--maxdeg 10")
    txt = open(log_path).read()
    m = re.search(r'Running unified Q63\b', txt)
    if not m:
        raise SystemExit(f"[gen_approx_exp] no 'Running unified Q63' block in {log_path}")
    seg = txt[m.end():]
    m2 = re.search(r'latex exact values poly:', seg)
    if not m2:
        raise SystemExit(f"[gen_approx_exp] no Q63 'latex exact values poly' in {log_path}")
    coeffs = {}
    # lines look like:  P &= <int> \cdot 2 ^ {-63} \\
    #                     &+ <int> \cdot 2 ^ {-63} \cdot x ^ {k}\\
    line_re = re.compile(r'(\d+)\s*\\cdot 2 \^ \{-63\}'
                         r'(?:\s*\\cdot x \^ \{(\d+)\})?')
    for line in seg[m2.end():].splitlines():
        if '====' in line or 'Running unified' in line:
            break
        mm = line_re.search(line)
        if mm:
            k = int(mm.group(2)) if mm.group(2) else 0
            coeffs[k] = int(mm.group(1))
    if not coeffs:
        raise SystemExit(f"[gen_approx_exp] Q63 poly parse found no coeffs in {log_path}")
    deg = max(coeffs)
    poly = [coeffs.get(k) for k in range(deg + 1)]
    if any(c is None for c in poly):
        raise SystemExit(f"[gen_approx_exp] Q63 poly parse incomplete: {coeffs}")
    return poly       # length deg+1; degree is len(poly)-1


# ============================================================
# Bit-exact Python mirrors of the C routines (validation)
# ============================================================
def mulh(x, y):
    return (x * y) >> 64


def approx_exp_q60(a_q60, POLY):
    """Mirror of the new C approx_exp (sampler scheme, Q60 in, Q63 out)."""
    if a_q60 < 0:
        a_q60 = 0
    V = mulh(a_q60, RCP)
    m = ((V + ((1 << 53) - 1)) >> 53) + (B - 1)
    n = m >> BETA
    j = m & (B - 1)
    R = m * LB - (a_q60 << (BETA + 4))
    X = (R >> BETA) & ((1 << 64) - 1)            # Q64 (mask to 64 bits)
    deg = len(POLY) - 1
    p = POLY[deg]
    for k in range(deg - 1, -1, -1):
        p = mulh(p, X) + POLY[k]
    s = TABLE[j]
    W = (mulh(p, s) << 1) & ((1 << 64) - 1)
    if n >= 64:
        return 0
    round_bit = (1 << n) >> 1
    return (W + round_bit) >> n


def exp_reduce_wide_q52(a_q52):
    if a_q52 < 0:
        a_q52 = 0
    V = mulh(a_q52, RCP)
    m = ((V + ((1 << Q52_CEIL_SHIFT) - 1)) >> Q52_CEIL_SHIFT) + (B - 1)
    R = m * LB - (a_q52 << Q52_ARG_SHIFT)
    low = (B - 1) * LB
    if R < low:
        m += 1
        R += LB
    n = m >> BETA
    j = m & (B - 1)
    X = (R >> BETA) & ((1 << 64) - 1)
    return X, n, j


def exp_reduce_small_signed_q52(d_q52):
    neg = 1 if d_q52 < 0 else 0
    ge_step = (neg ^ 1) & (1 if d_q52 >= STEP_Q52 else 0)
    m = (B - 1) + neg - ge_step
    R = m * LB + (d_q52 << Q52_ARG_SHIFT)
    n = m >> BETA
    j = m & (B - 1)
    X = (R >> BETA) & ((1 << 64) - 1)
    return X, n, j


def exp_irs_poly_q63(X, IRS_POLY):
    p = IRS_POLY[-1]
    for k in range(len(IRS_POLY) - 2, -1, -1):
        p = mulh(p, X) + IRS_POLY[k]
    return p


def approx_exp_wide_q52(a_q52, IRS_POLY):
    """Mirror of approx_exp_wide (table-reduced degree-7, Q52 in, Q63 out)."""
    X, n, j = exp_reduce_wide_q52(a_q52)
    p = exp_irs_poly_q63(X, IRS_POLY)
    W = (mulh(p, TABLE[j]) << 1) & ((1 << 64) - 1)
    if n >= 64:
        return 0
    round_bit = (1 << n) >> 1
    return (W + round_bit) >> n


def exp_small_signed_q62(d_q52, IRS_POLY):
    X, n, j = exp_reduce_small_signed_q52(d_q52)
    p = exp_irs_poly_q63(X, IRS_POLY)
    W = (mulh(p, TABLE[j]) << 1) & ((1 << 64) - 1)
    if n >= 64:
        W = 0
    else:
        W = (W + ((1 << n) >> 1)) >> n
    return W >> 1


# Module-level handles filled in main() (used by the mirrors above).
RCP = make_rcp()
LB = make_lb()
TABLE = make_table()


def main():
    print("gen_approx_exp.py - IRS fixed-point exp constants (reproducible)\n")

    poly = parse_poly_q63(POLY_LOG)
    irs_poly = parse_poly_q63(IRS_POLY_LOG)
    deg = len(poly) - 1
    irs_deg = len(irs_poly) - 1
    print(f"  sampler scheme: B={B} (BETA={BETA}), degree-{deg} Q63 minimax")
    print(f"    RCP = {RCP}  (bits {RCP.bit_length()})")
    print(f"    LB  = {LB}  (bits {LB.bit_length()})")
    print(f"    T[0]={TABLE[0]}  T[{B-1}]={TABLE[-1]}")
    print(f"    poly C0..C{deg} (from {os.path.basename(POLY_LOG)}):")
    for k, c in enumerate(poly):
        print(f"      C{k} = {c}  (= {c / (1 << POLY_Q):.12f}, bits {c.bit_length()})")
    print(f"\n  IRS wide/small scheme: B={B} (BETA={BETA}), degree-{irs_deg} Q63 minimax")
    print(f"    Q52_CEIL_SHIFT={Q52_CEIL_SHIFT}  Q52_ARG_SHIFT={Q52_ARG_SHIFT}  STEP_Q52={STEP_Q52}")
    print(f"    irs poly C0..C{irs_deg} (from {os.path.basename(IRS_POLY_LOG)}):")
    for k, c in enumerate(irs_poly):
        print(f"      I{k} = {c}  (= {c / (1 << POLY_Q):.12f}, bits {c.bit_length()})")

    # ---- validate the NEW sampler routine (Q60) ----
    N = 80000
    rel8_new = 0.0
    for i in range(N + 1):
        a = 8.0 * i / N
        got = approx_exp_q60(round(a * (1 << 60)), poly) / (1 << 63)
        want = math.exp(-a)
        rel8_new = max(rel8_new, abs(got - want) / want)
    print(f"\n  approx_exp (Q60) rel.err over sampler range [0,8]: {rel8_new:.3e} "
          f"(2^{math.log2(rel8_new):.1f})  [gate < 2^-44, need <= 2^-43]")
    ok = rel8_new < 2 ** -44

    # ---- validate the WIDE Taylor routine (Q52) ----
    rel8_w = 0.0
    for i in range(N + 1):
        a = 8.0 * i / N
        got = approx_exp_wide_q52(round(a * (1 << 52)), irs_poly) / (1 << 63)
        rel8_w = max(rel8_w, abs(got - math.exp(-a)) / math.exp(-a))
    abs40 = 0.0
    for i in range(N + 1):
        a = 40.0 * i / N
        got = approx_exp_wide_q52(round(a * (1 << 52)), irs_poly) / (1 << 63)
        abs40 = max(abs40, abs(got - math.exp(-a)))
    print(f"  approx_exp_wide (Q52) rel.err [0,8]: {rel8_w:.3e} (2^{math.log2(rel8_w):.1f}) "
          f"[gate < 2^-44]")
    print(f"  approx_exp_wide (Q52) abs.err [0,40]: {abs40:.3e} (2^{math.log2(abs40):.1f}) "
          f"[gate < 2^-50]")
    ok = ok and (rel8_w < 2 ** -44) and (abs40 < 2 ** -50)

    # ---- validate exp_small_signed ----
    rels = 0.0
    for i in range(N + 1):
        d = -0.05 + 0.10 * i / N
        got = exp_small_signed_q62(round(d * (1 << 52)), irs_poly) / (1 << 62)
        rels = max(rels, abs(got - math.exp(d)) / math.exp(d))
    print(f"  exp_small_signed rel.err [-0.05,0.05]: {rels:.3e} (2^{math.log2(rels):.1f}) "
          f"[gate < 2^-44]  deg={irs_deg}")
    ok = ok and (rels < 2 ** -44)

    small_max = max(BK_SQ[sec] / (MODES[sec] * MODES[sec]) for sec in MODES)
    rels_real = 0.0
    for i in range(N + 1):
        d = small_max * i / N
        got = exp_small_signed_q62(round(d * (1 << 52)), irs_poly) / (1 << 62)
        rels_real = max(rels_real, abs(got - math.exp(d)) / math.exp(d))
    print(f"  exp_small_signed rel.err actual IRS [0,{small_max:.6f}]: {rels_real:.3e} "
          f"(2^{math.log2(rels_real):.1f}) [gate < 2^-44]")
    ok = ok and (rels_real < 2 ** -44)

    recips = make_recips()
    print("\n  sampler reciprocal I80 = floor(2^80/(2 r^2))  (a_q60 = (num*I80)>>20):")
    for sec, (two_r2, I80, I_HI, I_LO) in recips.items():
        print(f"    IRS-{sec}: r={MODES[sec]} 2r^2={two_r2} I80={I80} HI={I_HI} LO={I_LO}")

    # ---- patch the @@AUTOGEN:approx_exp@@ region in approx_exp.h ----
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.normpath(os.path.join(here, "..", "ref", "approx_exp.h"))
    with io.StringIO() as f:
        # ---- sampler scheme ----
        f.write("/* ===== sampler scheme: approx_exp(a_q60) = round(exp(-a)*2^63) =====\n"
                " * exp(-a) = 2^{-n} * 2^{-j/B} * exp(r_l),  m=ceil(B*a/ln2)+(B-1),\n"
                f" * n=m>>BETA, j=m&(B-1), r_l in [(B-1)/B ln2, ln2) = [{B-1}/{B} ln2, ln2). */\n")
        f.write(f"#define APPROX_EXP_BETA   {BETA}\n")
        f.write("#define APPROX_EXP_B      (1u << APPROX_EXP_BETA)\n\n")
        f.write(f"/* RCP = floor((B/ln2) * 2^{RCP_SHIFT}); V = mulh64(a_q60, RCP) = "
                f"floor((B*a/ln2)*2^53). */\n")
        f.write(f"#define APPROX_EXP_RCP_SHIFT  {RCP_SHIFT}\n")
        f.write(f"#define APPROX_EXP_RCP        UINT64_C({RCP})\n\n")
        f.write(f"/* LB = round((ln2/B) * 2^(64+BETA)); R = m*LB - (a_q60<<(BETA+4)), "
                f"X = R>>BETA (Q64). */\n")
        f.write(f"#define APPROX_EXP_LB_SHIFT   {LB_SHIFT}\n")
        f.write(f"#define APPROX_EXP_LB         UINT64_C({LB})\n\n")
        f.write(f"/* Degree-{deg} Q63 minimax coeffs EXP_POLY[k], P(r_l) ~ exp(r_l)\n"
                " * (all positive -> unsigned Horner).  Source (bit-exact):\n"
                f" *   tools/LogPolyApprox/{_LOG_NAME}\n"
                " * (GALACTICS Sobolev fit + LLL Babai-rounding to Q63). */\n")
        f.write(f"#define EXP_POLY_DEG {deg}\n")
        f.write(f"static const uint64_t EXP_POLY[EXP_POLY_DEG + 1] = {{\n  ")
        f.write(", ".join(f"UINT64_C({c})" for c in poly))
        f.write("\n};\n\n")
        f.write("/* Q63 table T[j] = round(2^{-j/B} * 2^63), j in 0..B-1 (constant-time scan). */\n")
        f.write("static const uint64_t APPROX_EXP_TABLE_T[APPROX_EXP_B] = {\n")
        for i in range(0, B, 4):
            row = "  " + ", ".join(f"UINT64_C({v})" for v in TABLE[i:i + 4]) + ","
            f.write(row + "\n")
        f.write("};\n\n")

        # ---- IRS wide/small table-reduced scheme ----
        f.write("/* ===== IRS exp scheme: approx_exp_wide(a_q52), exp_small_signed(d_q52) =====\n"
                " * Same B=8 table reduction as the sampler, but Q52 input:\n"
                " * V = mulh64(a_q52, APPROX_EXP_RCP) = floor((B*a/ln2)*2^45);\n"
                " * X = R>>BETA is Q64.  The signed small path maps |d| < ln2/B into the\n"
                " * same residual interval with m=B-1 for d>=0 and m=B for d<0. */\n")
        f.write(f"#define APPROX_EXP_Q52_CEIL_SHIFT {Q52_CEIL_SHIFT}\n")
        f.write(f"#define APPROX_EXP_Q52_ARG_SHIFT  {Q52_ARG_SHIFT}\n")
        f.write(f"#define APPROX_EXP_STEP_Q52       UINT64_C({STEP_Q52})\n")
        f.write(f"#define EXP_IRS_POLY_DEG {irs_deg}\n")
        f.write(f"static const uint64_t EXP_IRS_POLY[EXP_IRS_POLY_DEG + 1] = {{\n  ")
        f.write(", ".join(f"UINT64_C({c})" for c in irs_poly))
        f.write("\n};\n\n")

        # ---- per-mode reciprocal ----
        f.write("/* ===== per-mode Q80 reciprocal of 2 r^2: a_q60 = (num*I80)>>20 ===== */\n")
        for sec, (two_r2, I80, I_HI, I_LO) in recips.items():
            f.write(f"#if LITHIUM_MODE == {sec}\n")
            f.write(f"/* 2 r^2 = {two_r2}, I80 = floor(2^80/2r^2) = {I80} */\n")
            f.write(f"#define SAMP_TWO_R2 UINT64_C({two_r2})\n")
            f.write(f"#define SAMP_I_HI   UINT64_C({I_HI})\n")
            f.write(f"#define SAMP_I_LO   UINT64_C({I_LO})\n")
            f.write(f"#endif\n")
        body = f.getvalue()
    autogen.patch_region(path, "approx_exp", body)

    print(f"\nPatched @@AUTOGEN:approx_exp@@ region in {path}")
    print("\n" + ("ALL OK" if ok else "PRECISION FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
