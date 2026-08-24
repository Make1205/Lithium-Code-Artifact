#!/usr/bin/env python3
"""irs_ref.py - executable spec + golden reference for the IRS exp-recurrence
rejection step (the IRS rejection-sampling design, three-region (f,g); Lithium.md Thm 1).

Patches ref/irs_rans.h with REJ_T and per-level M_Q62, then validates:
  * --grid : prints a FNV-1a hash of decide over a fixed (nrm,ip,r63) grid,
             which test_irs.c reproduces bit-for-bit (C<->Python lock);
  * default: the fixed-point three-region threshold matches the real-arithmetic
             threshold to within the Renyi budget at the true integer keygen
             lower gate BK_LOW_SQ <= ||s||^2, over the FULL proof window
             |u| <= U_CAP = 8 (dense window |u| <= 4 gamma via
             threshold_deviation + deep window 4 gamma < |u| <= U_CAP via
             deep_deviation with a cancellation-free small-side reference).
             U_CAP = 8 matches the validated ApproxExp domain a in [0,8]; the
             paper's capped-step hybrid (Lithium thm:irs-truncation) makes
             |u| > U_CAP a separately-budgeted tail event.

Constant-time fixed-point formulation (mirrored bit-for-bit in irs.c) --
TABLE + HORNER (see docs/design-notes/LITHIUM_x64_ApproxExp.md):
  S(a,gamma) = sum_{k>=0} (-1)^k e^{-ka-k^2 gamma} = sum_k c_k w^k, w = e^{-a},
  c_k = (-1)^k e^{-k^2 gamma}.  c_k depends ONLY on gamma (the per-signature
  norm ||s||), so precompute_c builds the c[] table once per signature; both the
  good series Sp = S(a,gamma) = Horner(c, e^{-|u|}) and the bad series
  Sb = S(-a_c,gamma) = Horner(c, e^{a_c}), a_c=min(|u|,2gamma), reuse that c[].
  sgn=(u<0); S^+ = sgn?Sb:Sp, S^- = sgn?Sp:Sb;
  thres_a=(u>=2g)?S^+:1-S^-; flag_minus = (r_rand*M < thres_a) -> z-v.
  This replaces the old forward term recurrence (5 muls/term, 3 exp/call) with
  two Horners (2 muls/term) + a per-signature table; per step only e^{-|u|}
  (wide) and e^{a_c} (small) are evaluated -- the e^{-2gamma}/e^{-gamma} bases
  fold into the once-per-signature table.

Truncation depth T: target-B params have alpha~4 (gamma~0.03); the alternating
series converges per term-pair as e^{-gamma}.  With T=33 the worst two-sided
relative kernel deviation on the sensitive window is 2^-43.28/2^-43.29/2^-43.83
per level (deep window ~2^-49), meeting the 2^-40.0 Renyi gate with >3.2 bits
of margin; T=32 gives only ~2^-40.6.  We use T=33: the minimal depth with the
intended working margin.

Run:  python3 irs_ref.py          # validate + emit
      python3 irs_ref.py --grid   # print C<->Python grid hash
"""
import os, sys, math
import autogen
import gen_approx_exp as approx_gen

T = 33
LN2 = math.log(2)
ONE62 = 1 << 62
IRS_POLY = approx_gen.parse_poly_q63(approx_gen.IRS_POLY_LOG)

# per level: (r, B_k, B_k')  full precision (param.md)
LEVELS = {
    "120": (182, 46.00893646172169, 44.5007315301564),
    "180": (204, 51.3038533806998,  49.889031744688346),
    "260": (257, 64.75410246578265, 63.245897534217356),
}

# Exact integer gates emitted by gen_bounds.py.  KeyGen accepts
# BK_LOW_SQ <= ||s||^2 <= BK_SQ, so these are the true gamma endpoints.
BK_LOW_SQ = {
    "120": 1981,
    "180": 2489,
    "260": 4001,
}
BK_SQ = {  # ceil(B_k^2) - 1
    "120": 2116,
    "180": 2632,
    "260": 4193,
}


# ---- fixed-point exp mirrors (identical to approx_exp.c) ----
def approx_exp_q52(a_q52):
    return approx_gen.approx_exp_wide_q52(a_q52, IRS_POLY)


def exp_small_signed_q62(d_q52):
    return approx_gen.exp_small_signed_q62(d_q52, IRS_POLY)


def mul62(x, y):
    return (x * y) >> 62


def Mrep(alpha):
    a = alpha
    num = 2 * a * math.sqrt(2 * math.pi) * math.exp(-(math.pi ** 2) * a * a / 2)
    den = math.exp(-1.0 / (2 * a * a)) * (1 - math.exp(-2 * (math.pi ** 2) * a * a))
    return 1.0 + num / den


def precompute_c(gamma_q52):
    """Per-signature coefficient table c_k = (-1)^k e^{-k^2 gamma} (Q62).
    Bases: Q=e^{-2gamma}, g1=e^{-gamma}; then e^{-k^2 gamma} via the ratio
    e^{-(2k-1)gamma} = g1 * Q^{k-1}.  Mirrors C precompute_c bit-for-bit."""
    Q = approx_exp_q52(2 * gamma_q52) >> 1       # e^{-2gamma} Q62
    g1 = approx_exp_q52(gamma_q52) >> 1          # e^{-gamma}  Q62
    c = [0] * T
    c[0] = ONE62
    E = ONE62; r = g1
    for k in range(1, T):
        E = mul62(E, r)                          # E = e^{-k^2 gamma}
        c[k] = -E if (k & 1) else E
        r = mul62(r, Q)                          # r = e^{-(2k+1)gamma}
    return c


def horner(c, w):
    """Sum_{k=0}^{T-1} c_k w^k via Horner (Q62, signed)."""
    acc = c[T - 1]
    for k in range(T - 2, -1, -1):
        acc = c[k] + mul62(acc, w)
    return acc


def step_fixed(ip, r, gamma_q52, M_q62, r_rand_63, c=None):
    """Bit-exact mirror of C decide internals (table + Horner). Returns
    (flag_minus, accepted, thres_a, thres_ab).  c may be precomputed (per
    signature); if None it is built from gamma_q52 (identical result)."""
    if c is None:
        c = precompute_c(gamma_q52)
    u_q52 = (ip << 52) // (r * r) if ip >= 0 else -((-ip << 52) // (r * r))
    a_q52 = abs(u_q52)
    sgn = 1 if u_q52 < 0 else 0
    two_g = 2 * gamma_q52
    a_c = min(a_q52, two_g)
    w = approx_exp_q52(a_q52) >> 1               # e^{-|u|}  Q62 (wide)
    wp = exp_small_signed_q62(a_c)               # e^{a_c}   Q62 (small, a_c>=0)
    Sp = horner(c, w)                            # S(a,gamma)
    Sb = horner(c, wp)                           # S(-a_c,gamma)
    S_pos = Sb if sgn else Sp
    S_neg = Sp if sgn else Sb
    is_ge = u_q52 >= two_g
    is_lt = u_q52 < -two_g
    is_II = (not is_ge) and (not is_lt)
    thres_a = S_pos if is_ge else (ONE62 - S_neg)
    thres_ab = (2 * ONE62 - S_pos - S_neg) if is_II else ONE62
    rM = (r_rand_63 * M_q62) >> 63
    flag_minus = 1 if rM < thres_a else 0
    accepted = 1 if rM < thres_ab else 0
    return flag_minus, accepted, thres_a, thres_ab


def fg_of(ip, r, gamma_q52, M_q62):
    """(f, g) = (a, b) of the rejection step (Q62 thresholds / M)."""
    _, _, thres_a, thres_ab = step_fixed(ip, r, gamma_q52, M_q62, 0)
    f = thres_a / M_q62
    g = thres_ab / M_q62 - f
    return f, g


def output_dist_deviation(r, v, gamma_q52, M_q62):
    """CONDITIONAL output distribution (given non-abort) vs ideal D_{Z,r}.
    Emitted mass P(z) = pi(z+v) f(z+v) + pi(z-v) g(z-v); by Lemma 1 with EXACT
    (f,g) this is pi(z)/M, so the normalised conditional output is exactly
    D_{Z,r}.  Residual = truncation + fixed-point error.  Returns (TV, max_rel)."""
    R = 12 * r
    Z = sum(math.exp(-(x * x) / (2 * r * r)) for x in range(-R, R + 1))
    def pi(x): return math.exp(-(x * x) / (2 * r * r)) / Z
    Pz = {}
    acc = 0.0
    for z in range(-R + v, R - v + 1):
        f_zp = fg_of((z + v) * v, r, gamma_q52, M_q62)[0]
        g_zm = fg_of((z - v) * v, r, gamma_q52, M_q62)[1]
        m = pi(z + v) * f_zp + pi(z - v) * g_zm
        Pz[z] = m
        acc += m
    tv = 0.0; mrel = 0.0
    for z, m in Pz.items():
        cond = m / acc
        tv += abs(cond - pi(z))
        if pi(z) > 1e-12:
            mrel = max(mrel, abs(cond - pi(z)) / pi(z))
    return tv / 2.0, mrel


def series_real(x, gamma):
    """High-precision floating reference for S(x,gamma)."""
    s = 0.0
    for k in range(1000):
        term = (-1.0 if (k & 1) else 1.0) * math.exp(-k * x - k * k * gamma)
        s += term
        if abs(term) < 1e-25:
            break
    return s


def exact_threshold(ip, r, gamma):
    """Exact real-arithmetic threshold for the same three-region decision."""
    u = ip / (r * r)
    sgn = u < 0.0
    a = abs(u)
    a_c = min(a, 2.0 * gamma)
    Sp = series_real(a, gamma)
    Sb = series_real(-a_c, gamma)
    S_pos = Sb if sgn else Sp
    S_neg = Sp if sgn else Sb
    return S_pos if u >= 2.0 * gamma else 1.0 - S_neg


def threshold_deviation(r, nrm_low, gamma_q52, M_q62):
    """Worst fixed-point threshold error at the true keygen lower gate.

    The old check used v=round(B_k') and gamma=v^2/(2r^2), which can sit above
    the actual accepted lower gate.  Here gamma is fixed from BK_LOW_SQ.
    The threshold is most sensitive around ip ~= ||s||^2, so we scan a dense
    symmetric window around that point.  The window also covers all three
    regions of the decision rule at this gamma.
    """
    gamma = nrm_low / (2.0 * r * r)
    c = precompute_c(gamma_q52)
    max_abs = 0.0
    max_rel = 0.0
    worst_ip = 0
    for ip in range(-2 * nrm_low, 2 * nrm_low + 1):
        fixed = step_fixed(ip, r, gamma_q52, M_q62, 0, c)[2] / M_q62
        exact = exact_threshold(ip, r, gamma)
        err = abs(fixed - exact)
        max_abs = max(max_abs, err)
        den = min(exact, 1.0 - exact)
        if den > 1e-12:
            rel = err / den
            if rel > max_rel:
                max_rel = rel
                worst_ip = ip
    return max_abs, max_rel, worst_ip


U_CAP = 8.0  # proof window |u| <= U_CAP; equals the validated ApproxExp domain


def series_one_minus(x, gamma):
    """Cancellation-free reference for 1 - S(x,gamma) = sum_{k>=1} (-1)^{k+1}
    e^{-kx-k^2 gamma}.  Factoring out the leading term e^{-x-gamma} keeps full
    double relative accuracy even when 1-S ~ e^{-x} is tiny (deep window)."""
    lead = math.exp(-x - gamma)
    s = 0.0
    for k in range(1, 2000):
        t = (-1.0) ** (k + 1) * math.exp(-(k - 1) * x - (k * k - 1) * gamma)
        s += t
        if abs(t) < 1e-22:
            break
    return lead * s


def deep_deviation(r, nrm_low, gamma_q52, M_q62):
    """Worst two-sided relative kernel deviation over the DEEP proof window
    4 gamma < |u| <= U_CAP (every integer ip; both signs).

    threshold_deviation covers |u| <= 4 gamma with den = min(theta, 1-theta);
    out there the small side ~ e^{-|u|} shrinks, so the reference must be
    computed cancellation-free (series_one_minus) and the fixed-point small
    side as an exact integer difference, else float64 subtraction noise
    (2^-53) swamps the ~2^-50 true deviation.  The float reference itself
    carries ~2^-50 relative noise, so the reported deep deviations are
    accurate only to ~1 bit -- immaterial against the 2^-40 gate.  The
    reference is the M=1 kernel; the exact-M kernel differs by M-1 ~ 2^-107,
    likewise immaterial."""
    gamma = nrm_low / (2.0 * r * r)
    c = precompute_c(gamma_q52)
    ip_hi = int(U_CAP * r * r)
    max_rel = 0.0
    worst_ip = 0
    for aip in range(2 * nrm_low + 1, ip_hi + 1):
        for ip in (aip, -aip):
            thres = step_fixed(ip, r, gamma_q52, M_q62, 0, c)[2]
            if ip > 0:
                fixed_small = (M_q62 - thres) / M_q62   # plus side 1-f
                fixed_big = thres / M_q62               # minus side f
            else:
                fixed_small = thres / M_q62             # minus side f
                fixed_big = (M_q62 - thres) / M_q62     # plus side 1-f
            small_exact = series_one_minus(abs(ip) / (r * r), gamma)
            rel = max(abs(fixed_small - small_exact) / small_exact,
                      abs(fixed_big - (1.0 - small_exact)) / (1.0 - small_exact))
            if rel > max_rel:
                max_rel = rel
                worst_ip = ip
    return max_rel, worst_ip


def decide_fixed(ip, r, gamma_q52, M_q62, r63):
    return step_fixed(ip, r, gamma_q52, M_q62, r63)[0]


def grid_hash(r, M_q62):
    """Match test_irs.c grid: nrms x ip x r63 -> FNV-1a 64 over decision bits."""
    h = 1469598103934665603
    MASK = (1 << 64) - 1
    minus = total = 0
    nrms = [500, 2000, 4000]
    r63s = [1 << 10, 1 << 40, 1 << 62, (1 << 63) - 1, 0x123456789]
    for nrm in nrms:
        g = (nrm << 52) // (2 * r * r)
        ip = -3000
        while ip <= 3000:
            for r63 in r63s:
                b = decide_fixed(ip, r, g, M_q62, r63) & 1
                h = ((h ^ b) * 1099511628211) & MASK
                minus += b; total += 1
            ip += 50
    return h, minus, total


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--grid":
        for sec, (r, Bk, Bkp) in LEVELS.items():
            M = Mrep(r / Bk); M_q62 = round(M * ONE62)
            h, mn, tot = grid_hash(r, M_q62)
            print(f"GRID mode={sec} hash={h} minus={mn} total={tot}")
        return 0

    print(f"irs_ref.py - IRS exp-recurrence reference (T={T}) + M / consts emission\n")
    consts = {}
    allok = True
    for sec, (r, Bk, Bkp) in LEVELS.items():
        alpha = r / Bk
        M = Mrep(alpha); M_q62 = round(M * ONE62)
        consts[sec] = (alpha, M, M_q62)
        # KeyGen accepts BK_LOW_SQ <= ||s||^2 <= BK_SQ.  The series tail is
        # worst at the smallest gamma while the theta mass peaks at the
        # largest, so the kernel gate is verified at BOTH integer endpoints.
        gate = 2 ** -40
        alpha_max = r / math.sqrt(BK_LOW_SQ[sec])
        print(f"--- IRS-{sec}: r={r} B_k={Bk:.4f} B_k'={Bkp:.4f} alpha_max=r/sqrt(BK_LOW_SQ)={alpha_max:.4f} ---")
        print(f"  M_alpha = {M:.12f}  M_q62 = {M_q62}")
        for tag, nrm in (("BK_LOW_SQ", BK_LOW_SQ[sec]), ("BK_SQ", BK_SQ[sec])):
            gamma_q52 = (nrm << 52) // (2 * r * r)
            gamma = nrm / (2.0 * r * r)
            abs_err, mrel, worst_ip = threshold_deviation(r, nrm, gamma_q52, M_q62)
            deep_rel, deep_ip = deep_deviation(r, nrm, gamma_q52, M_q62)
            ok = mrel < gate and deep_rel < gate
            allok &= ok
            abslog = math.log2(abs_err) if abs_err > 0 else -99.0
            mrlog = math.log2(mrel) if mrel > 0 else -99.0
            deeplog = math.log2(deep_rel) if deep_rel > 0 else -99.0
            print(f"  endpoint {tag}={nrm}: gamma={gamma:.5f}  gamma_q52={gamma_q52}")
            print(f"    threshold abs_err={abs_err:.2e} (2^{abslog:.1f})  max_rel={mrel:.2e} "
                  f"(2^{mrlog:.1f})  worst ip={worst_ip}  [gate rel < 2^-40] -> {'OK' if ok else 'FAIL'}")
            print(f"    deep window 4gamma<|u|<={U_CAP:.0f}: max_rel={deep_rel:.2e} "
                  f"(2^{deeplog:.1f})  worst ip={deep_ip}")

    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.normpath(os.path.join(here, "..", "ref", "irs_rans.h"))
    RSH = 52  # RECIP_SHIFT
    mask64 = (1 << 64) - 1
    L = []
    L.append("/* Lithium rejection-sampling constants (see tools/irs_ref.py).")
    L.append(" * REJ_T: series truncation depth.  M_Q62: per-level rep rate (Q62).")
    L.append(" * LITHIUM_*_RECIP_{HI,LO}: 128-bit reciprocals REC=floor(2^(52+SHIFT)/d) for the")
    L.append(" *   division-free u=ip*2^52/r^2 (d=r^2) and gamma=nrm*2^52/(2r^2) (d=2r^2). */")
    L.append(f"#define REJ_T {T}")
    L.append(f"#define RECIP_SHIFT {RSH}")
    for sec, (alpha, M, M_q62) in consts.items():
        r = LEVELS[sec][0]
        r2, tr2 = r * r, 2 * r * r
        rec_r2 = (1 << (52 + RSH)) // r2
        rec_2r2 = (1 << (52 + RSH)) // tr2
        L.append(f"#if LITHIUM_MODE == {sec}")
        L.append(f"/* alpha = {alpha:.6f}, M = {M:.12f}; r = {r}, r^2 = {r2} */")
        L.append(f"#define M_Q62 INT64_C({M_q62})")
        L.append(f"#define R2_RECIP_HI  UINT64_C({rec_r2 >> 64})")
        L.append(f"#define R2_RECIP_LO  UINT64_C({rec_r2 & mask64})")
        L.append(f"#define TWO_R2_RECIP_HI UINT64_C({rec_2r2 >> 64})")
        L.append(f"#define TWO_R2_RECIP_LO UINT64_C({rec_2r2 & mask64})")
        L.append(f"#endif")
    autogen.patch_region(path, "irs_consts", "\n".join(L) + "\n")
    print(f"\nPatched @@AUTOGEN:irs_consts@@ region in {path}")
    print("ALL OK" if allok else "FAILURES")
    return 0 if allok else 1


if __name__ == "__main__":
    raise SystemExit(main())
