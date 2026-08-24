/*
 * approx_exp.c - fixed-point exp(-a) for IRS.
 *
 * Two routines (constants in approx_exp.h, all reproduced by
 * tools/gen_approx_exp.py):
 *
 *   approx_exp(a_q60)       SAMPLER hot path.  round(exp(-a) * 2^63),
 *                           a = a_q60/2^60 >= 0, a in [0, ~8).  Scheme =
 * the lookup-table ApproxExp scheme (table size B = 2^BETA): a "exp(+)"
 * pipeline that turns exp(-a) into 2^{-n} * 2^{-j/B} * exp(r_l) with the
 * residual r_l kept in [(B-1)/B ln2, ln2) so a short Sobolev-optimal
 * Horner (Q63, all positive coeffs -> unsigned) reaches the 2^-43 target.
 * 2^{-j/B} comes from a B-entry constant-time table; 2^{-n} is a
 * round-to-nearest shift. B/degree are tuned by the (B,deg)x k-way sweep
 *                           (docs/design-notes/x64_ApproxExp.md):
 * B=8 (BETA=3), degree-6 -- the B-entry table scan dominates cost, so a
 * smaller table with a slightly higher degree is faster.  All constants
 * (incl. degree) come from approx_exp.h so this code is degree/B-agnostic.
 *
 *   approx_exp_wide(a_q52)  irs.c rejection recurrence.  round(exp(-a) *
 * 2^63), a in [0, 40].  It uses the same B=8 table range reduction as
 * approx_exp, but with Q52 input and a degree-7 IRS polynomial.
 *
 *   exp_small_signed(d_q52) Q62 table-reduced exp(d) for small SIGNED d
 * (|d| < ln2/8; irs.c uses d in [0,2gamma]).
 *
 * Constant time apart from the range-bounded variable shift by n.  n is
 * secret-derived (through a), so the security argument relies on
 * fixed-latency integer shifts plus the branchless n>=64 guard below, not
 * on publicness.
 */
#include "approx_exp.h"

/* ============================================================
 * Constant-time helpers (shared by approx_exp and approx_exp_4x)
 * ============================================================ */

/* Constant-time B-entry table lookup: scan the whole table, OR-masking the
 * matching slot.  j is treated as private (a deterministic function of the
 * secret seed via the CDT), so a data-dependent indexed load would leak
 * via the cache.  mask = -((i==idx)?1:0) without branching: diff==0 ->
 * (diff | -diff)>>63 == 0 -> mask = -1 diff!=0 -> (diff | -diff)>>63 == 1
 * -> mask =  0  */
static inline uint64_t ct_lookup_table(const uint64_t *table, uint64_t idx,
                                       uint64_t size)
{
    uint64_t r = 0;
    for (uint64_t i = 0; i < size; ++i) {
        uint64_t diff = i ^ idx;
        uint64_t mask = ((diff | (uint64_t)(-(int64_t)diff)) >> 63) - 1U;
        r |= table[i] & mask;
    }
    return r;
}

/* Step 5 tail: combine W (= 2^{-j/B}*exp(r_l) in Q63, value in [1,2)) with
 * the 2^{-n} factor by a round-to-nearest right shift.  n is
 * secret-derived but range-bounded; x86 shifts are fixed-latency.
 * Branchless n>=64 -> 0 guard (a underflows below the Q63 ULP); this
 * matches the original Taylor "(n<64)?...:0" semantics and keeps
 * approx_exp / approx_exp_4x bit-identical at the large-a tail. */
static inline uint64_t exp_combine_shift(uint64_t W, uint64_t n)
{
    uint64_t lt = (uint64_t)(-(int64_t)(n < 64)); /* all-ones iff n < 64 */
    unsigned sh = (unsigned)(n & 63);
    uint64_t round_bit = ((UINT64_C(1) << sh) >> 1) & lt;
    return (uint64_t)((((wide_uint128)W + round_bit) >> sh)) & lt;
}

/* ============================================================
 * Sampler exp: approx_exp(a_q60) = round(exp(-a) * 2^63)
 * ============================================================ */
uint64_t approx_exp(uint64_t a_q60)
{
    /* Step 1: m = ceil(B*a/ln2) + (B-1);  n = m>>BETA;  j = m & (B-1).
     * V = mulh64(a_q60, RCP) = floor((B*a/ln2)*2^53); +(2^53-1)>>53 makes
     * it ceil. */
    uint64_t V = mulh64(a_q60, APPROX_EXP_RCP);
    uint64_t m =
        ((V + ((UINT64_C(1) << 53) - 1U)) >> 53) + (APPROX_EXP_B - 1U);
    uint64_t n = m >> APPROX_EXP_BETA;
    uint64_t j = m & (APPROX_EXP_B - 1U);

    /* Step 2: r_l = (m/B)*ln2 - a in [(B-1)/B ln2, ln2).
     * rwide = m*LB - (a_q60 << (BETA+4)) = r_l * 2^(64+BETA);  X =
     * rwide>>BETA is Q64. */
    wide_uint128 rwide = (wide_uint128)m * APPROX_EXP_LB -
                         ((wide_uint128)a_q60 << (APPROX_EXP_BETA + 4));
    uint64_t X = (uint64_t)(rwide >> APPROX_EXP_BETA);

    /* Step 3: degree-D Horner P(r_l) ~ exp(r_l) in Q63.  Each step
     * p = mulh64(p,X)+C_k stays < 2^64 because exp(r_l) < exp(ln2) = 2. */
    uint64_t p = EXP_POLY[EXP_POLY_DEG];
    for (int k = EXP_POLY_DEG - 1; k >= 0; k--)
        p = mulh64(p, X) + EXP_POLY[k];

    /* Step 4: s = 2^{-j/B} in Q63 (constant-time). */
    uint64_t s = ct_lookup_table(APPROX_EXP_TABLE_T, j, APPROX_EXP_B);

    /* Step 5: W = p*s in Q63 (mulh gives Q62, <<1 -> Q63), then *2^{-n}.
     */
    uint64_t W = mulh64(p, s) << 1;
    return exp_combine_shift(W, n);
}

/* ============================================================
 * 4-way batched sampler exp: out[l] = approx_exp(a[l]) for l in 0..3.
 *
 * The four Gaussian samples are independent, so we batch four approx_exp
 * calls. Each lane's result is BIT-IDENTICAL to the scalar approx_exp
 * (same constants, same lookup masks, same shift) -- only the schedule
 * differs.  The win is microarchitectural:
 *   - Step 3 runs as a 4-wide Horner (step-major): the 4 mulh64 at each
 * Horner step are independent, so they fill the multiplier port at
 * ~1/cycle throughput instead of stalling on the scalar chain's ~3-cycle
 * mul latency.
 *   - Step 4 loads each table entry T[i] ONCE and applies 4 independent
 * masks, amortizing the B loads over 4 lanes and exposing 4 independent OR
 * chains that the out-of-order engine overlaps with the Horner multiplies.
 * `restrict` lets the compiler hoist the 4 loads and sink the 4 stores.
 * ============================================================ */
void approx_exp_4x(uint64_t *restrict out, const uint64_t *restrict a)
{
    uint64_t X[4], n[4], j[4], p[4], s[4];

    /* Steps 1+2: range reduction per lane (no cross-lane dependence). */
    for (int l = 0; l < 4; l++) {
        uint64_t V = mulh64(a[l], APPROX_EXP_RCP);
        uint64_t m =
            ((V + ((UINT64_C(1) << 53) - 1U)) >> 53) + (APPROX_EXP_B - 1U);
        n[l] = m >> APPROX_EXP_BETA;
        j[l] = m & (APPROX_EXP_B - 1U);
        wide_uint128 rwide = (wide_uint128)m * APPROX_EXP_LB -
                             ((wide_uint128)a[l] << (APPROX_EXP_BETA + 4));
        X[l] = (uint64_t)(rwide >> APPROX_EXP_BETA);
    }

    /* Step 3: degree-D Horner, step-major / lane-minor so the 4 multiplies
     * per step are independent (each lane is its own accumulator p[l]).
     * The outer loop walks the coefficients; the 4 inner multiplies of a
     * step have no mutual dependency, so they issue to the multiplier port
     * back-to-back at ~1/cycle instead of stalling on the ~3-cycle mul
     * latency. */
    for (int l = 0; l < 4; l++)
        p[l] = EXP_POLY[EXP_POLY_DEG];
    for (int k = EXP_POLY_DEG - 1; k >= 0; k--) {
        uint64_t c = EXP_POLY[k];
        p[0] = mulh64(p[0], X[0]) + c;
        p[1] = mulh64(p[1], X[1]) + c;
        p[2] = mulh64(p[2], X[2]) + c;
        p[3] = mulh64(p[3], X[3]) + c;
    }

    /* Step 4: 4-way constant-time table scan.  One load per entry, four
     * independent masks (same formula as ct_lookup_table -> bit-identical
     * s[l]). Full B-entry scan, no data-dependent branch/index -> constant
     * time. */
    s[0] = s[1] = s[2] = s[3] = 0;
    for (uint64_t i = 0; i < APPROX_EXP_B; i++) {
        uint64_t Ti = APPROX_EXP_TABLE_T[i];
        uint64_t d0 = i ^ j[0], d1 = i ^ j[1], d2 = i ^ j[2],
                 d3 = i ^ j[3];
        uint64_t m0 = ((d0 | (uint64_t)(-(int64_t)d0)) >> 63) - 1U;
        uint64_t m1 = ((d1 | (uint64_t)(-(int64_t)d1)) >> 63) - 1U;
        uint64_t m2 = ((d2 | (uint64_t)(-(int64_t)d2)) >> 63) - 1U;
        uint64_t m3 = ((d3 | (uint64_t)(-(int64_t)d3)) >> 63) - 1U;
        s[0] |= Ti & m0;
        s[1] |= Ti & m1;
        s[2] |= Ti & m2;
        s[3] |= Ti & m3;
    }

    /* Step 5: combine p*s (Q63) and apply 2^{-n} per lane. */
    for (int l = 0; l < 4; l++) {
        uint64_t W = mulh64(p[l], s[l]) << 1;
        out[l] = exp_combine_shift(W, n[l]);
    }
}

/* ============================================================
 * Batch exp dispatch: the single switch point shared by the scalar and the
 * AVX2 4-way samplers.  Default path is the integer Q63 backend
 * (approx_exp_4x in groups of 4 + scalar tail) and is BIT-IDENTICAL to a
 * per-attempt approx_exp -- so the avx2 build with AVX2EXP=0 stays
 * KAT-locked to ref.  With -DAVX2_EXP it routes to the float AVX2
 * kernel (own opt-in KAT).
 * ============================================================ */
void exp_batch(uint64_t *expv, const uint64_t *a, int n)
{
#if defined(AVX512IFMA_EXP)
    approx_exp_avx512ifma(expv, a, n);
#elif defined(AVX2_EXP)
    approx_exp_avx2(expv, a, n);
#else
    int jb = 0;
    for (; jb + 4 <= n; jb += 4)
        approx_exp_4x(&expv[jb], &a[jb]);
    for (; jb < n; jb++)
        expv[jb] = approx_exp(a[jb]);
#endif
}

static inline void exp_reduce_wide_q52(uint64_t a_q52, uint64_t *X,
                                       uint64_t *n, uint64_t *j)
{
    /* Q52 inputs cannot be promoted to the sampler's Q60 scale without
     * overflowing for a near 40.  Instead, reuse the same reciprocal in
     * its natural mulh scale: V has 45 fractional bits. */
    uint64_t V = mulh64(a_q52, APPROX_EXP_RCP);
    uint64_t m =
        ((V + ((UINT64_C(1) << APPROX_EXP_Q52_CEIL_SHIFT) - 1U)) >>
         APPROX_EXP_Q52_CEIL_SHIFT) +
        (APPROX_EXP_B - 1U);

    /* m should place the residual in [(B-1)/B ln2, ln2).  The reciprocal
     * can undershoot by one cell at boundaries, so bump once if the
     * residual would fall below the lower endpoint. */
    wide_uint128 rwide = (wide_uint128)m * APPROX_EXP_LB -
                         ((wide_uint128)a_q52 << APPROX_EXP_Q52_ARG_SHIFT);
    const wide_uint128 LOW =
        (wide_uint128)(APPROX_EXP_B - 1U) * APPROX_EXP_LB;
    uint64_t fix = (uint64_t)(rwide < LOW);
    m += fix;
    rwide += (wide_uint128)fix * APPROX_EXP_LB;

    *n = m >> APPROX_EXP_BETA;
    *j = m & (APPROX_EXP_B - 1U);
    *X = (uint64_t)(rwide >> APPROX_EXP_BETA);
}

static inline void exp_reduce_small_signed_q52(int64_t d_q52, uint64_t *X,
                                               uint64_t *n, uint64_t *j)
{
    /* exp_small_signed computes exp(+d).  Pick a public-form affine cell
     * so r = m/B*ln2 + d lands on the same positive residual interval used
     * by exp(-a).  The ge_step case covers the actual IRS max d slightly
     * beyond one table cell while staying below ln2/B. */
    uint64_t neg = (uint64_t)(d_q52 < 0);
    uint64_t nonneg = neg ^ 1U;
    uint64_t ge_step =
        nonneg & (uint64_t)((uint64_t)d_q52 >= APPROX_EXP_STEP_Q52);
    uint64_t m = (APPROX_EXP_B - 1U) + neg - ge_step;
    wide_int128 rwide =
        (wide_int128)((wide_uint128)m * APPROX_EXP_LB) +
        (wide_int128)d_q52 * ((wide_int128)1 << APPROX_EXP_Q52_ARG_SHIFT);

    *n = m >> APPROX_EXP_BETA;
    *j = m & (APPROX_EXP_B - 1U);
    *X = (uint64_t)((wide_uint128)rwide >> APPROX_EXP_BETA);
}

static inline uint64_t exp_irs_poly_q63(uint64_t X)
{
    uint64_t p = EXP_IRS_POLY[EXP_IRS_POLY_DEG];
    for (int k = EXP_IRS_POLY_DEG - 1; k >= 0; k--)
        p = mulh64(p, X) + EXP_IRS_POLY[k];
    return p;
}

static inline uint64_t exp_table_scale_q63(uint64_t p, uint64_t j)
{
    uint64_t s = ct_lookup_table(APPROX_EXP_TABLE_T, j, APPROX_EXP_B);
    return mulh64(p, s) << 1;
}

/* ============================================================
 * Wide exp: approx_exp_wide(a_q52) = round(exp(-a) * 2^63).  Used by
 * irs.c.  Q52 range reduction keeps the same residual interval as
 * approx_exp, so the polynomial degree drops from 18 to 7.
 * ============================================================ */
uint64_t approx_exp_wide(uint64_t a_q52)
{
    uint64_t X, n, j;
    exp_reduce_wide_q52(a_q52, &X, &n, &j);
    uint64_t W = exp_table_scale_q63(exp_irs_poly_q63(X), j);
    return exp_combine_shift(W, n);
}

int64_t exp_small_signed(int64_t d_q52)
{
    uint64_t X, n, j;
    exp_reduce_small_signed_q52(d_q52, &X, &n, &j);
    uint64_t W = exp_table_scale_q63(exp_irs_poly_q63(X), j);
    return (int64_t)(exp_combine_shift(W, n) >> 1);
}

void approx_exp_irs_pair(int64_t *w_q62, int64_t *wp_q62, uint64_t a_q52,
                         int64_t d_q52)
{
    uint64_t X0, n0, j0, X1, n1, j1;
    exp_reduce_wide_q52(a_q52, &X0, &n0, &j0);
    exp_reduce_small_signed_q52(d_q52, &X1, &n1, &j1);

    /* Step-major Horner exposes the two independent multiply chains from
     * decide(): w=exp(-|u|) and wp=exp(a_c). */
    uint64_t p0 = EXP_IRS_POLY[EXP_IRS_POLY_DEG];
    uint64_t p1 = p0;
    for (int k = EXP_IRS_POLY_DEG - 1; k >= 0; k--) {
        uint64_t c = EXP_IRS_POLY[k];
        p0 = mulh64(p0, X0) + c;
        p1 = mulh64(p1, X1) + c;
    }

    uint64_t s0 = 0;
    uint64_t s1 = 0;
    /* Full scan: j0/j1 are secret-derived, so do not index the table. */
    for (uint64_t i = 0; i < APPROX_EXP_B; i++) {
        uint64_t Ti = APPROX_EXP_TABLE_T[i];
        uint64_t d0 = i ^ j0;
        uint64_t d1 = i ^ j1;
        uint64_t m0 = ((d0 | (uint64_t)(-(int64_t)d0)) >> 63) - 1U;
        uint64_t m1 = ((d1 | (uint64_t)(-(int64_t)d1)) >> 63) - 1U;
        s0 |= Ti & m0;
        s1 |= Ti & m1;
    }

    uint64_t W0 = mulh64(p0, s0) << 1;
    uint64_t W1 = mulh64(p1, s1) << 1;
    *w_q62 = (int64_t)(exp_combine_shift(W0, n0) >> 1);
    *wp_q62 = (int64_t)(exp_combine_shift(W1, n1) >> 1);
}
