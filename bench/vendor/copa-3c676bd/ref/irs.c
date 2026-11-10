/*
 * irs.c - Iterative Rejection Sampling for IRS (exp-recurrence,
 * constant-time).
 *
 * See irs.h and tools/irs_ref.py (the bit-exact golden reference).  The
 * step is
 *
 *   a = |u|,  u = <z,v>/r^2,  gamma = ||s||^2/(2 r^2)
 *   Sp = S(a, gamma)      via base exp(-(a+gamma))          (approx_exp,
 * arg>=0) Sb = S(-a_c, gamma)   a_c = min(a, 2 gamma), base exp(a_c-gamma)
 *                                                            (exp_small_signed)
 *   S^+ = sgn? Sb : Sp ;  S^- = sgn? Sp : Sb     (sgn = u<0)
 *   thres_a  = (u>=2g)? S^+ : 1-S^-              (Q62)
 *   thres_ab = (|u|<2g)? 2-S^+-S^- : 1
 *   flag_minus = (r_rand*M < thres_a)            -> z-v  (else z+v)
 *
 * Reject (r_rand*M >= thres_ab) is ~2^-101 overall and ignored.
 */
#include "irs.h"

#include <stdint.h>
#include <string.h>

#include "approx_exp.h"
#include "params.h"
#include "polyvec.h"
#include "prof.h"
#include "reduce.h"
#include "symmetric.h"

/* ---- constant-time helpers (mask-based) ---- */
static inline int64_t ct_lt_s64(int64_t a, int64_t b)
{
    /* 1 if a < b, else 0; operands well within [-2^63/4, 2^63/4]. */
    return (int64_t)(((uint64_t)(a - b)) >> 63);
}
static inline int64_t ct_sel(int64_t cond, int64_t x, int64_t y)
{
    /* cond in {0,1}; m is all-zero or all-one. */
    int64_t m = -(cond & 1);
    return (x & m) | (y & ~m);
}
static inline int64_t mul62(int64_t x, int64_t y)
{
    return (int64_t)(((wide_int128)x * y) >> 62);
}
/* floor(x * 2^52 / d) for unsigned x, DIVISION-FREE / constant-time: a
 * 128-bit reciprocal  REC = floor(2^(52+RECIP_SHIFT)/d)
 * undershoots by at most 1 (x < 2^39 here vs the 2^128 fit bound), and one
 * branchless correction gives the exact floor.  Bit-identical to (x<<52)/d
 * -- but never emits a (variable- time) hardware/libgcc divide on
 * secret-derived x (and does not rely on the compiler strength-reducing
 * the 128-bit `/`). Mirrors tools/irs_ref.py. */
static inline uint64_t mul2p52_div(uint64_t x, uint64_t d,
                                   wide_uint128 REC)
{
    uint64_t q = (uint64_t)(((wide_uint128)x * REC) >> RECIP_SHIFT);
    uint64_t rem = (uint64_t)(((wide_uint128)x << 52) -
                              (wide_uint128)q * d); /* [0,2d) */
    q += (uint64_t)(rem >= d); /* branchless: 0 or 1 */
    return q;
}
static inline uint64_t load_le64(const uint8_t *p)
{
    return (uint64_t)p[0] | ((uint64_t)p[1] << 8) |
           ((uint64_t)p[2] << 16) | ((uint64_t)p[3] << 24) |
           ((uint64_t)p[4] << 32) | ((uint64_t)p[5] << 40) |
           ((uint64_t)p[6] << 48) | ((uint64_t)p[7] << 56);
}

/* v = x^j * in over Z[x]/(x^N+1): out[j..N-1] = in[0..N-1-j];
 * out[0..j-1] = -in[N-j..N-1] (negacyclic wrap).  j is the public
 * challenge index (no secret-dependent addressing). */
void poly_cyclic_shift(poly *out, const poly *in, unsigned int j)
{
    for (unsigned int k = 0; k < N; k++) {
        if (k >= j)
            out->coeffs[k] = in->coeffs[k - j];
        else
            out->coeffs[k] = -in->coeffs[k - j + N];
    }
}

int64_t compute_gamma_q52(const poly sk[K])
{
    int64_t nrm = 0;
    for (int c = 0; c < K; c++)
        nrm += poly_sqnorm(&sk[c]); /* ||s||^2 >= 0 */
    /* gamma_q52 = (nrm << 52) / (2 r^2), division-free (nrm >= 0). */
    const uint64_t d = 2u * (uint64_t)R * (uint64_t)R;
    const wide_uint128 REC =
        ((wide_uint128)TWO_R2_RECIP_HI << 64) | TWO_R2_RECIP_LO;
    return (int64_t)mul2p52_div((uint64_t)nrm, d, REC);
}

/* inner product <z, v> over the K-vector (exact int64). */
static int64_t inner(const poly z[K], const poly v[K])
{
    int64_t s = 0;
    for (int c = 0; c < K; c++)
        for (int i = 0; i < N; i++)
            s += (int64_t)z[c].coeffs[i] * v[c].coeffs[i];
    return s;
}

/* Per-signature coefficient table  ctab[k] = (-1)^k e^{-k^2 gamma}  (Q62)
 * -- the only gamma-dependent part of the series  S(.,gamma) = sum_k
 * ctab[k] w^k.  Built once per signature (gamma = ||s||^2/(2 r^2) is fixed
 * for all kappa steps) and reused by every decide.  Bases e^{-2g},
 * e^{-g} fold in here; the ratio e^{-(2k-1)g} = e^{-g} * (e^{-2g})^{k-1}
 * drives the k^2 recurrence.  Mirrors tools/irs_ref.py precompute_c
 * bit-for-bit. */
void precompute_c(int64_t *ctab, int64_t gamma_q52)
{
    const int64_t ONE62 = (int64_t)1 << 62;
    PROF_START(tx);
    int64_t qexp = (int64_t)(approx_exp_wide((uint64_t)(2 * gamma_q52)) >>
                             1); /* e^{-2g} */
    int64_t g1 =
        (int64_t)(approx_exp_wide((uint64_t)gamma_q52) >> 1); /* e^{-g}  */
    PROF_STOP(PT_EXP, tx);
    int64_t E = ONE62, r = g1;
    ctab[0] = ONE62; /* c_0 = 1 */
    for (int k = 1; k < REJ_T; k++) {
        E = mul62(E, r); /* E = e^{-k^2 g} */
        int64_t kneg = -(int64_t)(k & 1);
        ctab[k] = (E ^ kneg) - kneg; /* (-1)^k e^{-k^2 g} */
        r = mul62(r, qexp);          /* r = e^{-(2k+1)g} */
    }
}

/* Pure three-region decision from (ip = <z,v>, gamma_q52, precomputed
 * ctab, r63 uniform).  Returns flag_minus (1 => z-v).  Exposed for the
 * Python cross-check (tools/irs_ref.py).
 *
 * S(a,gamma) = sum_k ctab[k] w^k with w = e^{-a}: the good series uses
 * w = e^{-|u|} (wide range), the bad series w' = e^{a_c}, a_c =
 * min(|u|,2gamma) (small).  Each is one Horner over the per-signature ctab
 * (2 muls/term vs the old forward recurrence's 5), with one wide + one
 * small exp per call (vs 3). */
int decide(int64_t ip, int64_t gamma_q52, const int64_t *ctab,
           uint64_t r63)
{
    const int64_t ONE62 = (int64_t)1 << 62;

    /* u_q52 = trunc(ip * 2^52 / r^2), division-free / constant-time
     * (truncate toward zero: take |ip|, divide, reapply sign --
     * bit-identical to the signed 128-bit divide and to tools/irs_ref.py).
     */
    /* neg is all-ones iff ip < 0; use it to take |ip| branchlessly. */
    uint64_t neg = (uint64_t)0 - (uint64_t)(ip < 0);
    uint64_t aip = ((uint64_t)ip ^ neg) - neg;
    const uint64_t d = (uint64_t)R * (uint64_t)R;
    const wide_uint128 REC =
        ((wide_uint128)R2_RECIP_HI << 64) | R2_RECIP_LO;
    uint64_t uq = mul2p52_div(aip, d, REC);
    /* Reapply sign without branching. */
    int64_t u = (int64_t)((uq ^ neg) - neg);
    int64_t sgn = ct_lt_s64(u, 0);
    int64_t a = ct_sel(sgn, -u, u);
    int64_t two_g = 2 * gamma_q52;
    int64_t a_c = ct_sel(ct_lt_s64(a, two_g), a, two_g);

    PROF_START(tx);
    int64_t w, wp;
    /* The two exp inputs are independent once a and a_c are known; the
     * paired helper interleaves both Horner chains and uses one masked
     * table scan. */
    approx_exp_irs_pair(&w, &wp, (uint64_t)a, a_c);
    PROF_STOP(PT_EXP, tx);

    /* Sp = S(a,gamma) = Horner(ctab, w);  Sb = S(-a_c,gamma) =
     * Horner(ctab, w'). Intermediates stay < 2^62.x (alternating-series
     * cancellation) -> no int64 overflow (verified in tools/irs_ref.py).
     */
    int64_t Sp = ctab[REJ_T - 1];
    int64_t Sb = ctab[REJ_T - 1];
    for (int k = REJ_T - 2; k >= 0; k--) {
        Sp = ctab[k] + mul62(Sp, w);
        Sb = ctab[k] + mul62(Sb, wp);
    }

    int64_t S_pos = ct_sel(sgn, Sb, Sp);
    int64_t S_neg = ct_sel(sgn, Sp, Sb);

    int64_t is_ge = 1 - ct_lt_s64(u, two_g);
    int64_t thres_a = ct_sel(is_ge, S_pos, ONE62 - S_neg);
    int64_t rM = (int64_t)(((wide_uint128)r63 * (uint64_t)M_Q62) >> 63);
    return (int)(ct_lt_s64(rM, thres_a) & 1);
}

int step(poly z[K], const poly v[K], int64_t gamma_q52,
         const int64_t *ctab, keccak_state *rng)
{
    int64_t ip = inner(z, v);
    uint8_t rb[8];
    shake256_squeeze(rb, 8, rng);
    uint64_t r63 = load_le64(rb) >> 1;
    int flag_minus = decide(ip, gamma_q52, ctab, r63);

    /* m is -1 for the minus branch and 0 for the plus branch. */
    int32_t m = -(int32_t)(flag_minus & 1);
    for (int c = 0; c < K; c++)
        for (int i = 0; i < N; i++) {
            /* z += (m ? -v : +v): conditional-negate v then add.  m in
             * {0,-1}, so (v^m)-m == v (m=0) / -v (m=-1) -- bit-identical
             * to the old (minus&m)|(plus&~m) blend, fewer ops. */
            z[c].coeffs[i] += (v[c].coeffs[i] ^ m) - m;
        }
    return flag_minus;
}

void sign_irs(poly z[K], const poly *c, const poly sk[K],
              uint8_t *sign_bits, keccak_state *rng)
{
    int64_t gamma_q52 = compute_gamma_q52(sk);
    int64_t ctab[REJ_T];
    /* ctab depends only on this signature's gamma. */
    precompute_c(ctab, gamma_q52);
    poly v[K];
    unsigned cnt = 0;
    for (unsigned j = 0; j < N; j++) {
        if (c->coeffs[j] == 0)
            continue;
        for (int comp = 0; comp < K; comp++)
            poly_cyclic_shift(&v[comp], &sk[comp], j);
        int sb = step(z, v, gamma_q52, ctab, rng);
        if (sign_bits) {
            if (sb)
                sign_bits[cnt >> 3] |= (uint8_t)(1u << (cnt & 7));
            else
                sign_bits[cnt >> 3] &= (uint8_t) ~(1u << (cnt & 7));
        }
        cnt++;
    }
}
