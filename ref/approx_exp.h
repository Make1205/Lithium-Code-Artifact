/*
 * approx_exp.h - fixed-point exp(-a) for IRS (constants/validation in
 * tools/gen_approx_exp.py; every constant is reproduced there).
 *
 *   approx_exp(a_q60)      : exp(-a) * 2^63, a = a_q60/2^60 >= 0, a in [0,
 * ~8) (sampler Bernoulli accept).  [7/8 ln2, ln2) degree-6 minimax (Q63)
 * + 8-entry 2^{-j/8} table + 2^{-n} shift.  Q60 holds a < 16 in a
 * uint64. approx_exp_4x(out,a)   : 4 independent approx_exp calls, batched
 * for ILP (4-way Horner + 4-way constant-time lookup); each lane is
 * bit-identical to approx_exp(a[l]). approx_exp_wide(a_q52) : exp(-a) *
 * 2^63, a in [0, 40] (irs.c rejection recurrence).  Q52-aware table range
 * reduction + degree-7 minimax. exp_small_signed(d_q52): exp(d) * 2^62 for
 * small SIGNED d (|d| < ln2/8).
 *
 * Constant time: fixed-degree Horner + full-table-scan lookup, no
 * data-dependent branches. The final variable shift count n is
 * secret-derived but range-bounded; x86 integer shifts are fixed-latency
 * for these operands.
 */
#ifndef LITHIUM_APPROX_EXP_H
#define LITHIUM_APPROX_EXP_H

#include <stdint.h>

#include "config.h" /* LITHIUM_MODE for the per-mode constants below */

#if defined(__GNUC__) || defined(__clang__)
__extension__ typedef unsigned __int128 wide_uint128;
__extension__ typedef __int128 wide_int128;
static inline uint64_t mulh64(uint64_t a, uint64_t b)
{
    return (uint64_t)(((wide_uint128)a * b) >> 64);
}
static inline int64_t smulh64(int64_t a, int64_t b)
{
    return (int64_t)(((wide_int128)a * b) >> 64);
}
#elif defined(_MSC_VER)
#    include <intrin.h>
static inline uint64_t mulh64(uint64_t a, uint64_t b)
{
    return __umulh(a, b);
}
static inline int64_t smulh64(int64_t a, int64_t b)
{
    return __mulh(a, b);
}
#else
#    error "need 128-bit multiply"
#endif

/* Fixed-point exp constants (minimax/Taylor coeffs, range-reduction
 * reciprocals, per-mode 2r^2 reciprocal).  Every value is derived +
 * validated in tools/gen_approx_exp.py; regenerate via `make tables`.  Do
 * not edit by hand. */
/* @@AUTOGEN:approx_exp@@ BEGIN */
/* ===== sampler scheme: approx_exp(a_q60) = round(exp(-a)*2^63) =====
 * exp(-a) = 2^{-n} * 2^{-j/B} * exp(r_l),  m=ceil(B*a/ln2)+(B-1),
 * n=m>>BETA, j=m&(B-1), r_l in [(B-1)/B ln2, ln2) = [7/8 ln2, ln2). */
#define APPROX_EXP_BETA 3
#define APPROX_EXP_B (1u << APPROX_EXP_BETA)

/* RCP = floor((B/ln2) * 2^57); V = mulh64(a_q60, RCP) =
 * floor((B*a/ln2)*2^53). */
#define APPROX_EXP_RCP_SHIFT 57
#define APPROX_EXP_RCP UINT64_C(1663314137230540311)

/* LB = round((ln2/B) * 2^(64+BETA)); R = m*LB - (a_q60<<(BETA+4)), X =
 * R>>BETA (Q64). */
#define APPROX_EXP_LB_SHIFT 67
#define APPROX_EXP_LB UINT64_C(12786308645202655660)

/* Degree-6 Q63 minimax coeffs EXP_POLY[k], P(r_l) ~ exp(r_l)
 * (all positive -> unsigned Horner).  Source (bit-exact):
 *   tools/LogPolyApprox/expx_7div8mulln2_ln2_43_64.txt
 * (GALACTICS Sobolev fit + LLL Babai-rounding to Q63). */
#define EXP_POLY_DEG 6
static const uint64_t EXP_POLY[EXP_POLY_DEG + 1] = {
    UINT64_C(9223528939368774768), UINT64_C(9221695796025543192),
    UINT64_C(4619330203765820209), UINT64_C(1517980426284622595),
    UINT64_C(413086442818934603),  UINT64_C(51568330914993427),
    UINT64_C(24532161584429622)};

/* Q63 table T[j] = round(2^{-j/B} * 2^63), j in 0..B-1 (constant-time
 * scan). */
static const uint64_t APPROX_EXP_TABLE_T[APPROX_EXP_B] = {
    UINT64_C(9223372036854775808), UINT64_C(8457869449776733335),
    UINT64_C(7755900482342532474), UINT64_C(7112192101001162095),
    UINT64_C(6521908912666391106), UINT64_C(5980616842327661685),
    UINT64_C(5484249825272419512), UINT64_C(5029079263719320435),
};

/* ===== IRS exp scheme: approx_exp_wide(a_q52), exp_small_signed(d_q52)
 * ===== Same B=8 table reduction as the sampler, but Q52 input: V =
 * mulh64(a_q52, APPROX_EXP_RCP) = floor((B*a/ln2)*2^45); X = R>>BETA is
 * Q64.  The signed small path maps |d| < ln2/B into the same residual
 * interval with m=B-1 for d>=0 and m=B for d<0. */
#define APPROX_EXP_Q52_CEIL_SHIFT 45
#define APPROX_EXP_Q52_ARG_SHIFT 15
#define APPROX_EXP_STEP_Q52 UINT64_C(390207173010335)
#define EXP_IRS_POLY_DEG 7
static const uint64_t EXP_IRS_POLY[EXP_IRS_POLY_DEG + 1] = {
    UINT64_C(9223359195834938053), UINT64_C(9223529185475211575),
    UINT64_C(4610846930460523822), UINT64_C(1539778780641310092),
    UINT64_C(379492631815806242),  UINT64_C(82618948915657009),
    UINT64_C(8594191333134702),    UINT64_C(3504640452209541)};

/* ===== per-mode Q80 reciprocal of 2 r^2: a_q60 = (num*I80)>>20 ===== */
#if LITHIUM_MODE == 120
/* 2 r^2 = 66248, I80 = floor(2^80/2r^2) = 18248487797588292094 */
#    define SAMP_TWO_R2 UINT64_C(66248)
#    define SAMP_I_HI UINT64_C(0)
#    define SAMP_I_LO UINT64_C(18248487797588292094)
#endif
#if LITHIUM_MODE == 180
/* 2 r^2 = 83232, I80 = floor(2^80/2r^2) = 14524771958076571207 */
#    define SAMP_TWO_R2 UINT64_C(83232)
#    define SAMP_I_HI UINT64_C(0)
#    define SAMP_I_LO UINT64_C(14524771958076571207)
#endif
#if LITHIUM_MODE == 260
/* 2 r^2 = 132098, I80 = floor(2^80/2r^2) = 9151734466945973252 */
#    define SAMP_TWO_R2 UINT64_C(132098)
#    define SAMP_I_HI UINT64_C(0)
#    define SAMP_I_LO UINT64_C(9151734466945973252)
#endif
/* @@AUTOGEN:approx_exp@@ END */

/* Sampler exp: exp(-a) * 2^63, a = a_q60 / 2^60 >= 0, a in [0, ~8). */
uint64_t approx_exp(uint64_t a_q60);

/* 4-way batched approx_exp: out[l] = approx_exp(a[l]) for l in 0..3,
 * bit-identical to four scalar calls but structured for instruction-level
 * parallelism. */
void approx_exp_4x(uint64_t *out, const uint64_t *a);

/* Batch exp for the sampler mini-batch: fills expv[0..n) = exp(-a[i]) *
 * 2^63. Single switch point for the exp backend so the scalar
 * (single-lane) and the AVX2 4-way samplers share ONE implementation:
 *   - default  : integer Q63 approx_exp_4x + scalar tail (bit-identical to
 * ref, KAT-locked).
 *   - AVX2_EXP : float AVX2 kernel (separate, opt-in KAT) from
 *                approx_exp_avx2.c via approx_exp_avx2(). */
void exp_batch(uint64_t *expv, const uint64_t *a_q60, int n);

#ifdef AVX2_EXP
/* Float AVX2 exp backend (defined in avx2/approx_exp_avx2.c). Fills
 * expv[0..n) = round(exp(-a[i]) * 2^63), a[i] = a_q60[i]/2^60 in [0,~8).
 */
void approx_exp_avx2(uint64_t *expv, const uint64_t *a_q60, int n);
#endif

#ifdef AVX512IFMA_EXP
/* AVX-512 IFMA integer exp backend (defined in
 * avx512/approx_exp_avx512ifma.c). Fills expv[0..n) = round(exp(-a[i]) *
 * 2^63) via a Q51 deg-5 vpmadd52 Horner (8 lanes/zmm), max relative error
 * < 2^-43.  Opt-in (AVX512IFMAEXP=1) with its OWN KAT (the Q51 path is not
 * bit-identical to the scalar Q63 kernel; full reconciliation is a TODO).
 */
void approx_exp_avx512ifma(uint64_t *expv, const uint64_t *a_q60, int n);
#endif

/* Wide exp: exp(-a) * 2^63, a = a_q52 / 2^52 >= 0, valid a in [0, ~40]
 * (irs.c rejection recurrence; table-reduced degree-7). */
uint64_t approx_exp_wide(uint64_t a_q52);

/* exp(d) * 2^62 for small SIGNED d = d_q52 / 2^52, |d| < ln2/8.  Used by
 * irs.c for the clamped "bad" series base exp(a_c). */
int64_t exp_small_signed(int64_t d_q52);

/* Two-lane IRS exp schedule: w_q62=exp(-a_q52/2^52),
 * wp_q62=exp(d_q52/2^52). The two Horner chains and table scans are
 * interleaved to expose ILP. */
void approx_exp_irs_pair(int64_t *w_q62, int64_t *wp_q62, uint64_t a_q52,
                         int64_t d_q52);

#endif /* LITHIUM_APPROX_EXP_H */
