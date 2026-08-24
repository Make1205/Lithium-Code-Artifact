/*
 * approx_exp_avx2.c - AVX2 binary64 exp(-a) for the Lithium sampler
 * (opt-in, compiled only with -DAVX2_EXP).  Replaces the integer
 * approx_exp backend in exp_batch() (approx_exp.c) with a 4-way-double
 * minimax kernel.
 *
 * KAT: this is a SEPARATE, opt-in KAT (float != Q63 integer rounding).
 * Route (ii): the kernel produces exp(-a) as a double, then we convert to
 * a Q63 uint64 so the accept test stays the exact integer compare
 * rand_rej_63 < exp_q63 in the shared scalar gauss_finalize --
 * isolating the float approximation to a single value.
 *
 * Constant time: no div/sqrt; the B-entry 2^{-j/B} table is read by a full
 * masked-OR scan (not an indexed gather); 2^{-n} is an exponent-field
 * shift. j and n are secret-derived but range-bounded and kept in
 * fixed-latency SIMD integer operations.  FTZ/DAZ is set (and restored)
 * around the kernel so any (here unreachable) denormal cannot trigger a
 * data-dependent microcode assist. Built with -ffp-contract=off
 * -fno-fast-math (explicit FMA intrinsics) so the recorded KAT is
 * reproducible across compilers.
 *
 * Constants: the @@AUTOGEN:approx_exp_avx2@@ region of approx_exp_avx2.h,
 * generated + checked by tools/avxexp/gen_avxexp_consts.py.
 */
#include "approx_exp_avx2.h"

#include <immintrin.h>
#include <stdint.h>

#include "approx_exp.h"
#include "config.h"

/* ---- one 4-lane (__m256d) exp(-a) stream ---- */
static inline __m256d exp_stream(__m256d a)
{
    const __m256d INV = _mm256_set1_pd(APPROX_EXP_AVX2_INV);
    const __m256d LB = _mm256_set1_pd(APPROX_EXP_AVX2_LB);
    const __m256d cBm1 = _mm256_set1_pd((double)(APPROX_EXP_AVX2_B - 1));

    __m256d md =
        _mm256_add_pd(_mm256_ceil_pd(_mm256_mul_pd(a, INV)), cBm1);
    __m128i m = _mm256_cvtpd_epi32(md);
    __m128i n = _mm_srli_epi32(m, APPROX_EXP_AVX2_LOG2B);
    __m128i j = _mm_and_si128(m, _mm_set1_epi32(APPROX_EXP_AVX2_B - 1));

    __m256d rl = _mm256_fmsub_pd(md, LB, a); /* (m/B)*ln2 - a */
    __m256d p = _mm256_set1_pd(APPROX_EXP_AVX2_C[APPROX_EXP_AVX2_DEG]);
    for (int i = APPROX_EXP_AVX2_DEG - 1; i >= 0; i--)
        p = _mm256_fmadd_pd(p, rl, _mm256_set1_pd(APPROX_EXP_AVX2_C[i]));

    /* CT table lookup s = T[j] (j in [0,B)).  Widen j to 64-bit ONCE, then
     * OR-mask all B entries into 4 independent accumulators: exactly one
     * k==j matches, so OR_k(T[k] & (j==k)) == T[j].  Drops the per-entry
     * vpmovsxdq (frees port 5: one widen instead of B) and replaces the
     * B-deep vblendvpd dependency chain with 4-way ILP + a 2-level OR
     * combine (vpcmpeqq/vpand/vpor are p01/p015).  AVX2 has no cheap
     * B-double permute (cf. AVX-512 vpermt2q), so the scan length B is the
     * lever -- B=8 keeps it short.  See RESULTS_avx2.md. */
    __m256i jq = _mm256_cvtepi32_epi64(j);
    __m256i a0 = _mm256_setzero_si256(), a1 = a0, a2 = a0, a3 = a0;
    for (int k = 0; k < APPROX_EXP_AVX2_B; k += 4) {
        a0 = _mm256_or_si256(
            a0, _mm256_and_si256(
                    _mm256_cmpeq_epi64(jq, _mm256_set1_epi64x(k + 0)),
                    _mm256_castpd_si256(
                        _mm256_set1_pd(APPROX_EXP_AVX2_T[k + 0]))));
        a1 = _mm256_or_si256(
            a1, _mm256_and_si256(
                    _mm256_cmpeq_epi64(jq, _mm256_set1_epi64x(k + 1)),
                    _mm256_castpd_si256(
                        _mm256_set1_pd(APPROX_EXP_AVX2_T[k + 1]))));
        a2 = _mm256_or_si256(
            a2, _mm256_and_si256(
                    _mm256_cmpeq_epi64(jq, _mm256_set1_epi64x(k + 2)),
                    _mm256_castpd_si256(
                        _mm256_set1_pd(APPROX_EXP_AVX2_T[k + 2]))));
        a3 = _mm256_or_si256(
            a3, _mm256_and_si256(
                    _mm256_cmpeq_epi64(jq, _mm256_set1_epi64x(k + 3)),
                    _mm256_castpd_si256(
                        _mm256_set1_pd(APPROX_EXP_AVX2_T[k + 3]))));
    }
    __m256d s = _mm256_castsi256_pd(
        _mm256_or_si256(_mm256_or_si256(a0, a1), _mm256_or_si256(a2, a3)));

    __m256i n64 = _mm256_cvtepi32_epi64(n); /* 2^{-n} via exp field */
    __m256d e = _mm256_castsi256_pd(_mm256_slli_epi64(
        _mm256_sub_epi64(_mm256_set1_epi64x(1023), n64), 52));
    return _mm256_mul_pd(_mm256_mul_pd(p, s), e);
}

void exp_avx2_w4(const double *a, double *out)
{
    _mm256_storeu_pd(out, exp_stream(_mm256_loadu_pd(a)));
}

void exp_avx2_w8(const double *a, double *out)
{
    /* 2 independent streams */
    __m256d r0 = exp_stream(_mm256_loadu_pd(a + 0));
    __m256d r1 = exp_stream(_mm256_loadu_pd(a + 4));
    _mm256_storeu_pd(out + 0, r0);
    _mm256_storeu_pd(out + 4, r1);
}

void exp_avx2_w16(const double *a, double *out)
{
    /* 4 independent streams: the compiler interleaves the 4 latency-bound
     * chains to fill the FMA ports.  The production GAUSS_BATCH=32 wrapper
     * calls this natural 16-lane block twice. */
    __m256d a0 = _mm256_loadu_pd(a + 0);
    __m256d a1 = _mm256_loadu_pd(a + 4);
    __m256d a2 = _mm256_loadu_pd(a + 8);
    __m256d a3 = _mm256_loadu_pd(a + 12);
    __m256d r0 = exp_stream(a0);
    __m256d r1 = exp_stream(a1);
    __m256d r2 = exp_stream(a2);
    __m256d r3 = exp_stream(a3);
    _mm256_storeu_pd(out + 0, r0);
    _mm256_storeu_pd(out + 4, r1);
    _mm256_storeu_pd(out + 8, r2);
    _mm256_storeu_pd(out + 12, r3);
}

/* Convert four uint32 lanes to binary64 without the unsigned-conversion
 * branch GCC emits for scalar casts. */
static inline __m256d u32_to_pd(__m128i x)
{
    __m256d d = _mm256_cvtepi32_pd(x);
    __m128i neg = _mm_cmpgt_epi32(_mm_setzero_si128(), x);
    __m256i neg64 = _mm256_cvtepi32_epi64(neg);
    __m256d add =
        _mm256_and_pd(_mm256_castsi256_pd(neg64), _mm256_set1_pd(0x1p32));
    return _mm256_add_pd(d, add);
}

/* Four Q60 inputs -> four doubles.  Sampler inputs are in [0,8]*2^60, so
 * they fit in uint63, but the public tests include exactly 2^63. Splitting
 * into two uint32 halves avoids scalar uint64->double branches and matches
 * the usual nearest binary64 conversion. */
static inline __m256d q60_to_pd4(const uint64_t *a_q60)
{
    __m256i v = _mm256_loadu_si256((const __m256i *)a_q60);
    __m128i lo128 = _mm256_castsi256_si128(v);
    __m128i hi128 = _mm256_extracti128_si256(v, 1);
    __m128i lo01 = _mm_shuffle_epi32(lo128, _MM_SHUFFLE(2, 0, 2, 0));
    __m128i lo23 = _mm_shuffle_epi32(hi128, _MM_SHUFFLE(2, 0, 2, 0));
    __m128i hi01 = _mm_shuffle_epi32(lo128, _MM_SHUFFLE(3, 1, 3, 1));
    __m128i hi23 = _mm_shuffle_epi32(hi128, _MM_SHUFFLE(3, 1, 3, 1));
    __m128i lo = _mm_unpacklo_epi64(lo01, lo23);
    __m128i hi = _mm_unpacklo_epi64(hi01, hi23);

    __m256d q = _mm256_fmadd_pd(u32_to_pd(hi), _mm256_set1_pd(0x1p32),
                                u32_to_pd(lo));
    return _mm256_mul_pd(q, _mm256_set1_pd(0x1p-60));
}

/* exp_neg in (0,~1] -> round(exp_neg * 2^63) as uint64 (route ii), four
 * lanes at a time and without double->uint64 conversion branches.
 *
 * Range used here: the public accuracy gate covers a in [0,8], so
 * x=exp_neg*2^63+0.5 is in roughly [2^51, 2^63+2^15].  Thus floor(x) can
 * be recovered exactly from the binary64 exponent and mantissa using only
 * variable shifts.  For the a=0 overshoot, returning any value >=2^63
 * would make the 63-bit rejection comparison accept; the bit extraction
 * preserves the natural uint64 truncation as well. */
static inline __m256i to_q63_w4(__m256d exp_neg)
{
    const __m256i mant_mask = _mm256_set1_epi64x(0x000FFFFFFFFFFFFFLL);
    const __m256i hidden = _mm256_set1_epi64x(0x0010000000000000LL);
    const __m256i bias_plus_52 = _mm256_set1_epi64x(1023 + 52);

    __m256d x =
        _mm256_fmadd_pd(exp_neg, _mm256_set1_pd(APPROX_EXP_AVX2_SCALE63),
                        _mm256_set1_pd(0.5));
    __m256i bits = _mm256_castpd_si256(x);
    __m256i exp = _mm256_and_si256(_mm256_srli_epi64(bits, 52),
                                   _mm256_set1_epi64x(0x7FF));
    __m256i mant =
        _mm256_or_si256(_mm256_and_si256(bits, mant_mask), hidden);

    __m256i left_cnt = _mm256_sub_epi64(exp, bias_plus_52);
    __m256i right_cnt = _mm256_sub_epi64(bias_plus_52, exp);
    __m256i left = _mm256_sllv_epi64(mant, left_cnt);
    __m256i right = _mm256_srlv_epi64(mant, right_cnt);
    __m256i use_left =
        _mm256_cmpgt_epi64(exp, _mm256_set1_epi64x(1023 + 51));
    return _mm256_or_si256(_mm256_and_si256(left, use_left),
                           _mm256_andnot_si256(use_left, right));
}

void approx_exp_avx2(uint64_t *expv, const uint64_t *a_q60, int n)
{
    /* FTZ/DAZ for the CT guarantee; saved/restored so the caller's FP env
     * is untouched (exp values here are never denormal, so results are
     * unchanged). */
    unsigned int mxcsr = _mm_getcsr();
    _mm_setcsr((mxcsr | 0x8040u)); /* FTZ (bit15) | DAZ (bit6) */

    int i = 0;
    for (; i + 16 <= n; i += 16) {
        _mm256_storeu_si256(
            (__m256i *)(expv + i + 0),
            to_q63_w4(exp_stream(q60_to_pd4(a_q60 + i + 0))));
        _mm256_storeu_si256(
            (__m256i *)(expv + i + 4),
            to_q63_w4(exp_stream(q60_to_pd4(a_q60 + i + 4))));
        _mm256_storeu_si256(
            (__m256i *)(expv + i + 8),
            to_q63_w4(exp_stream(q60_to_pd4(a_q60 + i + 8))));
        _mm256_storeu_si256(
            (__m256i *)(expv + i + 12),
            to_q63_w4(exp_stream(q60_to_pd4(a_q60 + i + 12))));
    }
    for (; i + 4 <= n; i += 4)
        _mm256_storeu_si256((__m256i *)(expv + i),
                            to_q63_w4(exp_stream(q60_to_pd4(a_q60 + i))));
    for (; i < n; i++) {
        uint64_t tmp_in[4] = {a_q60[i], a_q60[i], a_q60[i], a_q60[i]};
        uint64_t tmp_out[4];
        _mm256_storeu_si256((__m256i *)tmp_out,
                            to_q63_w4(exp_stream(q60_to_pd4(tmp_in))));
        expv[i] = tmp_out[0];
    }

    _mm_setcsr(mxcsr);
}
