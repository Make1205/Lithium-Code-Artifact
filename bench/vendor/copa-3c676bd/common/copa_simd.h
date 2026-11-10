#ifndef COPA_SIMD_H
#define COPA_SIMD_H
#include <immintrin.h>
#include <stdint.h>
#if COPA_SIMD_BITS == 256
typedef __m256i copav_vec;
#    define COPAV_W 8
#    define cv_set(x) _mm256_set1_epi32((int32_t)(x))
#    define cv_set64(x) _mm256_set1_epi64x((int64_t)(x))
#    define cv_load(p) _mm256_loadu_si256((const __m256i *)(p))
#    define cv_store(p, x) _mm256_storeu_si256((__m256i *)(p), (x))
#    define cv_add(a, b) _mm256_add_epi32(a, b)
#    define cv_sub(a, b) _mm256_sub_epi32(a, b)
#    define cv_add64(a, b) _mm256_add_epi64(a, b)
#    define cv_sub64(a, b) _mm256_sub_epi64(a, b)
#    define cv_and(a, b) _mm256_and_si256(a, b)
#    define cv_or(a, b) _mm256_or_si256(a, b)
#    define cv_xor(a, b) _mm256_xor_si256(a, b)
#    define cv_gt(a, b) _mm256_cmpgt_epi32(a, b)
#    define cv_eq(a, b) _mm256_cmpeq_epi32(a, b)
#    define cv_gt64(a, b) _mm256_cmpgt_epi64(a, b)
#    define cv_mul32(a, b) _mm256_mullo_epi32(a, b)
#    define cv_mul_even(a, b) _mm256_mul_epu32(a, b)
#    define cv_mul_signed(a, b) _mm256_mul_epi32(a, b)
#    define cv_srl64(a, n) _mm256_srli_epi64(a, n)
#    define cv_sll64(a, n) _mm256_slli_epi64(a, n)
#    define cv_srl32(a, n) _mm256_srli_epi32(a, n)
#    define cv_sll32(a, n) _mm256_slli_epi32(a, n)
#    define cv_sra32(a, n) _mm256_srai_epi32(a, n)
#    define cv_odd(a, b) _mm256_blend_epi32(a, b, 0xaa)
#    define cv_perm(a, i) _mm256_permutevar8x32_epi32(a, i)
#else
typedef __m512i copav_vec;
#    define COPAV_W 16
#    define cv_set(x) _mm512_set1_epi32((int32_t)(x))
#    define cv_set64(x) _mm512_set1_epi64((int64_t)(x))
#    define cv_load(p) _mm512_loadu_si512((const void *)(p))
#    define cv_store(p, x) _mm512_storeu_si512((void *)(p), (x))
#    define cv_add(a, b) _mm512_add_epi32(a, b)
#    define cv_sub(a, b) _mm512_sub_epi32(a, b)
#    define cv_add64(a, b) _mm512_add_epi64(a, b)
#    define cv_sub64(a, b) _mm512_sub_epi64(a, b)
#    define cv_and(a, b) _mm512_and_si512(a, b)
#    define cv_or(a, b) _mm512_or_si512(a, b)
#    define cv_xor(a, b) _mm512_xor_si512(a, b)
#    define cv_gt(a, b) \
        _mm512_maskz_set1_epi32(_mm512_cmpgt_epi32_mask(a, b), -1)
#    define cv_eq(a, b) \
        _mm512_maskz_set1_epi32(_mm512_cmpeq_epi32_mask(a, b), -1)
#    define cv_gt64(a, b) \
        _mm512_maskz_set1_epi64(_mm512_cmpgt_epi64_mask(a, b), -1)
#    define cv_mul32(a, b) _mm512_mullo_epi32(a, b)
#    define cv_mul_even(a, b) _mm512_mul_epu32(a, b)
#    define cv_mul_signed(a, b) _mm512_mul_epi32(a, b)
#    define cv_srl64(a, n) _mm512_srli_epi64(a, n)
#    define cv_sll64(a, n) _mm512_slli_epi64(a, n)
#    define cv_srl32(a, n) _mm512_srli_epi32(a, n)
#    define cv_sll32(a, n) _mm512_slli_epi32(a, n)
#    define cv_sra32(a, n) _mm512_srai_epi32(a, n)
#    define cv_odd(a, b) _mm512_mask_blend_epi32(0xaaaa, a, b)
#    define cv_perm(a, i) _mm512_permutexvar_epi32(i, a)
#endif
static inline copav_vec cv_select(copav_vec a, copav_vec b, copav_vec mask)
{
    return cv_or(cv_and(a, mask), cv_and(b, cv_xor(mask, cv_set(-1))));
}
static inline copav_vec cv_reduce(copav_vec a, uint32_t p)
{
#if COPA_SIMD_BITS == 256
    return _mm256_min_epu32(a,cv_sub(a,cv_set(p)));
#else
    return _mm512_min_epu32(a,cv_sub(a,cv_set(p)));
#endif
}
static inline copav_vec cv_high32(copav_vec a, copav_vec b)
{
    copav_vec even = cv_mul_even(a, b),
              odd = cv_mul_even(cv_srl64(a, 32), cv_srl64(b, 32));
    return cv_odd(cv_srl64(even, 32), odd);
}
static inline copav_vec cv_shoup(copav_vec a, copav_vec w, copav_vec ws,
                                 uint32_t p)
{
    return cv_reduce(
        cv_sub(cv_mul32(a, w), cv_mul32(cv_high32(a, ws), cv_set(p))), p);
}
static inline copav_vec cv_mont(copav_vec a, copav_vec b, uint32_t p,
                                uint32_t qi)
{
    copav_vec ze = cv_mul_even(a, b),
              zo = cv_mul_even(cv_srl64(a, 32), cv_srl64(b, 32));
    copav_vec low = cv_odd(ze, cv_sll64(zo, 32)),
              t = cv_mul32(low, cv_set(qi));
    copav_vec even = cv_add64(ze, cv_mul_even(t, cv_set(p))),
              odd = cv_add64(zo, cv_mul_even(cv_srl64(t, 32), cv_set(p)));
    return cv_reduce(cv_odd(cv_srl64(even, 32), odd), p);
}
#endif
