#ifndef COPA_FP_VECTORS_H
#define COPA_FP_VECTORS_H
#include <immintrin.h>
#include <stdint.h>
/* Full 64x64 products from four unsigned 32-bit partial products. No
 * cross-sum overflows: mid has at most 34 bits. Round exactly as scalar
 * Q63. */
static inline __m256i copa_fp_mul4(__m256i a, __m256i b)
{
    const __m256i mask = _mm256_set1_epi64x(UINT64_C(0xffffffff));
    __m256i ah = _mm256_srli_epi64(a, 32), bh = _mm256_srli_epi64(b, 32);
    __m256i p00 = _mm256_mul_epu32(a, b), p01 = _mm256_mul_epu32(a, bh),
            p10 = _mm256_mul_epu32(ah, b), p11 = _mm256_mul_epu32(ah, bh);
    __m256i mid =
        _mm256_add_epi64(_mm256_srli_epi64(p00, 32),
                         _mm256_add_epi64(_mm256_and_si256(p01, mask),
                                          _mm256_and_si256(p10, mask)));
    __m256i high = _mm256_add_epi64(
        p11,
        _mm256_add_epi64(_mm256_srli_epi64(p01, 32),
                         _mm256_add_epi64(_mm256_srli_epi64(p10, 32),
                                          _mm256_srli_epi64(mid, 32))));
    __m256i low = _mm256_or_si256(_mm256_and_si256(p00, mask),
                                  _mm256_slli_epi64(mid, 32));
    return _mm256_add_epi64(_mm256_or_si256(_mm256_slli_epi64(high, 1),
                                            _mm256_srli_epi64(low, 63)),
                            _mm256_and_si256(_mm256_srli_epi64(low, 62),
                                             _mm256_set1_epi64x(1)));
}
static inline void copa_backend_mulq63x4(uint64_t out[4],
                                         const uint64_t a[4],
                                         const uint64_t b[4])
{
    _mm256_storeu_si256(
        (__m256i *)out,
        copa_fp_mul4(_mm256_loadu_si256((const __m256i *)a),
                     _mm256_loadu_si256((const __m256i *)b)));
}
#if COPA_SIMD_BITS == 512
static inline __m512i copa_fp_mul8(__m512i x,__m512i y){
    const __m512i mask = _mm512_set1_epi64(UINT64_C(0xffffffff));
    __m512i xh = _mm512_srli_epi64(x, 32), yh = _mm512_srli_epi64(y, 32);
    __m512i p00 = _mm512_mul_epu32(x, y), p01 = _mm512_mul_epu32(x, yh),
            p10 = _mm512_mul_epu32(xh, y), p11 = _mm512_mul_epu32(xh, yh);
    __m512i mid =
        _mm512_add_epi64(_mm512_srli_epi64(p00, 32),
                         _mm512_add_epi64(_mm512_and_si512(p01, mask),
                                          _mm512_and_si512(p10, mask)));
    __m512i high = _mm512_add_epi64(
        p11,
        _mm512_add_epi64(_mm512_srli_epi64(p01, 32),
                         _mm512_add_epi64(_mm512_srli_epi64(p10, 32),
                                          _mm512_srli_epi64(mid, 32))));
    __m512i low = _mm512_or_si512(_mm512_and_si512(p00, mask),
                                  _mm512_slli_epi64(mid, 32));
    return _mm512_add_epi64(_mm512_or_si512(_mm512_slli_epi64(high, 1),
                                              _mm512_srli_epi64(low, 63)),
                              _mm512_and_si512(_mm512_srli_epi64(low, 62),
                                               _mm512_set1_epi64(1)));
}
#endif
static inline void copa_backend_mulq63x8(uint64_t out[8],
                                         const uint64_t a[8],
                                         const uint64_t b[8])
{
#if COPA_SIMD_BITS == 512
    _mm512_storeu_si512(out,copa_fp_mul8(_mm512_loadu_si512(a),_mm512_loadu_si512(b)));
#else
    copa_backend_mulq63x4(out, a, b);
    copa_backend_mulq63x4(out + 4, a + 4, b + 4);
#endif
}
#endif
