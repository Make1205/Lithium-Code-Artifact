#ifndef COPA_POLY_VECTORS_H
#define COPA_POLY_VECTORS_H
#include <stddef.h>

#include "copa_params.h"
#include "copa_simd.h"

static inline copav_vec copa_poly_wrap(copav_vec x, uint32_t q)
{
    return cv_reduce(cv_add(x, cv_and(cv_set(q), cv_gt(cv_set(0), x))), q);
}
static inline copav_vec copa_poly_center(copav_vec x, uint32_t q)
{
    return cv_sub(x, cv_and(cv_set(q), cv_gt(x, cv_set(q / 2))));
}
static inline copav_vec copa_poly_high(copav_vec x)
{
    return cv_reduce(
        cv_srl32(cv_mul32(cv_add(x, cv_set(COPA_ALPHA >> 1)), cv_set(COPA_ALPHA_RECIP)), 24),
        COPA_Q / COPA_ALPHA);
}
static inline copav_vec copa_poly_square(copav_vec x)
{
    copav_vec odd = cv_srl64(x, 32);
    return cv_add64(cv_mul_signed(x, x), cv_mul_signed(odd, odd));
}
static inline uint64_t copa_poly_sum64(copav_vec x)
{
    uint64_t v[COPAV_W / 2], total = 0;
    cv_store(v, x);
    for (unsigned j = 0; j < COPAV_W / 2; ++j)
        total += v[j];
    return total;
}
static inline uint64_t copa_backend_norm2(const int32_t *x, unsigned n)
{
    copav_vec sum = cv_set(0);
    const unsigned end = n - n % COPAV_W;
    unsigned j = 0;
    for (; j < end; j += COPAV_W)
        sum = cv_add64(sum, copa_poly_square(cv_load(x + j)));
    uint64_t total = copa_poly_sum64(sum);
    for (; j < n; ++j)
        total += (uint64_t)((int64_t)x[j] * x[j]);
    return total;
}
static inline void copa_backend_rotate_add(int32_t *out, const int32_t *s,
                                           unsigned shift, int32_t mask)
{
    __asm__("" : "+r"(mask));
    copav_vec mv = cv_set(mask);
    unsigned j = shift;
    for (; j + COPAV_W <= COPA_N; j += COPAV_W) {
        copav_vec v = cv_sub(cv_xor(cv_load(s + j - shift), mv), mv);
        cv_store(out + j, cv_add(cv_load(out + j), v));
    }
    for (; j < COPA_N; ++j)
        out[j] += (s[j - shift] ^ mask) - mask;
    for (j = 0; j + COPAV_W <= shift; j += COPAV_W) {
        copav_vec v =
            cv_sub(cv_xor(cv_load(s + j + COPA_N - shift), mv), mv);
        cv_store(out + j, cv_sub(cv_load(out + j), v));
    }
    for (; j < shift; ++j)
        out[j] -= (s[j + COPA_N - shift] ^ mask) - mask;
}
/* Public challenge positions address a duplicated negacyclic secret. Keep
 * each destination vector in registers across all tau signed rotations. */
static inline void copa_backend_sparse_add(int32_t *out,const int32_t *secret,
        const unsigned *position,const int32_t *masks,unsigned count)
{
    int32_t extended[2*COPA_N];
    for(unsigned j=0;j<COPA_N;j+=COPAV_W){
        copav_vec v=cv_load(secret+j);
        cv_store(extended+j,cv_sub(cv_set(0),v));
        cv_store(extended+COPA_N+j,v);
    }
    /* Eight accumulators fit AVX2's register file. Amortize the public
     * index and sign broadcast across eight vectors without spilling all N. */
    for(unsigned j=0;j<COPA_N;j+=8*COPAV_W){
        copav_vec sum[8];
        for(unsigned v=0;v<8;++v)sum[v]=cv_load(out+j+v*COPAV_W);
        for(unsigned k=0;k<count;++k){
            int32_t mask=masks[k];__asm__("":"+r"(mask));
            copav_vec m=cv_set(mask);const int32_t *at=extended+COPA_N-position[k]+j;
            sum[0]=cv_add(sum[0],cv_sub(cv_xor(cv_load(at+0*COPAV_W),m),m));
            sum[1]=cv_add(sum[1],cv_sub(cv_xor(cv_load(at+1*COPAV_W),m),m));
            sum[2]=cv_add(sum[2],cv_sub(cv_xor(cv_load(at+2*COPAV_W),m),m));
            sum[3]=cv_add(sum[3],cv_sub(cv_xor(cv_load(at+3*COPAV_W),m),m));
            sum[4]=cv_add(sum[4],cv_sub(cv_xor(cv_load(at+4*COPAV_W),m),m));
            sum[5]=cv_add(sum[5],cv_sub(cv_xor(cv_load(at+5*COPAV_W),m),m));
            sum[6]=cv_add(sum[6],cv_sub(cv_xor(cv_load(at+6*COPAV_W),m),m));
            sum[7]=cv_add(sum[7],cv_sub(cv_xor(cv_load(at+7*COPAV_W),m),m));
        }
        for(unsigned v=0;v<8;++v)cv_store(out+j+v*COPAV_W,sum[v]);
    }
}
static inline void copa_backend_commit(int32_t *w, int32_t *w1,
                                       const int32_t *y2)
{
    for (unsigned j = 0; j < COPA_M * COPA_N; j += COPAV_W) {
        copav_vec v = copa_poly_wrap(
            cv_add(cv_load(w + j), cv_load(y2 + j)), COPA_Q);
        cv_store(w + j, v);
        cv_store(w1 + j, copa_poly_high(v));
    }
}
static inline uint64_t copa_backend_make_hint(int32_t *h, const int32_t *w,
                                              const int32_t *w1,
                                              const int32_t *z2)
{
    copav_vec sum = cv_set(0);
    for (unsigned j = 0; j < COPA_M * COPA_N; j += COPAV_W) {
        copav_vec u = copa_poly_wrap(
                      cv_sub(cv_load(w + j), cv_load(z2 + j)), COPA_Q),
                  high = cv_load(w1 + j);
        copav_vec hint = copa_poly_center(
            copa_poly_wrap(cv_sub(high, copa_poly_high(u)),
                           COPA_Q / COPA_ALPHA),
            COPA_Q / COPA_ALPHA);
        cv_store(h + j, hint);
        copav_vec zt = copa_poly_center(
            copa_poly_wrap(cv_sub(cv_mul32(high, cv_set(COPA_ALPHA)), u),
                           COPA_Q),
            COPA_Q);
        sum = cv_add64(sum, copa_poly_square(zt));
    }
    return copa_poly_sum64(sum);
}
static inline uint64_t copa_backend_reconstruct(int32_t *w1,
                                                const int32_t *u,
                                                const int32_t *h,
                                                const uint8_t c[COPA_N])
{
    copav_vec sum = cv_set(0);
    for (unsigned r = 0; r < COPA_M; ++r)
        for (unsigned j = 0; j < COPA_N; j += COPAV_W) {
            unsigned offset = r * COPA_N + j;
            copav_vec value = cv_load(u + offset);
            if (r == 0) {
#if COPA_SIMD_BITS == 256
                copav_vec ch = _mm256_cvtepu8_epi32(
                    _mm_loadl_epi64((const __m128i *)(c + j)));
#else
                copav_vec ch = _mm512_cvtepu8_epi32(
                    _mm_loadu_si128((const __m128i *)(c + j)));
#endif
                value = cv_sub(value, cv_mul32(ch, cv_set(COPA_Q / 2)));
            }
            value = copa_poly_wrap(value, COPA_Q);
            copav_vec high = copa_poly_wrap(
                cv_add(copa_poly_high(value), cv_load(h + offset)),
                COPA_Q / COPA_ALPHA);
            cv_store(w1 + offset, high);
            copav_vec zt = copa_poly_center(
                copa_poly_wrap(
                    cv_sub(cv_mul32(high, cv_set(COPA_ALPHA)), value),
                    COPA_Q),
                COPA_Q);
            sum = cv_add64(sum, copa_poly_square(zt));
        }
    return copa_poly_sum64(sum);
}
/* Byte-key decoding shares the signed-square kernel with norm checks.
 * The frozen signed16 little-endian wire format is unchanged. */
static inline uint64_t copa_backend_unpack_s16(int32_t *out,const uint8_t *in,unsigned count){
 copav_vec sum=cv_set(0);
 for(unsigned j=0;j<count;j+=COPAV_W){
#if COPA_SIMD_BITS==256
  copav_vec x=_mm256_cvtepi16_epi32(_mm_loadu_si128((const __m128i *)(in+2*j)));
#else
  copav_vec x=_mm512_cvtepi16_epi32(_mm256_loadu_si256((const __m256i *)(in+2*j)));
#endif
  cv_store(out+j,x);sum=cv_add64(sum,copa_poly_square(x));
 }
 return copa_poly_sum64(sum);
}

#endif
