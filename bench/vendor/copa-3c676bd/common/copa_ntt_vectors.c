/* Exact unsigned-Montgomery negacyclic NTTs, shared AVX2/AVX-512 schedule.
 */
#include <assert.h>
#include <string.h>

#include "copa_scalar64.h"
#include "copa_simd.h"
#include "ntt_crt.h"
#define P1 UINT32_C(8380417)
#define P2 UINT32_C(7340033)
typedef struct {
    uint32_t p, qinv;
    const uint32_t *twist, *post, *forward, *inverse, *twist_shoup,
        *post_shoup, *forward_shoup, *inverse_shoup, *zeta, *izeta,
        *zeta_shoup, *izeta_shoup, *scale, *scale_shoup;
} ntt_prime;
#include "copa_ntt_cache_tables.h"
#include "copa_simd_ntt_tables.h"
static inline uint32_t scalar_reduce(uint32_t a, uint32_t p)
{
    return a - (p & (uint32_t)copa64_ct_mask(a >= p));
}
static inline uint32_t scalar_shoup(uint32_t a, uint32_t w, uint32_t ws,
                                    uint32_t p)
{
    uint32_t q = (uint32_t)(((uint64_t)a * ws) >> 32);
    return scalar_reduce(a * w - q * p, p);
}
static inline uint32_t scalar_mont(uint32_t a, uint32_t b,
                                   const ntt_prime *p)
{
    uint64_t z = (uint64_t)a * b;
    uint32_t t = (uint32_t)z * p->qinv;
    return scalar_reduce((uint32_t)((z + (uint64_t)t * p->p) >> 32), p->p);
}
static void scalar_stages(uint32_t *x, size_t n, const ntt_prime *p,
                          unsigned inverse)
{
    for (size_t stage = 0; stage < (size_t)__builtin_ctz((unsigned)n);
         ++stage) {
        size_t len = inverse ? ((size_t)1 << stage) : n >> (stage + 1),
               root = n / (2 * len);
        for (size_t b = 0; b < n; b += 2 * len, ++root) {
            uint32_t w = (inverse ? p->izeta : p->zeta)[n - 1 + root],
                     ws = (inverse ? p->izeta_shoup
                                   : p->zeta_shoup)[n - 1 + root];
            for (size_t j = b; j < b + len; ++j) {
                uint32_t a = x[j], v = x[j + len];
                if (inverse) {
                    x[j] = scalar_reduce(a + v, p->p);
                    x[j + len] = scalar_shoup(
                        scalar_reduce(a + p->p - v, p->p), w, ws, p->p);
                } else {
                    v = scalar_shoup(v, w, ws, p->p);
                    x[j] = scalar_reduce(a + v, p->p);
                    x[j + len] = scalar_reduce(a + p->p - v, p->p);
                }
            }
        }
    }
}
/* For a<256*p<2^31, Shoup's unreduced remainder is <2*p.
 * Forward bounds grow by 2*p per level; inverse sums double per level. */
static inline copav_vec ntt_lazy_shoup(copav_vec a,copav_vec w,copav_vec ws,uint32_t p){return cv_sub(cv_mul32(a,w),cv_mul32(cv_high32(a,ws),cv_set(p)));}
static inline copav_vec ntt_canonical(copav_vec a,uint32_t p){copav_vec q=cv_high32(a,cv_set((uint32_t)(UINT64_C(0x100000000)/p)));return cv_reduce(cv_sub(a,cv_mul32(q,cv_set(p))),p);}
static void small_stage(uint32_t *x, size_t n, unsigned prime,
                        unsigned stage, unsigned inverse)
{
    unsigned len = 1u << stage;
    uint32_t index[COPAV_W] __attribute__((aligned(64)));
    uint32_t upper[COPAV_W] __attribute__((aligned(64)));
    for (unsigned j = 0; j < COPAV_W; ++j) {
        index[j] = j ^ len;
        upper[j] = 0u - (uint32_t)((j & len) != 0);
    }
    copav_vec ix = cv_load(index), mask = cv_load(upper);
    const ntt_prime *p = &ntt_primes[prime];
    const uint32_t *roots = (inverse ? ntv_inverse[prime][stage]
                                     : ntv_forward[prime][stage]) +
                            n - 1,
                   *shoup = (inverse ? ntv_inverse_shoup[prime][stage]
                                     : ntv_forward_shoup[prime][stage]) +
                            n - 1;
    for (size_t j = 0; j < n; j += COPAV_W) {
        copav_vec v = cv_load(x + j), partner = cv_perm(v, ix),
                  a = cv_select(partner, v, mask),
                  b = cv_select(v, partner, mask), w = cv_load(roots + j),
                  ws = cv_load(shoup + j), sum, diff;
        if (inverse) {
            sum = cv_add(a,b);
            diff = ntt_lazy_shoup(
                cv_sub(cv_add(a,cv_set((1u<<stage)*p->p)),b), w, ws,
                p->p);
        } else {
            b = ntt_lazy_shoup(b, w, ws, p->p);
            sum = cv_add(a,b);
            diff = cv_sub(cv_add(a,cv_set(2*p->p)),b);
        }
        cv_store(x + j, cv_select(diff, sum, mask));
    }
}
#include "copa_ntt_zip.h"
#include "copa_ntt_tiles.h"
static void vector_stages(uint32_t *x, size_t n, unsigned prime,
                          unsigned inverse)
{
    const ntt_prime *p = &ntt_primes[prime];
    unsigned logn = (unsigned)__builtin_ctz((unsigned)n);
    if (n < COPAV_W) {
        scalar_stages(x, n, p, inverse);
        return;
    }
    const size_t tile=4*COPAV_W;
    if(inverse&&n>=tile)for(size_t base=0;base<n;base+=tile)tile_inverse(x,n,base,prime);
    for (unsigned stage = 0; stage < logn; ++stage) {
        size_t len = inverse ? ((size_t)1 << stage) : n >> (stage + 1);
        if(n>=tile&&len<tile){if(!inverse){for(size_t base=0;base<n;base+=tile)tile_forward(x,n,base,prime);break;}continue;}
        if (len < COPAV_W) {
            small_stage(x, n, prime,
                        (unsigned)__builtin_ctz((unsigned)len), inverse);
            continue;
        }
        size_t root = n / (2 * len);
        for (size_t b = 0; b < n; b += 2 * len, ++root) {
            copav_vec w = cv_set(
                          (inverse ? p->izeta : p->zeta)[n - 1 + root]),
                      ws = cv_set((inverse ? p->izeta_shoup
                                           : p->zeta_shoup)[n - 1 + root]);
            for (size_t j = b; j < b + len; j += COPAV_W) {
                copav_vec a = cv_load(x + j), v = cv_load(x + j + len),
                          sum, diff;
                if (inverse) {
                    sum = cv_add(a,v);
                    diff = ntt_lazy_shoup(
                        cv_sub(cv_add(a,cv_set((1u<<stage)*p->p)),v),
                        w, ws, p->p);
                } else {
                    v = ntt_lazy_shoup(v, w, ws, p->p);
                    sum = cv_add(a,v);
                    diff = cv_sub(cv_add(a,cv_set(2*p->p)),v);
                }
                cv_store(x + j, sum);
                cv_store(x + j + len, diff);
            }
        }
    }
}
void copa_ntt_forward(uint32_t out[2][COPA_N], const int32_t *in, size_t n)
{
    assert(n && n <= COPA_N && !(n & (n - 1)));
    for (unsigned k = 0; k < 2; ++k) {
        const ntt_prime *p = &ntt_primes[k];
        size_t j = 0;
        for (; j + COPAV_W <= n; j += COPAV_W) {
            copav_vec v = cv_load(in + j);
            cv_store(out[k] + j,
                     cv_add(v, cv_and(cv_sra32(v, 31), cv_set(p->p))));
        }
        for (; j < n; ++j)
            out[k][j] = (uint32_t)in[j] +
                        (p->p & (uint32_t)copa64_ct_mask(in[j] < 0));
        vector_stages(out[k], n, k, 0);
        for(size_t t=0;t+COPAV_W<=n;t+=COPAV_W)cv_store(out[k]+t,ntt_canonical(cv_load(out[k]+t),p->p));
    }
}
static void dot(uint32_t out[2][COPA_N], const uint32_t a[][2][COPA_N],
                const uint32_t b[][2][COPA_N], unsigned count, size_t n,
                unsigned reversed)
{
    uint32_t reverse[COPAV_W] __attribute__((aligned(64)));
    for (unsigned j = 0; j < COPAV_W; ++j)
        reverse[j] = COPAV_W - 1 - j;
    copav_vec rev = cv_load(reverse);
    for (unsigned k = 0; k < 2; ++k) {
        const ntt_prime *p = &ntt_primes[k];
        size_t j = 0;
        for (; j + COPAV_W <= n; j += COPAV_W) {
            copav_vec sum = cv_set(0);
            for (unsigned r = 0; r < count; ++r) {
                copav_vec y =
                    reversed
                        ? cv_perm(cv_load(b[r][k] + n - COPAV_W - j), rev)
                        : cv_load(b[r][k] + j);
                sum = cv_reduce(cv_add(sum, cv_mont(cv_load(a[r][k] + j),
                                                    y, p->p, p->qinv)),
                                p->p);
            }
            cv_store(out[k] + j, sum);
        }
        for (; j < n; ++j) {
            uint32_t sum = 0;
            for (unsigned r = 0; r < count; ++r)
                sum = scalar_reduce(
                    sum + scalar_mont(a[r][k][j],
                                      b[r][k][reversed ? n - 1 - j : j],
                                      p),
                    p->p);
            out[k][j] = sum;
        }
    }
}
void copa_ntt_dot(uint32_t out[2][COPA_N], const uint32_t a[][2][COPA_N],
                  const uint32_t b[][2][COPA_N], unsigned count, size_t n)
{
    dot(out, a, b, count, n, 0);
}
void copa_ntt_dot_reversed(uint32_t out[2][COPA_N],
                           const uint32_t a[][2][COPA_N],
                           const uint32_t b[][2][COPA_N], unsigned count,
                           size_t n)
{
    dot(out, a, b, count, n, 1);
}
static void inverse_residues(uint32_t v[2][COPA_N],
                             const uint32_t in[2][COPA_N], size_t n)
{
    assert(n && n <= COPA_N && !(n & (n - 1)));
    unsigned logn = (unsigned)__builtin_ctz((unsigned)n);
    for (unsigned k = 0; k < 2; ++k) {
        const ntt_prime *p = &ntt_primes[k];
        memcpy(v[k], in[k], n * sizeof(uint32_t));
        vector_stages(v[k], n, k, 1);
        size_t j = 0;
        copav_vec w = cv_set(p->scale[logn]),
                  ws = cv_set(p->scale_shoup[logn]);
        for (; j + COPAV_W <= n; j += COPAV_W)
            cv_store(v[k] + j, cv_shoup(cv_load(v[k] + j), w, ws, p->p));
        for (; j < n; ++j)
            v[k][j] = scalar_shoup(v[k][j], p->scale[logn],
                                   p->scale_shoup[logn], p->p);
    }
}
static inline copav_vec crt_multiplier(copav_vec a, copav_vec b)
{
    copav_vec d = cv_sub(b, a);
    d = cv_add(d, cv_and(cv_sra32(d, 31), cv_set(2 * P2)));
    return cv_mont(cv_reduce(d, P2), cv_set(NTT_P1_INV_MONT_P2), P2,
                   ntt_primes[1].qinv);
}
static inline uint32_t scalar_crt_t(uint32_t a, uint32_t b)
{
    int64_t d = (int64_t)b - a;
    uint32_t u =
        (uint32_t)(d + (int64_t)((2 * P2) & copa64_ct_mask(d < 0)));
    return scalar_mont(scalar_reduce(u, P2), NTT_P1_INV_MONT_P2,
                       &ntt_primes[1]);
}
void copa_ntt_inverse_crt(int64_t *out, const uint32_t in[2][COPA_N],
                          size_t n)
{
    uint32_t v[2][COPA_N];
    inverse_residues(v, in, n);
    const uint64_t product = (uint64_t)P1 * P2;
    size_t j = 0;
    for (; j + COPAV_W <= n; j += COPAV_W) {
        copav_vec a = cv_load(v[0] + j),
                  t = crt_multiplier(a, cv_load(v[1] + j));
        copav_vec even =
                      cv_add64(cv_mul_even(t, cv_set(P1)),
                               cv_and(a, cv_set64(UINT64_C(0xffffffff)))),
                  odd = cv_add64(cv_mul_even(cv_srl64(t, 32), cv_set(P1)),
                                 cv_srl64(a, 32));
        even =
            cv_sub64(even, cv_and(cv_set64(product),
                                  cv_gt64(even, cv_set64(product / 2))));
        odd = cv_sub64(odd, cv_and(cv_set64(product),
                                   cv_gt64(odd, cv_set64(product / 2))));
#if COPA_SIMD_BITS == 256
        __m256i lo = _mm256_unpacklo_epi64(even, odd),
                hi = _mm256_unpackhi_epi64(even, odd);
        cv_store(out + j, _mm256_permute2x128_si256(lo, hi, 0x20));
        cv_store(out + j + 4, _mm256_permute2x128_si256(lo, hi, 0x31));
#else
        const __m512i lowindex =
                          _mm512_setr_epi64(0, 8, 1, 9, 2, 10, 3, 11),
                      highindex =
                          _mm512_setr_epi64(4, 12, 5, 13, 6, 14, 7, 15);
        cv_store(out + j, _mm512_permutex2var_epi64(even, lowindex, odd));
        cv_store(out + j + 8,
                 _mm512_permutex2var_epi64(even, highindex, odd));
#endif
    }
    for (; j < n; ++j) {
        uint64_t x =
            v[0][j] + (uint64_t)P1 * scalar_crt_t(v[0][j], v[1][j]);
        out[j] = (int64_t)x -
                 (int64_t)(product & copa64_ct_mask(x > product / 2));
    }
}
static inline copav_vec reduce_q(copav_vec a, uint32_t q,
                                 uint32_t reciprocal)
{
    return cv_reduce(
        cv_sub(a, cv_mul32(cv_high32(a, cv_set(reciprocal)), cv_set(q))),
        q);
}
void copa_ntt_inverse_modq(int32_t *out, const uint32_t in[2][COPA_N],
                           size_t n, uint32_t q)
{
    uint32_t v[2][COPA_N];
    inverse_residues(v, in, n);
    const uint64_t product = (uint64_t)P1 * P2;
    uint32_t reciprocal = (uint32_t)(UINT64_C(0x100000000) / q),
             p1q = P1 % q, pq = (uint32_t)(product % q);
    size_t j = 0;
    for (; j + COPAV_W <= n; j += COPAV_W) {
        copav_vec a = cv_load(v[0] + j),
                  t = crt_multiplier(a, cv_load(v[1] + j));
        copav_vec negative = cv_or(
            cv_gt(t, cv_set(P2 / 2)),
            cv_and(cv_eq(t, cv_set(P2 / 2)), cv_gt(a, cv_set(P1 / 2))));
        copav_vec r = reduce_q(
            cv_add(reduce_q(a, q, reciprocal),
                   cv_mul32(cv_set(p1q), reduce_q(t, q, reciprocal))),
            q, reciprocal);
        r = cv_sub(r, cv_and(cv_set(pq), negative));
        r = cv_add(r, cv_and(cv_set(q), cv_gt(cv_set(0), r)));
        cv_store(out + j, r);
    }
    for (; j < n; ++j) {
        uint32_t t = scalar_crt_t(v[0][j], v[1][j]);
        uint64_t negative = copa64_ct_mask(t > P2 / 2) |
                            (copa64_ct_mask(t == P2 / 2) &
                             copa64_ct_mask(v[0][j] > P1 / 2));
        uint32_t a = scalar_reduce(
                     v[0][j] -
                         (uint32_t)(((uint64_t)v[0][j] * reciprocal) >>
                                    32) *
                             q,
                     q),
                 b = scalar_reduce(
                     t - (uint32_t)(((uint64_t)t * reciprocal) >> 32) * q,
                     q),
                 r = a + p1q * b;
        r = scalar_reduce(
            r - (uint32_t)(((uint64_t)r * reciprocal) >> 32) * q, q);
        int32_t z = (int32_t)r - (int32_t)(pq & negative);
        out[j] = z + (int32_t)(q & (uint32_t)copa64_ct_mask(z < 0));
    }
}
void copa_ntt_crt_mul(int64_t *out, const int32_t *a, const int32_t *b,
                      size_t n)
{
    uint32_t aa[1][2][COPA_N], bb[1][2][COPA_N], cc[2][COPA_N];
    copa_ntt_forward(aa[0], a, n);
    copa_ntt_forward(bb[0], b, n);
    copa_ntt_dot(cc, aa, bb, 1, n);
    copa_ntt_inverse_crt(out, cc, n);
}
