#include "poly_inv.h"
#include "copa_scalar64.h"
#include <assert.h>

static void mul_mod(int32_t *out, const int32_t *a, const int32_t *b,
                    size_t n, uint32_t q)
{
    int64_t tmp[COPA_N];
    copa_ntt_crt_mul(tmp, a, b, n);
    for (size_t i = 0; i < n; ++i) out[i] = (int32_t)copa_mod_i64(tmp[i], q);
}
static int inverse(int32_t *out, const int32_t *f, size_t n, uint32_t q,
                   uint32_t phi)
{
    if (n == 1) {
        uint32_t a = copa_mod_i64(f[0], q), acc = 1;
        /* Euler exponent is public, multiplication inputs remain secret. */
        for (int bit = 31; bit >= 0; --bit) {
            acc = copa_mod_i64((int64_t)acc * acc, q);
            if (((phi-1) >> bit) & 1) acc = copa_mod_i64((int64_t)acc * a, q);
        }
        out[0] = (int32_t)acc;
        return copa_mod_i64((int64_t)a * acc, q) == 1;
    }
    size_t h = n/2;
    int32_t a[COPA_N/2] = {0}, b[COPA_N/2] = {0}, aa[COPA_N/2] = {0}, bb[COPA_N/2] = {0};
    int32_t norm[COPA_N/2], g[COPA_N/2];
    for (size_t i = 0; i < h; ++i) {
        a[i] = (int32_t)copa_mod_i64(f[2*i], q);
        b[i] = (int32_t)copa_mod_i64(f[2*i+1], q);
    }
    mul_mod(aa, a, a, h, q); mul_mod(bb, b, b, h, q);
    norm[0] = (int32_t)copa_mod_i64((int64_t)aa[0]+bb[h-1], q);
    for (size_t i = 1; i < h; ++i)
        norm[i] = (int32_t)copa_mod_i64((int64_t)aa[i]-bb[i-1], q);
    int ok = inverse(g, norm, h, q, phi);
    mul_mod(aa, a, g, h, q); mul_mod(bb, b, g, h, q);
    for (size_t i = 0; i < h; ++i) {
        out[2*i] = aa[i]; out[2*i+1] = (int32_t)copa_mod_i64(-(int64_t)bb[i], q);
    }
    return ok;
}
int copa_poly_inv(int32_t out[COPA_N], const int32_t f[COPA_N], uint32_t q)
{
    assert(q == 7680 || q == 12584 || q == 6048 || q == 8400);
    uint32_t phi = q == 7680 ? 2048 : q == 12584 ? 5280 : q == 6048 ? 1728 : 1920;
    int32_t result[COPA_N], check[COPA_N];
    uint64_t ok_mask = copa64_ct_mask((unsigned)inverse(result, f, COPA_N, q, phi));
    mul_mod(check, result, f, COPA_N, q);
    uint32_t bad = (uint32_t)(check[0] ^ 1);
    for (size_t i = 1; i < COPA_N; ++i) bad |= (uint32_t)check[i];
    uint32_t mask = (uint32_t)(ok_mask & copa64_ct_mask(bad == 0));
    for (size_t i = 0; i < COPA_N; ++i) out[i] = (int32_t)((uint32_t)result[i] & mask);
    return (int)(mask & 1);
}
