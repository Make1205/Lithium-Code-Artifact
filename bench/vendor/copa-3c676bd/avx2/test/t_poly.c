/* t_poly.c - gate for IRS poly layer: NTT roundtrip, NTT-multiply vs
 * schoolbook negacyclic, and ring inversion (poly_mul(f, f^{-1}) == 1).
 * Single PASS/FAIL line + exit status. Deterministic. */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "ntt.h"
#include "params.h"
#include "polyvec.h"
#include "reduce.h"

static uint64_t st = 0xC0FFEE123456789ULL;
static uint32_t xr(void)
{
    st ^= st << 13;
    st ^= st >> 7;
    st ^= st << 17;
    return (uint32_t)(st >> 27);
}
static int32_t cmod(long long a)
{
    a %= Q;
    if (a < 0)
        a += Q;
    if (a > Q / 2)
        a -= Q;
    return (int32_t)a;
}
static uint16_t umod(long long a)
{
    a %= Q;
    if (a < 0)
        a += Q;
    return (uint16_t)a;
}
static void rp(poly *p)
{
    for (int i = 0; i < N; i++)
        p->coeffs[i] = cmod((long long)(xr() % Q));
}

static void schoolbook(poly *r, const poly *a, const poly *b)
{
    long long t[2 * N];
    for (int i = 0; i < 2 * N; i++)
        t[i] = 0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            t[i + j] += (long long)a->coeffs[i] * b->coeffs[j];
    for (int i = 0; i < N; i++)
        r->coeffs[i] = cmod(t[i] - t[i + N]);
}

int main(void)
{
    int fails = 0, inv_ok = 0, inv_fail = 0;
    for (int tr = 0; tr < 1500; tr++) {
        poly a, b, c, ref;
        rp(&a);
        rp(&b);
        /* roundtrip: poly_ntt + poly_invntt_tomont leaves one Montgomery
         * factor R; strip it via montgomery_reduce16. */
        poly16 ah;
        poly cc;
        poly_ntt(&ah, &a);
        poly_invntt_tomont(&cc, &ah);
        for (int i = 0; i < N; i++)
            if (montgomery_reduce16((uint16_t)cc.coeffs[i]) !=
                umod(a.coeffs[i])) {
                fails++;
                break;
            }
        /* multiply */
        poly_mul(&c, &a, &b);
        schoolbook(&ref, &a, &b);
        for (int i = 0; i < N; i++)
            if (cmod(c.coeffs[i]) != cmod(ref.coeffs[i])) {
                fails++;
                break;
            }
        /* inversion */
        poly fi, prod;
        if (poly_invert(&fi, &a) == 0) {
            inv_ok++;
            poly_mul(&prod, &a, &fi);
            /* prod should be the constant poly 1 */
            for (int i = 0; i < N; i++) {
                int32_t want = (i == 0) ? 1 : 0;
                if (cmod(prod.coeffs[i]) != want) {
                    fails++;
                    break;
                }
            }
        } else
            inv_fail++;
    }
    if (fails == 0 && inv_ok > 0) {
        printf(
            "PASS t_poly mode=%d n=%d q=%d base_deg=%d inv_ok=%d "
            "inv_fail=%d\n",
            LITHIUM_MODE, N, Q, BASE_DEG, inv_ok, inv_fail);
        return 0;
    }
    printf("FAIL t_poly mode=%d fails=%d inv_ok=%d\n", LITHIUM_MODE, fails,
           inv_ok);
    return 1;
}
