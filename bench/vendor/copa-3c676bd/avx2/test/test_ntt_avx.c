/*
 * test_ntt_avx.c - direct correctness test for the AVX2 NTT path.
 *
 * Exercises poly_mul (poly_ntt -> poly_pointwise_montgomery/basemul ->
 * poly_invntt_tomont, all routed to the in-tree AVX2 assembly)
 * and poly_invert (which must stay on the scalar path).  Compares poly_mul
 * against an independent schoolbook negacyclic product mod q, and checks
 * a*a^{-1}==1.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#include "params.h"
#include "polyvec.h"

static void schoolbook(const poly *a, const poly *b, poly *c)
{
    long long t[2 * N];
    for (int i = 0; i < 2 * N; i++)
        t[i] = 0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            t[i + j] += (long long)a->coeffs[i] * (long long)b->coeffs[j];
    for (int i = 0; i < N; i++) {
        long long v = t[i] - t[i + N]; /* x^N = -1 */
        v %= Q;
        if (v < 0)
            v += Q;
        c->coeffs[i] = (int32_t)v;
    }
}

static int32_t canon(int32_t v)
{
    v %= Q;
    if (v < 0)
        v += Q;
    return v;
}

int main(void)
{
    srand(20260602);
    int mul_fail = 0, inv_ok = 0, inv_fail = 0;

    for (int it = 0; it < 3000; it++) {
        poly a, b, c, cs;
        for (int i = 0; i < N; i++) {
            a.coeffs[i] = rand() % Q;
            b.coeffs[i] = rand() % Q;
        }
        poly_mul(&c, &a, &b); /* AVX2 path */
        schoolbook(&a, &b, &cs);
        for (int i = 0; i < N; i++) {
            if (canon(c.coeffs[i]) != cs.coeffs[i]) {
                if (mul_fail < 5)
                    printf("  mul mismatch it=%d i=%d got=%d exp=%d\n", it,
                           i, canon(c.coeffs[i]), cs.coeffs[i]);
                mul_fail++;
                break;
            }
        }
    }

    for (int it = 0; it < 1000; it++) {
        poly a, ai, one;
        for (int i = 0; i < N; i++)
            a.coeffs[i] = rand() % Q;
        if (poly_invert(&ai, &a) != 0)
            continue; /* non-invertible: skip */
        poly_mul(&one, &a, &ai);
        int ok = 1;
        for (int i = 0; i < N; i++)
            if (canon(one.coeffs[i]) != (i == 0 ? 1 : 0)) {
                ok = 0;
                break;
            }
        if (ok)
            inv_ok++;
        else
            inv_fail++;
    }

    int fails = mul_fail + inv_fail;
    printf(
        "%s test_ntt_avx mode=%d  poly_mul_fails=%d  inv_ok=%d "
        "inv_fail=%d\n",
        fails ? "FAIL" : "PASS", LITHIUM_MODE, mul_fail, inv_ok, inv_fail);
    return fails ? 1 : 0;
}
