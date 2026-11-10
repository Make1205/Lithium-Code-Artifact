/* test_ntt.c - self-consistency gate for IRS NTT (full + incomplete),
 * uint16.
 *
 *   (1) invntt16_tomont(ntt16(a)) == a*R  (roundtrip; strip one R via
 * mreduce16) (2) ntt-based multiply == schoolbook negacyclic multiply mod
 * q
 *
 * Deterministic xorshift PRNG.  Prints a single PASS/FAIL line + exit
 * status.
 */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "ntt.h"
#include "params.h"
#include "reduce.h"

static uint64_t st = 0x123456789abcdef0ULL;
static uint32_t xr(void)
{
    st ^= st << 13;
    st ^= st >> 7;
    st ^= st << 17;
    return (uint32_t)(st >> 27);
}

static uint16_t umod(long long a)
{
    a %= Q;
    if (a < 0)
        a += Q;
    return (uint16_t)a;
}

static void rand_poly(uint16_t a[N])
{
    for (int i = 0; i < N; i++)
        a[i] = (uint16_t)(xr() % Q); /* [0,q) */
}

static void schoolbook(uint16_t r[N], const uint16_t a[N],
                       const uint16_t b[N])
{
    long long t[2 * N];
    for (int i = 0; i < 2 * N; i++)
        t[i] = 0;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            t[i + j] += (long long)a[i] * b[j];
    for (int i = 0; i < N; i++)
        r[i] = umod(t[i] - t[i + N]);
}

int main(void)
{
    int fails = 0;
    for (int trial = 0; trial < 2000; trial++) {
        uint16_t a[N], b[N], ah[N], bh[N], ch[N], ref[N];
        rand_poly(a);
        rand_poly(b);

        /* ntt16 + invntt16_tomont leaves one Montgomery factor R (a*R mod
         * q); strip it with montgomery_reduce16 to recover a. */
        memcpy(ah, a, sizeof a);
        ntt16(ah);
        invntt16_tomont(ah);
        for (int i = 0; i < N; i++)
            if (montgomery_reduce16(ah[i]) != a[i]) {
                fails++;
                break;
            }

        memcpy(ah, a, sizeof a);
        ntt16(ah);
        memcpy(bh, b, sizeof b);
        ntt16(bh);
#if BASE_DEG == 1
        for (int i = 0; i < N; i++)
            ch[i] = fqmul16(ah[i], bh[i]);
#else
        poly16_basemul_montgomery_native(ch, ah, bh);
#endif
        invntt16_tomont(ch);
        schoolbook(ref, a, b);
        for (int i = 0; i < N; i++)
            if (ch[i] != ref[i]) {
                fails++;
                break;
            }
    }
    if (fails == 0) {
        printf("PASS test_ntt mode=%d n=%d q=%d base_deg=%d\n",
               LITHIUM_MODE, N, Q, BASE_DEG);
        return 0;
    }
    printf("FAIL test_ntt mode=%d fails=%d\n", LITHIUM_MODE, fails);
    return 1;
}
