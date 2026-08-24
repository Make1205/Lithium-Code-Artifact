/*
 * speed_exp_avx512ifma.c - microbenchmark + correctness gate for the
 * AVX-512 IFMA integer exp(-a) kernel.
 *
 * Correctness (gate, return code): approx_exp_avx512ifma(a_q60) stays
 * below 2^-43 relative error versus libm exp(-a) on [0, AMAX].
 *
 * Timing (informational): median cycles per output for W8/W16/W32 batches
 * (1/2/4 independent 8-lane streams) plus the scalar Q63 baseline.  The
 * inputs are independent and drawn from a small L1-resident pool, matching
 * the style of ref/test/speed_approx_exp.c and avx2/test/speed_exp_avx2.c.
 *
 * Build: make AVX512IFMAEXP=1 TESTS=speed_exp_avx512ifma
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "approx_exp.h"
#include "config.h"
#include "cpucycles.h"

#define AMAX 7.369
#define NTESTS 2048
#define KBATCH 256
#define POOL 512

static uint64_t s_state;
static uint64_t sm64(void)
{
    uint64_t z = (s_state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static uint64_t rand_a_q60(void)
{
    const uint64_t max_q60 = (uint64_t)(AMAX * 0x1p60);
    return sm64() % (max_q60 + 1u);
}

static double accuracy_ifma(double *worst_a)
{
    const int niter = 40000;
    double worst = 0.0;
    *worst_a = 0.0;
    for (int i = 0; i <= niter; i++) {
        double a = AMAX * i / niter;
        uint64_t in[8], out[8];
        uint64_t q = (uint64_t)(a * 0x1p60 + 0.5);
        for (int k = 0; k < 8; k++)
            in[k] = q;
        approx_exp_avx512ifma(out, in, 8);
        double want = exp(-a);
        for (int k = 0; k < 8; k++) {
            double got = (double)out[k] / (double)(UINT64_C(1) << 63);
            double rel = fabs(got - want) / want;
            if (rel > worst) {
                worst = rel;
                *worst_a = a;
            }
        }
    }
    return worst;
}

static uint64_t bench_scalar(const uint64_t pool[POOL][32],
                             volatile uint64_t *sink)
{
    static uint64_t t[NTESTS];
    for (int s = 0; s < NTESTS; s++) {
        uint64_t a = cpucycles();
        for (int k = 0; k < KBATCH; k++) {
            const uint64_t *in = pool[(s * KBATCH + k) & (POOL - 1)];
            for (int i = 0; i < 8; i++)
                *sink ^= approx_exp(in[i]);
        }
        t[s] = cpucycles() - a;
    }
    return cpucycles_median(t, NTESTS);
}

static uint64_t bench_w8(const uint64_t pool[POOL][32],
                         volatile uint64_t *sink)
{
    static uint64_t t[NTESTS];
    for (int s = 0; s < NTESTS; s++) {
        uint64_t a = cpucycles();
        for (int k = 0; k < KBATCH; k++) {
            uint64_t out[8];
            approx_exp_avx512ifma(out, pool[(s * KBATCH + k) & (POOL - 1)],
                                  8);
            for (int i = 0; i < 8; i++)
                *sink ^= out[i];
        }
        t[s] = cpucycles() - a;
    }
    return cpucycles_median(t, NTESTS);
}

static uint64_t bench_w16(const uint64_t pool[POOL][32],
                          volatile uint64_t *sink)
{
    static uint64_t t[NTESTS];
    for (int s = 0; s < NTESTS; s++) {
        uint64_t a = cpucycles();
        for (int k = 0; k < KBATCH; k++) {
            uint64_t out[16];
            approx_exp_avx512ifma(out, pool[(s * KBATCH + k) & (POOL - 1)],
                                  16);
            for (int i = 0; i < 16; i++)
                *sink ^= out[i];
        }
        t[s] = cpucycles() - a;
    }
    return cpucycles_median(t, NTESTS);
}

static uint64_t bench_w32(const uint64_t pool[POOL][32],
                          volatile uint64_t *sink)
{
    static uint64_t t[NTESTS];
    for (int s = 0; s < NTESTS; s++) {
        uint64_t a = cpucycles();
        for (int k = 0; k < KBATCH; k++) {
            uint64_t out[32];
            approx_exp_avx512ifma(out, pool[(s * KBATCH + k) & (POOL - 1)],
                                  32);
            for (int i = 0; i < 32; i++)
                *sink ^= out[i];
        }
        t[s] = cpucycles() - a;
    }
    return cpucycles_median(t, NTESTS);
}

int main(void)
{
    static uint64_t pool[POOL][32];
    volatile uint64_t sink = 0;
    int ok = 1;

    double worst_a = 0.0;
    double rel = accuracy_ifma(&worst_a);
    int acc_ok = (rel < 0x1p-43);
    if (!acc_ok)
        ok = 0;

    s_state = 0x123456789ABCDEF0ULL;
    for (int i = 0; i < POOL; i++)
        for (int j = 0; j < 32; j++)
            pool[i][j] = rand_a_q60();

    for (int i = 0; i < POOL; i++) {
        uint64_t out[32];
        approx_exp_avx512ifma(out, pool[i], 32);
        sink ^= out[0];
    }

    uint64_t medScalar = bench_scalar(pool, &sink);
    uint64_t med8 = bench_w8(pool, &sink);
    uint64_t med16 = bench_w16(pool, &sink);
    uint64_t med32 = bench_w32(pool, &sink);

    double cScalar = (double)medScalar / KBATCH / 8.0;
    double c8 = (double)med8 / KBATCH / 8.0;
    double c16 = (double)med16 / KBATCH / 16.0;
    double c32 = (double)med32 / KBATCH / 32.0;
    printf("scalar Q63 approx_exp          : %.2f cyc/out\n", cScalar);
    printf("W8  AVX512IFMA exp (1 stream ) : %.2f cyc/out\n", c8);
    printf("W16 AVX512IFMA exp (2 streams) : %.2f cyc/out\n", c16);
    printf("W32 AVX512IFMA exp (4 streams) : %.2f cyc/out\n", c32);
    printf("speedup W32 vs scalar x%.2f   (W32 vs W8 x%.2f)\n",
           cScalar / c32, c8 / c32);
    printf(
        "%s speed_exp_avx512ifma mode=%d rel[0,%.3f]=%.3e at a=%.4f "
        "(gate 2^-43) sink=%llu\n",
        ok ? "PASS" : "FAIL", LITHIUM_MODE, AMAX, rel, worst_a,
        (unsigned long long)sink);
    return ok ? 0 : 1;
}
