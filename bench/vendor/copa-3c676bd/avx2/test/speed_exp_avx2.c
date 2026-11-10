/*
 * speed_exp_avx2.c - microbenchmark + correctness gate for the AVX2
 * binary64 exp(-a) kernels.
 *
 * Correctness (gate, return code): W4/W8/W16 all stay below 2^-43
 * relative error versus libm exp(-a) on [0, AMAX].  The production
 * Q60->Q63 wrapper is also checked against the raw W4 kernel's bit-level
 * conversion path.
 *
 * Timing (informational): median cycles per output for the three raw
 * kernels.  The inputs are independent and drawn from a small L1-resident
 * pool, matching the style of ref/test/speed_approx_exp.c.
 *
 * Build: make AVX2EXP=1 TESTS=speed_exp_avx2
 */
#include <immintrin.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "approx_exp_avx2.h"
#include "config.h"
#include "cpucycles.h"

#define AMAX 7.369
#define NCORR 1024
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

static double rand_a(void)
{
    return (double)(sm64() >> 11) * 0x1p-53 * AMAX;
}

static uint64_t rand_a_q60(void)
{
    const uint64_t max_q60 = (uint64_t)(AMAX * 0x1p60);
    return sm64() % (max_q60 + 1u);
}

static uint64_t q63_from_double_bits(double exp_neg)
{
    double x = exp_neg * APPROX_EXP_AVX2_SCALE63 + 0.5;
    uint64_t bits;
    memcpy(&bits, &x, sizeof(bits));
    uint64_t exp = (bits >> 52) & 0x7FFu;
    uint64_t mant = (bits & 0x000FFFFFFFFFFFFFULL) | 0x0010000000000000ULL;

    if (exp > 1023u + 51u)
        return mant << (exp - (1023u + 52u));
    return mant >> ((1023u + 52u) - exp);
}

static int check_q60_q63_wrapper(void)
{
    uint64_t in[NCORR], out[NCORR];

    s_state = 0x5152535455565758ULL;
    for (int i = 0; i < NCORR; i++)
        in[i] = rand_a_q60();
    in[0] = 0;
    in[1] = 1;
    in[2] = UINT64_C(1) << 60;
    in[3] = UINT64_C(1) << 63;

    approx_exp_avx2(out, in, NCORR);
    for (int i = 0; i < NCORR; i++) {
        double a = (double)in[i] * 0x1p-60;
        double aa[4] = {a, a, a, a};
        double e[4];
        exp_avx2_w4(aa, e);
        uint64_t want = q63_from_double_bits(e[0]);
        if (out[i] != want) {
            printf("WRAPPERDIFF i=%d in=%llu got=%llu want=%llu\n", i,
                   (unsigned long long)in[i], (unsigned long long)out[i],
                   (unsigned long long)want);
            return 0;
        }
    }
    return 1;
}

static double accuracy_width(void (*kernel)(const double *, double *),
                             int width)
{
    const int niter = 40000;
    double worst = 0.0;
    for (int i = 0; i <= niter; i++) {
        double a = AMAX * i / niter;
        double in[16], out[16];
        for (int k = 0; k < width; k++)
            in[k] = a;
        kernel(in, out);
        double want = exp(-a);
        for (int k = 0; k < width; k++) {
            double rel = fabs(out[k] - want) / want;
            if (rel > worst)
                worst = rel;
        }
    }
    return worst;
}

static uint64_t bench_w4(const double pool[POOL][16],
                         volatile double *sink)
{
    static uint64_t t[NTESTS];
    for (int s = 0; s < NTESTS; s++) {
        uint64_t a = cpucycles();
        for (int k = 0; k < KBATCH; k++) {
            double out[4];
            exp_avx2_w4(pool[(s * KBATCH + k) & (POOL - 1)], out);
            *sink += out[0] + out[1] + out[2] + out[3];
        }
        t[s] = cpucycles() - a;
    }
    return cpucycles_median(t, NTESTS);
}

static uint64_t bench_w8(const double pool[POOL][16],
                         volatile double *sink)
{
    static uint64_t t[NTESTS];
    for (int s = 0; s < NTESTS; s++) {
        uint64_t a = cpucycles();
        for (int k = 0; k < KBATCH; k++) {
            double out[8];
            exp_avx2_w8(pool[(s * KBATCH + k) & (POOL - 1)], out);
            *sink += out[0] + out[1] + out[2] + out[3] + out[4] + out[5] +
                     out[6] + out[7];
        }
        t[s] = cpucycles() - a;
    }
    return cpucycles_median(t, NTESTS);
}

static uint64_t bench_w16(const double pool[POOL][16],
                          volatile double *sink)
{
    static uint64_t t[NTESTS];
    for (int s = 0; s < NTESTS; s++) {
        uint64_t a = cpucycles();
        for (int k = 0; k < KBATCH; k++) {
            double out[16];
            exp_avx2_w16(pool[(s * KBATCH + k) & (POOL - 1)], out);
            *sink += out[0] + out[1] + out[2] + out[3] + out[4] + out[5] +
                     out[6] + out[7] + out[8] + out[9] + out[10] +
                     out[11] + out[12] + out[13] + out[14] + out[15];
        }
        t[s] = cpucycles() - a;
    }
    return cpucycles_median(t, NTESTS);
}

int main(void)
{
    static double pool[POOL][16];
    volatile double sink = 0.0;
    int ok = 1;

    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);

    double rel4 = accuracy_width(exp_avx2_w4, 4);
    double rel8 = accuracy_width(exp_avx2_w8, 8);
    double rel16 = accuracy_width(exp_avx2_w16, 16);
    int acc_ok = (rel4 < 0x1p-43) && (rel8 < 0x1p-43) && (rel16 < 0x1p-43);
    int wrap_ok = check_q60_q63_wrapper();
    if (!acc_ok || !wrap_ok)
        ok = 0;

    s_state = 0x123456789ABCDEF0ULL;
    for (int i = 0; i < POOL; i++)
        for (int j = 0; j < 16; j++)
            pool[i][j] = rand_a();

    for (int i = 0; i < POOL; i++) {
        double out[16];
        exp_avx2_w16(pool[i], out);
        sink += out[0];
    }

    uint64_t med4 = bench_w4(pool, &sink);
    uint64_t med8 = bench_w8(pool, &sink);
    uint64_t med16 = bench_w16(pool, &sink);

    double c4 = (double)med4 / KBATCH / 4.0;
    double c8 = (double)med8 / KBATCH / 8.0;
    double c16 = (double)med16 / KBATCH / 16.0;
    printf("W4  AVX2 exp (1 stream )        : %.2f cyc/out\n", c4);
    printf("W8  AVX2 exp (2 streams)        : %.2f cyc/out\n", c8);
    printf("W16 AVX2 exp (4 streams)        : %.2f cyc/out\n", c16);
    printf("speedup W16 vs W4 x%.2f   (W16 vs W8 x%.2f)\n", c4 / c16,
           c8 / c16);
    printf(
        "%s speed_exp_avx2 mode=%d rel[0,%.3f]={W4 %.3e, W8 %.3e, W16 "
        "%.3e} wrapper=%s (gate 2^-43) sink=%.3g\n",
        ok ? "PASS" : "FAIL", LITHIUM_MODE, AMAX, rel4, rel8, rel16,
        wrap_ok ? "OK" : "FAIL", (double)sink);
    return ok ? 0 : 1;
}
