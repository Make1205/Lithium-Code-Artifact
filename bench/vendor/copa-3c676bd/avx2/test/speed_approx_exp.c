/* speed_approx_exp.c - microbenchmark + correctness gate for
 * approx_exp_4x.
 *
 * Correctness (gate, return code):
 *   (1) bit-identity : approx_exp_4x lane[l] == approx_exp(a[l]) for a
 * large reproducible random set of a_q60 in [0, 15*2^60)
 * (Q60-representable); (2) accuracy     : approx_exp within 2^-43 relative
 * of libm exp(-a) over the sampler range [0, 8].
 *
 * Timing (informational): median cycles per output for
 *   - 1x latency-bound  : a serial chain (each call's input = previous
 * output), i.e. the cost of an ISOLATED approx_exp on the critical path --
 * the scalar Horner + 32-entry OR-scan run as one long dependency chain;
 *   - 1x throughput     : 4 *independent* scalar calls in a loop -- here
 * the compiler/out-of-order engine already interleaves them, so this is
 * close to the 4-way and is NOT the baseline the 4-way improves on;
 *   - 4x batched        : approx_exp_4x.
 * The meaningful number is 4x vs 1x-latency: the 4-way breaks the per-call
 * dependency chain (4 independent Horner muls / 4 independent OR-scans),
 * which is exactly what an isolated scalar call cannot do.  (4x ~
 * 1x-throughput is expected: the compiler already extracts ILP across
 * independent calls.)
 *
 * Build for meaningful cycle numbers (the library default is -O2; the
 * 4-way ILP win shows best at -O3 -march=native, which also emits
 * mulx/shrx):
 *   gcc -O3 -march=native -std=c99 -I. -DLITHIUM_MODE=120 \
 *       test/speed_approx_exp.c approx_exp.c test/cpucycles.c \
 *       -lm -o out/speed_approx_exp_120 (rdtsc measures
 * reference cycles; for retired core cycles build -DUSE_RDPMC after `echo
 * 2 > /sys/devices/cpu/rdpmc`.)
 */
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#include "approx_exp.h"
#include "cpucycles.h"

#define NCORR 4000000 /* bit-identity samples                    */
#define NTESTS 2048   /* timed samples                          */
#define KBATCH 256    /* work units per timed sample (amortize rdtsc) */
#define POOL 512      /* L1-resident input pool (512*4*8 = 16 KiB) */

static inline uint64_t mulhi(uint64_t x, uint64_t y)
{
    return mulh64(x, y);
}

static uint64_t falcon_approx_exp(uint64_t x)
{
    static const uint64_t EXPM_COEFFS[] = {
        0x00000004741183A3, 0x00000036548CFC06, 0x0000024FDCBF140A,
        0x0000171D939DE045, 0x0000D00CF58F6F84, 0x000680681CF796E3,
        0x002D82D8305B0FEA, 0x011111110E066FD0, 0x0555555555070F00,
        0x155555555581FF00, 0x400000000002B400, 0x7FFFFFFFFFFF4800,
        0x8000000000000000};
    uint64_t y = EXPM_COEFFS[0];
    for (size_t i = 1; i < (sizeof EXPM_COEFFS) / sizeof(uint64_t); i++) {
        uint64_t c = mulhi(x, y);
        y = EXPM_COEFFS[i] - c;
    }
    return y;
}

/* splitmix64 - deterministic, reproducible PRNG for the gate/pool. */
static uint64_t s_state;
static uint64_t sm64(void)
{
    uint64_t z = (s_state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

/* random a_q60 in [0, 15*2^60): a in [0,15), always < 2^64
 * (Q60-representable). */
static uint64_t rand_a_q60(void)
{
    return sm64() % (UINT64_C(15) << 60);
}

int main(void)
{
    static uint64_t t[NTESTS];
    volatile uint64_t sink = 0;
    int ok = 1;

    /* ---------------- correctness gate (1): bit-identity ----------------
     */
    s_state = 0xC0FFEEULL;
    long bad = 0;
    for (uint64_t i = 0; i < NCORR; i++) {
        uint64_t a4[4], o4[4];
        for (int l = 0; l < 4; l++)
            a4[l] = rand_a_q60();
        approx_exp_4x(o4, a4);
        for (int l = 0; l < 4; l++) {
            uint64_t ref = approx_exp(a4[l]);
            if (o4[l] != ref) {
                if (bad < 5)
                    printf("BITDIFF a=%llu 4x=%llu scalar=%llu\n",
                           (unsigned long long)a4[l],
                           (unsigned long long)o4[l],
                           (unsigned long long)ref);
                bad++;
            }
        }
    }
    if (bad)
        ok = 0;

    /* ---------------- correctness gate (2): accuracy vs libm
     * ---------------- */
    const int niter = 200000;
    double rel8 = 0.0;
    for (int i = 0; i <= niter; i++) {
        double a = 8.0 * i / niter;
        uint64_t q = (uint64_t)(a * (double)(UINT64_C(1) << 60) + 0.5);
        double got = (double)approx_exp(q) / (double)(UINT64_C(1) << 63);
        double want = exp(-a);
        double dr = fabs(got - want) / want;
        if (dr > rel8)
            rel8 = dr;
    }
    int acc_ok = (rel8 < 1.2e-13); /* 2^-43 ~ 1.137e-13 */
    if (!acc_ok)
        ok = 0;

    /* ---------------- timing pool ---------------- */
    static uint64_t pool[POOL][4];
    s_state = 0x1234567ULL;
    for (int i = 0; i < POOL; i++)
        for (int l = 0; l < 4; l++)
            pool[i][l] = rand_a_q60();

    /* warm up */
    for (int i = 0; i < POOL; i++) {
        uint64_t o[4];
        approx_exp_4x(o, pool[i]);
        sink ^= o[0];
    }

    /* 1x latency-bound: serial chain.  approx_exp output is in [0, 2^63]
     * which is a valid a_q60 (< 15*2^60), so feeding it straight back
     * makes each call depend on the previous -> exposes the full per-call
     * latency (Horner + OR-scan). */
    uint64_t acc = pool[0][0] | 1u;
    for (int s = 0; s < NTESTS; s++) {
        uint64_t a = cpucycles();
        for (int k = 0; k < KBATCH; k++)
            acc = approx_exp(acc);
        t[s] = cpucycles() - a;
    }
    sink ^= acc;
    uint64_t med1xL = cpucycles_median(t, NTESTS);

    /* 1x throughput: four INDEPENDENT scalar calls per unit (compiler
     * interleaves). */
    for (int s = 0; s < NTESTS; s++) {
        uint64_t a = cpucycles();
        for (int k = 0; k < KBATCH; k++) {
            const uint64_t *a4 = pool[(s * KBATCH + k) & (POOL - 1)];
            sink ^= approx_exp(a4[0]);
            sink ^= approx_exp(a4[1]);
            sink ^= approx_exp(a4[2]);
            sink ^= approx_exp(a4[3]);
        }
        t[s] = cpucycles() - a;
    }
    uint64_t med1xT = cpucycles_median(t, NTESTS);

    /* 4x batched: one approx_exp_4x per unit (independent across units).
     */
    for (int s = 0; s < NTESTS; s++) {
        uint64_t a = cpucycles();
        for (int k = 0; k < KBATCH; k++) {
            uint64_t o[4];
            approx_exp_4x(o, pool[(s * KBATCH + k) & (POOL - 1)]);
            sink ^= o[0] ^ o[1] ^ o[2] ^ o[3];
        }
        t[s] = cpucycles() - a;
    }
    uint64_t med4x = cpucycles_median(t, NTESTS);

    /* Falcon polynomial: same input pool, scalar 12th-degree Horner. */
    acc = pool[0][0] | 1u;
    for (int s = 0; s < NTESTS; s++) {
        uint64_t a = cpucycles();
        for (int k = 0; k < KBATCH; k++)
            acc = falcon_approx_exp(acc);
        t[s] = cpucycles() - a;
    }
    sink ^= acc;
    uint64_t medFalconL = cpucycles_median(t, NTESTS);

    for (int s = 0; s < NTESTS; s++) {
        uint64_t a = cpucycles();
        for (int k = 0; k < KBATCH; k++) {
            const uint64_t *a4 = pool[(s * KBATCH + k) & (POOL - 1)];
            sink ^= falcon_approx_exp(a4[0]);
            sink ^= falcon_approx_exp(a4[1]);
            sink ^= falcon_approx_exp(a4[2]);
            sink ^= falcon_approx_exp(a4[3]);
        }
        t[s] = cpucycles() - a;
    }
    uint64_t medFalconT = cpucycles_median(t, NTESTS);

    double cL =
        (double)med1xL / KBATCH; /* cyc per output, latency-bound */
    double cT =
        (double)med1xT / KBATCH / 4.0; /* cyc per output, throughput     */
    double c4 = (double)med4x / KBATCH / 4.0; /* cyc per output, 4x */
    double cFalconL = (double)medFalconL / KBATCH;
    double cFalconT = (double)medFalconT / KBATCH / 4.0;
    printf("1x latency-bound (serial chain) : %.2f cyc/out\n", cL);
    printf("1x throughput (4 indep calls)   : %.2f cyc/out\n", cT);
    printf("4x batched (approx_exp_4x)      : %.2f cyc/out\n", c4);
    printf("Falcon approx_exp serial chain  : %.2f cyc/out\n", cFalconL);
    // printf("Falcon approx_exp throughput    : %.2f cyc/out\n", cFalconT);
    printf(
        "speedup 4x vs 1x-latency  x%.2f   (4x vs 1x-throughput x%.2f)\n",
        cL / c4, cT / c4);
    printf(
        "%s speed_approx_exp mode=%d bitexact_bad=%ld rel[0,8]=%.3e (gate "
        "2^-43) sink=%llu\n",
        ok ? "PASS" : "FAIL", LITHIUM_MODE, bad, rel8,
        (unsigned long long)sink);
    return ok ? 0 : 1;
}
