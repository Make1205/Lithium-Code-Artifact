/* explore_exp.c - sweep harness for (B, degree, k-way) ApproxExp configs.
 * Build: gcc -O3 -march=native -DKWAY=k -I/tmp/exp -include cfg.h
 * explore_exp.c \ <ircode>/ref/test/cpucycles.c -lm
 * cfg.h defines
 * EB_BETA/EB_B/EB_DEG/EB_RCP/EB_LB/EB_POLY[]/EB_TABLE[].   */
#include <math.h>
#include <stdint.h>
#include <stdio.h>

#include "cfg.h"
#include "cpucycles.h"
#ifndef KWAY
#    define KWAY 4
#endif

static inline uint64_t mulh(uint64_t a, uint64_t b)
{
    return (uint64_t)(((unsigned __int128)a * b) >> 64);
}
static inline uint64_t combine(uint64_t W, uint64_t n)
{
    uint64_t lt = (uint64_t)(-(int64_t)(n < 64));
    unsigned sh = (unsigned)(n & 63);
    uint64_t rb = ((UINT64_C(1) << sh) >> 1) & lt;
    return (uint64_t)((((unsigned __int128)W + rb) >> sh)) & lt;
}

static uint64_t aexp(uint64_t a)
{
    uint64_t V = mulh(a, EB_RCP);
    uint64_t m = ((V + ((UINT64_C(1) << 53) - 1)) >> 53) + (EB_B - 1),
             n = m >> EB_BETA, j = m & (EB_B - 1);
    unsigned __int128 R = (unsigned __int128)m * EB_LB -
                          ((unsigned __int128)a << (EB_BETA + 4));
    uint64_t X = (uint64_t)(R >> EB_BETA), p = EB_POLY[EB_DEG];
    for (int k = EB_DEG - 1; k >= 0; k--)
        p = mulh(p, X) + EB_POLY[k];
    uint64_t s = 0;
    for (uint64_t i = 0; i < EB_B; i++) {
        uint64_t d = i ^ j, mk = ((d | (uint64_t)(-(int64_t)d)) >> 63) - 1;
        s |= EB_TABLE[i] & mk;
    }
    return combine(mulh(p, s) << 1, n);
}
static void aexp_kx(uint64_t *out, const uint64_t *a)
{
    uint64_t X[KWAY], n[KWAY], j[KWAY], p[KWAY], s[KWAY];
    for (int l = 0; l < KWAY; l++) {
        uint64_t V = mulh(a[l], EB_RCP);
        uint64_t m = ((V + ((UINT64_C(1) << 53) - 1)) >> 53) + (EB_B - 1);
        n[l] = m >> EB_BETA;
        j[l] = m & (EB_B - 1);
        unsigned __int128 R = (unsigned __int128)m * EB_LB -
                              ((unsigned __int128)a[l] << (EB_BETA + 4));
        X[l] = (uint64_t)(R >> EB_BETA);
    }
    for (int l = 0; l < KWAY; l++)
        p[l] = EB_POLY[EB_DEG];
    for (int k = EB_DEG - 1; k >= 0; k--) {
        uint64_t c = EB_POLY[k];
        for (int l = 0; l < KWAY; l++)
            p[l] = mulh(p[l], X[l]) + c;
    }
    for (int l = 0; l < KWAY; l++)
        s[l] = 0;
    for (uint64_t i = 0; i < EB_B; i++) {
        uint64_t Ti = EB_TABLE[i];
        for (int l = 0; l < KWAY; l++) {
            uint64_t d = i ^ j[l],
                     mk = ((d | (uint64_t)(-(int64_t)d)) >> 63) - 1;
            s[l] |= Ti & mk;
        }
    }
    for (int l = 0; l < KWAY; l++)
        out[l] = combine(mulh(p[l], s[l]) << 1, n[l]);
}

#define NT 2048
#define KB 256
#define POOL 4096
static uint64_t st = 0x9e37ULL;
static uint64_t sm(void)
{
    uint64_t z = (st += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static uint64_t ra(void)
{
    return sm() % (UINT64_C(15) << 60);
}

int main(void)
{
    static uint64_t t[NT];
    volatile uint64_t sink = 0;
    int ok = 1;
    long bad = 0;
    /* correctness: kx==scalar */
    for (int i = 0; i < 300000; i++) {
        uint64_t a[KWAY], o[KWAY];
        for (int l = 0; l < KWAY; l++)
            a[l] = ra();
        aexp_kx(o, a);
        for (int l = 0; l < KWAY; l++)
            if (o[l] != aexp(a[l]))
                bad++;
    }
    if (bad)
        ok = 0;
    double rel8 = 0;
    int M = 100000;
    for (int i = 0; i <= M; i++) {
        double x = 8.0 * i / M;
        uint64_t q = (uint64_t)(x * (double)(UINT64_C(1) << 60) + 0.5);
        double g = (double)aexp(q) / (double)(UINT64_C(1) << 63),
               w = exp(-x), d = fabs(g - w) / w;
        if (d > rel8)
            rel8 = d;
    }
    if (rel8 >= 1.2e-13)
        ok = 0;
    static uint64_t pool[POOL][KWAY];
    for (int i = 0; i < POOL; i++)
        for (int l = 0; l < KWAY; l++)
            pool[i][l] = ra();
    /* latency-bound scalar */
    uint64_t acc = pool[0][0] | 1u;
    for (int s = 0; s < NT; s++) {
        uint64_t a = cpucycles();
        for (int k = 0; k < KB; k++)
            acc = aexp(acc);
        t[s] = cpucycles() - a;
    }
    sink ^= acc;
    double lat = (double)cpucycles_median(t, NT) / KB;
    /* k-way throughput */
    for (int s = 0; s < NT; s++) {
        uint64_t a = cpucycles();
        for (int k = 0; k < KB; k++) {
            uint64_t o[KWAY];
            aexp_kx(o, pool[(s * KB + k) & (POOL - 1)]);
            for (int l = 0; l < KWAY; l++)
                sink ^= o[l];
        }
        t[s] = cpucycles() - a;
    }
    double thr = (double)cpucycles_median(t, NT) / KB / KWAY;
    printf(
        "B=%2d D=%d K=%d  lat=%5.2f  thru=%5.2f cyc/out  rel=2^%.1f  %s "
        "sink=%llu\n",
        EB_B, EB_DEG, KWAY, lat, thr, log2(rel8), ok ? "OK" : "BAD",
        (unsigned long long)sink);
    return ok ? 0 : 1;
}
