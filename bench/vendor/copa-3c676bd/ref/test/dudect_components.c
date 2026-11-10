/* dudect_components.c - optional component-level fixed-vs-random timing
 * smoke. */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "cpucycles.h"
#include "polyvec.h"
#include "rounding.h"

#ifndef DUDECT_COMPONENT_N
#    define DUDECT_COMPONENT_N 512
#endif

typedef struct dudect_stats {
    uint64_t n;
    double mean;
    double m2;
} dudect_stats;

typedef void (*component_fn)(unsigned cls, size_t idx);

static poly fixed_poly, random_poly[DUDECT_COMPONENT_N];
static poly16 ntt_buf;
static volatile int64_t sink64;
static volatile int32_t sink32;
static uint64_t prng_state = 0x6a09e667f3bcc909ULL;

static uint64_t xorshift64star(void)
{
    uint64_t x = prng_state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    prng_state = x;
    return x * 2685821657736338717ULL;
}

static void stats_push(dudect_stats *s, double x)
{
    double delta;
    s->n++;
    delta = x - s->mean;
    s->mean += delta / (double)s->n;
    s->m2 += delta * (x - s->mean);
}

static double stats_var(const dudect_stats *s)
{
    return s->n > 1 ? s->m2 / (double)(s->n - 1) : 0.0;
}

static double welch_t(const dudect_stats *a, const dudect_stats *b)
{
    double va = stats_var(a);
    double vb = stats_var(b);
    double denom = sqrt(va / (double)a->n + vb / (double)b->n);
    return denom == 0.0 ? 0.0 : (a->mean - b->mean) / denom;
}

static const poly *select_poly(unsigned cls, size_t idx)
{
    return cls == 0 ? &fixed_poly : &random_poly[idx % DUDECT_COMPONENT_N];
}

static void setup_inputs(void)
{
    memset(&fixed_poly, 0, sizeof fixed_poly);
    for (size_t i = 0; i < DUDECT_COMPONENT_N; i++) {
        for (int k = 0; k < N; k++) {
            uint64_t r = xorshift64star();
            int32_t v = (int32_t)(r % (uint64_t)(2 * Q - 1)) - (Q - 1);
            random_poly[i].coeffs[k] = v;
        }
    }
}

static void comp_sqnorm(unsigned cls, size_t idx)
{
    sink64 ^= poly_sqnorm(select_poly(cls, idx));
}

static void comp_ntt_roundtrip(unsigned cls, size_t idx)
{
    poly out;
    poly_ntt(&ntt_buf, select_poly(cls, idx));
    poly_invntt_tomont(&out, &ntt_buf);
    sink32 ^= out.coeffs[(idx * 13u) % N];
}

static void comp_rounding_batch(unsigned cls, size_t idx)
{
    const poly *p = select_poly(cls, idx);
    int32_t acc = 0;
    for (int k = 0; k < N; k++) {
        int32_t w = reduce2q(p->coeffs[k] + (int32_t)(idx & 255u));
        int32_t hb = highbits(w);
        acc ^= hbvalue(hb);
        acc ^= lsb(w);
        acc ^= makehint(w, w ^ 1);
        acc ^= usehint(w, makehint(w, w ^ 1));
    }
    sink32 ^= acc;
}

static int run_component(const char *name, component_fn fn)
{
    dudect_stats st[2] = {{0, 0.0, 0.0}, {0, 0.0, 0.0}};
    uint64_t overhead = cpucycles_overhead();

    for (size_t i = 0; i < 32; i++)
        fn((unsigned)(i & 1u), i);

    for (size_t i = 0; i < DUDECT_COMPONENT_N; i++) {
        unsigned cls = (unsigned)(xorshift64star() & 1u);
        uint64_t t0 = cpucycles();
        fn(cls, i);
        uint64_t t1 = cpucycles();
        uint64_t cyc = t1 - t0;
        if (cyc > overhead)
            cyc -= overhead;
        stats_push(&st[cls], (double)cyc);
    }

    double t = welch_t(&st[0], &st[1]);
    double at = fabs(t);
    const char *status = at < 4.5 ? "PASS" : at < 10.0 ? "WARN" : "FAIL";
    printf(
        "COMPONENT mode=%d name=%s traces=%d n0=%llu n1=%llu mean0=%.1f "
        "mean1=%.1f t=%.2f abs_t=%.2f status=%s\n",
        LITHIUM_MODE, name, DUDECT_COMPONENT_N,
        (unsigned long long)st[0].n, (unsigned long long)st[1].n,
        st[0].mean, st[1].mean, t, at, status);
    return status[0] == 'F';
}

int main(void)
{
    int fail = 0;
    setup_inputs();
    printf(
        "DUDECT_COMPONENTS mode=%d traces=%d "
        "scope=poly_sqnorm,ntt_roundtrip,rounding_batch\n",
        LITHIUM_MODE, DUDECT_COMPONENT_N);
    fail |= run_component("poly_sqnorm", comp_sqnorm);
    fail |= run_component("poly_ntt_invntt_roundtrip", comp_ntt_roundtrip);
    fail |= run_component("rounding_highbits_lowbits_lsb_batch",
                          comp_rounding_batch);
    printf("%s dudect_components mode=%d sink=%lld/%d\n",
           fail ? "WARN" : "PASS", LITHIUM_MODE, (long long)sink64,
           (int)sink32);
    return 0;
}
