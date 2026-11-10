/* dudect_sign.c - end-to-end fixed-secret vs random-secret timing test.
 *
 * This is a small dudect-style Welch t-test harness for
 * crypto_sign_signature. It deliberately provides its own randombytes():
 * key generation during setup uses a deterministic PRNG, while measured
 * signing calls consume the per-trace signing randomness prepared by this
 * harness.  The timed region contains only crypto_sign_signature(), not
 * key selection or input copying.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "cpucycles.h"
#include "params.h"
#include "sign.h"

#ifndef DUDECT_N
#    define DUDECT_N 2000
#endif

#ifndef DUDECT_SK_POOL
#    define DUDECT_SK_POOL 64
#endif

#define DUDECT_MSG_LEN 32
#define DUDECT_INPUT_LEN (DUDECT_MSG_LEN + RNDBYTES + 2)
#define DUDECT_SK_OFF_RHO 0
#define DUDECT_SK_OFF_TR (DUDECT_SK_OFF_RHO + SEEDBYTES)
#define DUDECT_SK_OFF_KEY (DUDECT_SK_OFF_TR + TRBYTES)
#define DUDECT_SK_OFF_NOISE (DUDECT_SK_OFF_KEY + SKSEEDBYTES)
#define DUDECT_SK_NOISE_BYTES (K * POLYNOISE_PACKEDBYTES)
#define DUDECT_SK_OFF_BN (DUDECT_SK_OFF_NOISE + DUDECT_SK_NOISE_BYTES)
#define DUDECT_SK_BN_BYTES (M * POLYPK_PACKEDBYTES)

static uint64_t prng_state = 0x9e3779b97f4a7c15ULL;
static uint8_t fixed_sk[CRYPTO_SECRETKEYBYTES];
static uint8_t random_sks[DUDECT_SK_POOL][CRYPTO_SECRETKEYBYTES];
static uint8_t active_sk[CRYPTO_SECRETKEYBYTES];
static uint8_t current_rnd[RNDBYTES];
static int in_measured_sign;

static uint64_t xorshift64star(void)
{
    uint64_t x = prng_state;

    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    prng_state = x;
    return x * 2685821657736338717ULL;
}

static void prng_bytes(uint8_t *out, size_t outlen)
{
    uint64_t x = 0;
    size_t left = 0;

    for (size_t i = 0; i < outlen; i++) {
        if (left == 0) {
            x = xorshift64star();
            left = 8;
        }
        out[i] = (uint8_t)x;
        x >>= 8;
        left--;
    }
}

void randombytes(uint8_t *out, size_t outlen)
{
    if (in_measured_sign) {
        if (outlen != RNDBYTES) {
            fprintf(stderr,
                    "unexpected randombytes() during measured sign: %zu\n",
                    outlen);
        }
        memcpy(out, current_rnd, outlen <= RNDBYTES ? outlen : RNDBYTES);
        if (outlen > RNDBYTES)
            memset(out + RNDBYTES, 0, outlen - RNDBYTES);
        return;
    }

    prng_bytes(out, outlen);
}

typedef struct dudect_stats {
    uint64_t n;
    double mean;
    double m2;
} dudect_stats;

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

static void copy_public_sk_fields(uint8_t *dst, const uint8_t *src)
{
    memcpy(dst + DUDECT_SK_OFF_RHO, src + DUDECT_SK_OFF_RHO, SEEDBYTES);
    memcpy(dst + DUDECT_SK_OFF_TR, src + DUDECT_SK_OFF_TR, TRBYTES);
    memcpy(dst + DUDECT_SK_OFF_BN, src + DUDECT_SK_OFF_BN,
           DUDECT_SK_BN_BYTES);
}

static void setup_keys(void)
{
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];

    if (crypto_sign_keypair(pk, fixed_sk) != 0) {
        fprintf(stderr, "fixed keypair failed\n");
        return;
    }

    for (size_t i = 0; i < DUDECT_SK_POOL; i++) {
        if (crypto_sign_keypair(pk, random_sks[i]) != 0) {
            fprintf(stderr, "random keypair %zu failed\n", i);
            return;
        }
        /* Keep the public transcript fixed across classes.  Class 1 keeps
         * its independently generated key/noise secret fields, while rho,
         * tr, and bn match class 0 so public A0 expansion does not
         * dominate the fixed-vs-random-secret timing test. */
        copy_public_sk_fields(random_sks[i], fixed_sk);
    }
}

static uint64_t measure_one(const uint8_t *msg, const uint8_t *rnd,
                            const uint8_t *sk)
{
    uint8_t sig[CRYPTO_BYTES];
    size_t siglen = 0;
    uint64_t t0, t1;
    int ret;

    memcpy(active_sk, sk, sizeof active_sk);
    memcpy(current_rnd, rnd, sizeof current_rnd);

    in_measured_sign = 1;
    t0 = cpucycles();
    ret = crypto_sign_signature(sig, &siglen, msg, DUDECT_MSG_LEN,
                                active_sk);
    t1 = cpucycles();
    in_measured_sign = 0;

    if (ret != 0 || siglen != CRYPTO_BYTES)
        fprintf(stderr, "sign failed: ret=%d siglen=%zu\n", ret, siglen);

    return t1 - t0;
}

int main(void)
{
    dudect_stats st[2] = {{0, 0.0, 0.0}, {0, 0.0, 0.0}};
    uint8_t input[DUDECT_INPUT_LEN];
    uint64_t overhead = cpucycles_overhead();

    setup_keys();

    for (size_t i = 0; i < 100; i++) {
        prng_bytes(input, sizeof input);
        (void)measure_one(input, input + DUDECT_MSG_LEN, fixed_sk);
    }

    for (size_t i = 0; i < DUDECT_N; i++) {
        const uint8_t *sk;
        uint64_t cycles;
        unsigned cls;

        prng_bytes(input, sizeof input);
        cls = input[DUDECT_INPUT_LEN - 1] & 1U;
        sk =
            cls == 0
                ? fixed_sk
                : random_sks[input[DUDECT_INPUT_LEN - 2] % DUDECT_SK_POOL];

        cycles = measure_one(input, input + DUDECT_MSG_LEN, sk);
        if (cycles > overhead)
            cycles -= overhead;
        stats_push(&st[cls], (double)cycles);
    }

    printf(
        "dudect-sign mode=%d traces=%d class0=fixed-sk "
        "class1=random-secret-same-public sk_pool=%d\n",
        LITHIUM_MODE, DUDECT_N, DUDECT_SK_POOL);
    printf("  n0=%llu mean0=%.1f var0=%.1f\n", (unsigned long long)st[0].n,
           st[0].mean, stats_var(&st[0]));
    printf("  n1=%llu mean1=%.1f var1=%.1f\n", (unsigned long long)st[1].n,
           st[1].mean, stats_var(&st[1]));
    printf("  Welch t=%.2f |t|=%.2f\n", welch_t(&st[0], &st[1]),
           fabs(welch_t(&st[0], &st[1])));

    return 0;
}
