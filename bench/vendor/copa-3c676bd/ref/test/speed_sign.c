/*
 * speed_sign.c - end-to-end cycle benchmark for keygen / sign / verify.
 *
 * keygen and sign use the tail benchmark because rejection sampling makes
 * their latency distributions wider. verify remains a median benchmark.
 */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "cpucycles.h"
#include "params.h"
#include "randombytes.h"
#include "sign.h"

#define NV 10000
#define MLEN 59

typedef struct keygen_bench_ctx {
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
} keygen_bench_ctx;

typedef struct sign_bench_ctx {
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t sig[CRYPTO_BYTES];
    size_t siglen;
} sign_bench_ctx;

static int keygen_once(void *vctx)
{
    keygen_bench_ctx *ctx = vctx;

    return crypto_sign_keypair(ctx->pk, ctx->sk);
}

static int sign_once(void *vctx)
{
    sign_bench_ctx *ctx = vctx;

    return crypto_sign_signature(ctx->sig, &ctx->siglen, ctx->sig, MLEN,
                                 ctx->sk);
}

int main(void)
{
    uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES],
        sig[CRYPTO_BYTES];
    uint8_t msg[MLEN];
    size_t siglen = 0;
    static uint64_t cv[NV];
    keygen_bench_ctx keygen_ctx;
    sign_bench_ctx sign_ctx;
    tail_sign_result keygen_result;
    tail_sign_result sign_result;

    randombytes(msg, sizeof msg);

    if (bench_tail_sign("keygen", &keygen_result, keygen_once,
                        &keygen_ctx) != 0)
        return 1;

    memcpy(pk, keygen_ctx.pk, sizeof pk);
    memcpy(sk, keygen_ctx.sk, sizeof sk);

    memset(&sign_ctx, 0, sizeof sign_ctx);
    memcpy(sign_ctx.sk, sk, sizeof sign_ctx.sk);

    if (bench_tail_sign("sign", &sign_result, sign_once, &sign_ctx) != 0)
        return 1;

    crypto_sign_signature(sig, &siglen, msg, sizeof msg, sk);
    for (int i = 0; i < NV; i++) {
        uint64_t a = cpucycles();
        (void)crypto_sign_verify(sig, siglen, msg, sizeof msg, pk);
        cv[i] = cpucycles() - a;
    }

    printf("mode=%d, pk=%dB, sig=%dB, pk+sig=%dB, mlen=%dB\n",
           LITHIUM_MODE, CRYPTO_PUBLICKEYBYTES, CRYPTO_BYTES,
           CRYPTO_PUBLICKEYBYTES + CRYPTO_BYTES, MLEN);
    bench_print_tail(&keygen_result);
    bench_print_tail(&sign_result);
    printf("verify median (%d runs): %llu cycles\n", NV,
           (unsigned long long)cpucycles_median(cv, NV));
    return 0;
}
