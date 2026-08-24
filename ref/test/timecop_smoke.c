/* timecop_smoke.c - ref mode-120 dynamic variable-latency smoke harness.
 */
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <valgrind/memcheck.h>
#include <valgrind/valgrind.h>

#include "fips202.h"
#include "params.h"
#include "sign.h"

static keccak_state drbg;
static int signing_phase;

void randombytes(uint8_t *out, size_t outlen)
{
    shake256_squeeze(out, outlen, &drbg);
    if (signing_phase && RUNNING_ON_VALGRIND)
        VALGRIND_MAKE_MEM_UNDEFINED(out, outlen);
}

static void drbg_init(void)
{
    uint8_t seed[32];
    for (size_t i = 0; i < sizeof seed; i++)
        seed[i] = (uint8_t)(0x5a ^ (uint8_t)(17 * i));
    shake256_init(&drbg);
    shake256_absorb(&drbg, seed, sizeof seed);
    shake256_finalize(&drbg);
}

static void mark_secret_key_material(uint8_t *sk)
{
    if (!RUNNING_ON_VALGRIND)
        return;
    size_t key_off = SEEDBYTES + TRBYTES;
    VALGRIND_MAKE_MEM_UNDEFINED(sk + key_off, SKSEEDBYTES);
    size_t noise_off = key_off + SKSEEDBYTES;
    VALGRIND_MAKE_MEM_UNDEFINED(sk + noise_off,
                                (size_t)K * POLYNOISE_PACKEDBYTES);
    size_t bn_off = noise_off + (size_t)K * POLYNOISE_PACKEDBYTES;
    VALGRIND_MAKE_MEM_DEFINED(sk, key_off);
    VALGRIND_MAKE_MEM_DEFINED(sk + bn_off, (size_t)M * POLYPK_PACKEDBYTES);
}

int main(void)
{
    drbg_init();
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t sig[CRYPTO_BYTES];
    uint8_t msg[33];
    size_t siglen = 0;

    for (size_t i = 0; i < sizeof msg; i++)
        msg[i] = (uint8_t)(0xa0u + i);

    signing_phase = 0;
    if (crypto_sign_keypair(pk, sk) != 0) {
        puts("FAIL timecop_smoke keypair");
        return 1;
    }

    mark_secret_key_material(sk);
    signing_phase = 1;
    if (crypto_sign_signature(sig, &siglen, msg, sizeof msg, sk) != 0) {
        puts("FAIL timecop_smoke sign");
        return 1;
    }
    signing_phase = 0;

    if (RUNNING_ON_VALGRIND) {
        VALGRIND_MAKE_MEM_DEFINED(sig, siglen);
        VALGRIND_MAKE_MEM_DEFINED(&siglen, sizeof siglen);
        VALGRIND_MAKE_MEM_DEFINED(pk, sizeof pk);
        VALGRIND_MAKE_MEM_DEFINED(msg, sizeof msg);
    }

    if (crypto_sign_verify(sig, siglen, msg, sizeof msg, pk) != 0) {
        puts("FAIL timecop_smoke verify");
        return 1;
    }

    printf("PASS timecop_smoke mode=%d siglen=%zu\n", LITHIUM_MODE,
           siglen);
    return 0;
}
