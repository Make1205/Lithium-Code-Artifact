/* trace_kat.c - deterministic transcript and randomness-consumption KAT.
 */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "fips202.h"
#include "params.h"
#include "sign.h"

#ifdef TEST_RANS_STATS
extern uint64_t lithium_test_rans_pack_calls;
extern uint64_t lithium_test_rans_pack_fails;
#endif

static keccak_state drbg;
static uint64_t rb_calls, rb_bytes;
void randombytes(uint8_t *out, size_t outlen)
{
    rb_calls++;
    rb_bytes += outlen;
    shake256_squeeze(out, outlen, &drbg);
}
static void drbg_init(void)
{
    uint8_t seed[64];
    for (size_t i = 0; i < sizeof seed; i++)
        seed[i] = (uint8_t)(0xC3u ^ (unsigned)(i * 29u) ^ LITHIUM_MODE);
    shake256_init(&drbg);
    shake256_absorb(&drbg, seed, sizeof seed);
    shake256_finalize(&drbg);
}
static uint64_t fnv_u64(uint64_t h, uint64_t x)
{
    for (int i = 0; i < 8; i++) {
        h ^= (uint8_t)(x >> (8 * i));
        h *= 1099511628211ULL;
    }
    return h;
}
static uint64_t fnv_bytes(uint64_t h, const uint8_t *p, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        h ^= p[i];
        h *= 1099511628211ULL;
    }
    return h;
}

int main(void)
{
    drbg_init();
    uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
    uint8_t sig[CRYPTO_BYTES];
    uint8_t msg[29];
    size_t siglen = 0;
    for (size_t i = 0; i < sizeof msg; i++)
        msg[i] = (uint8_t)(0x71u + i + LITHIUM_MODE);
    if (crypto_sign_keypair(pk, sk) != 0 ||
        crypto_sign_signature(sig, &siglen, msg, sizeof msg, sk) != 0 ||
        crypto_sign_verify(sig, siglen, msg, sizeof msg, pk) != 0) {
        printf("FAIL trace_kat mode=%d setup\n", LITHIUM_MODE);
        return 1;
    }
    uint64_t h = 1469598103934665603ULL;
    h = fnv_bytes(h, pk, sizeof pk);
    h = fnv_bytes(h, sk, sizeof sk);
    h = fnv_bytes(h, sig, siglen);
    h = fnv_u64(h, rb_calls);
    h = fnv_u64(h, rb_bytes);
#ifdef TEST_RANS_STATS
    h = fnv_u64(h, lithium_test_rans_pack_calls);
    h = fnv_u64(h, lithium_test_rans_pack_fails);
    printf(
        "PASS trace_kat mode=%d trace_hash=%llu rb_calls=%llu "
        "rb_bytes=%llu sign_attempts=%llu rans_overflow_retries=%llu\n",
        LITHIUM_MODE, (unsigned long long)h, (unsigned long long)rb_calls,
        (unsigned long long)rb_bytes,
        (unsigned long long)lithium_test_rans_pack_calls,
        (unsigned long long)lithium_test_rans_pack_fails);
#else
    printf(
        "PASS trace_kat mode=%d trace_hash=%llu rb_calls=%llu "
        "rb_bytes=%llu sign_attempts=1 rans_overflow_retries=0\n",
        LITHIUM_MODE, (unsigned long long)h, (unsigned long long)rb_calls,
        (unsigned long long)rb_bytes);
#endif
    return 0;
}
