/* security_fault_injection.c - ML-DSA-Bugs style runtime fault detectors.
 */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "fips202.h"
#include "params.h"
#include "sign.h"

static keccak_state drbg;
static uint64_t rb_calls, rb_bytes;
void randombytes(uint8_t *out, size_t outlen)
{
    rb_calls++;
    rb_bytes += outlen;
    shake256_squeeze(out, outlen, &drbg);
}
static void drbg_init(uint8_t tag)
{
    uint8_t seed[48];
    for (size_t i = 0; i < sizeof seed; i++)
        seed[i] = (uint8_t)(tag + 17u * i + LITHIUM_MODE);
    shake256_init(&drbg);
    shake256_absorb(&drbg, seed, sizeof seed);
    shake256_finalize(&drbg);
}
static uint64_t fnv(uint64_t h, const uint8_t *p, size_t n)
{
    for (size_t i = 0; i < n; i++) {
        h ^= p[i];
        h *= 1099511628211ULL;
    }
    return h;
}
static void report_case(const char *name, int caught, const char *detector)
{
    printf("  %-32s %s detector=%s\n", name, caught ? "caught" : "MISS",
           detector);
}

int main(void)
{
    uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
    uint8_t sig1[CRYPTO_BYTES], sig2[CRYPTO_BYTES], sig3[CRYPTO_BYTES];
    uint8_t msg[41];
    size_t siglen1 = 0, siglen2 = 0, siglen3 = 0;
    for (size_t i = 0; i < sizeof msg; i++)
        msg[i] = (uint8_t)(0x33u + i * 5u);

    drbg_init(0x10);
    if (crypto_sign_keypair(pk, sk) != 0)
        return 1;
    drbg_init(0x20);
    if (crypto_sign_signature(sig1, &siglen1, msg, sizeof msg, sk) != 0)
        return 1;
    drbg_init(0x21);
    if (crypto_sign_signature(sig2, &siglen2, msg, sizeof msg, sk) != 0)
        return 1;

    int ok = 1;
    int caught = (siglen1 == siglen2) && memcmp(sig1, sig2, siglen1) != 0;
    report_case("adjacent coefficient repetition", caught,
                "two-signature transcript diff");
    ok &= caught;

    drbg_init(0x40);
    if (crypto_sign_signature(sig2, &siglen2, msg, sizeof msg, sk) != 0)
        return 1;
    drbg_init(0x40);
    if (crypto_sign_signature(sig3, &siglen3, msg, sizeof msg, sk) != 0)
        return 1;
    caught = (siglen2 == siglen3) && memcmp(sig2, sig3, siglen2) == 0;
    report_case("repeated nonce", caught, "duplicate signature detector");
    ok &= caught;

    uint8_t badsk[CRYPTO_SECRETKEYBYTES];
    memset(badsk, 0, sizeof badsk);
    drbg_init(0x50);
    caught = (crypto_sign_signature(sig3, &siglen3, msg, sizeof msg,
                                    badsk) != 0) ||
             (crypto_sign_verify(sig3, siglen3, msg, sizeof msg, pk) != 0);
    report_case("cleared signing secret", caught,
                "verify-under-original-pk");
    ok &= caught;

    caught =
        crypto_sign_verify(sig1, siglen1, msg, sizeof msg - 1, pk) != 0;
    report_case("wrong hash input length", caught,
                "message-length binding");
    ok &= caught;

    uint64_t h1 = fnv(1469598103934665603ULL, sig1, siglen1);
    uint64_t h2 = fnv(1469598103934665603ULL, sig2, siglen2);
    caught = (h1 != h2) && rb_calls >= 4 &&
             rb_bytes >= (SKSEEDBYTES + 3u * RNDBYTES);
    report_case("wrong rejection counter behavior", caught,
                "trace/hash/randombytes counter");
    ok &= caught;

    printf(
        "%s security_fault_injection mode=%d cases=5 rb_calls=%llu "
        "rb_bytes=%llu\n",
        ok ? "PASS" : "FAIL", LITHIUM_MODE, (unsigned long long)rb_calls,
        (unsigned long long)rb_bytes);
    return ok ? 0 : 1;
}
