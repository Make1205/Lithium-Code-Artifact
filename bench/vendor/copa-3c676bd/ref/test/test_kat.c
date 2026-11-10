/* test_kat.c - deterministic Known-Answer Test for IRS.
 *
 * Provides its OWN deterministic randombytes() (SHAKE256-seeded DRBG) so
 * that keygen/sign are fully reproducible.  Generates N
 * keypairs+signatures from a fixed master seed, verifies each, and prints
 * an FNV-1a hash over the whole (pk || sk || sig) transcript.  The hash is
 * the KAT: any change in algorithm, constants, or byte layout changes it.
 * Re-running must reproduce it exactly, and verification must accept all
 * and reject tampered.
 *
 * Build WITHOUT ref/randombytes.c (this file supplies randombytes).
 */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "fips202.h"
#include "params.h"
#include "sign.h"

/* ---- deterministic DRBG: squeeze a single long SHAKE256 stream ---- */
static keccak_state drbg;
void randombytes(uint8_t *out, size_t outlen)
{
    shake256_squeeze(out, outlen, &drbg);
}
static void drbg_init(const uint8_t *seed, size_t n)
{
    shake256_init(&drbg);
    shake256_absorb(&drbg, seed, n);
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

int main(void)
{
    uint8_t master[32];
    for (int i = 0; i < 32; i++)
        master[i] = (uint8_t)(0xA5 ^ (i * 31 + LITHIUM_MODE));
    drbg_init(master, 32);

    const int nkat = 10;
    uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES],
        sig[CRYPTO_BYTES];
    size_t siglen;
    uint64_t h = 1469598103934665603ULL;
    int accept = 0, reject = 0;

    for (int t = 0; t < nkat; t++) {
        uint8_t msg[33];
        for (int i = 0; i < 33; i++)
            msg[i] = (uint8_t)(t * 7 + i);

        if (crypto_sign_keypair(pk, sk) != 0) {
            printf("FAIL test_kat keypair\n");
            return 1;
        }
        if (crypto_sign_signature(sig, &siglen, msg, sizeof msg, sk) !=
            0) {
            printf("FAIL test_kat sign\n");
            return 1;
        }

        h = fnv(h, pk, CRYPTO_PUBLICKEYBYTES);
        h = fnv(h, sk, CRYPTO_SECRETKEYBYTES);
        h = fnv(h, sig, siglen);

        if (crypto_sign_verify(sig, siglen, msg, sizeof msg, pk) == 0)
            accept++;
        sig[siglen / 2] ^= 0x40;
        if (crypto_sign_verify(sig, siglen, msg, sizeof msg, pk) != 0)
            reject++;
    }

    int ok = (accept == nkat) && (reject == nkat);
    printf(
        "%s test_kat mode=%d N=%d accept=%d reject=%d KAT_hash=%llu "
        "(pk=%d sk=%d sig=%zu)\n",
        ok ? "PASS" : "FAIL", LITHIUM_MODE, nkat, accept, reject,
        (unsigned long long)h, CRYPTO_PUBLICKEYBYTES,
        CRYPTO_SECRETKEYBYTES, siglen);
    return ok ? 0 : 1;
}
