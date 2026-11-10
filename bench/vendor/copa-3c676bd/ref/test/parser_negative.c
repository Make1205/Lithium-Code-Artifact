/* parser_negative.c - deterministic malformed-signature corpus. */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "fips202.h"
#include "params.h"
#include "sign.h"

static keccak_state drbg;
void randombytes(uint8_t *out, size_t outlen)
{
    shake256_squeeze(out, outlen, &drbg);
}
static void drbg_init(void)
{
    uint8_t seed[40];
    for (size_t i = 0; i < sizeof seed; i++)
        seed[i] = (uint8_t)(0x51u + 13u * i + LITHIUM_MODE);
    shake256_init(&drbg);
    shake256_absorb(&drbg, seed, sizeof seed);
    shake256_finalize(&drbg);
}

static int expect_reject(const char *name, const uint8_t *sig,
                         size_t siglen, const uint8_t *msg, size_t mlen,
                         const uint8_t *pk)
{
    int rc = crypto_sign_verify(sig, siglen, msg, mlen, pk);
    printf("  %-28s %s rc=%d\n", name, rc != 0 ? "reject" : "ACCEPT", rc);
    return rc != 0;
}

int main(void)
{
    drbg_init();
    uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
    uint8_t sig[CRYPTO_BYTES + 8], work[CRYPTO_BYTES + 8];
    uint8_t msg[37];
    size_t siglen = 0;
    for (size_t i = 0; i < sizeof msg; i++)
        msg[i] = (uint8_t)(0xA0u ^ (unsigned)i ^ LITHIUM_MODE);

    if (crypto_sign_keypair(pk, sk) != 0 ||
        crypto_sign_signature(sig, &siglen, msg, sizeof msg, sk) != 0 ||
        siglen != CRYPTO_BYTES ||
        crypto_sign_verify(sig, siglen, msg, sizeof msg, pk) != 0) {
        printf("FAIL parser_negative mode=%d setup\n", LITHIUM_MODE);
        return 1;
    }

    int ok = 1;
    const size_t lenoff = CHALLENGE_PACKEDBYTES;
    const size_t ransoff = lenoff + 2;
    size_t rlen = (size_t)sig[lenoff] | ((size_t)sig[lenoff + 1] << 8);
    const size_t zoff = ransoff + RANS_RESERVED_BYTES;

    memcpy(work, sig, sizeof work);
    ok &= expect_reject("signature length minus one", work,
                        CRYPTO_BYTES - 1, msg, sizeof msg, pk);

    memcpy(work, sig, sizeof work);
    work[CRYPTO_BYTES] = 0;
    ok &= expect_reject("signature length plus one", work,
                        CRYPTO_BYTES + 1, msg, sizeof msg, pk);

    memcpy(work, sig, sizeof work);
    work[lenoff] = (uint8_t)((RANS_RESERVED_BYTES + 1u) & 0xffu);
    work[lenoff + 1] = (uint8_t)((RANS_RESERVED_BYTES + 1u) >> 8);
    ok &= expect_reject("wrong reserve overflow", work, CRYPTO_BYTES, msg,
                        sizeof msg, pk);

    memcpy(work, sig, sizeof work);
    if (rlen < RANS_RESERVED_BYTES)
        work[ransoff + rlen] = 0x5a;
    else
        work[lenoff] ^= 1u;
    ok &= expect_reject("nonzero rANS padding", work, CRYPTO_BYTES, msg,
                        sizeof msg, pk);

    memcpy(work, sig, sizeof work);
    if (rlen > 0)
        work[ransoff] ^= 0x80u;
    ok &= expect_reject("noncanonical rANS stream", work, CRYPTO_BYTES,
                        msg, sizeof msg, pk);

    memcpy(work, sig, sizeof work);
    if (rlen > 0)
        work[ransoff + rlen - 1] ^= 0x01u;
    ok &= expect_reject("out-of-range hints", work, CRYPTO_BYTES, msg,
                        sizeof msg, pk);

    memcpy(work, sig, sizeof work);
    memset(work + zoff, 0xff, CRYPTO_BYTES - zoff);
    ok &= expect_reject("out-of-bound z", work, CRYPTO_BYTES, msg,
                        sizeof msg, pk);

    memcpy(work, sig, sizeof work);
    work[0] ^= 0x01u;
    work[CHALLENGE_PACKEDBYTES - 1] ^= 0x80u;
    ok &= expect_reject("malformed challenge packing", work, CRYPTO_BYTES,
                        msg, sizeof msg, pk);

    printf("%s parser_negative mode=%d cases=8\n", ok ? "PASS" : "FAIL",
           LITHIUM_MODE);
    return ok ? 0 : 1;
}
