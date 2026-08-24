/* diff_fuzz.c - deterministic bounded differential fuzz transcript. */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "fips202.h"
#include "packing.h"
#include "params.h"
#include "sign.h"

#ifndef DIFF_FUZZ_N
#    define DIFF_FUZZ_N 32
#endif

static keccak_state drbg;

void randombytes(uint8_t *out, size_t outlen)
{
    shake256_squeeze(out, outlen, &drbg);
}

static void drbg_init(void)
{
    uint8_t seed[48];
    for (size_t i = 0; i < sizeof seed; i++)
        seed[i] = (uint8_t)(0x5au ^ (unsigned)(i * 37u) ^ LITHIUM_MODE);
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

static uint64_t fnv_poly(uint64_t h, const poly *p)
{
    for (int i = 0; i < N; i++)
        h = fnv_u64(h, (uint32_t)p->coeffs[i]);
    return h;
}

static unsigned challenge_weight_from_sig(const uint8_t *sig)
{
    unsigned w = 0;
    for (int i = 0; i < N; i++)
        w += (unsigned)((sig[i >> 3] >> (i & 7)) & 1U);
    return w;
}

static unsigned padding_is_zero(const uint8_t *sig, size_t rlen)
{
    const uint8_t *rp = sig + CHALLENGE_PACKEDBYTES + 2;
    if (rlen > RANS_RESERVED_BYTES)
        return 0;
    for (size_t i = rlen; i < (size_t)RANS_RESERVED_BYTES; i++)
        if (rp[i] != 0)
            return 0;
    return 1;
}

static uint64_t unpack_summary(uint64_t h, const char *label,
                               const uint8_t *sig)
{
    poly c, hint[M], z1[ELL + 1];
    size_t off = CHALLENGE_PACKEDBYTES;
    size_t rlen = (size_t)sig[off] | ((size_t)sig[off + 1] << 8);
    unsigned pad = padding_is_zero(sig, rlen);
    unsigned cw = challenge_weight_from_sig(sig);
    int rc = unpack_sig(&c, hint, z1, sig);
    uint64_t ch = 1469598103934665603ULL;
    uint64_t hh = 1469598103934665603ULL;
    uint64_t zh = 1469598103934665603ULL;

    if (rc == 0) {
        ch = fnv_poly(ch, &c);
        for (int i = 0; i < M; i++)
            hh = fnv_poly(hh, &hint[i]);
        for (int i = 0; i < ELL + 1; i++)
            zh = fnv_poly(zh, &z1[i]);
    }
    h = fnv_u64(h, rlen);
    h = fnv_u64(h, pad);
    h = fnv_u64(h, cw);
    h = fnv_u64(h, (uint32_t)rc);
    h = fnv_u64(h, ch);
    h = fnv_u64(h, hh);
    h = fnv_u64(h, zh);
    printf(
        " %s_rans_len=%zu %s_pad_zero=%u %s_chal_w=%u %s_unpack_rc=%d "
        "%s_c_hash=%llu %s_h_hash=%llu %s_z1_hash=%llu",
        label, rlen, label, pad, label, cw, label, rc, label,
        (unsigned long long)ch, label, (unsigned long long)hh, label,
        (unsigned long long)zh);
    return h;
}

int main(void)
{
    uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
    uint8_t sig[CRYPTO_BYTES], malformed[CRYPTO_BYTES];
    uint8_t msg[96], badmsg[96];
    uint64_t transcript = 1469598103934665603ULL;
    int ok = 1;

    drbg_init();
    printf("DIFF_FUZZ mode=%d cases=%d pk=%d sk=%d sig=%d\n", LITHIUM_MODE,
           DIFF_FUZZ_N, CRYPTO_PUBLICKEYBYTES, CRYPTO_SECRETKEYBYTES,
           CRYPTO_BYTES);

    for (int t = 0; t < DIFF_FUZZ_N; t++) {
        size_t mlen = (size_t)((t * 17 + LITHIUM_MODE) % sizeof msg);
        size_t siglen = 0;
        unsigned bitmap = 0;
        uint64_t case_hash = 1469598103934665603ULL;

        for (size_t i = 0; i < mlen; i++)
            msg[i] =
                (uint8_t)(0x31u + (unsigned)i * 13u + (unsigned)t * 7u);
        memcpy(badmsg, msg, mlen);
        size_t bad_mlen = mlen;
        if (mlen > 0)
            badmsg[mlen / 2] ^= 0x80u;
        else {
            badmsg[0] = 0xa5u;
            bad_mlen = 1;
        }

        if (crypto_sign_keypair(pk, sk) != 0 ||
            crypto_sign_signature(sig, &siglen, msg, mlen, sk) != 0) {
            printf("FAIL setup case=%d\n", t);
            return 1;
        }
        if (crypto_sign_verify(sig, siglen, msg, mlen, pk) == 0)
            bitmap |= 1u;
        if (crypto_sign_verify(sig, siglen, badmsg, bad_mlen, pk) != 0)
            bitmap |= 2u;
        memcpy(malformed, sig, sizeof malformed);
        malformed[CHALLENGE_PACKEDBYTES] = 0xff;
        malformed[CHALLENGE_PACKEDBYTES + 1] = 0xff;
        if (crypto_sign_verify(malformed, sizeof malformed, msg, mlen,
                               pk) != 0)
            bitmap |= 4u;

        case_hash = fnv_bytes(case_hash, pk, sizeof pk);
        case_hash = fnv_bytes(case_hash, sk, sizeof sk);
        case_hash = fnv_bytes(case_hash, sig, siglen);
        case_hash = fnv_u64(case_hash, siglen);
        case_hash = fnv_u64(case_hash, bitmap);
        printf(
            "CASE mode=%d idx=%d mlen=%zu siglen=%zu bitmap=0x%x "
            "case_hash=%llu",
            LITHIUM_MODE, t, mlen, siglen, bitmap,
            (unsigned long long)case_hash);
        case_hash = unpack_summary(case_hash, "valid", sig);
        case_hash = unpack_summary(case_hash, "malformed", malformed);
        printf(" summary_hash=%llu\n", (unsigned long long)case_hash);

        transcript = fnv_u64(transcript, case_hash);
        if (siglen != CRYPTO_BYTES || bitmap != 0x7u)
            ok = 0;
    }

    printf("%s diff_fuzz mode=%d cases=%d transcript_hash=%llu\n",
           ok ? "PASS" : "FAIL", LITHIUM_MODE, DIFF_FUZZ_N,
           (unsigned long long)transcript);
    return ok ? 0 : 1;
}
