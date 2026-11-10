/* test_sign.c - end-to-end gate: keygen -> sign -> verify accept; tamper
 * -> reject. */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "packing.h"
#include "params.h"
#include "sign.h"

static int poly_eq(const poly *a, const poly *b)
{
    for (int k = 0; k < N; k++)
        if (a->coeffs[k] != b->coeffs[k])
            return 0;
    return 1;
}

static int compact_roundtrip_boundary(void)
{
    uint8_t sig[CRYPTO_BYTES];
    poly c, c2, h[M], h2[M];
    poly z1[ELL + 1], z2[ELL + 1];

    memset(&c, 0, sizeof c);
    memset(h, 0, sizeof h);
    memset(z1, 0, sizeof z1);
    for (int k = 0; k < N; k++)
        c.coeffs[k] = (k & 3) == 0;

    const int32_t heads[3] = {
        RANS_ZH_LO,
        0,
        RANS_ZH_LO + RANS_ZH_N - 1,
    };
    const int32_t lows[4] = {
        -(KSPLIT >> 1),
        -1,
        0,
        (KSPLIT >> 1) - 1,
    };
    const int32_t hints[3] = {
        RANS_HINT_LO,
        0,
        RANS_HINT_LO + RANS_HINT_N - 1,
    };
    for (int i = 0; i < ELL + 1; i++) {
        for (int j = 0; j < 12 && j < N; j++) {
            int32_t head = heads[j & 3 ? 1 : (j >> 2)];
            z1[i].coeffs[j] = head * KSPLIT + lows[j & 3];
        }
    }
    for (int i = 0; i < M; i++)
        for (int j = 0; j < 9 && j < N; j++)
            h[i].coeffs[j] = hints[j & 3 ? 1 : (j >> 2)];

    if (pack_sig(sig, &c, h, z1) != 0)
        return 0;
    memset(&c2, 0xA5, sizeof c2);
    memset(h2, 0xA5, sizeof h2);
    memset(z2, 0xA5, sizeof z2);
    if (unpack_sig(&c2, h2, z2, sig) != 0)
        return 0;
    if (!poly_eq(&c, &c2))
        return 0;
    for (int i = 0; i < M; i++)
        if (!poly_eq(&h[i], &h2[i]))
            return 0;
    for (int i = 0; i < ELL + 1; i++)
        if (!poly_eq(&z1[i], &z2[i]))
            return 0;
    return 1;
}

int main(void)
{
    uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
    uint8_t sig[CRYPTO_BYTES];
    uint8_t sm[CRYPTO_BYTES + sizeof("IRS end-to-end test message "
                                     "keygen/sign/verify path 0123")];
    uint8_t opened[sizeof(
        "IRS end-to-end test message keygen/sign/verify path 0123")];
    size_t siglen = 0;
    size_t smlen = 0, opened_len = 0;
    uint8_t msg[56] =
        "IRS end-to-end test message keygen/sign/verify path 0123";

    int ok_accept = 0, ok_reject_msg = 0, ok_reject_sig = 0;
    int ok_attached = 0, ok_reject_attached = 0;
    int ok_reject_rlen = 0, ok_reject_pad = 0, ok_zero_pad = 0;
    int rlen_trials = 0, pad_trials = 0, runs = 5;
    for (int t = 0; t < runs; t++) {
        if (crypto_sign_keypair(pk, sk) != 0) {
            printf("FAIL keypair\n");
            return 1;
        }
        memset(sig, 0xA5, sizeof sig);
        if (crypto_sign_signature(sig, &siglen, msg, sizeof msg, sk) !=
            0) {
            printf("FAIL sign\n");
            return 1;
        }
        if (crypto_sign_verify(sig, siglen, msg, sizeof msg, pk) == 0)
            ok_accept++;

        if (crypto_sign(sm, &smlen, msg, sizeof msg, sk) == 0 &&
            smlen == CRYPTO_BYTES + sizeof msg &&
            crypto_sign_open(opened, &opened_len, sm, smlen, pk) == 0 &&
            opened_len == sizeof msg &&
            memcmp(opened, msg, sizeof msg) == 0)
            ok_attached++;
        sm[0] ^= 0x01;
        if (crypto_sign_open(opened, &opened_len, sm, smlen, pk) != 0)
            ok_reject_attached++;
        sm[0] ^= 0x01;

        uint8_t saved = msg[0];
        msg[0] ^= 0x01;
        if (crypto_sign_verify(sig, siglen, msg, sizeof msg, pk) != 0)
            ok_reject_msg++;
        msg[0] = saved;

        /* flip a byte inside the rANS block (challenge + 2-byte len + a
         * few bytes) */
        size_t zoff = CHALLENGE_PACKEDBYTES + 2 + 6;
        sig[zoff] ^= 0x10;
        if (crypto_sign_verify(sig, siglen, msg, sizeof msg, pk) != 0)
            ok_reject_sig++;
        sig[zoff] ^= 0x10;

        size_t lenoff = CHALLENGE_PACKEDBYTES;
        size_t rlen = (size_t)sig[lenoff] | ((size_t)sig[lenoff + 1] << 8);
        int zero_pad = 1;
        for (size_t pad = rlen; pad < (size_t)RANS_RESERVED_BYTES; pad++)
            zero_pad &= (sig[lenoff + 2 + pad] == 0);
        if (zero_pad)
            ok_zero_pad++;
        if (rlen < RANS_RESERVED_BYTES) {
            size_t pad_off = lenoff + 2 + rlen;
            uint8_t saved0 = sig[lenoff];
            uint8_t saved1 = sig[lenoff + 1];
            uint8_t saved_pad = sig[pad_off];
            size_t new_rlen = rlen + 1;

            pad_trials++;
            sig[pad_off] = 0x42;
            if (crypto_sign_verify(sig, siglen, msg, sizeof msg, pk) != 0)
                ok_reject_pad++;
            sig[pad_off] = saved_pad;

            rlen_trials++;
            sig[lenoff] = (uint8_t)(new_rlen & 0xff);
            sig[lenoff + 1] = (uint8_t)((new_rlen >> 8) & 0xff);
            sig[pad_off] = 0x42;
            if (crypto_sign_verify(sig, siglen, msg, sizeof msg, pk) != 0)
                ok_reject_rlen++;
            sig[lenoff] = saved0;
            sig[lenoff + 1] = saved1;
            sig[pad_off] = saved_pad;
        }
    }

    int ok_roundtrip = compact_roundtrip_boundary();
    int ok = (ok_accept == runs) && (ok_reject_msg == runs) &&
             (ok_reject_sig == runs) && (ok_attached == runs) &&
             (ok_reject_attached == runs) && (rlen_trials > 0) &&
             (ok_reject_rlen == rlen_trials) && (pad_trials > 0) &&
             (ok_reject_pad == pad_trials) && (ok_zero_pad == runs) &&
             ok_roundtrip;
    printf(
        "%s test_sign mode=%d  accept=%d/%d reject_msg=%d/%d "
        "reject_sig=%d/%d attached=%d/%d reject_attached=%d/%d "
        "reject_rlen=%d/%d reject_pad=%d/%d zero_pad=%d/%d "
        "compact_roundtrip=%d  "
        "pk=%d sk=%d sig_raw=%zu (CRYPTO_BYTES target=%d)\n",
        ok ? "PASS" : "FAIL", LITHIUM_MODE, ok_accept, runs, ok_reject_msg,
        runs, ok_reject_sig, runs, ok_attached, runs, ok_reject_attached,
        runs, ok_reject_rlen, rlen_trials, ok_reject_pad, pad_trials,
        ok_zero_pad, runs, ok_roundtrip, CRYPTO_PUBLICKEYBYTES,
        CRYPTO_SECRETKEYBYTES, siglen, CRYPTO_BYTES);
    return ok ? 0 : 1;
}
