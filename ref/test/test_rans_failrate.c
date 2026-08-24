/* test_rans_failrate.c - statistical smoke test for rANS pack overflow.
 *
 * The 2^-100 target is proved by tools/SigSize.py.  This test catches
 * model drift at the implementation level: the old mode-180 table had an
 * about 11% pack overflow rate, which this test rejects with overwhelming
 * probability.
 */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "params.h"
#include "sign.h"

#ifndef RANS_STATS_N
#    define RANS_STATS_N 512
#endif

extern uint64_t lithium_test_rans_pack_calls;
extern uint64_t lithium_test_rans_pack_fails;

int main(void)
{
    uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
    uint8_t sig[CRYPTO_BYTES];
    uint8_t msg[64];
    size_t siglen = 0;

    memset(msg, 0, sizeof msg);
    if (crypto_sign_keypair(pk, sk) != 0) {
        printf("FAIL test_rans_failrate mode=%d keypair\n", LITHIUM_MODE);
        return 1;
    }

    lithium_test_rans_pack_calls = 0;
    lithium_test_rans_pack_fails = 0;
    for (unsigned i = 0; i < RANS_STATS_N; i++) {
        for (unsigned j = 0; j < 8; j++)
            msg[j] = (uint8_t)(i >> (8 * j));
        if (crypto_sign_signature(sig, &siglen, msg, sizeof msg, sk) !=
                0 ||
            siglen != CRYPTO_BYTES) {
            printf("FAIL test_rans_failrate mode=%d sign i=%u\n",
                   LITHIUM_MODE, i);
            return 1;
        }
        if (crypto_sign_verify(sig, siglen, msg, sizeof msg, pk) != 0) {
            printf("FAIL test_rans_failrate mode=%d verify i=%u\n",
                   LITHIUM_MODE, i);
            return 1;
        }
    }

    int ok = (lithium_test_rans_pack_fails == 0) &&
             (lithium_test_rans_pack_calls >= RANS_STATS_N);
    printf(
        "%s test_rans_failrate mode=%d signs=%d pack_calls=%llu "
        "pack_fails=%llu reserved=%d sig=%d\n",
        ok ? "PASS" : "FAIL", LITHIUM_MODE, RANS_STATS_N,
        (unsigned long long)lithium_test_rans_pack_calls,
        (unsigned long long)lithium_test_rans_pack_fails,
        RANS_RESERVED_BYTES, CRYPTO_BYTES);
    return ok ? 0 : 1;
}
