/*
 * speed_profile.c - per-component time + randomness breakdown of
 * crypto_sign.
 *
 * Build with the profiler enabled:  make profile   (or: make AVX2EXP=1
 * profile) Runs N real signatures (incl. rejection retries) and prints,
 * per signature:
 *   - cycles per component (setup / expand_A0 /
 * sample_y{shake,BaseSampler, ApproxExp,final} / commitment{ntt
 * fwd,pw,inv} / highbits / challenge / sign_irs{exp} / normcheck /
 * makehint / pack), as % of the sign total;
 *   - SHAKE bytes squeezed per context and the Gaussian consumption split.
 * See prof.h / prof.c.
 */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "params.h"
#include "prof.h"
#include "randombytes.h"
#include "sign.h"

#define NSIG 2000
#define NKG 500 /* keygen is costlier + uses randombytes; fewer iters */

int main(void)
{
    uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES],
        sig[CRYPTO_BYTES], msg[32];
    size_t siglen = 0;
    randombytes(msg, sizeof msg);
    if (crypto_sign_keypair(pk, sk) != 0) {
        printf("keypair FAIL\n");
        return 1;
    }

#if defined(PROF_TIME) || defined(PROF_RAND)
    char title[32];
    snprintf(title, sizeof title, "mode=%d%s", LITHIUM_MODE,
#    ifdef AVX2_EXP
             " (float exp)");
#    else
             "");
#    endif

    /* --- keygen --- (randombytes per call; rejection retries included) */
    for (int i = 0; i < 16; i++)
        crypto_sign_keypair(pk, sk);
    prof_reset();
    for (int i = 0; i < NKG; i++)
        crypto_sign_keypair(pk, sk);
    prof_report(title, NKG);

    /* --- sign --- (fresh key, then NSIG signatures incl. retries) */
    if (crypto_sign_keypair(pk, sk) != 0)
        return 1;
    for (int i = 0; i < 64; i++) /* warmup (caches, freq) */
        crypto_sign_signature(sig, &siglen, msg, sizeof msg, sk);
    prof_reset();
    for (int i = 0; i < NSIG; i++)
        crypto_sign_signature(sig, &siglen, msg, sizeof msg, sk);
    prof_report(title, NSIG);

    /* --- verify --- (one valid signature, verified NSIG times) */
    crypto_sign_signature(sig, &siglen, msg, sizeof msg, sk);
    for (int i = 0; i < 64; i++)
        crypto_sign_verify(sig, siglen, msg, sizeof msg, pk);
    prof_reset();
    for (int i = 0; i < NSIG; i++)
        crypto_sign_verify(sig, siglen, msg, sizeof msg, pk);
    prof_report(title, NSIG);
#else
    (void)siglen;
    printf(
        "speed_profile: build with `make profile` "
        "(PROF_TIME/PROF_RAND)\n");
#endif
    return 0;
}
