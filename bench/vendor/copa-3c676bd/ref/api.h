/*
 * api.h - NIST-style API header for Lithium signature scheme.
 *
 * This header exposes the standard SUPERCOP/NIST entry points
 * (crypto_sign_keypair, crypto_sign_signature, crypto_sign_verify,
 * crypto_sign, crypto_sign_open) with names mangled per LITHIUM_MODE,
 * so that liblithium_nist120_ref.so, liblithium_nist180_ref.so and
 * liblithium_nist260_ref.so can coexist in the same linker namespace.
 *
 * A caller compiled against mode 120 will have their `crypto_sign_keypair`
 * references rewritten by the preprocessor to
 * `lithium_nist120_ref_keypair`, which matches the symbol exported from
 * liblithium_nist120_ref.so.
 */

#ifndef LITHIUM_API_H
#define LITHIUM_API_H

#include <stddef.h>
#include <stdint.h>

#include "params.h"

/* CRYPTO_PUBLICKEYBYTES / CRYPTO_SECRETKEYBYTES / CRYPTO_BYTES come from
 * params.h (included above). */
/* CRYPTO_ALGNAME ("Lithium-120" / "Lithium-180" / "Lithium-260") comes
 * from config.h (included via params.h). */

/* NIST-style entry points; mangled to match the library symbols. */
int crypto_sign_keypair(uint8_t *pk, uint8_t *sk);

int crypto_sign_signature(uint8_t *sig, size_t *siglen, const uint8_t *m,
                          size_t mlen, const uint8_t *sk);

int crypto_sign_verify(const uint8_t *sig, size_t siglen, const uint8_t *m,
                       size_t mlen, const uint8_t *pk);

int crypto_sign(uint8_t *sm, size_t *smlen, const uint8_t *m, size_t mlen,
                const uint8_t *sk);

int crypto_sign_open(uint8_t *m, size_t *mlen, const uint8_t *sm,
                     size_t smlen, const uint8_t *pk);

#endif /* LITHIUM_API_H */
