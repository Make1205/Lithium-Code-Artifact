#ifndef COPA_NIST_RNG_H
#define COPA_NIST_RNG_H
#include <stddef.h>
#include <stdint.h>
/* Test-only AES-256 CTR_DRBG without derivation function, 48-byte entropy.
 * SP 800-90A Rev.1 sections 10.2.1.2, 10.2.1.3.1 and 10.2.1.5.1.
 * Matches the NIST PQC KAT rng.c initialization / generate convention.
 */
void copa_nist_seed(const uint8_t entropy[48]);
void randombytes(uint8_t *out,size_t length);
void copa_nist_selftest(void);
#endif
