#ifndef COPA_KEYGEN_H
#define COPA_KEYGEN_H
#include <stdint.h>
#include "copa_params.h"
#include "copa_matrix.h"
#include "copa_irs64.h"
typedef struct {
    uint8_t pk[COPA_PK_BYTES];
    uint8_t key[32];
    uint8_t hpk[64];
    int32_t s[COPA_K][COPA_N] __attribute__((aligned(64)));
    copa_irs64_context irs_context;
    uint32_t secret_ntt[COPA_K][2][COPA_N] __attribute__((aligned(64)));
    int64_t autocorrelation[COPA_N];
    copa_public_ntt public_ntt; /* runtime cache; absent from canonical KAT serialization */
} copa_secret_key;
typedef struct { uint64_t attempts,noninvertible,norm_reject,lower_norm_reject,packing_failure,noninvertible_even,noninvertible_odd,partial_norm_reject,partial_lower_norm_reject; } copa_keygen_stats;
void copa_expand_a(int32_t a[COPA_M][COPA_ELL][COPA_N],const uint8_t rho[32]);
void copa_expand_s(int32_t s[1+COPA_ELL][COPA_N],const uint8_t seed[32]);
int copa_keygen_raw(uint8_t pk[COPA_PK_BYTES],copa_secret_key *sk,copa_keygen_stats *stats);
int copa_keygen(uint8_t pk[COPA_PK_BYTES],copa_secret_key *sk,copa_keygen_stats *stats);
#endif
