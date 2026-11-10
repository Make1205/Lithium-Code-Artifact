#ifndef COPA_MATRIX_H
#define COPA_MATRIX_H
#include "copa_params.h"
#include <stdint.h>
typedef uint32_t copa_public_ntt[COPA_M][1+COPA_ELL][2][COPA_N] __attribute__((aligned(64)));
int copa_matrix_expand(copa_public_ntt matrix,const uint8_t pk[COPA_PK_BYTES]);
void copa_matrix_apply_ntt(int32_t out[COPA_M][COPA_N],const copa_public_ntt matrix,const uint32_t input[][2][COPA_N]);
void copa_matrix_apply(int32_t out[COPA_M][COPA_N],const copa_public_ntt matrix,const int32_t input[1+COPA_ELL][COPA_N]);
#ifdef COPA_VERIFY31
/* Prepared-Verify cache uses the same storage for a residue/Shoup pair. */
int copa_verify31_expand(copa_public_ntt matrix,const uint8_t pk[COPA_PK_BYTES]);
void copa_verify31_apply(int32_t out[COPA_M][COPA_N],const copa_public_ntt matrix,const int32_t z[1+COPA_ELL][COPA_N]);
#endif
#endif
