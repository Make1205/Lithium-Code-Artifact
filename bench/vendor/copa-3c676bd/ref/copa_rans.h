#ifndef COPA_RANS_H
#define COPA_RANS_H
#include "copa_params.h"
#include <stddef.h>
#include <stdint.h>
/* Frozen codec-v2 layout: state[4] || raw-low bits || rANS bytes || zero pad. */
#define COPA_RANS_BITS 20
#define COPA_RANS_LOWER (UINT32_C(1)<<23)
#define COPA_PAYLOAD (COPA_SIG_BYTES-32)
/* Returns 1 on success; overflow / invalid values return 0. Outputs are
 * fixed-capacity and zero-padded. used includes the four-byte coder state; no length field is stored. */
int copa_encode(uint8_t out[COPA_PAYLOAD],size_t *used,const int32_t z1[1+COPA_ELL][COPA_N],const int32_t h[COPA_M][COPA_N]);
int copa_decode(int32_t z1[1+COPA_ELL][COPA_N],int32_t h[COPA_M][COPA_N],const uint8_t in[COPA_PAYLOAD]);
#endif
