#ifndef COPA_POLY_INV_H
#define COPA_POLY_INV_H
#include <stdint.h>
#include "ntt_crt.h"
/* Invert in Z_q[X]/(X^256+1), with public q in {7680,8768,5120,6144}.
 * Inputs have |coefficient| <= 2^18; output is canonical modulo q.
 * Returns 1 on success, 0 and an all-zero output otherwise. Can alias input.
 * No retry/early exit depends on f. */
int copa_poly_inv(int32_t out[COPA_N], const int32_t f[COPA_N], uint32_t q);
#endif
