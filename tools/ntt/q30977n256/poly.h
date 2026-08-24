#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "ntt_ref.h"

typedef struct {
    int16_t coeffs[NTT_N];
} poly;

void poly_reduce(poly *r);

#endif
