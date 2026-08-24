#include "poly.h"

void poly_reduce(poly *r) {
    reduce_avx(r->coeffs);
}
