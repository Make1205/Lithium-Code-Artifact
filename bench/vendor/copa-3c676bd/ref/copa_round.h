#ifndef COPA_ROUND_H
#define COPA_ROUND_H
#include "copa_params.h"
#include <stdint.h>
/* Exact for x in [0,q): ceil(2^24/alpha) has enough reciprocal bits
 * for every installed alpha. Exhaustively checked for all four domains. */
static inline int32_t copa_highbits(int32_t x){
 uint32_t h=(((uint32_t)x+(COPA_ALPHA>>1))*COPA_ALPHA_RECIP)>>24;
 return (int32_t)(h-COPA_HINT_LEVELS*(h>=COPA_HINT_LEVELS));
}
#endif
