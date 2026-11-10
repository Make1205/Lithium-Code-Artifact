/*
 * rounding.h - IRS mod-2q HighBits / LSB / hint helpers (Lithium.md Fig.4
 * / A.1).
 *
 * Coefficient-level primitives; vector loops live in sign.c.
 *
 *   HighBits(w) bucket = round(w / tau) mod W1_MAX,  W1_MAX = 2(q-1)/tau.
 *   HighBits value     = bucket * tau.
 *   hint               = (bucket(w) - bucket(wtilde)) centered mod W1_MAX.
 *   UseHint            : bucket(w) = (bucket(wtilde) + hint) mod W1_MAX.
 */
#ifndef LITHIUM_ROUNDING_H
#define LITHIUM_ROUNDING_H

#include <stdint.h>

#include "params.h"

int32_t reduce2q(int32_t a); /* -> [0, 2q) */

int32_t highbits(int32_t w);         /* -> bucket in [0, W1_MAX) */
int32_t highbits_reduced(int32_t x); /* highbits for x already in [0,2q) */

int32_t hbvalue(int32_t bucket); /* -> bucket * tau */

int32_t lsb(int32_t a); /* LSB of (a mod 2q) */

int32_t makehint(int32_t w, int32_t wtilde); /* centered bucket diff */
int32_t makehint_hi(
    int32_t hi_w,
    int32_t wtilde); /* makehint w/ precomputed HighBits(w) */

int32_t usehint(int32_t wtilde, int32_t hint); /* -> recovered bucket */

#endif /* LITHIUM_ROUNDING_H */
