#ifndef LITHIUM_REDUCE_H
#define LITHIUM_REDUCE_H

#include <stdint.h>

#include "ntt.h" /* MONT / QINV / ... (NTT-domain Montgomery constants) */
#include "params.h"

/* ============================================================
 * uint16 Montgomery core, R = 2^16  (NTT engine; AVX2 16-bit lanes)
 * ============================================================
 * Per-mode constants come from the generated block in ntt.h:
 *
 *   MONT     = 2^16 mod q
 *   QINV     = -q^{-1} mod 2^16        (paired with the +m*q
 * reduce below) INVNTT_F = MONT^2 / zetas_size mod q   (final
 * invntt scale, fuses R^2)
 *
 * Bookkeeping (Kyber-style, see ntt.c): twiddles zetas[] are in the
 * Montgomery domain (zeta*R mod q).  fqmul16(x,c) = x*c*R^{-1} mod q, so a
 * Montgomery twiddle times a normal operand yields a normal product -> the
 * forward NTT is normal-domain in, normal-domain out.  Plain products
 * (basemul/pointwise) come out owing one R (value = product*R^{-1});
 * invntt_tomont's R^2 scale pays that R back so a*b lands in the normal
 * domain.
 *
 * Representation is canonical [0,q): EVERY add/sub is reduced
 * (addm16/subm16), because 2q can exceed 2^16 and a bare uint16 addition
 * would wrap.
 */

/* a*2^-16 mod q in [0,q), valid for 0 <= a < q*2^16 (e.g. a = product of
 * two coefficients in [0,q), since q < 2^16 => (q-1)^2 < q*2^16).  The a +
 * m*q sum can reach ~2qR < 2^33, so it is accumulated in uint64. */
static inline uint16_t montgomery_reduce16(uint32_t a)
{
    uint16_t m =
        (uint16_t)(a * (uint32_t)QINV); /* low 16 of a*(-q^{-1}) */
    uint32_t t =
        (uint32_t)((a + (uint64_t)m * Q) >> 16); /* exact; 0 <= t < 2q */
    return (uint16_t)(t >= (uint32_t)Q ? t - Q : t);
}

static inline uint16_t fqmul16(uint16_t a, uint16_t b)
{
    return montgomery_reduce16((uint32_t)a * b);
}

/* modular add/sub keeping the canonical [0,q) representation. */
static inline uint16_t addm16(uint16_t a, uint16_t b)
{
    uint32_t s =
        (uint32_t)a + b; /* up to 2q-2; widened to dodge uint16 wrap */
    return (uint16_t)(s >= (uint32_t)Q ? s - Q : s);
}
static inline uint16_t subm16(uint16_t a, uint16_t b)
{
    return (uint16_t)(a >= b
                          ? a - b
                          : a + Q - b); /* subm16(0,z) = q-z = -z mod q */
}

/* ============================================================
 * signed centered / mod-2q helpers (scheme domain; int32, unchanged)
 * ============================================================ */

int32_t reduce32(int32_t a); /* a mod q in (-q, q) */

int32_t caddq(int32_t a); /* conditional + q  */

int32_t freeze(int32_t a); /* standard representative in [0, q) */

int32_t caddq2(int32_t a); /* conditional + 2q */

int32_t reduce_mod_2q(int32_t a); /* canonical residue in [0, 2q) */

#endif
