/*
 * ntt.c - Forward / inverse NTT for IRS (unsigned 16-bit, R = 2^16).
 *
 * Coefficients are canonical uint16 in [0,q).  Because 2q can exceed 2^16,
 * NO lazy accumulation is possible: every butterfly add/sub is reduced
 * immediately (addm16/subm16).  (A future AVX2 port may reintroduce a lazy
 * [0,2q) schedule for q < 2^15 -- see tools/ntt -- but the reference stays
 * fully reduced.)
 *
 * Montgomery bookkeeping (Kyber-style; fqmul16(x,c) = x*c*R^{-1} mod q):
 *   - zetas[] are Montgomery-domain twiddles (zeta*R mod q).  A butterfly
 * does t = fqmul16(zeta_mont, a[j+len]) = zeta*a[j+len], so the FORWARD
 * NTT maps normal-domain input to normal-domain output (no Montgomery
 * factor added).
 *   - "Pointwise" multiplication (basemul / per-coeff) is a plain fqmul16
 *     product of two normal-domain values, so its output OWES one R
 *     (value = product*R^{-1}).
 *   - invntt_tomont's final scale INVNTT_F = MONT^2/zetas_size
 * fuses R^2 (not R): fqmul16(.,f) applies net *R/size, which pays back the
 * R the basemul owed AND normalizes by 1/size, so the product lands in the
 * normal domain with no leftover R.  (size = N full / N/2 incomplete =
 * 2^levels.) A bare ntt+invntt round-trip with no basemul ends one R high
 * = Montgomery domain -- hence the name "_tomont"; the tests strip one
 * montgomery_reduce16.
 *
 * AVX2 note: the last invntt level can fold f into that level's fqmul
 * twiddle so half the coeffs get the final scaling for free, leaving only
 * the pass-through half to an explicit fqmul16(.,f).  The scalar reference
 * keeps a clear uniform final loop (like Kyber's C reference); the fusion
 * is for the AVX2 port.
 *
 * Two flavours selected by BASE_DEG (full NTT vs incomplete,
 * basecase deg 2) share the single zetas[] table (size
 * NTT_ZETAS_SIZE): butterflies use zetas[1..ZETAS_SIZE-1]; basemul
 * uses zetas[ZETAS_SIZE/2
 * .. ZETAS_SIZE-1]
 * (+/- pairs, negation done as q-zeta in uint16).  Constants:
 * gen_zetas.py.
 */

#include "ntt.h"

#include <stdint.h>

#include "params.h"
#include "reduce.h"

/*************************************************
 * Name:        ntt16
 *
 * Description: In-place forward NTT, uint16 canonical [0,q). Normal-domain
 *in, normal-domain out.  Cooley-Tukey, bit-reversed output ordering.
 *              base_deg = 1: log2(N) layers (full NTT).
 *              base_deg = 2: log2(N)-1 layers (stop at len=1), output is
 *N/2 deg-1 polynomials in basecase rings Z_q[X]/(X^2 - zeta).
 *
 * Arguments:   - uint16_t a[N]: in/out coefficient array
 **************************************************/
void ntt16(uint16_t a[N])
{
    unsigned int len, start, j, k;
    uint16_t zeta, t;

    k = 0;
    for (len = N / 2; len >= (unsigned)BASE_DEG; len >>= 1) {
        for (start = 0; start < N; start = j + len) {
            zeta = zetas[++k];
            for (j = start; j < start + len; ++j) {
                t = fqmul16(zeta, a[j + len]);
                a[j + len] = subm16(a[j], t);
                a[j] = addm16(a[j], t);
            }
        }
    }
}

/*************************************************
 * Name:        invntt16_tomont
 *
 * Description: In-place inverse NTT + Montgomery scale, uint16 canonical
 *[0,q). Gentleman-Sande, layers matching the forward transform.  The final
 *              loop multiplies by INVNTT_F = MONT^2/zetas_size
 *(see file header: fuses R^2 so pointwise products return to normal
 *domain).
 *
 * Arguments:   - uint16_t a[N]: in/out coefficient array
 **************************************************/
void invntt16_tomont(uint16_t a[N])
{
    unsigned int start, len, j, k;
    uint16_t t, zeta;
    const uint16_t f = (uint16_t)INVNTT_F;

    k = NTT_ZETAS_SIZE;
    for (len = BASE_DEG; len < N; len <<= 1) {
        for (start = 0; start < N; start = j + len) {
            zeta = subm16(0, zetas[--k]); /* -zeta mod q */
            for (j = start; j < start + len; ++j) {
                t = a[j];
                a[j] = addm16(t, a[j + len]);
                a[j + len] = fqmul16(
                    zeta,
                    subm16(t, a[j + len])); /* (t - a[j+len])*(-zeta) */
            }
        }
    }

    for (j = 0; j < N; ++j)
        a[j] = fqmul16(f, a[j]);
}

#if BASE_DEG == 2
/*************************************************
 * Name:        basemul16
 *
 * Description: One basemul over Z_q[X]/(X^2 - zeta), uint16:
 *                r0 = a0*b0 + a1*b1*zeta
 *                r1 = a0*b1 + a1*b0
 *              All products via fqmul16; output owes one R (see file
 *header).
 *
 * Arguments:   - uint16_t r[2]:        output basecase pair
 *              - const uint16_t a[2]:  first input basecase pair
 *              - const uint16_t b[2]:  second input basecase pair
 *              - uint16_t zeta:        basecase modulus (Montgomery form)
 **************************************************/
static inline void basemul16(uint16_t r[2], const uint16_t a[2],
                             const uint16_t b[2], uint16_t zeta)
{
    uint16_t t;
    t = fqmul16(a[1], b[1]);
    t = fqmul16(t, zeta);
    t = addm16(t, fqmul16(a[0], b[0]));
    r[0] = t;
    t = fqmul16(a[0], b[1]);
    t = addm16(t, fqmul16(a[1], b[0]));
    r[1] = t;
}

/*************************************************
 * Name:        poly16_basemul_montgomery_native
 *
 * Description: basemul over all N/2 basecases.  Paired basecases share a
 *|zeta| with opposite signs: pair (4i,   4i+1) uses  zetas[ZETAS_SIZE/2 +
 *i] pair (4i+2, 4i+3) uses -zetas[ZETAS_SIZE/2 + i]  (= q - zeta)
 *
 * Arguments:   - uint16_t r[N]:        output (NTT domain, owes
 *one R)
 *              - const uint16_t a[N]:  first input (NTT domain)
 *              - const uint16_t b[N]:  second input (NTT domain)
 **************************************************/
void poly16_basemul_montgomery_native(uint16_t r[N], const uint16_t a[N],
                                      const uint16_t b[N])
{
    unsigned int i;
    const unsigned int half = NTT_ZETAS_SIZE / 2;
    for (i = 0; i < N / 4; ++i) {
        uint16_t zeta = zetas[half + i];
        basemul16(&r[4 * i], &a[4 * i], &b[4 * i], zeta);
        basemul16(&r[4 * i + 2], &a[4 * i + 2], &b[4 * i + 2],
                  subm16(0, zeta));
    }
}
#endif /* BASE_DEG == 2 */
