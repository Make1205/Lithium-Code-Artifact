/*
 * polyvec.h - polynomial arithmetic + NTWE matrix/vector operations for
 * IRS.
 *
 * The poly / poly16 types, ring-arithmetic declarations, and NTWE
 * matrix/vector declarations live here together.
 *
 * Two representations:
 *   poly   - signed int32 scheme domain (sampling, rejection, mod-2q lift,
 * norm, packing).  Coefficients may be negative / centered. poly16 -
 * unsigned uint16 NTT / mod-q domain, canonical [0,q).  This is the
 *            AVX2-aligned representation; the NTT engine (ntt.c) operates
 * on it. Conversion happens at the NTT boundary (poly_ntt /
 * poly_invntt_tomont).
 *
 * Public key: (rho, bn), where bn is the canonical serialized NTT-domain
 * vector bn = NTT_ref(b) and b = (A0 s0 + e) f^{-1} mod q in R_q^m.
 * The implicit verification matrix is A = [qj - 2b, 2 A0, 2 I_m]
 * (m x K), K = 1 + ell + m; coefficient-domain b is recovered only
 * conceptually as INTT_ref(bn).
 *
 * Canonical NTT serialization order: the on-wire order for any serialized
 * NTT-domain polynomial is exactly the coefficient order produced by the
 * scalar ref/ ntt16() routine after freeze-to-[0,q).  We choose ref/ order
 * as the wire invariant because it is portable, byte-stable, and
 * independent of AVX2/AVX-512 assembly slot layouts.  AVX backends must
 * import such a canonical polynomial with their nttunpack_* helper before
 * using it in assembly-domain products; nttunpack_* is the no-butterfly
 * forward shuffle ladder satisfying nttunpack_*(NTT_ref(a)) ==
 * NTT_avx*(a).
 *
 * Pre-NTT (Dilithium-style, see docs/design-notes/PreNTT.md): A0 is
 * uniform/public, so it is generated DIRECTLY in the NTT domain from rho
 * (expand_A0_ntt) and never forward-transformed; bn is already serialized
 * as b-hat = NTT_ref(b). The matrix/vector ops therefore take NTT-domain
 * A0h / bh (poly16) inputs, with AVX imports performed at unpack time.
 */
#ifndef LITHIUM_POLYVEC_H
#define LITHIUM_POLYVEC_H

#include <stdint.h>

#include "params.h"

/* ============================================================
 * poly types + ring arithmetic
 * ============================================================ */
typedef struct {
    int32_t coeffs[N];
} poly;

typedef struct {
    uint16_t coeffs[N];
} poly16;

void poly_reduce(poly *a); /* coeffs -> (-q,q) */

void poly_caddq(poly *a); /* conditional + q  */

void poly_freeze(poly *a); /* coeffs -> [0,q)  */

void poly_add(poly *r, const poly *a,
              const poly *b); /* normal domain (int32) */

void poly_sub(poly *r, const poly *a,
              const poly *b); /* normal domain (int32) */

/* signed poly -> backend-native NTT-domain poly16 (freezes to [0,q) then
 * forward NTT).  In ref/ this is the canonical NTT order; in AVX backends
 * this is the assembly slot layout. */
void poly_ntt(poly16 *out, const poly *in);

/* small signed poly -> backend-native NTT-domain poly16.  Precondition:
 * every coefficient is in [-(q-1), q-1].  The signed-to-[0,q)
 * conversion is written as explicit unsigned mask arithmetic (no C-level
 * branch/?:, no divide/Barrett) and is intended for Gaussian/noise/z1
 * inputs whose lifecycle proves this narrow range. */
void poly_ntt_small(poly16 *out, const poly *in);

/* signed poly -> canonical serialized NTT order, defined as the exact ref/
 * ntt16 output order.  Used to create bn = NTT_ref(b) for pk/sk bytes. */
void poly_ntt_canonical(poly16 *out, const poly *in);

/* canonical serialized NTT order -> backend-native NTT layout.  ref/ is
 * the identity; AVX2/AVX-512 call nttunpack_* before matrix products. */
void poly_ntt_import(poly16 *out, const poly16 *in);

/* Unpack M packed b-hat polynomials (QBITS-bit LSB-first fields) straight
 * into canonical uint16 [0,q) in ref-NTT order -- the fused replacement
 * for bitunpack()+freeze() in unpack_pk/unpack_sk.  Each packed
 * coefficient is in [0, 2^QBITS) and 2^QBITS < 2q, so a single conditional
 * subtract of q is byte-exact with freeze() over the whole field domain
 * (KAT-neutral). `packed` points at the first packed poly; stride
 * POLYPK_PACKEDBYTES.  No native NTT-slot shuffle (callers run
 * poly_ntt_import afterwards).  ref/ is scalar; AVX2/AVX-512 vectorize the
 * load + conditional subtract. */
void unpack_pk_bn(poly16 bn[M], const uint8_t *packed);

/* NTT-domain poly16 -> signed poly in [0,q) (inverse NTT + Montgomery
 * scale). */
void poly_invntt_tomont(poly *out, const poly16 *in);

/* NTT-domain pointwise product (per-coeff full NTT / basemul incomplete
 * NTT). Output owes one Montgomery R (paid back by a later
 * poly_invntt_tomont). */
void poly_pointwise_montgomery(poly16 *r, const poly16 *a,
                               const poly16 *b);

/* NTT-domain modular add (canonical [0,q)). */
void poly16_add(poly16 *r, const poly16 *a, const poly16 *b);

/* Multiply in R_q (normal domain): r = a*b mod (x^n+1, q), output in
 * [0,q). */
void poly_mul(poly *r, const poly *a, const poly *b);

/* Uniform sample mod q from SHAKE-128(seed || nonce). */
void poly_uniform(poly *a, const uint8_t seed[SEEDBYTES], uint16_t nonce);

/* 4-way uniform expansion (canonical layout for AVX2 KAT-compat): all four
 * lanes are real (caller handles remainders with scalar poly_uniform). */
void poly_uniform_4x(poly *a0, poly *a1, poly *a2, poly *a3,
                     const uint8_t seed[SEEDBYTES], uint16_t nonce0,
                     uint16_t nonce1, uint16_t nonce2, uint16_t nonce3);

/* 8-way uniform expansion (AVX-512 keccak8x; KAT-neutral wrt the 4->8
 * widening -- each lane reproduces poly_uniform(seed, nonce_j)). */
void poly_uniform_8x(poly *a0, poly *a1, poly *a2, poly *a3, poly *a4,
                     poly *a5, poly *a6, poly *a7,
                     const uint8_t seed[SEEDBYTES], uint16_t nonce0,
                     uint16_t nonce1, uint16_t nonce2, uint16_t nonce3,
                     uint16_t nonce4, uint16_t nonce5, uint16_t nonce6,
                     uint16_t nonce7);

/* A0 expansion as 8 EQUAL coefficient chunks: lane j rejection-samples
 * per_lane uniform [0,q) coeffs from SHAKE128(seed||j) into
 * out[j*per_lane .. (j+1)*per_lane).  Every q < 2^16, so the accepted
 * samples are written straight as uint16_t (no int32 widening).  Used by
 * expand_A0_ntt so the 8-way Keccak engine is ALWAYS fully utilised
 * (M*ELL*N is a multiple of 8 in every mode), with no idle lanes / scalar
 * tail for any M*ELL.  Backends: ref = 8 scalar SHAKE128, avx2 = 2x
 * keccak4x + SIMD reject, avx512 = 1x keccak8x + SIMD reject. */
void poly_uniform_chunkx8(uint16_t *out, int per_lane,
                          const uint8_t seed[SEEDBYTES]);

/* Invert a in R_q.  Returns 0 on success (r = a^{-1}), -1 if a is not
 * invertible (some NTT residue / quadratic norm is zero). */
int poly_invert(poly *r, const poly *a);

/* Centered squared L2 norm (sum of (coeff mod^pm q)^2). */
int64_t poly_sqnorm(const poly *a);

/* ============================================================
 * NTWE matrix/vector operations
 * ============================================================ */
/* Expand A0h in [m][ell] DIRECTLY in the NTT domain (poly16, [0,q)) from
 * rho -- "seed-to-NTT": the uniform [0,q) samples ARE the forward-NTT
 * representation (no ntt16), so A0 := INTT(A0h).  Uniformity is preserved
 * because NTT is a bijection R_q <-> T_q.  ref/ ntt16 output order is the
 * canonical NTT order; AVX2/AVX-512 permute that order with nttunpack_*
 * into their assembly slot layout before multiplication. */
void expand_A0_ntt(poly16 A0h[M][ELL], const uint8_t rho[SEEDBYTES]);

/* b = (A0 s0 + e) * f^{-1} mod q.  Returns 0 on success, -1 if f not
 * invertible. A0h is NTT-domain (from expand_A0_ntt); f is the
 * normal-domain secret poly (= 2 f0 + 1).  The public-key line format
 * stores bn = NTT_ref(b) in canonical ref/ ntt16 order; backends import bn
 * to their native bh layout before calling mat_mul_*(). */
int compute_b(poly b[M], const poly16 A0h[M][ELL], const poly s0[ELL],
              const poly e[M], const poly *f);

/* Commitment w = A x mod 2q  (m polys, coeffs in [0,2q)).
 * x is the length-K vector [x0, xs(ell), xe(m)]; A0h and bh are the
 * NTT-domain public matrix / vector (bh = NTT(b), cached by the caller).
 */
void mat_mul_2q(poly w[M], const poly16 A0h[M][ELL], const poly16 bh[M],
                const poly x[K]);

/* Verifier recompute wtilde' = A1 z1 - q c j mod 2q  (m polys, [0,2q)).
 * z1 = [z0, zs(ell)] (length ell+1); c is the challenge poly; A0h/bh are
 * NTT-domain. */
void mat_mul_z1_2q(poly wt[M], const poly16 A0h[M][ELL],
                   const poly16 bh[M], const poly z1[ELL + 1],
                   const poly *c);

#endif /* LITHIUM_POLYVEC_H */
