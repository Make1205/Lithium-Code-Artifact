/*
 * irs.h - Iterative Rejection Sampling for IRS (exp-recurrence,
 * constant-time).
 *
 * Per Lithium.md (RejectSample) + the IRS rejection-sampling design: for
 * each nonzero coefficient x^j of the challenge c, the secret column v =
 * x^j * s is added to z with a sign chosen by the three-region (f,g)
 * bimodal rule, so that the emitted z = y + s c' follows D_{Z,r}
 * independent of s. The overall reject probability is ~2^-101, so the
 * abort branch is dropped (always proceed).
 */
#ifndef LITHIUM_H
#define LITHIUM_H

#include <stdint.h>

#include "params.h"
#include "polyvec.h"
#include "symmetric.h"

/* v = in * x^j in R = Z[x]/(x^n+1).  Exposed for testing. */
void poly_cyclic_shift(poly *out, const poly *in, unsigned int j);

/* Per-key gamma = ||s||^2 / (2 r^2) in Q52 fixed point. */
int64_t compute_gamma_q52(const poly sk[K]);

/* Per-signature series-coefficient table ctab[k] = (-1)^k e^{-k^2 gamma}
 * (Q62), built once from gamma_q52 and reused by decide (table +
 * Horner form of the series S).  ctab must hold REJ_T entries. */
void precompute_c(int64_t *ctab, int64_t gamma_q52);

/* Pure three-region decision from ip=<z,v>, per-key gamma_q52, the
 * per-signature ctab, and uniform r63.  Returns flag_minus (1 => z-v).
 * Exposed for the Python cross-check (tools/irs_ref.py). */
int decide(int64_t ip, int64_t gamma_q52, const int64_t *ctab,
           uint64_t r63);

/* One rejection step on the K-vectors z (in/out) and v (= x^j s).
 *   gamma_q52 : per-key constant; ctab : per-signature table; rng : SHAKE
 * state. Returns the sign decision (1 => z-v branch, 0 => z+v branch);
 * updates z. */
int step(poly z[K], const poly v[K], int64_t gamma_q52,
         const int64_t *ctab, keccak_state *rng);

/* Full signing core: z := y + s c' for the challenge c (binary,
 * weight<=kappa). z   : in = y (K polys), out = response. c   : challenge
 * poly (coeffs in {0,1}). sk  : secret [f, s0, e] (K polys). rng : SHAKE
 * state providing per-step randomness. sign_bits (optional, may be NULL)
 * receives one decision bit per nonzero c. */
void sign_irs(poly z[K], const poly *c, const poly sk[K],
              uint8_t *sign_bits, keccak_state *rng);

#endif /* LITHIUM_H */
