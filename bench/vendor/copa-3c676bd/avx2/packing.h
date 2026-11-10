/*
 * packing.h - serialization for IRS keys and signatures.
 *
 * Compact signature (the on-wire format, fixed size CRYPTO_BYTES):
 *   challenge(n-bit) || uint16 rans_len || rANS(z-head ++ hint) padded to
 *   RANS_RESERVED_BYTES || z-lo (ell+1 components, t bits/coef).
 *   z-head[i][k] = round(z1[i][k] / 2^t); z-lo = z1 - 2^t * z-head (t-bit
 * signed); hint = the m hint polynomials.  pack_sig returns -2 on rANS
 * overflow (signer retries).
 *
 * RAW packing (pack_sig_raw) is kept for module-level validation.
 *
 * Key format note: pk/sk carry bn = NTT_ref(b), not backend-native bh
 * bytes. "NTT_ref" means the exact coefficient order produced by ref/
 * ntt16() after reducing b to [0,q).  This canonical NTT order is the wire
 * format because it is independent of AVX2/AVX-512 assembly layouts.  AVX
 * backends unpack bn and then call nttunpack_avx / nttunpack_avx512 to
 * obtain their native NTT slot order before matrix products.
 */
#ifndef LITHIUM_PACKING_H
#define LITHIUM_PACKING_H

#include <stdint.h>

#include "params.h"
#include "polyvec.h"

/* ---- RAW (validation) signature ---- */
#define HINT_PACKEDBYTES ((N * W1_BITS + 7) / 8)
#define Z1_PACKEDBYTES (N * 2)
#define SIG_RAW_BYTES                               \
    (CHALLENGE_PACKEDBYTES + M * HINT_PACKEDBYTES + \
     (ELL + 1) * Z1_PACKEDBYTES)

/* bn must be canonical NTT wire data: exact ref/ ntt16 order and each
 * coefficient in [0,q).  Do not pass backend-native AVX NTT output here:
 * some mode-180 AVX paths are lazy [0,2q), and AVX slot order is not the
 * serialized order. */
void pack_pk(uint8_t *pk, const uint8_t rho[SEEDBYTES],
             const poly16 bn[M]);
void unpack_pk(uint8_t rho[SEEDBYTES], poly16 bn[M], const uint8_t *pk);

/* Same bn contract as pack_pk: callers must provide poly_ntt_canonical()
 * output, not poly_ntt()/poly_ntt_small() backend-native output. */
void pack_sk(uint8_t *sk, const uint8_t rho[SEEDBYTES],
             const uint8_t tr[TRBYTES], const uint8_t key[SKSEEDBYTES],
             const poly *f0, const poly s0[ELL], const poly e[M],
             const poly16 bn[M]);
void unpack_sk(uint8_t rho[SEEDBYTES], uint8_t tr[TRBYTES],
               uint8_t key[SKSEEDBYTES], poly *f0, poly s0[ELL], poly e[M],
               poly16 bn[M], const uint8_t *sk);

void pack_sig_raw(uint8_t *sig, const poly *c, const poly h[M],
                  const poly z1[ELL + 1]);
void unpack_sig_raw(poly *c, poly h[M], poly z1[ELL + 1],
                    const uint8_t *sig);

/* ---- compact (on-wire) signature ---- */
/* pack_sig reads exactly M hint polynomials and the first
 * ELL+1 z polynomials; callers may pass a read-only prefix view.
 */
int pack_sig(uint8_t sig[CRYPTO_BYTES], const poly *c, const poly h[M],
             const poly z1[ELL + 1]); /* 0 ok, -2 rANS overflow */
int unpack_sig(poly *c, poly h[M], poly z1[ELL + 1],
               const uint8_t sig[CRYPTO_BYTES]);

#endif /* LITHIUM_PACKING_H */
