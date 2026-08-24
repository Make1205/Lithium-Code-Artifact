#ifndef NTT_REF_H
#define NTT_REF_H

#include <stdint.h>

#define NTT_Q 30977
#define NTT_N 256
#define NTT_LEVELS 7          /* INCOMPLETE NTT: stop at len=2, then basemul */

extern uint16_t ntt_ref_R2;   /* R^2 mod q */

void     ntt_ref_init(void);
uint16_t ntt_ref_fqmul(uint16_t a, uint16_t b);

void ntt_ref(uint16_t r[256]);                    /* forward, in place */
void invntt_tomont_ref(uint16_t r[256]);          /* inverse + to-Montgomery, in place */
/* incomplete NTT: 128 degree-2 residue products mod (x^2 - gamma_i) */
void basemul_ref(const uint16_t a[256], const uint16_t b[256], uint16_t c[256]);
void polymul_schoolbook(const uint16_t a[256], const uint16_t b[256], uint16_t c[256]);

/* AVX2 assembly entry points + constant tables.
 * UNSIGNED [0,q) representation with LAZY reduction (csub schedule); the
 * incomplete NTT stops at len=2 and a separate basemul finishes the product.
 * qdata is {q, 2^16-q, 0}.  Twiddles use the Seiler precompute: each butterfly
 * consumes a (zl,zh) pair, so the zeta tables hold 2x16-lane vectors per
 * butterfly; z0/z0inv/scale are 32-int16 (zl[16] then zh[16]).  The inverse
 * scale is 128^-1*R^2, so invntt_tomont maps a bare roundtrip to a*R. */
extern const int16_t ntt_qdata[48];
extern const int16_t ntt_qinv[16];
extern const int16_t ntt_r2[32];
extern const int16_t ntt_z0[32], ntt_z0inv[32], ntt_scale[32];
extern const int16_t ntt_zetas_fwd[96][16];
extern const int16_t ntt_zetas_inv[96][16];
extern const int16_t ntt_basemul_gamma[32][16];
void ntt_avx    (int16_t *poly, const int16_t *qdata, const int16_t (*ztab)[16], const int16_t *z0,    const int16_t *scale);
void invntt_tomont_avx(int16_t *poly, const int16_t *qdata, const int16_t (*ztab)[16], const int16_t *z0inv, const int16_t *scale);
void basemul_avx(int16_t *c, const int16_t *a, const int16_t *b, const int16_t *qdata, const int16_t (*gammas)[16]);
void reduce_avx(int16_t *poly);
/* nttunpack: standard-order int32 samples (in [0,q)) -> avx NTT slot layout
 * (Kyber poly_nttunpack analogue: forward shuffle network, no butterflies). */
void nttunpack_avx(int16_t *dst, const int32_t *src);

/* Opt-in AVX-512BW full-ZMM path.  Built only by the AVX512=1 Makefile targets. */
extern const int16_t ntt512_qdata[64];
extern const int16_t ntt512_scale[64];
extern const int16_t ntt512_qinv[32];
extern const int16_t ntt512_r2[64];
extern const int16_t ntt512_zetas_fwd[56][32];
extern const int16_t ntt512_zetas_inv[56][32];
extern const int16_t ntt512_basemul_gamma[16][32];
void ntt_avx512    (int16_t *poly, const int16_t *qdata, const int16_t (*ztab)[32], const int16_t *scale);
void invntt_tomont_avx512(int16_t *poly, const int16_t *qdata, const int16_t (*ztab)[32], const int16_t *scale);
void basemul_avx512(int16_t *c, const int16_t *a, const int16_t *b, const int16_t *qdata, const int16_t (*gammas)[32]);
void nttunpack_avx512(int16_t *dst, const int32_t *src);

#endif
