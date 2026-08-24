#ifndef NTT_REF_H
#define NTT_REF_H

#include <stdint.h>

#define NTT_Q 64513
#define NTT_N 512

extern uint16_t ntt_ref_R2;                       /* R^2 mod q */

void     ntt_ref_init(void);
uint16_t ntt_ref_fqmul(uint16_t a, uint16_t b);

void ntt_ref(uint16_t r[NTT_N]);                    /* forward, in place */
void invntt_tomont_ref(uint16_t r[NTT_N]);          /* inverse + to-Montgomery, in place */
void pointwise_ref(const uint16_t a[NTT_N], const uint16_t b[NTT_N], uint16_t c[NTT_N]);
void polymul_schoolbook(const uint16_t a[NTT_N], const uint16_t b[NTT_N], uint16_t c[NTT_N]);

/* the AVX2 assembly entry points and their constant tables.
 * Twiddles use the Seiler precompute: each butterfly consumes a (zl,zh) pair,
 * so the zeta tables hold 2x16-lane vectors per butterfly and ninv is 32-int16
 * (zl[16] then zh[16]).  qdata is {q, 2^16-q, 0} (qinv folded into zl).
 * cross_fwd = level0(2 pairs)+level1(2 pairs); cross_inv = level1+level0. */
extern const int16_t ntt_qdata[48];
extern const int16_t ntt_ninv[32];
extern const int16_t ntt_cross_fwd[8][16];
extern const int16_t ntt_cross_inv[8][16];
extern const int16_t ntt_zetas_fwd[224][16];
extern const int16_t ntt_zetas_inv[224][16];
void ntt_avx   (int16_t *poly, const int16_t *qdata, const int16_t (*ztab)[16], const int16_t (*cross)[16], const int16_t *ninv);
void invntt_tomont_avx(int16_t *poly, const int16_t *qdata, const int16_t (*ztab)[16], const int16_t (*cross)[16], const int16_t *ninv);
void pointwise_avx(int16_t *c, const int16_t *a, const int16_t *b, const int16_t *qdata);
/* nttunpack: standard-order int32 samples (in [0,q)) -> avx NTT slot layout
 * (Kyber poly_nttunpack analogue: forward shuffle network, no butterflies). */
void nttunpack_avx(int16_t *dst, const int32_t *src);

/* Opt-in AVX-512BW full-ZMM path.  Built only by the AVX512=1 Makefile targets. */
extern const int16_t ntt512_qdata[64];
extern const int16_t ntt512_scale[64];
extern const int16_t ntt512_qinv[32];
extern const int16_t ntt512_r2[64];
extern const int16_t ntt512_zetas_fwd[144][32];
extern const int16_t ntt512_zetas_inv[144][32];
void ntt_avx512      (int16_t *poly, const int16_t *qdata, const int16_t (*ztab)[32], const int16_t *scale);
void invntt_tomont_avx512(int16_t *poly, const int16_t *qdata, const int16_t (*ztab)[32], const int16_t *scale);
void pointwise_avx512(int16_t *c, const int16_t *a, const int16_t *b, const int16_t *qdata);
void nttunpack_avx512(int16_t *dst, const int32_t *src);

#endif
