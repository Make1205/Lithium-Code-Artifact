#ifndef NTT_REF_H
#define NTT_REF_H

#include <stdint.h>

#define NTT_Q 45569
#define NTT_N 256

extern uint16_t ntt_ref_R2;                       /* R^2 mod q */

void     ntt_ref_init(void);
uint16_t ntt_ref_fqmul(uint16_t a, uint16_t b);

void ntt_ref(uint16_t r[256]);                    /* forward, in place */
void invntt_tomont_ref(uint16_t r[256]);          /* inverse + to-Montgomery, in place */
void pointwise_ref(const uint16_t a[256], const uint16_t b[256], uint16_t c[256]);
void polymul_schoolbook(const uint16_t a[256], const uint16_t b[256], uint16_t c[256]);

/* the AVX2 assembly entry points and their constant tables.
 * Twiddles use the Seiler precompute: each butterfly consumes a (zl,zh) pair,
 * so the zeta tables hold 2x16-lane vectors per butterfly and z0/z0inv/ninv are
 * 32-int16 (zl[16] then zh[16]).  qdata is {q, 2^16-q, 0} (qinv folded into zl). */
extern const int16_t ntt_qdata[48];
extern const int16_t ntt_qinv[16];
extern const int16_t ntt_z0[32], ntt_z0inv[32], ntt_ninv[32];
extern const int16_t ntt_zetas_fwd[112][16];
extern const int16_t ntt_zetas_inv[112][16];
void ntt_avx   (int16_t *poly, const int16_t *qdata, const int16_t (*ztab)[16], const int16_t *z0,    const int16_t *ninv);
void invntt_tomont_avx(int16_t *poly, const int16_t *qdata, const int16_t (*ztab)[16], const int16_t *z0inv, const int16_t *ninv);
void pointwise_avx(int16_t *c, const int16_t *a, const int16_t *b, const int16_t *qdata);
/* nttunpack: standard-order int32 samples (in [0,q)) -> avx NTT slot layout
 * (Kyber poly_nttunpack analogue: forward shuffle network, no butterflies). */
void nttunpack_avx(int16_t *dst, const int32_t *src);

#endif
