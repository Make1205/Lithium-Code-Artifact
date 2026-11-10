#ifndef COPA_NTT_CRT_H
#define COPA_NTT_CRT_H
#include <stddef.h>
#include <stdint.h>
#define COPA_N 256
/* Exact negacyclic product over Z[X]/(X^degree+1), degree a power of two
 * <= 256. Inputs must have |coefficient| <= 2^18. The resulting absolute
 * coefficient bound is < P1*P2/2. Outputs may alias neither input.
 * The transform schedule and reduction schedule depend only on degree.
 */
void copa_ntt_crt_mul(int64_t *out, const int32_t *a, const int32_t *b,
                      size_t degree);
/* Prepared bit-reversed negacyclic transforms for public/key caches and matrix products.
 * Forward residues are ordinary; pointwise Montgomery multiplication adds
 * R^-1, compensated by the inverse scale. The layout is shared by all backends; in/out degree remains public. */
void copa_ntt_forward(uint32_t out[2][COPA_N],const int32_t *in,size_t degree);
void copa_ntt_dot(uint32_t out[2][COPA_N],const uint32_t a[][2][COPA_N],const uint32_t b[][2][COPA_N],unsigned count,size_t degree);
void copa_ntt_dot_reversed(uint32_t out[2][COPA_N],const uint32_t a[][2][COPA_N],const uint32_t b[][2][COPA_N],unsigned count,size_t degree);
void copa_ntt_inverse_crt(int64_t *out,const uint32_t in[2][COPA_N],size_t degree);
void copa_ntt_inverse_modq(int32_t *out,const uint32_t in[2][COPA_N],size_t degree,uint32_t modulus);
/* Canonical reduction, fixed schedule for a public modulus >= 2. */
uint32_t copa_mod_i64(int64_t value, uint32_t modulus);
#endif
