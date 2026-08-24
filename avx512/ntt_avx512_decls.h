/*
 * ntt_avx512_decls.h - extern declarations for the in-tree AVX-512 NTT
 * assembly and its constant tables, switched by LITHIUM_MODE.
 *
 * The assembly + constants are pulled (by explicit path in the Makefile)
 * from tools/ntt/q57601n256/  (mode 120, incomplete base_deg=2, strict)
 *      tools/ntt/q30977n256/  (mode 180, incomplete base_deg=2, basemul
 *                              SELF-canonicalises to [0,q): no reduce
 * site) tools/ntt/q64513n512/  (mode 260, complete base_deg=1, strict) The
 * generated assembly and constant tables are namespaced with
 * LITHIUM_NAMESPACE so multiple Lithium modes/implementations can coexist
 * in one link unit.  These prototypes are copied from the matching
 * tools/ntt/<qset>/ntt_ref.h (the ntt512_* block) so the compiler
 * type-checks the (mode-specific!) basemul argument order.
 *
 * Vs the AVX2 decls (ntt_avx_decls.h): the AVX-512 forward/inverse
 * transforms take the SAME 4-arg signature (poly, qdata, ztab, scale) for
 * ALL modes (no z0 / cross args), the twiddle/const tables are [...][32]
 * (full-ZMM, 32x int16), and there is no reduce_avx512 (mode 180 basemul
 * self-reduces).
 */
#ifndef LITHIUM_NTT_AVX512_DECLS_H
#define LITHIUM_NTT_AVX512_DECLS_H

#include <stdint.h>

#include "config.h"

#define ntt_avx512 LITHIUM_NAMESPACE(ntt_avx512)
#define invntt_tomont_avx512 LITHIUM_NAMESPACE(invntt_tomont_avx512)
#define basemul_avx512 LITHIUM_NAMESPACE(basemul_avx512)
#define pointwise_avx512 LITHIUM_NAMESPACE(pointwise_avx512)
#define nttunpack_avx512 LITHIUM_NAMESPACE(nttunpack_avx512)
#define nttunpack16_avx512 LITHIUM_NAMESPACE(nttunpack16_avx512)
#define poly_freeze_to_u16_avx512 \
    LITHIUM_NAMESPACE(poly_freeze_to_u16_avx512)
#define poly_freeze_small_to_u16_avx512 \
    LITHIUM_NAMESPACE(poly_freeze_small_to_u16_avx512)
#define poly_cast_modq_avx512 LITHIUM_NAMESPACE(poly_cast_modq_avx512)
#define ntt512_qdata LITHIUM_NAMESPACE(ntt512_qdata)
#define ntt512_scale LITHIUM_NAMESPACE(ntt512_scale)
#define ntt512_qinv LITHIUM_NAMESPACE(ntt512_qinv)
#define ntt512_r2 LITHIUM_NAMESPACE(ntt512_r2)
#define ntt512_zetas_fwd LITHIUM_NAMESPACE(ntt512_zetas_fwd)
#define ntt512_zetas_inv LITHIUM_NAMESPACE(ntt512_zetas_inv)
#define ntt512_basemul_gamma LITHIUM_NAMESPACE(ntt512_basemul_gamma)

/* nttunpack: canonical ref/ NTT order -> AVX-512 NTT slot layout.  The
 * canonical order is exactly the coefficient order produced by ref/
 * ntt16() after freezing inputs to [0,q); this is also the pk/sk wire
 * order for serialized bn = NTT_ref(b).  AVX-512 assembly uses a different
 * slot layout, so unpacking pk/sk bn or seed-to-NTT A0h must replay the
 * forward shuffle network with NO butterflies before basemul.  Equivalence
 * is tested as nttunpack_avx512(NTT_ref(a)) == NTT_avx512(a).  The
 * signature is mode-independent (n is baked into the linked ntt_avx512.S).
 */
void nttunpack_avx512(int16_t *dst, const int32_t *src);

/* nttunpack16: like nttunpack_avx512 but the source is ALREADY
 * canonical-order uint16 (in [0,q)), so it skips the int32->uint16 narrow.
 * Same slot-layout result; lets poly_ntt_import unpack a poly16 in place
 * (no int32 staging). */
void nttunpack16_avx512(int16_t *dst, const int16_t *src);

/* Standard-order NTT input conversion helpers.  They only convert the
 * coefficient representation; unlike nttunpack_avx512(), they do not
 * apply the AVX-512 NTT slot shuffle. */
void poly_freeze_to_u16_avx512(int16_t *dst, const int32_t *src);
void poly_freeze_small_to_u16_avx512(int16_t *dst, const int32_t *src);
void poly_cast_modq_avx512(int16_t *dst, const int32_t *src);

#if LITHIUM_MODE == 260
/* ---- complete NTT, n=512, q=64513 ---- */
extern const int16_t ntt512_qdata[64];
extern const int16_t ntt512_scale[64];
extern const int16_t ntt512_qinv[32];
extern const int16_t ntt512_r2[64];
extern const int16_t ntt512_zetas_fwd[144][32];
extern const int16_t ntt512_zetas_inv[144][32];
void ntt_avx512(int16_t *poly, const int16_t *qdata,
                const int16_t (*ztab)[32], const int16_t *scale);
void invntt_tomont_avx512(int16_t *poly, const int16_t *qdata,
                          const int16_t (*ztab)[32], const int16_t *scale);
void pointwise_avx512(int16_t *c, const int16_t *a, const int16_t *b,
                      const int16_t *qdata);

#else
/* ---- incomplete NTT, n=256 (120: q57601 strict, 180: q30977 self-reduce)
 */
extern const int16_t ntt512_qdata[64];
extern const int16_t ntt512_scale[64];
extern const int16_t ntt512_qinv[32];
extern const int16_t ntt512_r2[64];
extern const int16_t ntt512_zetas_fwd[56][32];
extern const int16_t ntt512_zetas_inv[56][32];
extern const int16_t ntt512_basemul_gamma[16][32];
void ntt_avx512(int16_t *poly, const int16_t *qdata,
                const int16_t (*ztab)[32], const int16_t *scale);
void invntt_tomont_avx512(int16_t *poly, const int16_t *qdata,
                          const int16_t (*ztab)[32], const int16_t *scale);

#    if LITHIUM_MODE == 120
/* q57601: basemul args (c, qdata, a, b, gamma); output strict [0,q). */
void basemul_avx512(int16_t *c, const int16_t *qdata, const int16_t *a,
                    const int16_t *b, const int16_t (*gamma)[32]);
#    else /* LITHIUM_MODE == 180 */
/* q30977: basemul args (c, a, b, qdata, gamma) -- DIFFERENT order.  The
 * AVX-512 basemul SELF-canonicalises to [0,q) (its montmul512/addmod512/
 * dmul512 macros each end in a branchless conditional subtract), so --
 * unlike the AVX2 path -- there is NO reduce_avx512 call after it. */
void basemul_avx512(int16_t *c, const int16_t *a, const int16_t *b,
                    const int16_t *qdata, const int16_t (*gamma)[32]);
#    endif

#endif /* LITHIUM_MODE */

#endif /* LITHIUM_NTT_AVX512_DECLS_H */
