/*
 * ntt_avx_decls.h - extern declarations for the in-tree AVX2 NTT assembly
 * and its constant tables, switched by LITHIUM_MODE.
 *
 * The assembly + constants are pulled (by explicit path in the Makefile)
 * from tools/ntt/q57601n256/  (mode 120, incomplete base_deg=2, strict)
 *   tools/ntt/q30977n256/  (mode 180, incomplete base_deg=2, LAZY ->
 * reduce_avx) tools/ntt/q64513n512/  (mode 260, complete base_deg=1,
 * strict) The generated assembly and constant tables are namespaced with
 * LITHIUM_NAMESPACE so multiple Lithium modes/implementations can coexist
 * in one link unit.  These prototypes are copied from the matching
 * tools/ntt/<qset>/ntt_ref.h so the compiler type-checks the (mode-
 * specific!) basemul argument order -- see the AVX2 polyvec routing.
 */
#ifndef LITHIUM_NTT_AVX_DECLS_H
#define LITHIUM_NTT_AVX_DECLS_H

#include <stdint.h>

#include "config.h"

#define ntt_avx LITHIUM_NAMESPACE(ntt_avx)
#define invntt_tomont_avx LITHIUM_NAMESPACE(invntt_tomont_avx)
#define basemul_avx LITHIUM_NAMESPACE(basemul_avx)
#define pointwise_avx LITHIUM_NAMESPACE(pointwise_avx)
#define reduce_avx LITHIUM_NAMESPACE(reduce_avx)
#define nttunpack_avx LITHIUM_NAMESPACE(nttunpack_avx)
#define nttunpack16_avx LITHIUM_NAMESPACE(nttunpack16_avx)
#define poly_freeze_to_u16_avx LITHIUM_NAMESPACE(poly_freeze_to_u16_avx)
#define poly_freeze_small_to_u16_avx \
    LITHIUM_NAMESPACE(poly_freeze_small_to_u16_avx)
#define poly_cast_modq_avx LITHIUM_NAMESPACE(poly_cast_modq_avx)
#define ntt_qdata LITHIUM_NAMESPACE(ntt_qdata)
#define ntt_ninv LITHIUM_NAMESPACE(ntt_ninv)
#define ntt_cross_fwd LITHIUM_NAMESPACE(ntt_cross_fwd)
#define ntt_cross_inv LITHIUM_NAMESPACE(ntt_cross_inv)
#define ntt_zetas_fwd LITHIUM_NAMESPACE(ntt_zetas_fwd)
#define ntt_zetas_inv LITHIUM_NAMESPACE(ntt_zetas_inv)
#define ntt_z0 LITHIUM_NAMESPACE(ntt_z0)
#define ntt_z0inv LITHIUM_NAMESPACE(ntt_z0inv)
#define ntt_scale LITHIUM_NAMESPACE(ntt_scale)
#define ntt_basemul_gamma LITHIUM_NAMESPACE(ntt_basemul_gamma)

/* nttunpack: canonical ref/ NTT order -> AVX2 NTT slot layout.  The
 * canonical order is exactly the coefficient order produced by ref/
 * ntt16() after freezing inputs to [0,q); this is also the pk/sk wire
 * order for serialized bn = NTT_ref(b).  AVX2 assembly uses a different
 * slot layout, so unpacking pk/sk bn or seed-to-NTT A0h must replay the
 * forward shuffle network with NO butterflies before basemul.  Equivalence
 * is tested as nttunpack_avx(NTT_ref(a)) == NTT_avx(a).  The signature is
 * mode-independent (n is baked into the linked ntt.S). */
void nttunpack_avx(int16_t *dst, const int32_t *src);

/* nttunpack16: like nttunpack_avx but the source is ALREADY
 * canonical-order uint16 (in [0,q)), so it skips the int32->uint16 narrow.
 * Same slot-layout result; lets poly_ntt_import unpack a poly16 in place
 * (no int32 staging). */
void nttunpack16_avx(int16_t *dst, const int16_t *src);

/* Standard-order NTT input conversion helpers.  They only convert the
 * coefficient representation; unlike nttunpack_avx(), they do not apply
 * the AVX NTT slot shuffle. */
void poly_freeze_to_u16_avx(int16_t *dst, const int32_t *src);
void poly_freeze_small_to_u16_avx(int16_t *dst, const int32_t *src);
void poly_cast_modq_avx(int16_t *dst, const int32_t *src);

#if LITHIUM_MODE == 260
/* ---- complete NTT, n=512, q=64513 ---- */
extern const int16_t ntt_qdata[48];
extern const int16_t ntt_ninv[32];
extern const int16_t ntt_cross_fwd[8][16];
extern const int16_t ntt_cross_inv[8][16];
extern const int16_t ntt_zetas_fwd[224][16];
extern const int16_t ntt_zetas_inv[224][16];
void ntt_avx(int16_t *poly, const int16_t *qdata,
             const int16_t (*ztab)[16], const int16_t (*cross)[16],
             const int16_t *ninv);
void invntt_tomont_avx(int16_t *poly, const int16_t *qdata,
                       const int16_t (*ztab)[16],
                       const int16_t (*cross)[16], const int16_t *ninv);
void pointwise_avx(int16_t *c, const int16_t *a, const int16_t *b,
                   const int16_t *qdata);

#else
/* ---- incomplete NTT, n=256 (120: q57601 strict, 180: q30977 lazy) ----
 */
extern const int16_t ntt_qdata[48];
extern const int16_t ntt_z0[32], ntt_z0inv[32], ntt_scale[32];
extern const int16_t ntt_zetas_fwd[96][16];
extern const int16_t ntt_zetas_inv[96][16];
extern const int16_t ntt_basemul_gamma[32][16];
void ntt_avx(int16_t *poly, const int16_t *qdata,
             const int16_t (*ztab)[16], const int16_t *z0,
             const int16_t *scale);
void invntt_tomont_avx(int16_t *poly, const int16_t *qdata,
                       const int16_t (*ztab)[16], const int16_t *z0inv,
                       const int16_t *scale);

#    if LITHIUM_MODE == 120
/* q57601: basemul args (c, qdata, a, b, gamma) */
void basemul_avx(int16_t *c, const int16_t *qdata, const int16_t *a,
                 const int16_t *b, const int16_t (*gamma)[16]);
#    else /* LITHIUM_MODE == 180 */
/* q30977: basemul args (c, a, b, qdata, gamma) -- DIFFERENT order -- plus
 * a lazy-reduction canonicaliser reduce_avx([0,2q) -> [0,q)). */
void basemul_avx(int16_t *c, const int16_t *a, const int16_t *b,
                 const int16_t *qdata, const int16_t (*gamma)[16]);
void reduce_avx(int16_t *poly);
#    endif

#endif /* LITHIUM_MODE */

#endif /* LITHIUM_NTT_AVX_DECLS_H */
