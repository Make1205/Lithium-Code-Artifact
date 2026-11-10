/*
 * namespace.h - central symbol namespacing for IRS.
 *
 * Every IRS function is renamed to LITHIUM_NAMESPACE(x) =
 * lithium_nist<mode>_<ref|avx2>_x (see config.h) so the three parameter
 * sets
 * -- and the ref vs avx2 libraries -- can coexist in one link unit without
 * symbol clashes.  All the renames live here (instead of being scattered
 * through the per-module headers) so they can be managed and toggled in
 * one place.
 *
 * Toggle: define DISABLE_NAMESPACE=1 (e.g. -DDISABLE_NAMESPACE=1, or in
 * your IDE) to keep the plain, un-mangled function names -- convenient for
 * debugging, indexing, and single-library builds.  Leave it UNDEFINED for
 * the normal namespaced build; do NOT set its value in source.
 *
 * config.h includes this at its end, AFTER defining LITHIUM_NAMESPACE, so
 * the renames are visible everywhere config.h reaches (i.e. the whole
 * tree, via params.h).
 */
#ifndef LITHIUM_NAMESPACE_H
#define LITHIUM_NAMESPACE_H

#if !DISABLE_NAMESPACE

/* ---- top-level API (api.h, sign.h) ---- */
#    define crypto_sign_keypair LITHIUM_NAMESPACE(keypair)
#    define crypto_sign_signature LITHIUM_NAMESPACE(signature)
#    define crypto_sign_verify LITHIUM_NAMESPACE(verify)
#    define crypto_sign LITHIUM_NAMESPACE(sign)
#    define crypto_sign_open LITHIUM_NAMESPACE(open)

/* ---- reduce.h ---- */
#    define reduce32 LITHIUM_NAMESPACE(reduce32)
#    define caddq LITHIUM_NAMESPACE(caddq)
#    define freeze LITHIUM_NAMESPACE(freeze)
#    define caddq2 LITHIUM_NAMESPACE(caddq2)
#    define reduce_mod_2q LITHIUM_NAMESPACE(reduce_mod_2q)

/* ---- ntt.h ---- */
#    define ntt16 LITHIUM_NAMESPACE(ntt16)
#    define invntt16_tomont LITHIUM_NAMESPACE(invntt16_tomont)
#    define poly16_basemul_montgomery_native \
        LITHIUM_NAMESPACE(poly16_basemul_montgomery_native)

/* ---- rounding.h ---- */
#    define reduce2q LITHIUM_NAMESPACE(reduce2q)
#    define highbits LITHIUM_NAMESPACE(highbits)
#    define hbvalue LITHIUM_NAMESPACE(hbvalue)
#    define lsb LITHIUM_NAMESPACE(lsb)
#    define makehint LITHIUM_NAMESPACE(makehint)
#    define usehint LITHIUM_NAMESPACE(usehint)

/* ---- poly.h / polyvec.h ---- */
#    define poly_reduce LITHIUM_NAMESPACE(poly_reduce)
#    define poly_caddq LITHIUM_NAMESPACE(poly_caddq)
#    define poly_freeze LITHIUM_NAMESPACE(poly_freeze)
#    define poly_add LITHIUM_NAMESPACE(poly_add)
#    define poly_sub LITHIUM_NAMESPACE(poly_sub)
#    define poly_ntt LITHIUM_NAMESPACE(poly_ntt)
#    define poly_ntt_small LITHIUM_NAMESPACE(poly_ntt_small)
#    define poly_ntt_canonical LITHIUM_NAMESPACE(poly_ntt_canonical)
#    define poly_ntt_import LITHIUM_NAMESPACE(poly_ntt_import)
#    define unpack_pk_bn LITHIUM_NAMESPACE(unpack_pk_bn)
#    define poly_invntt_tomont LITHIUM_NAMESPACE(poly_invntt_tomont)
#    define poly_pointwise_montgomery \
        LITHIUM_NAMESPACE(poly_pointwise_montgomery)
#    define poly16_add LITHIUM_NAMESPACE(poly16_add)
#    define poly_mul LITHIUM_NAMESPACE(poly_mul)
#    define poly_uniform LITHIUM_NAMESPACE(poly_uniform)
#    define poly_uniform_4x LITHIUM_NAMESPACE(poly_uniform_4x)
#    define poly_invert LITHIUM_NAMESPACE(poly_invert)
#    define poly_sqnorm LITHIUM_NAMESPACE(poly_sqnorm)
#    define expand_A0_ntt LITHIUM_NAMESPACE(expand_A0_ntt)
#    define compute_b LITHIUM_NAMESPACE(compute_b)
#    define mat_mul_2q LITHIUM_NAMESPACE(mat_mul_2q)
#    define mat_mul_z1_2q LITHIUM_NAMESPACE(mat_mul_z1_2q)

/* ---- packing.h ---- */
#    define pack_pk LITHIUM_NAMESPACE(pack_pk)
#    define unpack_pk LITHIUM_NAMESPACE(unpack_pk)
#    define pack_sk LITHIUM_NAMESPACE(pack_sk)
#    define unpack_sk LITHIUM_NAMESPACE(unpack_sk)
#    define pack_sig_raw LITHIUM_NAMESPACE(pack_sig_raw)
#    define unpack_sig_raw LITHIUM_NAMESPACE(unpack_sig_raw)
#    define pack_sig LITHIUM_NAMESPACE(pack_sig)
#    define unpack_sig LITHIUM_NAMESPACE(unpack_sig)

/* ---- irs_rans.h ---- */
#    define rans_encode LITHIUM_NAMESPACE(rans_encode)
#    define rans_decode LITHIUM_NAMESPACE(rans_decode)

/* ---- sampler.c RCDT tables (external so the avx2 8-way sampler can share
 * them) ---- */
#    define RCDT_Z LITHIUM_NAMESPACE(RCDT_Z)
#    define RCDT_NOISE LITHIUM_NAMESPACE(RCDT_NOISE)

/* ---- irs.h ---- */
#    define poly_cyclic_shift LITHIUM_NAMESPACE(poly_cyclic_shift)
#    define compute_gamma_q52 LITHIUM_NAMESPACE(compute_gamma_q52)
#    define precompute_c LITHIUM_NAMESPACE(precompute_c)
#    define decide LITHIUM_NAMESPACE(decide)
#    define step LITHIUM_NAMESPACE(step)
#    define sign_irs LITHIUM_NAMESPACE(sign_irs)

#endif /* !DISABLE_NAMESPACE */
#endif /* LITHIUM_NAMESPACE_H */
