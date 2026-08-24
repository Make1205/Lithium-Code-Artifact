/*
 * params.h - Parameters for the Lithium signature scheme.
 *
 * Three parameter sets:
 *
 * Lithium-120: n=256, q=57601, ell=2, m=2, kappa=30,
 * r=182, incomplete NTT (base_deg=2)
 *
 * Lithium-180: n=256, q=30977, ell=3, m=3, kappa=58,
 * r=204, incomplete NTT (base_deg=2)
 *
 * Lithium-260: n=512, q=64513, ell=2, m=2, kappa=58,
 * r=257, full NTT (base_deg=1)
 *
 * Design-level constants are generated into @@AUTOGEN@@ regions inside the
 * hand-written headers (DO NOT hand-edit a region; regenerate via `make
 * tables`): bounds here (gen_bounds.py); NTT/Montgomery in ntt.h
 * (gen_zetas.py); exp in approx_exp.h (gen_approx_exp.py); RCDT in
 * sampler.c (gen_rcdt.py); rANS tables
 * + reserve + irs constants in irs_rans.h
 * (gen_rans_tables.py/SigSize.py/irs_ref.py).
 *
 * LITHIUM_MODE is selected via config.h.
 */

#ifndef LITHIUM_PARAMS_H
#define LITHIUM_PARAMS_H

#include <stdint.h>

#include "config.h"

/* ---- shared (hash / seed sizes) ---- */
#define SEEDBYTES 32
#define SKSEEDBYTES 64
#define CRHBYTES 64
#define TRBYTES 64
#define RNDBYTES 32

/* ---- mode-specific ---- */
#if LITHIUM_MODE == 120
#    define N 256
#    define Q 57601
#    define QBITS 16
#    define BASE_DEG 2
#    define ELL 2
#    define M 2
#    define KAPPA 30
#    define TAU 512
#    define R 182
#    define T 7
#    define CHASH_BYTES 32
#    define RCDT_ENTRIES 14
#    define W1_BITS 8 /* ceil(log2(2(q-1)/tau = 225)) */

#elif LITHIUM_MODE == 180
#    define N 256
#    define Q 30977
#    define QBITS 15
#    define BASE_DEG 2
#    define ELL 3
#    define M 3
#    define KAPPA 58
#    define TAU 512
#    define R 204
#    define T 7
#    define CHASH_BYTES 48
#    define RCDT_ENTRIES 16
#    define W1_BITS 7 /* ceil(log2(2(q-1)/tau = 121)) */

#elif LITHIUM_MODE == 260
#    define N 512
#    define Q 64513
#    define QBITS 16
#    define BASE_DEG 1
#    define ELL 2
#    define M 2
#    define KAPPA 58
#    define TAU 1024
#    define R 257
#    define T 8
#    define CHASH_BYTES 64
#    define RCDT_ENTRIES 20
#    define W1_BITS 7 /* ceil(log2(2(q-1)/tau = 126)) */

#else
#    error "Unsupported LITHIUM_MODE (expected 120, 180, or 260)"
#endif

/* ---- generated design-level constants ---- */
/* NTT/Montgomery constants now live in ntt.h (reduce.h pulls them in);
 * they are not re-exported here to avoid a params<->ntt include cycle. */
#include "irs_rans.h" /* RANS_RESERVED_BYTES (+ rANS tables, irs_consts) */

/* @@AUTOGEN:bounds@@ BEGIN */
/* EXACT integer squared norm thresholds (see tools/gen_bounds.py).
 *   KeyGen accept iff BK_LOW_SQ <= norm_sq <= BK_SQ (B_k' <= ||s|| < B_k)
 *   Sign   accept iff norm_sq <= BS_SQ;  Verify accept iff norm_sq <=
 * BV_SQ IRS-120: BK_LOW_SQ=1981 <= ||s||^2 <= BK_SQ=2116, BS_SQ=65448100,
 * BV_SQ=102137895 IRS-180: BK_LOW_SQ=2489 <= ||s||^2 <= BK_SQ=2632,
 * BS_SQ=108076816, BV_SQ=163678800 IRS-260: BK_LOW_SQ=4001 <= ||s||^2 <=
 * BK_SQ=4193, BS_SQ=231344100, BV_SQ=426258202
 */
#if LITHIUM_MODE == 120
/* B_k = 46.00893646172169  B_k' = 44.5007315301564  B_s = 8090  B_v =
 * 10106.329462945841 */
#    define BK_SQ INT64_C(2116)
#    define BK_LOW_SQ INT64_C(1981)
#    define BS_SQ INT64_C(65448100)
#    define BV_SQ INT64_C(102137895)
#endif

#if LITHIUM_MODE == 180
/* B_k = 51.3038533806998  B_k' = 49.889031744688346  B_s = 10396  B_v =
 * 12793.701585698811 */
#    define BK_SQ INT64_C(2632)
#    define BK_LOW_SQ INT64_C(2489)
#    define BS_SQ INT64_C(108076816)
#    define BV_SQ INT64_C(163678800)
#endif

#if LITHIUM_MODE == 260
/* B_k = 64.75410246578265  B_k' = 63.245897534217356  B_s = 15210  B_v =
 * 20646.021465119095 */
#    define BK_SQ INT64_C(4193)
#    define BK_LOW_SQ INT64_C(4001)
#    define BS_SQ INT64_C(231344100)
#    define BV_SQ INT64_C(426258202)
#endif
/* @@AUTOGEN:bounds@@ END */

/* ---- derived structural sizes ---- */
#define K (1 + ELL + M) /* z parts */
#define KSPLIT (1 << T) /* head/tail step */
#define CHALLENGE_PACKEDBYTES ((N + 7) / 8)

/* ---- mod 2q (rounding.c) ---- */
#define DQ (2 * Q)
#define HINT_MOD (2 * (Q - 1))
#define W1_MAX                                   \
    (HINT_MOD / TAU) /* tau-buckets: 225/121/126 \
                      */

/* ---- Gaussian sampler ---- */
#define GAUSS_K 128
#define K_BITS 7
#define Y_BITS 7
#define TWO_K_BITS 8
#define RCDT_LIMB_BITS 31
#define RCDT_BITS 93
#define NOISE_RCDT_ENTRIES 10
#define GAUSS_XMAX (RCDT_ENTRIES - 1)
#define NOISE_BITS 5

/* ---- packing sizes ---- */
#define POLYPK_PACKEDBYTES ((N * QBITS + 7) / 8)
#define POLYNOISE_PACKEDBYTES ((N * NOISE_BITS + 7) / 8)
#define POLYZ_LO_PACKEDBYTES ((N * T + 7) / 8)
#define RANS_BLOCK_BYTES (2 + RANS_RESERVED_BYTES)

/* ---- key + signature sizes ---- */
#define CRYPTO_PUBLICKEYBYTES (SEEDBYTES + M * POLYPK_PACKEDBYTES)
#define CRYPTO_SECRETKEYBYTES                                        \
    (SEEDBYTES + TRBYTES + SKSEEDBYTES + K * POLYNOISE_PACKEDBYTES + \
     M * POLYPK_PACKEDBYTES)
#define CRYPTO_BYTES                            \
    (CHALLENGE_PACKEDBYTES + RANS_BLOCK_BYTES + \
     (ELL + 1) * POLYZ_LO_PACKEDBYTES)

#endif /* LITHIUM_PARAMS_H */
