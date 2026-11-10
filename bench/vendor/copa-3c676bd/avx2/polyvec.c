/*
 * polyvec.c - polynomial arithmetic + NTWE matrix/vector operations for
 * IRS.
 *
 * This file contains ring arithmetic, the signed<->uint16 NTT boundary,
 * ring inversion, squared norms, and the NTWE matrix/vector operations.
 *
 * Boundary between the signed scheme domain (poly, int32) and the uint16
 * NTT / mod-q domain (poly16, [0,q)): poly_ntt freezes + forward-NTTs into
 * a poly16; poly_invntt_tomont inverse-NTTs a poly16 back to a signed poly
 * in [0,q).  The Montgomery bookkeeping (R = 2^16) is documented in ntt.c
 * / reduce.h.
 *
 * All ring products go through the uint16 NTT.  The mod-2q lift is exactly
 * the structure of A = [qj - 2b, 2 A0, 2 I_m]; see polyvec.h.  Montgomery
 * convention (validated by t_poly): pointwise_montgomery (owes one R)
 * followed by invntt_tomont (R^2 scale, pays it back) gives a*b in [0,q)
 * with no leftover R.
 */
#include "polyvec.h"

#include <stdint.h>
#include <string.h>

#include "config.h" /* namespace for the merged AVX2 NTT symbols */
#include "fips202.h"
#include "fips202x4.h" /* shake128x4 for the merged 4-way uniform expansion */
#include "ntt.h"
#include "ntt_avx_decls.h" /* nttunpack_avx: standard -> avx NTT layout ladder */
#include "prof.h"
#include "reduce.h"
#ifdef __AVX2__
#    include <immintrin.h> /* AVX2 poly_sqnorm fast path (scalar #else is the
                          * portable reference; exact integer -> KAT-neutral) */
#endif

/* ---- small modular helpers (q is a 16-bit prime), values in [0,q) ----
 */
#if Q >= 65536
#    error "polyvec mod-q helpers assume a 16-bit modulus"
#endif

/* Reduce a product of two mod-q residues.  Since q < 2^16, every caller
 * has a < q^2 < 2^32.  REC=floor(2^32/q) makes the reciprocal quotient
 * undershoot by at most one, so one masked subtract gives the exact
 * residue without a runtime division or C modulus operation. */
#define MODQ_RECIP32 (((uint64_t)1 << 32) / (uint64_t)Q)
static inline uint16_t modq_reduce_product(uint32_t a)
{
    uint32_t qh = (uint32_t)(((uint64_t)a * MODQ_RECIP32) >> 32);
    uint32_t r = a - qh * (uint32_t)Q; /* r < 2q */
    uint32_t ge = (uint32_t)0 - (uint32_t)(r >= (uint32_t)Q);
    return (uint16_t)(r - (ge & (uint32_t)Q));
}

static inline uint16_t modq_mul_u16(uint16_t a, uint16_t b)
{
    return modq_reduce_product((uint32_t)a * b);
}

static inline uint16_t modq_sqr_u16(uint16_t a)
{
    return modq_reduce_product((uint32_t)a * a);
}

static inline uint16_t modq_sub_u16(uint16_t a, uint16_t b)
{
    uint32_t r = (uint32_t)a + (uint32_t)Q - (uint32_t)b;
    uint32_t ge = (uint32_t)0 - (uint32_t)(r >= (uint32_t)Q);
    return (uint16_t)(r - (ge & (uint32_t)Q));
}

static inline uint16_t small_to_u16_modq(int32_t x)
{
    /* Precondition: -(Q - 1) <= x <= Q - 1.  The high bit
     * of the modulo-2^32 unsigned view selects +q for negative inputs;
     * unsigned arithmetic makes the mask explicit and avoids
     * implementation-defined signed shifts. */
    uint32_t ux = (uint32_t)x;
    uint32_t neg = ux >> 31;
    uint32_t mask = 0u - neg;
    return (uint16_t)(ux + (mask & (uint32_t)Q));
}

static uint16_t modq_pow(uint16_t base, uint32_t e)
{
    uint16_t b = base;
    uint16_t r = 1;
    while (e) {
        if (e & 1u)
            r = modq_mul_u16(r, b);
        b = modq_sqr_u16(b);
        e >>= 1;
    }
    return r;
}
static uint16_t modq_inv(uint16_t a)
{
    return modq_pow(a, (uint32_t)(Q - 2));
}
/* Batch modular inverse (Montgomery's trick): inv[i] = v[i]^{-1} mod q for
 * all i, with ONE Fermat modexp + ~3n modmuls instead of n modexps -- the
 * dominant cost of poly_invert.  v[i] in [0,q); returns -1 if any v[i]==0
 * (non-invertible).  inv and v MUST NOT alias.  Each result equals the
 * per-element modq_inv (modular inverse is unique) -> bit-identical / KAT-
 * neutral.  No data-dependent branch on the value (the only early-out is
 * the structural non-invertibility reject, as before). */
static int batch_modq_inv(uint16_t *inv, const uint16_t *v, int n)
{
    uint16_t acc = 1;
    int zero = 0;
    for (int i = 0; i < n; i++) {
        inv[i] = (uint16_t)acc; /* prefix product v[0..i-1] */
        zero |= (v[i] == 0);
        acc = modq_mul_u16(acc, v[i]);
    }
    if (zero)
        return -1; /* some element non-invertible -> resample in the caller
                    */
    acc = modq_inv(acc); /* (prod v)^{-1}, the single Fermat */
    for (int i = n - 1; i >= 0; i--) {
        inv[i] = modq_mul_u16(acc, inv[i]); /* = v[i]^{-1} */
        acc = modq_mul_u16(acc, v[i]);
    }
    return 0;
}

void poly_reduce(poly *a)
{
    for (int i = 0; i < N; i++)
        a->coeffs[i] = reduce32(a->coeffs[i]);
}
void poly_caddq(poly *a)
{
    for (int i = 0; i < N; i++)
        a->coeffs[i] = caddq(a->coeffs[i]);
}
void poly_freeze(poly *a)
{
    for (int i = 0; i < N; i++)
        a->coeffs[i] = freeze(a->coeffs[i]);
}

void poly_add(poly *r, const poly *a, const poly *b)
{
    for (int i = 0; i < N; i++)
        r->coeffs[i] = a->coeffs[i] + b->coeffs[i];
}
void poly_sub(poly *r, const poly *a, const poly *b)
{
    for (int i = 0; i < N; i++)
        r->coeffs[i] = a->coeffs[i] - b->coeffs[i];
}

/* signed poly -> [0,q) uint16 -> forward NTT (normal-domain in/out).
 * In the AVX2 build (USE_AVX2_NTT) the three NTT-domain primitives
 * below are provided below in this file, where they route to
 * the in-tree AVX2 NTT assembly.  poly_invert keeps its OWN scalar
 * transforms (see below) so it never depends on the AVX shuffle order. */
#ifndef USE_AVX2_NTT
void poly_ntt(poly16 *out, const poly *in)
{
    for (int i = 0; i < N; i++)
        out->coeffs[i] = (uint16_t)freeze(in->coeffs[i]);
    ntt16(out->coeffs);
}

void poly_ntt_small(poly16 *out, const poly *in)
{
    for (int i = 0; i < N; i++)
        out->coeffs[i] = small_to_u16_modq(in->coeffs[i]);
    ntt16(out->coeffs);
}

/* inverse NTT (+ R^2 scale) of a NTT-domain poly16 -> signed poly in
 * [0,q). */
void poly_invntt_tomont(poly *out, const poly16 *in)
{
    uint16_t t[N];
    memcpy(t, in->coeffs, sizeof t);
    invntt16_tomont(t);
    for (int i = 0; i < N; i++)
        out->coeffs[i] = (int32_t)t[i];
}
#endif /* !USE_AVX2_NTT */

static void poly_ntt_cast(poly16 *out, const poly *in)
{
#ifdef USE_AVX2_NTT
    /* in is lifecycle-proven canonical [0,q); only narrow to uint16 before
     * the native forward NTT. */
    poly_cast_modq_avx((int16_t *)out->coeffs, in->coeffs);
#    if LITHIUM_MODE == 260
    ntt_avx((int16_t *)out->coeffs, ntt_qdata, ntt_zetas_fwd,
            ntt_cross_fwd, ntt_ninv);
#    else
    ntt_avx((int16_t *)out->coeffs, ntt_qdata, ntt_zetas_fwd, ntt_z0,
            ntt_scale);
#    endif
#else
    for (int i = 0; i < N; i++)
        out->coeffs[i] = (uint16_t)in->coeffs[i];
    ntt16(out->coeffs);
#endif
}

void poly16_add(poly16 *r, const poly16 *a, const poly16 *b)
{
    for (int i = 0; i < N; i++)
        r->coeffs[i] = addm16(a->coeffs[i], b->coeffs[i]);
}

#ifndef USE_AVX2_NTT
void poly_pointwise_montgomery(poly16 *r, const poly16 *a, const poly16 *b)
{
#    if BASE_DEG == 1
    for (int i = 0; i < N; i++)
        r->coeffs[i] = fqmul16(a->coeffs[i], b->coeffs[i]);
#    else
    poly16_basemul_montgomery_native(r->coeffs, a->coeffs, b->coeffs);
#    endif
}
#endif /* !USE_AVX2_NTT */

void poly_mul(poly *r, const poly *a, const poly *b)
{
    poly16 ah, bh, rh;
    poly_ntt(&ah, a);
    poly_ntt(&bh, b);
    poly_pointwise_montgomery(&rh, &ah, &bh);
    poly_invntt_tomont(r, &rh); /* == a*b in [0,q), no leftover R */
}

void poly_uniform(poly *a, const uint8_t seed[SEEDBYTES], uint16_t nonce)
{
    uint8_t ext[SEEDBYTES + 2];
    uint8_t blk[SHAKE128_RATE];
    keccak_state st;
    const uint32_t mask = (uint32_t)((1u << QBITS) - 1u);
    unsigned ctr = 0;

    memcpy(ext, seed, SEEDBYTES);
    ext[SEEDBYTES] = (uint8_t)(nonce & 0xff);
    ext[SEEDBYTES + 1] = (uint8_t)(nonce >> 8);
    shake128_init(&st);
    shake128_absorb(&st, ext, sizeof ext);
    shake128_finalize(&st);

    while (ctr < (unsigned)N) {
        shake128_squeezeblocks(blk, 1, &st);
        for (unsigned i = 0; i + 2 <= SHAKE128_RATE && ctr < (unsigned)N;
             i += 2) {
            uint32_t t = (uint32_t)blk[i] | ((uint32_t)blk[i + 1] << 8);
            t &= mask;
            if (t < (uint32_t)Q)
                a->coeffs[ctr++] = (int32_t)t;
        }
    }
}

/* 4-way uniform expansion (canonical KAT layout for AVX2): 4 independent
 * SHAKE-128 streams (seed, nonce_j).  Sequential here; the AVX2 build
 * replaces it with the shake128x4 version below, which is
 * byte-identical per lane.  All four lanes are real (expand_A0 handles
 * any remainder with scalar poly_uniform). */
#ifndef USE_AVX2_SHAKE4X
void poly_uniform_4x(poly *a0, poly *a1, poly *a2, poly *a3,
                     const uint8_t seed[SEEDBYTES], uint16_t nonce0,
                     uint16_t nonce1, uint16_t nonce2, uint16_t nonce3)
{
    poly_uniform(a0, seed, nonce0);
    poly_uniform(a1, seed, nonce1);
    poly_uniform(a2, seed, nonce2);
    poly_uniform(a3, seed, nonce3);
}

/* scalar A0 chunk fallback (the avx2 build uses 2x keccak4x in
 * the shake128x4 path below instead). */
void poly_uniform_chunkx8(int32_t *flat, int per_lane,
                          const uint8_t seed[SEEDBYTES])
{
    const uint32_t mask = (uint32_t)((1u << QBITS) - 1u);
    for (int j = 0; j < 8; j++) {
        uint8_t ext[SEEDBYTES + 2];
        uint8_t blk[SHAKE128_RATE];
        keccak_state st;
        memcpy(ext, seed, SEEDBYTES);
        ext[SEEDBYTES] = (uint8_t)j;
        ext[SEEDBYTES + 1] = 0;
        shake128_init(&st);
        shake128_absorb(&st, ext, sizeof ext);
        shake128_finalize(&st);
        int ctr = 0;
        int32_t *o = flat + (size_t)j * per_lane;
        while (ctr < per_lane) {
            shake128_squeezeblocks(blk, 1, &st);
            for (unsigned i = 0; i + 2 <= SHAKE128_RATE && ctr < per_lane;
                 i += 2) {
                uint32_t t =
                    (uint32_t)blk[i] | ((uint32_t)blk[i + 1] << 8);
                t &= mask;
                if (t < (uint32_t)Q)
                    o[ctr++] = (int32_t)t;
            }
        }
    }
}
#endif /* !USE_AVX2_SHAKE4X */

int poly_invert(poly *r, const poly *a)
{
    poly16 th, rh;
    /* Scalar forward NTT (inlined): residue inversion below indexes
     * coefficients in the scalar Cooley-Tukey order (zetas[half+i],
     * adjacent quadratic pairs), so poly_invert MUST NOT use the AVX2
     * shuffle-ordered transform.  This is byte-identical to the ref
     * poly_ntt body. */
    for (int i = 0; i < N; i++)
        th.coeffs[i] = (uint16_t)freeze(a->coeffs[i]);
    ntt16(th.coeffs); /* exact residues of a in [0,q) */
#if BASE_DEG == 1
    /* full NTT: invert every residue, batched (one Fermat for all
     * N).
     */
    if (batch_modq_inv(rh.coeffs, th.coeffs, N) != 0)
        return -1;
#else
    /* incomplete NTT: invert (a0 + a1 X) in F_q[X]/(X^2 - rr) per
     * quadratic pair.  Pair p (p=0..N/2-1) is at coeffs[2p,2p+1]; p even
     * uses root zetas[half + p/2], p odd uses -root.  D = a0^2 - rr*a1^2
     * is the determinant; batch-invert all D (one Fermat), then finish
     * each pair out0 = a0/D, out1 = -a1/D using the same
     * two-coefficient inverse formula. */
    {
        const unsigned half = NTT_ZETAS_SIZE / 2;
        const unsigned npair = (unsigned)N / 2;
        uint16_t Dv[N / 2], Di[N / 2];
        for (unsigned p = 0; p < npair; p++) {
            uint16_t a0 = th.coeffs[2 * p], a1 = th.coeffs[2 * p + 1];
            uint16_t root = montgomery_reduce16(zetas[half + p / 2]);
            uint16_t rr = (p & 1u) ? modq_sub_u16(0, root) : root;
            uint16_t a0sq = modq_sqr_u16(a0);
            uint16_t rr_a1sq = modq_mul_u16(rr, modq_sqr_u16(a1));
            Dv[p] = modq_sub_u16(a0sq, rr_a1sq);
        }
        if (batch_modq_inv(Di, Dv, (int)npair) != 0)
            return -1;
        for (unsigned p = 0; p < npair; p++) {
            uint16_t a0 = th.coeffs[2 * p], a1 = th.coeffs[2 * p + 1];
            uint16_t di = Di[p];
            rh.coeffs[2 * p] = modq_mul_u16(a0, di);
            rh.coeffs[2 * p + 1] = modq_mul_u16(modq_sub_u16(0, a1), di);
        }
    }
#endif
    /* exact residues -> coeff domain: invntt_tomont({Y}) =
     * invNTT_plain({Y*R}), so feed residue * R^{-1} =
     * montgomery_reduce16(residue). */
    for (int i = 0; i < N; i++)
        rh.coeffs[i] = montgomery_reduce16(rh.coeffs[i]);
    /* Scalar inverse NTT (inlined), byte-identical to the ref
     * poly_invntt_tomont body; see the forward note above for why
     * poly_invert stays scalar. */
    {
        uint16_t t[N];
        memcpy(t, rh.coeffs, sizeof t);
        invntt16_tomont(t);
        for (int i = 0; i < N; i++)
            r->coeffs[i] = (int32_t)t[i];
    }
    poly_freeze(r);
    return 0;
}

int64_t poly_sqnorm(const poly *a)
{
#ifdef __AVX2__
    /* Norm inputs are centered with |coeff| < q (z/s/z1 are far below q),
     * so reduce32 is the identity and freeze+center is just the masked
     * centering below -- no Barrett needed.  Exact int64 sum (vpmuldq) ->
     * KAT-neutral. */
    const __m256i qvec = _mm256_set1_epi32(Q);
    const __m256i QH = _mm256_set1_epi32(Q / 2);
    __m256i acc = _mm256_setzero_si256();
    for (int i = 0; i < N; i += 8) {
        __m256i x = _mm256_loadu_si256((const __m256i *)(a->coeffs + i));
        __m256i f =
            _mm256_add_epi32(x, _mm256_and_si256(_mm256_srai_epi32(x, 31),
                                                 qvec)); /* x<0 -> +q */
        __m256i c =
            _mm256_sub_epi32(f, _mm256_and_si256(_mm256_cmpgt_epi32(f, QH),
                                                 qvec)); /* center */
        acc = _mm256_add_epi64(acc, _mm256_mul_epi32(c, c));
        acc = _mm256_add_epi64(acc,
                               _mm256_mul_epi32(_mm256_srli_epi64(c, 32),
                                                _mm256_srli_epi64(c, 32)));
    }
    int64_t t[4];
    _mm256_storeu_si256((__m256i *)t, acc);
    return t[0] + t[1] + t[2] + t[3];
#else
    int64_t s = 0;
    for (int i = 0; i < N; i++) {
        int32_t c = freeze(a->coeffs[i]);
        int32_t gt = ((Q / 2) - c) >> 31;
        c -= gt & Q;
        s += (int64_t)c * c;
    }
    return s;
#endif
}

/* ============================================================
 * NTWE matrix/vector operations
 * ============================================================ */
void expand_A0_ntt(poly16 A0h[M][ELL], const uint8_t rho[SEEDBYTES])
{
    /* Seed-to-NTT (Pre-NTT): sample the whole A0 (M*ELL*N uniform [0,q)
     * coeffs) as 8 equal chunks from SHAKE128(rho||j), j=0..7 -- byte-
     * identical per lane to the ref scalar, full keccak4x utilisation for
     * every mode (no idle lanes).  Every q < 2^16, so the standard-order
     * samples are written straight into (uint16_t*)A0h (poly p = i*ELL+j
     * -> the slice [p*N ..)).  poly_ntt_import then lays each poly, IN
     * PLACE, into the avx ntt.S slot order (transpose-equivalent root, so
     * the avx A0 matches ref's at the scheme level -> shared KAT). */
    _Static_assert(sizeof(poly16) == 2 * (size_t)N,
                   "(uint16_t*)A0h cast needs poly16 == N uint16_t");
    poly_uniform_chunkx8((uint16_t *)A0h, (M * ELL * N) / 8, rho);

    for (int i = 0; i < M; i++)
        for (int j = 0; j < ELL; j++)
            poly_ntt_import(&A0h[i][j], &A0h[i][j]);
}

int compute_b(poly b[M], const poly16 A0h[M][ELL], const poly s0[ELL],
              const poly e[M], const poly *f)
{
    poly finv;
    if (poly_invert(&finv, f) != 0)
        return -1; /* f not invertible */

    poly16 s0h[ELL], finv_ntt, acc_h, prod, eh;
    poly acc_n;
    poly_ntt_cast(&finv_ntt, &finv);
    for (int j = 0; j < ELL; j++)
        poly_ntt_small(&s0h[j], &s0[j]);
    /* A0h is already NTT-domain (seed-to-NTT): no forward transform.
     * s0 is key noise, hence small signed; use the branchless narrow NTT
     * boundary instead of the generic freeze/Barrett path. */

    for (int i = 0; i < M; i++) {
        /* acc_h = sum_j A0[i][j] * s0[j]  (NTT-domain accumulate, then
         * invNTT) */
        poly_pointwise_montgomery(&acc_h, &A0h[i][0], &s0h[0]);
        for (int j = 1; j < ELL; j++) {
            poly_pointwise_montgomery(&prod, &A0h[i][j], &s0h[j]);
            poly16_add(&acc_h, &acc_h, &prod);
        }
        poly_invntt_tomont(&acc_n,
                           &acc_h); /* == A0 s0 (no leftover R), [0,q) */
        poly_add(&acc_n, &acc_n, &e[i]); /* (A0 s0 + e), normal domain */
        poly_reduce(&acc_n);
        /* b_i = (A0 s0 + e) * f^{-1} */
        poly_ntt_small(&eh, &acc_n);
        poly_pointwise_montgomery(&acc_h, &eh, &finv_ntt);
        poly_invntt_tomont(&b[i], &acc_h); /* == (A0 s0 + e) f^{-1} */
        poly_freeze(&b[i]);                /* b_i in [0,q) */
    }
    return 0;
}

/* t_i = -b_i x0 + (A0 xs)_i (mod q), normal domain in [0,q).  Inputs
 * NTT-domain. */
static void compute_t(poly *t, const poly16 A0h[M][ELL],
                      const poly16 *bh_i, const poly16 *x0h,
                      const poly16 xsh[ELL], int i)
{
    poly16 acc, prod;
    PROF_START(tp);
    poly_pointwise_montgomery(&acc, bh_i, x0h);
    PROF_STOP(PT_NTT_PW, tp);
    for (int k = 0; k < N; k++)
        acc.coeffs[k] = subm16(0, acc.coeffs[k]); /* -b_i x0 */
    for (int j = 0; j < ELL; j++) {
        PROF_START(tp2);
        poly_pointwise_montgomery(&prod, &A0h[i][j], &xsh[j]);
        PROF_STOP(PT_NTT_PW, tp2);
        poly16_add(&acc, &acc, &prod);
    }
    PROF_START(ti);
    poly_invntt_tomont(t, &acc);
    PROF_STOP(PT_NTT_INV, ti); /* -b_i x0 + (A0 xs)_i */
}

void mat_mul_2q(poly w[M], const poly16 A0h[M][ELL], const poly16 bh[M],
                const poly x[K])
{
    /* A0h / bh are already NTT-domain (Pre-NTT): only the fresh vector x
     * is forward-transformed here.  x is the Gaussian y vector; the
     * sampler bound 128*RCDT_ENTRIES+127 is < q in every mode. */
    poly16 x0h, xsh[ELL];
    PROF_START(tfwd);
    poly_ntt_small(&x0h, &x[0]);
    for (int j = 0; j < ELL; j++)
        poly_ntt_small(&xsh[j], &x[1 + j]);
    PROF_STOP(PT_NTT_FWD, tfwd);

    for (int i = 0; i < M; i++) {
        poly t;
        compute_t(&t, A0h, &bh[i], &x0h, xsh, i);
        poly_add(&t, &t, &x[1 + ELL + i]); /* + xe_i */
        poly_freeze(&t);                   /* t in [0,q) */
        /* v = 2t + q*(x0 mod 2); w = reduce_mod_2q(v).  After freeze
         * t in [0,q) so v in [0,3q) (always >= 0), and reduce_mod_2q
         * collapses to a single conditional subtract of 2q -> branchless
         * SIMD (bit-identical to the scalar reduce_mod_2q over this
         * range). The parity uses the TRUE coefficient parity (raw x0 &
         * 1), not freeze(x0)&1 (freeze flips the parity of negatives). */
#ifdef __AVX2__
        {
            const __m256i one = _mm256_set1_epi32(1);
            const __m256i vq = _mm256_set1_epi32(Q);
            const __m256i v2q = _mm256_set1_epi32(DQ);
            const __m256i v2qm1 = _mm256_set1_epi32(DQ - 1);
            const __m256i zero = _mm256_setzero_si256();
            for (int k = 0; k < N; k += 8) {
                __m256i v = _mm256_slli_epi32(
                    _mm256_loadu_si256((const __m256i *)&t.coeffs[k]), 1);
                if (i == 0) {
                    __m256i x0 = _mm256_loadu_si256(
                        (const __m256i *)&x[0].coeffs[k]);
                    __m256i m = _mm256_sub_epi32(
                        zero, _mm256_and_si256(x0, one)); /* 0 or -1 */
                    v = _mm256_add_epi32(v, _mm256_and_si256(m, vq));
                }
                __m256i ge = _mm256_cmpgt_epi32(v, v2qm1); /* v >= 2q */
                v = _mm256_sub_epi32(v, _mm256_and_si256(ge, v2q));
                _mm256_storeu_si256((__m256i *)&w[i].coeffs[k], v);
            }
        }
#else
        for (int k = 0; k < N; k++) {
            int32_t v = 2 * t.coeffs[k];
            if (i == 0)
                v += (int32_t)((int64_t)Q * (x[0].coeffs[k] & 1));
            w[i].coeffs[k] = reduce_mod_2q(v);
        }
#endif
    }
}

void mat_mul_z1_2q(poly wt[M], const poly16 A0h[M][ELL],
                   const poly16 bh[M], const poly z1[ELL + 1],
                   const poly *c)
{
    /* A0h / bh are already NTT-domain (Pre-NTT): only z1 is transformed.
     * Verifier z1 is pre-bounded by rANS z-head range plus t-bit low part
     * before the later norm check, and the resulting maxima are < q. */
    poly16 z0h, zsh[ELL];
    poly_ntt_small(&z0h, &z1[0]);
    for (int j = 0; j < ELL; j++)
        poly_ntt_small(&zsh[j], &z1[1 + j]);

    for (int i = 0; i < M; i++) {
        poly t;
        compute_t(&t, A0h, &bh[i], &z0h, zsh,
                  i); /* t = -b_i z0 + (A0 zs)_i, already [0,q) */
        /* v = 2t + q*(z0 mod 2) - q*(c mod 2); w = reduce_mod_2q(v).
         * compute_t ends in poly_invntt_tomont so t is already [0,q) (no
         * intervening poly_add, unlike mat_mul_2q), so the freeze was a
         * no-op -- dropped.  v in [-q,3q): reduce_mod_2q collapses to a
         * single +2q (if v<0) and a single -2q (if v>=2q), branchless SIMD
         * (bit-identical to scalar reduce_mod_2q over this range).  Raw
         * coefficient parities (not freeze()&1). */
#ifdef __AVX2__
        {
            const __m256i one = _mm256_set1_epi32(1);
            const __m256i vq = _mm256_set1_epi32(Q);
            const __m256i v2q = _mm256_set1_epi32(DQ);
            const __m256i v2qm1 = _mm256_set1_epi32(DQ - 1);
            const __m256i zero = _mm256_setzero_si256();
            for (int k = 0; k < N; k += 8) {
                __m256i v = _mm256_slli_epi32(
                    _mm256_loadu_si256((const __m256i *)&t.coeffs[k]), 1);
                if (i == 0) {
                    __m256i z0 = _mm256_loadu_si256(
                        (const __m256i *)&z1[0].coeffs[k]);
                    __m256i cc =
                        _mm256_loadu_si256((const __m256i *)&c->coeffs[k]);
                    __m256i zm = _mm256_sub_epi32(
                        zero, _mm256_and_si256(z0, one)); /* 0 or -1 */
                    __m256i cm =
                        _mm256_sub_epi32(zero, _mm256_and_si256(cc, one));
                    v = _mm256_add_epi32(v, _mm256_and_si256(zm, vq));
                    v = _mm256_sub_epi32(v, _mm256_and_si256(cm, vq));
                }
                __m256i neg = _mm256_cmpgt_epi32(zero, v); /* v < 0   */
                __m256i ge = _mm256_cmpgt_epi32(v, v2qm1); /* v >= 2q */
                v = _mm256_add_epi32(v, _mm256_and_si256(neg, v2q));
                v = _mm256_sub_epi32(v, _mm256_and_si256(ge, v2q));
                _mm256_storeu_si256((__m256i *)&wt[i].coeffs[k], v);
            }
        }
#else
        for (int k = 0; k < N; k++) {
            int32_t v = 2 * t.coeffs[k];
            if (i == 0) {
                int32_t ck = c->coeffs[k] & 1;
                v += (int32_t)((int64_t)Q * (z1[0].coeffs[k] & 1));
                v -= (int32_t)((int64_t)Q * ck);
            }
            wt[i].coeffs[k] = reduce_mod_2q(v);
        }
#endif
    }
}

/* ==========================================================================
 * AVX2 NTT-domain primitives + 4-way uniform expansion.  poly_ntt /
 * poly_invntt_tomont / poly_pointwise_montgomery route to the vendored NTT
 * assembly; poly_uniform_4x / poly_uniform_chunkx8 use shake128x4.  The
 * scalar counterparts above are guarded out under USE_AVX2_NTT /
 * _SHAKE4X.
 * ==========================================================================
 */
/* Canonical pk/sk wire order is ref/ ntt16 order, even in AVX builds. */
void poly_ntt_canonical(poly16 *out, const poly *in)
{
    poly_freeze_to_u16_avx((int16_t *)out->coeffs, in->coeffs);
    ntt16(out->coeffs);
}

void poly_ntt_import(poly16 *out, const poly16 *in)
{
    /* Canonical uint16 [0,q) -> avx NTT slot layout, NO int32 staging: the
     * 16-bit unpack consumes the samples directly.  out may alias in
     * (expand_A0_ntt unpacks A0h in place); the shuffle is in-place safe
     * -- each 128-coeff half loads all its registers before any store. */
    nttunpack16_avx((int16_t *)out->coeffs, (const int16_t *)in->coeffs);
}

/* Fused bitunpack+freeze of the M packed b-hat polys (see polyvec.h).  The
 * conditional subtract of q is byte-exact with freeze() over [0, 2^QBITS)
 * (KAT-neutral).  QBITS==16: each coeff is a plain LE uint16, so the whole
 * job is load + masked-subtract.  QBITS<16: SIMD 15-bit field extraction
 * (pshufb spread + per-lane variable shift), 8 coeffs/iteration. */
void unpack_pk_bn(poly16 bn[M], const uint8_t *packed)
{
    for (int i = 0; i < M; i++) {
        const uint8_t *p = packed + (size_t)i * POLYPK_PACKEDBYTES;
        int k = 0;
#if QBITS == 16
#    ifdef __AVX2__
        /* cond-sub q via the min idiom: for c in [0,2^16), min(c, c-q mod
         * 2^16) == (c>=q ? c-q : c) == freeze(c).  N is a multiple of 16.
         */
        const __m256i vq = _mm256_set1_epi16((short)Q);
        for (; k + 16 <= N; k += 16) {
            __m256i c = _mm256_loadu_si256((const __m256i *)(p + 2 * k));
            c = _mm256_min_epu16(c, _mm256_sub_epi16(c, vq));
            _mm256_storeu_si256((__m256i *)(bn[i].coeffs + k), c);
        }
#    endif
        for (; k < N; k++) {
            uint32_t c =
                (uint32_t)p[2 * k] | ((uint32_t)p[2 * k + 1] << 8);
            uint32_t ge = (uint32_t)0 - (uint32_t)(c >= (uint32_t)Q);
            bn[i].coeffs[k] = (uint16_t)(c - (ge & (uint32_t)Q));
        }
#else /* QBITS == 15 (mode 180): straddled fields, real extraction */
#    ifdef __AVX2__
        /* 8 coeffs occupy exactly 15 bytes (15*8=120 bits).  Spread each
         * field's <=3 source bytes into a 32-bit lane, shift right by the
         * field's start-bit-in-byte, mask 15 bits, cond-sub q, narrow to
         * uint16.  Last group is scalar so the 16-byte load never
         * over-reads the final poly. */
        const __m256i shuf = _mm256_setr_epi8(
            0, 1, 2, -1, 1, 2, 3, -1, 3, 4, 5, -1, 5, 6, 7, -1, 7, 8, 9,
            -1, 9, 10, 11, -1, 11, 12, 13, -1, 13, 14, 15, -1);
        const __m256i shifts = _mm256_setr_epi32(0, 7, 6, 5, 4, 3, 2, 1);
        const __m256i mask15 = _mm256_set1_epi32(0x7FFF);
        const __m256i vqd = _mm256_set1_epi32(Q);
        const __m256i vqm1d = _mm256_set1_epi32(Q - 1);
        for (; k + 8 <= N - 8; k += 8) {
            __m128i raw = _mm_loadu_si128(
                (const __m128i *)(p + (size_t)(k >> 3) * 15));
            __m256i s = _mm256_shuffle_epi8(
                _mm256_broadcastsi128_si256(raw), shuf);
            s = _mm256_and_si256(_mm256_srlv_epi32(s, shifts), mask15);
            __m256i ge = _mm256_cmpgt_epi32(s, vqm1d); /* s >= q */
            s = _mm256_sub_epi32(s, _mm256_and_si256(ge, vqd));
            /* 8x32 -> 8x16 contiguous */
            __m256i pk16 = _mm256_packus_epi32(s, _mm256_setzero_si256());
            pk16 = _mm256_permute4x64_epi64(pk16, 0x08);
            _mm_storeu_si128((__m128i *)(bn[i].coeffs + k),
                             _mm256_castsi256_si128(pk16));
        }
#    endif
        /* scalar tail (also the whole loop when __AVX2__ is off): byte-
         * aligned because 8 coeffs == 15 bytes, so it resumes at
         * p[15*k/8]. */
        {
            uint32_t acc = 0;
            int accbits = 0;
            size_t inpos = (size_t)(k >> 3) * 15;
            for (; k < N; k++) {
                while (accbits < QBITS) {
                    acc |= (uint32_t)p[inpos++] << accbits;
                    accbits += 8;
                }
                uint32_t c = acc & ((1u << QBITS) - 1u);
                acc >>= QBITS;
                accbits -= QBITS;
                uint32_t ge = (uint32_t)0 - (uint32_t)(c >= (uint32_t)Q);
                bn[i].coeffs[k] = (uint16_t)(c - (ge & (uint32_t)Q));
            }
        }
#endif
    }
}

/* Lock the reinterpret-cast contract used below. */
_Static_assert(sizeof(poly16) == 2 * (size_t)N,
               "poly16 must be exactly N uint16_t for the int16_t* cast");

void poly_ntt(poly16 *out, const poly *in)
{
    poly_freeze_to_u16_avx((int16_t *)out->coeffs, in->coeffs);
#if LITHIUM_MODE == 260
    ntt_avx((int16_t *)out->coeffs, ntt_qdata, ntt_zetas_fwd,
            ntt_cross_fwd, ntt_ninv);
#else
    ntt_avx((int16_t *)out->coeffs, ntt_qdata, ntt_zetas_fwd, ntt_z0,
            ntt_scale);
#endif
}

void poly_ntt_small(poly16 *out, const poly *in)
{
    poly_freeze_small_to_u16_avx((int16_t *)out->coeffs, in->coeffs);
#if LITHIUM_MODE == 260
    ntt_avx((int16_t *)out->coeffs, ntt_qdata, ntt_zetas_fwd,
            ntt_cross_fwd, ntt_ninv);
#else
    ntt_avx((int16_t *)out->coeffs, ntt_qdata, ntt_zetas_fwd, ntt_z0,
            ntt_scale);
#endif
}

void poly_invntt_tomont(poly *out, const poly16 *in)
{
    poly16 t;
    memcpy(&t, in, sizeof t);
#if LITHIUM_MODE == 260
    invntt_tomont_avx((int16_t *)t.coeffs, ntt_qdata, ntt_zetas_inv,
                      ntt_cross_inv, ntt_ninv);
#else
    invntt_tomont_avx((int16_t *)t.coeffs, ntt_qdata, ntt_zetas_inv,
                      ntt_z0inv, ntt_scale);
#endif
    for (int i = 0; i < N; i++)
        out->coeffs[i] = (int32_t)t.coeffs[i];
}

void poly_pointwise_montgomery(poly16 *r, const poly16 *a, const poly16 *b)
{
#if LITHIUM_MODE == 260
    pointwise_avx((int16_t *)r->coeffs, (const int16_t *)a->coeffs,
                  (const int16_t *)b->coeffs, ntt_qdata);
#elif LITHIUM_MODE == 120
    basemul_avx((int16_t *)r->coeffs, ntt_qdata,
                (const int16_t *)a->coeffs, (const int16_t *)b->coeffs,
                ntt_basemul_gamma);
#else /* LITHIUM_MODE == 180 (q30977, lazy) */
    basemul_avx((int16_t *)r->coeffs, (const int16_t *)a->coeffs,
                (const int16_t *)b->coeffs, ntt_qdata, ntt_basemul_gamma);
    reduce_avx(
        (int16_t *)r->coeffs); /* [0,2q) -> [0,q): the ONE reduce site */
#endif
}

void poly_uniform_4x(poly *a0, poly *a1, poly *a2, poly *a3,
                     const uint8_t seed[SEEDBYTES], uint16_t nonce0,
                     uint16_t nonce1, uint16_t nonce2, uint16_t nonce3)
{
    poly *a[4] = {a0, a1, a2, a3};
    uint16_t nonce[4] = {nonce0, nonce1, nonce2, nonce3};

    uint8_t inbuf[4][SEEDBYTES + 2];
    for (int j = 0; j < 4; j++) {
        memcpy(inbuf[j], seed, SEEDBYTES);
        inbuf[j][SEEDBYTES] = (uint8_t)(nonce[j] & 0xff);
        inbuf[j][SEEDBYTES + 1] = (uint8_t)(nonce[j] >> 8);
    }

    keccakx4_state st;
    shake128x4_absorb_once(&st, inbuf[0], inbuf[1], inbuf[2], inbuf[3],
                           SEEDBYTES + 2);

    const uint32_t mask = (uint32_t)((1u << QBITS) - 1u);
    unsigned ctr[4] = {0, 0, 0, 0};
    uint8_t blk[4][SHAKE128_RATE];

    while (ctr[0] < (unsigned)N || ctr[1] < (unsigned)N ||
           ctr[2] < (unsigned)N || ctr[3] < (unsigned)N) {
        shake128x4_squeezeblocks(blk[0], blk[1], blk[2], blk[3], 1, &st);
        for (int j = 0; j < 4; j++) {
            for (unsigned i = 0;
                 i + 2 <= SHAKE128_RATE && ctr[j] < (unsigned)N; i += 2) {
                uint32_t t =
                    (uint32_t)blk[j][i] | ((uint32_t)blk[j][i + 1] << 8);
                t &= mask;
                if (t < (uint32_t)Q)
                    a[j]->coeffs[ctr[j]++] = (int32_t)t;
            }
        }
    }
}

/* ---- vectorized rejection sampling for the A0 expansion ----------------
 * Kyber rej_uniform_avx-style SIMD rejection (the build has -mbmi2, so the
 * compaction is the table-free pdep/pext path -- no idx[256][8]).  Simpler
 * than Kyber: our candidates are already 16-bit aligned, so we load uint16
 * directly and skip the 12-bit unpack.  Bit-identical to the ref scalar
 * reject (same candidate order, same unsigned q-test, same per-block
 * stop), which keeps A0h -- hence the top-level KAT -- unchanged. */
#ifdef __AVX2__
static inline int compact8(uint16_t *o, __m128i f, unsigned m)
{
    /* m: 8-bit accept mask (bit k = candidate k of this 128-bit half).
     * Gather the even byte-offsets (2k) of accepted lanes into a pshufb
     * control, then compact the surviving uint16 to the front.  Stores a
     * full 8 u16; trailing lanes are garbage the caller overwrites (the
     * bulk loop guarantees >= 8 slots of headroom).  Returns popcount(m).
     */
    uint64_t e = _pdep_u64(m, 0x0101010101010101ULL); /* byte k <- bit k */
    e = (e << 8) - e; /* -> 0x00/0xFF per byte */
    uint64_t off =
        _pext_u64(0x0E0C0A0806040200ULL, e); /* accepted offsets */
    __m128i lo = _mm_cvtsi64_si128((long long)off);
    __m128i hi =
        _mm_add_epi8(lo, _mm_set1_epi8(1));   /* high byte of the u16 */
    __m128i ctrl = _mm_unpacklo_epi8(lo, hi); /* off0,off0+1,off1,... */
    _mm_storeu_si128((__m128i *)o, _mm_shuffle_epi8(f, ctrl));
    return _mm_popcnt_u32(m);
}
#endif

/* Rejection-sample one SHAKE128 block into o[ctr..cap), appending accepted
 * uint16 [0,q).  Returns the new ctr. */
static int reject_block_avx2(uint16_t *o, int ctr, int cap,
                             const uint8_t *blk)
{
    const uint32_t qmask = (uint32_t)((1u << QBITS) - 1u);
    unsigned pos = 0;
#ifdef __AVX2__
    const __m256i maskv = _mm256_set1_epi16((short)qmask);
    const __m256i qm1 = _mm256_set1_epi16((short)(Q - 1));
    /* 16-wide bulk while a full 8-u16 compact store stays inside [0,cap).
     */
    while (pos + 32 <= SHAKE128_RATE && ctr <= cap - 16) {
        __m256i f = _mm256_and_si256(
            _mm256_loadu_si256((const __m256i *)(blk + pos)), maskv);
        /* unsigned (f < Q)  ==  (min_epu16(f, Q-1) == f) */
        __m256i acc = _mm256_cmpeq_epi16(_mm256_min_epu16(f, qm1), f);
        uint32_t good = _pext_u32((uint32_t)_mm256_movemask_epi8(acc),
                                  0x55555555u); /* 1 bit per coeff */
        pos += 32;
        ctr += compact8(o + ctr, _mm256_castsi256_si128(f), good & 0xFFu);
        ctr += compact8(o + ctr, _mm256_extracti128_si256(f, 1),
                        (good >> 8) & 0xFFu);
    }
#endif
    /* scalar tail: the <16 remainder of the block and the final fill to
     * cap
     * -- identical to the ref reject, so the accepted stream is unchanged.
     */
    for (; pos + 2 <= SHAKE128_RATE && ctr < cap; pos += 2) {
        uint32_t t =
            ((uint32_t)blk[pos] | ((uint32_t)blk[pos + 1] << 8)) & qmask;
        if (t < (uint32_t)Q)
            o[ctr++] = (uint16_t)t;
    }
    return ctr;
}

/* A0 in 8 equal coefficient chunks (see polyvec.h): 8 sub-streams
 * SHAKE128(seed||j) via 2x keccak4x (lanes 0-3, 4-7), each lane rejection-
 * sampling per_lane coeffs into out[j*per_lane ..).  Byte-identical per
 * lane to the ref scalar; full 8-way utilisation for every mode. */
void poly_uniform_chunkx8(uint16_t *out, int per_lane,
                          const uint8_t seed[SEEDBYTES])
{
    uint8_t inbuf[8][SEEDBYTES + 2];
    for (int j = 0; j < 8; j++) {
        memcpy(inbuf[j], seed, SEEDBYTES);
        inbuf[j][SEEDBYTES] = (uint8_t)j;
        inbuf[j][SEEDBYTES + 1] = 0;
    }
    keccakx4_state lo, hi;
    shake128x4_absorb_once(&lo, inbuf[0], inbuf[1], inbuf[2], inbuf[3],
                           SEEDBYTES + 2);
    shake128x4_absorb_once(&hi, inbuf[4], inbuf[5], inbuf[6], inbuf[7],
                           SEEDBYTES + 2);

    int ctr[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    uint8_t blk[8][SHAKE128_RATE];
    int done;
    do {
        shake128x4_squeezeblocks(blk[0], blk[1], blk[2], blk[3], 1, &lo);
        shake128x4_squeezeblocks(blk[4], blk[5], blk[6], blk[7], 1, &hi);
        done = 1;
        for (int j = 0; j < 8; j++) {
            ctr[j] = reject_block_avx2(out + (size_t)j * per_lane, ctr[j],
                                       per_lane, blk[j]);
            if (ctr[j] < per_lane)
                done = 0;
        }
    } while (!done);
}
