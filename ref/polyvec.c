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

#include "fips202.h"
#include "ntt.h"
#include "prof.h"
#include "reduce.h"

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

/* signed poly -> [0,q) uint16 -> forward NTT (normal-domain in/out).  This
 * is the scalar reference; the avx2 build maintains its own polyvec.c that
 * routes these three NTT-domain primitives to the in-tree NTT assembly.
 * poly_invert keeps its OWN scalar transforms (see below) so it never
 * depends on the NTT shuffle order. */
void poly_ntt(poly16 *out, const poly *in)
{
    for (int i = 0; i < N; i++)
        out->coeffs[i] = (uint16_t)freeze(in->coeffs[i]);
    ntt16(out->coeffs);
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

void poly_ntt_small(poly16 *out, const poly *in)
{
    for (int i = 0; i < N; i++)
        out->coeffs[i] = small_to_u16_modq(in->coeffs[i]);
    ntt16(out->coeffs);
}

void poly_ntt_canonical(poly16 *out, const poly *in)
{
    poly_ntt(out, in);
}

void poly_ntt_import(poly16 *out, const poly16 *in)
{
    memcpy(out, in, sizeof *out);
}

/* Fused bitunpack+freeze for the M packed b-hat polynomials (see
 * polyvec.h). Scalar reference; AVX2/AVX-512 override this with a vector
 * load + masked subtract.  The conditional subtract of q is byte-exact
 * with freeze() for every packed value in [0, 2^QBITS) because 2^QBITS <
 * 2q -> KAT-neutral. */
void unpack_pk_bn(poly16 bn[M], const uint8_t *packed)
{
    for (int i = 0; i < M; i++) {
        const uint8_t *p = packed + (size_t)i * POLYPK_PACKEDBYTES;
#if QBITS == 16
        /* Each coefficient occupies exactly two bytes (LSB-first), so the
         * packed field is a plain little-endian uint16 -- no bit straddle.
         */
        for (int k = 0; k < N; k++) {
            uint32_t c =
                (uint32_t)p[2 * k] | ((uint32_t)p[2 * k + 1] << 8);
            uint32_t ge = (uint32_t)0 - (uint32_t)(c >= (uint32_t)Q);
            bn[i].coeffs[k] = (uint16_t)(c - (ge & (uint32_t)Q));
        }
#else
        /* QBITS < 16: fields straddle bytes.  Fused LSB-first 15-bit
         * reader straight to uint16 (7+15 < 32, so a uint32 accumulator
         * suffices -- no uint64, no int32 tmp[N], no freeze()). */
        uint32_t acc = 0;
        int accbits = 0;
        size_t inpos = 0;
        for (int k = 0; k < N; k++) {
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
#endif
    }
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

void poly16_add(poly16 *r, const poly16 *a, const poly16 *b)
{
    for (int i = 0; i < N; i++)
        r->coeffs[i] = addm16(a->coeffs[i], b->coeffs[i]);
}

void poly_pointwise_montgomery(poly16 *r, const poly16 *a, const poly16 *b)
{
#if BASE_DEG == 1
    for (int i = 0; i < N; i++)
        r->coeffs[i] = fqmul16(a->coeffs[i], b->coeffs[i]);
#else
    poly16_basemul_montgomery_native(r->coeffs, a->coeffs, b->coeffs);
#endif
}

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
 * replaces it with an optimized multi-lane version that is
 * byte-identical per lane.  All four lanes are real (expand_A0 handles
 * any remainder with scalar poly_uniform). */
void poly_uniform_4x(poly *a0, poly *a1, poly *a2, poly *a3,
                     const uint8_t seed[SEEDBYTES], uint16_t nonce0,
                     uint16_t nonce1, uint16_t nonce2, uint16_t nonce3)
{
    poly_uniform(a0, seed, nonce0);
    poly_uniform(a1, seed, nonce1);
    poly_uniform(a2, seed, nonce2);
    poly_uniform(a3, seed, nonce3);
}

/* A0 in 8 equal coefficient chunks (see polyvec.h).  Scalar: 8 SHAKE128
 * sub-streams; the avx2/avx512 builds replace this with 2x keccak4x /
 * 1x keccak8x (byte-identical per lane). */
void poly_uniform_chunkx8(uint16_t *out, int per_lane,
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
        uint16_t *o = out + (size_t)j * per_lane;
        while (ctr < per_lane) {
            shake128_squeezeblocks(blk, 1, &st);
            for (unsigned i = 0; i + 2 <= SHAKE128_RATE && ctr < per_lane;
                 i += 2) {
                uint32_t t =
                    (uint32_t)blk[i] | ((uint32_t)blk[i + 1] << 8);
                t &= mask;
                if (t < (uint32_t)Q)
                    o[ctr++] = (uint16_t)t;
            }
        }
    }
}

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
    int64_t s = 0;
    for (int i = 0; i < N; i++) {
        int32_t c = freeze(a->coeffs[i]);
        int32_t gt = ((Q / 2) - c) >> 31;
        c -= gt & Q;
        s += (int64_t)c * c;
    }
    return s;
}

/* ============================================================
 * NTWE matrix/vector operations
 * ============================================================ */
void expand_A0_ntt(poly16 A0h[M][ELL], const uint8_t rho[SEEDBYTES])
{
    /* Seed-to-NTT (Dilithium-style, docs/design-notes/PreNTT.md): sample
     * uniform [0,q) and store the samples DIRECTLY as the NTT-domain
     * poly16
     * -- no ntt16.  The samples ARE A0h; coefficient-domain A0 :=
     * INTT(A0h) is never materialized (every consumer wants A0h).
     * Uniformity is preserved because NTT is a bijection R_q <-> T_q.
     *
     * A0h[M][ELL] is one contiguous uint16_t[M*ELL*N] block and every
     * q < 2^16, so sample the 16-bit coeffs straight into it (poly
     * p = i*ELL+j occupies the slice [p*N .. (p+1)*N)).  ref's ntt16 order
     * IS the canonical/standard order, so no slot permutation is needed.
     */
    _Static_assert(sizeof(poly16) == 2 * (size_t)N,
                   "(uint16_t*)A0h cast needs poly16 == N uint16_t");
    poly_uniform_chunkx8((uint16_t *)A0h, (M * ELL * N) / 8, rho);
}

int compute_b(poly b[M], const poly16 A0h[M][ELL], const poly s0[ELL],
              const poly e[M], const poly *f)
{
    poly finv;
    if (poly_invert(&finv, f) != 0)
        return -1; /* f not invertible */

    poly16 s0h[ELL], finv_ntt, acc_h, prod, eh;
    poly acc_n;
    poly_ntt(&finv_ntt, &finv);
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
        poly_ntt(&eh, &acc_n);
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
        for (int k = 0; k < N; k++) {
            int32_t v = 2 * t.coeffs[k];
            if (i == 0) {
                /* q*x0 mod 2q == q*(x0 mod 2): use the TRUE coefficient
                 * parity, not freeze(x0)&1 (freeze adds odd q to
                 * negatives, flipping parity). */
                v += (int32_t)((int64_t)Q * (x[0].coeffs[k] & 1));
            }
            /* v in [0,3q): mod-2q is one conditional subtract (matches the
             * avx2/avx512 lift; no general Barrett reduce_mod_2q). */
            v -= ((DQ - 1 - v) >> 31) & DQ;
            w[i].coeffs[k] = v;
        }
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
        /* compute_t ends in poly_invntt_tomont, so t is already in [0,q);
         * unlike mat_mul_2q there is no intervening poly_add, so the
         * freeze was a no-op -- dropped. */
        for (int k = 0; k < N; k++) {
            int32_t v = 2 * t.coeffs[k];
            if (i == 0) {
                /* true parities (raw coeffs), not freeze()&1 */
                int32_t ck = c->coeffs[k] & 1;
                v += (int32_t)((int64_t)Q * (z1[0].coeffs[k] & 1));
                v -= (int32_t)((int64_t)Q * ck);
            }
            /* v in [-q,3q): mod-2q via one conditional add + one
             * conditional subtract (matches avx2/avx512; no Barrett). */
            v += (v >> 31) & DQ;
            v -= ((DQ - 1 - v) >> 31) & DQ;
            wt[i].coeffs[k] = v;
        }
    }
}
