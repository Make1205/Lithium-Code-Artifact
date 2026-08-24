/*
 * sign.c - IRS KeyGen / Sign / Verify (compressed variant, Lithium.md
 * Fig.4).
 *
 * NTWE keys: s = [f, s0, e], f = 2 f0 + 1, b = (A0 s0 + e) f^{-1};
 *            A = [qj - 2b, 2 A0, 2 I] implicit, A s = q j (mod 2q).
 * CSign: y <- D_{Z,r}^K; w = A y mod 2q; c = H(HighBits(w), LSB(w[0]),
 * mu); z = y + s c' via iterative rejection sampling (irs.c);
 *        ||z|| <= B_s; split z=(z1,z2); h = MakeHint(w, w-2 z2);
 * sig=(c,h,z1). CVerify: wtilde' = A1 z1 - q c j mod 2q; w' =
 * UseHint(wtilde', h); c' = H(w', LSB(z0-c), mu); accept iff c'==c and
 * ||[z1,z2']|| <= B_v.
 *
 * RAW signature packing (packing.c) for validation; rANS replaces it
 * later.
 */
#include "sign.h"

#include <stdint.h>
#include <string.h>

#include "fips202.h"
#include "irs.h"
#include "packing.h"
#include "params.h"
#include "polyvec.h"
#include "prof.h"
#include "randombytes.h"
#include "reduce.h"
#include "rounding.h"
#include "sampler.h"
/* AVX-512/AVX2 highbits/makehint fast paths. The scalar #else path is the
 * portable reference. Inputs are pre-bounded, so Barrett reduction is a
 * no-op here and the vector path reduces to exact shifts/masks. */
#if defined(__AVX512F__) || defined(__AVX2__)
#    include <immintrin.h>
#    if TAU == 1024
#        define TAU_LOG2 10
#    elif TAU == 512
#        define TAU_LOG2 9
#    else
#        error "TAU must be a power of two (512/1024)"
#    endif
#endif
#if defined(__AVX512F__)
/* HighBits of 16 coeffs w in [0,2q): (w+tau/2)>>log2(tau), then mod
 * W1_MAX. */
static inline __m512i highbits16(__m512i w)
{
    __m512i b = _mm512_srli_epi32(
        _mm512_add_epi32(w, _mm512_set1_epi32(TAU / 2)), TAU_LOG2);
    __mmask16 ge =
        _mm512_cmpgt_epi32_mask(b, _mm512_set1_epi32(W1_MAX - 1));
    return _mm512_mask_sub_epi32(b, ge, b, _mm512_set1_epi32(W1_MAX));
}
#elif defined(__AVX2__)
/* HighBits of 8 coeffs w in [0,2q): (w+tau/2)>>log2(tau), then mod W1_MAX.
 */
static inline __m256i highbits8(__m256i w)
{
    __m256i b = _mm256_srli_epi32(
        _mm256_add_epi32(w, _mm256_set1_epi32(TAU / 2)), TAU_LOG2);
    __m256i ge = _mm256_cmpgt_epi32(b, _mm256_set1_epi32(W1_MAX - 1));
    return _mm256_sub_epi32(
        b, _mm256_and_si256(ge, _mm256_set1_epi32(W1_MAX)));
}
#endif

#ifdef TEST_RANS_STATS
uint64_t lithium_test_rans_pack_calls = 0;
uint64_t lithium_test_rans_pack_fails = 0;
#endif

/* challenge: kappa ones, Dilithium-style Fisher-Yates from SHAKE256(seed).
 */
static void poly_challenge(poly *c, const uint8_t seed[CHASH_BYTES])
{
    uint8_t buf[SHAKE256_RATE];
    keccak_state st;
    shake256_init(&st);
    shake256_absorb(&st, seed, CHASH_BYTES);
    shake256_finalize(&st);
    shake256_squeezeblocks(buf, 1, &st);
    unsigned pos = 0;
    for (int k = 0; k < N; k++)
        c->coeffs[k] = 0;
    for (int i = N - KAPPA; i < N; i++) {
        unsigned b;
        do {
            if (pos >= SHAKE256_RATE) {
                shake256_squeezeblocks(buf, 1, &st);
                pos = 0;
            }
            b = buf[pos++];
        } while (b > (unsigned)i);
        c->coeffs[i] = c->coeffs[b];
        c->coeffs[b] = 1;
    }
}

static void hash_challenge(uint8_t seed[CHASH_BYTES],
                           const int32_t buckets[M][N],
                           const uint8_t w0bits[(N + 7) / 8],
                           const uint8_t mu[CRHBYTES])
{
    keccak_state st;
    uint8_t bbuf[M * ((N * W1_BITS + 7) / 8)];
    /* Pack M*N buckets, W1_BITS each, LSB-first into one
     * continuous bitstream (N*W1_BITS is a multiple of 8 in every
     * mode, so this is byte-identical to the per-coeff bit loop --
     * KAT-neutral).  A bit-accumulator writer replaces the old per-bit `if
     * (x & (1<<b))`: its byte-emit schedule depends only on the iteration
     * count (not the data), so it has NO data-dependent branches -> no
     * misprediction stall, and no memset.  Portable scalar win (helps ref
     * too). */
    _Static_assert(W1_BITS <= 8, "bit-writer assumes <=8 bits/coeff");
    _Static_assert((N * W1_BITS) % 8 == 0, "row must be byte-aligned");
    {
        const uint32_t mask = (1u << W1_BITS) - 1u;
        uint64_t acc = 0;
        int accbits = 0;
        size_t outpos = 0;
        for (int i = 0; i < M; i++)
            for (int k = 0; k < N; k++) {
                acc |= (uint64_t)((uint32_t)buckets[i][k] & mask)
                       << accbits;
                accbits += W1_BITS;
                while (accbits >= 8) {
                    bbuf[outpos++] = (uint8_t)acc;
                    acc >>= 8;
                    accbits -= 8;
                }
            }
    }
    shake256_init(&st);
    shake256_absorb(&st, bbuf, sizeof bbuf);
    shake256_absorb(&st, w0bits, (N + 7) / 8);
    shake256_absorb(&st, mu, CRHBYTES);
    shake256_finalize(&st);
    shake256_squeeze(seed, CHASH_BYTES, &st);
}

/* Gaussian y-vector: ALL K polys sampled from ONE unified keccak4x
 * combined stream (4 sub-streams of a single base nonce, consumed
 * back-to-back) -- full keccak4x throughput with near-zero waste (only the
 * final partial 4-block unit is unused).  ref and avx2 produce identical
 * output. */
static void sample_gauss_polys(poly *fse, int kpolys,
                               const uint8_t seed[SKSEEDBYTES],
                               uint64_t base)
{
    int32_t *rp[K];
    for (int c = 0; c < kpolys; c++)
        rp[c] = fse[c].coeffs; /* sample int32 straight into poly.coeffs */
    sample_gauss_vec(rp, kpolys, seed, base);
}

/* Key-noise vector: ALL K polys sampled from ONE unified combined stream
 * (GAUSS_COMB_LANES sub-streams of a single base nonce, consumed
 * back-to-back) -- full Keccak throughput with no idle lanes (vs the old
 * per-poly sample_polys_4x, which padded the K-mod-lanes remainder).  ref
 * and avx2/avx512 produce identical output. */
static void sample_noise_polys(poly *fse, int kpolys,
                               const uint8_t seed[SKSEEDBYTES],
                               uint64_t base)
{
    int32_t *rp[K];
    for (int c = 0; c < kpolys; c++)
        rp[c] = fse[c].coeffs; /* sample int32 straight into poly.coeffs */
    sample_noise_vec(rp, kpolys, seed, base);
}

_Static_assert(K == 1 + ELL + M, "secret vector layout mismatch");

int crypto_sign_keypair(uint8_t *pk, uint8_t *sk)
{
    /* expanded = rho (SEEDBYTES) || rhoprime (SKSEEDBYTES) || key
     * (SKSEEDBYTES). The master seed is SEEDBYTES of randomness with the
     * (ELL, M) parameter pair appended so distinct parameter sets derive
     * independent keys from the same RNG draw. rhoprime is passed by
     * pointer to the noise sampler, which copies it internally. */
    uint8_t seedbuf[SEEDBYTES + 2];
    uint8_t expanded[SEEDBYTES + SKSEEDBYTES + SKSEEDBYTES];
    randombytes(seedbuf, SEEDBYTES);
    seedbuf[SEEDBYTES + 0] = (uint8_t)ELL;
    seedbuf[SEEDBYTES + 1] = (uint8_t)M;
    shake256(expanded, sizeof expanded, seedbuf, SEEDBYTES + 2);
    const uint8_t *rho = expanded, *rhoprime = expanded + SEEDBYTES,
                  *key = expanded + SEEDBYTES + SKSEEDBYTES;

    poly16 A0h[M][ELL];
    PROF_START(t_ka0);
    expand_A0_ntt(A0h, rho); /* A0 generated directly in NTT domain */
    PROF_STOP(PT_KG_EXPAND_A0, t_ka0);

    /* Noise vector [f0, s0(ell), e(m)] held in ONE contiguous poly[K] so
     * the sampler writes int32 straight into it (no int16 staging);
     * f0/s0/e are readability aliases into it. */
    poly fse[K], f, b[M], s[K];
    poly *f0 = &fse[0], *s0 = &fse[1], *e = &fse[1 + ELL];
    /* Keygen retries are a one-time provisioning cost in deployments such
     * as modern certificate issuance and TLS 1.3 server authentication:
     * the generated key is expected to be reused long-term, so this loop
     * is not a per-message leakage channel. */
    for (unsigned attempt = 0;; attempt++) {
        uint64_t base = (uint64_t)attempt << 8;
        /* Noise vector [f0, s0(ell), e(m)] sampled from ONE unified
         * combined stream (no idle Keccak lanes; see sample_noise_polys).
         */
        PROF_START(t_kn);
        sample_noise_polys(fse, K, rhoprime, base);
        PROF_STOP(PT_KG_NOISE, t_kn);
        for (int k = 0; k < N; k++)
            f.coeffs[k] = 2 * f0->coeffs[k] + (k == 0 ? 1 : 0);

        s[0] = f;
        for (int j = 0; j < ELL; j++)
            s[1 + j] = s0[j];
        for (int i = 0; i < M; i++)
            s[1 + ELL + i] = e[i];
        PROF_START(t_knm);
        int64_t nrm = 0;
        for (int c = 0; c < K; c++)
            nrm += poly_sqnorm(&s[c]);
        PROF_STOP(PT_KG_NORM, t_knm);
        if (nrm < BK_LOW_SQ || nrm > BK_SQ)
            continue; /* keygen norm window */
        PROF_START(t_kcb);
        int binv = compute_b(b, A0h, s0, e, &f);
        PROF_STOP(PT_KG_COMPUTEB, t_kcb);
        if (binv != 0)
            continue; /* f not invertible */
        break;
    }

    PROF_START(t_kpk);
    poly16 bn[M];
    for (int i = 0; i < M; i++)
        poly_ntt_canonical(&bn[i], &b[i]);
    pack_pk(pk, rho, bn);
    uint8_t tr[TRBYTES];
    shake256(tr, TRBYTES, pk, CRYPTO_PUBLICKEYBYTES);
    pack_sk(sk, rho, tr, key, f0, s0, e, bn);
    PROF_STOP(PT_KG_PACK, t_kpk);
    return 0;
}

int crypto_sign_signature(uint8_t *sig, size_t *siglen, const uint8_t *m,
                          size_t mlen, const uint8_t *sk)
{
    uint8_t rho[SEEDBYTES], tr[TRBYTES], key[SKSEEDBYTES];
    poly s[K];
    poly16 bn[M];
    unpack_sk(rho, tr, key, &s[0], &s[1], &s[1 + ELL], bn, sk);
    for (int k = 0; k < N; k++)
        s[0].coeffs[k] = 2 * s[0].coeffs[k] + (k == 0 ? 1 : 0);

    PROF_CTX(PC_SETUP);
    PROF_START(t_setup);
    uint8_t mu[CRHBYTES];
    {
        keccak_state st;
        shake256_init(&st);
        shake256_absorb(&st, tr, TRBYTES);
        shake256_absorb(&st, m, mlen);
        shake256_finalize(&st);
        shake256_squeeze(mu, CRHBYTES, &st);
    }

    uint8_t rnd[RNDBYTES], rpbuf[SKSEEDBYTES + RNDBYTES + CRHBYTES];
    randombytes(rnd, RNDBYTES);
    memcpy(rpbuf, key, SKSEEDBYTES);
    memcpy(rpbuf + SKSEEDBYTES, rnd, RNDBYTES);
    memcpy(rpbuf + SKSEEDBYTES + RNDBYTES, mu, CRHBYTES);
    uint8_t rhoprime[SKSEEDBYTES];
    shake256(rhoprime, SKSEEDBYTES, rpbuf, sizeof rpbuf);
    PROF_STOP(PT_SETUP, t_setup);

    poly16 A0h[M][ELL], bh[M];
    PROF_CTX(PC_A0);
    PROF_START(t_a0);
    expand_A0_ntt(A0h, rho); /* A0 directly in NTT domain (Pre-NTT) */
    for (int i = 0; i < M; i++)
        poly_ntt_import(&bh[i],
                        &bn[i]); /* import canonical bn to native bh */
    PROF_STOP(PT_EXPAND_A0, t_a0);

    for (unsigned attempt = 0;; attempt++) {
        poly z[K], w[M];
        /* y is only read before sign_irs mutates z.  Make y a readability
         * alias for z, sample directly into z, and avoid one full
         * poly[K] memcpy per signing attempt. */
        poly *y = z;
        uint64_t base = (uint64_t)attempt * K;
        /* Gaussian y vector with nonces base+0..base+K-1, batched 4-way
         * (KAT-identical to K sequential sample_gauss_N calls).  z is
         * contiguous, so sample straight into it (y is its alias). */
        PROF_CTX(PC_GAUSS);
        PROF_START(t_y);
        sample_gauss_polys(z, K, rhoprime, base);
        PROF_STOP(PT_SAMPLE_Y, t_y);

        PROF_START(t_commit);
        mat_mul_2q(w, A0h, bh, y);
        PROF_STOP(PT_COMMIT, t_commit);

        uint8_t seedc[CHASH_BYTES];
        poly c;
        {
            int32_t buckets[M][N];
            uint8_t w0bits[(N + 7) / 8];
            PROF_START(t_hb);
            memset(w0bits, 0, sizeof w0bits);
#if defined(__AVX512F__)
            for (int i = 0; i < M; i++)
                for (int k = 0; k < N; k += 16)
                    _mm512_storeu_si512(
                        (__m512i *)&buckets[i][k],
                        highbits16(_mm512_loadu_si512(
                            (const __m512i *)&w[i].coeffs[k])));
#elif defined(__AVX2__)
            for (int i = 0; i < M; i++)
                for (int k = 0; k < N; k += 8)
                    _mm256_storeu_si256(
                        (__m256i *)&buckets[i][k],
                        highbits8(_mm256_loadu_si256(
                            (const __m256i *)&w[i].coeffs[k])));
#else
            for (int i = 0; i < M; i++)
                for (int k = 0; k < N; k++)
                    buckets[i][k] = highbits(w[i].coeffs[k]);
#endif
            for (int k = 0; k < N; k++) {
                /* w[0] is already in [0,2q) (mat_mul_2q stored it via
                 * reduce_mod_2q), so lsb()==&1 -- skip the redundant
                 * Barrett (KAT-neutral; the verify side's lsb(z0-c) is NOT
                 * pre-bounded, so it keeps lsb()). */
                uint8_t bit = (uint8_t)(1u << (k & 7));
                uint8_t mask =
                    (uint8_t)(0u - ((uint32_t)w[0].coeffs[k] & 1u));
                w0bits[k >> 3] |= bit & mask;
            }
            PROF_STOP(PT_HIGHBITS, t_hb);

            PROF_CTX(PC_CHALLENGE);
            PROF_START(t_ch);
            hash_challenge(seedc, buckets, w0bits, mu);
            poly_challenge(&c, seedc);
            PROF_STOP(PT_CHALLENGE, t_ch);
        }

        PROF_CTX(PC_IRS);
        PROF_START(t_irs);
        {
            keccak_state rng;
            {
                uint8_t rbuf[SKSEEDBYTES + CRHBYTES + CHASH_BYTES];
                memcpy(rbuf, key, SKSEEDBYTES);
                memcpy(rbuf + SKSEEDBYTES, mu, CRHBYTES);
                memcpy(rbuf + SKSEEDBYTES + CRHBYTES, seedc, CHASH_BYTES);
                shake256_init(&rng);
                shake256_absorb(&rng, rbuf, sizeof rbuf);
                shake256_finalize(&rng);
            }
            sign_irs(z, &c, s, NULL, &rng);
        }
        PROF_STOP(PT_SIGN, t_irs);

        int64_t nrm = 0;
        PROF_START(t_nc);
        for (int cc = 0; cc < K; cc++)
            nrm += poly_sqnorm(&z[cc]);
        PROF_STOP(PT_NORMCHECK, t_nc);
        if (nrm > BS_SQ)
            continue; /* ||z|| > B_s; near-2^-100 retry by parameters */

        const poly *z1 = z; /* pack_sig reads only the z prefix. */
        poly *h = w;
        PROF_START(t_mh);
#if defined(__AVX512F__)
        {
            const __m512i dqv = _mm512_set1_epi32(DQ);
            const __m512i W1 = _mm512_set1_epi32(W1_MAX);
            const __m512i W1H = _mm512_set1_epi32(W1_MAX / 2);
            const __m512i zero = _mm512_setzero_si512();
            for (int i = 0; i < M; i++) {
                poly *z2i = &z[1 + ELL + i];
                for (int k = 0; k < N; k += 16) {
                    __m512i wv = _mm512_loadu_si512(
                        (const __m512i *)&w[i].coeffs[k]);
                    __m512i z2v = _mm512_loadu_si512(
                        (const __m512i *)&z2i->coeffs[k]);
                    /* wt = reduce_mod_2q(w - 2*z2), |w-2z2|<2q -> two
                     * masked corrections cover [2q,4q) and (-2q,0). */
                    __m512i d =
                        _mm512_sub_epi32(wv, _mm512_add_epi32(z2v, z2v));
                    d = _mm512_mask_sub_epi32(
                        d,
                        _mm512_cmpgt_epi32_mask(d,
                                                _mm512_set1_epi32(DQ - 1)),
                        d, dqv);
                    d = _mm512_mask_add_epi32(
                        d, _mm512_cmpgt_epi32_mask(zero, d), d,
                        dqv); /* d<0 */
                    /* makehint = centermod(highbits(w)-highbits(wt),
                     * W1_MAX) */
                    __m512i diff =
                        _mm512_sub_epi32(highbits16(wv), highbits16(d));
                    diff = _mm512_mask_add_epi32(
                        diff, _mm512_cmpgt_epi32_mask(zero, diff), diff,
                        W1);
                    diff = _mm512_mask_sub_epi32(
                        diff, _mm512_cmpgt_epi32_mask(diff, W1H), diff,
                        W1);
                    _mm512_storeu_si512((__m512i *)&h[i].coeffs[k], diff);
                }
            }
        }
#elif defined(__AVX2__)
        {
            const __m256i dqv = _mm256_set1_epi32(DQ);
            const __m256i W1 = _mm256_set1_epi32(W1_MAX);
            const __m256i W1H = _mm256_set1_epi32(W1_MAX / 2);
            for (int i = 0; i < M; i++) {
                poly *z2i = &z[1 + ELL + i];
                for (int k = 0; k < N; k += 8) {
                    __m256i wv = _mm256_loadu_si256(
                        (const __m256i *)&w[i].coeffs[k]);
                    __m256i z2v = _mm256_loadu_si256(
                        (const __m256i *)&z2i->coeffs[k]);
                    /* wt = reduce_mod_2q(w - 2*z2), |w-2z2| < 2q (z2
                     * small) -> two masked corrections cover [2q,4q) and
                     * (-2q,0). */
                    __m256i d =
                        _mm256_sub_epi32(wv, _mm256_add_epi32(z2v, z2v));
                    d = _mm256_sub_epi32(
                        d,
                        _mm256_and_si256(_mm256_cmpgt_epi32(
                                             d, _mm256_set1_epi32(DQ - 1)),
                                         dqv));
                    d = _mm256_add_epi32(
                        d,
                        _mm256_and_si256(_mm256_srai_epi32(d, 31), dqv));
                    /* makehint = centermod(highbits(w) - highbits(wt),
                     * W1_MAX) */
                    __m256i diff =
                        _mm256_sub_epi32(highbits8(wv), highbits8(d));
                    diff = _mm256_add_epi32(
                        diff,
                        _mm256_and_si256(_mm256_srai_epi32(diff, 31), W1));
                    diff = _mm256_sub_epi32(
                        diff, _mm256_and_si256(
                                  _mm256_cmpgt_epi32(diff, W1H), W1));
                    _mm256_storeu_si256((__m256i *)&h[i].coeffs[k], diff);
                }
            }
        }
#else
        for (int i = 0; i < M; i++) {
            poly *z2i = &z[1 + ELL + i];
            for (int k = 0; k < N; k++) {
                int32_t wi = w[i].coeffs[k];
                int32_t z2c = z2i->coeffs[k]; /* small centered */
                int32_t wt = reduce2q(wi - 2 * z2c);
                h[i].coeffs[k] = makehint(wi, wt);
            }
        }
#endif
        PROF_STOP(PT_MAKEHINT, t_mh);
        PROF_START(t_pk);
#ifdef TEST_RANS_STATS
        lithium_test_rans_pack_calls++;
#endif
        int packrc = pack_sig(sig, &c, h, z1);
        PROF_STOP(PT_PACK, t_pk);
        if (packrc != 0) {
#ifdef TEST_RANS_STATS
            lithium_test_rans_pack_fails++;
#endif
            continue; /* rANS overflow: near-2^-100 retry by parameters */
        }
        *siglen = CRYPTO_BYTES;
        return 0;
    }
}

int crypto_sign_verify(const uint8_t *sig, size_t siglen, const uint8_t *m,
                       size_t mlen, const uint8_t *pk)
{
    if (siglen != CRYPTO_BYTES)
        return -1;
    PROF_START(t_vun);
    uint8_t rho[SEEDBYTES];
    poly16 bn[M];
    unpack_pk(rho, bn, pk);

    poly c, h[M], z1[ELL + 1];
    if (unpack_sig(&c, h, z1, sig) != 0)
        return -1;
    PROF_STOP(PT_VF_UNPACK, t_vun);

    PROF_START(t_vsu);
    uint8_t tr[TRBYTES], mu[CRHBYTES];
    shake256(tr, TRBYTES, pk, CRYPTO_PUBLICKEYBYTES);
    {
        keccak_state st;
        shake256_init(&st);
        shake256_absorb(&st, tr, TRBYTES);
        shake256_absorb(&st, m, mlen);
        shake256_finalize(&st);
        shake256_squeeze(mu, CRHBYTES, &st);
    }
    PROF_STOP(PT_VF_SETUP, t_vsu);

    PROF_START(t_va0);
    poly16 A0h[M][ELL], bh[M];
    expand_A0_ntt(A0h, rho); /* A0 directly in NTT domain (Pre-NTT) */
    for (int i = 0; i < M; i++)
        poly_ntt_import(&bh[i],
                        &bn[i]); /* import canonical bn to native bh */
    PROF_STOP(PT_VF_A0, t_va0);

    PROF_START(t_vmm);
    poly wt[M];
    mat_mul_z1_2q(wt, A0h, bh, z1, &c);
    PROF_STOP(PT_VF_MATMUL, t_vmm);

    PROF_START(t_vhb);
    int32_t buckets[M][N];
    uint8_t w0bits[(N + 7) / 8];
    memset(w0bits, 0, sizeof w0bits);
    for (int i = 0; i < M; i++)
        for (int k = 0; k < N; k++)
            buckets[i][k] = usehint(wt[i].coeffs[k], h[i].coeffs[k]);
    for (int k = 0; k < N; k++) {
        int32_t z0 = z1[0].coeffs[k];
        int32_t ck = c.coeffs[k];
        if (lsb(z0 - ck))
            w0bits[k >> 3] |= (uint8_t)(1u << (k & 7));
    }
    PROF_STOP(PT_VF_HINT, t_vhb);

    PROF_START(t_vch);
    uint8_t seedc[CHASH_BYTES];
    poly cprime;
    hash_challenge(seedc, buckets, w0bits, mu);
    poly_challenge(&cprime, seedc);
    for (int k = 0; k < N; k++)
        if (cprime.coeffs[k] != c.coeffs[k]) {
            PROF_STOP(PT_VF_CHALLENGE, t_vch);
            return -2; /* challenge mismatch */
        }
    PROF_STOP(PT_VF_CHALLENGE, t_vch);

    PROF_START(t_vnm);
    int64_t nrm = 0;
    for (int c2 = 0; c2 < ELL + 1; c2++)
        nrm += poly_sqnorm(&z1[c2]);
    for (int i = 0; i < M; i++) {
        for (int k = 0; k < N; k++) {
            int32_t wval = hbvalue(buckets[i][k]); /* bucket*tau */
            int32_t w0c = (i == 0) ? ((w0bits[k >> 3] >> (k & 7)) & 1) : 0;
            int32_t diff = wval - wt[i].coeffs[k] + w0c;
            /* diff in (-2q, 2q): wval = bucket*tau in [0, 2q-tau], wt in
             * [0,2q), w0c in {0,1}.  So mod-2q is one conditional add --
             * no full Barrett reduce2q (per-coeff on the verify path). */
            int32_t d = diff + ((diff >> 31) & DQ);
            int32_t z2 = d >> 1;
            if (z2 > Q / 2)
                z2 -= Q; /* mod^pm q */
            nrm += (int64_t)z2 * z2;
        }
    }
    PROF_STOP(PT_VF_NORM, t_vnm);
    if (nrm > BV_SQ)
        return -3;
    return 0;
}

int crypto_sign(uint8_t *sm, size_t *smlen, const uint8_t *m, size_t mlen,
                const uint8_t *sk)
{
    size_t siglen = 0;

    memmove(sm + CRYPTO_BYTES, m, mlen);
    if (crypto_sign_signature(sm, &siglen, sm + CRYPTO_BYTES, mlen, sk) !=
        0)
        return -1;
    if (siglen != CRYPTO_BYTES)
        return -1;
    *smlen = siglen + mlen;
    return 0;
}

int crypto_sign_open(uint8_t *m, size_t *mlen, const uint8_t *sm,
                     size_t smlen, const uint8_t *pk)
{
    if (smlen < CRYPTO_BYTES)
        return -1;
    if (crypto_sign_verify(sm, CRYPTO_BYTES, sm + CRYPTO_BYTES,
                           smlen - CRYPTO_BYTES, pk) != 0)
        return -1;
    memmove(m, sm + CRYPTO_BYTES, smlen - CRYPTO_BYTES);
    *mlen = smlen - CRYPTO_BYTES;
    return 0;
}
