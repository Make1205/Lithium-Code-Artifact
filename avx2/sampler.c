/*
 * sampler.c - discrete Gaussian samplers for IRS (ref).
 *
 * BLISS base sampler: k = 128, sigma2 = r/128, 93-bit RCDT (3x31 limbs)
 * from the in-file RCDT tables.  Bernoulli accept exp(-y(y+256x)/(2 r^2))
 * via approx_exp(a_q60), a_q60 = (num * I80) >> 20 with I80 =
 * floor(2^80/(2 r^2)) (gen_approx_exp.py).
 *
 * Noise sampler: D_{Z,1} via the shared 10-entry RCDT, no ApproxExp;
 * magnitude
 * + sign + z=0 half-reject (folds the two zero-signs into one).
 *
 * Unified SHAKE-256 stream per call; mini-batch layout identical between
 * ref and a future AVX2 4x (KAT-compatible).  The GAUSS_BATCH accept tests
 * in a mini-batch are independent, so their exp() values are evaluated
 * 4-way via approx_exp_4x (instruction-level parallelism) -- bit-identical
 * to a per- attempt approx_exp, so sample_gauss_N_4x stays KAT-locked to
 * sample_gauss_N.
 */
#include "sampler.h"

#include <immintrin.h>
#include <string.h>

#include "approx_exp.h"
#include "fips202x4.h" /* keccak4x backend for the merged 4-way samplers */
#include "prof.h"

/* 93-bit reverse-CDT tables.  External `const` (declared in sampler.h) so
 * optimized samplers share the canonical tables.
 * Generated + validated by tools/gen_rcdt.py; regenerate via `make
 * tables`.  Do not edit. */
/* @@AUTOGEN:rcdt@@ BEGIN */
/* 93-bit reverse-CDT tables, 3 x 31-bit limbs (low, mid, high); entry i
 * encodes Pr[X >= i+1] * 2^93 (half-discrete-gaussian).  See
 * tools/gen_rcdt.py. The generator checks byte-exact table regeneration
 * and Renyi bounds. RCDT_Z[RCDT_ENTRIES][3] : BLISS base sampler, sigma2 =
 * r/128 RCDT_NOISE[10][3]           : key-noise sampler, sigma = 1
 * (shared) */
#if LITHIUM_MODE == 120
const uint32_t RCDT_Z[RCDT_ENTRIES][3] = {
    {0x7A81A80CU, 0x1C543995U, 0x47E904C9U},
    {0x526DB044U, 0x0EB15A73U, 0x1C1C254BU},
    {0x66D499C6U, 0x3EEB9B72U, 0x0740B175U},
    {0x6A76244DU, 0x28DEE8C3U, 0x01323963U},
    {0x722B35B5U, 0x5CA79891U, 0x001FAD77U},
    {0x6862CF2DU, 0x0122DAD4U, 0x00020847U},
    {0x0DA12494U, 0x1F1A637AU, 0x00001490U},
    {0x4AC0D4AAU, 0x54FA2079U, 0x0000007FU},
    {0x20D97E2FU, 0x72907420U, 0x00000001U},
    {0x6BFAF40CU, 0x02335EAEU, 0x00000000U},
    {0x61F1476AU, 0x00031F01U, 0x00000000U},
    {0x4E2B723FU, 0x000002B3U, 0x00000000U},
    {0x369F3BB7U, 0x00000001U, 0x00000000U},
    {0x003AC8E8U, 0x00000000U, 0x00000000U},
};
#endif

#if LITHIUM_MODE == 180
const uint32_t RCDT_Z[RCDT_ENTRIES][3] = {
    {0x3E0EF881U, 0x33B014DCU, 0x4CBF84DAU},
    {0x1154B70DU, 0x67DBF313U, 0x22A77549U},
    {0x0F30D703U, 0x5237B27AU, 0x0B553B20U},
    {0x4A623BA8U, 0x7C42BC28U, 0x029E0304U},
    {0x08D8DB0FU, 0x29AEC850U, 0x006B85CBU},
    {0x4D5B61AAU, 0x269E0E6DU, 0x000BDE2CU},
    {0x503F9525U, 0x1FA37291U, 0x0000E51FU},
    {0x40FB6AC4U, 0x2C8F7338U, 0x00000BC1U},
    {0x53D2D77EU, 0x5EF14273U, 0x00000068U},
    {0x69A49D58U, 0x3BF8203AU, 0x00000002U},
    {0x70E8EE5AU, 0x05093725U, 0x00000000U},
    {0x7B4E0978U, 0x000DE26CU, 0x00000000U},
    {0x0D770E6DU, 0x000019DAU, 0x00000000U},
    {0x3F5868BAU, 0x00000020U, 0x00000000U},
    {0x0DC88430U, 0x00000000U, 0x00000000U},
    {0x0007E0E1U, 0x00000000U, 0x00000000U},
};
#endif

#if LITHIUM_MODE == 260
const uint32_t RCDT_Z[RCDT_ENTRIES][3] = {
    {0x135CE1E7U, 0x0CE9F6A3U, 0x5590CA61U},
    {0x51107A58U, 0x63FA476DU, 0x3014BAD6U},
    {0x78F93DEAU, 0x75915F0DU, 0x163E3905U},
    {0x66F3987BU, 0x1F70F6E1U, 0x085880F0U},
    {0x29576A36U, 0x240B0CBDU, 0x02834FEFU},
    {0x3F9AD1B7U, 0x0AB033B3U, 0x009A4A22U},
    {0x76ECDA55U, 0x6751D5D7U, 0x001D5157U},
    {0x265706D2U, 0x7A9C4D13U, 0x00046592U},
    {0x1BB09EFCU, 0x1EDF58E5U, 0x000084DAU},
    {0x64C7828FU, 0x6964CD5DU, 0x00000C50U},
    {0x0FBF9931U, 0x199CB932U, 0x000000E5U},
    {0x62A8475DU, 0x068F7A91U, 0x0000000DU},
    {0x7237064DU, 0x4A72297EU, 0x00000000U},
    {0x4E9A9404U, 0x029837FCU, 0x00000000U},
    {0x562C8A76U, 0x00121871U, 0x00000000U},
    {0x456780E5U, 0x0000629CU, 0x00000000U},
    {0x628E5D0AU, 0x000001A3U, 0x00000000U},
    {0x39BBD11DU, 0x00000005U, 0x00000000U},
    {0x071304D0U, 0x00000000U, 0x00000000U},
    {0x000E3F77U, 0x00000000U, 0x00000000U},
};
#endif

const uint32_t RCDT_NOISE[NOISE_RCDT_ENTRIES][3] = {
    {0x057FB7ABU, 0x2DF4C635U, 0x36FED273U},
    {0x055455E2U, 0x2E17C41BU, 0x0AB7433DU},
    {0x22A8A96CU, 0x45896D14U, 0x00D5F51CU},
    {0x44D9083FU, 0x399BDCFBU, 0x000656E6U},
    {0x36EE7761U, 0x11C429AFU, 0x000011E7U},
    {0x61F254F6U, 0x574A372EU, 0x00000012U},
    {0x7B1A38BCU, 0x03977C96U, 0x00000000U},
    {0x76A6C503U, 0x00008224U, 0x00000000U},
    {0x6397C3ABU, 0x00000006U, 0x00000000U},
    {0x00109FA6U, 0x00000000U, 0x00000000U},
};
/* @@AUTOGEN:rcdt@@ END */

/* ---- BLISS base CDT (sigma2 = r/128), GAUSS_BATCH magnitudes ----
 * The AVX2 build replaces this scalar kernel with the 8-way SIMD version
 * below, which is bit-identical to this borrow chain; it is guarded out
 * here when USE_AVX2_SAMPLER is defined. */
#ifndef USE_AVX2_SAMPLER
int sampler_sigma2(int32_t *z_out, const uint8_t *rand)
{
    for (int s = 0; s < GAUSS_BATCH; s++) {
        int group = s >> 3, lane = s & 7;
        int byte_base = group * 96 + lane * 4;
        uint32_t v0 = load_le32(rand + byte_base) & 0x7FFFFFFFU;
        uint32_t v1 = load_le32(rand + byte_base + 32) & 0x7FFFFFFFU;
        uint32_t v2 = load_le32(rand + byte_base + 64) & 0x7FFFFFFFU;
        int32_t z = 0;
#    if defined(__clang__)
#        pragma clang loop vectorize(disable) interleave(disable)
#    endif
        for (int i = 0; i < RCDT_ENTRIES; i++) {
            uint32_t cc;
            /* exact 31-bit-limb borrow chain: borrow0 from ct_lt_u32 (no
             * incoming borrow), then sign bit of (v - Z - borrow) at each
             * higher limb -- the operands are <= 31 bits so v-Z-borrow
             * lands in [-2^31, 2^31) and bit 31 is the true borrow-out. */
            cc = ct_lt_u32(v0, RCDT_Z[i][0]);
            cc = ((uint32_t)(v1 - RCDT_Z[i][1] - cc)) >> 31;
            cc = ((uint32_t)(v2 - RCDT_Z[i][2] - cc)) >> 31;
            z += (int32_t)cc;
        }
        z_out[s] = z;
    }
    return GAUSS_BATCH;
}
#endif /* !USE_AVX2_SAMPLER */

/* ---- batched uniform 7-bit sampler (AVX2 vpshufb/vpsrlvd unpack) ----
 * Extracts GAUSS_BATCH 7-bit fields, LSB-first, from the bit-packed rand.
 * Closed form for field i: y[i] = (LE32(rand[B_i..]) >> S_i) & 0x7F with
 * B_i=floor(7i/8), S_i=(7i)&7.  Process 8 fields/group (i=8g..8g+7): a
 * 16-byte window at rand+7g (since 7*8g is byte-aligned), broadcast to
 * both 128-lanes; vpshufb gathers each field's 4 source bytes into its
 * 32-bit lane (AVX2 vpshufb is per-128-lane, hence the broadcast), vpsrlvd
 * applies the per-lane shift (AVX2 has vpsrlvd, not the 16-bit vpsrlvw),
 * vpand masks to 7 bits; a second vpshufb packs the 8 low bytes for two
 * 32-bit stores.  Bit-identical to the scalar bit-accumulator (same
 * LSB-first stream); the constants are the closed form, not magic.
 * Requires GAUSS_BATCH==32 && Y_BITS==7 (every mode); the 16-byte
 * window over-reads into the >=28-byte gauss buffer. Non-static: the 4-way
 * sampler reuses it. */
typedef char
    static_assert_y32x7_avx2[(GAUSS_BATCH == 32 && Y_BITS == 7) ? 1 : -1];
void sampler_y(uint8_t y_out[GAUSS_BATCH], const uint8_t *rand)
{
    /* per 32-bit lane l in [0,8): the 4 source bytes [floor(7l/8)+0..3]
     * within the group's 16-byte window (low-128 = l 0..3, high-128 =
     * l 4..7). */
    static const uint8_t WIN_IDX[32] = {0, 1, 2, 3, 0, 1, 2, 3, 1, 2, 3,
                                        4, 2, 3, 4, 5, 3, 4, 5, 6, 4, 5,
                                        6, 7, 5, 6, 7, 8, 6, 7, 8, 9};
    /* gather byte 0 of each u32 lane to the low 4 bytes of each 128-lane.
     */
    static const uint8_t LOW_IDX[32] = {
        0,    4,    8,    12,   0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
        0x80, 0x80, 0x80, 0x80, 0x80, 0,    4,    8,    12,   0x80, 0x80,
        0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80};
    const __m256i win = _mm256_loadu_si256((const __m256i *)WIN_IDX);
    const __m256i low = _mm256_loadu_si256((const __m256i *)LOW_IDX);
    const __m256i sh =
        _mm256_setr_epi32(0, 7, 6, 5, 4, 3, 2, 1); /* (7l)&7 */
    const __m256i m7f = _mm256_set1_epi32(0x7F);
    for (int g = 0; g < GAUSS_BATCH / 8; g++) {
        __m256i wb = _mm256_broadcastsi128_si256(
            _mm_loadu_si128((const __m128i *)(rand + 7 * g)));
        __m256i v = _mm256_and_si256(
            _mm256_srlv_epi32(_mm256_shuffle_epi8(wb, win), sh), m7f);
        __m256i bb = _mm256_shuffle_epi8(v, low);
        uint32_t lo =
            (uint32_t)_mm_cvtsi128_si32(_mm256_castsi256_si128(bb));
        uint32_t hi =
            (uint32_t)_mm_cvtsi128_si32(_mm256_extracti128_si256(bb, 1));
        memcpy(y_out + 8 * g + 0, &lo, 4);
        memcpy(y_out + 8 * g + 4, &hi, 4);
    }
}

/* ---- one BLISS attempt, split into (a) the approx_exp input and (b) the
 * accept/reject + write tail, so the GAUSS_BATCH independent attempts of a
 * mini-batch can have their exp() evaluated 4-way (see sample_gauss_N).
 * ---- */

/* a_q60 = (num * I80) >> 20, num = y*(y+2kx), I80 = floor(2^80/(2 r^2)).
 * Depends only on (x,y), so it can be precomputed for the whole
 * mini-batch. Non-static so the AVX2 4-way sampler can reuse the exact Q60
 * input. */
uint64_t gauss_a_q60(uint32_t x, uint32_t y)
{
    uint32_t t = y + (x << TWO_K_BITS); /* y + 2kx */
    uint64_t num = (uint64_t)y * (uint64_t)t;
    wide_uint128 mwide = (wide_uint128)num * SAMP_I_LO +
                         ((wide_uint128)((uint64_t)num * SAMP_I_HI) << 64);
    return (uint64_t)(mwide >> 20);
}

/* Batched a_q60, 4-wide (AVX2).  a_q60 = (num*I_LO) >> 20, num = y*(y+2kx)
 * < 2^20, I_LO is 64-bit (SAMP_I_HI==0 every mode).  Split I_LO =
 * L1*2^32 + L0: num*I_LO = num*L1*2^32 + num*L0, and 2^32 = 0 (mod 2^20),
 * so a_q60 = (num*L1) << 12  +  (num*L0) >> 20 is EXACT (a 64-bit add;
 * both 32x32->64 products via vpmuludq).  Bit-identical to the scalar
 * 128-bit multiply -> KAT-neutral.  AVX2 has no 52x52 IFMA (the avx512
 * build's path); this u32-limb form is the AVX2 equivalent. */
void gauss_a_q60_batch(uint64_t *out, const int32_t *z, const uint8_t *y,
                       int n)
{
    const __m256i L0 =
        _mm256_set1_epi64x((long long)(SAMP_I_LO & 0xFFFFFFFFu));
    const __m256i L1 = _mm256_set1_epi64x((long long)(SAMP_I_LO >> 32));
    int j = 0;
    for (; j + 4 <= n; j += 4) {
        __m256i xv = _mm256_cvtepi32_epi64(
            _mm_loadu_si128((const __m128i *)(z + j))); /* 4x x (>=0) */
        __m256i yv =
            _mm256_cvtepu8_epi64(_mm_cvtsi32_si128((int)load_le32(y + j)));
        __m256i t =
            _mm256_add_epi64(yv, _mm256_slli_epi64(xv, TWO_K_BITS));
        __m256i num = _mm256_mul_epu32(yv, t); /* y*t < 2^20 */
        __m256i lo = _mm256_mul_epu32(num, L0);
        __m256i hi = _mm256_mul_epu32(num, L1);
        __m256i aq = _mm256_add_epi64(_mm256_slli_epi64(hi, 12),
                                      _mm256_srli_epi64(lo, 20));
        _mm256_storeu_si256((__m256i *)(out + j), aq);
    }
    for (; j < n; j++) /* tail (GAUSS_BATCH=32 is a multiple of 4) */
        out[j] = gauss_a_q60((uint32_t)z[j], (uint32_t)y[j]);
}

/* Accept/reject tail given the precomputed exp_q63 = approx_exp(a_q60).
 * Writes the signed candidate to *r and returns 1 if accepted.  Non-static
 * so the AVX2 4-way sampler reuses the exact accept/zero-fold/sign logic
 * (KAT-critical). */
int gauss_finalize(int32_t *r, uint32_t x, uint32_t y, uint64_t exp_q63,
                   const uint8_t *rej_rand, const uint8_t *signs,
                   size_t idx)
{
    int32_t candidate = ((int32_t)x << K_BITS) | (int32_t)y;
    uint64_t rand_tail = load_le64(rej_rand);
    uint8_t sign_r0 = (uint8_t)(rand_tail & 1U);
    uint64_t rand_rej_63 = rand_tail >> 1;

    int accepted = (int)ct_lt_u64(rand_rej_63, exp_q63);

    /* candidate==0: accept only on one of the two signs (fold). */
    uint32_t r_is_zero = ct_is_zero_u32((uint32_t)candidate);
    accepted &= (int)(1U ^ (r_is_zero & (1U ^ (uint32_t)sign_r0)));

    uint8_t sign = (signs[idx >> 3] >> (idx & 7)) & 1;
    *r = ct_sel_i32(sign, -candidate, candidate);
    return accepted;
}

void sample_gauss_N(int32_t *r, const uint8_t seed[SEEDBYTES],
                    uint64_t nonce, size_t len)
{
    if (len > (size_t)N)
        return; /* defensive API boundary; in-tree callers pass <= N */

    uint8_t buf[GAUSS_BUF_SIZE];
    stream256_state state;
    stream256_init(&state, seed, nonce);

    size_t sign_bytes = (len + 7) / 8;
    size_t init_needed = sign_bytes + MINIBATCH_RAND_BYTES;
    size_t init_nblocks =
        (init_needed + STREAM256_BLOCKBYTES - 1) / STREAM256_BLOCKBYTES;
    {
        PROF_START(ts);
        stream256_squeezeblocks(buf, init_nblocks, &state);
        PROF_STOP(PT_G_SHAKE, ts);
    }

    uint8_t signs[(N + 7) / 8];
    memcpy(signs, buf, sign_bytes);
    PROF_USE(PU_SIGNS, sign_bytes);
    size_t pos = sign_bytes;
    size_t avail = init_nblocks * STREAM256_BLOCKBYTES - sign_bytes;

    size_t coefcnt = 0;
    int32_t z[GAUSS_BATCH];
    uint8_t y[GAUSS_BATCH];

    /* Isochronous Gaussian sampling accepts runtime variance from public
     * rejection progress; the accepted values are protected by branchless
     * finalize logic and output-indexed signs. */
    while (coefcnt < len) {
        if (avail < (size_t)MINIBATCH_RAND_BYTES) {
            if (avail > 0)
                memmove(buf, buf + pos, avail);
            pos = 0;
            PROF_START(ts);
            stream256_squeezeblocks(buf + avail, SQUEEZE_NBLOCKS, &state);
            PROF_STOP(PT_G_SHAKE, ts);
            avail += SQUEEZE_NBLOCKS * STREAM256_BLOCKBYTES;
        }
        PROF_START(tb);
        sampler_sigma2(z, buf + pos);
        PROF_STOP(PT_G_BASESAMP, tb);
        PROF_USE(PU_SIGMA2, SIGMA2_RAND_BYTES);
        pos += SIGMA2_RAND_BYTES;
        avail -= SIGMA2_RAND_BYTES;

        PROF_START(tf);
        sampler_y(y, buf + pos);
        PROF_USE(PU_Y, Y_RAND_BYTES);
        pos += Y_RAND_BYTES;
        avail -= Y_RAND_BYTES;

        /* The GAUSS_BATCH accept tests are independent (each a_q60 depends
         * only on its own (z[j],y[j])), so evaluate exp() as a batch via
         * the shared dispatch (integer approx_exp_4x by default --
         * bit-identical to a per- attempt approx_exp; float AVX2 kernel
         * under -DAVX2_EXP). */
        uint64_t a_q60[GAUSS_BATCH], expv[GAUSS_BATCH];
        gauss_a_q60_batch(a_q60, z, y, GAUSS_BATCH);
        PROF_STOP(PT_G_FINAL, tf);
        PROF_START(te);
        exp_batch(expv, a_q60, GAUSS_BATCH);
        PROF_STOP(PT_G_APPROXEXP, te);

        PROF_START(tf2);
        for (int j = 0; j < GAUSS_BATCH && coefcnt < len; j++) {
            int acc =
                gauss_finalize(&r[coefcnt], (uint32_t)z[j], (uint32_t)y[j],
                               expv[j], buf + pos, signs, coefcnt);
            PROF_USE(PU_GTAIL, GAUSS_RAND_BYTES);
            pos += GAUSS_RAND_BYTES;
            avail -= GAUSS_RAND_BYTES;
            if (acc)
                coefcnt++;
        }
        PROF_STOP(PT_G_FINAL, tf2);
    }
}

/* ---- shared unified combined-stream core (KAT-critical; identical math
 * to sample_gauss_N, only the byte source differs).  Samples ONE poly from
 * the persistent stream gs, refilling one round-robin 4-block unit at a
 * time; the leftover carries to the NEXT poly so the whole y-vector wastes
 * only the final partial unit.  Shared by the ref (4 scalar SHAKE) and
 * avx2 (keccak4x) backends -> bit-identical. ---- */
static void gs_ensure(gauss_stream *gs, size_t need)
{
    while (gs->avail < need) {
        if (gs->avail && gs->pos)
            memmove(gs->buf, gs->buf + gs->pos, gs->avail);
        gs->pos = 0;
        PROF_START(ts);
        gs->squeeze(gs->st, gs->buf + gs->avail);
        PROF_STOP(PT_G_SHAKE, ts);
        gs->avail += GAUSS_COMB_LANES * STREAM256_BLOCKBYTES;
    }
}

void gauss_stream_poly(gauss_stream *gs, int32_t *r)
{
    uint8_t signs[(N + 7) / 8];
    const size_t sign_bytes = (N + 7) / 8;
    size_t coefcnt = 0;
    int32_t z[GAUSS_BATCH];
    uint8_t y[GAUSS_BATCH];

    gs_ensure(gs, sign_bytes);
    memcpy(signs, gs->buf + gs->pos, sign_bytes);
    PROF_USE(PU_SIGNS, sign_bytes);
    gs->pos += sign_bytes;
    gs->avail -= sign_bytes;

    while (coefcnt < (size_t)N) {
        gs_ensure(gs, (size_t)MINIBATCH_RAND_BYTES);
        PROF_START(tb);
        sampler_sigma2(z, gs->buf + gs->pos);
        PROF_STOP(PT_G_BASESAMP, tb);
        PROF_USE(PU_SIGMA2, SIGMA2_RAND_BYTES);
        gs->pos += SIGMA2_RAND_BYTES;
        gs->avail -= SIGMA2_RAND_BYTES;

        PROF_START(tf);
        sampler_y(y, gs->buf + gs->pos);
        PROF_USE(PU_Y, Y_RAND_BYTES);
        gs->pos += Y_RAND_BYTES;
        gs->avail -= Y_RAND_BYTES;
        uint64_t a_q60[GAUSS_BATCH], expv[GAUSS_BATCH];
        gauss_a_q60_batch(a_q60, z, y, GAUSS_BATCH);
        PROF_STOP(PT_G_FINAL, tf);
        PROF_START(te);
        exp_batch(expv, a_q60, GAUSS_BATCH);
        PROF_STOP(PT_G_APPROXEXP, te);

        PROF_START(tf2);
        for (int j = 0; j < GAUSS_BATCH && coefcnt < (size_t)N; j++) {
            int acc =
                gauss_finalize(&r[coefcnt], (uint32_t)z[j], (uint32_t)y[j],
                               expv[j], gs->buf + gs->pos, signs, coefcnt);
            PROF_USE(PU_GTAIL, GAUSS_RAND_BYTES);
            gs->pos += GAUSS_RAND_BYTES;
            gs->avail -= GAUSS_RAND_BYTES;
            if (acc)
                coefcnt++;
        }
        PROF_STOP(PT_G_FINAL, tf2);
    }
}

/* 4-way batch wrappers.  These sequential versions are the canonical KAT
 * definition (4 independent (seed,nonce) SHAKE streams == 4 single calls);
 * the AVX2 build replaces them with the keccak4x versions below
 * that produce byte-identical per-lane output, so they are guarded out
 * here when USE_AVX2_SHAKE4X is defined. */
#ifndef USE_AVX2_SHAKE4X
void sample_gauss_N_4x(int32_t *r0, int32_t *r1, int32_t *r2, int32_t *r3,
                       const uint8_t seed[SEEDBYTES], uint64_t nonce0,
                       uint64_t nonce1, uint64_t nonce2, uint64_t nonce3,
                       size_t len0, size_t len1, size_t len2, size_t len3)
{
    if (len0)
        sample_gauss_N(r0, seed, nonce0, len0);
    if (len1)
        sample_gauss_N(r1, seed, nonce1, len1);
    if (len2)
        sample_gauss_N(r2, seed, nonce2, len2);
    if (len3)
        sample_gauss_N(r3, seed, nonce3, len3);
}

void sample_noise_N_4x(int32_t *r0, int32_t *r1, int32_t *r2, int32_t *r3,
                       const uint8_t seed[SEEDBYTES], uint64_t nonce0,
                       uint64_t nonce1, uint64_t nonce2, uint64_t nonce3,
                       size_t len0, size_t len1, size_t len2, size_t len3)
{
    if (len0)
        sample_noise_N(r0, seed, nonce0, len0);
    if (len1)
        sample_noise_N(r1, seed, nonce1, len1);
    if (len2)
        sample_noise_N(r2, seed, nonce2, len2);
    if (len3)
        sample_noise_N(r3, seed, nonce3, len3);
}

/* ref combined backend: GAUSS_COMB_LANES scalar SHAKE-256 sub-streams
 * (seed||nonce||j), one block each per round-robin unit.  The avx2 build
 * replaces this with 2x keccak4x below, producing byte-identical
 * units; this scalar copy is excluded under USE_AVX2_SHAKE4X. */
static void comb_ref_squeeze8(void *v, uint8_t *out)
{
    keccak_state *s = (keccak_state *)v;
    for (int j = 0; j < GAUSS_COMB_LANES; j++)
        shake256_squeezeblocks(out + (size_t)j * STREAM256_BLOCKBYTES, 1,
                               &s[j]);
}

void sample_gauss_vec(int32_t *const r[], int cnt,
                      const uint8_t seed[SKSEEDBYTES], uint64_t nonce)
{
    keccak_state s[GAUSS_COMB_LANES];
    uint8_t in[SKSEEDBYTES + 9];
    memcpy(in, seed, SKSEEDBYTES);
    for (int i = 0; i < 8; i++)
        in[SKSEEDBYTES + i] = (uint8_t)(nonce >> (8 * i));
    for (int j = 0; j < GAUSS_COMB_LANES; j++) {
        in[SKSEEDBYTES + 8] = (uint8_t)j;
        shake256_absorb_once(&s[j], in, SKSEEDBYTES + 9);
    }
    gauss_stream gs;
    gs.st = s;
    gs.squeeze = comb_ref_squeeze8;
    gs.pos = 0;
    gs.avail = 0;
    for (int i = 0; i < cnt; i++)
        gauss_stream_poly(&gs, r[i]);
}

void sample_noise_vec(int32_t *const r[], int cnt,
                      const uint8_t seed[SKSEEDBYTES], uint64_t nonce)
{
    keccak_state s[GAUSS_COMB_LANES];
    uint8_t in[SKSEEDBYTES + 9];
    memcpy(in, seed, SKSEEDBYTES);
    for (int i = 0; i < 8; i++)
        in[SKSEEDBYTES + i] = (uint8_t)(nonce >> (8 * i));
    for (int j = 0; j < GAUSS_COMB_LANES; j++) {
        in[SKSEEDBYTES + 8] = (uint8_t)j;
        shake256_absorb_once(&s[j], in, SKSEEDBYTES + 9);
    }
    gauss_stream gs;
    gs.st = s;
    gs.squeeze = comb_ref_squeeze8;
    gs.pos = 0;
    gs.avail = 0;
    for (int i = 0; i < cnt; i++)
        noise_stream_poly(&gs, r[i]);
}
#endif /* !USE_AVX2_SHAKE4X */

/* ---- unified key-noise vector core (grouped, SIMD-friendly) -- see ref/
 * sampler.c.  The scalar noise_magnitude_batch is guarded out here when
 * USE_AVX2_SAMPLER is defined: the AVX2 build supplies the 8-way
 * SIMD version below, mirroring how sampler_sigma2 is replaced. ----
 */
#ifndef USE_AVX2_SAMPLER
void noise_magnitude_batch(int32_t *m_out, const uint8_t *rand)
{
    for (int s = 0; s < NOISE_BATCH; s++) {
        int group = s >> 3, lane = s & 7;
        int byte_base = group * 96 + lane * 4;
        uint32_t v0 = load_le32(rand + byte_base) & 0x7FFFFFFFU;
        uint32_t v1 = load_le32(rand + byte_base + 32) & 0x7FFFFFFFU;
        uint32_t v2 = load_le32(rand + byte_base + 64) & 0x7FFFFFFFU;
        int32_t m = 0;
#    if defined(__clang__)
#        pragma clang loop vectorize(disable) interleave(disable)
#    endif
        for (int i = 0; i < NOISE_RCDT_ENTRIES; i++) {
            uint32_t cc;
            cc = ct_lt_u32(v0, RCDT_NOISE[i][0]);
            cc = ((uint32_t)(v1 - RCDT_NOISE[i][1] - cc)) >> 31;
            cc = ((uint32_t)(v2 - RCDT_NOISE[i][2] - cc)) >> 31;
            m += (int32_t)cc;
        }
        m_out[s] = m;
    }
}
#endif /* !USE_AVX2_SAMPLER */

void noise_stream_poly(gauss_stream *gs, int32_t *r)
{
    size_t cnt = 0;
    int32_t m[NOISE_BATCH];
    while (cnt < (size_t)N) {
        gs_ensure(gs, (size_t)NOISE_MINIBATCH_RAND_BYTES);
        PROF_START(tb);
        noise_magnitude_batch(m, gs->buf + gs->pos);
        PROF_STOP(PT_G_BASESAMP, tb);
        PROF_USE(PU_SIGMA2, NOISE_CDT_BYTES);
        gs->pos += NOISE_CDT_BYTES;
        gs->avail -= NOISE_CDT_BYTES;
        /* One tail byte packs 4 candidates (2 bits each: bit0 sign, bit1
         * zero-fold), LSB-first.  Consume the whole NOISE_BATCH/4 tail
         * block up front so the cursor advance is independent of the early
         * cnt==N break -- deterministic and identical across backends. */
        const uint8_t *tailp = gs->buf + gs->pos;
        gs->pos += NOISE_BATCH / 4;
        gs->avail -= NOISE_BATCH / 4;
        PROF_USE(PU_GTAIL, NOISE_BATCH / 4);
        for (int j = 0; j < NOISE_BATCH && cnt < (size_t)N; j++) {
            uint8_t f = (uint8_t)((tailp[j >> 2] >> (2 * (j & 3))) & 3u);
            uint8_t sign = f & 1;
            uint8_t zero_fold = (f >> 1) & 1;
            uint32_t reject =
                ct_is_zero_u32((uint32_t)m[j]) & (uint32_t)zero_fold;
            int32_t v = m[j];
            r[cnt] = ct_sel_i32(sign, -v, v);
            cnt += (size_t)(1U ^ reject);
        }
    }
}

/* ============================================================
 * sigma = 1 noise sampler (keygen): magnitude via 10-entry RCDT, sign,
 * fold 0. Per coefficient: 12 bytes CDT (3x31) + 1 byte (bit0 sign, bit1
 * zero-fold).
 * ============================================================ */
/* Non-static so the AVX2 4-way noise sampler reuses the exact RCDT scan.
 */
int noise_magnitude(const uint8_t *rand12)
{
    uint32_t v0 = load_le32(rand12) & 0x7FFFFFFFU;
    uint32_t v1 = load_le32(rand12 + 4) & 0x7FFFFFFFU;
    uint32_t v2 = load_le32(rand12 + 8) & 0x7FFFFFFFU;
    int32_t m = 0;
#if defined(__clang__)
#    pragma clang loop vectorize(disable) interleave(disable)
#endif
    for (int i = 0; i < NOISE_RCDT_ENTRIES; i++) {
        uint32_t cc;
        /* exact 31-bit-limb borrow chain (see sampler_sigma2). */
        cc = ct_lt_u32(v0, RCDT_NOISE[i][0]);
        cc = ((uint32_t)(v1 - RCDT_NOISE[i][1] - cc)) >> 31;
        cc = ((uint32_t)(v2 - RCDT_NOISE[i][2] - cc)) >> 31;
        m += (int32_t)cc;
    }
    return m;
}

void sample_noise_N(int32_t *r, const uint8_t seed[SEEDBYTES],
                    uint64_t nonce, size_t len)
{
    /* 13 bytes/attempt; squeeze NOISE_NBLOCKS blocks per refill. Buffer
     * holds the (<need) residual plus one full squeeze. */
    enum {
        NOISE_NBLOCKS = 4,
        NOISE_BUF = NOISE_NBLOCKS * STREAM256_BLOCKBYTES + 16
    };
    uint8_t buf[NOISE_BUF];
    stream256_state state;
    stream256_init(&state, seed, nonce);
    size_t pos = 0, avail = 0;
    size_t cnt = 0;
    const size_t need = 13;
    while (cnt < len) {
        if (avail < need) {
            if (avail > 0)
                memmove(buf, buf + pos, avail);
            pos = 0;
            stream256_squeezeblocks(buf + avail, NOISE_NBLOCKS, &state);
            avail += NOISE_NBLOCKS * STREAM256_BLOCKBYTES;
        }
        int m = noise_magnitude(buf + pos);
        uint8_t tail = buf[pos + 12];
        pos += need;
        avail -= need;
        uint8_t sign = tail & 1;
        uint8_t zero_fold = (tail >> 1) & 1;
        /* Zero-fold reject WITHOUT branching on the secret magnitude m (m
         * is a key-noise coefficient): build a mask, always write the
         * candidate, and advance only when accepted -- same shape as
         * gauss_finalize.  (The loop iteration count still varies, but
         * that is the usual rejection-sampling timing on the public SHAKE
         * stream, not a branch on the secret value.) */
        uint32_t reject =
            ct_is_zero_u32((uint32_t)m) & (uint32_t)zero_fold;
        int32_t v = m;
        r[cnt] = ct_sel_i32(sign, -v, v);
        cnt += (size_t)(1U ^ reject); /* accepted -> advance */
    }
}

/* ==========================================================================
 * AVX2 vectorized samplers.  This block provides the 8-way sigma2
 * sampler (sampler_sigma2), the 8-way key-noise CDT
 * (noise_magnitude_batch), keccak4x 4-way Gaussian/noise wrappers
 * (sample_gauss_N_4x / sample_noise_N_4x), and the unified
 * 8-substream y/noise vector backends (sample_gauss_vec /
 * sample_noise_vec). The scalar counterparts above are guarded out under
 * USE_AVX2_*; all symbols are prototyped in sampler.h.
 * ==========================================================================
 */
/* GAUSS_BATCH must be a whole number of 16-sample (2x8) units.  With the
 * canonical GAUSS_BATCH=32 the loop runs twice -- the user's "wrap an
 * outer loop, no extra unroll": the 16-sample (two-groups-of-8) core is
 * unchanged, just applied to each 192-byte half of the mini-batch. */
typedef char static_assert_batch16x_avx[(GAUSS_BATCH % 16 == 0) ? 1 : -1];

int sampler_sigma2(int32_t *z_out, const uint8_t *rand)
{
    const __m256i mask31 = _mm256_set1_epi32(0x7FFFFFFF);

    for (int t = 0; t < GAUSS_BATCH / 16; t++) {
        const uint8_t *rb =
            rand + (size_t)t * 192; /* 16 samples = 2x96 B */
        int32_t *zo = z_out + t * 16;

        /* group a = samples 0..7 (bytes 0..95), group b = samples 8..15
         * (96..191) */
        __m256i v0a = _mm256_and_si256(
            _mm256_loadu_si256((const __m256i *)(rb + 0)), mask31);
        __m256i v1a = _mm256_and_si256(
            _mm256_loadu_si256((const __m256i *)(rb + 32)), mask31);
        __m256i v2a = _mm256_and_si256(
            _mm256_loadu_si256((const __m256i *)(rb + 64)), mask31);
        __m256i v0b = _mm256_and_si256(
            _mm256_loadu_si256((const __m256i *)(rb + 96)), mask31);
        __m256i v1b = _mm256_and_si256(
            _mm256_loadu_si256((const __m256i *)(rb + 128)), mask31);
        __m256i v2b = _mm256_and_si256(
            _mm256_loadu_si256((const __m256i *)(rb + 160)), mask31);

        __m256i za = _mm256_setzero_si256();
        __m256i zb = _mm256_setzero_si256();

#if defined(__clang__)
#    pragma clang loop vectorize(disable) interleave(disable)
#endif
        for (int i = 0; i < RCDT_ENTRIES; i++) {
            __m256i Z0 = _mm256_set1_epi32((int)RCDT_Z[i][0]);
            __m256i Z1 = _mm256_set1_epi32((int)RCDT_Z[i][1]);
            __m256i Z2 = _mm256_set1_epi32((int)RCDT_Z[i][2]);

            /* exact 31-bit-limb borrow chain, two independent 8-lane
             * groups for ILP (dual accumulator).  limb0: cmpgt(Z0,v0) =
             * -borrow0 in {0,-1}; higher limbs: borrow = bit 31 of (v - Z
             * - borrow_in). */
            __m256i ca = _mm256_cmpgt_epi32(Z0, v0a); /* -borrow0 */
            __m256i cb = _mm256_cmpgt_epi32(Z0, v0b);
            /* borrow1 = (v1 + (-borrow0) - Z1) >> 31 */
            ca = _mm256_srli_epi32(
                _mm256_sub_epi32(_mm256_add_epi32(v1a, ca), Z1), 31);
            cb = _mm256_srli_epi32(
                _mm256_sub_epi32(_mm256_add_epi32(v1b, cb), Z1), 31);
            /* borrow2 = (v2 - borrow1 - Z2) >> 31 */
            ca = _mm256_srli_epi32(
                _mm256_sub_epi32(_mm256_sub_epi32(v2a, ca), Z2), 31);
            cb = _mm256_srli_epi32(
                _mm256_sub_epi32(_mm256_sub_epi32(v2b, cb), Z2), 31);
            za = _mm256_add_epi32(za, ca);
            zb = _mm256_add_epi32(zb, cb);
        }

        /* za = samples t*16+0..7, zb = +8..15: store the 32-bit lanes
         * directly in sample order (int32 output needs no pack/permute --
         * the old packs+permute4x64 only undid the int16 in-lane
         * interleave). */
        _mm256_storeu_si256((__m256i *)zo, za);
        _mm256_storeu_si256((__m256i *)(zo + 8), zb);
    }
    return GAUSS_BATCH;
}

/* AVX2 8-way key-noise magnitude scan (noise_magnitude_batch) -- the same
 * 8-way borrow chain as sampler_sigma2 above, just over the 10-entry
 * RCDT_NOISE table and NOISE_BATCH samples (magnitude only, no y).
 * Replaces the scalar noise_magnitude_batch (guarded out of sampler.c by
 * USE_AVX2_SAMPLER); reads the identical grouped layout (group*96
 * + lane*4; v0/v1/v2 at +0/+32/+64), so it is bit-identical to the scalar
 * borrow chain -> KAT-neutral.  Vectorizes the previously-scalar keygen
 * noise sampler (cold-but-now-faster path). */
typedef char static_assert_nbatch16x_avx[(NOISE_BATCH % 16 == 0) ? 1 : -1];

void noise_magnitude_batch(int32_t *m_out, const uint8_t *rand)
{
    const __m256i mask31 = _mm256_set1_epi32(0x7FFFFFFF);

    for (int t = 0; t < NOISE_BATCH / 16; t++) {
        const uint8_t *rb =
            rand + (size_t)t * 192; /* 16 samples = 2x96 B */
        int32_t *mo = m_out + t * 16;

        __m256i v0a = _mm256_and_si256(
            _mm256_loadu_si256((const __m256i *)(rb + 0)), mask31);
        __m256i v1a = _mm256_and_si256(
            _mm256_loadu_si256((const __m256i *)(rb + 32)), mask31);
        __m256i v2a = _mm256_and_si256(
            _mm256_loadu_si256((const __m256i *)(rb + 64)), mask31);
        __m256i v0b = _mm256_and_si256(
            _mm256_loadu_si256((const __m256i *)(rb + 96)), mask31);
        __m256i v1b = _mm256_and_si256(
            _mm256_loadu_si256((const __m256i *)(rb + 128)), mask31);
        __m256i v2b = _mm256_and_si256(
            _mm256_loadu_si256((const __m256i *)(rb + 160)), mask31);

        __m256i za = _mm256_setzero_si256();
        __m256i zb = _mm256_setzero_si256();

#if defined(__clang__)
#    pragma clang loop vectorize(disable) interleave(disable)
#endif
        for (int i = 0; i < NOISE_RCDT_ENTRIES; i++) {
            __m256i Z0 = _mm256_set1_epi32((int)RCDT_NOISE[i][0]);
            __m256i Z1 = _mm256_set1_epi32((int)RCDT_NOISE[i][1]);
            __m256i Z2 = _mm256_set1_epi32((int)RCDT_NOISE[i][2]);

            __m256i ca = _mm256_cmpgt_epi32(Z0, v0a); /* -borrow0 */
            __m256i cb = _mm256_cmpgt_epi32(Z0, v0b);
            ca = _mm256_srli_epi32(
                _mm256_sub_epi32(_mm256_add_epi32(v1a, ca), Z1), 31);
            cb = _mm256_srli_epi32(
                _mm256_sub_epi32(_mm256_add_epi32(v1b, cb), Z1), 31);
            ca = _mm256_srli_epi32(
                _mm256_sub_epi32(_mm256_sub_epi32(v2a, ca), Z2), 31);
            cb = _mm256_srli_epi32(
                _mm256_sub_epi32(_mm256_sub_epi32(v2b, cb), Z2), 31);
            za = _mm256_add_epi32(za, ca);
            zb = _mm256_add_epi32(zb, cb);
        }

        /* int32 output: store za (samples t*16+0..7) and zb (+8..15)
         * directly; no pack/permute. */
        _mm256_storeu_si256((__m256i *)mo, za);
        _mm256_storeu_si256((__m256i *)(mo + 8), zb);
    }
}
#define BLK ((size_t)STREAM256_BLOCKBYTES) /* 136 */

/* ---- Gaussian (D_{Z,r}) ----
 * GAUSS4X_INIT sizes the initial 4-way squeeze so the common case needs no
 * continuation squeeze; the refill path still runs (and is verified
 * KAT-stable) when a rejection-heavy lane needs more.  Overridable
 * (-DGAUSS4X_INIT=n) so a tiny value can force-exercise the refill path in
 * testing.  The buffer must hold BOTH the initial fill (INIT*BLK) AND the
 * worst-case refill occupancy: a refill memmoves `avail` (< MINIBATCH, the
 * trigger) to the front then squeezes SQUEEZE_NBLOCKS*BLK more, so
 * occupancy stays < MINIBATCH_RAND_BYTES + SQUEEZE_NBLOCKS*BLK.  Hence the
 * buffer is sized to max(INIT*BLK, MINIBATCH) + SQUEEZE_NBLOCKS*BLK (+
 * margin); the direct (INIT+SQUEEZE)*BLK form under-sizes it once
 * GAUSS_BATCH grew 16->32 pushed MINIBATCH past INIT*BLK for tiny test
 * INITs.
 */
#ifndef GAUSS4X_INIT
/* Per-mode initial 4-way squeeze, sized just above the measured per-poly
 * consumption (mode120≈51, mode180≈62, mode260≈95 blocks/poly) so the
 * common case needs no continuation squeeze; the refill path (verified by
 * `make check-refill`) handles rejection-heavy outliers. */
#    if LITHIUM_MODE == 120
#        define GAUSS4X_INIT 56
#    elif LITHIUM_MODE == 180
#        define GAUSS4X_INIT 68
#    else /* 260 */
#        define GAUSS4X_INIT 104
#    endif
#endif
#define GAUSS4X_INIT_BYTES ((size_t)GAUSS4X_INIT * BLK)
#define GAUSS4X_BUF                                     \
    ((GAUSS4X_INIT_BYTES > (size_t)MINIBATCH_RAND_BYTES \
          ? GAUSS4X_INIT_BYTES                          \
          : (size_t)MINIBATCH_RAND_BYTES) +             \
     (size_t)SQUEEZE_NBLOCKS * BLK + 32)

void sample_gauss_N_4x(int32_t *r0, int32_t *r1, int32_t *r2, int32_t *r3,
                       const uint8_t seed[SEEDBYTES], uint64_t nonce0,
                       uint64_t nonce1, uint64_t nonce2, uint64_t nonce3,
                       size_t len0, size_t len1, size_t len2, size_t len3)
{
    if (len0 > (size_t)N || len1 > (size_t)N || len2 > (size_t)N ||
        len3 > (size_t)N)
        return; /* defensive API boundary; in-tree callers pass <= N */

    int32_t *r[4] = {r0, r1, r2, r3};
    uint64_t non[4] = {nonce0, nonce1, nonce2, nonce3};
    size_t len[4] = {len0, len1, len2, len3};

    uint8_t inbuf[4][SEEDBYTES + 8];
    for (int j = 0; j < 4; j++) {
        memcpy(inbuf[j], seed, SEEDBYTES);
        for (int i = 0; i < 8; i++)
            inbuf[j][SEEDBYTES + i] = (uint8_t)(non[j] >> (8 * i));
    }

    keccakx4_state st;
    PROF_START(ts0);
    shake256x4_absorb_once(&st, inbuf[0], inbuf[1], inbuf[2], inbuf[3],
                           SEEDBYTES + 8);

    uint8_t buf[4][GAUSS4X_BUF];
    shake256x4_squeezeblocks(buf[0], buf[1], buf[2], buf[3], GAUSS4X_INIT,
                             &st);
    PROF_STOP(PT_G_SHAKE, ts0);

    uint8_t signs[4][(N + 7) / 8];
    size_t pos[4], avail[4], cnt[4];
    for (int j = 0; j < 4; j++) {
        size_t sb = (len[j] + 7) / 8; /* matches ref sample_gauss_N */
        if (len[j]) {
            memcpy(signs[j], buf[j], sb);
            PROF_USE(PU_SIGNS, sb);
        }
        pos[j] = sb;
        avail[j] = (size_t)GAUSS4X_INIT * BLK - sb;
        cnt[j] = 0;
    }

    while (cnt[0] < len[0] || cnt[1] < len[1] || cnt[2] < len[2] ||
           cnt[3] < len[3]) {
        int need = 0;
        for (int j = 0; j < 4; j++)
            if (cnt[j] < len[j] && avail[j] < (size_t)MINIBATCH_RAND_BYTES)
                need = 1;
        if (need) {
            for (int j = 0; j < 4; j++) {
                if (cnt[j] >= len[j]) {
                    pos[j] = 0;
                    avail[j] = 0;
                    continue;
                }
                if (avail[j] && pos[j])
                    memmove(buf[j], buf[j] + pos[j], avail[j]);
                pos[j] = 0;
            }
            PROF_START(ts);
            shake256x4_squeezeblocks(buf[0] + avail[0], buf[1] + avail[1],
                                     buf[2] + avail[2], buf[3] + avail[3],
                                     SQUEEZE_NBLOCKS, &st);
            PROF_STOP(PT_G_SHAKE, ts);
            for (int j = 0; j < 4; j++)
                avail[j] += (size_t)SQUEEZE_NBLOCKS * BLK;
        }

        for (int j = 0; j < 4; j++) {
            if (cnt[j] >= len[j])
                continue;
            if (avail[j] < (size_t)MINIBATCH_RAND_BYTES)
                continue;

            int32_t z[GAUSS_BATCH];
            uint8_t y[GAUSS_BATCH];
            PROF_START(tb);
            sampler_sigma2(z, buf[j] + pos[j]);
            PROF_STOP(PT_G_BASESAMP, tb);
            PROF_USE(PU_SIGMA2, SIGMA2_RAND_BYTES);
            pos[j] += SIGMA2_RAND_BYTES;
            avail[j] -= SIGMA2_RAND_BYTES;

            PROF_START(tf);
            sampler_y(y, buf[j] + pos[j]);
            PROF_USE(PU_Y, Y_RAND_BYTES);
            pos[j] += Y_RAND_BYTES;
            avail[j] -= Y_RAND_BYTES;

            uint64_t a_q60[GAUSS_BATCH], expv[GAUSS_BATCH];
            gauss_a_q60_batch(a_q60, z, y, GAUSS_BATCH);
            PROF_STOP(PT_G_FINAL, tf);
            PROF_START(te);
            exp_batch(expv, a_q60, GAUSS_BATCH);
            PROF_STOP(PT_G_APPROXEXP, te);

            PROF_START(tf2);
            for (int k = 0; k < GAUSS_BATCH && cnt[j] < len[j]; k++) {
                int acc = gauss_finalize(
                    &r[j][cnt[j]], (uint32_t)z[k], (uint32_t)y[k], expv[k],
                    buf[j] + pos[j], signs[j], cnt[j]);
                PROF_USE(PU_GTAIL, GAUSS_RAND_BYTES);
                pos[j] += GAUSS_RAND_BYTES;
                avail[j] -= GAUSS_RAND_BYTES;
                if (acc)
                    cnt[j]++;
            }
            PROF_STOP(PT_G_FINAL, tf2);
        }
    }
}

/* ---- noise (D_{Z,1}) ---- */
#define NOISE_NEED 13     /* 12 RCDT + 1 tail, per ref */
#define NOISE4X_NBLOCKS 4 /* matches ref sample_noise_N refill */
#ifndef NOISE4X_INIT
#    if N <= 256
#        define NOISE4X_INIT 32
#    else
#        define NOISE4X_INIT 60
#    endif
#endif
#define NOISE4X_BUF ((size_t)(NOISE4X_INIT + NOISE4X_NBLOCKS) * BLK + 16)

void sample_noise_N_4x(int32_t *r0, int32_t *r1, int32_t *r2, int32_t *r3,
                       const uint8_t seed[SEEDBYTES], uint64_t nonce0,
                       uint64_t nonce1, uint64_t nonce2, uint64_t nonce3,
                       size_t len0, size_t len1, size_t len2, size_t len3)
{
    int32_t *r[4] = {r0, r1, r2, r3};
    uint64_t non[4] = {nonce0, nonce1, nonce2, nonce3};
    size_t len[4] = {len0, len1, len2, len3};

    uint8_t inbuf[4][SEEDBYTES + 8];
    for (int j = 0; j < 4; j++) {
        memcpy(inbuf[j], seed, SEEDBYTES);
        for (int i = 0; i < 8; i++)
            inbuf[j][SEEDBYTES + i] = (uint8_t)(non[j] >> (8 * i));
    }

    keccakx4_state st;
    shake256x4_absorb_once(&st, inbuf[0], inbuf[1], inbuf[2], inbuf[3],
                           SEEDBYTES + 8);

    uint8_t buf[4][NOISE4X_BUF];
    shake256x4_squeezeblocks(buf[0], buf[1], buf[2], buf[3], NOISE4X_INIT,
                             &st);

    size_t pos[4], avail[4], cnt[4];
    for (int j = 0; j < 4; j++) {
        pos[j] = 0;
        avail[j] = (size_t)NOISE4X_INIT * BLK;
        cnt[j] = 0;
    }

    while (cnt[0] < len[0] || cnt[1] < len[1] || cnt[2] < len[2] ||
           cnt[3] < len[3]) {
        int need = 0;
        for (int j = 0; j < 4; j++)
            if (cnt[j] < len[j] && avail[j] < (size_t)NOISE_NEED)
                need = 1;
        if (need) {
            for (int j = 0; j < 4; j++) {
                if (cnt[j] >= len[j]) {
                    pos[j] = 0;
                    avail[j] = 0;
                    continue;
                }
                if (avail[j] && pos[j])
                    memmove(buf[j], buf[j] + pos[j], avail[j]);
                pos[j] = 0;
            }
            shake256x4_squeezeblocks(buf[0] + avail[0], buf[1] + avail[1],
                                     buf[2] + avail[2], buf[3] + avail[3],
                                     NOISE4X_NBLOCKS, &st);
            for (int j = 0; j < 4; j++)
                avail[j] += (size_t)NOISE4X_NBLOCKS * BLK;
        }

        for (int j = 0; j < 4; j++) {
            if (cnt[j] >= len[j])
                continue;
            if (avail[j] < (size_t)NOISE_NEED)
                continue;

            int m = noise_magnitude(buf[j] + pos[j]);
            uint8_t tail = buf[j][pos[j] + 12];
            pos[j] += NOISE_NEED;
            avail[j] -= NOISE_NEED;
            uint8_t sign = tail & 1;
            uint8_t zero_fold = (tail >> 1) & 1;
            uint32_t reject =
                ct_is_zero_u32((uint32_t)m) & (uint32_t)zero_fold;
            int32_t v = m;
            r[j][cnt[j]] = ct_sel_i32(sign, -v, v);
            cnt[j] += (size_t)(1U ^ reject);
        }
    }
}

/* ---- unified combined GAUSS_COMB_LANES-lane Gaussian y-vector, keccak4x
 * backend ---- GAUSS_COMB_LANES (=8) domain-separated sub-streams
 * SHAKE256(seed||nonce||j) squeezed one block each per round-robin unit,
 * emitted via TWO keccak4x sponges (lanes 0-3, 4-7).  The unit
 * [lane0][lane1]...[lane7] is byte-identical to the ref 8-scalar-SHAKE and
 * the avx512 1x-keccak8x backends; the shared core (sampler.c
 * gauss_stream_poly) samples all cnt polys from one persistent stream.
 * avx2 stays perf-neutral: 2 keccak4x per 8-block unit == same per-block
 * cost as the old 1 keccak4x per 4-block unit. */
typedef struct {
    keccakx4_state lo, hi;
} comb8_avx_state;

static void comb_avx_squeeze8(void *v, uint8_t *out)
{
    comb8_avx_state *c = (comb8_avx_state *)v;
    shake256x4_squeezeblocks(out + 0 * BLK, out + 1 * BLK, out + 2 * BLK,
                             out + 3 * BLK, 1, &c->lo);
    shake256x4_squeezeblocks(out + 4 * BLK, out + 5 * BLK, out + 6 * BLK,
                             out + 7 * BLK, 1, &c->hi);
}

void sample_gauss_vec(int32_t *const r[], int cnt,
                      const uint8_t seed[SKSEEDBYTES], uint64_t nonce)
{
    uint8_t in[GAUSS_COMB_LANES][SKSEEDBYTES + 9];
    for (int j = 0; j < GAUSS_COMB_LANES; j++) {
        memcpy(in[j], seed, SKSEEDBYTES);
        for (int i = 0; i < 8; i++)
            in[j][SKSEEDBYTES + i] = (uint8_t)(nonce >> (8 * i));
        in[j][SKSEEDBYTES + 8] = (uint8_t)j;
    }
    comb8_avx_state c;
    shake256x4_absorb_once(&c.lo, in[0], in[1], in[2], in[3],
                           SKSEEDBYTES + 9);
    shake256x4_absorb_once(&c.hi, in[4], in[5], in[6], in[7],
                           SKSEEDBYTES + 9);
    gauss_stream gs;
    gs.st = &c;
    gs.squeeze = comb_avx_squeeze8;
    gs.pos = 0;
    gs.avail = 0;
    for (int i = 0; i < cnt; i++)
        gauss_stream_poly(&gs, r[i]);
}

/* unified key-noise vector, same 8-substream (2x keccak4x) backend as the
 * Gaussian; the shared core (sampler.c noise_stream_poly) carves all cnt
 * noise polys from one persistent stream -> byte-identical to ref. */
void sample_noise_vec(int32_t *const r[], int cnt,
                      const uint8_t seed[SKSEEDBYTES], uint64_t nonce)
{
    uint8_t in[GAUSS_COMB_LANES][SKSEEDBYTES + 9];
    for (int j = 0; j < GAUSS_COMB_LANES; j++) {
        memcpy(in[j], seed, SKSEEDBYTES);
        for (int i = 0; i < 8; i++)
            in[j][SKSEEDBYTES + i] = (uint8_t)(nonce >> (8 * i));
        in[j][SKSEEDBYTES + 8] = (uint8_t)j;
    }
    comb8_avx_state c;
    shake256x4_absorb_once(&c.lo, in[0], in[1], in[2], in[3],
                           SKSEEDBYTES + 9);
    shake256x4_absorb_once(&c.hi, in[4], in[5], in[6], in[7],
                           SKSEEDBYTES + 9);
    gauss_stream gs;
    gs.st = &c;
    gs.squeeze = comb_avx_squeeze8;
    gs.pos = 0;
    gs.avail = 0;
    for (int i = 0; i < cnt; i++)
        noise_stream_poly(&gs, r[i]);
}
