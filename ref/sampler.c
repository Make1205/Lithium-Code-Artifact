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
 * ref and optimized backends (KAT-compatible).  The GAUSS_BATCH accept
 * tests in a mini-batch are independent, so their exp() values are
 * evaluated 4-way via approx_exp_4x (instruction-level parallelism) --
 * bit-identical to a per- attempt approx_exp, so sample_gauss_N_4x stays
 * KAT-locked to sample_gauss_N.  The Gaussian sampler follows the
 * isochronous Gaussian-sampling discipline: runtime varies with rejection
 * progress, but not in a way that exposes attacker-usable information.
 */
#include "sampler.h"

#include <string.h>

#include "approx_exp.h"
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

/* ---- BLISS base CDT (sigma2 = r/128), GAUSS_BATCH magnitudes ---- */
int sampler_sigma2(int32_t *z_out, const uint8_t *rand)
{
    for (int s = 0; s < GAUSS_BATCH; s++) {
        int group = s >> 3, lane = s & 7;
        int byte_base = group * 96 + lane * 4;
        uint32_t v0 = load_le32(rand + byte_base) & 0x7FFFFFFFU;
        uint32_t v1 = load_le32(rand + byte_base + 32) & 0x7FFFFFFFU;
        uint32_t v2 = load_le32(rand + byte_base + 64) & 0x7FFFFFFFU;
        int32_t z = 0;
#if defined(__clang__)
#    pragma clang loop vectorize(disable) interleave(disable)
#endif
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

/* ---- batched uniform Y_BITS-bit sampler ----
 * Non-static so optimized samplers can reuse the exact bit-packing
 * (KAT-critical). */
void sampler_y(uint8_t y_out[GAUSS_BATCH], const uint8_t *rand)
{
    const uint32_t y_mask = (1U << Y_BITS) - 1U;
    uint32_t acc = 0;
    unsigned acc_bits = 0, byte_idx = 0;
    for (int i = 0; i < GAUSS_BATCH; ++i) {
        while (acc_bits < (unsigned)Y_BITS) {
            acc |= (uint32_t)rand[byte_idx++] << acc_bits;
            acc_bits += 8;
        }
        y_out[i] = (uint8_t)(acc & y_mask);
        acc >>= Y_BITS;
        acc_bits -= Y_BITS;
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

/* Batched a_q60 (scalar; the AVX-512 build replaces this with an 8-lane
 * IFMA kernel in sampler_avx512.c). */
void gauss_a_q60_batch(uint64_t *out, const int32_t *z, const uint8_t *y,
                       int n)
{
    for (int j = 0; j < n; j++)
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
 * definition (4 independent (seed,nonce) SHAKE streams == 4 single calls).
 */
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
 * replaces this with 2x keccak4x and the avx512 build with 1x keccak8x,
 * both producing byte-identical units. */
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

/* ---- unified key-noise vector (no-waste, SIMD-friendly grouped layout)
 * ---- Mirrors the Gaussian unified stream; see sampler.h.
 * noise_magnitude_batch is the grouped (group*96+lane*4) 10-entry RCDT
 * scan -- structurally identical to the scalar sampler_sigma2, just the
 * RCDT_NOISE table -- so the AVX-512 build can replace it with the
 * same 16-lane SIMD kernel.  The borrow chain is the exact 31-bit-limb (v
 * - Z - borrow)>>31 (see sampler_sigma2). */
void noise_magnitude_batch(int32_t *m_out, const uint8_t *rand)
{
    for (int s = 0; s < NOISE_BATCH; s++) {
        int group = s >> 3, lane = s & 7;
        int byte_base = group * 96 + lane * 4;
        uint32_t v0 = load_le32(rand + byte_base) & 0x7FFFFFFFU;
        uint32_t v1 = load_le32(rand + byte_base + 32) & 0x7FFFFFFFU;
        uint32_t v2 = load_le32(rand + byte_base + 64) & 0x7FFFFFFFU;
        int32_t m = 0;
#if defined(__clang__)
#    pragma clang loop vectorize(disable) interleave(disable)
#endif
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

/* Shared noise core (KAT-critical; identical math regardless of byte
 * source). Carves ONE noise poly from the persistent stream gs: per
 * mini-batch read a grouped CDT block (NOISE_BATCH magnitudes), then
 * per-candidate one tail byte (bit0 sign, bit1 zero-fold), branchless
 * accept (zero-fold only on m==0).  Leftover carries to the next poly. */
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

/* ============================================================
 * sigma = 1 noise sampler (keygen): magnitude via 10-entry RCDT, sign,
 * fold 0. Per coefficient: 12 bytes CDT (3x31) + 1 byte (bit0 sign, bit1
 * zero-fold).  (The per-attempt single-stream path below is retained for
 * tests; keygen uses the unified sample_noise_vec above.)
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
