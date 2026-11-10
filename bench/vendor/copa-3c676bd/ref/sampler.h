/*
 * sampler.h - discrete Gaussian samplers for IRS.
 *
 * Two samplers:
 *
 *  1) sample_gauss_N / sample_gauss_N_4x : BLISS-style sampler for D_{Z,
 * r} (r = 182/204/257).  k = 128 convolution: x ~ D_{Z+, r/128} via a
 * 93-bit RCDT (RCDT_Z), y ~ U{0..127}, candidate 128*x + y,
 * accepted with prob exp(-y(y+256x)/(2 r^2)) via approx_exp, then signed;
 * z=0 folded. The 4x variant defines the canonical stream layout for AVX2
 * KAT-compat.
 *
 *  2) sample_noise_N : D_{Z, 1} via the shared 10-entry RCDT
 * (RCDT_NOISE), no ApproxExp.  Used for the NTWE key noise f0, s0,
 * e.
 *
 * Both use one SHAKE-256 stream per call (seed + nonce).
 */
#ifndef LITHIUM_SAMPLER_H
#define LITHIUM_SAMPLER_H

#include <stddef.h>
#include <stdint.h>

#include "params.h"
#include "symmetric.h"

/* Gaussian mini-batch size.  32 (= two 16-lane AVX-512 sigma2 streams, the
 * IPC-optimal width) is the canonical default; the ref scalar loop and the
 * AVX2 kernel (which wraps its 16-sample core in a GAUSS_BATCH/16 outer
 * loop) adapt to it.  The mini-batch byte grouping depends on GAUSS_BATCH,
 * so changing it re-records the KAT (done deliberately,
 * ref==avx2==avx512). */
#ifndef GAUSS_BATCH
#    define GAUSS_BATCH 32
#endif

/* 93-bit reverse-CDT tables (3 x 31-bit limbs), defined in sampler.c (the
 * @@AUTOGEN:rcdt@@ region) and shared with optimized samplers via these
 * extern decls.  Generated + validated by tools/gen_rcdt.py. */
extern const uint32_t RCDT_Z[RCDT_ENTRIES]
                            [3]; /* BLISS base, sigma2=r/128 */
extern const uint32_t RCDT_NOISE[NOISE_RCDT_ENTRIES]
                                [3]; /* key noise, sigma=1 */

/* CDT random bytes per mini-batch: 12 bytes/sample, 3x31 AVX2 layout.
 * The canonical GAUSS_BATCH=32 is two 16-sample units, each 192 bytes. */
#define SIGMA2_RAND_BYTES (GAUSS_BATCH * 12)

/* y: GAUSS_BATCH * Y_BITS bits, bit-packed.  For GAUSS_BATCH=32 and
 * Y_BITS=7 this is 28 bytes. */
#define Y_RAND_BYTES ((GAUSS_BATCH * Y_BITS + 7) / 8)

/* Per-attempt tail: 8 bytes = 1 sign-for-zero bit + 63 rejection bits. */
#define GAUSS_RAND_BYTES 8

#define MINIBATCH_RAND_BYTES \
    (SIGMA2_RAND_BYTES + Y_RAND_BYTES + GAUSS_BATCH * GAUSS_RAND_BYTES)

#define SQUEEZE_NBLOCKS                                  \
    ((MINIBATCH_RAND_BYTES + STREAM256_BLOCKBYTES - 1) / \
     STREAM256_BLOCKBYTES)

#define GAUSS_BUF_SIZE \
    (MINIBATCH_RAND_BYTES + SQUEEZE_NBLOCKS * STREAM256_BLOCKBYTES)

static inline uint64_t load_le64(const uint8_t *p)
{
    return (uint64_t)p[0] | ((uint64_t)p[1] << 8) |
           ((uint64_t)p[2] << 16) | ((uint64_t)p[3] << 24) |
           ((uint64_t)p[4] << 32) | ((uint64_t)p[5] << 40) |
           ((uint64_t)p[6] << 48) | ((uint64_t)p[7] << 56);
}
static inline uint32_t load_le32(const uint8_t *p)
{
    return (uint32_t)p[0] | ((uint32_t)p[1] << 8) |
           ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}
static inline uint32_t ct_lt_u32(uint32_t a, uint32_t b)
{
    return (uint32_t)(((uint64_t)a - (uint64_t)b) >> 32) & 1U;
}
static inline uint64_t ct_lt_u64(uint64_t a, uint64_t b)
{
    /* Inputs used here are below 2^63, so the high bit of a-b is borrow.
     */
    return (a - b) >> 63;
}
static inline uint32_t ct_is_zero_u32(uint32_t x)
{
    return 1U ^ ((x | (0U - x)) >> 31);
}
static inline int16_t ct_sel_i16(uint32_t cond, int16_t x, int16_t y)
{
    uint32_t m = 0U - (cond & 1U);
    return (int16_t)(((uint32_t)(uint16_t)x & m) |
                     ((uint32_t)(uint16_t)y & ~m));
}
static inline int32_t ct_sel_i32(uint32_t cond, int32_t x, int32_t y)
{
    uint32_t m = 0U - (cond & 1U);
    return (int32_t)(((uint32_t)x & m) | ((uint32_t)y & ~m));
}

/* BLISS base CDT: GAUSS_BATCH magnitudes in [0, RCDT_ENTRIES] into
 * z_out. Scalar in ref; replaced by an 8-way SIMD kernel in the AVX2
 * build. */
int sampler_sigma2(int32_t *z_out, const uint8_t *rand);

/* Mini-batch helpers exposed for the AVX2 4-way sampler (KAT-critical: the
 * 4-way path reuses the exact y bit-packing, Q60 input, and
 * accept/sign/fold logic). */
void sampler_y(uint8_t y_out[GAUSS_BATCH], const uint8_t *rand);
uint64_t gauss_a_q60(uint32_t x, uint32_t y);
/* Batched a_q60 over the mini-batch (out[i] = gauss_a_q60(z[i], y[i])).
 * Scalar in ref/avx2; an 8-lane vpmadd52 (IFMA) kernel in the AVX-512
 * build (guarded out by USE_AVX512_SAMPLER) -- exact, KAT-neutral.
 */
void gauss_a_q60_batch(uint64_t *out, const int32_t *z, const uint8_t *y,
                       int n);
int noise_magnitude(
    const uint8_t *rand12); /* 10-entry RCDT scan, [0,10] */
int gauss_finalize(int32_t *r, uint32_t x, uint32_t y, uint64_t exp_q63,
                   const uint8_t *rej_rand, const uint8_t *signs,
                   size_t idx);

/* D_{Z,r}: len signed samples (|r[i]| <= 128*ENTRIES + 127). */
void sample_gauss_N(int32_t *r, const uint8_t seed[SEEDBYTES],
                    uint64_t nonce, size_t len);

/* 4-way (canonical layout for AVX2 KAT-compat). Unused lanes: len=0. */
void sample_gauss_N_4x(int32_t *r0, int32_t *r1, int32_t *r2, int32_t *r3,
                       const uint8_t seed[SEEDBYTES], uint64_t nonce0,
                       uint64_t nonce1, uint64_t nonce2, uint64_t nonce3,
                       size_t len0, size_t len1, size_t len2, size_t len3);

/* ---- unified "combined N-lane" Gaussian y-vector sampler ----
 * The whole y-vector (cnt polys) is sampled from ONE continuous
 * round-robin stream of GAUSS_COMB_LANES domain-separated sub-streams
 * SHAKE256(seed || nonce_LE8 || j), j=0..GAUSS_COMB_LANES-1 (41-byte
 * inputs, so no collision with the 40-byte single-nonce streams).  Polys
 * are consumed back-to-back from one persistent buffer (leftover carries
 * across polys), so the Keccak engine is fully utilised and the ONLY waste
 * is the final partial unit.
 *
 * GAUSS_COMB_LANES = 8 matches the 8-way AVX-512 Keccak (keccak8x): the
 * dominant Gaussian SHAKE cost runs at full x8 width.  The backends emit
 * the IDENTICAL unit ([lane0 blk][lane1]...[lane7]) -- ref = 8 scalar
 * SHAKE, avx2 = 2x keccak4x, avx512 = 1x keccak8x -- so all three are
 * bit-identical; the shared sampling loop (gauss_stream_poly) is
 * lane-count-agnostic.  (This 4->8 widening, with GAUSS_BATCH=32, is the
 * deliberate KAT re-record.) */
#define GAUSS_COMB_LANES 8
#define GAUSS_COMB_BUF \
    (MINIBATCH_RAND_BYTES + GAUSS_COMB_LANES * STREAM256_BLOCKBYTES + 64)
typedef struct {
    void *st; /* backend keccak state(s)              */
    void (*squeeze)(void *st,
                    uint8_t *out); /* append GAUSS_COMB_LANES blocks */
    size_t pos, avail;             /* cursor into buf                */
    uint8_t buf[GAUSS_COMB_BUF];
} gauss_stream;

void gauss_stream_poly(gauss_stream *gs,
                       int32_t *r); /* shared gauss core */

/* Backend (ref: 8 scalar SHAKE; avx2: 2x keccak4x; avx512: 1x keccak8x) --
 * selected like the _4x/_8x fns.  Samples cnt polys of N coeffs
 * each into r[0..cnt-1] from one unified stream. */
void sample_gauss_vec(int32_t *const r[], int cnt,
                      const uint8_t seed[SKSEEDBYTES], uint64_t nonce);

/* ---- unified "combined N-lane" key-noise (D_{Z,1}) vector sampler ----
 * Mirrors the Gaussian unified stream: the whole noise vector (cnt =
 * K polys: f0, s0(ell), e(m)) is carved back-to-back from ONE
 * round-robin of GAUSS_COMB_LANES sub-streams -- no per-poly partitioning,
 * so the 8-way Keccak never idles a lane (the only waste is the final
 * partial unit).
 *
 * The noise mini-batch reuses the sigma2 GROUPED byte layout so its
 * magnitude scan can be the same 16-lane SIMD kernel: NOISE_BATCH
 * candidates laid out group*96 + lane*4 (v0/v1/v2 at +0/+32/+64), 12 B
 * each, followed by a block of NOISE_BATCH per-candidate tail bytes (bit0
 * sign, bit1 zero-fold).  noise_magnitude_batch computes all NOISE_BATCH
 * magnitudes; the per-candidate sign/zero-fold/accept then runs serially
 * (branchless). This regroups the noise bytes vs the per-attempt
 * sample_noise_N, so it re-records the KAT (bundled with the
 * GAUSS_BATCH/8-substream change). */
#ifndef NOISE_BATCH
#    define NOISE_BATCH 32
#endif
#define NOISE_CDT_BYTES (NOISE_BATCH * 12) /* grouped 3x31-bit limbs */
/* Per-candidate tail is just 2 bits (bit0 sign, bit1 zero-fold), so 4
 * candidates share one tail byte: NOISE_BATCH/4 tail bytes per mini-batch
 * (NOISE_BATCH is a multiple of 4 in every mode). */
#define NOISE_MINIBATCH_RAND_BYTES (NOISE_CDT_BYTES + NOISE_BATCH / 4)

/* NOISE_BATCH magnitudes (each in [0, NOISE_RCDT_ENTRIES]) from a
 * grouped CDT block.  Scalar in ref/avx2; replaced by a 16-lane SIMD
 * kernel in the AVX-512 build (guarded out by USE_AVX512_SAMPLER).
 */
void noise_magnitude_batch(int32_t *m_out, const uint8_t *cdt);
void noise_stream_poly(gauss_stream *gs,
                       int32_t *r); /* shared noise core */
void sample_noise_vec(int32_t *const r[], int cnt,
                      const uint8_t seed[SKSEEDBYTES], uint64_t nonce);

/* D_{Z,1} key noise: len signed samples in roughly [-9, 9]. */
void sample_noise_N(int32_t *r, const uint8_t seed[SEEDBYTES],
                    uint64_t nonce, size_t len);

/* 4-way noise (canonical layout for AVX2 KAT-compat). Unused lanes: len=0.
 */
void sample_noise_N_4x(int32_t *r0, int32_t *r1, int32_t *r2, int32_t *r3,
                       const uint8_t seed[SEEDBYTES], uint64_t nonce0,
                       uint64_t nonce1, uint64_t nonce2, uint64_t nonce3,
                       size_t len0, size_t len1, size_t len2, size_t len3);

#endif /* LITHIUM_SAMPLER_H */
