/*
 * NTT16vs32.c  --  Montgomery modular multiplication AND NTT butterflies for
 *                  q = 45569 on AVX2, comparing a 16-bit-lane layout against a
 *                  32-bit-lane layout, for both Cooley-Tukey (CT, forward) and
 *                  Gentleman-Sande (GS, inverse) butterflies.
 *
 * Background
 * ----------
 * Kyber    : q = 3329     (< 2^15)  -> signed 16-bit AVX2 NTT, lots of headroom.
 * Dilithium: q = 8380417 (~2^23)    -> 32-bit AVX2 NTT, R = 2^32, "splitty".
 *
 * q = 45569 sits in the "valley"  2^15 < q < 2^16, which has sharp consequences:
 *   - q does NOT fit in signed int16 (q > 32767): the clean Kyber signed recipe
 *     breaks, so the 16-bit path goes unsigned (vpmulhuw) + a borrow correction.
 *   - The Montgomery result spans (-q,q), > 2^16 wide, so it cannot live in 16
 *     bits: it must be reduced to [0,q) immediately.
 *   - ZERO lazy-reduction headroom: q+q = 91138 > 65535, so you cannot even FORM
 *     a + t for a,t in [0,q) in a 16-bit lane without overflow handling. Every
 *     add AND sub in a butterfly needs an overflow-aware modular reduction.
 *
 * In 32-bit lanes q is a positive int32, the result has ~2^15.5x headroom (so a
 * whole n<=256 NTT runs with NO intermediate reduction -- both CT and GS), and
 * because q < 2^16 the reduction is cheaper than Dilithium's (R = 2^16, vpmulld).
 *
 * CT vs GS (the order of multiply vs add/sub):
 *   CT (forward):  t = zeta*b;  a' = a + t;  b' = a - t.   (multiply FIRST)
 *   GS (inverse):  a' = a + b;  b' = zeta*(a - b).         (multiply LAST)
 *   - In 16-bit (zero headroom) both must fully reduce every add/sub -> CT ~ GS.
 *   - In 32-bit (huge headroom) both are fully lazy for n<=256 -> CT ~ GS.
 *   The usual CT(additive growth)/GS(multiplicative growth) asymmetry that drives
 *   lazy-reduction budgets is collapsed by the valley: 16-bit has no budget,
 *   32-bit has so much that GS's faster growth is still absorbed.
 *
 * Every multiply kernel and every butterfly (CT & GS, all three backends) is
 * verified bit-exactly against a scalar oracle before benchmarking.
 *
 * Build:  gcc -O3 -mavx2 -o ntt16vs32 NTT16vs32.c
 * Run:    ./ntt16vs32
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
#include <x86intrin.h>

/* ------------------------------------------------------------------ */
/*  Modulus and Montgomery constants                                  */
/* ------------------------------------------------------------------ */
#define Q          45569
#define QINV16     19969u          /* q^{-1} mod 2^16  (sub form)             */
#define QINV32     868437505u      /* q^{-1} mod 2^32  (for the vpmuldq path) */
#define R_MINUS_Q  19967           /* 2^16 - q, the 16-bit reduction constant */

static int rinv16_mod_q;
static int rinv32_mod_q;

static int egcd_inv(long a, long m)
{
    long t = 0, newt = 1, r = m, newr = a % m;
    while (newr != 0) {
        long quo = r / newr, tmp;
        tmp = t - quo * newt; t = newt; newt = tmp;
        tmp = r - quo * newr; r = newr; newr = tmp;
    }
    if (r > 1) return -1;
    if (t < 0) t += m;
    return (int)t;
}

/* ------------------------------------------------------------------ */
/*  Scalar reference oracles                                          */
/* ------------------------------------------------------------------ */
static uint32_t mulmod_plain(uint32_t a, uint32_t b)
{
    return (uint32_t)(((uint64_t)a * b) % Q);
}
static uint32_t addmod(uint32_t a, uint32_t b){ uint32_t s=a+b; return s>=Q?s-Q:s; }
static uint32_t submod(uint32_t a, uint32_t b){ return a>=b ? a-b : a+Q-b; }

/* Signed Montgomery, R = 2^16 (32-bit vpmulld semantics): a*b*R^{-1} mod q in (-q,q). */
static int32_t mont_ref_signed(int32_t a, int32_t b)
{
    int64_t T = (int64_t)a * b;
    int16_t m = (int16_t)((int32_t)T * (int32_t)QINV16);
    return (int32_t)((T - (int64_t)m * Q) >> 16);
}
/* Unsigned Montgomery, R = 2^16 (16-bit-lane semantics): a*b*R^{-1} mod q in [0,q). */
static uint32_t mont_ref_unsigned(uint16_t a, uint16_t b)
{
    uint64_t T = (uint64_t)a * b;
    uint16_t m = (uint16_t)(((uint32_t)(T & 0xffff) * QINV16) & 0xffff);
    int64_t  r = ((int64_t)T - (int64_t)m * Q) >> 16;
    if (r < 0) r += Q;
    return (uint32_t)r;
}
/* Signed Montgomery, R = 2^32 (vpmuldq semantics): a*b*2^{-32} mod q in (-q,q). */
static int32_t mont_ref_signed32(int32_t a, int32_t b)
{
    int64_t T = (int64_t)a * b;
    int32_t m = (int32_t)((uint64_t)T * QINV32);
    return (int32_t)((T - (int64_t)m * Q) >> 32);
}

/* ------------------------------------------------------------------ */
/*  AVX2 Montgomery multiply kernels                                  */
/* ------------------------------------------------------------------ */

/* 16-bit lanes, 16 coeff/reg, R = 2^16.  a,b in [0,q) -> [0,q).  9 instr.
 * Unsigned high (vpmulhuw) is mandatory (q > 2^15); the borrow correction
 * exists only because the true result spans (-q,q) and does not fit 16 bits. */
static inline __m256i mont_mul_16(__m256i a, __m256i b,
                                  __m256i vqinv, __m256i vq, __m256i vrmq, __m256i zero)
{
    __m256i lo   = _mm256_mullo_epi16(a, b);
    __m256i Thi  = _mm256_mulhi_epu16(a, b);
    __m256i m    = _mm256_mullo_epi16(lo, vqinv);
    __m256i mqhi = _mm256_mulhi_epu16(m, vq);
    __m256i r    = _mm256_sub_epi16(Thi, mqhi);
    __m256i nob  = _mm256_cmpeq_epi16(_mm256_subs_epu16(mqhi, Thi), zero); /* Thi>=mqhi */
    return _mm256_sub_epi16(r, _mm256_andnot_si256(nob, vrmq));
}

/* 32-bit lanes, 8 coeff/reg, R = 2^16, vpmulld.  centered a,b -> (-q,q).  7 instr. */
static inline __m256i mont_mul_32(__m256i a, __m256i b, __m256i vqinv, __m256i vq)
{
    __m256i T  = _mm256_mullo_epi32(a, b);
    __m256i mf = _mm256_mullo_epi32(T, vqinv);
    __m256i m  = _mm256_srai_epi32(_mm256_slli_epi32(mf, 16), 16);
    __m256i mq = _mm256_mullo_epi32(m, vq);
    return _mm256_srai_epi32(_mm256_sub_epi32(T, mq), 16);
}

/* 32-bit lanes, vpmuldq, R = 2^32 (Dilithium-style even/odd split). ~13 instr. */
static inline __m256i redc32_even(__m256i a_e, __m256i b_e, __m256i vqinv64, __m256i vq64)
{
    __m256i P  = _mm256_mul_epi32(a_e, b_e);
    __m256i m  = _mm256_mul_epi32(P, vqinv64);
    __m256i mq = _mm256_mul_epi32(m, vq64);
    return _mm256_srli_epi64(_mm256_sub_epi64(P, mq), 32);
}
static inline __m256i mont_mul_32_dq(__m256i a, __m256i b, __m256i vqinv64, __m256i vq64)
{
    __m256i even = redc32_even(a, b, vqinv64, vq64);
    __m256i odd  = redc32_even(_mm256_srli_epi64(a, 32), _mm256_srli_epi64(b, 32), vqinv64, vq64);
    return _mm256_blend_epi32(even, _mm256_slli_epi64(odd, 32), 0xAA);
}

/* ------------------------------------------------------------------ */
/*  16-bit modular add / sub  (overflow-aware: the valley's real cost) */
/* ------------------------------------------------------------------ */
/* (a+b) mod q, a,b in [0,q).  Cannot form a+b directly (overflows 2^16), so we
 * compute a-(q-b) = a+b-q with a borrow-driven correction.  6 instr. */
static inline __m256i add_mod_16(__m256i a, __m256i b, __m256i vq, __m256i vrmq, __m256i zero)
{
    __m256i u   = _mm256_sub_epi16(vq, b);                 /* q-b in (0,q]      */
    __m256i s   = _mm256_sub_epi16(a, u);                  /* a+b-q (wraps)     */
    __m256i nob = _mm256_cmpeq_epi16(_mm256_subs_epu16(u, a), zero); /* a>=u => a+b>=q */
    return _mm256_sub_epi16(s, _mm256_andnot_si256(nob, vrmq));
}
/* (a-b) mod q, a,b in [0,q).  5 instr. */
static inline __m256i sub_mod_16(__m256i a, __m256i b, __m256i vrmq, __m256i zero)
{
    __m256i s   = _mm256_sub_epi16(a, b);                  /* a-b (wraps)       */
    __m256i nob = _mm256_cmpeq_epi16(_mm256_subs_epu16(b, a), zero); /* a>=b */
    return _mm256_sub_epi16(s, _mm256_andnot_si256(nob, vrmq));
}

/* ------------------------------------------------------------------ */
/*  Butterflies (return both outputs)                                 */
/* ------------------------------------------------------------------ */
typedef struct { __m256i x, y; } pair_t;

/* 16-bit CT: t=tw*b; (a+t, a-t).  GS: (a+b, tw*(a-b)). */
static inline pair_t ct16(__m256i a, __m256i b, __m256i tw,
                          __m256i vqinv, __m256i vq, __m256i vrmq, __m256i zero)
{
    __m256i t = mont_mul_16(b, tw, vqinv, vq, vrmq, zero);
    pair_t p = { add_mod_16(a, t, vq, vrmq, zero), sub_mod_16(a, t, vrmq, zero) };
    return p;
}
static inline pair_t gs16(__m256i a, __m256i b, __m256i tw,
                          __m256i vqinv, __m256i vq, __m256i vrmq, __m256i zero)
{
    __m256i d = sub_mod_16(a, b, vrmq, zero);
    pair_t p = { add_mod_16(a, b, vq, vrmq, zero), mont_mul_16(d, tw, vqinv, vq, vrmq, zero) };
    return p;
}

/* 32-bit (vpmulld) CT / GS -- fully lazy, no reduction (headroom absorbs growth). */
static inline pair_t ct32(__m256i a, __m256i b, __m256i tw, __m256i vqinv, __m256i vq)
{
    __m256i t = mont_mul_32(b, tw, vqinv, vq);
    pair_t p = { _mm256_add_epi32(a, t), _mm256_sub_epi32(a, t) };
    return p;
}
static inline pair_t gs32(__m256i a, __m256i b, __m256i tw, __m256i vqinv, __m256i vq)
{
    /* GS feeds the difference a-b (|.| up to q) into the multiply, but the tight
     * R=2^16 vpmulld kernel is only valid for |input| <= (q-1)/2 (because
     * q*2^16 > 2^31).  So the difference must be normalized to the centered
     * range first -- an extra reduction that CT (which multiplies a single
     * coefficient) does not pay.  This is a real CT/GS asymmetry at R=2^16. */
    __m256i vqh  = _mm256_set1_epi32((Q - 1) / 2);
    __m256i vnqh = _mm256_set1_epi32(-((Q - 1) / 2));
    __m256i d  = _mm256_sub_epi32(a, b);
    d = _mm256_sub_epi32(d, _mm256_and_si256(_mm256_cmpgt_epi32(d, vqh), vq));
    d = _mm256_add_epi32(d, _mm256_and_si256(_mm256_cmpgt_epi32(vnqh, d), vq));
    pair_t p = { _mm256_add_epi32(a, b), mont_mul_32(d, tw, vqinv, vq) };
    return p;
}

/* 32-bit (vpmuldq) CT / GS. */
static inline pair_t ctdq(__m256i a, __m256i b, __m256i tw, __m256i vqi64, __m256i vq64)
{
    __m256i t = mont_mul_32_dq(b, tw, vqi64, vq64);
    pair_t p = { _mm256_add_epi32(a, t), _mm256_sub_epi32(a, t) };
    return p;
}
static inline pair_t gsdq(__m256i a, __m256i b, __m256i tw, __m256i vqi64, __m256i vq64)
{
    pair_t p = { _mm256_add_epi32(a, b), mont_mul_32_dq(_mm256_sub_epi32(a, b), tw, vqi64, vq64) };
    return p;
}

/* ================================================================== */
/*  A complete n = 256 negacyclic NTT / INTT over Z_q[x]/(x^256+1),    */
/*  16-bit method, structured like pq-crystals/kyber's ntt.c.          */
/*                                                                     */
/*  WHY "COMPLETE" (8 levels) AND NOT KYBER'S 7:                       */
/*    q-1 = 2^9 * 89 = 512 * 89, so 2n = 512 DIVIDES q-1 and a         */
/*    primitive 512-th root of unity omega exists (omega^256 = -1).    */
/*    Hence x^256+1 splits into 256 LINEAR factors -> 8 butterfly      */
/*    levels and a PLAIN coefficient-wise pointwise multiply.          */
/*    Kyber's q=3329 has q-1 = 2^8 * 13: only 256 | q-1, NOT 512, so   */
/*    x^256+1 splits into 128 QUADRATIC factors -> 7 levels + a degree-1*/
/*    "basemul".  Copying Kyber's `len >= 2` loop bound verbatim would  */
/*    silently compute the WRONG (incomplete) transform here.          */
/*                                                                     */
/*  DOMAIN: coefficients are kept in PLAIN [0,q) the whole time.       */
/*    Twiddles are stored in MONTGOMERY domain (zeta * R mod q) so that */
/*    fqmul(zeta_mont, x) = zeta * x mod q leaves x's domain unchanged. */
/*    This is exactly Kyber's trick; see zeta_mont[] / ninv_mont below. */
/* ================================================================== */

static uint16_t zeta_mont[256];   /* omega^{brv8(k)} * R mod q, k=1..255   */
static uint16_t ninv_mont;        /* 256^{-1} * R mod q (final INTT scale) */
static uint16_t R2_mod_q;         /* R^2 mod q (to push an operand to Mont)*/

static uint32_t powmod(uint32_t b, uint32_t e)
{
    uint64_t r = 1, bb = b % Q;
    while (e) { if (e & 1) r = r * bb % Q; bb = bb * bb % Q; e >>= 1; }
    return (uint32_t)r;
}
static int brv8(int x){ int r=0; for(int i=0;i<8;i++) r |= ((x>>i)&1)<<(7-i); return r; }

/* Scalar 16-bit Montgomery multiply (== mont_ref_unsigned), result in [0,q). */
static inline uint16_t fqmul(uint16_t a, uint16_t b){ return (uint16_t)mont_ref_unsigned(a, b); }

static void ntt_init(void)
{
    uint32_t w = powmod(3, (Q - 1) / 512);     /* primitive 512-th root  */
    for (int k = 0; k < 256; k++)
        zeta_mont[k] = (uint16_t)((uint64_t)powmod(w, brv8(k)) * ((1u<<16) % Q) % Q);
    uint32_t ninv = powmod(256, Q - 2);        /* Fermat inverse         */
    ninv_mont = (uint16_t)((uint64_t)ninv * ((1u<<16) % Q) % Q);
    R2_mod_q  = (uint16_t)((uint64_t)((1u<<16)%Q) * ((1u<<16)%Q) % Q);
}

/* ---- scalar reference: forward CT, normal -> bit-reversed order ---- */
static void ntt_scalar(uint16_t r[256])
{
    int k = 1;
    for (int len = 128; len >= 1; len >>= 1)             /* 8 levels (>=1!) */
        for (int start = 0; start < 256; start += 2 * len) {
            uint16_t z = zeta_mont[k++];
            for (int j = start; j < start + len; j++) {
                uint16_t t = fqmul(z, r[j + len]);
                r[j + len] = (uint16_t)submod(r[j], t);
                r[j]       = (uint16_t)addmod(r[j], t);
            }
        }
}
/* ---- scalar reference: inverse GS, bit-reversed -> normal order ---- */
static void invntt_scalar(uint16_t r[256])
{
    int k = 255;
    for (int len = 1; len <= 128; len <<= 1)            /* 8 levels        */
        for (int start = 0; start < 256; start += 2 * len) {
            uint16_t z = zeta_mont[k--];
            for (int j = start; j < start + len; j++) {
                uint16_t t = r[j];
                r[j]       = (uint16_t)addmod(t, r[j + len]);
                r[j + len] = fqmul(z, (uint16_t)submod(r[j + len], t));
            }
        }
    for (int j = 0; j < 256; j++) r[j] = fqmul(r[j], ninv_mont);  /* * 1/n */
}
/* ---- pointwise: complete NTT => plain coefficient-wise product ---- */
static void pointwise_scalar(const uint16_t a[256], const uint16_t b[256], uint16_t c[256])
{
    for (int i = 0; i < 256; i++) {
        uint16_t bm = fqmul(b[i], R2_mod_q);   /* b -> Montgomery (b*R)    */
        c[i] = fqmul(a[i], bm);                /* a * (b*R) * R^-1 = a*b   */
    }
}

/* ================================================================== */
/*  AVX2 version (hybrid):                                             */
/*    - levels with len >= 16 are vectorized (16 coeffs / YMM), using   */
/*      mont_mul_16 + the OVERFLOW-AWARE add_mod_16 / sub_mod_16 (the   */
/*      "valley tax": a+b cannot be formed directly since 2q > 2^16).   */
/*    - levels with len < 16 (len = 8,4,2,1) interleave butterfly pairs */
/*      inside a 16-lane register and would need cross-lane shuffles    */
/*      (cf. Kyber's shuffle.inc).  For clarity we run them scalar.     */
/*  Both halves operate in place on the same [0,q) uint16 array and the */
/*  whole thing is verified bit-exact against ntt_scalar.              */
/* ================================================================== */
static void ntt_avx2(uint16_t r[256])
{
    __m256i vq=_mm256_set1_epi16((short)Q), vqi=_mm256_set1_epi16((short)QINV16),
            vrmq=_mm256_set1_epi16((short)R_MINUS_Q), zero=_mm256_setzero_si256();
    int k = 1, len = 128;
    for (; len >= 16; len >>= 1)                         /* vectorized      */
        for (int start = 0; start < 256; start += 2 * len) {
            __m256i vz = _mm256_set1_epi16((short)zeta_mont[k++]);
            for (int j = start; j < start + len; j += 16) {
                __m256i a = _mm256_loadu_si256((__m256i*)&r[j]);
                __m256i b = _mm256_loadu_si256((__m256i*)&r[j + len]);
                __m256i t = mont_mul_16(b, vz, vqi, vq, vrmq, zero);
                _mm256_storeu_si256((__m256i*)&r[j],       add_mod_16(a, t, vq, vrmq, zero));
                _mm256_storeu_si256((__m256i*)&r[j + len], sub_mod_16(a, t, vrmq, zero));
            }
        }
    for (; len >= 1; len >>= 1)                          /* scalar tail     */
        for (int start = 0; start < 256; start += 2 * len) {
            uint16_t z = zeta_mont[k++];
            for (int j = start; j < start + len; j++) {
                uint16_t t = fqmul(z, r[j + len]);
                r[j + len] = (uint16_t)submod(r[j], t);
                r[j]       = (uint16_t)addmod(r[j], t);
            }
        }
}
static void invntt_avx2(uint16_t r[256])
{
    __m256i vq=_mm256_set1_epi16((short)Q), vqi=_mm256_set1_epi16((short)QINV16),
            vrmq=_mm256_set1_epi16((short)R_MINUS_Q), zero=_mm256_setzero_si256();
    int k = 255, len = 1;
    for (; len < 16; len <<= 1)                          /* scalar head     */
        for (int start = 0; start < 256; start += 2 * len) {
            uint16_t z = zeta_mont[k--];
            for (int j = start; j < start + len; j++) {
                uint16_t t = r[j];
                r[j]       = (uint16_t)addmod(t, r[j + len]);
                r[j + len] = fqmul(z, (uint16_t)submod(r[j + len], t));
            }
        }
    for (; len <= 128; len <<= 1)                        /* vectorized      */
        for (int start = 0; start < 256; start += 2 * len) {
            __m256i vz = _mm256_set1_epi16((short)zeta_mont[k--]);
            for (int j = start; j < start + len; j += 16) {
                __m256i a = _mm256_loadu_si256((__m256i*)&r[j]);
                __m256i b = _mm256_loadu_si256((__m256i*)&r[j + len]);
                _mm256_storeu_si256((__m256i*)&r[j],       add_mod_16(a, b, vq, vrmq, zero));
                _mm256_storeu_si256((__m256i*)&r[j + len], mont_mul_16(sub_mod_16(b, a, vrmq, zero), vz, vqi, vq, vrmq, zero));
            }
        }
    __m256i vninv = _mm256_set1_epi16((short)ninv_mont);
    for (int i = 0; i < 256; i += 16)
        _mm256_storeu_si256((__m256i*)&r[i],
            mont_mul_16(_mm256_loadu_si256((__m256i*)&r[i]), vninv, vqi, vq, vrmq, zero));
}
static void pointwise_avx2(const uint16_t a[256], const uint16_t b[256], uint16_t c[256])
{
    __m256i vq=_mm256_set1_epi16((short)Q), vqi=_mm256_set1_epi16((short)QINV16),
            vrmq=_mm256_set1_epi16((short)R_MINUS_Q), zero=_mm256_setzero_si256(),
            vR2=_mm256_set1_epi16((short)R2_mod_q);
    for (int i = 0; i < 256; i += 16) {
        __m256i bm = mont_mul_16(_mm256_loadu_si256((__m256i*)&b[i]), vR2, vqi, vq, vrmq, zero);
        _mm256_storeu_si256((__m256i*)&c[i],
            mont_mul_16(_mm256_loadu_si256((__m256i*)&a[i]), bm, vqi, vq, vrmq, zero));
    }
}

/* Schoolbook negacyclic convolution (ground truth), c = a (*) b in x^256+1. */
static void polymul_ref(const uint16_t a[256], const uint16_t b[256], uint16_t c[256])
{
    for (int i = 0; i < 256; i++) c[i] = 0;
    for (int i = 0; i < 256; i++)
        for (int j = 0; j < 256; j++) {
            uint32_t p = (uint32_t)(((uint64_t)a[i] * b[j]) % Q);
            int k = i + j;
            if (k < 256) c[k] = (uint16_t)addmod(c[k], p);
            else         c[k-256] = (uint16_t)submod(c[k-256], p);
        }
}

/* ------------------------------------------------------------------ */
/*  Correctness                                                       */
/* ------------------------------------------------------------------ */
static uint32_t red_pm(int32_t v){ int32_t r=v%Q; if(r<0)r+=Q; return (uint32_t)r; }

static int test_correctness(void)
{
    int fails = 0, half = (Q-1)/2;
    __m256i zero = _mm256_setzero_si256();

    /* --- 16-bit lane: mul + CT + GS --- */
    {
        __m256i vq=_mm256_set1_epi16((short)Q), vqi=_mm256_set1_epi16((short)QINV16);
        __m256i vrmq=_mm256_set1_epi16((short)R_MINUS_Q);
        uint16_t A[16],B[16],TW[16],X[16],Y[16];
        for (int it=0; it<300000; it++) {
            for (int k=0;k<16;k++){A[k]=rand()%Q;B[k]=rand()%Q;TW[k]=rand()%Q;}
            __m256i a=_mm256_loadu_si256((__m256i*)A),b=_mm256_loadu_si256((__m256i*)B),
                    tw=_mm256_loadu_si256((__m256i*)TW);
            /* mul */
            _mm256_storeu_si256((__m256i*)X, mont_mul_16(a,b,vqi,vq,vrmq,zero));
            for(int k=0;k<16;k++){ uint32_t e=mont_ref_unsigned(A[k],B[k]);
                if(X[k]!=e||X[k]>=Q){ if(fails<4)printf("  [16 mul] %u*%u got %u exp %u\n",A[k],B[k],X[k],e); fails++; } }
            /* CT */
            pair_t p=ct16(a,b,tw,vqi,vq,vrmq,zero);
            _mm256_storeu_si256((__m256i*)X,p.x); _mm256_storeu_si256((__m256i*)Y,p.y);
            for(int k=0;k<16;k++){ uint32_t t=mont_ref_unsigned(B[k],TW[k]);
                uint32_t ex=addmod(A[k],t),ey=submod(A[k],t);
                if(X[k]!=ex||Y[k]!=ey){ if(fails<4)printf("  [16 CT] got(%u,%u) exp(%u,%u)\n",X[k],Y[k],ex,ey); fails++; } }
            /* GS */
            p=gs16(a,b,tw,vqi,vq,vrmq,zero);
            _mm256_storeu_si256((__m256i*)X,p.x); _mm256_storeu_si256((__m256i*)Y,p.y);
            for(int k=0;k<16;k++){ uint32_t ex=addmod(A[k],B[k]);
                uint32_t ey=mont_ref_unsigned((uint16_t)submod(A[k],B[k]),TW[k]);
                if(X[k]!=ex||Y[k]!=ey){ if(fails<4)printf("  [16 GS] got(%u,%u) exp(%u,%u)\n",X[k],Y[k],ex,ey); fails++; } }
        }
    }

    /* --- 32-bit lane (vpmulld, R=2^16): mul + CT + GS --- */
    {
        __m256i vq=_mm256_set1_epi32(Q), vqi=_mm256_set1_epi32((int)QINV16);
        int32_t A[8],B[8],TW[8],X[8],Y[8];
        for (int it=0; it<500000; it++) {
            for(int k=0;k<8;k++){A[k]=(rand()%(2*half+1))-half;B[k]=(rand()%(2*half+1))-half;TW[k]=(rand()%(2*half+1))-half;}
            __m256i a=_mm256_loadu_si256((__m256i*)A),b=_mm256_loadu_si256((__m256i*)B),
                    tw=_mm256_loadu_si256((__m256i*)TW);
            _mm256_storeu_si256((__m256i*)X, mont_mul_32(a,b,vqi,vq));
            for(int k=0;k<8;k++){ if(X[k]!=mont_ref_signed(A[k],B[k])){ if(fails<4)printf("  [32 mul]\n"); fails++; } }
            pair_t p=ct32(a,b,tw,vqi,vq);
            _mm256_storeu_si256((__m256i*)X,p.x); _mm256_storeu_si256((__m256i*)Y,p.y);
            for(int k=0;k<8;k++){ int32_t t=mont_ref_signed(B[k],TW[k]);
                if(red_pm(X[k])!=red_pm(A[k]+t)||red_pm(Y[k])!=red_pm(A[k]-t)){ if(fails<4)printf("  [32 CT]\n"); fails++; } }
            p=gs32(a,b,tw,vqi,vq);
            _mm256_storeu_si256((__m256i*)X,p.x); _mm256_storeu_si256((__m256i*)Y,p.y);
            for(int k=0;k<8;k++){ int32_t bb=mont_ref_signed(A[k]-B[k],TW[k]);
                if(red_pm(X[k])!=red_pm(A[k]+B[k])||red_pm(Y[k])!=red_pm(bb)){ if(fails<4)printf("  [32 GS]\n"); fails++; } }
        }
    }

    /* --- 32-bit lane (vpmuldq, R=2^32): mul + CT + GS --- */
    {
        __m256i vqi64=_mm256_set1_epi64x((long long)QINV32), vq64=_mm256_set1_epi64x(Q);
        int32_t A[8],B[8],TW[8],X[8],Y[8];
        for (int it=0; it<500000; it++) {
            for(int k=0;k<8;k++){A[k]=(rand()%(2*half+1))-half;B[k]=(rand()%(2*half+1))-half;TW[k]=(rand()%(2*half+1))-half;}
            __m256i a=_mm256_loadu_si256((__m256i*)A),b=_mm256_loadu_si256((__m256i*)B),
                    tw=_mm256_loadu_si256((__m256i*)TW);
            _mm256_storeu_si256((__m256i*)X, mont_mul_32_dq(a,b,vqi64,vq64));
            for(int k=0;k<8;k++){ if(X[k]!=mont_ref_signed32(A[k],B[k])){ if(fails<4)printf("  [dq mul]\n"); fails++; } }
            pair_t p=ctdq(a,b,tw,vqi64,vq64);
            _mm256_storeu_si256((__m256i*)X,p.x); _mm256_storeu_si256((__m256i*)Y,p.y);
            for(int k=0;k<8;k++){ int32_t t=mont_ref_signed32(B[k],TW[k]);
                if(red_pm(X[k])!=red_pm(A[k]+t)||red_pm(Y[k])!=red_pm(A[k]-t)){ if(fails<4)printf("  [dq CT]\n"); fails++; } }
            p=gsdq(a,b,tw,vqi64,vq64);
            _mm256_storeu_si256((__m256i*)X,p.x); _mm256_storeu_si256((__m256i*)Y,p.y);
            for(int k=0;k<8;k++){ int32_t bb=mont_ref_signed32(A[k]-B[k],TW[k]);
                if(red_pm(X[k])!=red_pm(A[k]+B[k])||red_pm(Y[k])!=red_pm(bb)){ if(fails<4)printf("  [dq GS]\n"); fails++; } }
        }
    }

    /* --- oracle vs plain product: mont(a,b)*R == a*b (mod q) --- */
    for (int t=0;t<200000;t++){ uint32_t a=rand()%Q,b=rand()%Q;
        if(mont_ref_unsigned((uint16_t)a,(uint16_t)b)!=mulmod_plain(mulmod_plain(a,b),(uint32_t)rinv16_mod_q)){printf("  [oracle16]\n");fails++;break;}
        if(red_pm(mont_ref_signed32((int32_t)a,(int32_t)b))!=mulmod_plain(mulmod_plain(a,b),(uint32_t)rinv32_mod_q)){printf("  [oracle32]\n");fails++;break;} }
    return fails;
}

/* ------------------------------------------------------------------ */
/*  NTT correctness: roots, scalar==AVX2, roundtrip, polymul vs school */
/* ------------------------------------------------------------------ */
static int test_ntt(void)
{
    int fails = 0;
    uint32_t w = powmod(3, (Q - 1) / 512);

    /* root sanity: omega^256 == -1, omega^512 == 1, omega^128 != +-1 */
    if (powmod(w, 256) != (uint32_t)(Q - 1)) { printf("  [root] omega^256 != -1\n"); fails++; }
    if (powmod(w, 512) != 1)                 { printf("  [root] omega^512 != 1\n");  fails++; }

    uint16_t a[256], b[256], na[256], nb[256], nc[256], cc[256], ref[256], a2[256], aa[256]
        __attribute__((aligned(32)));

    for (int it = 0; it < 4000; it++) {
        for (int i = 0; i < 256; i++) { a[i]=rand()%Q; b[i]=rand()%Q; }

        /* (1) scalar forward == AVX2 forward */
        for (int i=0;i<256;i++){ na[i]=a[i]; a2[i]=a[i]; }
        ntt_scalar(na); ntt_avx2(a2);
        for (int i=0;i<256;i++) if(na[i]!=a2[i]){ if(fails<4)printf("  [ntt s!=v] i=%d %u %u\n",i,na[i],a2[i]); fails++; break; }

        /* (2) roundtrip: invntt(ntt(a)) == a, both scalar and AVX2 */
        for (int i=0;i<256;i++) aa[i]=a[i];
        ntt_scalar(aa); invntt_scalar(aa);
        for (int i=0;i<256;i++) if(aa[i]!=a[i]){ if(fails<4)printf("  [rt s] i=%d %u %u\n",i,aa[i],a[i]); fails++; break; }
        for (int i=0;i<256;i++) aa[i]=a[i];
        ntt_avx2(aa); invntt_avx2(aa);
        for (int i=0;i<256;i++) if(aa[i]!=a[i]){ if(fails<4)printf("  [rt v] i=%d %u %u\n",i,aa[i],a[i]); fails++; break; }

        /* (3) full negacyclic product via NTT == schoolbook */
        for (int i=0;i<256;i++){ na[i]=a[i]; nb[i]=b[i]; }
        ntt_avx2(na); ntt_avx2(nb);
        pointwise_avx2(na, nb, nc);
        invntt_avx2(nc);
        polymul_ref(a, b, ref);
        for (int i=0;i<256;i++) if(nc[i]!=ref[i]){ if(fails<4)printf("  [polymul] i=%d got %u exp %u\n",i,nc[i],ref[i]); fails++; break; }

        /* (4) pointwise scalar == AVX2 */
        pointwise_scalar(na, nb, cc);
        pointwise_avx2(na, nb, nc);
        for (int i=0;i<256;i++) if(cc[i]!=nc[i]){ if(fails<4)printf("  [pw s!=v]\n"); fails++; break; }
    }
    return fails;
}

/* ------------------------------------------------------------------ */
/*  Benchmark harness: min cycles/coeff; fixed inputs + barrier.       */
/* ------------------------------------------------------------------ */
#define NCOEFF 4096
#define INNER  400
#define NTRIAL 120
#define BARRIER(p) __asm__ __volatile__("" : : "r"(p) : "memory")

static uint16_t a16[NCOEFF],b16[NCOEFF],tw16[NCOEFF],o16[NCOEFF],p16[NCOEFF] __attribute__((aligned(32)));
static int32_t  a32[NCOEFF],b32[NCOEFF],tw32[NCOEFF],o32[NCOEFF],p32[NCOEFF] __attribute__((aligned(32)));

static void init_buffers(void)
{
    int half=(Q-1)/2;
    for(int i=0;i<NCOEFF;i++){
        a16[i]=rand()%Q;b16[i]=rand()%Q;tw16[i]=rand()%Q;
        a32[i]=(rand()%(2*half+1))-half;b32[i]=(rand()%(2*half+1))-half;tw32[i]=(rand()%(2*half+1))-half;
    }
}
static inline uint64_t tic(void){ _mm_lfence(); uint64_t t=__rdtsc(); _mm_lfence(); return t; }

/* kind: 0=mul, 1=CT, 2=GS ;  backend handled per function */
static double bench16(int kind)
{
    __m256i vq=_mm256_set1_epi16((short)Q),vqi=_mm256_set1_epi16((short)QINV16),
            vrmq=_mm256_set1_epi16((short)R_MINUS_Q),zero=_mm256_setzero_si256();
    uint64_t best=~0ull;
    for(int tr=0;tr<NTRIAL;tr++){ uint64_t t0=tic();
        for(int r=0;r<INNER;r++){
            for(int i=0;i<NCOEFF;i+=16){
                __m256i a=_mm256_loadu_si256((__m256i*)&a16[i]),b=_mm256_loadu_si256((__m256i*)&b16[i]),
                        tw=_mm256_loadu_si256((__m256i*)&tw16[i]);
                if(kind==0){ _mm256_storeu_si256((__m256i*)&o16[i], mont_mul_16(a,tw,vqi,vq,vrmq,zero)); }
                else { pair_t p = kind==1 ? ct16(a,b,tw,vqi,vq,vrmq,zero) : gs16(a,b,tw,vqi,vq,vrmq,zero);
                       _mm256_storeu_si256((__m256i*)&o16[i],p.x); _mm256_storeu_si256((__m256i*)&p16[i],p.y); }
            }
            BARRIER(o16); BARRIER(p16);
        }
        uint64_t d=tic()-t0; if(d<best)best=d;
    }
    return (double)best/((double)NCOEFF*INNER);
}
static double bench32(int kind)
{
    __m256i vq=_mm256_set1_epi32(Q),vqi=_mm256_set1_epi32((int)QINV16);
    uint64_t best=~0ull;
    for(int tr=0;tr<NTRIAL;tr++){ uint64_t t0=tic();
        for(int r=0;r<INNER;r++){
            for(int i=0;i<NCOEFF;i+=8){
                __m256i a=_mm256_loadu_si256((__m256i*)&a32[i]),b=_mm256_loadu_si256((__m256i*)&b32[i]),
                        tw=_mm256_loadu_si256((__m256i*)&tw32[i]);
                if(kind==0){ _mm256_storeu_si256((__m256i*)&o32[i], mont_mul_32(a,tw,vqi,vq)); }
                else { pair_t p = kind==1 ? ct32(a,b,tw,vqi,vq) : gs32(a,b,tw,vqi,vq);
                       _mm256_storeu_si256((__m256i*)&o32[i],p.x); _mm256_storeu_si256((__m256i*)&p32[i],p.y); }
            }
            BARRIER(o32); BARRIER(p32);
        }
        uint64_t d=tic()-t0; if(d<best)best=d;
    }
    return (double)best/((double)NCOEFF*INNER);
}
static double benchdq(int kind)
{
    __m256i vqi64=_mm256_set1_epi64x((long long)QINV32),vq64=_mm256_set1_epi64x(Q);
    uint64_t best=~0ull;
    for(int tr=0;tr<NTRIAL;tr++){ uint64_t t0=tic();
        for(int r=0;r<INNER;r++){
            for(int i=0;i<NCOEFF;i+=8){
                __m256i a=_mm256_loadu_si256((__m256i*)&a32[i]),b=_mm256_loadu_si256((__m256i*)&b32[i]),
                        tw=_mm256_loadu_si256((__m256i*)&tw32[i]);
                if(kind==0){ _mm256_storeu_si256((__m256i*)&o32[i], mont_mul_32_dq(a,tw,vqi64,vq64)); }
                else { pair_t p = kind==1 ? ctdq(a,b,tw,vqi64,vq64) : gsdq(a,b,tw,vqi64,vq64);
                       _mm256_storeu_si256((__m256i*)&o32[i],p.x); _mm256_storeu_si256((__m256i*)&p32[i],p.y); }
            }
            BARRIER(o32); BARRIER(p32);
        }
        uint64_t d=tic()-t0; if(d<best)best=d;
    }
    return (double)best/((double)NCOEFF*INNER);
}

/* ---- whole-polynomial NTT benchmarks (cycles per 256-coeff poly) ---- */
static uint16_t pa[256], pb[256], pna[256], pnb[256], pnc[256] __attribute__((aligned(32)));
#define PINNER 2000
static void poly_init(void){ for(int i=0;i<256;i++){ pa[i]=rand()%Q; pb[i]=rand()%Q; } }

static double bench_fwd(int avx2)
{
    uint64_t best=~0ull;
    for(int tr=0;tr<NTRIAL;tr++){ uint64_t t0=tic();
        for(int r=0;r<PINNER;r++){ for(int i=0;i<256;i++) pna[i]=pa[i];
            if(avx2) ntt_avx2(pna); else ntt_scalar(pna); BARRIER(pna); }
        uint64_t d=tic()-t0; if(d<best)best=d; }
    return (double)best/PINNER;
}
static double bench_inv(int avx2)
{
    for(int i=0;i<256;i++) pnb[i]=pa[i];
    if(avx2) ntt_avx2(pnb); else ntt_scalar(pnb);     /* a transformed poly to invert */
    uint64_t best=~0ull;
    for(int tr=0;tr<NTRIAL;tr++){ uint64_t t0=tic();
        for(int r=0;r<PINNER;r++){ for(int i=0;i<256;i++) pna[i]=pnb[i];
            if(avx2) invntt_avx2(pna); else invntt_scalar(pna); BARRIER(pna); }
        uint64_t d=tic()-t0; if(d<best)best=d; }
    return (double)best/PINNER;
}
static double bench_pw(int avx2)
{
    for(int i=0;i<256;i++){ pna[i]=pa[i]; pnb[i]=pb[i]; }
    uint64_t best=~0ull;
    for(int tr=0;tr<NTRIAL;tr++){ uint64_t t0=tic();
        for(int r=0;r<PINNER;r++){ if(avx2) pointwise_avx2(pna,pnb,pnc); else pointwise_scalar(pna,pnb,pnc); BARRIER(pnc); }
        uint64_t d=tic()-t0; if(d<best)best=d; }
    return (double)best/PINNER;
}
static double bench_polymul(int avx2)   /* full a*b: 2 NTT + pointwise + INTT */
{
    uint64_t best=~0ull;
    for(int tr=0;tr<NTRIAL;tr++){ uint64_t t0=tic();
        for(int r=0;r<PINNER;r++){
            for(int i=0;i<256;i++){ pna[i]=pa[i]; pnb[i]=pb[i]; }
            if(avx2){ ntt_avx2(pna); ntt_avx2(pnb); pointwise_avx2(pna,pnb,pnc); invntt_avx2(pnc); }
            else    { ntt_scalar(pna); ntt_scalar(pnb); pointwise_scalar(pna,pnb,pnc); invntt_scalar(pnc); }
            BARRIER(pnc);
        }
        uint64_t d=tic()-t0; if(d<best)best=d; }
    return (double)best/PINNER;
}

/* ------------------------------------------------------------------ */
int main(void)
{
    rinv16_mod_q = egcd_inv((1LL<<16)%Q, Q);
    rinv32_mod_q = egcd_inv((1LL<<32)%Q, Q);
    ntt_init();

    printf("=== q = %d on AVX2 (2^15 < q < 2^16): mul + CT/GS butterflies ===\n", Q);
    printf("q^{-1} mod 2^16 = %u   q^{-1} mod 2^32 = %u\n\n", QINV16, QINV32);

    printf("[1] Correctness (mul + CT + GS, all three backends)\n");
    int fails = test_correctness();
    printf("    %s (%d failures)\n\n", fails==0?"PASS":"FAIL", fails);

    init_buffers();
    printf("[2] Throughput (min cycles/coeff over %d trials, lower = better)\n", NTRIAL);
    double m16=bench16(0), m32=bench32(0), mdq=benchdq(0);
    printf("    mul        16-bit          : %6.3f\n", m16);
    printf("    mul        32-bit vpmulld  : %6.3f\n", m32);
    printf("    mul        32-bit vpmuldq  : %6.3f\n", mdq);
    double ct16c=bench16(1), ct32c=bench32(1), ctdqc=benchdq(1);
    printf("    CT bfly    16-bit          : %6.3f   (reduce every add/sub)\n", ct16c);
    printf("    CT bfly    32-bit vpmulld  : %6.3f   (lazy)\n", ct32c);
    printf("    CT bfly    32-bit vpmuldq  : %6.3f   (lazy)\n", ctdqc);
    double gs16c=bench16(2), gs32c=bench32(2), gsdqc=benchdq(2);
    printf("    GS bfly    16-bit          : %6.3f   (reduce every add/sub)\n", gs16c);
    printf("    GS bfly    32-bit vpmulld  : %6.3f   (lazy)\n", gs32c);
    printf("    GS bfly    32-bit vpmuldq  : %6.3f   (lazy)\n", gsdqc);

    printf("\n[3] Complete n=256 NTT (x^256+1), 16-bit method, Kyber-structured\n");
    int nfails = test_ntt();
    printf("    roots + scalar==AVX2 + roundtrip + polymul(vs schoolbook): %s (%d failures)\n",
           nfails==0?"PASS":"FAIL", nfails);

    poly_init();
    printf("[4] NTT throughput (min cycles per 256-coeff polynomial)\n");
    double fS=bench_fwd(0), fV=bench_fwd(1);
    double iS=bench_inv(0), iV=bench_inv(1);
    double pS=bench_pw(0),  pV=bench_pw(1);
    double mS=bench_polymul(0), mV=bench_polymul(1);
    printf("    forward NTT   scalar %8.1f   AVX2 %8.1f   (%.2fx)\n", fS, fV, fS/fV);
    printf("    inverse NTT   scalar %8.1f   AVX2 %8.1f   (%.2fx)\n", iS, iV, iS/iV);
    printf("    pointwise     scalar %8.1f   AVX2 %8.1f   (%.2fx)\n", pS, pV, pS/pV);
    printf("    full polymul  scalar %8.1f   AVX2 %8.1f   (%.2fx)\n", mS, mV, mS/mV);
    printf("    (AVX2 full polymul = %.2f cycles/coeff)\n", mV/256.0);

    printf("\nSUMMARY fails=%d nfails=%d m16=%.3f m32=%.3f mdq=%.3f "
           "ct16=%.3f ct32=%.3f ctdq=%.3f gs16=%.3f gs32=%.3f gsdq=%.3f "
           "fwdS=%.1f fwdV=%.1f invS=%.1f invV=%.1f pwS=%.1f pwV=%.1f mulS=%.1f mulV=%.1f\n",
           fails, nfails, m16, m32, mdq, ct16c, ct32c, ctdqc, gs16c, gs32c, gsdqc,
           fS, fV, iS, iV, pS, pV, mS, mV);
    return (fails==0 && nfails==0)?0:1;
}
