/* bench.c -- historical B=32 correctness + throughput + constant-time check
 * for three exp(-a) kernels used in docs/design-notes/AVXExp.md: (1) scalar integer
 * ApproxExp   (mulh64, Q51)        -- the AVX2-can't-vectorize baseline
 *   (2) AVX2 double, 4-wide        (binary64 minimax)   -- the FP route
 * (Algorithm 1) (3) AVX-512 IFMA integer, 8-wide (vpmadd52, Q51)    -- the
 * integer-SIMD route (Section 8)
 *
 * All three use the SAME away-from-0 range reduction (B=32, interval
 * [31/32 ln2, ln2]) so they compute the same exp(-a); only the inner
 * arithmetic differs.  This remains useful for the IFMA path and historical
 * comparison; the production AVX2Exp code has since moved to B=8/degree-6.
 *
 * Build:  gcc -O2 -mavx2 -mfma -mavx512f -mavx512dq -mavx512ifma -o bench
 * bench.c -lm Run  :  ./bench Target: any CPU with AVX-512 IFMA (e.g.
 * Intel the recorded CPU generation / Ice Lake, AMD Zen4).
 */
#include <immintrin.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ---- constants (from tools/run_polyapprox.sh & run_avxexp.sh) ---- */
static const uint64_t CQ[6] = {/* Q51 integer coeffs, true = CQ[i]*2^-51 */
                               2251231781277982ULL, 2256720101499690ULL,
                               1108248280973626ULL, 408721085595457ULL,
                               58973132430121ULL,   37125483605256ULL};
static const double CD[6] = {/* binary64 minimax coeffs */
                             0.999747725807801690, 1.002185155330213240,
                             0.492160834701707961, 0.181508945025961682,
                             0.026189149824962118, 0.016487064382500742};
static uint64_t TQ[32]; /* TQ[j] = round(2^{-j/32} * 2^51) */
static double TD[32];   /* TD[j] = 2^{-j/32}              */
static const double LN2 = 0.69314718055994530942;
static void init_tables(void)
{
    for (int j = 0; j < 32; j++) {
        TD[j] = pow(2.0, -(double)j / 32.0);
        TQ[j] = (uint64_t)llround(ldexp(TD[j], 51));
    }
}

/* ===== (1) scalar integer (mulh64, Q51) ===== */
static inline uint64_t mulsh(uint64_t a, uint64_t b, int s)
{
    return (uint64_t)(((unsigned __int128)a * b) >> s);
}
static double exp_scalar(double a)
{
    const int B = 32, Q = 51;
    long m = (long)ceil(B * a / LN2) + (B - 1);
    long n = m >> 5, j = m & 31;
    double rl = ((double)m / B) * LN2 - a;
    uint64_t X = (uint64_t)llround(ldexp(rl, Q));
    uint64_t p = CQ[5];
    for (int i = 4; i >= 0; i--)
        p = mulsh(p, X, Q) + CQ[i];
    uint64_t ps = mulsh(p, TQ[j], Q);
    return ldexp((double)ps, -Q - (int)n);
}

/* ===== (2) AVX2 double, 4-wide ===== */
static void exp_avx2_4(const double *a, double *out)
{
    const double LB = LN2 / 32.0, INV = 32.0 / LN2;
    __m256d va = _mm256_loadu_pd(a);
    __m256d md = _mm256_add_pd(
        _mm256_ceil_pd(_mm256_mul_pd(va, _mm256_set1_pd(INV))),
        _mm256_set1_pd(31.0));
    __m128i m = _mm256_cvtpd_epi32(md);
    __m128i n = _mm_srli_epi32(m, 5),
            j = _mm_and_si128(m, _mm_set1_epi32(31));
    __m256d rl = _mm256_fmsub_pd(md, _mm256_set1_pd(LB), va);
    __m256d p = _mm256_set1_pd(CD[5]);
    for (int i = 4; i >= 0; i--)
        p = _mm256_fmadd_pd(p, rl, _mm256_set1_pd(CD[i]));
    __m256d s = _mm256_setzero_pd();
    for (int k = 0; k < 32; k++) {
        __m128i cmp = _mm_cmpeq_epi32(j, _mm_set1_epi32(k));
        __m256d mask = _mm256_castsi256_pd(_mm256_cvtepi32_epi64(cmp));
        s = _mm256_blendv_pd(s, _mm256_set1_pd(TD[k]), mask);
    }
    __m256i n64 = _mm256_cvtepi32_epi64(n);
    __m256d e = _mm256_castsi256_pd(_mm256_slli_epi64(
        _mm256_sub_epi64(_mm256_set1_epi64x(1023), n64), 52));
    _mm256_storeu_pd(out, _mm256_mul_pd(_mm256_mul_pd(p, s), e));
}

/* ===== (3) AVX-512 IFMA integer, 8-wide ===== */
static inline __m512i mulsh52(__m512i p, __m512i x)
{ /* (p*x) >> 51, p,x < 2^52 */
    __m512i lo = _mm512_madd52lo_epu64(_mm512_setzero_si512(), p, x);
    __m512i hi = _mm512_madd52hi_epu64(_mm512_setzero_si512(), p, x);
    return _mm512_or_si512(_mm512_slli_epi64(hi, 1),
                           _mm512_srli_epi64(lo, 51));
}
static void exp_ifma8(const double *a, double *out)
{
    const double LB = LN2 / 32.0, INV = 32.0 / LN2,
                 SC = (double)(1ULL << 51);
    __m512d va = _mm512_loadu_pd(a);
    __m512d md = _mm512_add_pd(
        _mm512_ceil_pd(_mm512_mul_pd(va, _mm512_set1_pd(INV))),
        _mm512_set1_pd(31.0));
    __m512i m = _mm512_cvtpd_epi64(md);
    __m512i n = _mm512_srli_epi64(m, 5),
            j = _mm512_and_si512(m, _mm512_set1_epi64(31));
    __m512d rl = _mm512_fmsub_pd(md, _mm512_set1_pd(LB), va);
    __m512i X = _mm512_cvtpd_epi64(_mm512_mul_pd(rl, _mm512_set1_pd(SC)));
    __m512i p = _mm512_set1_epi64(CQ[5]);
    for (int i = 4; i >= 0; i--)
        p = _mm512_add_epi64(mulsh52(p, X), _mm512_set1_epi64(CQ[i]));
    __m512i s = _mm512_setzero_si512();
    for (int k = 0; k < 32; k++) {
        __mmask8 mk = _mm512_cmpeq_epi64_mask(j, _mm512_set1_epi64(k));
        s = _mm512_mask_blend_epi64(mk, s, _mm512_set1_epi64(TQ[k]));
    }
    __m512i ps = mulsh52(p, s);
    __m512d e = _mm512_castsi512_pd(_mm512_slli_epi64(
        _mm512_sub_epi64(_mm512_set1_epi64(1023 - 51), n), 52));
    _mm512_storeu_pd(out, _mm512_mul_pd(_mm512_cvtepi64_pd(ps), e));
}

/* ---- helpers ---- */
static inline uint64_t rdtsc(void)
{
    _mm_lfence();
    unsigned x;
    uint64_t t = __rdtscp(&x);
    _mm_lfence();
    return t;
}
static double now_ns(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec * 1e9 + t.tv_nsec;
}
static int cmp_u64(const void *x, const void *y)
{
    uint64_t a = *(const uint64_t *)x, b = *(const uint64_t *)y;
    return a < b ? -1 : a > b;
}

#define AMAX 7.369
static double frand_a(void)
{
    return (double)rand() / RAND_MAX * AMAX;
}

int main(void)
{
    init_tables();
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);

    /* ---------- correctness vs libm ---------- */
    double ws = 0, wa = 0, wi = 0, wsi = 0;
    for (double a = 0; a <= AMAX; a += 0.0002) {
        double ref = exp(-a);
        double in4[4], o4[4], in8[8], o8[8];
        for (int k = 0; k < 4; k++)
            in4[k] = a;
        for (int k = 0; k < 8; k++)
            in8[k] = a;
        exp_avx2_4(in4, o4);
        exp_ifma8(in8, o8);
        double sc = exp_scalar(a);
        double rs = fabs((sc - ref) / ref), ra = fabs((o4[0] - ref) / ref),
               ri = fabs((o8[0] - ref) / ref);
        if (rs > ws)
            ws = rs;
        if (ra > wa)
            wa = ra;
        if (ri > wi)
            wi = ri;
        /* IFMA 8-wide vs scalar reference: same Q51 fixed-point core;
           differ only by the last-bit rounding of the (double) range
           reduction => agree to ~2^-49, not bit-exact */
        double rsi = fabs((o8[0] - sc) / sc);
        if (rsi > wsi)
            wsi = rsi;
    }
    printf(
        "== correctness (max relative error over a in [0,%.3f] vs libm "
        "exp) ==\n",
        AMAX);
    printf("  scalar int (mulh64 Q51) : 2^%.2f\n", log2(ws));
    printf("  AVX2 double  (4-wide)   : 2^%.2f\n", log2(wa));
    printf("  AVX-512 IFMA (8-wide)   : 2^%.2f\n", log2(wi));
    printf("  target                  : 2^-43  (all must beat this)\n");
    printf(
        "  IFMA(8) vs scalar ref   : max rel diff 2^%.1f (same Q51 core; "
        "differ only by range-redux rounding)\n\n",
        log2(wsi));

    /* ---------- throughput ---------- */
    const int N = 8192, K = 2000, R = 21;
    double *A = aligned_alloc(64, N * sizeof(double)),
           *O = aligned_alloc(64, N * sizeof(double));
    for (int i = 0; i < N; i++)
        A[i] = frand_a();
    volatile double sink = 0;
    uint64_t cyc[3][R];
    double ns[3];
    for (int kern = 0; kern < 3; kern++) {
        double best_ns = 1e30;
        for (int r = 0; r < R; r++) {
            uint64_t c0 = rdtsc();
            double t0 = now_ns();
            for (int it = 0; it < K; it++) {
                if (kern == 0) {
                    for (int i = 0; i < N; i++)
                        O[i] = exp_scalar(A[i]);
                } else if (kern == 1) {
                    for (int i = 0; i < N; i += 4)
                        exp_avx2_4(A + i, O + i);
                } else {
                    for (int i = 0; i < N; i += 8)
                        exp_ifma8(A + i, O + i);
                }
                sink += O[it & (N - 1)];
            }
            double t1 = now_ns();
            uint64_t c1 = rdtsc();
            cyc[kern][r] = c1 - c0;
            double e = (t1 - t0) / ((double)N * K);
            if (e < best_ns)
                best_ns = e;
        }
        ns[kern] = best_ns;
    }
    for (int kern = 0; kern < 3; kern++)
        qsort(cyc[kern], R, sizeof(uint64_t), cmp_u64);
    const char *nm[3] = {"scalar int (1x lane)", "AVX2 double (4 lanes)",
                         "AVX-512 IFMA (8 lanes)"};
    printf(
        "== throughput (best of %d runs, %d elems x %d iters; FTZ/DAZ on) "
        "==\n",
        R, N, K);
    printf("  %-24s %10s %10s %10s %8s\n", "kernel", "ns/exp", "Mexp/s",
           "tsc/exp", "speedup");
    double base_ns = ns[0];
    for (int kern = 0; kern < 3; kern++) {
        double tscpe = (double)cyc[kern][R / 2] / ((double)N * K);
        printf("  %-24s %10.3f %10.1f %10.2f %7.2fx\n", nm[kern], ns[kern],
               1e3 / ns[kern], tscpe, base_ns / ns[kern]);
    }
    printf(
        "  (tsc/exp = TSC reference cycles ~3.6GHz; core cycles are fewer "
        "at boost.)\n\n");

    /* ---------- dudect-style constant-time test (fixed vs random inputs)
     * ---------- */
    /* time a batch of BATCH elems; class0 = all a=a_max (fixed), class1 =
     * random.   */
    /* No data-dependent branch/memory => expect |t| small. */
    const int BATCH = 512, M = 300000;
    /* Pre-fill BOTH input buffers ONCE (no writes inside the timed loop,
     * so the two classes are perfectly symmetric apart from the operand
     * values themselves).   */
    double *bf = aligned_alloc(64, BATCH * sizeof(double)),
           *br = aligned_alloc(64, BATCH * sizeof(double)),
           *bo = aligned_alloc(64, BATCH * sizeof(double));
    for (int i = 0; i < BATCH; i++) {
        bf[i] = AMAX;
        br[i] = frand_a();
    }
    for (int tcase = 1; tcase < 3; tcase++) { /* 1=AVX2, 2=IFMA */
        /* warmup */
        for (int w = 0; w < 2000; w++) {
            if (tcase == 1) {
                for (int i = 0; i < BATCH; i += 4)
                    exp_avx2_4(bf + i, bo + i);
            } else {
                for (int i = 0; i < BATCH; i += 8)
                    exp_ifma8(bf + i, bo + i);
            }
        }
        uint64_t *t0s = malloc(M * sizeof(uint64_t)),
                 *t1s = malloc(M * sizeof(uint64_t));
        long c0 = 0, c1 = 0;
        for (int it = 0; it < M; it++) {
            int cls = rand() & 1;
            double *B =
                cls ? br : bf; /* only a pointer select, no writes */
            uint64_t a = rdtsc();
            if (tcase == 1) {
                for (int i = 0; i < BATCH; i += 4)
                    exp_avx2_4(B + i, bo + i);
            } else {
                for (int i = 0; i < BATCH; i += 8)
                    exp_ifma8(B + i, bo + i);
            }
            uint64_t d = rdtsc() - a;
            sink += bo[0];
            if (cls)
                t1s[c1++] = d;
            else
                t0s[c0++] = d;
        }
        /* crop to the interrupt-free lower 50% of each class, then Welch t
         */
        qsort(t0s, c0, sizeof(uint64_t), cmp_u64);
        qsort(t1s, c1, sizeof(uint64_t), cmp_u64);
        long u0 = c0 / 2, u1 = c1 / 2;
        double m0 = 0, m1 = 0, v0 = 0, v1 = 0;
        for (long i = 0; i < u0; i++)
            m0 += t0s[i];
        m0 /= u0;
        for (long i = 0; i < u1; i++)
            m1 += t1s[i];
        m1 /= u1;
        for (long i = 0; i < u0; i++) {
            double e = t0s[i] - m0;
            v0 += e * e;
        }
        v0 /= u0;
        for (long i = 0; i < u1; i++) {
            double e = t1s[i] - m1;
            v1 += e * e;
        }
        v1 /= u1;
        double t = (m0 - m1) / sqrt(v0 / u0 + v1 / u1);
        printf(
            "== constant-time (dudect-style, %s; %d batches of %d, "
            "lower-50%% crop) ==\n",
            tcase == 1 ? "AVX2 double" : "AVX-512 IFMA", M, BATCH);
        printf(
            "  fixed mean=%.1f  random mean=%.1f  (tsc/batch)   Welch |t| "
            "= %.2f  -> %s\n",
            m0, m1, fabs(t),
            fabs(t) < 4.5 ? "no leakage detected (|t|<4.5)"
                          : "|t|>=4.5 (see note)");
        free(t0s);
        free(t1s);
    }
    printf("\n(sink=%.3g)\n", (double)sink);
    free(A);
    free(O);
    free(bf);
    free(br);
    free(bo);
    return 0;
}
