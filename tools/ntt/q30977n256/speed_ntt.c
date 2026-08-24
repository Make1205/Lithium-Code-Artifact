/*
 * speed_ntt.c - performance benchmarks for the INCOMPLETE n=256 NTT/INTT/basemul
 * (q=30977, 16-bit UNSIGNED [0,q) + lazy), in the style of the IRS speed-sign harness.
 *
 * IMPORTANT: the timed region contains ONLY the transform under test.  All
 * inputs (and, for the inverse, the pre-transform) are built BEFORE the timing
 * loop -- we pre-fill NTESTS independent buffers and time one transform each, so
 * each sample is pure compute (no memcpy / data construction inside the window).
 *
 * Scalar (ntt_ref) and vector (ntt_avx) paths are benchmarked SEPARATELY.
 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "ntt_ref.h"
#include "cpucycles.h"
#include "speed_print.h"

#define NTESTS 10000

static uint16_t in_u16[NTESTS][NTT_N];
static int16_t  in_i16[NTESTS][NTT_N];
static uint16_t fwd_u16[NTESTS][NTT_N];
static int16_t  fwd_i16[NTESTS][NTT_N];
static int16_t  bm_a[NTESTS][NTT_N], bm_b[NTESTS][NTT_N], bm_c[NTESTS][NTT_N];
#ifdef TEST_AVX512
static int16_t  in512_i16[NTESTS][NTT_N];
static int16_t  fwd512_i16[NTESTS][NTT_N];
static int16_t  bm512_a[NTESTS][NTT_N], bm512_b[NTESTS][NTT_N], bm512_c[NTESTS][NTT_N];
#endif

static uint64_t t[NTESTS];

static void avx_fwd(int16_t *p){ ntt_avx   (p, ntt_qdata, ntt_zetas_fwd, ntt_z0,    ntt_scale); }
static void avx_inv(int16_t *p){ invntt_tomont_avx(p, ntt_qdata, ntt_zetas_inv, ntt_z0inv, ntt_scale); }
#ifdef TEST_AVX512
static void avx512_fwd(int16_t *p){ ntt_avx512   (p, ntt512_qdata, ntt512_zetas_fwd, ntt512_scale); }
static void avx512_inv(int16_t *p){ invntt_tomont_avx512(p, ntt512_qdata, ntt512_zetas_inv, ntt512_scale); }
#endif

int main(void){
    unsigned int i;
    ntt_ref_init();

    srand(20260601);
    for(i=0;i<NTESTS;i++)
        for(int j=0;j<NTT_N;j++){
            int v = rand()%NTT_Q;
            in_u16[i][j]=(uint16_t)v;
            in_i16[i][j]=(int16_t)v;
#ifdef TEST_AVX512
            in512_i16[i][j]=(int16_t)v;
#endif
        }
    for(i=0;i<NTESTS;i++){
        memcpy(fwd_u16[i], in_u16[i], sizeof in_u16[i]); ntt_ref(fwd_u16[i]);
        memcpy(fwd_i16[i], in_i16[i], sizeof in_i16[i]); avx_fwd(fwd_i16[i]);
#ifdef TEST_AVX512
        memcpy(fwd512_i16[i], in512_i16[i], sizeof in512_i16[i]); avx512_fwd(fwd512_i16[i]);
#endif
        /* two pre-transformed operands for the basemul benchmark */
        for(int j=0;j<NTT_N;j++){
            int va = rand()%NTT_Q, vb = rand()%NTT_Q;
            bm_a[i][j]=(int16_t)va; bm_b[i][j]=(int16_t)vb;
#ifdef TEST_AVX512
            bm512_a[i][j]=(int16_t)va; bm512_b[i][j]=(int16_t)vb;
#endif
        }
        avx_fwd(bm_a[i]); avx_fwd(bm_b[i]);
#ifdef TEST_AVX512
        avx512_fwd(bm512_a[i]); avx512_fwd(bm512_b[i]);
#endif
    }

    printf("=== INCOMPLETE NTT/INTT/basemul benchmarks  (n=%d, q=%d, 16-bit unsigned[0,q)+lazy AVX2 vs scalar) ===\n", NTT_N, NTT_Q);
    printf("NTESTS=%d, cycles via rdtsc (overhead-subtracted median/average)\n\n", NTESTS);

    /* ============================ SCALAR (ntt_ref) ============================ */
    for(i=0;i<NTESTS;i++){ t[i]=cpucycles(); ntt_ref(in_u16[i]); }
    print_results("scalar  forward NTT  (ntt_ref):", t, NTESTS);

    for(i=0;i<NTESTS;i++){ t[i]=cpucycles(); invntt_tomont_ref(fwd_u16[i]); }
    print_results("scalar  inverse NTT  (invntt_tomont_ref):", t, NTESTS);

    /* ============================ VECTOR (ntt_avx) ============================ */
    for(i=0;i<NTESTS;i++){ t[i]=cpucycles(); avx_fwd(in_i16[i]); }
    print_results("vector  forward NTT  (ntt_avx):", t, NTESTS);

    for(i=0;i<NTESTS;i++){ t[i]=cpucycles(); avx_inv(fwd_i16[i]); }
    print_results("vector  inverse NTT  (invntt_tomont_avx):", t, NTESTS);

    for(i=0;i<NTESTS;i++){ t[i]=cpucycles(); basemul_avx(bm_c[i], bm_a[i], bm_b[i], ntt_qdata, ntt_basemul_gamma); }
    print_results("vector  basemul      (basemul_avx):", t, NTESTS);

#ifdef TEST_AVX512
    /* ============================ VECTOR (AVX-512 full ZMM) ============================ */
    for(i=0;i<NTESTS;i++){ t[i]=cpucycles(); avx512_fwd(in512_i16[i]); }
    print_results("avx512  forward NTT  (ntt_avx512):", t, NTESTS);

    for(i=0;i<NTESTS;i++){ t[i]=cpucycles(); avx512_inv(fwd512_i16[i]); }
    print_results("avx512  inverse NTT  (invntt_tomont_avx512):", t, NTESTS);

    for(i=0;i<NTESTS;i++){ t[i]=cpucycles(); basemul_avx512(bm512_c[i], bm512_a[i], bm512_b[i], ntt512_qdata, ntt512_basemul_gamma); }
    print_results("avx512  basemul      (basemul_avx512):", t, NTESTS);
#endif

    return 0;
}
