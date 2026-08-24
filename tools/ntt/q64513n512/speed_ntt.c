/*
 * speed_ntt.c - performance benchmarks for the complete n=512 NTT/INTT
 * (q=64513, 16-bit), in the the IRS speed-sign harness style.
 *
 * The timed region contains ONLY the transform under test.  All inputs (and the
 * pre-transforms feeding the inverse benchmarks) are constructed BEFORE every
 * timing loop, so each sample is pure compute (no memcpy / no data build).
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

static uint16_t in_u16[NTESTS][NTT_N];   /* scalar path  */
static int16_t  in_i16[NTESTS][NTT_N];   /* vector path  */
static uint16_t fwd_u16[NTESTS][NTT_N];  /* pre-transformed (scalar) for invntt_tomont_ref */
static int16_t  fwd_i16[NTESTS][NTT_N];  /* pre-transformed (vector) for invntt_tomont_avx */
static int16_t  pw_a[NTESTS][NTT_N], pw_b[NTESTS][NTT_N], pw_c[NTESTS][NTT_N];
#ifdef TEST_AVX512
static int16_t  in512_i16[NTESTS][NTT_N];
static int16_t  fwd512_i16[NTESTS][NTT_N];
static int16_t  pw512_a[NTESTS][NTT_N], pw512_b[NTESTS][NTT_N], pw512_c[NTESTS][NTT_N];
#endif

static uint64_t t[NTESTS];

static void avx_fwd(int16_t *p){ ntt_avx   (p, ntt_qdata, ntt_zetas_fwd, ntt_cross_fwd, ntt_ninv); }
static void avx_inv(int16_t *p){ invntt_tomont_avx(p, ntt_qdata, ntt_zetas_inv, ntt_cross_inv, ntt_ninv); }
static void avx_pointwise(int16_t *c, const int16_t *a, const int16_t *b){ pointwise_avx(c, a, b, ntt_qdata); }
#ifdef TEST_AVX512
static void avx512_fwd(int16_t *p){ ntt_avx512   (p, ntt512_qdata, ntt512_zetas_fwd, ntt512_scale); }
static void avx512_inv(int16_t *p){ invntt_tomont_avx512(p, ntt512_qdata, ntt512_zetas_inv, ntt512_scale); }
#endif

int main(void){
    unsigned int i;
    ntt_ref_init();

    /* ---------- build ALL inputs up front (outside every timed region) ---------- */
    srand(20260531);
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
        for(int j=0;j<NTT_N;j++){
            int va = rand()%NTT_Q, vb = rand()%NTT_Q;
            pw_a[i][j]=(int16_t)va; pw_b[i][j]=(int16_t)vb;
#ifdef TEST_AVX512
            pw512_a[i][j]=(int16_t)va; pw512_b[i][j]=(int16_t)vb;
#endif
        }
        avx_fwd(pw_a[i]); avx_fwd(pw_b[i]);
#ifdef TEST_AVX512
        memcpy(fwd512_i16[i], in512_i16[i], sizeof in512_i16[i]); avx512_fwd(fwd512_i16[i]);
        avx512_fwd(pw512_a[i]); avx512_fwd(pw512_b[i]);
#endif
    }

    printf("=== NTT/INTT benchmarks  (n=%d, q=%d, 16-bit AVX2 vs scalar) ===\n", NTT_N, NTT_Q);
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

    for(i=0;i<NTESTS;i++){ t[i]=cpucycles(); avx_pointwise(pw_c[i], pw_a[i], pw_b[i]); }
    print_results("vector  pointwise    (pointwise_avx):", t, NTESTS);

#ifdef TEST_AVX512
    /* ============================ VECTOR (AVX-512 full ZMM) ============================ */
    for(i=0;i<NTESTS;i++){ t[i]=cpucycles(); avx512_fwd(in512_i16[i]); }
    print_results("avx512  forward NTT  (ntt_avx512):", t, NTESTS);

    for(i=0;i<NTESTS;i++){ t[i]=cpucycles(); avx512_inv(fwd512_i16[i]); }
    print_results("avx512  inverse NTT  (invntt_tomont_avx512):", t, NTESTS);

    for(i=0;i<NTESTS;i++){ t[i]=cpucycles(); pointwise_avx512(pw512_c[i], pw512_a[i], pw512_b[i], ntt512_qdata); }
    print_results("avx512  pointwise    (pointwise_avx512):", t, NTESTS);
#endif

    return 0;
}
