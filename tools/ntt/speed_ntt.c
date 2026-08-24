/*
 * speed_ntt.c - performance benchmarks for the complete n=256 NTT/INTT
 * (q=45569, 16-bit), in the style of the IRS speed-sign harness.
 *
 * Conventions copied from that harness:
 *   - cpucycles()/print_results() infrastructure (median + average, minus the
 *     measured rdtsc overhead);
 *   - NTESTS samples per routine.
 *
 * IMPORTANT: the timed region contains ONLY the transform under test.  All
 * inputs (and, for routines that consume a transformed polynomial, the
 * pre-transform) are constructed BEFORE the timing loop.  Each sample copies a
 * pristine input into the work buffer outside cpucycles()? -- no: even the
 * memcpy is kept out of the timed window by transforming in place over a buffer
 * that is restored between *batches*, not between samples.  To keep each sample
 * pure-compute while still feeding the routine a valid (untransformed) input,
 * we pre-fill NTESTS independent input buffers and time one transform each.
 *
 * Scalar (ntt_ref) and vector (ntt_avx) paths are benchmarked SEPARATELY.
 *
 * Build: see Makefile target `speed_ntt`.
 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "ntt_ref.h"
#include "cpucycles.h"
#include "speed_print.h"

#define NTESTS 10000

/* NTESTS independent input polynomials, built once, never inside a timed loop. */
static uint16_t in_u16[NTESTS][NTT_N];   /* for the scalar path  */
static int16_t  in_i16[NTESTS][NTT_N];   /* for the vector path  */
static uint16_t fwd_u16[NTESTS][NTT_N];  /* pre-transformed (scalar) inputs for invntt_tomont_ref */
static int16_t  fwd_i16[NTESTS][NTT_N];  /* pre-transformed (vector) inputs for invntt_tomont_avx */
static int16_t  pw_a[NTESTS][NTT_N], pw_b[NTESTS][NTT_N], pw_c[NTESTS][NTT_N];

static uint64_t t[NTESTS];

static void avx_fwd(int16_t *p){ ntt_avx   (p, ntt_qdata, ntt_zetas_fwd, ntt_z0,    ntt_ninv); }
static void avx_inv(int16_t *p){ invntt_tomont_avx(p, ntt_qdata, ntt_zetas_inv, ntt_z0inv, ntt_ninv); }
static void avx_pointwise(int16_t *c, const int16_t *a, const int16_t *b){ pointwise_avx(c, a, b, ntt_qdata); }

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
        }
    /* pre-transform copies so the inverse benchmarks time ONLY the inverse */
    for(i=0;i<NTESTS;i++){
        memcpy(fwd_u16[i], in_u16[i], sizeof in_u16[i]); ntt_ref(fwd_u16[i]);
        memcpy(fwd_i16[i], in_i16[i], sizeof in_i16[i]); avx_fwd(fwd_i16[i]);
        for(int j=0;j<NTT_N;j++){
            int va = rand()%NTT_Q, vb = rand()%NTT_Q;
            pw_a[i][j]=(int16_t)va; pw_b[i][j]=(int16_t)vb;
        }
        avx_fwd(pw_a[i]); avx_fwd(pw_b[i]);
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

    return 0;
}
