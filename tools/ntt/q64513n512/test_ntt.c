/*
 * test_ntt.c -- CORRECTNESS tests for the AVX2 assembly NTT/INTT (ntt.S) of the
 * complete n=512 negacyclic transform over Z_q[x]/(x^512+1), q=64513, 16-bit.
 *
 * Correctness is pinned end-to-end, order-agnostically, by:
 *   [1] inverse round-trip:  invntt_tomont(ntt(a)) == a*R
 *   [2] full negacyclic product via the AVX NTT == schoolbook convolution
 *
 * THE q>2^15 TRAP: coefficients live in [0,q) but are stored in int16_t, so
 * values in (2^15, q) are NEGATIVE as int16.  Because 2^16 is NOT congruent to
 * 0 mod q, read them back as UNSIGNED (uint16_t); reducing the signed value mod
 * q gives a wrong answer.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ntt_ref.h"

#define Q NTT_Q
#define N NTT_N
#define MONT ((uint16_t)((1u << 16) % Q))
#define U(x) ((uint16_t)(x))          /* read an int16 lane as its [0,q) value */

#ifdef TEST_AVX512
#define VEC_NAME "AVX-512"
static void AVX_NTT(int16_t *p){ ntt_avx512   (p, ntt512_qdata, ntt512_zetas_fwd, ntt512_scale); }
static void AVX_INV(int16_t *p){ invntt_tomont_avx512(p, ntt512_qdata, ntt512_zetas_inv, ntt512_scale); }
static void AVX_POINTWISE(int16_t *c, const int16_t *a, const int16_t *b){ pointwise_avx512(c, a, b, ntt512_qdata); }
#define AVX_UNPACK nttunpack_avx512
#else
#define VEC_NAME "AVX2"
static void AVX_NTT(int16_t *p){ ntt_avx   (p, ntt_qdata, ntt_zetas_fwd, ntt_cross_fwd, ntt_ninv); }
static void AVX_INV(int16_t *p){ invntt_tomont_avx(p, ntt_qdata, ntt_zetas_inv, ntt_cross_inv, ntt_ninv); }
static void AVX_POINTWISE(int16_t *c, const int16_t *a, const int16_t *b){ pointwise_avx(c, a, b, ntt_qdata); }
#define AVX_UNPACK nttunpack_avx
#endif

int main(void){
    ntt_ref_init();
    srand(20260531);

    /* [1] round trip: invntt_tomont(ntt(a)) == a*R */
    int f1=0;
    for(int it=0; it<5000; it++){
        int16_t a[N]; uint16_t ref[N];
        for(int i=0;i<N;i++){ int v=rand()%Q; a[i]=(int16_t)v; ref[i]=(uint16_t)((uint32_t)v*MONT%Q); }
        AVX_NTT(a); AVX_INV(a);
        for(int i=0;i<N;i++) if(U(a[i])!=ref[i]){ if(f1<6)printf("  [rt] it=%d i=%d got=%u exp=%u\n",it,i,U(a[i]),ref[i]); f1++; break; }
    }
    printf("[1] %s round trip invntt_tomont(ntt(a)) == a*R, 5000 polys: %s (%d fails)\n", VEC_NAME, f1?"FAIL":"PASS", f1);

    /* [2] full negacyclic product via AVX NTT == schoolbook */
    int f2=0;
    for(int it=0; it<2000; it++){
        uint16_t ar[N], br[N], cref[N];
        int16_t na[N], nb[N], nc[N];
        for(int i=0;i<N;i++){ ar[i]=rand()%Q; br[i]=rand()%Q; na[i]=(int16_t)ar[i]; nb[i]=(int16_t)br[i]; }
        AVX_NTT(na); AVX_NTT(nb);
        AVX_POINTWISE(nc, na, nb);
        AVX_INV(nc);
        polymul_schoolbook(ar,br,cref);
        for(int i=0;i<N;i++) if(U(nc[i])!=cref[i]){ if(f2<6)printf("  [mul] it=%d i=%d got=%u exp=%u\n",it,i,U(nc[i]),cref[i]); f2++; break; }
    }
    printf("[2] negacyclic product (%s NTT+pointwise+invntt_tomont) == schoolbook, 2000 pairs: %s (%d fails)\n", VEC_NAME, f2?"FAIL":"PASS", f2);

    /* [3] nttunpack: standard-order samples -> avx layout matches AVX ntt.
     * A0 seed-to-NTT relies on nttunpack(ntt_ref(p)) == ntt_avx(p) (mod q): the
     * scalar ntt_ref output is the standard-order NTT rep, and nttunpack lays it
     * out in the avx NTT slot order (Kyber gen_matrix / poly_nttunpack). */
    int f3=0;
    for(int it=0; it<5000; it++){
        uint16_t p[N]; int32_t src[N]; int16_t want[N], got[N];
        for(int i=0;i<N;i++){ uint16_t v=rand()%Q; p[i]=v; want[i]=(int16_t)v; }
        ntt_ref(p);
        for(int i=0;i<N;i++) src[i]=p[i];
        AVX_NTT(want);
        AVX_UNPACK(got, src);
        for(int i=0;i<N;i++) if(U(got[i])%Q!=U(want[i])%Q){ if(f3<6)printf("  [unpack] it=%d i=%d got=%u exp=%u\n",it,i,U(got[i])%Q,U(want[i])%Q); f3++; break; }
    }
    printf("[3] %s nttunpack(ntt_ref(p)) == ntt(p) (mod q), 5000 polys: %s (%d fails)\n", VEC_NAME, f3?"FAIL":"PASS", f3);

    int fails = f1 + f2 + f3;
    printf("\nSUMMARY ntt_asm fails=%d\n", fails);
    return fails?1:0;
}
