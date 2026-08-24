/*
 * test_ntt.c -- CORRECTNESS tests for the AVX2 assembly NTT/INTT/basemul (ntt.S)
 * of the INCOMPLETE n=256 negacyclic transform over Z_q[x]/(x^256+1), q=57601.
 *
 * INCOMPLETE: the NTT stops at len=2; the full negacyclic product is computed as
 *   a*b  ==  invntt_tomont( basemul( ntt(a), ntt(b) ) )
 * where basemul multiplies per quadratic factor (x^2 - gamma_i), NOT pointwise.
 *
 * Correctness is pinned end-to-end, order-agnostically, by:
 *   [0] scalar-ref self-consistency: invntt_tomont_ref(ntt_ref(a)) == a*R, and the
 *       scalar basemul convolution == schoolbook (validates the oracle).
 *   [1] AVX round-trip:  invntt_tomont_avx(ntt_avx(a)) == a*R
 *   [2] full negacyclic product via the AVX path == schoolbook convolution.
 *
 * THE q>2^15 TRAP: coefficients live in [0,q) but are stored in int16_t, so
 * values in (2^15, q) are NEGATIVE as int16.  Because 2^16 is NOT congruent to 0
 * mod q, read them back as UNSIGNED (uint16_t).
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
static void AVX_BASEMUL(int16_t *c, const int16_t *a, const int16_t *b){
    basemul_avx512(c, ntt512_qdata, a, b, ntt512_basemul_gamma);
}
#define AVX_UNPACK nttunpack_avx512
#else
#define VEC_NAME "AVX2"
static void AVX_NTT(int16_t *p){ ntt_avx   (p, ntt_qdata, ntt_zetas_fwd, ntt_z0,    ntt_scale); }
static void AVX_INV(int16_t *p){ invntt_tomont_avx(p, ntt_qdata, ntt_zetas_inv, ntt_z0inv, ntt_scale); }
static void AVX_BASEMUL(int16_t *c, const int16_t *a, const int16_t *b){
    basemul_avx(c, ntt_qdata, a, b, ntt_basemul_gamma);
}
#define AVX_UNPACK nttunpack_avx
#endif

int main(void){
    ntt_ref_init();
    srand(20260601);

    /* [0] scalar reference self-check (oracle validation) */
    int f0=0;
    for(int it=0; it<2000; it++){
        uint16_t a[N], a0[N], ar[N], br[N], na[N], nb[N], nc[N], cref[N];
        for(int i=0;i<N;i++){ a[i]=rand()%Q; a0[i]=a[i]; ar[i]=rand()%Q; br[i]=rand()%Q; }
        memcpy(na,a0,sizeof na); ntt_ref(na); memcpy(a,na,sizeof a); invntt_tomont_ref(a);
        for(int i=0;i<N;i++) if(a[i]!=(uint16_t)((uint32_t)a0[i]*MONT%Q)){ f0++; break; }
        memcpy(na,ar,sizeof na); memcpy(nb,br,sizeof nb); ntt_ref(na); ntt_ref(nb);
        basemul_ref(na,nb,nc); invntt_tomont_ref(nc);
        polymul_schoolbook(ar,br,cref);
        for(int i=0;i<N;i++) if(nc[i]!=cref[i]){ f0++; break; }
    }
    printf("[0] scalar ref self-check (rt-to-mont + basemul==schoolbook), 2000: %s (%d fails)\n", f0?"FAIL":"PASS", f0);

    /* [1] AVX round trip: invntt_tomont(ntt(a)) == a*R */
    int f1=0;
    for(int it=0; it<5000; it++){
        int16_t a[N]; uint16_t ref[N];
        for(int i=0;i<N;i++){ int v=rand()%Q; a[i]=(int16_t)v; ref[i]=(uint16_t)((uint32_t)v*MONT%Q); }
        AVX_NTT(a); AVX_INV(a);
        for(int i=0;i<N;i++) if(U(a[i])!=ref[i]){ if(f1<6)printf("  [rt] it=%d i=%d got=%u exp=%u\n",it,i,U(a[i]),ref[i]); f1++; break; }
    }
    printf("[1] %s round trip invntt_tomont(ntt(a)) == a*R, 5000 polys: %s (%d fails)\n", VEC_NAME, f1?"FAIL":"PASS", f1);

    /* [2] full negacyclic product via AVX ntt->basemul->invntt_tomont == schoolbook */
    int f2=0;
    for(int it=0; it<2000; it++){
        uint16_t ar[N], br[N], cref[N];
        int16_t na[N], nb[N], nc[N];
        for(int i=0;i<N;i++){ ar[i]=rand()%Q; br[i]=rand()%Q; na[i]=(int16_t)ar[i]; nb[i]=(int16_t)br[i]; }
        AVX_NTT(na); AVX_NTT(nb);
        AVX_BASEMUL(nc, na, nb);
        AVX_INV(nc);
        polymul_schoolbook(ar,br,cref);
        for(int i=0;i<N;i++) if(U(nc[i])!=cref[i]){ if(f2<6)printf("  [mul] it=%d i=%d got=%u exp=%u\n",it,i,U(nc[i]),cref[i]); f2++; break; }
    }
    printf("[2] negacyclic product (%s ntt->basemul->invntt_tomont) == schoolbook, 2000 pairs: %s (%d fails)\n", VEC_NAME, f2?"FAIL":"PASS", f2);

    /* [3] nttunpack: standard-order samples -> avx layout matches AVX ntt.
     * A0 seed-to-NTT relies on nttunpack(ntt_ref(p)) == ntt_avx(p): the scalar
     * ntt_ref output is the standard-order NTT representation, and nttunpack
     * lays it out in the avx NTT slot order (Kyber gen_matrix / poly_nttunpack). */
    int f3=0;
    for(int it=0; it<5000; it++){
        uint16_t p[N]; int32_t src[N]; int16_t want[N], got[N];
        for(int i=0;i<N;i++){ uint16_t v=rand()%Q; p[i]=v; want[i]=(int16_t)v; }
        ntt_ref(p);                              /* standard-order NTT representation */
        for(int i=0;i<N;i++) src[i]=p[i];        /* samples as int32 in [0,q) */
        AVX_NTT(want);                           /* avx-layout NTT representation (may be lazy [0,2q)) */
        AVX_UNPACK(got, src);                    /* standard -> avx layout */
        for(int i=0;i<N;i++) if(U(got[i])%Q!=U(want[i])%Q){ if(f3<6)printf("  [unpack] it=%d i=%d got=%u exp=%u\n",it,i,U(got[i])%Q,U(want[i])%Q); f3++; break; }
    }
    printf("[3] %s nttunpack(ntt_ref(p)) == ntt(p), 5000 polys: %s (%d fails)\n", VEC_NAME, f3?"FAIL":"PASS", f3);

    int fails = f0 + f1 + f2 + f3;
    printf("\nSUMMARY ntt_asm fails=%d\n", fails);
    return fails?1:0;
}
