/*
 * test_ntt.c -- CORRECTNESS tests for the AVX2 assembly NTT/INTT + BASEMUL
 * (ntt.S) of the INCOMPLETE n=256 negacyclic transform over Z_q[x]/(x^256+1),
 * q=30977, 16-bit UNSIGNED [0,q) + LAZY reduction (csub schedule).
 *
 * (Performance is measured separately in speed_ntt.c.)
 *
 * Correctness is pinned end-to-end, order-agnostically, by:
 *   [1] inverse round-trip to Montgomery: invntt_tomont(ntt(a)) == a*R
 *   [2] full negacyclic product via the AVX NTT+BASEMUL == schoolbook
 * Both validate the merged/shuffled data path without assuming any particular
 * mid-transform output permutation.  The scalar ntt_ref provides the schoolbook
 * oracle (polymul_schoolbook) and the constants.
 *
 * READBACK: coefficients live in UNSIGNED [0,2q) (lazy) -> at most 2q-1 = 61953,
 * which exceeds int16 range, so a lane can come back as a NEGATIVE int16 whose
 * true residue is the value mod q.  Because 2^16 is NOT congruent to 0 mod q,
 * reading a lane and taking %q must use the SIGNED value: U(x) below maps the
 * int16 bit pattern back to its residue in [0,q).
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ntt_ref.h"
#include "poly.h"

#define Q NTT_Q
#define N NTT_N
#define MONT ((uint16_t)((1u << 16) % Q))
static inline uint16_t U(int16_t x){ int r = (int)x % Q; if(r<0) r+=Q; return (uint16_t)r; }
static inline uint16_t URAW(int16_t x){ return (uint16_t)x; }

#ifdef TEST_AVX512
#define VEC_NAME "AVX-512"
static void AVX_NTT(int16_t *p){ ntt_avx512   (p, ntt512_qdata, ntt512_zetas_fwd, ntt512_scale); }
static void AVX_INV(int16_t *p){ invntt_tomont_avx512(p, ntt512_qdata, ntt512_zetas_inv, ntt512_scale); }
static void AVX_BMUL(int16_t *c, const int16_t *a, const int16_t *b){ basemul_avx512(c, a, b, ntt512_qdata, ntt512_basemul_gamma); }
#define AVX_UNPACK nttunpack_avx512
#else
#define VEC_NAME "AVX2"
static void AVX_NTT(int16_t *p){ ntt_avx   (p, ntt_qdata, ntt_zetas_fwd, ntt_z0,    ntt_scale); }
static void AVX_INV(int16_t *p){ invntt_tomont_avx(p, ntt_qdata, ntt_zetas_inv, ntt_z0inv, ntt_scale); }
static void AVX_BMUL(int16_t *c, const int16_t *a, const int16_t *b){ basemul_avx(c, a, b, ntt_qdata, ntt_basemul_gamma); }
#define AVX_UNPACK nttunpack_avx
#endif

int main(void){
    ntt_ref_init();
    srand(20260601);

    /* [0] scalar self-consistency: invntt_tomont_ref(ntt_ref(a))==a*R, and the scalar
       basemul reproduces schoolbook -- guards the oracle itself. */
    int f0=0;
    {
        for(int it=0; it<2000; it++){
            uint16_t a[N], aa[N], b[N], nc[N], cref[N];
            for(int i=0;i<N;i++){ a[i]=rand()%Q; b[i]=rand()%Q; }
            memcpy(aa,a,sizeof a); ntt_ref(aa); invntt_tomont_ref(aa);
            for(int i=0;i<N;i++) if(aa[i]!=(uint16_t)((uint32_t)a[i]*MONT%Q)){ f0++; break; }
            uint16_t na[N], nb[N];
            memcpy(na,a,sizeof a); ntt_ref(na);
            memcpy(nb,b,sizeof b); ntt_ref(nb);
            basemul_ref(na,nb,nc); invntt_tomont_ref(nc);
            polymul_schoolbook(a,b,cref);
            for(int i=0;i<N;i++) if(nc[i]!=cref[i]){ f0++; break; }
        }
        printf("[0] scalar ref self-check (rt-to-mont + basemul vs schoolbook): %s (%d fails)\n", f0?"FAIL":"PASS", f0);
    }

    /* [1] round trip: invntt_tomont(ntt(a)) == a*R (AVX) */
    int f1=0;
    for(int it=0; it<5000; it++){
        int16_t a[N]; uint16_t ref[N];
        for(int i=0;i<N;i++){ int v=rand()%Q; a[i]=(int16_t)v; ref[i]=(uint16_t)((uint32_t)v*MONT%Q); }
        AVX_NTT(a); AVX_INV(a);
        for(int i=0;i<N;i++) if(U(a[i])!=ref[i]){ if(f1<6)printf("  [rt] it=%d i=%d got=%u exp=%u\n",it,i,U(a[i]),ref[i]); f1++; break; }
    }
    printf("[1] %s round trip invntt_tomont(ntt(a)) == a*R, 5000 polys: %s (%d fails)\n", VEC_NAME, f1?"FAIL":"PASS", f1);

    /* [2] full negacyclic product via AVX NTT + BASEMUL == schoolbook */
    int f2=0;
    for(int it=0; it<2000; it++){
        uint16_t ar[N], br[N], cref[N];
        int16_t na[N], nb[N], nc[N];
        for(int i=0;i<N;i++){ ar[i]=rand()%Q; br[i]=rand()%Q; na[i]=(int16_t)ar[i]; nb[i]=(int16_t)br[i]; }
        AVX_NTT(na); AVX_NTT(nb);
        AVX_BMUL(nc, na, nb);
        AVX_INV(nc);
        polymul_schoolbook(ar,br,cref);
        for(int i=0;i<N;i++) if(U(nc[i])!=cref[i]){ if(f2<6)printf("  [mul] it=%d i=%d got=%u exp=%u\n",it,i,U(nc[i]),cref[i]); f2++; break; }
    }
    printf("[2] negacyclic product (%s NTT+basemul+invntt_tomont) == schoolbook, 2000 pairs: %s (%d fails)\n", VEC_NAME, f2?"FAIL":"PASS", f2);

    int f3=0;
    {
        static const uint16_t edges[] = {0, Q-1, Q, 2*Q-1, 2*Q, 65535};
        poly p;
        for(size_t i=0;i<sizeof(edges)/sizeof(edges[0]);i++) {
            for(int j=0;j<N;j++) p.coeffs[j] = (int16_t)edges[i];
            poly_reduce(&p);
            uint16_t exp = (uint16_t)(edges[i] % Q);
            for(int j=0;j<N;j++) if(URAW(p.coeffs[j]) != exp){ f3++; break; }
        }
        for(int it=0; it<2000; it++) {
            uint16_t ref[N];
            for(int i=0;i<N;i++){ uint16_t v=(uint16_t)rand(); p.coeffs[i]=(int16_t)v; ref[i]=(uint16_t)(v%Q); }
            poly_reduce(&p);
            for(int i=0;i<N;i++) if(URAW(p.coeffs[i]) != ref[i]){ if(f3<6)printf("  [reduce] it=%d i=%d got=%u exp=%u\n",it,i,URAW(p.coeffs[i]),ref[i]); f3++; break; }
        }
    }
    printf("[3] poly_reduce/reduce_avx canonicalizes uint16 coeffs: %s (%d fails)\n", f3?"FAIL":"PASS", f3);

    /* [4] nttunpack: standard-order samples -> avx layout matches AVX ntt.
     * A0 seed-to-NTT relies on nttunpack(ntt_ref(p)) == ntt_avx(p) (mod q): the
     * scalar ntt_ref output is the standard-order NTT rep, and nttunpack lays it
     * out in the avx NTT slot order (Kyber gen_matrix / poly_nttunpack). */
    int f4=0;
    for(int it=0; it<5000; it++){
        uint16_t p[N]; int32_t src[N]; int16_t want[N], got[N];
        for(int i=0;i<N;i++){ uint16_t v=rand()%Q; p[i]=v; want[i]=(int16_t)v; }
        ntt_ref(p);
        for(int i=0;i<N;i++) src[i]=p[i];
        AVX_NTT(want);                           /* avx-layout NTT rep (LAZY [0,2q): read UNSIGNED) */
        AVX_UNPACK(got, src);                    /* standard -> avx layout */
        for(int i=0;i<N;i++) if(URAW(got[i])%Q!=URAW(want[i])%Q){ if(f4<6)printf("  [unpack] it=%d i=%d got=%u exp=%u\n",it,i,URAW(got[i])%Q,URAW(want[i])%Q); f4++; break; }
    }
    printf("[4] %s nttunpack(ntt_ref(p)) == ntt(p) (mod q), 5000 polys: %s (%d fails)\n", VEC_NAME, f4?"FAIL":"PASS", f4);

    int fails = f0 + f1 + f2 + f3 + f4;
    printf("\nSUMMARY ntt_asm fails=%d\n", fails);
    return fails?1:0;
}
