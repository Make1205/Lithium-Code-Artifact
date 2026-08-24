/*
 * ntt_ref.c -- portable SCALAR reference implementation of the complete n=256
 * negacyclic NTT / INTT over Z_q[x]/(x^256+1), q=45569, 16-bit Montgomery.
 *
 * This is the textbook (non-vectorized) transform: it is used both as the
 * correctness oracle for the AVX2 assembly and as the scalar baseline in the
 * speed benchmark.  It mirrors gen_ntt.py exactly (same twiddle convention).
 *
 * NOTE: the scalar transform leaves outputs in standard bit-reversed order,
 * which is NOT the same storage permutation as the merged/shuffled AVX2
 * routine.  The two agree on the end-to-end results that matter (round-trip
 * and convolution), which is what the test harness checks; they are not
 * expected to match coefficient-by-coefficient mid-transform.
 */
#include "ntt_ref.h"

#define R 65536u

static uint32_t RINV;            /* (2^16)^-1 mod q */
static uint16_t ZMONT[256];      /* omega^{brv8(k)} * R mod q */
static uint16_t NINVTOMONT;      /* 256^-1 * R^2 mod q */
uint16_t ntt_ref_R2;             /* R^2 mod q (exported for diagnostics) */

static uint32_t powmod(uint32_t b, uint32_t e){
    uint64_t r=1,x=b%NTT_Q; while(e){ if(e&1)r=r*x%NTT_Q; x=x*x%NTT_Q; e>>=1; } return (uint32_t)r;
}
static int brv8(int x){ int r=0; for(int i=0;i<8;i++) r|=((x>>i)&1)<<(7-i); return r; }

uint16_t ntt_ref_fqmul(uint16_t a, uint16_t b){
    return (uint16_t)((uint64_t)a*b%NTT_Q*RINV%NTT_Q);
}
static uint16_t addm(uint16_t a, uint16_t b){ uint32_t s=a+b; return (uint16_t)(s>=NTT_Q?s-NTT_Q:s); }
static uint16_t subm(uint16_t a, uint16_t b){ return (uint16_t)(a>=b?a-b:a+NTT_Q-b); }

void ntt_ref_init(void){
    RINV = powmod(R%NTT_Q, NTT_Q-2);
    uint32_t w = powmod(3,(NTT_Q-1)/512);
    for(int k=0;k<256;k++) ZMONT[k]=(uint16_t)((uint64_t)powmod(w,brv8(k))*(R%NTT_Q)%NTT_Q);
    ntt_ref_R2 = (uint16_t)((uint64_t)(R%NTT_Q)*(R%NTT_Q)%NTT_Q);
    NINVTOMONT = (uint16_t)((uint64_t)powmod(256,NTT_Q-2)*ntt_ref_R2%NTT_Q);
}

void ntt_ref(uint16_t r[256]){
    int k=1;
    for(int len=128; len>=1; len>>=1)
        for(int s=0; s<256; s+=2*len){ uint16_t z=ZMONT[k++];
            for(int j=s;j<s+len;j++){ uint16_t t=ntt_ref_fqmul(z,r[j+len]); r[j+len]=subm(r[j],t); r[j]=addm(r[j],t); } }
}
void invntt_tomont_ref(uint16_t r[256]){
    int k=255;
    for(int len=1; len<=128; len<<=1)
        for(int s=0; s<256; s+=2*len){ uint16_t z=ZMONT[k--];
            for(int j=s;j<s+len;j++){ uint16_t t=r[j]; r[j]=addm(t,r[j+len]); r[j+len]=ntt_ref_fqmul(z,subm(r[j+len],t)); } }
    for(int j=0;j<256;j++) r[j]=ntt_ref_fqmul(r[j],NINVTOMONT);
}
void pointwise_ref(const uint16_t a[256], const uint16_t b[256], uint16_t c[256]){
    for(int i=0;i<256;i++) c[i]=ntt_ref_fqmul(a[i], b[i]);
}
void polymul_schoolbook(const uint16_t a[256], const uint16_t b[256], uint16_t c[256]){
    for(int i=0;i<256;i++) c[i]=0;
    for(int i=0;i<256;i++) for(int j=0;j<256;j++){ uint16_t p=(uint16_t)((uint32_t)a[i]*b[j]%NTT_Q); int k=i+j;
        if(k<256) c[k]=addm(c[k],p); else c[k-256]=subm(c[k-256],p); }
}
