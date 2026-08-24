/*
 * ntt_ref.c -- portable SCALAR reference for the INCOMPLETE n=256 negacyclic
 * NTT / INTT over Z_q[x]/(x^256+1), q=57601, 16-bit Montgomery.
 *
 * INCOMPLETE: q-1 = 2^8*3^2*5^2; n=256 | q-1 but 2n=512 does NOT, so x^256+1
 * splits into 128 quadratic factors (x^2 - gamma_i).  The forward transform runs
 * 7 CT levels (len 128..2) and STOPS at len=2; the pointwise step is a basemul
 * (per quadratic factor), not a coefficient-wise product.
 *
 * This is the textbook (non-vectorized) transform: oracle for the AVX2 assembly
 * and scalar baseline for the speed benchmark.  It mirrors gen_ntt.py exactly
 * (same psi, brv7, gamma, scale = 1/128 conventions).
 *
 * NOTE: the scalar transform leaves outputs in standard bit-reversed order,
 * which is NOT the same storage permutation as the merged/shuffled AVX2 routine.
 * The two agree on the end-to-end results that matter (round-trip and
 * convolution), which is what the test harness checks.
 */
#include "ntt_ref.h"

#define R 65536u

static uint32_t RINV;            /* (2^16)^-1 mod q */
static uint16_t ZMONT[128];      /* psi^{brv7(k)} * R mod q,  k=1..127 */
static uint16_t GAMMA[128];      /* psi^(2*brv7(i)+1) mod q (plain) */
static uint16_t SCALETOMONT;     /* 128^-1 * R^2 mod q */
uint16_t ntt_ref_R2;             /* R^2 mod q (exported for callers) */

static uint32_t powmod(uint32_t b, uint32_t e){
    uint64_t r=1,x=b%NTT_Q; while(e){ if(e&1)r=r*x%NTT_Q; x=x*x%NTT_Q; e>>=1; } return (uint32_t)r;
}
/* multiplicative order of g mod q (factoring q-1), then the SMALLEST primitive
 * `order`-th root -- the same rule gen_zetas.py / gen_ntt.py use, so this oracle
 * matches the avx assembly's (ref-aligned) zeta layout. */
static uint32_t mod_order(uint32_t g){
    uint32_t n=NTT_Q-1, m=NTT_Q-1;
    for(uint32_t p=2; (uint64_t)p*p<=m; p++) if(m%p==0){
        while(m%p==0) m/=p;
        while(n%p==0 && powmod(g,n/p)==1) n/=p;
    }
    if(m>1){ while(n%m==0 && powmod(g,n/m)==1) n/=m; }
    return n;
}
static uint32_t prim_root(uint32_t order){
    for(uint32_t g=2; g<NTT_Q; g++) if(mod_order(g)==order) return g;
    return 0;
}
static int brv7(int x){ int r=0; for(int i=0;i<7;i++) r|=((x>>i)&1)<<(6-i); return r; }

uint16_t ntt_ref_fqmul(uint16_t a, uint16_t b){
    return (uint16_t)((uint64_t)a*b%NTT_Q*RINV%NTT_Q);
}
static uint16_t addm(uint16_t a, uint16_t b){ uint32_t s=a+b; return (uint16_t)(s>=NTT_Q?s-NTT_Q:s); }
static uint16_t subm(uint16_t a, uint16_t b){ return (uint16_t)(a>=b?a-b:a+NTT_Q-b); }

void ntt_ref_init(void){
    RINV = powmod(R%NTT_Q, NTT_Q-2);
    uint32_t psi = prim_root(256);                    /* smallest primitive 256th root (=68), psi^128=-1 */
    for(int k=0;k<128;k++) ZMONT[k]=(uint16_t)((uint64_t)powmod(psi,brv7(k))*(R%NTT_Q)%NTT_Q);
    for(int i=0;i<128;i++) GAMMA[i]=(uint16_t)powmod(psi, 2*brv7(i)+1);
    ntt_ref_R2 = (uint16_t)((uint64_t)(R%NTT_Q)*(R%NTT_Q)%NTT_Q);
    SCALETOMONT = (uint16_t)((uint64_t)powmod(128,NTT_Q-2)*ntt_ref_R2%NTT_Q);
}

/* forward: 7 CT levels, len 128..2 (STOP at len=2) */
void ntt_ref(uint16_t r[256]){
    int k=1;
    for(int len=128; len>=2; len>>=1)
        for(int s=0; s<256; s+=2*len){ uint16_t z=ZMONT[k++];
            for(int j=s;j<s+len;j++){ uint16_t t=ntt_ref_fqmul(z,r[j+len]); r[j+len]=subm(r[j],t); r[j]=addm(r[j],t); } }
}
/* inverse: 7 GS levels, len 2..128, then multiply by 1/128 */
void invntt_tomont_ref(uint16_t r[256]){
    int k=127;
    for(int len=2; len<=128; len<<=1)
        for(int s=0; s<256; s+=2*len){ uint16_t z=ZMONT[k--];
            for(int j=s;j<s+len;j++){ uint16_t t=r[j]; r[j]=addm(t,r[j+len]); r[j+len]=ntt_ref_fqmul(z,subm(r[j+len],t)); } }
    for(int j=0;j<256;j++) r[j]=ntt_ref_fqmul(r[j],SCALETOMONT);
}
/* basemul: per quadratic factor (x^2 - gamma_i), standard layout pair i=(2i,2i+1).
 *   c0 = a0*b0 + gamma_i*(a1*b1);  c1 = a0*b1 + a1*b0   (all plain mod q). */
void basemul_ref(const uint16_t a[256], const uint16_t b[256], uint16_t c[256]){
    for(int i=0;i<128;i++){
        uint16_t a0=a[2*i],a1=a[2*i+1],b0=b[2*i],b1=b[2*i+1],g=GAMMA[i];
        uint16_t a0b0=ntt_ref_fqmul(a0,b0);
        uint16_t a1b1=ntt_ref_fqmul(a1,b1);
        uint16_t a0b1=ntt_ref_fqmul(a0,b1);
        uint16_t a1b0=ntt_ref_fqmul(a1,b0);
        c[2*i]   = addm(a0b0, (uint16_t)((uint32_t)g*a1b1%NTT_Q));
        c[2*i+1] = addm(a0b1, a1b0);
    }
}
void polymul_schoolbook(const uint16_t a[256], const uint16_t b[256], uint16_t c[256]){
    for(int i=0;i<256;i++) c[i]=0;
    for(int i=0;i<256;i++) for(int j=0;j<256;j++){ uint16_t p=(uint16_t)((uint32_t)a[i]*b[j]%NTT_Q); int k=i+j;
        if(k<256) c[k]=addm(c[k],p); else c[k-256]=subm(c[k-256],p); }
}
