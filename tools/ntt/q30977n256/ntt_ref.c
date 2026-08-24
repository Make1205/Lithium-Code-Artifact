/*
 * ntt_ref.c -- portable SCALAR reference for the INCOMPLETE n=256 negacyclic
 * NTT / INTT over Z_q[x]/(x^256+1), q=30977, 16-bit Montgomery.
 *
 * q-1 = 2^8 * 11^2, so 2n=512 does NOT divide q-1 (only n=256 does): the NTT is
 * INCOMPLETE -- 7 butterfly levels (len 128..2), leaving 128 degree-2 residues
 * mod (x^2 - gamma_i), gamma_i = psi^(2*brv7(i)+1).  The pointwise step is a
 * BASEMUL of those 128 pairs, not a coefficient-wise multiply.  Final inverse
 * scale = (2^7)^{-1} = 128^{-1} mod q.
 *
 * This scalar transform is the correctness oracle AND the scalar speed baseline.
 * It is SELF-CONSISTENT (scalar basemul on scalar ntt output) and agrees with
 * the AVX path on the end-to-end results the harness checks (round-trip and full
 * negacyclic product); the two are NOT expected to match coefficient-by-
 * coefficient mid-transform (different storage permutations).
 *
 * Because q < 2^15 every coefficient in [0,q) fits a positive int16, so the
 * usual Montgomery arithmetic over uint16 works without the q>2^15 "valley" tax.
 */
#include "ntt_ref.h"

#define R 65536u

static uint32_t RINV;            /* (2^16)^-1 mod q */
static uint16_t ZMONT[128];      /* psi^{brv7(k)} * R mod q, k=0..127 */
static uint16_t ZINV[128];       /* (ZMONT[k])^{-1} in Montgomery (for GS) */
static uint16_t GAMMA[128];      /* psi^{2*brv7(i)+1} mod q (basemul gamma) */
static uint16_t SCALETOMONT;     /* (2^7)^-1 * R^2 mod q */
uint16_t ntt_ref_R2;             /* R^2 mod q */

static uint32_t powmod(uint32_t b, uint32_t e){
    uint64_t r=1,x=b%NTT_Q; while(e){ if(e&1)r=r*x%NTT_Q; x=x*x%NTT_Q; e>>=1; } return (uint32_t)r;
}
/* SMALLEST primitive order-th root (same rule as gen_zetas.py / gen_ntt.py), so
 * this oracle matches the avx assembly's (ref-aligned) zeta layout. */
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
    uint32_t psi = prim_root(256);                          /* smallest primitive 256th root (=430) */
    for(int k=0;k<128;k++){
        uint32_t zplain = powmod(psi,brv7(k));
        ZMONT[k]=(uint16_t)((uint64_t)zplain*(R%NTT_Q)%NTT_Q);
        ZINV[k] =(uint16_t)((uint64_t)powmod(zplain,NTT_Q-2)*(R%NTT_Q)%NTT_Q); /* z^{-1} in Montgomery */
    }
    for(int i=0;i<128;i++) GAMMA[i]=(uint16_t)powmod(psi, 2*brv7(i)+1);
    ntt_ref_R2 = (uint16_t)((uint64_t)(R%NTT_Q)*(R%NTT_Q)%NTT_Q);
    SCALETOMONT = (uint16_t)((uint64_t)powmod(1u<<NTT_LEVELS,NTT_Q-2)*ntt_ref_R2%NTT_Q);
}

/* forward: 7 CT levels, len 128 -> 2 (stop before len=1).  Twiddle index k runs
 * 1.. in the (len desc, s asc) order; the inverse below reuses the SAME k per
 * (len,s) block so the two are exact inverses. */
void ntt_ref(uint16_t r[256]){
    int k=1;
    for(int len=128; len>=2; len>>=1)
        for(int s=0; s<256; s+=2*len){ uint16_t z=ZMONT[k++];
            for(int j=s;j<s+len;j++){ uint16_t t=ntt_ref_fqmul(z,r[j+len]); r[j+len]=subm(r[j],t); r[j]=addm(r[j],t); } }
}

/* basemul: 128 products mod (x^2 - gamma_i), pairs (r[2i], r[2i+1]) */
void basemul_ref(const uint16_t a[256], const uint16_t b[256], uint16_t c[256]){
    for(int i=0;i<128;i++){
        uint16_t a0=a[2*i], a1=a[2*i+1], b0=b[2*i], b1=b[2*i+1];
        /* (a0+a1 x)(b0+b1 x) mod (x^2-gamma): c0=a0b0+gamma a1b1, c1=a0b1+a1b0 */
        uint16_t a0b0=ntt_ref_fqmul(a0,b0);
        uint16_t a1b1=ntt_ref_fqmul(a1,b1);
        uint16_t a0b1=ntt_ref_fqmul(a0,b1);
        uint16_t a1b0=ntt_ref_fqmul(a1,b0);
        uint16_t gab =(uint16_t)((uint32_t)GAMMA[i]*a1b1%NTT_Q);   /* plain gamma*(a1*b1) */
        c[2*i]   = addm(a0b0, gab);
        c[2*i+1] = addm(a0b1, a1b0);
    }
}

/* inverse: reverse the 7 CT levels as GS (len 2 -> 128), reusing the SAME
 * twiddle index k = base(len) + s/(2*len) the forward used for block (len,s);
 * then scale by 1/2^7.  base(len) = 256/(2*len) is the count of all blocks at
 * the (strictly larger) lengths that came before this one in the forward order. */
void invntt_tomont_ref(uint16_t r[256]){
    for(int len=2; len<=128; len<<=1){
        int base = 256/(2*len);                       /* 64,32,16,8,4,2,1 */
        for(int s=0; s<256; s+=2*len){ uint16_t zi=ZINV[base + s/(2*len)];
            /* GS: a=A'+B'; b=z^{-1}(A'-B'); the 2^7 factor folds into the 1/2^7 scale */
            for(int j=s;j<s+len;j++){ uint16_t t=r[j]; r[j]=addm(t,r[j+len]); r[j+len]=ntt_ref_fqmul(zi,subm(t,r[j+len])); } }
    }
    for(int j=0;j<256;j++) r[j]=ntt_ref_fqmul(r[j],SCALETOMONT);
}

void polymul_schoolbook(const uint16_t a[256], const uint16_t b[256], uint16_t c[256]){
    for(int i=0;i<256;i++) c[i]=0;
    for(int i=0;i<256;i++) for(int j=0;j<256;j++){ uint16_t p=(uint16_t)((uint32_t)a[i]*b[j]%NTT_Q); int k=i+j;
        if(k<256) c[k]=addm(c[k],p); else c[k-256]=subm(c[k-256],p); }
}
