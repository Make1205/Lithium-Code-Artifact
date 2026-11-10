#include "copa_prof.h"
#include "copa_keygen.h"
#include "quant.h"
#include "pack_pk.h"
#include "poly_inv.h"
#include "fips202.h"
#include "copa_xofbuf.h"
#include "randombytes.h"
#include <string.h>
#include "copa_noise_tables.h"
#ifdef COPA_VEC_POLY_HEADER
#include COPA_VEC_POLY_HEADER
#endif
#ifdef COPA_VEC_LATTICE
#include "../common/copa_cdt_vectors.h"
#if COPA_SIMD_BITS==256
#include "../common/copa_compact8.h"
#endif
#endif
#ifdef COPA_PARALLEL_XOF
#include "copa_xof8.h"
#endif

void copa_expand_a(int32_t a[COPA_M][COPA_ELL][COPA_N],const uint8_t rho[32])
{
    int32_t *flat=&a[0][0][0];unsigned per=COPA_M*COPA_ELL*COPA_N/8;
#if defined(COPA_PARALLEL_XOF) && defined(COPA_VEC_LATTICE) && COPA_M*COPA_ELL>=12
    {
#define A_PER (COPA_M*COPA_ELL*COPA_N/8)
#define A_THRESHOLD (65536/COPA_Q*COPA_Q)
#ifndef COPA_A_TEST_PREFIX_BLOCKS
#define A_BLOCKS ((2*A_PER*65536+A_THRESHOLD*136-1)/(A_THRESHOLD*136))
#else
#define A_BLOCKS COPA_A_TEST_PREFIX_BLOCKS
#endif
     uint8_t raw[8][A_BLOCKS*136];uint8_t *ptr[8];copa_xof8 batch;unsigned failed=0;
     for(unsigned lane=0;lane<8;++lane)ptr[lane]=raw[lane];
     copa_xof8_seed_init(&batch,rho,0x11);copa_xof8_prefix(&batch,ptr,A_BLOCKS*136);
     for(unsigned lane=0;lane<8;++lane){unsigned count=0,j=0;
      for(;j+2*COPAV_W<=A_BLOCKS*136&&count+COPAV_W<=per;j+=2*COPAV_W){
#if COPA_SIMD_BITS==256
       copav_vec v=_mm256_cvtepu16_epi32(_mm_loadu_si128((const __m128i *)(raw[lane]+j)));
#else
       copav_vec v=_mm512_cvtepu16_epi32(_mm256_loadu_si256((const __m256i *)(raw[lane]+j)));
#endif
       /* v<2^16: floor(2^20/q) underestimates the quotient by <=1. */
       copav_vec quot=cv_srl32(cv_mul32(v,cv_set((1u<<20)/COPA_Q)),20),rem=cv_reduce(cv_sub(v,cv_mul32(quot,cv_set(COPA_Q))),COPA_Q);
#if COPA_SIMD_BITS==256
       unsigned mask=(unsigned)_mm256_movemask_ps(_mm256_castsi256_ps(cv_gt(cv_set(A_THRESHOLD),v)));
       cv_store(flat+lane*per+count,cv_perm(rem,cv_load(copa_compact8[mask])));
#else
       unsigned mask=_mm512_cmpgt_epi32_mask(cv_set(A_THRESHOLD),v);_mm512_mask_compressstoreu_epi32(flat+lane*per+count,(__mmask16)mask,rem);
#endif
       count+=(unsigned)__builtin_popcount(mask);
      }
      for(;j+2<=A_BLOCKS*136&&count<per;j+=2){unsigned v=raw[lane][j]+256u*raw[lane][j+1];if(v<A_THRESHOLD)flat[lane*per+count++]=(int32_t)(v%COPA_Q);}
      failed|=count<per;
     }
     if(!failed)return;
     /* Public-seed rare fallback: restart this local XOF, overwrite all A. */
#undef A_PER
#undef A_THRESHOLD
#undef A_BLOCKS
    }
#endif

#ifdef COPA_PARALLEL_XOF
    copa_xof8 engine;uint8_t blocks[8][136];unsigned count[8]={0},done=0;
    copa_xof8_seed_init(&engine,rho,0x11);
    while(done<8*per){copa_xof8_block(&engine,blocks);for(unsigned lane=0;lane<8;++lane){unsigned j=0;
#ifdef COPA_VEC_LATTICE
     for(;j+2*COPAV_W<=136&&count[lane]+COPAV_W<=per;j+=2*COPAV_W){
#if COPA_SIMD_BITS==256
      copav_vec v=_mm256_cvtepu16_epi32(_mm_loadu_si128((const __m128i *)(blocks[lane]+j)));
#else
      copav_vec v=_mm512_cvtepu16_epi32(_mm256_loadu_si256((const __m256i *)(blocks[lane]+j)));
#endif
      copav_vec reduced=cv_sub(v,cv_mul32(qv_div_positive(v,COPA_Q),cv_set(COPA_Q)));
#if COPA_SIMD_BITS==256
      unsigned mask=(unsigned)_mm256_movemask_ps(_mm256_castsi256_ps(cv_gt(cv_set(65536/COPA_Q*COPA_Q),v)));
      cv_store(flat+lane*per+count[lane],cv_perm(reduced,cv_load(copa_compact8[mask])));
#else
      unsigned mask=_mm512_cmpgt_epi32_mask(cv_set(65536/COPA_Q*COPA_Q),v);
      _mm512_mask_compressstoreu_epi32(flat+lane*per+count[lane],(__mmask16)mask,reduced);
#endif
      unsigned accepted=(unsigned)__builtin_popcount(mask);count[lane]+=accepted;done+=accepted;
     }
#endif
     for(;j<136&&count[lane]<per;j+=2){unsigned v=blocks[lane][j]+256u*blocks[lane][j+1];if(v<65536/COPA_Q*COPA_Q){flat[lane*per+count[lane]++]=(int32_t)(v%COPA_Q);++done;}}
    }}
#else
    for(unsigned lane=0;lane<8;++lane){
        uint8_t ext[35];copa_xofbuf st;memcpy(ext,rho,32);ext[32]=0x11;ext[33]=(uint8_t)lane;ext[34]=0;
        copa_xof_init(&st,ext,sizeof(ext));
        for(unsigned j=0;j<per;++j){unsigned v;do{v=copa_xof_u16(&st);}while(v>=65536/COPA_Q*COPA_Q);flat[lane*per+j]=(int32_t)(v%COPA_Q);}
    }
#endif
}
void copa_expand_s(int32_t s[1+COPA_ELL][COPA_N],const uint8_t seed[32])
{
    const uint64_t (*cdt)[3]=COPA_SIGMA100==150?noise_s150:noise_s070;
    unsigned tail=COPA_SIGMA100==150?18:9,per=(1+COPA_ELL)*COPA_N/8;
    int32_t *flat=&s[0][0];
#ifdef COPA_PARALLEL_XOF
    copa_xof8 engine;uint8_t raw[8][((1+COPA_ELL)*COPA_N/8)*24];uint8_t *out[8];
    for(unsigned lane=0;lane<8;++lane)out[lane]=raw[lane];
    copa_xof8_seed_init(&engine,seed,0x12);copa_xof8_prefix(&engine,out,per*24);
#endif
    for(unsigned lane=0;lane<8;++lane){
#ifndef COPA_PARALLEL_XOF
        uint8_t ext[35];copa_xofbuf st;memcpy(ext,seed,32);ext[32]=0x12;ext[33]=(uint8_t)lane;ext[34]=0;
        copa_xof_init(&st,ext,sizeof(ext));
#endif
        for(unsigned group=0;group<per;group+=32){
            uint8_t rb[32*24];uint64_t u[3][32];int32_t value[32];
#ifdef COPA_PARALLEL_XOF
            memcpy(rb,raw[lane]+24*group,sizeof(rb));
#else
            copa_xof_read(&st,rb,sizeof(rb));
#endif
            for(unsigned j=0;j<32;++j){value[j]=-(int32_t)tail;for(unsigned limb=0;limb<3;++limb){uint64_t v=0;for(unsigned b=0;b<8;++b)v|=(uint64_t)rb[24*j+8*limb+b]<<(8*b);u[limb][j]=v;}}
#ifdef COPA_VEC_LATTICE
            uint32_t counts[32];copa_cdt_vectors(counts,u,cdt,2*tail);for(unsigned j=0;j<32;++j)value[j]=(int32_t)tail-(int32_t)counts[j];
#else
            for(unsigned k=0;k<2*tail;++k){uint64_t t0=cdt[k][0],t1=cdt[k][1],t2=cdt[k][2];
                for(unsigned j=0;j<32;++j){unsigned less=(u[2][j]<t2)|((u[2][j]==t2)&((u[1][j]<t1)|((u[1][j]==t1)&(u[0][j]<t0))));value[j]+=(int32_t)(1-less);}
            }
#endif
            for(unsigned j=0;j<32;++j)flat[lane*per+group+j]=value[j];
        }
    }
}
static int copa_keygen_internal(uint8_t pk[COPA_PK_BYTES],copa_secret_key *sk,copa_keygen_stats *stats,unsigned prepare)
{
    if(stats)memset(stats,0,sizeof(*stats));
    for(;;){
        uint8_t zeta[32],seeds[96],in[33];int32_t finv[COPA_N];
        int32_t a[COPA_M][COPA_ELL][COPA_N],d[COPA_M][COPA_N],x[COPA_M][COPA_N],b[COPA_M][COPA_N];
        if(stats)++stats->attempts;
        randombytes(zeta,32);CP_BEGIN(phash);memcpy(in,zeta,32);in[32]=0x01;shake256(seeds,96,in,sizeof(in));CP_END(CP_KG_HASH,phash);
        CP_BEGIN(ps);copa_expand_s((int32_t (*)[COPA_N])sk->s,seeds+32);
        CP_END(CP_S,ps);
        uint32_t parity=0;for(unsigned i=0;i<COPA_N;++i)parity^=(uint32_t)sk->s[0][i];
#ifndef COPA_KEYGEN_FULL_INVERSE
        /* For even q, an even f(1) shares X+1 with X^256+1 modulo 2.
         * This is the existing nonunit rejection, before its expensive
         * norm-tower inversion. Discarded attempts consume identical RNG. */
        if(!(parity&1u)){if(stats){++stats->noninvertible;++stats->noninvertible_even;}continue;}
#endif
        CP_BEGIN(pearly);uint64_t partial_norm;
#ifdef COPA_VEC_POLY_HEADER
        partial_norm=copa_backend_norm2(&sk->s[0][0],(1+COPA_ELL)*COPA_N);
#else
        partial_norm=0;for(unsigned r=0;r<1+COPA_ELL;++r)for(unsigned i=0;i<COPA_N;++i)partial_norm+=(uint64_t)((int64_t)sk->s[r][i]*sk->s[r][i]);
#endif
        CP_END(CP_KG_NORM,pearly);
#if COPA_PARTIAL_MIN_NORM>0
        if(partial_norm<COPA_PARTIAL_MIN_NORM){if(stats){++stats->lower_norm_reject;++stats->partial_lower_norm_reject;}continue;}
#endif
        if(partial_norm>COPA_BK_SQ){if(stats){++stats->norm_reject;++stats->partial_norm_reject;}continue;}
        CP_BEGIN(pa);copa_expand_a(a,seeds);CP_END(CP_A,pa);CP_BEGIN(pd);copa_expand_d(&d[0][0],seeds,COPA_MODE);CP_END(CP_D,pd);
        CP_BEGIN(pkgmat);uint32_t fg[1+COPA_ELL][2][COPA_N];
        for(unsigned j=0;j<1+COPA_ELL;++j)copa_ntt_forward(fg[j],sk->s[j],COPA_N);
        for(unsigned r=0;r<COPA_M;++r){int32_t first[COPA_N];
            for(unsigned i=0;i<COPA_N;++i)first[i]=-d[r][i];
            copa_ntt_forward(sk->public_ntt[r][0],first,COPA_N);
            for(unsigned j=0;j<COPA_ELL;++j)copa_ntt_forward(sk->public_ntt[r][1+j],a[r][j],COPA_N);
        }
        copa_matrix_apply_ntt(x,sk->public_ntt,fg);CP_END(CP_KG_MATRIX,pkgmat);
        CP_BEGIN(pq);
        int32_t cvec[COPA_M][COPA_N];
#ifdef COPA_VEC_LATTICE
        copa_quant_all(&cvec[0][0],&sk->s[1+COPA_ELL][0],&x[0][0]);
#else
        unsigned width=copa_quant_width(COPA_MODE);
        for(unsigned k=0;k<64;++k){
            int32_t block[16],quant[16];
            for(unsigned r=0;r<width;++r)block[r]=x[r/4][k+(r%4)*64];
            copa_quant_block(quant,block,COPA_MODE);
            for(unsigned r=0;r<width;++r){unsigned poly=r/4,idx=k+(r%4)*64;cvec[poly][idx]=quant[r];sk->s[1+COPA_ELL+poly][idx]=quant[r]-block[r];}
        }
#endif
        CP_END(CP_QUANT,pq);CP_BEGIN(pnorm);int64_t norm=0;for(unsigned r=0;r<COPA_K;++r)for(unsigned i=0;i<COPA_N;++i)norm+=(int64_t)sk->s[r][i]*sk->s[r][i];
        CP_END(CP_KG_NORM,pnorm);if(norm<COPA_MIN_NORM){if(stats)++stats->lower_norm_reject;continue;}
        if(norm>COPA_BK_SQ){if(stats)++stats->norm_reject;continue;}
        CP_BEGIN(pinv);int unit=copa_poly_inv(finv,sk->s[0],COPA_Q);CP_END(CP_INVERSE,pinv);if(!unit){if(stats){++stats->noninvertible;if(parity&1u)++stats->noninvertible_odd;else ++stats->noninvertible_even;}continue;}
        CP_BEGIN(pcompb);cvec[0][0]+=COPA_Q/2;
        uint32_t inv_ntt[1][2][COPA_N];copa_ntt_forward(inv_ntt[0],finv,COPA_N);
        for(unsigned r=0;r<COPA_M;++r){
            uint32_t cv[1][2][COPA_N],prod[2][COPA_N];int64_t product[COPA_N];
            copa_ntt_forward(cv[0],cvec[r],COPA_N);copa_ntt_dot(prod,cv,inv_ntt,1,COPA_N);copa_ntt_inverse_crt(product,prod,COPA_N);
            int32_t first[COPA_N];
            for(unsigned i=0;i<COPA_N;++i){b[r][i]=(int32_t)copa_mod_i64(product[i],COPA_Q);first[i]=-(int32_t)copa_mod_i64((int64_t)b[r][i]+d[r][i],COPA_Q);}
            if(prepare)copa_ntt_forward(sk->public_ntt[r][0],first,COPA_N);
        }
        CP_END(CP_KG_B,pcompb);CP_BEGIN(ppack);if(!copa_pack_pk(pk,COPA_PK_BYTES,seeds,&b[0][0],COPA_MODE)){if(stats)++stats->packing_failure;return -1;}
        memcpy(sk->pk,pk,COPA_PK_BYTES);memcpy(sk->key,seeds+64,32);
        keccak_state st;uint8_t domain=0x02;shake256_init(&st);shake256_absorb(&st,pk,COPA_PK_BYTES);shake256_absorb(&st,&domain,1);shake256_finalize(&st);shake256_squeeze(sk->hpk,64,&st);
        CP_END(CP_PACK,ppack);if(prepare){
        memcpy(sk->secret_ntt,fg,sizeof(fg));
        for(unsigned r=1+COPA_ELL;r<COPA_K;++r)copa_ntt_forward(sk->secret_ntt[r],sk->s[r],COPA_N);
        uint32_t correlation[2][COPA_N];copa_ntt_dot_reversed(correlation,sk->secret_ntt,sk->secret_ntt,COPA_K,COPA_N);
        copa_ntt_inverse_crt(sk->autocorrelation,correlation,COPA_N);
        copa_irs64_prepare(&sk->irs_context,(uint64_t)norm);
        }
        return 0;
    }
}

int copa_keygen(uint8_t pk[COPA_PK_BYTES],copa_secret_key *sk,copa_keygen_stats *stats){return copa_keygen_internal(pk,sk,stats,1);}
int copa_keygen_raw(uint8_t pk[COPA_PK_BYTES],copa_secret_key *sk,copa_keygen_stats *stats){return copa_keygen_internal(pk,sk,stats,0);}
