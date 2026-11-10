#include "copa_prof.h"
#include "copa_api.h"
#include "copa_sampler.h"
#include "copa_rans.h"
#include "quant.h"
#include "pack_pk.h"
#include "ntt_crt.h"
#include "fips202.h"
#include "copa_xofbuf.h"
#include "randombytes.h"
#ifdef COPA_VEC_POLY_HEADER
#include COPA_VEC_POLY_HEADER
#endif
#include <string.h>
#ifdef COPA_BENCH_ITERATIONS
#include "../common/copa_bench_ticks.h"
uint64_t copa_bench_iteration_cycles,copa_bench_attempts;
#endif
static void message_hash(uint8_t mu[64],const uint8_t hpk[64],const uint8_t *m,size_t len)
{
    keccak_state st;uint8_t dom=3;shake256_init(&st);shake256_absorb(&st,hpk,64);shake256_absorb(&st,m,len);shake256_absorb(&st,&dom,1);shake256_finalize(&st);shake256_squeeze(mu,64,&st);
}
static void challenge_hash(uint8_t seed[32],const int32_t w1[COPA_M][COPA_N],const uint8_t mu[64])
{
    keccak_state st;uint8_t dom=5;shake256_init(&st);
    uint8_t packed[2*COPA_M*COPA_N];for(unsigned r=0;r<COPA_M;++r)for(unsigned i=0;i<COPA_N;++i){packed[2*(r*COPA_N+i)]=(uint8_t)w1[r][i];packed[2*(r*COPA_N+i)+1]=(uint8_t)(w1[r][i]>>8);}
    shake256_absorb(&st,packed,sizeof(packed));
    shake256_absorb(&st,mu,64);shake256_absorb(&st,&dom,1);shake256_finalize(&st);shake256_squeeze(seed,32,&st);
}
static void ball(uint8_t c[COPA_N],const uint8_t seed[32])
{
    uint8_t in[33],rb;copa_xofbuf st;memcpy(in,seed,32);in[32]=6;copa_xof_init(&st,in,sizeof(in));memset(c,0,COPA_N);
    for(unsigned i=COPA_N-COPA_TAU;i<COPA_N;++i){do{rb=copa_xof_u8(&st);}while(rb>i);c[i]=c[rb];c[rb]=1;}
}
#include "copa_round.h"
#define highbits copa_highbits
static int32_t wrap_small(int32_t x,unsigned q){x+=(int32_t)q*(x<0);return x-(int32_t)q*((uint32_t)x>=q);}
static int32_t centered(int32_t x,unsigned q){return x-(int32_t)q*(int32_t)((unsigned)x>q/2);}
int copa_sign(uint8_t sig[COPA_SIG_BYTES],const uint8_t *m,size_t len,const copa_secret_key *sk,copa_sign_stats *stats)
{
    if(stats)memset(stats,0,sizeof(*stats));
    CP_BEGIN(pset);uint64_t norm=0;for(unsigned r=0;r<COPA_K;++r)for(unsigned i=0;i<COPA_N;++i)norm+=(uint64_t)((int64_t)sk->s[r][i]*sk->s[r][i]);
    if(norm<COPA_MIN_NORM||norm>COPA_BK_SQ)return -1;

    uint8_t mu[64],seed[32],in[129];message_hash(mu,sk->hpk,m,len);
    memcpy(in,sk->key,32);randombytes(in+32,32);memcpy(in+64,mu,64);in[128]=4;shake256(seed,32,in,sizeof(in));
    if(norm!=sk->irs_context.norm)return -1;
    CP_END(CP_SIGN_SETUP,pset);
#ifdef COPA_USE_IRS_POOL
    copa_irs_pool irs_pool;irs_pool.ready=0;
#endif
#ifdef COPA_BENCH_ITERATIONS
    uint64_t iteration_start=copa_bench_ticks();
#endif
    for(uint64_t attempt=0;;++attempt){
#ifdef COPA_BENCH_ITERATIONS
        ++copa_bench_attempts;
#endif
        if(stats)++stats->attempts;
        int32_t z[COPA_K][COPA_N],w[COPA_M][COPA_N],w1[COPA_M][COPA_N],h[COPA_M][COPA_N];
        uint8_t c[COPA_N];uint32_t z_ntt[COPA_K][2][COPA_N];CP_BEGIN(py);copa_expand_y(z,seed,attempt);CP_END(CP_Y,py);CP_BEGIN(pay);
        for(unsigned r=0;r<COPA_K;++r)copa_ntt_forward(z_ntt[r],z[r],COPA_N);
        copa_matrix_apply_ntt(w,sk->public_ntt,z_ntt);CP_END(CP_AY,pay);CP_BEGIN(ph);
#ifdef COPA_VEC_POLY_HEADER
        copa_backend_commit(&w[0][0],&w1[0][0],&z[1+COPA_ELL][0]);
#else
        for(unsigned r=0;r<COPA_M;++r)for(unsigned i=0;i<COPA_N;++i){w[r][i]=wrap_small(w[r][i]+z[1+COPA_ELL+r][i],COPA_Q);w1[r][i]=highbits(w[r][i]);}
#endif
        challenge_hash(sig,w1,mu);ball(c,sig);CP_END(CP_HIGH_CH,ph);CP_BEGIN(pi);
        unsigned steps;int irs_ok;
#ifdef COPA_USE_IRS_POOL
        irs_ok=copa_irs_with_pool(z,c,sk,z_ntt,seed,attempt,&steps,&irs_pool);
#else
        irs_ok=copa_irs(z,c,sk,z_ntt,seed,attempt,&steps);
#endif
        if(stats)stats->irs_steps+=steps;
        CP_END(CP_IRS,pi);if(!irs_ok)continue;
        if(stats)++stats->irs_accepted;
        CP_BEGIN(pbs);int64_t bs=0;
#ifdef COPA_VEC_POLY_HEADER
        bs=(int64_t)copa_backend_norm2(&z[0][0],COPA_K*COPA_N);
#else
        for(unsigned r=0;r<COPA_K;++r)for(unsigned i=0;i<COPA_N;++i)bs+=(int64_t)z[r][i]*z[r][i];
#endif
        CP_END(CP_BS,pbs);if(bs>(int64_t)COPA_BS*COPA_BS){if(stats)++stats->bs_restarts;continue;}
        CP_BEGIN(pbv);int64_t bv=0;
#ifdef COPA_VEC_POLY_HEADER
        bv=(int64_t)copa_backend_norm2(&z[0][0],(1+COPA_ELL)*COPA_N);
        bv+=(int64_t)copa_backend_make_hint(&h[0][0],&w[0][0],&w1[0][0],&z[1+COPA_ELL][0]);
#else
for(unsigned r=0;r<1+COPA_ELL;++r)for(unsigned i=0;i<COPA_N;++i)bv+=(int64_t)z[r][i]*z[r][i];
        for(unsigned r=0;r<COPA_M;++r)for(unsigned i=0;i<COPA_N;++i){
            int32_t u=wrap_small(w[r][i]-z[1+COPA_ELL+r][i],COPA_Q);
            int32_t diff=wrap_small(w1[r][i]-highbits(u),COPA_Q/COPA_ALPHA);
            h[r][i]=centered(diff,COPA_Q/COPA_ALPHA);
            int32_t zt=centered(wrap_small(COPA_ALPHA*w1[r][i]-u,COPA_Q),COPA_Q);bv+=(int64_t)zt*zt;
        }
#endif
        CP_END(CP_HINT_BV,pbv);if(bv>(int64_t)COPA_BV*COPA_BV){if(stats)++stats->bv_restarts;continue;}
        CP_BEGIN(prans);size_t used;int enc_ok=copa_encode(sig+32,&used,(const int32_t (*)[COPA_N])z,h);CP_END(CP_RANS,prans);if(!enc_ok){if(stats)++stats->encode_overflows;continue;}
        if(stats)stats->payload_used=used;
#ifdef COPA_BENCH_ITERATIONS
        copa_bench_iteration_cycles+=copa_bench_ticks()-iteration_start;
#endif
        return 0;
    }
}
int copa_verify_prepared(const uint8_t *sig,size_t siglen,const uint8_t *m,size_t len,const copa_public_key *key)
{
    if(siglen!=COPA_SIG_BYTES)return -1;
    int32_t z1[1+COPA_ELL][COPA_N],h[COPA_M][COPA_N],u[COPA_M][COPA_N],w1[COPA_M][COPA_N];
    CP_BEGIN(pd);int decoded=copa_decode(z1,h,sig+32);CP_END(CP_DECODE,pd);if(!decoded)return -1;
    CP_BEGIN(pnorm);uint64_t z1_norm=0;
#ifdef COPA_VEC_POLY_HEADER
    z1_norm=copa_backend_norm2(&z1[0][0],(1+COPA_ELL)*COPA_N);
#else
    for(unsigned r=0;r<1+COPA_ELL;++r)for(unsigned i=0;i<COPA_N;++i)z1_norm+=(uint64_t)((int64_t)z1[r][i]*z1[r][i]);
#endif
    CP_END(CP_VF_NORM,pnorm);if(z1_norm>(uint64_t)COPA_BV*COPA_BV)return -1;
    uint8_t mu[64],c[COPA_N],seed[32];
    CP_BEGIN(pc);message_hash(mu,key->hpk,m,len);ball(c,sig);CP_END(CP_CHALLENGE,pc);CP_BEGIN(pntt);
#ifdef COPA_VERIFY31
    copa_verify31_apply(u,key->matrix,z1);
#else
    copa_matrix_apply(u,key->matrix,z1);
#endif
    CP_END(CP_VF_NTT,pntt);CP_BEGIN(prec);int64_t bv=(int64_t)z1_norm;
#ifdef COPA_VEC_POLY_HEADER
    bv+=copa_backend_reconstruct(&w1[0][0],&u[0][0],&h[0][0],c);
#else
    for(unsigned r=0;r<COPA_M;++r)for(unsigned i=0;i<COPA_N;++i){
        int32_t value=wrap_small(u[r][i]-(r==0?COPA_Q/2*c[i]:0),COPA_Q);
        w1[r][i]=wrap_small(highbits(value)+h[r][i],COPA_Q/COPA_ALPHA);
        int32_t zt=centered(wrap_small(COPA_ALPHA*w1[r][i]-value,COPA_Q),COPA_Q);bv+=(int64_t)zt*zt;
    }
#endif
    CP_END(CP_RECONSTRUCT,prec);if(bv>(int64_t)COPA_BV*COPA_BV)return -1;
    CP_BEGIN(pc2);challenge_hash(seed,w1,mu);CP_END(CP_CHALLENGE,pc2);unsigned bad=0;for(unsigned i=0;i<32;++i)bad|=seed[i]^sig[i];return bad?-1:0;
}

int copa_verify(const uint8_t *sig,size_t siglen,const uint8_t *m,size_t len,const uint8_t *pk,size_t pklen){
 if(siglen!=COPA_SIG_BYTES||pklen!=COPA_PK_BYTES)return -1;
 copa_public_key key;if(copa_prepare_public(&key,pk,pklen))return -1;
 return copa_verify_prepared(sig,siglen,m,len,&key);
}
