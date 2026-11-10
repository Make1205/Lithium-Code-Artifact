#include "copa_sampler_profile.h"
#ifdef COPA_IRS_CT_CHECK
#include <valgrind/memcheck.h>
#endif
#ifdef COPA_SAMPLER_EXP
#include "copa_sampler_legacy_experiment.inc"
#else
#include "copa_prof.h"
#include "copa_sampler.h"
#include "fips202.h"
#include "copa_xofbuf.h"
#include <string.h>
#include "copa_sampler_tables.h"
#include "copa_scalar64_tables.h"
#include "copa_gaussian_fp.h"
#include "ntt_crt.h"
#ifdef COPA_VEC_POLY_HEADER
#include COPA_VEC_POLY_HEADER
#endif
#ifdef COPA_PARALLEL_XOF
#include "copa_xof8.h"
#endif
#ifdef COPA_VEC_SAMPLER_HEADER
#include COPA_VEC_SAMPLER_HEADER
#endif
static void sample_batch_body(uint32_t candidate[32],uint32_t accepted[32],const uint64_t u[][32],const uint8_t low[32],const uint8_t tail[256]){
#ifdef COPA_VEC_SAMPLER_HEADER
 copa_backend_gaussian_batch(candidate,accepted,u,low,tail);
#else
 uint32_t x[32]={0};
 for(unsigned k=0;k<COPA_BASE_TAIL;++k)for(unsigned i=0;i<32;++i){
  uint64_t borrow=0;
  for(unsigned limb=0;limb<COPA_BASE_LIMBS32;++limb)
   borrow=((uint64_t)(uint32_t)u[limb][i]-copa_base_rcdt[k][limb]-borrow)>>63;
  x[i]+=(uint32_t)borrow;
 }
 for(unsigned j=0;j<32;++j){
  uint64_t num=(uint64_t)low[j]*(low[j]+(2u<<COPA_BASE_BITS)*x[j])*10000,rand=0;
  for(unsigned b=0;b<8;++b)rand|=(uint64_t)tail[8*j+b]<<(8*b);
  candidate[j]=(1u<<COPA_BASE_BITS)*x[j]+low[j];
  accepted[j]=((rand&UINT64_C(0x00ffffffffffffff))<copa_gaussian_threshold56(num))&((candidate[j]!=0)|(((unsigned)(rand>>56)&1u)^1u));
 }
#endif
}
void copa_sample_batch(uint32_t candidate[32],uint32_t accepted[32],const uint64_t u[][32],const uint8_t low[32],const uint8_t tail[256]){copa_gaussian_env env=copa_gaussian_enter();sample_batch_body(candidate,accepted,u,low,tail);copa_gaussian_leave(env);}
#define BASE_BYTES (32*4*COPA_BASE_LIMBS32)
#define Y_BYTES (4*COPA_BASE_BITS)
#define TAIL_BYTES 7
#define BATCH_BYTES (BASE_BYTES+Y_BYTES+32*TAIL_BYTES+4)
#define UNIT_BYTES (8*SHAKE256_RATE)
#ifdef COPA_SAMPLER_STATS
uint64_t copa_y_generated,copa_y_consumed,copa_y_candidates;
#define Y_STAT(x) do{x;}while(0)
#else
#define Y_STAT(x) ((void)0)
#endif
typedef struct {
#ifdef COPA_PARALLEL_XOF
 copa_xof8 engine;
#else
 keccak_state lane[8];
#endif
 uint8_t buf[BATCH_BYTES+UNIT_BYTES];size_t pos,avail;} stream;
static uint64_t le64(const uint8_t *b){uint64_t v=0;for(unsigned j=0;j<8;++j)v|=(uint64_t)b[j]<<(8*j);return v;}
static void ensure(stream *s,size_t n)
{
    while(s->avail<n){CP_BEGIN(pyshake);if(s->avail)memmove(s->buf,s->buf+s->pos,s->avail);s->pos=0;
#ifdef COPA_PARALLEL_XOF
 copa_xof8_block(&s->engine,(uint8_t (*)[SHAKE256_RATE])(s->buf+s->avail));
#else
 for(unsigned j=0;j<8;++j)shake256_squeezeblocks(s->buf+s->avail+j*SHAKE256_RATE,1,&s->lane[j]);
#endif
 CP_END(CP_Y_SHAKE,pyshake);Y_STAT(copa_y_generated+=UNIT_BYTES);s->avail+=UNIT_BYTES;}
}
void copa_expand_y(int32_t y[COPA_K][COPA_N],const uint8_t seed[32],uint64_t attempt)
{
    copa_gaussian_env fp_env=copa_gaussian_enter();
    stream st;st.pos=st.avail=0;uint8_t ext[42];memcpy(ext,seed,32);ext[32]=0x13;
    for(unsigned j=0;j<8;++j)ext[33+j]=(uint8_t)(attempt>>(8*j));
    CP_BEGIN(pabsorb);
#ifdef COPA_PARALLEL_XOF
    uint8_t inputs[8][42];const uint8_t *ptr[8];for(unsigned j=0;j<8;++j){memcpy(inputs[j],ext,41);inputs[j][41]=(uint8_t)j;ptr[j]=inputs[j];}copa_xof8_init(&st.engine,ptr,42);
#else
    for(unsigned j=0;j<8;++j){ext[41]=(uint8_t)j;shake256_absorb_once(&st.lane[j],ext,sizeof(ext));}
#endif
    CP_END(CP_Y_SHAKE,pabsorb);
    for(unsigned poly=0;poly<COPA_K;++poly){
        unsigned count=0;
        while(count<COPA_N){
            ensure(&st,BATCH_BYTES);uint8_t low[32];
            const uint8_t *base=st.buf+st.pos;
            uint64_t u[COPA_BASE_LIMBS32][32];
            for(unsigned limb=0;limb<COPA_BASE_LIMBS32;++limb)for(unsigned i=0;i<32;++i){
                const uint8_t *p=base+128*limb+4*i;
                u[limb][i]=(uint32_t)p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);
            }
            st.pos+=BASE_BYTES;st.avail-=BASE_BYTES;
            uint32_t acc=0;unsigned bits=0,idx=0;
            for(unsigned j=0;j<32;++j){while(bits<COPA_BASE_BITS){acc|=(uint32_t)st.buf[st.pos+idx++]<<bits;bits+=8;}low[j]=(uint8_t)(acc&((1u<<COPA_BASE_BITS)-1));acc>>=COPA_BASE_BITS;bits-=COPA_BASE_BITS;}
            st.pos+=Y_BYTES;st.avail-=Y_BYTES;
            uint32_t candidate[32],accepted[32];uint8_t tail[256];
            const uint8_t *packed=st.buf+st.pos;
            const uint8_t *signs=packed+32*TAIL_BYTES;
            for(unsigned j=0;j<32;++j){memcpy(tail+8*j,packed+7*j,7);tail[8*j+7]=(uint8_t)((signs[j/8]>>(j%8))&1);}
            sample_batch_body(candidate,accepted,u,low,tail);
            st.pos+=32*TAIL_BYTES+4;st.avail-=32*TAIL_BYTES+4;
            Y_STAT(copa_y_consumed+=BATCH_BYTES;copa_y_candidates+=32);
            for(unsigned j=0;j<32&&count<COPA_N;++j){unsigned sign=tail[8*j+7]&1;int32_t mask=-(int32_t)sign;
                y[poly][count]=((int32_t)candidate[j]^mask)-mask;if(accepted[j])++count;}


        }
    }
    copa_gaussian_leave(fp_env);
}


#endif
#ifdef COPA_IRS_EXP
#include "copa_irs_experiment.inc"
#else
static int irs_run(int32_t z[COPA_K][COPA_N],const uint8_t challenge[COPA_N],
             const copa_secret_key *sk,const uint32_t z_ntt[COPA_K][2][COPA_N],
             const uint8_t seed[32],uint64_t attempt,unsigned *steps
#ifdef COPA_USE_IRS_POOL
,copa_irs_pool *pool
#endif
)
{
    CP_BEGIN(pproj);uint32_t correlation[2][COPA_N];int64_t projection[COPA_N];
    copa_ntt_dot_reversed(correlation,z_ntt,sk->secret_ntt,COPA_K,COPA_N);copa_ntt_inverse_crt(projection,correlation,COPA_N);CP_END(CP_IRS_PROJECT,pproj);
    unsigned position[COPA_TAU],count=0;int32_t masks[COPA_TAU];
    for(unsigned j=0;j<COPA_N;++j)if(challenge[j])position[count++]=j;
    uint8_t ext[41];memcpy(ext,seed,32);ext[32]=0x15;for(unsigned j=0;j<8;++j)ext[33+j]=(uint8_t)(attempt>>(8*j));
    copa_xofbuf rng;
#ifdef COPA_USE_IRS_POOL
    if(pool){uint64_t group=attempt/COPA_IRS_POOL_LANES;if(!pool->ready||pool->group!=group){uint8_t inputs[COPA_IRS_POOL_LANES][41];const uint8_t *ptr[COPA_IRS_POOL_LANES];for(unsigned lane=0;lane<COPA_IRS_POOL_LANES;++lane){memcpy(inputs[lane],ext,33);uint64_t nonce=group*COPA_IRS_POOL_LANES+lane;for(unsigned j=0;j<8;++j)inputs[lane][33+j]=(uint8_t)(nonce>>(8*j));ptr[lane]=inputs[lane];}
#if COPA_IRS_POOL_LANES==4
     copa_xof4_init(&pool->engine,ptr,41);
#else
     copa_xof8_init(&pool->engine,ptr,41);
#endif
     pool->ready=1;pool->group=group;pool->filled=0;}}
    else
#endif
    copa_xof_init(&rng,ext,sizeof(ext));
    for(unsigned k=0;k<count;++k){
        unsigned j=position[k];uint8_t rb[48];uint64_t words[6];
        /* Step k is the fixed slice [48*k,48*(k+1)) of this attempt's
         * domain-15 stream; another attempt never shares its XOF state. */
#ifdef COPA_USE_IRS_POOL
        if(pool){while(pool->filled<48*(k+1)){uint8_t blocks[COPA_IRS_POOL_LANES][136];
#if COPA_IRS_POOL_LANES==4
         copa_xof4_block(&pool->engine,blocks);
#else
         copa_xof8_block(&pool->engine,blocks);
#endif
         for(unsigned lane=0;lane<COPA_IRS_POOL_LANES;++lane)memcpy(pool->data[lane]+pool->filled,blocks[lane],136);
         pool->filled+=136;}
         memcpy(rb,pool->data[attempt%COPA_IRS_POOL_LANES]+48*k,48);
        }else
#endif
        {copa_xof_read(&rng,rb,sizeof(rb));}
        for(unsigned b=0;b<6;++b)words[b]=le64(rb+8*b);
        CP_BEGIN(ptrans);unsigned minus,accepted=copa_irs64_step(&minus,projection[j],&sk->irs_context,words);CP_END(CP_IRS_TRANSITION,ptrans);*steps=k+1;
        if(!accepted)return 0; /* permitted step-boundary early rejection */
        masks[k]=-(int32_t)minus;int64_t mask=-(int64_t)minus;
        for(unsigned next=k+1;next<count;++next){unsigned target=position[next];int64_t c=sk->autocorrelation[target-j];projection[target]+=(c^mask)-mask;}
    }
    CP_BEGIN(prot);
    /* All updates commute. Rejected candidates need no vector writes. */
#ifdef COPA_VEC_POLY_HEADER
    for(unsigned r=0;r<COPA_K;++r)
        copa_backend_sparse_add(z[r],sk->s[r],position,masks,count);
#else
    for(unsigned r=0;r<COPA_K;++r)for(unsigned k=0;k<count;++k){unsigned j=position[k];int32_t mask=masks[k];
        for(unsigned i=j;i<COPA_N;++i)z[r][i]+=(sk->s[r][i-j]^mask)-mask;
        for(unsigned i=0;i<j;++i)z[r][i]-=(sk->s[r][i+COPA_N-j]^mask)-mask;
    }
#endif
    CP_END(CP_IRS_ROTATE,prot);return 1;
}

#endif

int copa_irs(int32_t z[COPA_K][COPA_N],const uint8_t c[COPA_N],const copa_secret_key *sk,const uint32_t zn[COPA_K][2][COPA_N],const uint8_t seed[32],uint64_t attempt,unsigned *steps){return irs_run(z,c,sk,zn,seed,attempt,steps
#ifdef COPA_USE_IRS_POOL
,NULL
#endif
);}
#ifdef COPA_USE_IRS_POOL
int copa_irs_with_pool(int32_t z[COPA_K][COPA_N],const uint8_t c[COPA_N],const copa_secret_key *sk,const uint32_t zn[COPA_K][2][COPA_N],const uint8_t seed[32],uint64_t attempt,unsigned *steps,copa_irs_pool *pool){return irs_run(z,c,sk,zn,seed,attempt,steps,pool);}
#endif
