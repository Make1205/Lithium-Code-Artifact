/* Historical opt-in k experiments; scalar proposal oracle on every backend. */
#ifdef COPA_VEC_SAMPLER_HEADER
#undef COPA_VEC_SAMPLER_HEADER
#endif
#include "copa_prof.h"
#ifdef COPA_IRS_CT_CHECK
#include <valgrind/memcheck.h>
#endif
#include "copa_sampler_profile.h"
#include "copa_sampler.h"
#include "fips202.h"
#include "copa_xofbuf.h"
#include <string.h>
#include "copa_sampler_tables.h"
#ifdef COPA_SAMPLER_EXP
#include "copa_sampler_exp.h"
#endif
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
 SP_BEGIN(pb);uint32_t x[32]={0};
#ifdef COPA_SAMPLER_EXP
 for(unsigned k=0;k<COPA_BASE_TAIL;++k)for(unsigned i=0;i<32;++i){unsigned lt=0;for(unsigned l=0;l<COPA_BASE_LIMBS;++l)lt=(u[l][i]<copa_base_rcdt[k][l])|((u[l][i]==copa_base_rcdt[k][l])&lt);x[i]+=lt;}
#else
 for(unsigned k=0;k<COPA_BASE_TAIL;++k){uint64_t t0=copa_base_rcdt[k][0],t1=copa_base_rcdt[k][1],t2=copa_base_rcdt[k][2];for(unsigned i=0;i<32;++i)x[i]+=(u[2][i]<t2)|((u[2][i]==t2)&((u[1][i]<t1)|((u[1][i]==t1)&(u[0][i]<t0))));}
#endif
 SP_END(SP_BASE,pb);SP_BEGIN(pe);
 for(unsigned j=0;j<32;++j){uint64_t num=(uint64_t)low[j]*(low[j]+(2u<<COPA_BASE_BITS)*x[j])*10000;
 candidate[j]=(1u<<COPA_BASE_BITS)*x[j]+low[j];
#ifdef COPA_SAMPLER_EXP
 uint64_t words[3];for(unsigned w=0;w<3;++w){words[w]=0;for(unsigned b=0;b<8;++b)words[w]|=(uint64_t)tail[24*j+8*w+b]<<(8*b);}
 accepted[j]=copa_exp_accept(copa_exp_probability_bits(num),words)&((candidate[j]!=0)|(unsigned)(words[2]&1));
#else
 uint64_t rand=0;for(unsigned b=0;b<8;++b)rand|=(uint64_t)tail[8*j+b]<<(8*b);accepted[j]=((rand>>1)<copa_gaussian_threshold(num))&((candidate[j]!=0)|(unsigned)(rand&1));
#endif
 }SP_END(SP_EXP,pe);
#endif
}
void copa_sample_batch(uint32_t candidate[32],uint32_t accepted[32],const uint64_t u[][32],const uint8_t low[32],const uint8_t tail[256]){copa_gaussian_env env=copa_gaussian_enter();sample_batch_body(candidate,accepted,u,low,tail);copa_gaussian_leave(env);}
#define BASE_BYTES (32*8*COPA_BASE_LIMBS)
#define Y_BYTES (4*COPA_BASE_BITS)
#ifdef COPA_SAMPLER_EXP
#define TAIL_BYTES 24
#else
#define TAIL_BYTES 8
#endif
#define BATCH_BYTES (BASE_BYTES+Y_BYTES+32*TAIL_BYTES)
#define UNIT_BYTES (8*SHAKE256_RATE)
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
    while(s->avail<n){if(s->avail)memmove(s->buf,s->buf+s->pos,s->avail);s->pos=0;
SP_BEGIN(px);
#ifdef COPA_PARALLEL_XOF
 copa_xof8_block(&s->engine,(uint8_t (*)[SHAKE256_RATE])(s->buf+s->avail));
#else
 for(unsigned j=0;j<8;++j)shake256_squeezeblocks(s->buf+s->avail+j*SHAKE256_RATE,1,&s->lane[j]);
#endif
 SP_END(SP_Y_SHAKE,px);SC_ADD(SC_Y_SQUEEZED,UNIT_BYTES);s->avail+=UNIT_BYTES;}
}
void copa_expand_y(int32_t y[COPA_K][COPA_N],const uint8_t seed[32],uint64_t attempt)
{
    copa_gaussian_env fp_env=copa_gaussian_enter();
    SC_ADD(SC_Y_CALLS,1);SP_BEGIN(pinit);stream st;st.pos=st.avail=0;uint8_t ext[42];memcpy(ext,seed,32);ext[32]=0x13;
    for(unsigned j=0;j<8;++j)ext[33+j]=(uint8_t)(attempt>>(8*j));
#ifdef COPA_PARALLEL_XOF
    uint8_t inputs[8][42];const uint8_t *ptr[8];for(unsigned j=0;j<8;++j){memcpy(inputs[j],ext,41);inputs[j][41]=(uint8_t)j;ptr[j]=inputs[j];}copa_xof8_init(&st.engine,ptr,42);
#else
    for(unsigned j=0;j<8;++j){ext[41]=(uint8_t)j;shake256_absorb_once(&st.lane[j],ext,sizeof(ext));}
#endif
    SP_END(SP_Y_SHAKE,pinit);
#ifdef COPA_SAMPLER_EXP
    unsigned complete=0;
    do {
     uint8_t signs[COPA_K*32];ensure(&st,sizeof(signs));memcpy(signs,st.buf+st.pos,sizeof(signs));st.pos+=sizeof(signs);st.avail-=sizeof(signs);SC_ADD(SC_Y_USED,sizeof(signs));
     unsigned count=0;
     for(unsigned block=0;block<COPA_EXP_CANDIDATES/32;++block){
      ensure(&st,BATCH_BYTES);const uint8_t *base=st.buf+st.pos;uint64_t u[COPA_BASE_LIMBS][32];uint8_t low[32];
      for(unsigned i=0;i<32;++i){unsigned off=(i/8)*(64*COPA_BASE_LIMBS)+(i%8)*8;for(unsigned limb=0;limb<COPA_BASE_LIMBS;++limb)u[limb][i]=le64(base+off+64*limb);}
      st.pos+=BASE_BYTES;st.avail-=BASE_BYTES;uint32_t acc=0;unsigned bits=0,idx=0;
      for(unsigned j=0;j<32;++j){while(bits<COPA_BASE_BITS){acc|=(uint32_t)st.buf[st.pos+idx++]<<bits;bits+=8;}low[j]=(uint8_t)(acc&((1u<<COPA_BASE_BITS)-1));acc>>=COPA_BASE_BITS;bits-=COPA_BASE_BITS;}
      st.pos+=Y_BYTES;st.avail-=Y_BYTES;uint32_t candidate[32],accepted[32];sample_batch_body(candidate,accepted,u,low,st.buf+st.pos);st.pos+=32*TAIL_BYTES;st.avail-=32*TAIL_BYTES;
      SC_ADD(SC_CANDIDATES,32);SC_ADD(SC_Y_USED,BATCH_BYTES);SP_BEGIN(psel);
      for(unsigned j=0;j<32;++j){SC_ADD(SC_ACCEPTED,accepted[j]);SC_ADD(SC_EXP_ACCEPTED,accepted[j]+((candidate[j]==0)&!accepted[j]));if(count<COPA_K*COPA_N){unsigned sign=(signs[count/8]>>(count%8))&1;int32_t mask=-(int32_t)sign;(&y[0][0])[count]=((int32_t)candidate[j]^mask)-mask;if(accepted[j])++count;}}
      SP_END(SP_SELECT,psel);
     }
     complete=count==COPA_K*COPA_N;SC_ADD(SC_WHOLE_RETRIES,!complete);
    }while(!complete);
#else
    for(unsigned poly=0;poly<COPA_K;++poly){
        uint8_t signs[32];ensure(&st,32);memcpy(signs,st.buf+st.pos,32);st.pos+=32;st.avail-=32;SC_ADD(SC_Y_USED,32);
        unsigned count=0;
        while(count<COPA_N){
            ensure(&st,BATCH_BYTES);uint8_t low[32];
            const uint8_t *base=st.buf+st.pos;
            uint64_t u[COPA_BASE_LIMBS][32];
            for(unsigned i=0;i<32;++i){unsigned offset=(i/8)*(64*COPA_BASE_LIMBS)+(i%8)*8;for(unsigned limb=0;limb<COPA_BASE_LIMBS;++limb)u[limb][i]=le64(base+offset+64*limb);}
            st.pos+=BASE_BYTES;st.avail-=BASE_BYTES;
            uint32_t acc=0;unsigned bits=0,idx=0;
            for(unsigned j=0;j<32;++j){while(bits<COPA_BASE_BITS){acc|=(uint32_t)st.buf[st.pos+idx++]<<bits;bits+=8;}low[j]=(uint8_t)(acc&((1u<<COPA_BASE_BITS)-1));acc>>=COPA_BASE_BITS;bits-=COPA_BASE_BITS;}
            st.pos+=Y_BYTES;st.avail-=Y_BYTES;
            uint32_t candidate[32],accepted[32];const uint8_t *tail=st.buf+st.pos;
            sample_batch_body(candidate,accepted,u,low,tail);
            /* Complete the fixed 32-proposal block before selecting outputs.
             * The following polynomial starts at the next full block boundary. */
            st.pos+=32*TAIL_BYTES;st.avail-=32*TAIL_BYTES;
            SC_ADD(SC_CANDIDATES,32);SC_ADD(SC_Y_USED,BATCH_BYTES);
#ifdef COPA_SAMPLER_PROFILE
            for(unsigned j=0;j<32;++j){SC_ADD(SC_ACCEPTED,accepted[j]);SC_ADD(SC_EXP_ACCEPTED,accepted[j]+((candidate[j]==0)&!accepted[j]));}
#endif
            SP_BEGIN(psel);for(unsigned j=0;j<32&&count<COPA_N;++j){unsigned sign=(signs[count/8]>>(count%8))&1;int32_t mask=-(int32_t)sign;y[poly][count]=((int32_t)candidate[j]^mask)-mask;if(accepted[j])++count;}SP_END(SP_SELECT,psel);

        }
    }
#endif
    copa_gaussian_leave(fp_env);
}
