#include "copa_q120_sampler.h"
#include "fips202.h"
#include <string.h>
#include "copa_q120_base_tables.h"
#define BASE_BYTES (32*24)
#define Y_BYTES 28
#define TAIL_BYTES 16
#define BATCH_BYTES (BASE_BYTES+Y_BYTES+32*TAIL_BYTES)
#define UNIT_BYTES (8*SHAKE256_RATE)
typedef struct {keccak_state lane[8];uint8_t buf[BATCH_BYTES+UNIT_BYTES];size_t pos,avail;} stream;
static uint64_t le64(const uint8_t *b){uint64_t v=0;for(unsigned j=0;j<8;++j)v|=(uint64_t)b[j]<<(8*j);return v;}
static copa_u128 le128(const uint8_t *b){return (copa_u128)le64(b)|((copa_u128)le64(b+8)<<64);}
static void ensure(stream *s,size_t n)
{
    while(s->avail<n){if(s->avail)memmove(s->buf,s->buf+s->pos,s->avail);s->pos=0;for(unsigned j=0;j<8;++j)shake256_squeezeblocks(s->buf+s->avail+j*SHAKE256_RATE,1,&s->lane[j]);s->avail+=UNIT_BYTES;}
}
static unsigned lt192(const uint64_t a[3],const uint64_t b[3])
{unsigned lt=0,eq=1;for(int j=2;j>=0;--j){lt|=eq&(a[j]<b[j]);eq&=a[j]==b[j];}return lt;}
void copa_expand_y(int32_t y[COPA_K][COPA_N],const uint8_t seed[32],uint64_t attempt)
{
    stream st;st.pos=st.avail=0;uint8_t ext[42];memcpy(ext,seed,32);ext[32]=0x13;
    for(unsigned j=0;j<8;++j)ext[33+j]=(uint8_t)(attempt>>(8*j));
    for(unsigned j=0;j<8;++j){ext[41]=(uint8_t)j;shake256_absorb_once(&st.lane[j],ext,sizeof(ext));}
    for(unsigned poly=0;poly<COPA_K;++poly){
        uint8_t signs[32];ensure(&st,32);memcpy(signs,st.buf+st.pos,32);st.pos+=32;st.avail-=32;
        unsigned count=0;
        while(count<COPA_N){
            ensure(&st,BATCH_BYTES);unsigned x[32];uint8_t low[32];
            const uint8_t *base=st.buf+st.pos;
            for(unsigned i=0;i<32;++i){
                unsigned offset=(i/8)*192+(i%8)*8;uint64_t u[3];
                for(unsigned limb=0;limb<3;++limb)u[limb]=le64(base+offset+64*limb);
                x[i]=0;for(unsigned k=0;k<COPA_BASE_TAIL;++k)x[i]+=lt192(u,copa_base_rcdt[k]);
            }
            st.pos+=BASE_BYTES;st.avail-=BASE_BYTES;
            uint32_t acc=0;unsigned bits=0,idx=0;
            for(unsigned j=0;j<32;++j){while(bits<7){acc|=(uint32_t)st.buf[st.pos+idx++]<<bits;bits+=8;}low[j]=(uint8_t)(acc&127);acc>>=7;bits-=7;}
            st.pos+=Y_BYTES;st.avail-=Y_BYTES;
            for(unsigned j=0;j<32&&count<COPA_N;++j){
                uint64_t num=(uint64_t)low[j]*(low[j]+256*x[j])*10000;
                copa_u128 a=copa_ratio_fixed(num,2*COPA_R100*COPA_R100),p=copa_exp_negative(a);
                copa_u128 rand=le128(st.buf+st.pos);st.pos+=TAIL_BYTES;st.avail-=TAIL_BYTES;
                unsigned candidate=128*x[j]+low[j];
                unsigned accept=((rand>>1)&(COPA_ONE-1))<p;
                accept&=(candidate!=0)|(unsigned)(rand&1);
                unsigned sign=(signs[count/8]>>(count%8))&1;
                int32_t mask=-(int32_t)sign;y[poly][count]=((int32_t)candidate^mask)-mask;
                if(accept)++count; /* isochronous rejection-progress discipline inherited from Lithium */
            }
        }
    }
}
static copa_i128 select128(copa_i128 a,copa_i128 b,unsigned bit)
{copa_u128 m=0-(copa_u128)bit;return (copa_i128)(((copa_u128)a&m)|((copa_u128)b&~m));}
void copa_irs_table(copa_i128 table[COPA_IRS_TERMS],uint64_t norm_sq)
{
    for(unsigned k=0;k<COPA_IRS_TERMS;++k){
        copa_u128 a=copa_ratio_fixed(norm_sq*10000*k*k,2*COPA_R100*COPA_R100);
        copa_i128 v=(copa_i128)copa_exp_negative(a);
        table[k]=select128(-v,v,k&1);
    }
}
void copa_irs_thresholds(copa_i128 *minus,copa_i128 *accept,int64_t dot,uint64_t norm_sq,const copa_i128 table[COPA_IRS_TERMS])
{
    copa_i128 twice_gamma=(copa_i128)copa_ratio_fixed(norm_sq*10000,COPA_R100*COPA_R100);
    uint64_t neg=0-(uint64_t)(dot<0),absdot=((uint64_t)dot^neg)-neg;
    copa_i128 a=(copa_i128)copa_ratio_fixed(absdot*10000,COPA_R100*COPA_R100);
    copa_i128 ac=select128(a,twice_gamma,a<twice_gamma);
    copa_i128 w=(copa_i128)copa_exp_negative((copa_u128)a),wp=(copa_i128)copa_exp_positive((copa_u128)ac);
    copa_i128 sp=table[COPA_IRS_TERMS-1],sb=sp;
    for(int k=COPA_IRS_TERMS-2;k>=0;--k){sp=table[k]+copa_mul_fixed(sp,w);sb=table[k]+copa_mul_fixed(sb,wp);}
    copa_i128 spos=select128(sb,sp,dot<0),sneg=select128(sp,sb,dot<0);
    copa_i128 ta=select128(spos,(copa_i128)COPA_ONE-sneg,dot>=(int64_t)norm_sq);
    copa_i128 tab=select128(2*(copa_i128)COPA_ONE-spos-sneg,(copa_i128)COPA_ONE,(dot>=-(int64_t)norm_sq)&(dot<(int64_t)norm_sq));
    *minus=select128(0,ta,ta<0);*accept=select128(0,tab,tab<0);
}
int copa_irs(int32_t z[COPA_K][COPA_N],const uint8_t challenge[COPA_N],
             const int32_t s[COPA_K][COPA_N],uint64_t norm_sq,
             const copa_i128 table[COPA_IRS_TERMS],const uint8_t seed[32],uint64_t attempt)
{
    uint8_t ext[41];memcpy(ext,seed,32);ext[32]=0x15;for(unsigned j=0;j<8;++j)ext[33+j]=(uint8_t)(attempt>>(8*j));
    keccak_state rng;shake256_absorb_once(&rng,ext,sizeof(ext));
    unsigned accepted=1;
    for(unsigned j=0;j<COPA_N;++j){
        if(!challenge[j])continue; /* public challenge support */
        int32_t v[COPA_K][COPA_N];int64_t dot=0;
        for(unsigned r=0;r<COPA_K;++r)for(unsigned i=0;i<COPA_N;++i){v[r][i]=i>=j?s[r][i-j]:-s[r][i+COPA_N-j];dot+=(int64_t)z[r][i]*v[r][i];}
        copa_i128 ta,tab;copa_irs_thresholds(&ta,&tab,dot,norm_sq,table);
        uint8_t rb[16];shake256_squeeze(rb,sizeof(rb),&rng);
        copa_i128 rM=copa_mul_fixed((copa_i128)(le128(rb)&(COPA_ONE-1)),(copa_i128)COPA_IRS_M);
        unsigned minus=rM<ta;accepted&=rM<tab;int32_t mask=-(int32_t)minus;
        for(unsigned r=0;r<COPA_K;++r)for(unsigned i=0;i<COPA_N;++i)z[r][i]+=(v[r][i]^mask)-mask;
    }
    return (int)accepted;
}
