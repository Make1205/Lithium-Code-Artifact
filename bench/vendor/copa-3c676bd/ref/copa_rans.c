#include "copa_rans.h"
#ifdef COPA_VEC_LATTICE
#include "../common/copa_simd.h"
#endif
#include <string.h>
#include "copa_rans_tables.h"
#include "copa_rans_lookup.h"
#ifndef COPA_LEGACY_RANS
#include "copa_rans_fast.h"
#endif
#define LOWER COPA_RANS_LOWER
#define SCALE (UINT32_C(1)<<COPA_RANS_BITS)
#define STATE_LIMIT ((uint64_t)LOWER<<8)
#define ZCOUNT ((1+COPA_ELL)*COPA_N)
#define HCOUNT (COPA_M*COPA_N)
#define LOW_BYTES (ZCOUNT*COPA_BLO/8)
_Static_assert(COPA_RANS_BITS==COPA_RANS_TABLE_BITS,"frequency precision mismatch");
_Static_assert(COPA_BLO==COPA_RANS_TABLE_LOW_BITS,"raw-low precision mismatch");
_Static_assert(SCALE<=LOWER,"normalization bound must cover frequency scale");
_Static_assert(STATE_LIMIT<=UINT32_MAX,"normalized states fit uint32_t");
_Static_assert(20+LOW_BYTES<=COPA_PAYLOAD,"fixed fields fit payload");
static int high(int32_t v){return v>=0?v/(1<<COPA_BLO):-((-v+(1<<COPA_BLO)-1)/(1<<COPA_BLO));}
static void putlow(uint8_t *buf,size_t *bit,uint32_t value){
 unsigned shift=(unsigned)(*bit&7);size_t byte=*bit>>3;uint32_t v=(value&((1u<<COPA_BLO)-1))<<shift;
 buf[byte]|=(uint8_t)v;if(shift+COPA_BLO>8)buf[byte+1]|=(uint8_t)(v>>8);*bit+=COPA_BLO;
}
static uint32_t getlow(const uint8_t *buf,size_t *bit){
 unsigned shift=(unsigned)(*bit&7);size_t byte=*bit>>3;uint32_t v=buf[byte]|((uint32_t)buf[byte+1]<<8);
 *bit+=COPA_BLO;return (v>>shift)&((1u<<COPA_BLO)-1);
}

int copa_encode(uint8_t out[COPA_PAYLOAD],size_t *used,const int32_t z1[1+COPA_ELL][COPA_N],const int32_t h[COPA_M][COPA_N])
{
    for(unsigned i=0;i<ZCOUNT;++i)if((&z1[0][0])[i]<-COPA_BS||(&z1[0][0])[i]>COPA_BS)return 0;
    uint8_t bytes[COPA_PAYLOAD];uint8_t *ptr=bytes+sizeof(bytes);uint32_t state=LOWER;
    for(int i=ZCOUNT+HCOUNT-1;i>=0;--i){
        const uint32_t *cdf;int index;
        if(i>=ZCOUNT){index=(&h[0][0])[i-ZCOUNT]-COPA_HINT_MIN;if(index<0||index>=COPA_HINT_COUNT)return 0;cdf=copa_hint_cdf;}
        else{index=high((&z1[0][0])[i])-COPA_HIGH_MIN;if(index<0||index>=COPA_HIGH_COUNT)return 0;cdf=copa_high_cdf;}
        uint32_t start=cdf[index],freq=cdf[index+1]-start;
        uint64_t maximum=((uint64_t)(LOWER>>COPA_RANS_BITS)<<8)*freq;
        while(state>=maximum){if(ptr==bytes)return 0;*--ptr=(uint8_t)state;state>>=8;}
#ifdef COPA_LEGACY_RANS
        uint64_t next=((uint64_t)(state/freq)<<COPA_RANS_BITS)+(state%freq)+start;
#else
        const uint64_t *reciprocal=i>=ZCOUNT?copa_hint_reciprocal:copa_high_reciprocal;
        uint32_t quotient=(uint32_t)(((uint64_t)state*reciprocal[index])>>32),rem=state-quotient*freq,correction=rem>=freq;
        quotient+=correction;rem-=correction*freq;
        uint64_t next=((uint64_t)quotient<<COPA_RANS_BITS)+rem+start;
#endif
        if(next<LOWER||next>=STATE_LIMIT)return 0;
        state=(uint32_t)next;
    }
    size_t n=(size_t)(bytes+sizeof(bytes)-ptr),len=4+LOW_BYTES+n;
    if(len>COPA_PAYLOAD)return 0;
    memset(out,0,COPA_PAYLOAD);
    for(unsigned j=0;j<4;++j)out[j]=(uint8_t)(state>>(8*j));
    size_t bit=0;for(unsigned i=0;i<ZCOUNT;++i)putlow(out+4,&bit,(uint32_t)(&z1[0][0])[i]);
    memcpy(out+4+LOW_BYTES,ptr,n);if(used)*used=len;return 1;
}
#ifdef COPA_LEGACY_RANS
static unsigned symbol(const uint32_t *cdf,const uint16_t *lookup,uint32_t value)
{
    unsigned s=lookup[value>>(COPA_RANS_BITS-COPA_RANS_LOOKUP_BITS)];
    while(cdf[s+1]<=value)++s;
    return s;
}
#endif
/* The zero hint is the common case; its fixed interval needs no lookup.
 * Separate coefficient/hint loops let the compiler specialize both models. */
static inline __attribute__((always_inline)) int decode_one(uint32_t *state,size_t *pos,unsigned *symbol_out,const uint8_t *in,const uint32_t *cdf,const void *table,unsigned hint){
 uint32_t slot=*state&(SCALE-1),start,freq;unsigned sym;
 if(hint && (uint32_t)(slot-copa_hint_cdf[-COPA_HINT_MIN])<copa_hint_cdf[1-COPA_HINT_MIN]-copa_hint_cdf[-COPA_HINT_MIN]){
  sym=(unsigned)-COPA_HINT_MIN;start=copa_hint_cdf[sym];freq=copa_hint_cdf[sym+1]-start;
 }else{
#ifdef COPA_LEGACY_RANS
  sym=symbol(cdf,(const uint16_t *)table,slot);start=cdf[sym];freq=cdf[sym+1]-start;
#else
  const uint64_t *fast=(const uint64_t *)table;uint64_t entry=fast[slot>>(COPA_RANS_BITS-COPA_FAST_LOOKUP_BITS)];sym=(unsigned)entry&511;start=(uint32_t)(entry>>9)&(SCALE-1);freq=(uint32_t)(entry>>29);
  while(start+freq<=slot){++sym;start=cdf[sym];freq=cdf[sym+1]-start;}
#endif
 }
 /* slot-start<freq<=SCALE and normalized state<STATE_LIMIT imply
  * next<freq*(state/SCALE+1)<=STATE_LIMIT. No 64-bit product is needed. */
 uint32_t next=freq*(*state>>COPA_RANS_BITS)+slot-start;
 while(next<LOWER){if(*pos>=COPA_PAYLOAD)return 0;next=(next<<8)|in[(*pos)++];}
 *state=next;*symbol_out=sym;return 1;
}
static int decode_low_fields(int32_t *z,const uint8_t *in){
#ifdef COPA_VEC_LATTICE
 uint8_t shuffle[4*COPAV_W];uint32_t shifts[COPAV_W];
 for(unsigned j=0;j<COPAV_W;++j){shuffle[4*j]=(uint8_t)(j*COPA_BLO/8);shuffle[4*j+1]=shuffle[4*j]+1;shuffle[4*j+2]=shuffle[4*j+3]=0x80;shifts[j]=j*COPA_BLO%8;}
 copav_vec bad=cv_set(0),indices=cv_load(shuffle),shift=cv_load(shifts);
 for(unsigned i=0;i<ZCOUNT;i+=COPAV_W){__m128i raw=_mm_loadu_si128((const __m128i *)(in+i*COPA_BLO/8));copav_vec low;
#if COPA_SIMD_BITS==256
  low=_mm256_shuffle_epi8(_mm256_broadcastsi128_si256(raw),indices);low=_mm256_srlv_epi32(low,shift);
#else
  low=_mm512_shuffle_epi8(_mm512_broadcast_i32x4(raw),indices);low=_mm512_srlv_epi32(low,shift);
#endif
  copav_vec value=cv_add(cv_load(z+i),cv_and(low,cv_set((1u<<COPA_BLO)-1)));bad=cv_or(bad,cv_or(cv_gt(value,cv_set(COPA_BS)),cv_gt(cv_set(-COPA_BS),value)));cv_store(z+i,value);
 }
#if COPA_SIMD_BITS==256
 return _mm256_movemask_epi8(bad)==0;
#else
 return _mm512_cmpneq_epi32_mask(bad,cv_set(0))==0;
#endif
#else
 size_t bit=0;for(unsigned i=0;i<ZCOUNT;++i){int32_t value=z[i]+(int32_t)getlow(in,&bit);if(value<-COPA_BS||value>COPA_BS)return 0;z[i]=value;}return 1;
#endif
}
int copa_decode(int32_t z1[1+COPA_ELL][COPA_N],int32_t h[COPA_M][COPA_N],const uint8_t in[COPA_PAYLOAD]){
 uint32_t state=0;for(unsigned j=0;j<4;++j)state|=(uint32_t)in[j]<<(8*j);
 if(state<LOWER||state>=STATE_LIMIT)return 0;
 size_t pos=4+LOW_BYTES;unsigned sym;
#ifdef COPA_LEGACY_RANS
 const void *high_table=copa_high_lookup,*hint_table=copa_hint_lookup;
#else
 const void *high_table=copa_high_fast,*hint_table=copa_hint_fast;
#endif
 for(unsigned i=0;i<ZCOUNT;++i){
  if(!decode_one(&state,&pos,&sym,in,copa_high_cdf,high_table,0))return 0;
  (&z1[0][0])[i]=((int32_t)sym+COPA_HIGH_MIN)*(1<<COPA_BLO);
 }
 for(unsigned i=0;i<HCOUNT;++i){if(!decode_one(&state,&pos,&sym,in,copa_hint_cdf,hint_table,1))return 0;(&h[0][0])[i]=(int32_t)sym+COPA_HINT_MIN;}
 if(state!=LOWER||!decode_low_fields(&z1[0][0],in+4))return 0;
 for(size_t i=pos;i<COPA_PAYLOAD;++i)if(in[i])return 0;
 return 1;
}
