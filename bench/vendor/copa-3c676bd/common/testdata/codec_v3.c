/* Frozen byte-codec reference from bc73ff8. */
#include "copa_rans.h"
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
_Static_assert(4+LOW_BYTES<=COPA_PAYLOAD,"fixed fields fit payload");
static int high(int32_t v){return v>=0?v/(1<<COPA_BLO):-((-v+(1<<COPA_BLO)-1)/(1<<COPA_BLO));}
static void putlow(uint8_t *buf,size_t *bit,uint32_t value){
 unsigned shift=(unsigned)(*bit&7);size_t byte=*bit>>3;uint32_t v=(value&((1u<<COPA_BLO)-1))<<shift;
 buf[byte]|=(uint8_t)v;if(shift+COPA_BLO>8)buf[byte+1]|=(uint8_t)(v>>8);*bit+=COPA_BLO;
}
static uint32_t getlow(const uint8_t *buf,size_t *bit){
 unsigned shift=(unsigned)(*bit&7);size_t byte=*bit>>3;uint32_t v=buf[byte];
 if(shift+COPA_BLO>8)v|=(uint32_t)buf[byte+1]<<8;
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
int copa_decode(int32_t z1[1+COPA_ELL][COPA_N],int32_t h[COPA_M][COPA_N],const uint8_t in[COPA_PAYLOAD])
{
    uint32_t state=0;for(unsigned j=0;j<4;++j)state|=(uint32_t)in[j]<<(8*j);
    if(state<LOWER||state>=STATE_LIMIT)return 0;
    size_t pos=4+LOW_BYTES,bit=0;
    /* The symbol count is fixed. Reaching LOWER in the middle is not EOF. */
    for(unsigned i=0;i<ZCOUNT+HCOUNT;++i){
        const uint32_t *cdf=i<ZCOUNT?copa_high_cdf:copa_hint_cdf;
        const uint16_t *lookup=i<ZCOUNT?copa_high_lookup:copa_hint_lookup;
        uint32_t slot=state&(SCALE-1);unsigned s;uint32_t start,freq;
#ifdef COPA_LEGACY_RANS
        s=symbol(cdf,lookup,slot);start=cdf[s];freq=cdf[s+1]-start;
#else
        (void)lookup;const uint64_t *fast=i<ZCOUNT?copa_high_fast:copa_hint_fast;
        uint64_t entry=fast[slot>>(COPA_RANS_BITS-COPA_FAST_LOOKUP_BITS)];s=(unsigned)entry&511;start=(uint32_t)(entry>>9)&(SCALE-1);freq=(uint32_t)(entry>>29);
        while(start+freq<=slot){++s;start=cdf[s];freq=cdf[s+1]-start;}
#endif
        if(i<ZCOUNT){
            int32_t v=((int32_t)s+COPA_HIGH_MIN)*(1<<COPA_BLO)+(int32_t)getlow(in+4,&bit);
            if(v<-COPA_BS||v>COPA_BS)return 0;
            (&z1[0][0])[i]=v;
        }else (&h[0][0])[i-ZCOUNT]=(int32_t)s+COPA_HINT_MIN;
        uint64_t next=(uint64_t)freq*(state>>COPA_RANS_BITS)+slot-start;
        if(next>=STATE_LIMIT)return 0;
        state=(uint32_t)next;
        while(state<LOWER){if(pos>=COPA_PAYLOAD)return 0;state=(state<<8)|in[pos++];}
    }
    if(state!=LOWER)return 0;
    for(size_t i=pos;i<COPA_PAYLOAD;++i)if(in[i])return 0;
    return 1;
}
