#include "copa_rans.h"
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#ifdef NDEBUG
#error "Codec checks require assertions"
#endif
static uint64_t seed=1;
static uint32_t rnd(void){seed^=seed<<13;seed^=seed>>7;seed^=seed<<17;return (uint32_t)seed;}
int main(void)
{
    int32_t z[1+COPA_ELL][COPA_N],h[COPA_M][COPA_N],dz[1+COPA_ELL][COPA_N],dh[COPA_M][COPA_N];
    uint8_t packet[COPA_PAYLOAD],copy[COPA_PAYLOAD];size_t used;unsigned pass=0,overflow=0,zero_final_byte=0,canonical_fuzz=0;
    for(unsigned t=0;t<2000;++t){
        for(unsigned r=0;r<1+COPA_ELL;++r)for(unsigned i=0;i<COPA_N;++i)z[r][i]=t%4==0?0:(int32_t)(rnd()%101)-50;
        for(unsigned r=0;r<COPA_M;++r)for(unsigned i=0;i<COPA_N;++i)h[r][i]=t%4==0?0:((rnd()%8)<2?((rnd()&1)?1:-1):0);
        if(t%4==0)h[COPA_M-1][COPA_N-1]=COPA_Q/COPA_ALPHA/2; /* Force required terminal zero renormalization bytes. */
        if(!copa_encode(packet,&used,z,h)){++overflow;continue;}
        ++pass;assert(copa_decode(dz,dh,packet));assert(!memcmp(z,dz,sizeof(z))&&!memcmp(h,dh,sizeof(h)));
        size_t again;assert(copa_encode(copy,&again,dz,dh)&&again==used&&!memcmp(packet,copy,sizeof(copy)));
        zero_final_byte+=packet[used-1]==0;
        if(used<COPA_PAYLOAD){memcpy(copy,packet,sizeof(copy));copy[COPA_PAYLOAD-1]=1;assert(!copa_decode(dz,dh,copy));}
        memcpy(copy,packet,sizeof(copy));memset(copy,0,4);assert(!copa_decode(dz,dh,copy));
        memcpy(copy,packet,sizeof(copy));copy[3]|=0x80;assert(!copa_decode(dz,dh,copy));
        for(unsigned k=0;k<4;++k){
            memcpy(copy,packet,sizeof(copy));unsigned j=rnd()%COPA_PAYLOAD;copy[j]^=(uint8_t)(1u<<(rnd()%8));
            if(copa_decode(dz,dh,copy)){uint8_t recoded[COPA_PAYLOAD];assert(copa_encode(recoded,&again,dz,dh));assert(!memcmp(copy,recoded,sizeof(copy)));++canonical_fuzz;}
        }
    }
    assert(pass&&zero_final_byte);
    memset(z,0,sizeof(z));memset(h,0,sizeof(h));
    z[0][0]=COPA_BS;z[0][1]=-COPA_BS;
    assert(copa_encode(packet,&used,z,h)&&copa_decode(dz,dh,packet));
    assert(!memcmp(z,dz,sizeof(z))&&!memcmp(h,dh,sizeof(h)));
    z[0][0]=COPA_BS+1;assert(!copa_encode(packet,&used,z,h));
    z[0][0]=-COPA_BS-1;assert(!copa_encode(packet,&used,z,h));
    z[0][0]=0;h[0][0]=COPA_Q/COPA_ALPHA/2+1;assert(!copa_encode(packet,&used,z,h));
    printf("PASS codec-v2 mode=%d roundtrips=%u overflow_inputs=%u required_final_zero_cases=%u canonical_mutations=%u\n",COPA_MODE,pass,overflow,zero_final_byte,canonical_fuzz);
    return 0;
}
