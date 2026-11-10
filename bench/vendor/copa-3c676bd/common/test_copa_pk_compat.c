#include "pack_pk.h"
#include "quant.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>
int old_copa_pack_pk(uint8_t *,size_t,const uint8_t *,const int32_t *,unsigned);
int old_copa_unpack_pk(uint8_t *,int32_t *,const uint8_t *,size_t,unsigned);
static uint64_t seed=UINT64_C(0x9d3b41f52768ace1);
static unsigned next(void){seed^=seed<<13;seed^=seed>>7;seed^=seed<<17;return (unsigned)seed;}
int main(void){
 const unsigned modes[]={2,3,5,6};
 for(unsigned mi=0;mi<4;++mi){unsigned mode=modes[mi],q=copa_modulus(mode),n=copa_label_width(mode),size=(unsigned)copa_pk_bytes(mode),rows=copa_rows(mode);
  for(unsigned test=0;test<10000;++test){uint8_t rho[32],r0[32],r1[32],old[1480],now[1480];int32_t b[1024],b0[1024],b1[1024];for(unsigned j=0;j<32;++j)rho[j]=(uint8_t)next();
   for(unsigned k=0;k<64;++k){uint32_t label[16];int32_t x[16];for(unsigned j=0;j<n;++j)label[j]=next()%copa_label_radix(mode,j);assert(copa_unlab(x,label,mode));for(unsigned j=0;j<n;++j)b[(j/4)*256+(j%4)*64+k]=x[j];}
   if(mode==6)for(unsigned j=0;j<256;++j)b[512+j]=3*(int32_t)(next()%(q/3));
   assert(copa_pack_pk(now,size,rho,b,mode));assert(old_copa_pack_pk(old,size,rho,b,mode));assert(!memcmp(old,now,size));
   assert(copa_unpack_pk(r0,b0,now,size,mode)&&old_copa_unpack_pk(r1,b1,now,size,mode));assert(!memcmp(r0,r1,32)&&!memcmp(b0,b1,rows*256*4)&&!memcmp(b,b0,rows*256*4));
   now[next()%size]^=(uint8_t)(1u<<(next()%8));int a=copa_unpack_pk(r0,b0,now,size,mode),c=old_copa_unpack_pk(r1,b1,now,size,mode);assert(a==c);if(a)assert(!memcmp(r0,r1,32)&&!memcmp(b0,b1,rows*256*4));
  }
  printf("PASS PK mode=%u: 10000 exact old/new wires and canonical/mutated parses\n",mode);
 }
 return 0;
}
