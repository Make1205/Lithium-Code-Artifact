#include "copa_sampler.h"
#include "copa_sampler_tables.h"
#include "copa_gaussian_fp.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>
int main(void){unsigned old=_mm_getcsr(),points=0;
 for(unsigned x=0;x<=COPA_BASE_TAIL;++x)for(unsigned y=0;y<(1u<<COPA_BASE_BITS);++y){
  uint64_t u[COPA_BASE_LIMBS32][32];uint8_t low[32],tail[256];uint32_t expected[32],ea[32],actual[32],aa[32];
  _mm_setcsr(0x1f80);uint64_t A=copa_gaussian_threshold56((uint64_t)y*(y+(2u<<COPA_BASE_BITS)*x)*10000);
  for(unsigned j=0;j<32;++j){for(unsigned l=0;l<COPA_BASE_LIMBS32;++l)u[l][j]=x<COPA_BASE_TAIL?copa_base_rcdt[x][l]:0;low[j]=(uint8_t)y;uint64_t r=(A+j%3-1)&UINT64_C(0x00ffffffffffffff);for(unsigned k=0;k<7;++k)tail[8*j+k]=(uint8_t)(r>>(8*k));tail[8*j+7]=(uint8_t)(j&1);}
  copa_sample_batch(expected,ea,u,low,tail);
  for(unsigned j=0;j<32;++j)assert(expected[j]==(1u<<COPA_BASE_BITS)*x+y);
  for(unsigned rc=0;rc<4;++rc)for(unsigned ftz=0;ftz<2;++ftz)for(unsigned daz=0;daz<2;++daz){
   unsigned requested=0x1f80|(rc<<13)|(ftz<<15)|(daz<<6);_mm_setcsr(requested);copa_sample_batch(actual,aa,u,low,tail);
   assert(_mm_getcsr()==requested);assert(!memcmp(expected,actual,sizeof(actual)));assert(!memcmp(ea,aa,sizeof(aa)));++points;
  }
 }
 _mm_setcsr(old);printf("PASS mode=%d actual sampler batches: %u proposal/environment pairs, restore exact\n",COPA_MODE,points);
}
