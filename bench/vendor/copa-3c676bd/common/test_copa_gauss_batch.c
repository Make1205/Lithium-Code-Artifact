#include "copa_sampler.h"
#include "copa_sampler_tables.h"
#include "copa_gauss_simd.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>
void copa_ref_sample_batch(uint32_t[32],uint32_t[32],const uint64_t[][32],const uint8_t[32],const uint8_t[256]);
static uint64_t state=1;
static uint64_t rng(void){state^=state<<13;state^=state>>7;state^=state<<17;return state;}
int main(void){
 for(unsigned t=0;t<10000;++t){uint64_t u[COPA_BASE_LIMBS32][32];uint8_t low[32],tail[256];uint32_t ref[32],vec[32],ra[32],va[32];
  for(unsigned j=0;j<32;++j){for(unsigned k=0;k<COPA_BASE_LIMBS32;++k)u[k][j]=t==0?0:t==1?UINT64_MAX:rng();low[j]=(uint8_t)((t<2?UINT64_MAX:rng())&((1u<<COPA_BASE_BITS)-1));}
  for(unsigned j=0;j<256;++j)tail[j]=(uint8_t)(t==0?0:t==1?255:rng());
  copa_ref_sample_batch(ref,ra,u,low,tail);copa_gaussian_env env=copa_gaussian_enter();copa_backend_gaussian_batch(vec,va,u,low,tail);copa_gaussian_leave(env);
  assert(!memcmp(ref,vec,sizeof(ref)));assert(!memcmp(ra,va,sizeof(ra)));
 }
 printf("PASS mode=%d 320000 Gaussian proposals and accept bits match scalar oracle\n",COPA_MODE);return 0;
}
