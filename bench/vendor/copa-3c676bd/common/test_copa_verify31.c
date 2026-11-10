#include "copa_params.h"
#include "copa_prof.h"
#define COPA_VERIFY31_CORE_ONLY
#include "copa_verify31.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>
static uint64_t rng=UINT64_C(0x123456789abcdef);
static uint32_t next(void){rng^=rng<<13;rng^=rng>>7;rng^=rng<<17;return (uint32_t)rng;}
int main(void){
 for(unsigned t=0;t<1000;++t){int32_t a[256],b[256];uint32_t ca[2][256],cb[256],v[256];
  for(unsigned i=0;i<256;++i){a[i]=(int32_t)(next()%COPA_Q)-(t%2?COPA_Q/2:0);b[i]=(int32_t)(next()%33)-16;}
  q31_cache_impl(ca,a,t%2?0:1);q31_forward_ntt(cb,b);
  for(unsigned i=0;i<256;++i){assert(ca[0][i]<Q31_P);assert(ca[1][i]==((uint64_t)ca[0][i]<<32)/Q31_P);}
  for(unsigned i=0;i<256;++i){uint32_t x=cb[i],w=ca[0][i],ws=ca[1][i];uint32_t r=(uint32_t)(x*w-(uint32_t)(((uint64_t)x*ws)>>32)*Q31_P);v[i]=r-Q31_P*(r>=Q31_P);}
  q31_stages(v,1);
  for(unsigned i=0;i<256;++i){uint32_t x=v[i],w=q31_post[i],ws=q31_post_shoup[i],r=(uint32_t)(x*w-(uint32_t)(((uint64_t)x*ws)>>32)*Q31_P);r-=Q31_P*(r>=Q31_P);int64_t got=(int64_t)r-(int64_t)Q31_P*(r>Q31_P/2),want=0;
   for(unsigned j=0;j<=i;++j)want+=(int64_t)a[j]*b[i-j];
   for(unsigned j=i+1;j<256;++j)want-=(int64_t)a[j]*b[i+256-j];
   assert(got==want);
  }
 }
 puts("PASS 1000 Verify31 schoolbook products and 256000 cache reciprocals");return 0;
}
