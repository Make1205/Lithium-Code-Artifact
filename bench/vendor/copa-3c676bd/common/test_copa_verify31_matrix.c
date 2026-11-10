#include "copa_api.h"
#define COPA_VERIFY31_CORE_ONLY
#include "copa_verify31.h"
#include <assert.h>
#include <stdio.h>
static uint64_t seed=UINT64_C(0x671943bc5d28e0af);
static unsigned next(void){seed^=seed<<13;seed^=seed>>7;seed^=seed<<17;return (unsigned)seed;}
int main(void){
 copa_public_ntt cache;int32_t a[COPA_M][1+COPA_ELL][COPA_N],z[1+COPA_ELL][COPA_N],out[COPA_M][COPA_N];
 for(unsigned test=0;test<100;++test){
  for(unsigned k=0;k<1+COPA_ELL;++k)for(unsigned i=0;i<COPA_N;++i)z[k][i]=test<2?0:(int32_t)(next()%17)-8;
  if(test<2)z[test][255]=(test?1:-1)*COPA_BV;
  for(unsigned r=0;r<COPA_M;++r)for(unsigned k=0;k<1+COPA_ELL;++k){
   for(unsigned i=0;i<COPA_N;++i)a[r][k][i]=(test<2?COPA_Q-1:(int32_t)(next()%COPA_Q))*(k?1:-1);
   q31_cache_impl(cache[r][k],a[r][k],k!=0);
  }
  copa_verify31_apply(out,cache,z);
  for(unsigned r=0;r<COPA_M;++r)for(unsigned i=0;i<COPA_N;++i){int64_t want=0;
   for(unsigned k=0;k<1+COPA_ELL;++k){for(unsigned j=0;j<=i;++j)want+=(int64_t)a[r][k][j]*z[k][i-j];for(unsigned j=i+1;j<COPA_N;++j)want-=(int64_t)a[r][k][j]*z[k][i+COPA_N-j];}
   want%=COPA_Q;if(want<0)want+=COPA_Q;assert(want==out[r][i]);
  }
 }
 printf("PASS Verify31 matrix mode=%d: 100 schoolbook matrices including sparse norm boundary\n",COPA_MODE);return 0;
}
