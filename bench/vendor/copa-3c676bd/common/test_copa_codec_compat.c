#include "copa_rans.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>
int old_copa_encode(uint8_t *,size_t *,const int32_t (*)[COPA_N],const int32_t (*)[COPA_N]);
int old_copa_decode(int32_t (*)[COPA_N],int32_t (*)[COPA_N],const uint8_t *);
static uint64_t seed=UINT64_C(0x51c92e784abd630f);
static unsigned next(void){seed^=seed<<13;seed^=seed>>7;seed^=seed<<17;return (unsigned)seed;}
int main(void){
 int32_t z[1+COPA_ELL][COPA_N],h[COPA_M][COPA_N],a[1+COPA_ELL][COPA_N],b[1+COPA_ELL][COPA_N],ha[COPA_M][COPA_N],hb[COPA_M][COPA_N];uint8_t in[COPA_PAYLOAD],old[COPA_PAYLOAD],base[COPA_PAYLOAD];unsigned accepted=0,have_base=0;
 for(unsigned test=0;test<100000;++test){
  if(test%100==0){for(unsigned i=0;i<(1+COPA_ELL)*COPA_N;++i)(&z[0][0])[i]=(int32_t)(next()%41)-20;memset(h,0,sizeof(h));size_t n0=0,n1=0;int p=copa_encode(in,&n0,z,h),q=old_copa_encode(old,&n1,z,h);assert(p==q);if(p){assert(n0==n1&&!memcmp(in,old,sizeof(in)));memcpy(base,in,sizeof(in));have_base=1;}else continue;}
  else{if(!have_base)continue;memcpy(in,base,sizeof(in));in[next()%sizeof(in)]^=(uint8_t)(1u<<(next()%8));}
  int x=copa_decode(a,ha,in),y=old_copa_decode(b,hb,in);assert(x==y);if(x){++accepted;assert(!memcmp(a,b,sizeof(a))&&!memcmp(ha,hb,sizeof(ha)));}
 }
 assert(accepted>=1000);printf("PASS codec mode=%d: 100000 valid/mutated old/new parses; accepted=%u\n",COPA_MODE,accepted);return 0;
}
