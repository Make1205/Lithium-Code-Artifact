#include "ntt_crt.h"
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <x86intrin.h>
void copa_aux16_mul(int64_t *,const int32_t *,const int32_t *);
static uint64_t state=1;static uint32_t rng(void){state^=state<<13;state^=state>>7;state^=state<<17;return(uint32_t)state;}
static uint64_t ticks(void){unsigned aux;_mm_lfence();uint64_t t=__rdtscp(&aux);_mm_lfence();return t;}
static int cmp(const void *a,const void *b){uint64_t x=*(const uint64_t*)a,y=*(const uint64_t*)b;return(x>y)-(x<y);}
static volatile uint64_t sink;
int main(void){int32_t a[32][256],b[32][256];int64_t x[256],y[256];
 for(unsigned k=0;k<32;++k){for(unsigned j=0;j<256;++j){a[k][j]=(int32_t)(rng()%18001)-9000;b[k][j]=(int32_t)(rng()%6001)-3000;}if(k<2)for(unsigned j=0;j<256;++j){a[k][j]=k?-(1<<18):(1<<18);b[k][j]=1<<18;}copa_aux16_mul(x,a[k],b[k]);copa_ntt_crt_mul(y,a[k],b[k],256);assert(!memcmp(x,y,sizeof(x)));}
 uint64_t large[101],small[101];
 for(unsigned r=0;r<121;++r){uint64_t start=ticks();for(unsigned i=0;i<256;++i){copa_ntt_crt_mul(x,a[i&31],b[i&31],256);sink^=(uint64_t)x[i];}uint64_t t1=ticks()-start;
  start=ticks();for(unsigned i=0;i<256;++i){copa_aux16_mul(x,a[i&31],b[i&31]);sink^=(uint64_t)x[i];}uint64_t t2=ticks()-start;if(r>=20){large[r-20]=t1;small[r-20]=t2;}}
 qsort(large,101,8,cmp);qsort(small,101,8,cmp);printf("{\"degree\":256,\"samples\":101,\"products_per_sample\":256,\"two_large_cycles_per_product\":%.6f,\"three_16bit_cycles_per_product\":%.6f,\"three_over_two\":%.9f,\"checked_products\":32}\n",(double)large[50]/256,(double)small[50]/256,(double)small[50]/large[50]);return 0;
}
