#include "copa_fp_simd.h"
#include "copa_scalar64.h"
#include <immintrin.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static uint64_t state=1;static uint64_t next(void){state^=state<<13;state^=state>>7;state^=state<<17;return state;}
static uint64_t ticks(void){unsigned aux;_mm_lfence();return __rdtscp(&aux);}
static __attribute__((noinline)) void hardware(uint64_t out[8],const uint64_t a[8],const uint64_t b[8]){for(unsigned i=0;i<8;++i){unsigned long long high;uint64_t lo=_mulx_u64(a[i],b[i],&high);out[i]=((uint64_t)high<<1|lo>>63)+((lo>>62)&1);}}
static __attribute__((noinline)) void vectors(uint64_t out[8],const uint64_t a[8],const uint64_t b[8]){copa_backend_mulq63x8(out,a,b);}
static int cmp(const void *a,const void *b){uint64_t x=*(const uint64_t *)a,y=*(const uint64_t *)b;return(x>y)-(x<y);}
static volatile uint64_t sink;
int main(void){uint64_t a[8],b[8],x[8],y[8];for(unsigned t=0;t<100000;++t){for(unsigned j=0;j<8;++j){a[j]=next();b[j]=next();}vectors(x,a,b);hardware(y,a,b);assert(!memcmp(x,y,sizeof(x)));}
 uint64_t timing[2][101];for(unsigned sample=0;sample<121;++sample)for(unsigned method=0;method<2;++method){for(unsigned j=0;j<8;++j){a[j]=next()>>1;b[j]=(UINT64_C(1)<<62)+(UINT64_C(1)<<61);}uint64_t start=ticks();for(unsigned t=0;t<10000;++t){if(method)hardware(x,a,b);else vectors(x,a,b);for(unsigned j=0;j<8;++j)a[j]=x[j];}uint64_t elapsed=ticks()-start;sink^=x[0];if(sample>=20)timing[method][sample-20]=elapsed;}
 qsort(timing[0],101,8,cmp);qsort(timing[1],101,8,cmp);printf("{\"vector_cycles_per_eight\":%.6f,\"hardware_mulx_cycles_per_eight\":%.6f,\"hardware_over_vector\":%.9f,\"random_full_width_comparisons\":800000}\n",timing[0][50]/10000.,timing[1][50]/10000.,(double)timing[1][50]/timing[0][50]);return 0;}
