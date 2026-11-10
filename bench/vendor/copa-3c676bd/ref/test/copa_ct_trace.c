/* A deterministic instruction/cache/branch trace region, independent of TSC. */
#include "ntt_crt.h"
#include <valgrind/callgrind.h>
#include <stdint.h>
#include <stdlib.h>
static uint64_t state=1;static volatile uint64_t sink;
static uint32_t rng(void){state^=state<<13;state^=state>>7;state^=state<<17;return(uint32_t)state;}
int main(int argc,char **argv){
 unsigned cls=argc>1?(unsigned)strtoul(argv[1],0,10):0;uint32_t mask=0u-cls;
 int32_t a[256],b[256];int64_t out[256];
 for(unsigned j=0;j<256;++j){int32_t x=(int32_t)(rng()%7)-3,y=(int32_t)(rng()%11)-5;a[j]=(int32_t)(((uint32_t)x&mask)|((uint32_t)(j==0)&~mask));b[j]=(int32_t)(((uint32_t)y&mask)|((uint32_t)(j==0)&~mask));}
 CALLGRIND_START_INSTRUMENTATION;CALLGRIND_TOGGLE_COLLECT;
 for(unsigned r=0;r<100;++r){copa_ntt_crt_mul(out,a,b,256);sink^=(uint64_t)out[r&255];}
 CALLGRIND_TOGGLE_COLLECT;CALLGRIND_STOP_INSTRUMENTATION;return 0;
}
