/* Deterministic SHAKE corpus; histogram collection is test-only. */
#include "copa_sampler.h"
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
extern uint64_t copa_y_generated,copa_y_consumed,copa_y_candidates;
int main(int argc,char **argv){
 uint64_t n=argc>1?strtoull(argv[1],0,10):UINT64_C(100000000),done=0,calls=0;
 unsigned cap=(16*COPA_R100+99)/100+1;uint64_t *hist=calloc(2*cap+1,sizeof(*hist));assert(hist);
 unsigned seed_tag=argc>2?(unsigned)strtoul(argv[2],0,10):0;
 uint8_t seed[32];for(unsigned j=0;j<32;++j)seed[j]=(uint8_t)(j*13+7)^(uint8_t)(seed_tag?(0x5a+7*j):0);
 while(done<n){int32_t y[COPA_K][COPA_N];copa_expand_y(y,seed,calls++);
  for(unsigned j=0;j<COPA_K*COPA_N&&done<n;++j){int32_t v=(&y[0][0])[j];assert(v>=-(int32_t)cap&&v<=(int32_t)cap);++hist[v+cap];++done;}}
 printf("{\"seed_tag\":%u,\"mode\":%d,\"samples\":%"PRIu64",\"calls\":%"PRIu64",\"squeezed_bytes\":%"PRIu64",\"consumed_bytes\":%"PRIu64",\"candidates\":%"PRIu64",\"offset\":%u,\"histogram\":[",seed_tag,COPA_MODE,n,calls,copa_y_generated,copa_y_consumed,copa_y_candidates,cap);
 for(unsigned j=0;j<2*cap+1;++j)printf("%s%"PRIu64,j?",":"",hist[j]);
 puts("]}");free(hist);
}
