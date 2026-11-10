#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "../ref/quant.c"
static uint64_t rng=41;
static uint32_t next(void){rng^=rng<<13;rng^=rng>>7;rng^=rng<<17;return(uint32_t)rng;}
int main(void){int32_t in[COPA_M*256],out[COPA_M*256],error[COPA_M*256];
 for(unsigned trial=0;trial<1000;++trial){for(unsigned j=0;j<COPA_M*256;++j)in[j]=trial<4?(int32_t)trial:(int32_t)(next()%COPA_Q);copa_quant_all(out,error,in);
  for(unsigned k=0;k<64;++k){int32_t x[16],y[16];for(unsigned j=0;j<COPA_M*4;++j)x[j]=in[(j/4)*256+(j%4)*64+k];copa_quant_block(y,x,COPA_MODE);for(unsigned j=0;j<COPA_M*4;++j){unsigned idx=(j/4)*256+(j%4)*64+k;if(y[j]!=out[idx]){printf("FAIL mode=%d trial=%u block=%u coord=%u scalar=%d vector=%d\n",COPA_MODE,trial,k,j,y[j],out[idx]);return 1;}assert(error[idx]==out[idx]-in[idx]);}}
 }
 printf("PASS mode=%d lattice 64000 blocks, all ties and codewords\n",COPA_MODE);return 0;
}
