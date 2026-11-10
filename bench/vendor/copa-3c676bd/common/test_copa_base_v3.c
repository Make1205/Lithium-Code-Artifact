/* Wide integers are an independent test oracle, not sampler arithmetic. */
#include "copa_sampler_tables.h"
#include "copa_sampler.h"
#ifdef COPA_SIMD_BITS
#include "copa_gauss_simd.h"
#endif
#include <assert.h>
#include <stdio.h>
#include <string.h>
typedef unsigned __int128 wide;
static uint64_t state=1;
static uint32_t rng(void){state^=state<<13;state^=state>>7;state^=state<<17;return (uint32_t)state;}
static wide entry(unsigned row){wide v=0;for(unsigned l=0;l<COPA_BASE_LIMBS32;++l)v|=(wide)copa_base_rcdt[row][l]<<(32*l);return v;}
static void check(uint64_t u[][32]){
 uint32_t out[32];
#ifdef COPA_SIMD_BITS
 copa_base96_vectors(out,u);
#else
 uint32_t accepted[32];uint8_t low[32]={0},tail[256]={0};copa_sample_batch(out,accepted,u,low,tail);for(unsigned j=0;j<32;++j)out[j]>>=COPA_BASE_BITS;
#endif
 for(unsigned j=0;j<32;++j){wide v=0;for(unsigned l=0;l<COPA_BASE_LIMBS32;++l)v|=(wide)(uint32_t)u[l][j]<<(32*l);unsigned x=0;for(unsigned k=0;k<COPA_BASE_TAIL;++k)x+=v<entry(k);assert(out[j]==x);}
}
int main(void){uint64_t u[COPA_BASE_LIMBS32][32];
 for(unsigned k=0;k<COPA_BASE_TAIL;++k){for(unsigned j=0;j<32;++j){wide v=entry(k)+(int)(j%3)-1;for(unsigned l=0;l<COPA_BASE_LIMBS32;++l)u[l][j]=(uint32_t)(v>>(32*l));}check(u);}
 for(unsigned t=0;t<1000;++t){for(unsigned l=0;l<COPA_BASE_LIMBS32;++l)for(unsigned j=0;j<32;++j)u[l][j]=t<2?(t?UINT32_MAX:0):rng();check(u);}
 printf("PASS mode=%d: all RCDT threshold +/-1 boundaries, zero/max and 32000 random vectors\n",COPA_MODE);}
