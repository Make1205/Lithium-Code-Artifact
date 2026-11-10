/* Fixed-schedule kernel taint test; sampler rejection progress is outside it. */
#include "copa_sampler.h"
#include "copa_sampler_tables.h"
#include "copa_irs64.h"
#include <valgrind/memcheck.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
int main(int argc,char **argv){
 if(argc>1&&!strcmp(argv[1],"irs")){
  copa_irs64_context ctx;copa_irs64_prepare(&ctx,COPA_MIN_NORM+1);int64_t dot=3;
  uint64_t words[6]={1,2,3,4,5,6};unsigned minus,ok;
  VALGRIND_MAKE_MEM_UNDEFINED(&ctx,sizeof(ctx));VALGRIND_MAKE_MEM_UNDEFINED(&dot,sizeof(dot));VALGRIND_MAKE_MEM_UNDEFINED(words,sizeof(words));
  ok=copa_irs64_step(&minus,dot,&ctx,words);
  VALGRIND_MAKE_MEM_DEFINED(&ok,sizeof(ok));VALGRIND_MAKE_MEM_DEFINED(&minus,sizeof(minus));
 }else{
  uint64_t u[COPA_BASE_LIMBS32][32];uint8_t low[32],tail[256];uint32_t candidate[32],accepted[32];
  for(unsigned t=0;t<4;++t){
   for(unsigned l=0;l<COPA_BASE_LIMBS32;++l)for(unsigned j=0;j<32;++j)u[l][j]=(uint32_t)(j*193+l*4729+t*0x5fffffffU);
   for(unsigned j=0;j<32;++j)low[j]=(uint8_t)(j&((1u<<COPA_BASE_BITS)-1));
   memset(tail,(int)t*73,sizeof(tail));
   VALGRIND_MAKE_MEM_UNDEFINED(u,sizeof(u));VALGRIND_MAKE_MEM_UNDEFINED(low,sizeof(low));VALGRIND_MAKE_MEM_UNDEFINED(tail,sizeof(tail));
   for(unsigned j=0;j<32;++j)low[j]&=(1u<<COPA_BASE_BITS)-1;
   copa_sample_batch(candidate,accepted,u,low,tail);
   VALGRIND_MAKE_MEM_DEFINED(candidate,sizeof(candidate));VALGRIND_MAKE_MEM_DEFINED(accepted,sizeof(accepted));
  }
 }
 puts("PASS fixed-kernel taint trace complete");return 0;
}
