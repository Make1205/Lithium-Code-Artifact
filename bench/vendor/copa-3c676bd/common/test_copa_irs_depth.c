#include "copa_irs64.h"
#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
void old_prepare(void *,uint64_t);
void old_prob(copa_fp63 *,copa_fp63 *,int64_t,const void *);
static unsigned at(unsigned j,unsigned terms){
#ifdef COPA_VEC_FP_HEADER
 return j/2+(j%2)*(terms/2);
#else
 (void)terms;return j;
#endif
}
int main(void){unsigned long pairs=0,zeros=0;
 for(unsigned norm=COPA_MIN_NORM;norm<=COPA_BK_SQ;++norm){
  uint8_t old[32768] __attribute__((aligned(64)));copa_irs64_context now;old_prepare(old,norm);copa_irs64_prepare(&now,norm);
  const uint64_t *c=(const uint64_t *)(old+offsetof(copa_irs64_context,c));
  for(unsigned j=0;j<COPA64_IRS_TERMS;++j)assert(now.c[at(j,COPA64_IRS_TERMS)]==c[at(j,OLD_T)]);
  for(unsigned j=COPA64_IRS_TERMS;j<OLD_T;++j){assert(c[at(j,OLD_T)]==0);++zeros;}
  int64_t dots[]={-(int64_t)norm-1,-(int64_t)norm,1-(int64_t)norm,-2,-1,0,1,2,(int64_t)norm-1,norm,(int64_t)norm+1};
  for(unsigned j=0;j<sizeof(dots)/sizeof(*dots);++j){copa_fp63 a,b,x,y;old_prob(&a,&b,dots[j],old);copa_irs64_probabilities(&x,&y,dots[j],&now);assert(a.m==x.m&&a.e==x.e&&b.m==y.m&&b.e==y.e);++pairs;}
 }
 printf("PASS IRS mode=%d all norms [%d,%d], %lu zero omitted coefficients, %lu identical transition pairs\n",COPA_MODE,COPA_MIN_NORM,COPA_BK_SQ,zeros,pairs);return 0;
}
