/* Positive-pair lower bound on S(a+2 gamma), all a>=0, all integer norms.
 * Coefficient enclosures and roundoff padding are deliberately much wider
 * than the analytic implementation errors. This is an offline certificate. */
#include "copa_irs64.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
int main(void){
 long double global=1;unsigned worst=0;long double max_coefficient_error=0;unsigned long boxes=0;
 for(unsigned norm=COPA_MIN_NORM;norm<=COPA_BK_SQ;++norm){
  copa_irs64_context ctx;copa_irs64_prepare(&ctx,norm);
  long double lo[COPA64_IRS_TERMS],hi[COPA64_IRS_TERMS];
  long double coefficient_error=0;
  for(unsigned j=0;j<COPA64_IRS_TERMS;++j){long double v=ldexpl((long double)ctx.c[j],-63);long double exact=expl(-(long double)norm*10000*j*j/(2*(long double)COPA_R100*COPA_R100));coefficient_error+=fabsl(v-exact)+0x1p-61L*exact;lo[j]=fmaxl(0,v-0x1p-45L);hi[j]=v+0x1p-45L;}
  max_coefficient_error=fmaxl(max_coefficient_error,coefficient_error);
  unsigned parts=(unsigned)ceill(2048/sqrtl(norm));
  for(unsigned b=0;b<parts;++b){long double wl=(long double)b/parts,wh=(long double)(b+1)/parts,power=1,lower=0;
   for(unsigned j=0;j<COPA64_IRS_TERMS;j+=2){long double term=fmaxl(0,lo[j]-hi[j+1]*wh-0x1p-50L);lower+=power*term;power*=wl*wl;}
   lower-=0x1p-30L;
   if(lower<global){global=lower;worst=norm;}++boxes;
  }
 }
 printf("{\"mode\":%d,\"norms\":%d,\"minimum_squared_norm\":%d,\"boxes\":%lu,\"S_uniform_lower_bound\":%.18Lf,\"worst_norm\":%u,\"coefficient_absolute_error_sum_bound\":%.21Le}\n",COPA_MODE,COPA_BK_SQ-COPA_MIN_NORM+1,COPA_MIN_NORM,boxes,global,worst,max_coefficient_error);assert(global>0.4L);return 0;
}
