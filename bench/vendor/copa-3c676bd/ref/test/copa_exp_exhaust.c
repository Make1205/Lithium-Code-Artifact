/* Every possible nonzero rational exponent of the integer-input IRS and
 * Gaussian kernels: num is always 10000 times an integer. The long-double
 * oracle has a separately reserved 2^-61 relative uncertainty. */
#include "copa_scalar64.h"
#include "copa_scalar64_tables.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
int main(void){
 uint64_t count=((uint64_t)copa64_r2.den*88+9999)/10000;long double worst=0,absolute=0;uint64_t point=0;
 for(uint64_t k=0;k<count;++k){copa_fp63 q=copa64_exp_ratio(k*10000,&copa64_r2);long double a=ldexpl((long double)q.m,q.e-63),p=expl(-(long double)(k*10000)/copa64_r2.den),err=fabsl(a/p-1);
  absolute=fmaxl(absolute,fabsl(a-p));if(err>worst){worst=err;point=k;}
 }
 printf("{\"mode\":%d,\"integer_exponents\":%llu,\"max_oracle_relative_error\":%.21Le,\"max_absolute_error\":%.21Le,\"oracle_uncertainty\":\"2^-61\",\"worst_integer\":%llu}\n",COPA_MODE,(unsigned long long)count,worst,absolute,(unsigned long long)point);
 assert(worst+0x1p-61L<0x1p-50L);return 0;
}
