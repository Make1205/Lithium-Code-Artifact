/* Include the implementation to test its complete 128-bit quotient domain. */
#include "../ref/pack_pk.c"
#include <stdio.h>
static uint64_t seed=UINT64_C(0x724acfed91b38e05);
static uint64_t next(void){seed^=seed<<13;seed^=seed>>7;seed^=seed<<17;return seed;}
int main(void){const unsigned radix[]={3,1573,189,525};
 for(unsigned k=0;k<4;++k){uint64_t d=(uint64_t)radix[k]*radix[k]*radix[k]*radix[k];
  for(unsigned j=0;j<100000;++j){copa_pk_wide x=((copa_pk_wide)next()<<64)|next();
   if(j<64)x=j;
   if(j==64)x=~(copa_pk_wide)0;
   assert(divide_four(x,radix[k])==x/d);
  }
 }
 puts("PASS 400000 full-width mixed-radix reciprocal quotients");return 0;
}
