/* Optional Valgrind check: installed HighBits domain, secret input taint. */
#include "copa_round.h"
#include <valgrind/memcheck.h>
static volatile int32_t sink;
int main(void){
 for(int32_t v=0;v<COPA_Q;++v){int32_t x=v;
  VALGRIND_MAKE_MEM_UNDEFINED(&x,sizeof(x));int32_t y=copa_highbits(x);
  VALGRIND_MAKE_MEM_DEFINED(&y,sizeof(y));sink=y;
 }
 return 0;
}
