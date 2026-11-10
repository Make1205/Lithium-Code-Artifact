#include "fips202.h"
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
extern int __real_PQCLEAN_randombytes(uint8_t *,size_t);
int __wrap_PQCLEAN_randombytes(uint8_t *out,size_t len){
 static uint64_t counter;const char *seed=getenv("BENCH_SEED");
 if(!seed)return __real_PQCLEAN_randombytes(out,len);
 size_t n=strlen(seed);uint8_t input[n+8];memcpy(input,seed,n);
 for(unsigned j=0;j<8;++j)input[n+j]=(uint8_t)(counter>>(8*j));++counter;
 shake256(out,len,input,n+8);return 0;
}
