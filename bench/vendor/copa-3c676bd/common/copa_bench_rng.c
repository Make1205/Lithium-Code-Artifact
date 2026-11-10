#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "fips202.h"
extern void __real_randombytes(uint8_t *,size_t);
void __wrap_randombytes(uint8_t *out,size_t len){static unsigned init;static const char *seed;static uint64_t counter; if(!init){seed=getenv("BENCH_SEED");init=1;}if(!seed){__real_randombytes(out,len);return;}keccak_state st;uint8_t count[8];for(unsigned j=0;j<8;++j)count[j]=(uint8_t)(counter>>(8*j));++counter;shake256_init(&st);shake256_absorb(&st,(const uint8_t *)seed,strlen(seed));shake256_absorb(&st,count,8);shake256_finalize(&st);shake256_squeeze(out,len,&st);}
