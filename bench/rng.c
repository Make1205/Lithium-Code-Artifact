#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "fips202.h"
extern int bench_testing;
void bench_fill(void *,size_t);
void bench_native_randombytes(uint8_t *,size_t);
void randombytes(uint8_t *out,size_t n){
 if(bench_testing){bench_fill(out,n);return;}
 /* Optional calibration entropy exactly follows COPA 3c676bd's wrapper.
    No algorithm call site or upstream provider is edited. */
 const char *seed=getenv("BENCH_SEED");
 if(!seed){bench_native_randombytes(out,n);return;}
 static uint64_t counter;uint8_t c[8];for(unsigned j=0;j<8;j++)c[j]=(uint8_t)(counter>>(8*j));counter++;
 keccak_state st;shake256_init(&st);shake256_absorb(&st,(const uint8_t*)seed,strlen(seed));shake256_absorb(&st,c,8);shake256_finalize(&st);shake256_squeeze(out,n,&st);
}
