#ifdef NDEBUG
#error Benchmarks require active correctness assertions
#endif
#include "api.h"
#include "prof.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
int main(int argc,char **argv){unsigned n=argc>1?(unsigned)strtoul(argv[1],0,10):1001;uint8_t pk[CRYPTO_PUBLICKEYBYTES],sk[CRYPTO_SECRETKEYBYTES],sig[CRYPTO_BYTES],m[59]={0};size_t len;
 assert(!crypto_sign_keypair(pk,sk));prof_reset();for(unsigned i=0;i<n;++i)assert(!crypto_sign_signature(sig,&len,m,sizeof(m),sk));
 unsigned ids[]={PT_SAMPLE_Y,PT_COMMIT,PT_HIGHBITS,PT_CHALLENGE,PT_SIGN,PT_NORMCHECK,PT_MAKEHINT,PT_PACK};uint64_t cycles=0;for(unsigned j=0;j<sizeof(ids)/sizeof(*ids);++j)cycles+=prof_cyc[ids[j]];
 printf("{\"scheme\":\"LITHIUM-%d\",\"samples\":%u,\"iterations\":%llu,\"iteration_cycles\":%llu,\"fixed_cycles\":%llu}\n",LITHIUM_MODE,n,(unsigned long long)prof_cnt[PT_SAMPLE_Y],(unsigned long long)cycles,(unsigned long long)(prof_cyc[PT_SETUP]+prof_cyc[PT_EXPAND_A0]));return 0;}
