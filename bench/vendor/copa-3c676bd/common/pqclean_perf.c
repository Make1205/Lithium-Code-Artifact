#ifdef NDEBUG
#error Benchmarks require active correctness assertions
#endif
/* Byte-key, serialized-TSC benchmark adapter for PQClean clean signatures. */
#include "api.h"
#include "copa_bench_ticks.h"
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#define CAT_(a,b) a##b
#define CAT(a,b) CAT_(a,b)
#define API(n) CAT(PQ_PREFIX,n)
static int cmp(const void *a,const void *b){uint64_t x=*(const uint64_t *)a,y=*(const uint64_t *)b;return (x>y)-(x<y);}
int main(int argc,char **argv){unsigned n=argc>1?(unsigned)strtoul(argv[1],0,10):10001;assert(n>=101);unsigned warm=20;
 uint8_t pk[API(CRYPTO_PUBLICKEYBYTES)],sk[API(CRYPTO_SECRETKEYBYTES)],sig[API(CRYPTO_BYTES)],msg[59]={0};size_t len;
 uint64_t *kg=calloc(n,8),*sg=calloc(n,8),*vf=calloc(n,8);assert(kg&&sg&&vf);
 for(unsigned i=0;i<n+warm;++i){uint64_t t=copa_bench_ticks();assert(!API(crypto_sign_keypair)(pk,sk));uint64_t a=copa_bench_ticks()-t;
  t=copa_bench_ticks();assert(!API(crypto_sign_signature)(sig,&len,msg,sizeof(msg),sk));uint64_t b=copa_bench_ticks()-t;
  t=copa_bench_ticks();assert(!API(crypto_sign_verify)(sig,len,msg,sizeof(msg),pk));uint64_t c=copa_bench_ticks()-t;
  if(i>=warm){kg[i-warm]=a;sg[i-warm]=b;vf[i-warm]=c;}}
 qsort(kg,n,8,cmp);qsort(sg,n,8,cmp);qsort(vf,n,8,cmp);
 printf("{\"scheme\":\"%s\",\"key_input\":\"bytes\",\"samples\":%u,\"warmup\":%u,\"keygen\":%"PRIu64",\"sign\":%"PRIu64",\"verify\":%"PRIu64"}\n",API(CRYPTO_ALGNAME),n,warm,kg[n/2],sg[n/2],vf[n/2]);free(kg);free(sg);free(vf);return 0;}
