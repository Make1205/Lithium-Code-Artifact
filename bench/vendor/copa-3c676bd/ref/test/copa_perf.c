/* Common measurement protocol for Copa and upstream LITHIUM backends. */
#include <assert.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fips202.h"
extern void __real_randombytes(uint8_t *,size_t);
void __wrap_randombytes(uint8_t *out,size_t len){static unsigned init;static const char *seed;static uint64_t counter; if(!init){seed=getenv("BENCH_SEED");init=1;}if(!seed){__real_randombytes(out,len);return;}keccak_state st;uint8_t count[8];for(unsigned j=0;j<8;++j)count[j]=(uint8_t)(counter>>(8*j));++counter;shake256_init(&st);shake256_absorb(&st,(const uint8_t *)seed,strlen(seed));shake256_absorb(&st,count,8);shake256_finalize(&st);shake256_squeeze(out,len,&st);}
#include <x86intrin.h>
#ifdef BENCH_LITHIUM
#include "api.h"
extern uint64_t lithium_bench_sign_attempts,lithium_bench_iteration_cycles;
#else
#include "copa_api.h"
#ifdef COPA_BENCH_ITERATIONS
extern uint64_t copa_bench_iteration_cycles,copa_bench_attempts;
#endif
#endif
#ifdef NDEBUG
#error Benchmark requires correctness assertions
#endif
static uint64_t ticks(void){unsigned aux;_mm_lfence();uint64_t t=__rdtscp(&aux);_mm_lfence();return t;}
static int cmp(const void *a,const void *b){uint64_t x=*(const uint64_t*)a,y=*(const uint64_t*)b;return(x>y)-(x<y);}
int main(int argc,char **argv){
 unsigned n=argc>1?(unsigned)strtoul(argv[1],NULL,10):1001,warm=20;assert(n>=101);
 uint64_t *kg=calloc(n,sizeof(*kg)),*sg=calloc(n,sizeof(*sg)),*vf=calloc(n,sizeof(*vf)),*cached_sg=calloc(n,8),*cached_vf=calloc(n,8);assert(kg&&sg&&vf&&cached_sg&&cached_vf);
 uint8_t msg[59]={0};uint64_t actual_iteration_cycles=0;uint64_t attempts=0,steps=0,sign_cycles=0,keygen_attempts=0,verify_cycles=0,cached_sign_cycles=0,cached_verify_cycles=0,cached_attempts=0;
#ifdef BENCH_LITHIUM
 uint8_t pk[CRYPTO_PUBLICKEYBYTES],sk[CRYPTO_SECRETKEYBYTES],sig[CRYPTO_BYTES];size_t len;
#else
 uint8_t pk[COPA_PK_BYTES],sig[COPA_SIG_BYTES],skbytes[COPA_SK_BYTES];copa_expanded_sk sk;copa_expanded_pk public_key;size_t len;
#endif
 for(unsigned i=0;i<n+warm;++i){uint64_t start=ticks();
#ifdef BENCH_LITHIUM
  assert(!crypto_sign_keypair(pk,sk));uint64_t a=ticks()-start;
  uint64_t before=lithium_bench_sign_attempts,iteration_before=lithium_bench_iteration_cycles;
  start=ticks();assert(!crypto_sign_signature(sig,&len,msg,sizeof(msg),sk));uint64_t b=ticks()-start;uint64_t actual_it=lithium_bench_iteration_cycles-iteration_before;
  start=ticks();assert(!crypto_sign_verify(sig,len,msg,sizeof(msg),pk));uint64_t c=ticks()-start;
  uint64_t tries=lithium_bench_sign_attempts-before,stps=tries*KAPPA,kgtries=0,cs=0,cv=0,ct=0;
#else
  copa_sign_stats stats;copa_keygen_stats kgstats;assert(!copa_keypair_bytes(pk,skbytes,&kgstats));uint64_t a=ticks()-start;
#ifdef COPA_BENCH_ITERATIONS
  uint64_t iteration_before=copa_bench_iteration_cycles;
#endif
  start=ticks();assert(!copa_sign_signature(sig,&len,msg,sizeof(msg),skbytes,sizeof(skbytes),&stats));uint64_t b=ticks()-start;
#ifdef COPA_BENCH_ITERATIONS
  uint64_t actual_it=copa_bench_iteration_cycles-iteration_before;
#else
  uint64_t actual_it=0;
#endif
  start=ticks();assert(!copa_verify(sig,sizeof(sig),msg,sizeof(msg),pk,sizeof(pk)));uint64_t c=ticks()-start;
  uint64_t tries=stats.attempts,stps=stats.irs_steps?stats.irs_steps:tries*COPA_TAU,kgtries=kgstats.attempts;
  assert(!crypto_sign_expand_sk(&sk,skbytes));
  assert(!crypto_sign_expand_pk(&public_key,pk));
#ifdef COPA_BENCH_ITERATIONS
  uint64_t cached_before=copa_bench_attempts;
#endif
  start=ticks();assert(!crypto_sign_signature_expanded(sig,&len,msg,sizeof(msg),&sk));uint64_t cs=ticks()-start;
  start=ticks();assert(!crypto_sign_verify_expanded(sig,sizeof(sig),msg,sizeof(msg),&public_key));uint64_t cv=ticks()-start;
#ifdef COPA_BENCH_ITERATIONS
  uint64_t ct=copa_bench_attempts-cached_before;
#else
  uint64_t ct=0; /* iteration statistics require benchmark instrumentation */
#endif
  crypto_sign_clear_expanded_sk(&sk);crypto_sign_clear_expanded_pk(&public_key);
#endif
  if(i>=warm){kg[i-warm]=a;sg[i-warm]=b;vf[i-warm]=c;attempts+=tries;actual_iteration_cycles+=actual_it;steps+=stps;sign_cycles+=b;keygen_attempts+=kgtries;verify_cycles+=c;cached_sg[i-warm]=cs;cached_vf[i-warm]=cv;cached_sign_cycles+=cs;cached_verify_cycles+=cv;cached_attempts+=ct;}
 }
 qsort(cached_sg,n,8,cmp);qsort(cached_vf,n,8,cmp);qsort(kg,n,sizeof(*kg),cmp);qsort(sg,n,sizeof(*sg),cmp);qsort(vf,n,sizeof(*vf),cmp);unsigned p99=(99*n+99)/100-1;
#ifdef BENCH_LITHIUM
 printf("{\"expanded_api\":false,\"scheme\":\"LITHIUM-%d\",",LITHIUM_MODE);
#else
 printf("{\"expanded_api\":true,\"scheme\":\"Copa-%s\",",COPA_MODE==2?"2":COPA_MODE==3?"3":COPA_MODE==5?"5a":"5b");
#endif
 printf("\"sign_actual_mean_cycles_per_iteration\":%.9f,\"fixed_mean_cycles\":%.9f,",attempts?(double)actual_iteration_cycles/attempts:0.0,(double)(sign_cycles-actual_iteration_cycles)/n);
 printf("\"key_input\":\"bytes\",\"includes_key_parsing_and_cache\":true,\"sign_mean\":%.9f,\"verify_mean\":%.9f,\"cached_sign\":%"PRIu64",\"cached_verify\":%"PRIu64",\"cached_sign_mean\":%.9f,\"cached_verify_mean\":%.9f,\"cached_mean_iterations\":%.9f,\"cached_mean_cycles_per_iteration\":%.9f,",(double)sign_cycles/n,(double)verify_cycles/n,cached_sg[n/2],cached_vf[n/2],(double)cached_sign_cycles/n,(double)cached_verify_cycles/n,(double)cached_attempts/n,cached_attempts?(double)cached_sign_cycles/cached_attempts:0.0);
 printf("\"verify_p99_p50\":%.9f,",(double)vf[p99]/vf[n/2]);
 printf("\"samples\":%u,\"warmup\":%u,\"message_bytes\":59,\"keygen\":%"PRIu64",\"sign\":%"PRIu64",\"verify\":%"PRIu64",\"keygen_p99_p50\":%.9f,\"sign_p99_p50\":%.9f,\"mean_iterations\":%.9f,\"irs_mean_steps_per_iteration\":%.9f,\"irs_mean_steps_per_signature\":%.9f,\"sign_amortized_median_cycles_per_iteration\":%.9f,\"sign_amortized_mean_cycles_per_iteration\":%.9f,\"keygen_mean_attempts\":%.9f}\n",n,warm,kg[n/2],sg[n/2],vf[n/2],(double)kg[p99]/kg[n/2],(double)sg[p99]/sg[n/2],(double)attempts/n,(double)steps/attempts,(double)steps/n,(double)sg[n/2]*n/attempts,(double)sign_cycles/attempts,(double)keygen_attempts/n);free(kg);free(sg);free(vf);free(cached_sg);free(cached_vf);return 0;
}
