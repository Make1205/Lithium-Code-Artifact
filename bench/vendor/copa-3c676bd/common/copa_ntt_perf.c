#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <x86intrin.h>
#ifdef BENCH_LITHIUM
#include "ntt.h"
#include "ntt_avx_decls.h"
#if LITHIUM_MODE==260
#define BENCH_FWD(a) ntt_avx((int16_t *)(a),ntt_qdata,ntt_zetas_fwd,ntt_cross_fwd,ntt_ninv)
#define BENCH_INV(a) invntt_tomont_avx((int16_t *)(a),ntt_qdata,ntt_zetas_inv,ntt_cross_inv,ntt_ninv)
#else
#define BENCH_FWD(a) ntt_avx((int16_t *)(a),ntt_qdata,ntt_zetas_fwd,ntt_z0,ntt_scale)
#define BENCH_INV(a) invntt_tomont_avx((int16_t *)(a),ntt_qdata,ntt_zetas_inv,ntt_z0inv,ntt_scale)
#endif
#else
#include "copa_ntt_vectors.c"
#endif
static uint64_t tick(void){unsigned aux;_mm_lfence();uint64_t t=__rdtscp(&aux);_mm_lfence();return t;}
static int compare(const void *a,const void *b){uint64_t x=*(const uint64_t*)a,y=*(const uint64_t*)b;return(x>y)-(x<y);}
static volatile uint64_t sink;
int main(void){enum{COUNT=10001};uint64_t f[COUNT],v[COUNT],c[COUNT],empty[COUNT];
#ifdef BENCH_LITHIUM
 uint16_t input[N] __attribute__((aligned(64))),transformed[N] __attribute__((aligned(64))),buffer[N] __attribute__((aligned(64)));for(unsigned j=0;j<N;++j)input[j]=(j*177+29)%Q;memcpy(transformed,input,sizeof(input));BENCH_FWD(transformed);
#else
 int32_t input[256];uint32_t transformed[2][256],buffer[2][256];int64_t output[256];for(unsigned j=0;j<256;++j)input[j]=(j*177+29)%6144;copa_ntt_forward(transformed,input,256);
#endif
 for(unsigned j=0;j<COUNT+100;++j){uint64_t a=tick(),e=tick()-a;
#ifdef BENCH_LITHIUM
  memcpy(buffer,input,sizeof(buffer));a=tick();BENCH_FWD(buffer);uint64_t x=tick()-a;sink+=buffer[0];memcpy(buffer,transformed,sizeof(buffer));a=tick();BENCH_INV(buffer);uint64_t y=tick()-a,z=0;sink+=buffer[0];
#else
  a=tick();copa_ntt_forward(buffer,input,256);uint64_t x=tick()-a;sink+=buffer[0][0];a=tick();inverse_residues(buffer,transformed,256);uint64_t y=tick()-a;sink+=buffer[0][0];a=tick();copa_ntt_inverse_crt(output,transformed,256);uint64_t z=tick()-a;sink+=(uint64_t)output[0];
#endif
  if(j>=100){f[j-100]=x;v[j-100]=y;c[j-100]=z;empty[j-100]=e;}
 }
 qsort(f,COUNT,8,compare);qsort(v,COUNT,8,compare);qsort(c,COUNT,8,compare);qsort(empty,COUNT,8,compare);
#ifdef BENCH_LITHIUM
 printf("{\"scheme\":\"LITHIUM-%d\",\"n\":%d,\"primes\":1,",LITHIUM_MODE,N);
#else
 printf("{\"scheme\":\"Copa\",\"n\":256,\"primes\":2,");
#endif
 printf("\"samples\":%d,\"timestamp_overhead\":%llu,\"forward_raw\":%llu,\"inverse_raw\":%llu,\"inverse_and_CRT_raw\":%llu,\"forward_net\":%llu,\"inverse_net\":%llu}\n",COUNT,(unsigned long long)empty[COUNT/2],(unsigned long long)f[COUNT/2],(unsigned long long)v[COUNT/2],(unsigned long long)c[COUNT/2],(unsigned long long)(f[COUNT/2]-empty[COUNT/2]),(unsigned long long)(v[COUNT/2]-empty[COUNT/2]));return 0;
}
