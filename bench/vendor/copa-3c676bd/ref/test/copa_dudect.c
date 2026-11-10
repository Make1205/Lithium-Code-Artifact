#define _GNU_SOURCE
#include "copa_keygen.h"
#include "copa_sampler.h"
#include "copa_irs64.h"
#include "copa_scalar64_tables.h"
#include "ntt_crt.h"
#include "poly_inv.h"
#include "quant.h"
#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <x86intrin.h>
static uint64_t state=UINT64_C(0x123456789abcdef1);
static uint64_t rng(void){state^=state<<13;state^=state>>7;state^=state<<17;return state;}
static uint64_t ticks(void){unsigned aux;_mm_lfence();uint64_t t=__rdtscp(&aux);_mm_lfence();return t;}
static volatile uint64_t sink;
typedef struct{long double mean,m2;uint64_t n;} moments;
static void record(moments *s,uint64_t x){++s->n;long double delta=x-s->mean;s->mean+=delta/s->n;s->m2+=delta*(x-s->mean);}
static int cmp(const void *a,const void *b){uint64_t x=*(const uint64_t *)a,y=*(const uint64_t *)b;return (x>y)-(x<y);}
typedef struct {
 uint8_t seed[32];int32_t x[COPA_N],g[COPA_N],out[(1+COPA_ELL)*COPA_N];
 int64_t product[COPA_N],dot;uint64_t random[6],nums[4];copa_irs64_context irs;
 unsigned minus,status;copa_fp63 exp[4];
} ct_input;
#define NOINLINE __attribute__((noinline))
static NOINLINE void run_quant(ct_input *v){copa_quant_block(v->out,v->x,COPA_MODE);}
static NOINLINE void run_noise(ct_input *v){copa_expand_s((int32_t (*)[COPA_N])v->out,v->seed);}
static NOINLINE void run_inverse(ct_input *v){v->status=(unsigned)copa_poly_inv(v->out,v->x,COPA_Q);}
static NOINLINE void run_irs(ct_input *v){v->status=copa_irs64_step(&v->minus,v->dot,&v->irs,v->random);}
static NOINLINE void run_ntt(ct_input *v){copa_ntt_crt_mul(v->product,v->x,v->g,COPA_N);}
static NOINLINE void run_exp(ct_input *v){copa64_exp_ratio4(v->exp,v->nums,&copa64_r2);}
int main(int argc,char **argv)
{
 const char *component=argc>1?argv[1]:"quant";unsigned samples=argc>2?(unsigned)strtoul(argv[2],NULL,10):20000;
 void (*worker)(ct_input *)=NULL;
 if(!strcmp(component,"quant"))worker=run_quant;else if(!strcmp(component,"noise"))worker=run_noise;
 else if(!strcmp(component,"inverse"))worker=run_inverse;else if(!strcmp(component,"irs"))worker=run_irs;
 else if(!strcmp(component,"ntt"))worker=run_ntt;else if(!strcmp(component,"exp"))worker=run_exp;else return 2;
 moments m[2]={{0}};ct_input input;memset(&input,0,sizeof(input));
 copa_irs64_context contexts[2];memset(contexts,0,sizeof(contexts));copa_irs64_prepare(&contexts[0],1);copa_irs64_prepare(&contexts[1],64);
 uint64_t *times=calloc(samples,sizeof(*times)),*sorted=calloc(samples,sizeof(*sorted));unsigned char *classes=calloc(samples,1);assert(times&&sorted&&classes);
 for(unsigned trial=0;trial<samples+200;++trial){
  unsigned cls=(unsigned)(rng()&1);uint32_t mask=0u-cls;
  for(unsigned j=0;j<32;++j)input.seed[j]=(uint8_t)rng()&(uint8_t)mask;
  for(unsigned j=0;j<COPA_N;++j){int32_t x=(int32_t)(rng()%7)-3,g=(int32_t)(rng()%11)-5;input.x[j]=(int32_t)(((uint32_t)x&mask)|((uint32_t)(j==0)&~mask));input.g[j]=(int32_t)(((uint32_t)g&mask)|((uint32_t)(j==0)&~mask));}
  if(worker==run_inverse){
   /* Compare units with the same public success status. For odd j,
    * 1+2*X^j is a unit modulo every supported q: a scalar root over
    * each odd prime factor (<512) cannot have order 512. */
   unsigned index=1+2*(rng()%128);memset(input.x,0,sizeof(input.x));input.x[0]=1;input.x[index]=(int32_t)(2*cls);
  }
  /* Prepare both contexts once, and copy fixed-size inputs outside timing.
   * This avoids class-dependent coefficient generation before each sample. */
  memcpy(&input.irs,&contexts[cls],sizeof(input.irs));input.dot=-(int64_t)((1+63*cls)*cls);
  for(unsigned j=0;j<6;++j)input.random[j]=rng();
  for(unsigned j=0;j<4;++j)input.nums[j]=10000*(1+cls*(100+j));
  __asm__ volatile(""::"m"(input):"memory");
  uint64_t start=ticks();worker(&input);uint64_t elapsed=ticks()-start;
  sink+=(uint32_t)input.out[0]+input.status+input.minus;
  if(worker==run_inverse)assert(input.status==1);
  if(trial>=200){unsigned j=trial-200;times[j]=sorted[j]=elapsed;classes[j]=(unsigned char)cls;record(&m[cls],elapsed);}
 }
    long double variance=m[0].m2/(m[0].n-1)/m[0].n+m[1].m2/(m[1].n-1)/m[1].n;
    long double t=(m[0].mean-m[1].mean)/sqrtl(variance);
    qsort(sorted,samples,sizeof(*sorted),cmp);
    long double max_t=fabsl(t);
    for(unsigned percent=90;percent<=99;percent+=3){
        uint64_t cutoff=sorted[(size_t)samples*percent/100];moments crop[2]={{0}};
        for(unsigned i=0;i<samples;++i)if(times[i]<=cutoff)record(&crop[classes[i]],times[i]);
        long double v=crop[0].m2/(crop[0].n-1)/crop[0].n+crop[1].m2/(crop[1].n-1)/crop[1].n;
        long double tc=(crop[0].mean-crop[1].mean)/sqrtl(v);if(fabsl(tc)>max_t)max_t=fabsl(tc);
    }
    free(times);free(sorted);free(classes);
    printf("{\"mode\":%d,\"component\":\"%s\",\"samples\":%u,\"n0\":%"PRIu64",\"n1\":%"PRIu64",\"mean0\":%.3Lf,\"mean1\":%.3Lf,\"welch_t\":%.6Lf,\"max_abs_t_raw_and_90_93_96_99_percent_crops\":%.6Lf,\"threshold\":4.5,\"result\":\"%s\"}\n",COPA_MODE,component,samples,m[0].n,m[1].n,m[0].mean,m[1].mean,t,max_t,max_t<4.5?"PASS":"DETECTED");
    return max_t<4.5?0:1;
}
