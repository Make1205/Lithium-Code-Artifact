/* Per-call stage medians. Nested buckets are explicitly not additive. */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <x86intrin.h>
#ifdef BENCH_LITHIUM
#include "api.h"
#include "prof.h"
#define NB PT_NBUCKETS
static const char *names[]={"PT_SETUP","PT_EXPAND_A0","PT_SAMPLE_Y","PT_G_SHAKE","PT_G_BASESAMP","PT_G_APPROXEXP","PT_G_FINAL","PT_COMMIT","PT_NTT_FWD","PT_NTT_PW","PT_NTT_INV","PT_HIGHBITS","PT_CHALLENGE","PT_SIGN","PT_EXP","PT_NORMCHECK","PT_MAKEHINT","PT_PACK","PT_KG_EXPAND_A0","PT_KG_NOISE","PT_KG_NORM","PT_KG_COMPUTEB","PT_KG_PACK","PT_VF_UNPACK","PT_VF_SETUP","PT_VF_A0","PT_VF_MATMUL","PT_VF_HINT","PT_VF_CHALLENGE","PT_VF_NORM"};
#ifdef PROFILE_ENABLED
#define COUNTERS prof_cyc
#define RESET() prof_reset()
#else
static uint64_t empty[NB];
#define COUNTERS empty
#define RESET() memset(empty,0,sizeof(empty))
#endif
#else
#include "copa_api.h"
#include "copa_prof.h"
#define NB CP_COUNT
uint64_t copa_profile_cycles[CP_COUNT],copa_profile_calls[CP_COUNT];
#define COUNTERS copa_profile_cycles
#define RESET() do{memset(COUNTERS,0,sizeof(COUNTERS));memset(copa_profile_calls,0,sizeof(copa_profile_calls));}while(0)
static const char *names[]={"CP_PARSE_SK","CP_Y","CP_AY","CP_HIGH_CH","CP_IRS","CP_BS","CP_HINT_BV","CP_RANS","CP_PK_UNLAB","CP_A","CP_D","CP_MATRIX_NTT","CP_DECODE","CP_VF_NTT","CP_CHALLENGE","CP_RECONSTRUCT","CP_S","CP_INVERSE","CP_QUANT","CP_PACK","CP_NTT_FORWARD","CP_NTT_POINTWISE","CP_NTT_INVERSE","CP_SK_DECODE","CP_SK_HASH","CP_SK_NTT","CP_IRS_PREP","CP_SIGN_SETUP","CP_B_REBUILD","CP_Y_SHAKE","CP_IRS_PROJECT","CP_IRS_TRANSITION","CP_IRS_ROTATE","CP_UNLAB","CP_B_NTT","CP_A_NTT","CP_VF_NORM","CP_KG_HASH","CP_KG_NORM","CP_KG_B","CP_KG_MATRIX"};
#endif
static uint64_t tick(void){unsigned a;_mm_lfence();uint64_t t=__rdtscp(&a);_mm_lfence();return t;}
static int cmp(const void *a,const void *b){uint64_t x=*(const uint64_t*)a,y=*(const uint64_t*)b;return(x>y)-(x<y);}
int main(int argc,char **argv){assert(argc==3);unsigned n=(unsigned)strtoul(argv[2],0,10);assert(n>=101);const char *op=argv[1];
 uint64_t *v=calloc((size_t)(NB+2)*n,sizeof(*v));assert(v);uint8_t m[59]={0};size_t len=0;
 uint64_t attempts=0,even=0,odd=0,norm=0,other=0,partial=0;
#ifdef BENCH_LITHIUM
 uint8_t pk[CRYPTO_PUBLICKEYBYTES],sk[CRYPTO_SECRETKEYBYTES],sig[CRYPTO_BYTES];
 assert(!crypto_sign_keypair(pk,sk));assert(!crypto_sign_signature(sig,&len,m,sizeof(m),sk));
#else
 uint8_t pk[COPA_PK_BYTES],sk[COPA_SK_BYTES],sig[COPA_SIG_BYTES];
 assert(!copa_keypair_bytes(pk,sk,0));assert(!copa_sign_signature(sig,&len,m,sizeof(m),sk,sizeof(sk),0));
#endif
 for(unsigned sample=0;sample<n+20;++sample){RESET();uint64_t start=tick();
#ifdef BENCH_LITHIUM
  if(!strcmp(op,"keygen")){assert(!crypto_sign_keypair(pk,sk));}
  else if(!strcmp(op,"sign")){assert(!crypto_sign_signature(sig,&len,m,sizeof(m),sk));}
  else{assert(!strcmp(op,"verify"));assert(!crypto_sign_verify(sig,len,m,sizeof(m),pk));}
#else
  copa_keygen_stats st={0};
  if(!strcmp(op,"keygen")){assert(!copa_keypair_bytes(pk,sk,&st));}
  else if(!strcmp(op,"sign")){assert(!copa_sign_signature(sig,&len,m,sizeof(m),sk,sizeof(sk),0));}
  else{assert(!strcmp(op,"verify"));assert(!copa_verify(sig,len,m,sizeof(m),pk,sizeof(pk)));}
#endif
  uint64_t elapsed=tick()-start;
  uint64_t accounted=0;
#ifdef BENCH_LITHIUM
 const unsigned kg[]={PT_KG_EXPAND_A0,PT_KG_NOISE,PT_KG_NORM,PT_KG_COMPUTEB,PT_KG_PACK};
 const unsigned sg[]={PT_SETUP,PT_EXPAND_A0,PT_SAMPLE_Y,PT_COMMIT,PT_HIGHBITS,PT_CHALLENGE,PT_SIGN,PT_NORMCHECK,PT_MAKEHINT,PT_PACK};
 const unsigned vf[]={PT_VF_UNPACK,PT_VF_SETUP,PT_VF_A0,PT_VF_MATMUL,PT_VF_HINT,PT_VF_CHALLENGE,PT_VF_NORM};
#else
 const unsigned kg[]={CP_KG_HASH,CP_S,CP_INVERSE,CP_A,CP_D,CP_QUANT,CP_KG_NORM,CP_KG_B,CP_PACK,CP_KG_MATRIX};
 const unsigned sg[]={CP_PARSE_SK,CP_SIGN_SETUP,CP_Y,CP_AY,CP_HIGH_CH,CP_IRS,CP_BS,CP_HINT_BV,CP_RANS};
 const unsigned vf[]={CP_PK_UNLAB,CP_A,CP_D,CP_MATRIX_NTT,CP_DECODE,CP_VF_NTT,CP_CHALLENGE,CP_RECONSTRUCT,CP_VF_NORM};
#endif
 const unsigned *ids=!strcmp(op,"keygen")?kg:!strcmp(op,"sign")?sg:vf;
 unsigned groups=!strcmp(op,"keygen")?sizeof(kg)/sizeof(*kg):!strcmp(op,"sign")?sizeof(sg)/sizeof(*sg):sizeof(vf)/sizeof(*vf);
 for(unsigned j=0;j<groups;++j)accounted+=COUNTERS[ids[j]];
 assert(accounted<=elapsed);
 if(sample>=20){unsigned i=sample-20;v[(NB+1)*n+i]=elapsed-accounted;v[NB*n+i]=elapsed;for(unsigned j=0;j<NB;++j)v[j*n+i]=COUNTERS[j];
#ifndef BENCH_LITHIUM
   attempts+=st.attempts;even+=st.noninvertible_even;odd+=st.noninvertible_odd;norm+=st.norm_reject;other+=st.lower_norm_reject+st.packing_failure;partial+=st.partial_norm_reject;
#endif
  }
 }
 for(unsigned j=0;j<=NB+1;++j)qsort(v+j*n,n,sizeof(*v),cmp);
#ifdef BENCH_LITHIUM
 printf("{\"scheme\":\"LITHIUM-%d\",",LITHIUM_MODE);
#else
 printf("{\"scheme\":\"Copa-%s\",",COPA_MODE==2?"2":COPA_MODE==3?"3":COPA_MODE==5?"5a":"5b");
#endif
 printf("\"partial_upper_norm\":%llu,",(unsigned long long)partial);
 printf("\"operation\":\"%s\",\"samples\":%u,\"median_cycles\":%llu,\"keygen_attempts\":%llu,\"nonunit_even\":%llu,\"nonunit_odd\":%llu,\"norm_reject\":%llu,\"other\":%llu,\"stages\":{",op,n,(unsigned long long)v[NB*n+n/2],(unsigned long long)attempts,(unsigned long long)even,(unsigned long long)odd,(unsigned long long)norm,(unsigned long long)other);
 for(unsigned j=0;j<NB;++j){printf("%s\"%s\":%llu",j?",":"",names[j],(unsigned long long)v[j*n+n/2]);}printf(",\"OTHER\":%llu}}\n",(unsigned long long)v[(NB+1)*n+n/2]);free(v);return 0;
}
