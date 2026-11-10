#include "copa_sign.h"
#include "api.h"
#include "sign.h"
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <x86intrin.h>
static uint64_t ticks(void){unsigned aux;_mm_lfence();uint64_t t=__rdtscp(&aux);_mm_lfence();return t;}
static int cmp(const void *a,const void *b){uint64_t x=*(const uint64_t*)a,y=*(const uint64_t*)b;return (x>y)-(x<y);}
static uint64_t median(uint64_t *x){qsort(x,101,sizeof(*x),cmp);return x[50];}
int main(void)
{
    uint8_t pk[COPA_PK_BYTES],sig[COPA_SIG_BYTES],msg[59]={0};copa_secret_key sk;
    uint8_t lpk[CRYPTO_PUBLICKEYBYTES],lsk[CRYPTO_SECRETKEYBYTES],lsig[CRYPTO_BYTES];size_t llen=0;
    uint64_t ckg[101],cs[101],cv[101],lkg[101],ls[101],lv[101],start;
    for(unsigned j=0;j<121;++j){
        start=ticks();assert(!copa_keygen(pk,&sk,NULL));uint64_t a=ticks()-start;
        start=ticks();assert(!copa_sign(sig,msg,sizeof(msg),&sk,NULL));uint64_t b=ticks()-start;
        start=ticks();assert(!copa_verify(sig,sizeof(sig),msg,sizeof(msg),pk,sizeof(pk)));uint64_t c=ticks()-start;
        start=ticks();assert(!crypto_sign_keypair(lpk,lsk));uint64_t d=ticks()-start;
        start=ticks();assert(!crypto_sign_signature(lsig,&llen,msg,sizeof(msg),lsk));uint64_t e=ticks()-start;
        start=ticks();assert(!crypto_sign_verify(lsig,llen,msg,sizeof(msg),lpk));uint64_t f=ticks()-start;
        if(j>=20){unsigned k=j-20;ckg[k]=a;cs[k]=b;cv[k]=c;lkg[k]=d;ls[k]=e;lv[k]=f;}
    }
    printf("{\"mode\":%d,\"samples\":101,\"warmup\":20,\"message_bytes\":59,\"unit\":\"serialized TSC ticks\",\"Copa_keygen\":%"PRIu64",\"Copa_sign\":%"PRIu64",\"Copa_verify\":%"PRIu64",\"Lithium120_keygen\":%"PRIu64",\"Lithium120_sign\":%"PRIu64",\"Lithium120_verify\":%"PRIu64"}\n",COPA_MODE,median(ckg),median(cs),median(cv),median(lkg),median(ls),median(lv));
    return 0;
}
