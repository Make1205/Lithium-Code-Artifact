#include "copa_api.h"
#include "copa_kat/nist_rng.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static void zero(const void *p,size_t n){const uint8_t *b=p;for(size_t j=0;j<n;++j)assert(b[j]==0);}
int main(int argc,char **argv){
 unsigned samples=argc>1?(unsigned)strtoul(argv[1],0,10):100;assert(samples&&samples<=10000);
 copa_nist_selftest();
 uint8_t pk[COPA_PK_BYTES],sk[COPA_SK_BYTES],a[COPA_SIG_BYTES],b[COPA_SIG_BYTES],seed[48],message[65];
 copa_expanded_sk es;copa_expanded_pk ep;
 for(unsigned key=0;key<10;++key){for(unsigned j=0;j<48;++j)seed[j]=(uint8_t)(j+key*19);copa_nist_seed(seed);assert(!copa_keypair_bytes(pk,sk,0));
  assert(!crypto_sign_expand_sk(&es,sk));assert(!crypto_sign_expand_pk(&ep,pk));
  for(unsigned sample=0;sample<samples;++sample){for(unsigned j=0;j<65;++j)message[j]=(uint8_t)(key+sample+j);for(unsigned j=0;j<48;++j)seed[j]=(uint8_t)(sample*17+j+key);size_t na=0,nb=0,mlen=sample%66;
   copa_nist_seed(seed);assert(!copa_sign_signature(a,&na,message,mlen,sk,sizeof(sk),0));
   copa_nist_seed(seed);assert(!crypto_sign_signature_expanded(b,&nb,message,mlen,&es));
   assert(na==COPA_SIG_BYTES&&nb==na&&!memcmp(a,b,na));assert(!copa_verify(a,na,message,mlen,pk,sizeof(pk)));assert(!crypto_sign_verify_expanded(b,nb,message,mlen,&ep));
   b[0]^=1;assert(crypto_sign_verify_expanded(b,nb,message,mlen,&ep));assert(copa_verify(b,nb,message,mlen,pk,sizeof(pk)));
  }
  crypto_sign_clear_expanded_sk(&es);crypto_sign_clear_expanded_pk(&ep);zero(&es,sizeof(es));zero(&ep,sizeof(ep));
 }
 memset(sk+COPA_PK_BYTES+96,0,sizeof(sk)-COPA_PK_BYTES-96);memset(&es,0xa5,sizeof(es));assert(crypto_sign_expand_sk(&es,sk));zero(&es,sizeof(es));
 memset(pk+32,255,sizeof(pk)-32);memset(&ep,0xa5,sizeof(ep));assert(crypto_sign_expand_pk(&ep,pk));zero(&ep,sizeof(ep));
 printf("PASS mode=%d expanded API: 10 keys, %u byte-identical signatures, mutual verification, malformed keys and complete wipes\n",COPA_MODE,10*samples);return 0;
}
