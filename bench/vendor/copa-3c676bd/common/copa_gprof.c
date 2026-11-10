#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef BENCH_LITHIUM
#include "api.h"
#else
#include "copa_api.h"
#endif
extern void moncontrol(int);
int main(int argc,char **argv){assert(argc==3);moncontrol(0);unsigned n=(unsigned)strtoul(argv[2],0,10);uint8_t m[59]={0};size_t len;
#ifdef BENCH_LITHIUM
 uint8_t pk[CRYPTO_PUBLICKEYBYTES],sk[CRYPTO_SECRETKEYBYTES],sig[CRYPTO_BYTES];assert(!crypto_sign_keypair(pk,sk));assert(!crypto_sign_signature(sig,&len,m,sizeof(m),sk));
#else
 uint8_t pk[COPA_PK_BYTES],sk[COPA_SK_BYTES],sig[COPA_SIG_BYTES];assert(!copa_keypair_bytes(pk,sk,0));assert(!copa_sign_signature(sig,&len,m,sizeof(m),sk,sizeof(sk),0));
#endif
 moncontrol(1);for(unsigned i=0;i<n;++i){
#ifdef BENCH_LITHIUM
 if(!strcmp(argv[1],"keygen"))assert(!crypto_sign_keypair(pk,sk));else if(!strcmp(argv[1],"sign"))assert(!crypto_sign_signature(sig,&len,m,sizeof(m),sk));else assert(!crypto_sign_verify(sig,len,m,sizeof(m),pk));
#else
 if(!strcmp(argv[1],"keygen"))assert(!copa_keypair_bytes(pk,sk,0));else if(!strcmp(argv[1],"sign"))assert(!copa_sign_signature(sig,&len,m,sizeof(m),sk,sizeof(sk),0));else assert(!copa_verify(sig,len,m,sizeof(m),pk,sizeof(pk)));
#endif
 }moncontrol(0);printf("PASS %s samples=%u\n",argv[1],n);return 0;}
