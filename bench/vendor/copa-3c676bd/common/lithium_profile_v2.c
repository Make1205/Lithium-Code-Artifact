#include "api.h"
#include "prof.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc,char **argv){assert(argc==3);unsigned n=(unsigned)strtoul(argv[2],0,10);uint8_t pk[CRYPTO_PUBLICKEYBYTES],sk[CRYPTO_SECRETKEYBYTES],sig[CRYPTO_BYTES],m[59]={0};size_t len;assert(!crypto_sign_keypair(pk,sk));assert(!crypto_sign_signature(sig,&len,m,sizeof(m),sk));prof_reset();for(unsigned i=0;i<n;++i){if(!strcmp(argv[1],"keygen"))assert(!crypto_sign_keypair(pk,sk));else if(!strcmp(argv[1],"sign"))assert(!crypto_sign_signature(sig,&len,m,sizeof(m),sk));else assert(!crypto_sign_verify(sig,len,m,sizeof(m),pk));}printf("mode=%d op=%s samples=%u\n",LITHIUM_MODE,argv[1],n);prof_report("default AVX2",n);return 0;}
