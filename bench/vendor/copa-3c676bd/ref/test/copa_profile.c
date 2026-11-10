/* gprof entry point: one operation per process, no timing claims from -pg. */
#include "copa_sign.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc,char **argv){
 assert(argc==3);unsigned n=(unsigned)strtoul(argv[2],0,10);assert(n);
 uint8_t pk[COPA_PK_BYTES],sig[COPA_SIG_BYTES],msg[59]={0};copa_secret_key sk;
 assert(!copa_keygen(pk,&sk,0));assert(!copa_sign(sig,msg,sizeof(msg),&sk,0));
 for(unsigned i=0;i<n;++i){
  if(!strcmp(argv[1],"keygen")){assert(!copa_keygen(pk,&sk,0));}
  else if(!strcmp(argv[1],"sign")){assert(!copa_sign(sig,msg,sizeof(msg),&sk,0));}
  else {assert(!strcmp(argv[1],"verify"));assert(!copa_verify(sig,sizeof(sig),msg,sizeof(msg),pk,sizeof(pk)));}
 }
 printf("mode=%d op=%s samples=%u\n",COPA_MODE,argv[1],n);return 0;
}
