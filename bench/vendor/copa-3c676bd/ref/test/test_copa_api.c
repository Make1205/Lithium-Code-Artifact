#include "copa_api.h"
#include "copa_kat/nist_rng.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>
int main(void){uint8_t pk[COPA_PK_BYTES],bytes[COPA_SK_BYTES],other[COPA_SK_BYTES],a[COPA_SIG_BYTES],b[COPA_SIG_BYTES],seed[48]={0};copa_secret_key sk,decoded;copa_public_key public_key;size_t len;
 copa_nist_seed(seed);assert(!copa_keygen(pk,&sk,0));assert(!copa_pack_sk(bytes,&sk));assert(!copa_unpack_sk(&decoded,bytes,sizeof(bytes)));assert(!copa_pack_sk(other,&decoded));assert(!memcmp(bytes,other,sizeof(bytes)));assert(!copa_prepare_public(&public_key,pk,sizeof(pk)));
 for(unsigned j=0;j<20;++j){seed[0]=(uint8_t)j;copa_nist_seed(seed);assert(!copa_sign(a,seed,48,&sk,0));copa_nist_seed(seed);assert(!copa_sign_signature(b,&len,seed,48,bytes,sizeof(bytes),0));assert(len==sizeof(b)&&!memcmp(a,b,sizeof(a)));assert(!copa_verify(b,len,seed,48,pk,sizeof(pk)));assert(!copa_verify_prepared(b,len,seed,48,&public_key));}
 assert(copa_unpack_sk(&decoded,bytes,sizeof(bytes)-1));bytes[COPA_PK_BYTES+32]^=1;assert(copa_unpack_sk(&decoded,bytes,sizeof(bytes)));printf("PASS mode=%d canonical byte API and cached/byte-identical signatures\n",COPA_MODE);return 0;}
