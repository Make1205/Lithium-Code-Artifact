#include "nist_rng.h"
#ifdef COPA_KAT_OPENSSL
#include <openssl/evp.h>
#else
#include "nist_aes.h"
#endif
#include <stdlib.h>
#include <string.h>
#include <assert.h>
static uint8_t key[32],v[16];
static void aes(const uint8_t k[32],const uint8_t in[16],uint8_t out[16])
{
#ifdef COPA_KAT_OPENSSL
    EVP_CIPHER_CTX *ctx=EVP_CIPHER_CTX_new();int n=0,final=0;
    if(!ctx||!EVP_EncryptInit_ex(ctx,EVP_aes_256_ecb(),NULL,k,NULL)||!EVP_CIPHER_CTX_set_padding(ctx,0)||!EVP_EncryptUpdate(ctx,out,&n,in,16)||!EVP_EncryptFinal_ex(ctx,out+n,&final)||n+final!=16)abort();
    EVP_CIPHER_CTX_free(ctx);
#else
    copa_test_aes256(k,in,out);
#endif
}
static void increment(void){for(int j=15;j>=0;--j)if(++v[j])break;}
static void update(const uint8_t *provided)
{
    uint8_t tmp[48];for(unsigned j=0;j<3;++j){increment();aes(key,v,tmp+16*j);}
    if(provided)for(unsigned j=0;j<48;++j)tmp[j]^=provided[j];
    memcpy(key,tmp,32);memcpy(v,tmp+32,16);
}
void copa_nist_seed(const uint8_t entropy[48]){memset(key,0,32);memset(v,0,16);update(entropy);}
void randombytes(uint8_t *out,size_t length)
{
    while(length){uint8_t block[16];increment();aes(key,v,block);size_t n=length<16?length:16;memcpy(out,block,n);out+=n;length-=n;}
    update(NULL);
}
void copa_nist_selftest(void)
{
    static const uint8_t expected[32]={0x06,0x15,0x50,0x23,0x4d,0x15,0x8c,0x5e,0xc9,0x55,0x95,0xfe,0x04,0xef,0x7a,0x25,0x76,0x7f,0x2e,0x24,0xcc,0x2b,0xc4,0x79,0xd0,0x9d,0x86,0xdc,0x9a,0xbc,0xfd,0xe7};
    uint8_t entropy[48],out[32];for(unsigned j=0;j<48;++j)entropy[j]=(uint8_t)j;
    copa_nist_seed(entropy);randombytes(out,sizeof(out));assert(!memcmp(out,expected,32));
}
