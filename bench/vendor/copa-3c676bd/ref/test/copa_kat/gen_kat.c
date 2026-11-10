#include "nist_rng.h"
#include "copa_sign.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>
static void hex(const char *label,const uint8_t *p,size_t n){printf("%s = ",label);for(size_t i=0;i<n;++i)printf("%02X",p[i]);putchar('\n');}
int main(void)
{
    copa_nist_selftest();uint8_t entropy[48],seeds[10][48],messages[10][330];
    for(unsigned j=0;j<48;++j)entropy[j]=(uint8_t)j;
    copa_nist_seed(entropy);
    for(unsigned j=0;j<10;++j){randombytes(seeds[j],48);randombytes(messages[j],33*(j+1));}
    printf("# Copa mode %d, AES-256 CTR_DRBG, seed/master bytes 00..2F\n",COPA_MODE);
    for(unsigned j=0;j<10;++j){
        copa_secret_key sk;uint8_t pk[COPA_PK_BYTES],sig[COPA_SIG_BYTES];copa_nist_seed(seeds[j]);
        copa_keygen_stats stats;assert(!copa_keygen(pk,&sk,&stats));
        assert(stats.partial_lower_norm_reject==0);assert(!copa_sign(sig,messages[j],33*(j+1),&sk,NULL));
        assert(!copa_verify(sig,sizeof(sig),messages[j],33*(j+1),pk,sizeof(pk)));
        printf("\ncount = %u\n",j);hex("seed",seeds[j],48);printf("mlen = %u\n",33*(j+1));hex("msg",messages[j],33*(j+1));hex("pk",pk,sizeof(pk));
        /* Portable secret-key serialization: pk || K || H(pk) || signed16le(s). */
        uint8_t serialized[COPA_PK_BYTES+32+64+2*COPA_K*COPA_N];size_t pos=0;
        memcpy(serialized+pos,sk.pk,COPA_PK_BYTES);pos+=COPA_PK_BYTES;memcpy(serialized+pos,sk.key,32);pos+=32;memcpy(serialized+pos,sk.hpk,64);pos+=64;
        for(unsigned r=0;r<COPA_K;++r)for(unsigned i=0;i<COPA_N;++i){uint16_t x=(uint16_t)sk.s[r][i];serialized[pos++]=(uint8_t)x;serialized[pos++]=(uint8_t)(x>>8);}
        hex("sk",serialized,pos);printf("siglen = %u\n",COPA_SIG_BYTES);hex("sig",sig,sizeof(sig));
    }
    return 0;
}
