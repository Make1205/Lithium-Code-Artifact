#include "copa_api.h"
#ifdef COPA_SIMD_BITS
#define COPA_EXP_BACKEND COPA_SIMD_BITS
#else
#define COPA_EXP_BACKEND 0
#endif
#define EXP_SK_VALID (UINT64_C(0x4350000000000000)^((uint64_t)COPA_MODE<<40)^((uint64_t)COPA_EXP_BACKEND<<24)^sizeof(copa_secret_key))
#ifdef COPA_VERIFY31
#define EXP_PK_CACHE_KIND 31
#else
#define EXP_PK_CACHE_KIND 23
#endif
#define EXP_PK_VALID (UINT64_C(0x4351000000000000)^((uint64_t)COPA_MODE<<40)^((uint64_t)COPA_EXP_BACKEND<<24)^((uint64_t)EXP_PK_CACHE_KIND<<34)^sizeof(copa_public_key))
/* Volatile writes clear every cache byte and padding, including the marker.
 * Expanded objects are implementation-local, never serialized wire keys. */
static void copa_expanded_zero(void *p,size_t n){volatile uint8_t *v=(volatile uint8_t *)p;while(n--)*v++=0;}
void crypto_sign_clear_expanded_sk(copa_expanded_sk *expanded){if(expanded)copa_expanded_zero(expanded,sizeof(*expanded));}
void crypto_sign_clear_expanded_pk(copa_expanded_pk *expanded){if(expanded)copa_expanded_zero(expanded,sizeof(*expanded));}
int crypto_sign_expand_sk(copa_expanded_sk *expanded,const uint8_t sk[COPA_SK_BYTES]){
 if(!expanded)return -1;
 crypto_sign_clear_expanded_sk(expanded);
 if(!sk||copa_unpack_sk(&expanded->key,sk,COPA_SK_BYTES)){crypto_sign_clear_expanded_sk(expanded);return -1;}
 expanded->initialized=EXP_SK_VALID;return 0;
}
int crypto_sign_expand_pk(copa_expanded_pk *expanded,const uint8_t pk[COPA_PK_BYTES]){
 if(!expanded)return -1;
 crypto_sign_clear_expanded_pk(expanded);
 if(!pk||copa_prepare_public(&expanded->key,pk,COPA_PK_BYTES)){crypto_sign_clear_expanded_pk(expanded);return -1;}
 expanded->initialized=EXP_PK_VALID;return 0;
}
int crypto_sign_signature_expanded(uint8_t sig[COPA_SIG_BYTES],size_t *siglen,const uint8_t *m,size_t mlen,const copa_expanded_sk *expanded){
 if(!siglen)return -1;
 *siglen=0;
 if(!sig||!expanded||expanded->initialized!=EXP_SK_VALID||(!m&&mlen))return -1;
 int result=copa_sign(sig,m,mlen,&expanded->key,0);if(!result)*siglen=COPA_SIG_BYTES;return result;
}
int crypto_sign_verify_expanded(const uint8_t *sig,size_t siglen,const uint8_t *m,size_t mlen,const copa_expanded_pk *expanded){
 if(!sig||!expanded||expanded->initialized!=EXP_PK_VALID||(!m&&mlen))return -1;
 return copa_verify_prepared(sig,siglen,m,mlen,&expanded->key);
}
