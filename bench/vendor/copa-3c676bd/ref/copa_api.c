#include "copa_prof.h"
#include "copa_api.h"
#include "fips202.h"
#include "ntt_crt.h"
#include <string.h>
#ifdef COPA_VEC_POLY_HEADER
#include COPA_VEC_POLY_HEADER
#endif
static void copa_hash_public(uint8_t out[64],const uint8_t pk[COPA_PK_BYTES]){
 keccak_state st;uint8_t d=2;shake256_init(&st);shake256_absorb(&st,pk,COPA_PK_BYTES);shake256_absorb(&st,&d,1);shake256_finalize(&st);shake256_squeeze(out,64,&st);
}
int copa_prepare_public(copa_public_key *key,const uint8_t *pk,size_t length){
 if(length!=COPA_PK_BYTES)return -1;
#ifdef COPA_VERIFY31
 if(!copa_verify31_expand(key->matrix,pk))return -1;
#else
 if(!copa_matrix_expand(key->matrix,pk))return -1;
#endif
 copa_hash_public(key->hpk,pk);return 0;
}
int copa_pack_sk(uint8_t out[COPA_SK_BYTES],const copa_secret_key *key){
 size_t p=0;memcpy(out,key->pk,COPA_PK_BYTES);p+=COPA_PK_BYTES;memcpy(out+p,key->key,32);p+=32;memcpy(out+p,key->hpk,64);p+=64;
 for(unsigned r=0;r<COPA_K;++r)for(unsigned i=0;i<COPA_N;++i){int32_t v=key->s[r][i];if(v < -32768||v>32767)return -1;out[p++]=(uint8_t)v;out[p++]=(uint8_t)((uint32_t)v>>8);}
 return 0;
}
int copa_unpack_sk(copa_secret_key *key,const uint8_t *in,size_t length){
 if(length!=COPA_SK_BYTES)return -1;
 CP_BEGIN(psd);size_t p=0;memcpy(key->pk,in,COPA_PK_BYTES);p+=COPA_PK_BYTES;memcpy(key->key,in+p,32);p+=32;memcpy(key->hpk,in+p,64);p+=64;
 uint64_t norm;
#ifdef COPA_VEC_POLY_HEADER
 norm=copa_backend_unpack_s16(&key->s[0][0],in+p,COPA_K*COPA_N);p+=2*COPA_K*COPA_N;
#else
 norm=0;for(unsigned r=0;r<COPA_K;++r)for(unsigned i=0;i<COPA_N;++i){uint32_t v=(uint32_t)in[p]|((uint32_t)in[p+1]<<8);p+=2;key->s[r][i]=(int32_t)(v^32768)-32768;norm+=(uint64_t)((int64_t)key->s[r][i]*key->s[r][i]);}
#endif
 if(norm<COPA_MIN_NORM||norm>COPA_BK_SQ){return -1;}
 CP_END(CP_SK_DECODE,psd);
 CP_BEGIN(psh);uint8_t hash[64];copa_hash_public(hash,key->pk);unsigned bad=0;for(unsigned i=0;i<64;++i)bad|=hash[i]^key->hpk[i];if(bad)return -1;CP_END(CP_SK_HASH,psh);
 if(!copa_matrix_expand(key->public_ntt,key->pk))return -1;
 CP_BEGIN(psn);for(unsigned r=0;r<COPA_K;++r)copa_ntt_forward(key->secret_ntt[r],key->s[r],COPA_N);
 uint32_t corr[2][COPA_N];copa_ntt_dot_reversed(corr,key->secret_ntt,key->secret_ntt,COPA_K,COPA_N);copa_ntt_inverse_crt(key->autocorrelation,corr,COPA_N);CP_END(CP_SK_NTT,psn);
 CP_BEGIN(pip);copa_irs64_prepare(&key->irs_context,norm);CP_END(CP_IRS_PREP,pip);return 0;
}
int copa_keypair_bytes(uint8_t pk[COPA_PK_BYTES],uint8_t sk[COPA_SK_BYTES],copa_keygen_stats *stats){copa_secret_key key;if(copa_keygen_raw(pk,&key,stats))return -1;CP_BEGIN(ppacksk);int packed=copa_pack_sk(sk,&key);CP_END(CP_PACK,ppacksk);return packed;}
int copa_sign_signature(uint8_t sig[COPA_SIG_BYTES],size_t *siglen,const uint8_t *m,size_t mlen,const uint8_t *sk,size_t sklen,copa_sign_stats *stats){copa_secret_key key;CP_BEGIN(pp);int parsed=copa_unpack_sk(&key,sk,sklen);CP_END(CP_PARSE_SK,pp);if(parsed)return -1;int result=copa_sign(sig,m,mlen,&key,stats);if(!result)*siglen=COPA_SIG_BYTES;return result;}
