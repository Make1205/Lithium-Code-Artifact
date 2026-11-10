#ifndef COPA_API_H
#define COPA_API_H
#include "copa_sign.h"
#define COPA_SK_BYTES (COPA_PK_BYTES+32+64+2*COPA_K*COPA_N)
/* Existing canonical KAT format: pk || key || H(pk) || signed16le(s). */
int copa_pack_sk(uint8_t out[COPA_SK_BYTES],const copa_secret_key *key);
int copa_unpack_sk(copa_secret_key *key,const uint8_t *in,size_t length);
int copa_keypair_bytes(uint8_t pk[COPA_PK_BYTES],uint8_t sk[COPA_SK_BYTES],copa_keygen_stats *stats);
int copa_sign_signature(uint8_t sig[COPA_SIG_BYTES],size_t *siglen,const uint8_t *m,size_t mlen,const uint8_t *sk,size_t sklen,copa_sign_stats *stats);
typedef struct {copa_public_ntt matrix;uint8_t hpk[64];} copa_public_key;
int copa_prepare_public(copa_public_key *key,const uint8_t *pk,size_t length);
int copa_verify_prepared(const uint8_t *sig,size_t siglen,const uint8_t *m,size_t mlen,const copa_public_key *key);
/* Runtime objects are backend/version specific and must not overlap input.
 * They contain caches, not a new serialized key format. Clear after use. */
typedef struct {uint64_t initialized;copa_secret_key key;} copa_expanded_sk;
typedef struct {uint64_t initialized;copa_public_key key;} copa_expanded_pk;
int crypto_sign_expand_sk(copa_expanded_sk *expanded,const uint8_t sk[COPA_SK_BYTES]);
int crypto_sign_signature_expanded(uint8_t sig[COPA_SIG_BYTES],size_t *siglen,const uint8_t *m,size_t mlen,const copa_expanded_sk *expanded);
int crypto_sign_expand_pk(copa_expanded_pk *expanded,const uint8_t pk[COPA_PK_BYTES]);
int crypto_sign_verify_expanded(const uint8_t *sig,size_t siglen,const uint8_t *m,size_t mlen,const copa_expanded_pk *expanded);
void crypto_sign_clear_expanded_sk(copa_expanded_sk *expanded);
void crypto_sign_clear_expanded_pk(copa_expanded_pk *expanded);
#endif
