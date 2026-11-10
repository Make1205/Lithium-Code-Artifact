#ifndef COPA_SIGN_H
#define COPA_SIGN_H
#include "copa_keygen.h"
#include <stddef.h>
typedef struct {uint64_t attempts,irs_accepted,bs_restarts,bv_restarts,encode_overflows,payload_used,irs_steps;} copa_sign_stats;
int copa_sign(uint8_t sig[COPA_SIG_BYTES],const uint8_t *message,size_t length,const copa_secret_key *sk,copa_sign_stats *stats);
int copa_verify(const uint8_t *sig,size_t siglen,const uint8_t *message,size_t length,const uint8_t *pk,size_t pklen);
#endif
