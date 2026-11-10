#ifndef BENCH_ADAPTER_H
#define BENCH_ADAPTER_H
#include "api.h"
#include "irs_rans.h"
#define B_N N
#define B_Q Q
#define B_PK CRYPTO_PUBLICKEYBYTES
#define B_SK CRYPTO_SECRETKEYBYTES
#define B_SIG CRYPTO_BYTES
#define b_keygen crypto_sign_keypair
#define b_sign crypto_sign_signature
#define b_verify crypto_sign_verify
static size_t b_unpadded(const uint8_t *p,size_t n){size_t r=p[CHALLENGE_PACKEDBYTES]|((size_t)p[CHALLENGE_PACKEDBYTES+1]<<8);return n-RANS_RESERVED_BYTES+r;}
#endif
