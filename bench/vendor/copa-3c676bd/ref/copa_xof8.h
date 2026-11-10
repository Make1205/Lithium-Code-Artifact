#ifndef COPA_XOF8_H
#define COPA_XOF8_H
#include <stddef.h>
#include <stdint.h>
#define COPA_XOF8_RATE 136
/* Opaque storage: scalar eight-state, two Keccak4x, or one Keccak8x engine. */
typedef struct __attribute__((aligned(64))) {uint64_t opaque[208];} copa_xof8;
void copa_xof4_init(copa_xof8 *s,const uint8_t *inputs[4],size_t len);
void copa_xof4_block(copa_xof8 *s,uint8_t out[4][COPA_XOF8_RATE]);
/* Exact SHAKE256(seed || domain || lane || 0), eight fixed lanes. */
void copa_xof8_seed_init(copa_xof8 *,const uint8_t seed[32],uint8_t domain);
void copa_xof8_init(copa_xof8 *s,const uint8_t *inputs[8],size_t len);
void copa_xof8_block(copa_xof8 *s,uint8_t out[8][COPA_XOF8_RATE]);
static inline void copa_xof8_prefix(copa_xof8 *s,uint8_t *out[8],size_t len){
 uint8_t block[8][COPA_XOF8_RATE];size_t done=0;
 while(done<len){copa_xof8_block(s,block);size_t n=len-done;if(n>COPA_XOF8_RATE)n=COPA_XOF8_RATE;
  for(unsigned j=0;j<8;++j)for(size_t k=0;k<n;++k)out[j][done+k]=block[j][k];
  done+=n;
 }
}
#endif
