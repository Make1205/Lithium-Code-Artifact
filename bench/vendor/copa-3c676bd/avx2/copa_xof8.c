#include "copa_xof8.h"
#include "fips202x4.h"
_Static_assert(2*sizeof(keccakx4_state)<=sizeof(copa_xof8),"two Keccak4x states fit");
void copa_keccak4_absorb(void *,const uint8_t *,const uint8_t *,const uint8_t *,const uint8_t *,size_t);
void copa_keccak4_block(void *,uint8_t *,uint8_t *,uint8_t *,uint8_t *);
void copa_xof8_init(copa_xof8 *s,const uint8_t *in[8],size_t n){copa_keccak4_absorb(s,in[0],in[1],in[2],in[3],n);copa_keccak4_absorb((uint8_t *)s+sizeof(keccakx4_state),in[4],in[5],in[6],in[7],n);}
void copa_xof8_block(copa_xof8 *s,uint8_t out[8][136]){copa_keccak4_block(s,out[0],out[1],out[2],out[3]);copa_keccak4_block((uint8_t *)s+sizeof(keccakx4_state),out[4],out[5],out[6],out[7]);}

void copa_xof4_init(copa_xof8 *s,const uint8_t *in[4],size_t n){copa_keccak4_absorb(s,in[0],in[1],in[2],in[3],n);}
void copa_xof4_block(copa_xof8 *s,uint8_t out[4][136]){copa_keccak4_block(s,out[0],out[1],out[2],out[3]);}

/* Broadcast the shared seed directly into the padded sponge states.
 * Squeeze performs the first permutation, exactly as generic absorb_once. */
void copa_xof8_seed_init(copa_xof8 *s,const uint8_t seed[32],uint8_t domain){
 uint64_t words[4]={0};for(unsigned i=0;i<4;++i)for(unsigned j=0;j<8;++j)words[i]|=(uint64_t)seed[8*i+j]<<(8*j);
 for(unsigned group=0;group<2;++group){__m256i *v=(__m256i *)((uint8_t *)s+group*sizeof(keccakx4_state));
  for(unsigned j=0;j<25;++j)v[j]=_mm256_setzero_si256();
  for(unsigned j=0;j<4;++j)v[j]=_mm256_set1_epi64x((int64_t)words[j]);
  uint64_t base=(uint64_t)domain|UINT64_C(0x1f000000);
  v[4]=_mm256_setr_epi64x(base|((4*group+0)<<8),base|((4*group+1)<<8),base|((4*group+2)<<8),base|((4*group+3)<<8));v[16]=_mm256_set1_epi64x(INT64_MIN);
 }
}
