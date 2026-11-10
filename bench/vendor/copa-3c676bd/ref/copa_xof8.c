#include "copa_xof8.h"
#include "fips202.h"
#include <string.h>
_Static_assert(sizeof(keccak_state)*8<=sizeof(copa_xof8),"eight scalar states fit opaque storage");
void copa_xof8_init(copa_xof8 *s,const uint8_t *inputs[8],size_t len){
 keccak_state states[8];for(unsigned j=0;j<8;++j)shake256_absorb_once(&states[j],inputs[j],len);memcpy(s,states,sizeof(states));
}
void copa_xof8_block(copa_xof8 *s,uint8_t out[8][COPA_XOF8_RATE]){
 keccak_state states[8];memcpy(states,s,sizeof(states));for(unsigned j=0;j<8;++j)shake256_squeezeblocks(out[j],1,&states[j]);memcpy(s,states,sizeof(states));
}

void copa_xof4_init(copa_xof8 *s,const uint8_t *inputs[4],size_t len){keccak_state states[4];for(unsigned j=0;j<4;++j)shake256_absorb_once(&states[j],inputs[j],len);memcpy(s,states,sizeof(states));}
void copa_xof4_block(copa_xof8 *s,uint8_t out[4][136]){keccak_state states[4];memcpy(states,s,sizeof(states));for(unsigned j=0;j<4;++j)shake256_squeezeblocks(out[j],1,&states[j]);memcpy(s,states,sizeof(states));}

void copa_xof8_seed_init(copa_xof8 *s,const uint8_t seed[32],uint8_t domain){
 uint8_t inputs[8][35];const uint8_t *ptr[8];
 for(unsigned lane=0;lane<8;++lane){memcpy(inputs[lane],seed,32);inputs[lane][32]=domain;inputs[lane][33]=(uint8_t)lane;inputs[lane][34]=0;ptr[lane]=inputs[lane];}
 copa_xof8_init(s,ptr,35);
}
