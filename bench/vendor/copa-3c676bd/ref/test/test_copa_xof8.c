#include "copa_xof8.h"
#include "fips202.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>
int main(void){uint8_t inputs[8][136],expected[8][816],block[8][136];const uint8_t *p[8];unsigned sizes[]={0,1,7,8,35,42,64,135};
 for(unsigned lane=0;lane<8;++lane){p[lane]=inputs[lane];for(unsigned j=0;j<136;++j)inputs[lane][j]=(uint8_t)(j*19+lane*37);}
 for(unsigned t=0;t<8;++t){copa_xof8 s;copa_xof8_init(&s,p,sizes[t]);for(unsigned lane=0;lane<8;++lane)shake256(expected[lane],816,p[lane],sizes[t]);
  for(unsigned b=0;b<6;++b){copa_xof8_block(&s,block);for(unsigned lane=0;lane<8;++lane)assert(!memcmp(block[lane],expected[lane]+136*b,136));}
 }
 uint8_t seed[32];for(unsigned j=0;j<32;++j)seed[j]=(uint8_t)(j*17+3);
 for(unsigned dom=0;dom<256;++dom){copa_xof8 s;copa_xof8_seed_init(&s,seed,(uint8_t)dom);
  for(unsigned lane=0;lane<8;++lane){memcpy(inputs[lane],seed,32);inputs[lane][32]=(uint8_t)dom;inputs[lane][33]=(uint8_t)lane;inputs[lane][34]=0;shake256(expected[lane],272,inputs[lane],35);}
  for(unsigned b=0;b<2;++b){copa_xof8_block(&s,block);for(unsigned lane=0;lane<8;++lane)assert(!memcmp(block[lane],expected[lane]+136*b,136));}
 }
 puts("PASS eight-lane SHAKE256: eight input lengths, eight lanes, six blocks each");return 0;
}
