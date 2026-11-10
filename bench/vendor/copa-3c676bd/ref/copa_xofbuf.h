#ifndef COPA_XOFBUF_H
#define COPA_XOFBUF_H
#include "fips202.h"
#include <string.h>
typedef struct {keccak_state state;uint8_t block[SHAKE256_RATE];unsigned pos;} copa_xofbuf;
static inline void copa_xof_init(copa_xofbuf *s,const uint8_t *in,size_t len){shake256_absorb_once(&s->state,in,len);s->pos=SHAKE256_RATE;}
static inline void copa_xof_read(copa_xofbuf *s,uint8_t *out,size_t n){
 while(n){if(s->pos==SHAKE256_RATE){shake256_squeezeblocks(s->block,1,&s->state);s->pos=0;}
  size_t take=SHAKE256_RATE-s->pos;if(take>n)take=n;memcpy(out,s->block+s->pos,take);out+=take;s->pos+=(unsigned)take;n-=take;
 }
}
static inline unsigned copa_xof_u16(copa_xofbuf *s){
 /* Used by streams consuming only two-byte candidates. */
 if(s->pos==SHAKE256_RATE){shake256_squeezeblocks(s->block,1,&s->state);s->pos=0;}
 unsigned v=s->block[s->pos]+256u*s->block[s->pos+1];s->pos+=2;return v;
}
static inline uint8_t copa_xof_u8(copa_xofbuf *s){
 if(s->pos==SHAKE256_RATE){shake256_squeezeblocks(s->block,1,&s->state);s->pos=0;}return s->block[s->pos++];
}
#endif
