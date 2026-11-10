#include "copa_fixed.h"
void ratio(uint64_t out[2],uint64_t n,uint32_t d){copa_u128 v=copa_ratio_fixed(n,d);out[0]=(uint64_t)v;out[1]=(uint64_t)(v>>64);}
void expq(uint64_t out[2],uint64_t lo,uint64_t hi,unsigned positive){copa_u128 a=(copa_u128)lo|((copa_u128)hi<<64),v=positive?copa_exp_positive(a):copa_exp_negative(a);out[0]=(uint64_t)v;out[1]=(uint64_t)(v>>64);}
void mulq(uint64_t out[2],uint64_t al,uint64_t ah,uint64_t bl,uint64_t bh){copa_i128 a=(copa_i128)((copa_u128)al|((copa_u128)ah<<64)),b=(copa_i128)((copa_u128)bl|((copa_u128)bh<<64));copa_u128 v=(copa_u128)copa_mul_fixed(a,b);out[0]=(uint64_t)v;out[1]=(uint64_t)(v>>64);}
