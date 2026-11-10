#ifndef COPA_SAMPLER_EXP_H
#define COPA_SAMPLER_EXP_H
#include "copa_gaussian_fp.h"
#include "copa_scalar64.h"
/* IEEE encoding of the same compensated FMA polynomial, without Q63
 * quantization. Normalize its significand for an exact 192-bit Bernoulli. */
static inline uint64_t copa_exp_probability_bits(uint64_t num){
 double value=(double)num,x=copa_gaussian_rounded(value*CG_INV_HI),xe=fma(value,CG_INV_LO,fma(value,CG_INV_HI,-x));int n=(int)(copa_gaussian_rounded(x*CG_LOG2E)+.5);double prod=copa_gaussian_rounded(n*CG_LN2_HI),pe=fma((double)n,CG_LN2_HI,-prod),r=x-prod,z=r-x,re=(x-(r-z))+(-prod-z),rl=fma(-(double)n,CG_LN2_LO,xe-pe)+re;
 double p=cg_exp_coeff[11];for(int k=10;k>=0;--k)p=fma(p,r,cg_exp_coeff[k]);p=fma(-rl,p,p);uint64_t bits;memcpy(&bits,&p,8);return bits-((uint64_t)n<<52);
}
static inline unsigned copa_exp_accept(uint64_t bits,const uint64_t u[3]){
 copa_fp63 p={((bits&UINT64_C(0xfffffffffffff))|UINT64_C(0x10000000000000))<<10,(int32_t)(bits>>52)-1022};return copa64_bernoulli(p,u);
}
#endif
