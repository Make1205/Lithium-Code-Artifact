#ifndef COPA_GAUSSIAN_FP_H
#define COPA_GAUSSIAN_FP_H
#include "copa_params.h"
#include "copa_gaussian_fp_tables.h"
#include <math.h>
#include <stdint.h>
#include <string.h>
#if defined(__x86_64__)
#include <xmmintrin.h>
typedef unsigned copa_gaussian_env;
static inline copa_gaussian_env copa_gaussian_enter(void){unsigned old=_mm_getcsr();unsigned safe=0x1f80;__asm__ volatile("ldmxcsr %0"::"m"(safe):"memory");return old;}
static inline void copa_gaussian_leave(copa_gaussian_env old){__asm__ volatile("ldmxcsr %0"::"m"(old):"memory");}
#else
#include <fenv.h>
typedef fenv_t copa_gaussian_env;
static inline copa_gaussian_env copa_gaussian_enter(void){fenv_t old;feholdexcept(&old);fesetround(FE_TONEAREST);return old;}
static inline void copa_gaussian_leave(copa_gaussian_env old){fesetenv(&old);}
#endif
static inline double copa_gaussian_rounded(double x){
#if defined(__x86_64__)
 __asm__("":"+x"(x));
#elif defined(__aarch64__)
 __asm__("":"+w"(x));
#else
 volatile double y=x;return y;
#endif
 return x;
}
/* num is an exact integer <2^31. Compensated reciprocal and ln(2)
 * reduction avoid division. Every intermediate is normal or zero.
 * IEEE FMA semantics define the same polynomial on scalar/SIMD backends. */
static inline uint64_t copa_gaussian_threshold(uint64_t num){
 double value=(double)(uint32_t)num,x=copa_gaussian_rounded(value*CG_INV_HI),xe=fma(value,CG_INV_LO,fma(value,CG_INV_HI,-x));int n=(int)(copa_gaussian_rounded(x*CG_LOG2E)+.5);double prod=copa_gaussian_rounded(n*CG_LN2_HI),pe=fma((double)n,CG_LN2_HI,-prod),r=x-prod,z=r-x,re=(x-(r-z))+(-prod-z),rl=fma(-(double)n,CG_LN2_LO,xe-pe)+re;
 double p=cg_exp_coeff[11];for(int k=10;k>=0;--k)p=fma(p,r,cg_exp_coeff[k]);p=fma(-rl,p,p);uint64_t bits;memcpy(&bits,&p,8);bits-=(uint64_t)n<<52;unsigned exponent=(unsigned)(bits>>52)&2047;return((bits&UINT64_C(0xfffffffffffff))|UINT64_C(0x10000000000000))<<(exponent-1012);
}
/* V3 compares an independent 56-bit uniform with this actual threshold.
 * Round to nearest, ties up, after the existing binary64/Q63 evaluation. */
#include "copa_scalar64.h"
#include "copa_scalar64_tables.h"
static inline uint64_t copa_gaussian_threshold56(uint64_t num){
#ifdef COPA_GAUSS_INTEGER
 return (copa64_to_q63(copa64_exp_ratio(num,&copa64_r2))+64)>>7;
#else
 return (copa_gaussian_threshold(num)+64)>>7;
#endif
}
#endif
