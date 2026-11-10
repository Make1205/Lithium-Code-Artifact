#ifndef COPA_SCALAR64_H
#define COPA_SCALAR64_H
#include "copa_params.h"
#include <stdint.h>
#if defined(__GNUC__) || defined(__clang__)
#define COPA_FLATTEN __attribute__((flatten))
#else
#define COPA_FLATTEN
#endif
/* Value m*2^(e-63), m==0 or 2^62<=m<2^63. Horner arithmetic is Q63.
 * A separate binary scale prevents small nonzero probabilities losing bits. */
typedef struct {uint64_t m;int32_t e;} copa_fp63;
typedef struct {uint32_t den,shift;uint64_t reciprocal,fraction_reciprocal,log_reciprocal;} copa64_ratio_ctx;
/* One hardware 64x64 product; no __int128 type or software wide arithmetic. */
static inline uint64_t copa64_product(uint64_t a,uint64_t b,uint64_t *lo){
#if defined(__x86_64__)
 uint64_t hi;__asm__("mulq %3":"=a"(*lo),"=d"(hi):"a"(a),"r"(b):"cc");return hi;
#elif defined(__aarch64__)
 uint64_t hi;*lo=a*b;__asm__("umulh %0,%1,%2":"=r"(hi):"r"(a),"r"(b));return hi;
#else
#error "A 64-bit hardware multiply-high instruction is required"
#endif
}
static inline uint64_t copa64_mulq63(uint64_t a,uint64_t b){uint64_t lo,hi=copa64_product(a,b,&lo);return((hi<<1)|(lo>>63))+((lo>>62)&1);}
static inline uint64_t copa64_mask(unsigned b){
#if defined(__clang__)
 /* An opaque loop-invariant zero blocks LLVM's mask-to-branch rewrite. */
 uint64_t zero=0;__asm__("":"+r"(zero));return zero-(uint64_t)b;
#else
 return UINT64_C(0)-(uint64_t)b;
#endif
}
static inline uint64_t copa64_ct_mask(unsigned b){
 uint64_t mask=UINT64_C(0)-(uint64_t)b;
 /* Value barrier: the optimizer must not turn a mask into a secret branch. */
#if defined(__GNUC__) || defined(__clang__)
 __asm__("":"+r"(mask));
#endif
 return mask;
}
/* Avoid a legacy bit-scan instruction in baseline x86 builds. Native LZCNT
 * and AArch64 CLZ are fixed-latency; the portable fallback has six stages. */
static inline unsigned copa64_clz_nonzero(uint64_t x){
#if defined(__LZCNT__) || defined(__aarch64__)
 return (unsigned)__builtin_clzll(x);
#else
 unsigned n=0;for(unsigned shift=32;shift;shift>>=1){unsigned amount=shift&(unsigned)copa64_ct_mask((x>>(64-shift))==0);x<<=amount;n+=amount;}return n;
#endif
}
static inline uint64_t copa64_select(uint64_t a,uint64_t b,unsigned bit){uint64_t m=copa64_ct_mask(bit);return(a&m)|(b&~m);}
/* Optional backend header supplies bit-exact parallel Q63 products. */
#ifdef COPA_VEC_FP_HEADER
#include COPA_VEC_FP_HEADER
#endif
copa_fp63 copa64_exp_ratio(uint64_t num,const copa64_ratio_ctx *ctx);
void copa64_exp_ratio4(copa_fp63 out[4],const uint64_t num[4],const copa64_ratio_ctx *ctx);
#ifdef COPA_VEC_FP_HEADER
void copa64_irs_polynomials(copa_fp63 *factor,uint64_t out[3],uint64_t num,uint64_t tq,uint64_t cq,const copa64_ratio_ctx *ctx);
#endif
uint64_t copa64_fraction(uint64_t rem,const copa64_ratio_ctx *ctx);
copa_fp63 copa64_normalize(uint64_t m,int32_t e);
copa_fp63 copa64_multiply(copa_fp63 a,copa_fp63 b);
copa_fp63 copa64_add(copa_fp63 a,copa_fp63 b);
uint64_t copa64_to_q63(copa_fp63 a);
unsigned copa64_bernoulli(copa_fp63 p,const uint64_t uniform[3]);
#endif
