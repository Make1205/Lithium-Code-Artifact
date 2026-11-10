#ifndef COPA_GAUSS_VECTORS_H
#define COPA_GAUSS_VECTORS_H
#include "copa_cdt_vectors.h"
#include "copa_gaussian_fp.h"
#include "copa_sampler_tables.h"
/* Full table scan; carry a 32-bit-limb borrow in 64-bit vector lanes. */
static inline void copa_base96_vectors(uint32_t out[32],const uint64_t u[][32]){
 for(unsigned g=0;g<32;g+=CG_LANES){
  copav_vec words[COPA_BASE_LIMBS32],sum=cv_set64(0);
  for(unsigned l=0;l<COPA_BASE_LIMBS32;++l)words[l]=cv_and(cv_load(u[l]+g),cv_set64(UINT64_C(0xffffffff)));
  for(unsigned k=0;k<COPA_BASE_TAIL;++k){
   copav_vec borrow=cv_set64(0);
   for(unsigned l=0;l<COPA_BASE_LIMBS32;++l)
    borrow=cv_srl64(cv_sub64(cv_sub64(words[l],cv_set64(copa_base_rcdt[k][l])),borrow),63);
   sum=cv_add64(sum,borrow);
  }
  uint64_t x[CG_LANES];cv_store(x,sum);for(unsigned j=0;j<CG_LANES;++j)out[g+j]=(uint32_t)x[j];
 }
}
#if COPA_SIMD_BITS==256
static inline __m256d cg_round256(__m256d x){__asm__("":"+x"(x));return x;}
typedef __m256d cg_fp;
typedef struct{__m256d r,rl;__m128i ni;} cg_exp_state;
static inline cg_exp_state cg_exp_prepare(copav_vec nums){
 __m128i input=_mm256_castsi256_si128(nums);__m128i high=_mm256_extracti128_si256(nums,1);input=_mm_unpacklo_epi64(_mm_shuffle_epi32(input,0x88),_mm_shuffle_epi32(high,0x88));
 __m256d value=_mm256_cvtepi32_pd(input),x=cg_round256(_mm256_mul_pd(value,_mm256_set1_pd(CG_INV_HI)));
 __m256d xe=_mm256_fmadd_pd(value,_mm256_set1_pd(CG_INV_LO),_mm256_fmadd_pd(value,_mm256_set1_pd(CG_INV_HI),_mm256_sub_pd(_mm256_setzero_pd(),x)));
 __m128i ni=_mm256_cvttpd_epi32(_mm256_add_pd(cg_round256(_mm256_mul_pd(x,_mm256_set1_pd(CG_LOG2E))),_mm256_set1_pd(.5)));__m256d n=_mm256_cvtepi32_pd(ni),prod=cg_round256(_mm256_mul_pd(n,_mm256_set1_pd(CG_LN2_HI))),pe=_mm256_fmadd_pd(n,_mm256_set1_pd(CG_LN2_HI),_mm256_sub_pd(_mm256_setzero_pd(),prod));
 __m256d r=_mm256_sub_pd(x,prod),z=_mm256_sub_pd(r,x),re=_mm256_add_pd(_mm256_sub_pd(x,_mm256_sub_pd(r,z)),_mm256_sub_pd(_mm256_sub_pd(_mm256_setzero_pd(),prod),z)),rl=_mm256_add_pd(_mm256_fmadd_pd(_mm256_sub_pd(_mm256_setzero_pd(),n),_mm256_set1_pd(CG_LN2_LO),_mm256_sub_pd(xe,pe)),re);

 cg_exp_state st={r,rl,ni};return st;
}
static inline copav_vec cg_exp_finish(cg_fp p,cg_exp_state st){p=_mm256_fmadd_pd(_mm256_sub_pd(_mm256_setzero_pd(),st.rl),p,p);
 copav_vec bits=_mm256_sub_epi64(_mm256_castpd_si256(p),_mm256_slli_epi64(_mm256_cvtepi32_epi64(st.ni),52)),exponent=_mm256_sub_epi64(_mm256_and_si256(_mm256_srli_epi64(bits,52),_mm256_set1_epi64x(2047)),_mm256_set1_epi64x(1012));return _mm256_sllv_epi64(_mm256_or_si256(_mm256_and_si256(bits,_mm256_set1_epi64x(UINT64_C(0xfffffffffffff))),_mm256_set1_epi64x(UINT64_C(0x10000000000000))),exponent);

}
#else
static inline __m512d cg_round512(__m512d x){__asm__("":"+v"(x));return x;}
typedef __m512d cg_fp;
typedef struct{__m512d r,rl;__m256i ni;} cg_exp_state;
static inline cg_exp_state cg_exp_prepare(copav_vec nums){
 __m256i input=_mm512_cvtepi64_epi32(nums);__m512d value=_mm512_cvtepi32_pd(input),x=cg_round512(_mm512_mul_pd(value,_mm512_set1_pd(CG_INV_HI)));__m512d xe=_mm512_fmadd_pd(value,_mm512_set1_pd(CG_INV_LO),_mm512_fmadd_pd(value,_mm512_set1_pd(CG_INV_HI),_mm512_sub_pd(_mm512_setzero_pd(),x)));
 __m256i ni=_mm512_cvttpd_epi32(_mm512_add_pd(cg_round512(_mm512_mul_pd(x,_mm512_set1_pd(CG_LOG2E))),_mm512_set1_pd(.5)));__m512d n=_mm512_cvtepi32_pd(ni),prod=cg_round512(_mm512_mul_pd(n,_mm512_set1_pd(CG_LN2_HI))),pe=_mm512_fmadd_pd(n,_mm512_set1_pd(CG_LN2_HI),_mm512_sub_pd(_mm512_setzero_pd(),prod));__m512d r=_mm512_sub_pd(x,prod),z=_mm512_sub_pd(r,x),re=_mm512_add_pd(_mm512_sub_pd(x,_mm512_sub_pd(r,z)),_mm512_sub_pd(_mm512_sub_pd(_mm512_setzero_pd(),prod),z)),rl=_mm512_add_pd(_mm512_fmadd_pd(_mm512_sub_pd(_mm512_setzero_pd(),n),_mm512_set1_pd(CG_LN2_LO),_mm512_sub_pd(xe,pe)),re);

 cg_exp_state st={r,rl,ni};return st;
}
static inline copav_vec cg_exp_finish(cg_fp p,cg_exp_state st){p=_mm512_fmadd_pd(_mm512_sub_pd(_mm512_setzero_pd(),st.rl),p,p);copav_vec bits=_mm512_sub_epi64(_mm512_castpd_si512(p),_mm512_slli_epi64(_mm512_cvtepi32_epi64(st.ni),52)),exponent=_mm512_sub_epi64(_mm512_and_si512(_mm512_srli_epi64(bits,52),_mm512_set1_epi64(2047)),_mm512_set1_epi64(1012));return _mm512_sllv_epi64(_mm512_or_si512(_mm512_and_si512(bits,_mm512_set1_epi64(UINT64_C(0xfffffffffffff))),_mm512_set1_epi64(UINT64_C(0x10000000000000))),exponent);

}
#endif
#ifdef COPA_GAUSS_INTEGER
static inline copav_vec copa_gauss_exp_vector(copav_vec nums){
 uint64_t in[CG_LANES],out[CG_LANES];cv_store(in,nums);
 for(unsigned j=0;j<CG_LANES;j+=4){copa_fp63 p[4];copa64_exp_ratio4(p,in+j,&copa64_r2);for(unsigned k=0;k<4;++k)out[j+k]=copa64_to_q63(p[k]);}
 return cv_load(out);
}
static inline void copa_gauss_exp_pair(copav_vec out[2],const copav_vec nums[2]){out[0]=copa_gauss_exp_vector(nums[0]);out[1]=copa_gauss_exp_vector(nums[1]);}
#else
static inline copav_vec copa_gauss_exp_vector(copav_vec nums){cg_exp_state st=cg_exp_prepare(nums);
#if COPA_SIMD_BITS==256
 cg_fp p=_mm256_set1_pd(cg_exp_coeff[11]);for(int k=10;k>=0;--k)p=_mm256_fmadd_pd(p,st.r,_mm256_set1_pd(cg_exp_coeff[k]));
#else
 cg_fp p=_mm512_set1_pd(cg_exp_coeff[11]);for(int k=10;k>=0;--k)p=_mm512_fmadd_pd(p,st.r,_mm512_set1_pd(cg_exp_coeff[k]));
#endif
 return cg_exp_finish(p,st);
}
static inline void copa_gauss_exp_pair(copav_vec out[2],const copav_vec nums[2]){cg_exp_state a=cg_exp_prepare(nums[0]),b=cg_exp_prepare(nums[1]);
#if COPA_SIMD_BITS==256
 cg_fp p=_mm256_set1_pd(cg_exp_coeff[11]),q=p;for(int k=10;k>=0;--k){cg_fp c=_mm256_set1_pd(cg_exp_coeff[k]);p=_mm256_fmadd_pd(p,a.r,c);q=_mm256_fmadd_pd(q,b.r,c);}
#else
 cg_fp p=_mm512_set1_pd(cg_exp_coeff[11]),q=p;for(int k=10;k>=0;--k){cg_fp c=_mm512_set1_pd(cg_exp_coeff[k]);p=_mm512_fmadd_pd(p,a.r,c);q=_mm512_fmadd_pd(q,b.r,c);}
#endif
 out[0]=cg_exp_finish(p,a);out[1]=cg_exp_finish(q,b);
}
#endif
static inline void copa_backend_gaussian_batch(uint32_t candidate[32],uint32_t accepted[32],const uint64_t u[][32],const uint8_t low[32],const uint8_t tail[256]){
 uint32_t base[32];copa_base96_vectors(base,u);
 _Static_assert(((1ULL<<COPA_BASE_BITS)-1)*(((1ULL<<COPA_BASE_BITS)-1)+(2ULL<<COPA_BASE_BITS)*COPA_BASE_TAIL)*10000 < (1ULL<<31),"Gaussian exponents fit signed32 conversion");
 for(unsigned g=0;g<32;g+=COPAV_W){
  copav_vec b=cv_load(base+g),low32;
#if COPA_SIMD_BITS==256
  low32=_mm256_cvtepu8_epi32(_mm_loadl_epi64((const __m128i *)(low+g)));
#else
  low32=_mm512_cvtepu8_epi32(_mm_loadu_si128((const __m128i *)(low+g)));
#endif
  copav_vec values=cv_add(cv_sll32(b,COPA_BASE_BITS),low32);
  copav_vec nums=cv_mul32(cv_mul32(low32,cv_add(low32,cv_sll32(b,COPA_BASE_BITS+1))),cv_set(10000));
  cv_store(candidate+g,values);copav_vec inputs[2],v64[2],p[2];
#if COPA_SIMD_BITS==256
  inputs[0]=_mm256_cvtepu32_epi64(_mm256_castsi256_si128(nums));inputs[1]=_mm256_cvtepu32_epi64(_mm256_extracti128_si256(nums,1));
  v64[0]=_mm256_cvtepu32_epi64(_mm256_castsi256_si128(values));v64[1]=_mm256_cvtepu32_epi64(_mm256_extracti128_si256(values,1));
#else
  inputs[0]=_mm512_cvtepu32_epi64(_mm512_castsi512_si256(nums));inputs[1]=_mm512_cvtepu32_epi64(_mm512_extracti64x4_epi64(nums,1));
  v64[0]=_mm512_cvtepu32_epi64(_mm512_castsi512_si256(values));v64[1]=_mm512_cvtepu32_epi64(_mm512_extracti64x4_epi64(values,1));
#endif
  copa_gauss_exp_pair(p,inputs);
  for(unsigned lane=0;lane<2;++lane){copav_vec rd=cv_load(tail+8*(g+lane*CG_LANES));
   copav_vec threshold=cv_srl64(cv_add64(p[lane],cv_set64(64)),7);
   copav_vec sign=cv_and(cv_srl64(rd,56),cv_set64(1));
   copav_vec keep=cv_and(cg_ult(cv_and(rd,cv_set64(UINT64_C(0x00ffffffffffffff))),threshold),cv_or(cv_gt64(v64[lane],cv_set(0)),cg_eq(sign,cv_set64(0))));
   uint64_t yes[CG_LANES];cv_store(yes,keep);for(unsigned j=0;j<CG_LANES;++j)accepted[g+lane*CG_LANES+j]=(unsigned)(yes[j]&1);
  }
 }
}
#endif
