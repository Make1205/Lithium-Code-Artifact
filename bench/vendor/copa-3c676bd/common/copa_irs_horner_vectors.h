#ifndef COPA_IRS_HORNER_VECTORS_H
#define COPA_IRS_HORNER_VECTORS_H
static inline void copa_backend_positive_horner(uint64_t out[4],const uint64_t *coeff,unsigned terms,uint64_t w,uint64_t w8){
 __m256i sums=_mm256_setzero_si256();
#if COPA_SIMD_BITS==512
 __m512i factors=_mm512_inserti64x4(_mm512_castsi256_si512(_mm256_set1_epi64x((int64_t)w)),_mm256_set1_epi64x((int64_t)w8),1);
 for(int k=(int)terms/8-1;k>=0;--k){__m256i even=_mm256_loadu_si256((const __m256i *)(coeff+4*k)),odd=_mm256_loadu_si256((const __m256i *)(coeff+terms/2+4*k));
  __m512i paired=_mm512_inserti64x4(_mm512_castsi256_si512(odd),sums,1),products=copa_fp_mul8(paired,factors);sums=_mm256_add_epi64(_mm256_sub_epi64(even,_mm512_castsi512_si256(products)),_mm512_extracti64x4_epi64(products,1));
 }
#else
 __m256i factor=_mm256_set1_epi64x((int64_t)w),power=_mm256_set1_epi64x((int64_t)w8);
 for(int k=(int)terms/8-1;k>=0;--k){__m256i even=_mm256_loadu_si256((const __m256i *)(coeff+4*k)),odd=_mm256_loadu_si256((const __m256i *)(coeff+terms/2+4*k));sums=_mm256_add_epi64(_mm256_sub_epi64(even,copa_fp_mul4(odd,factor)),copa_fp_mul4(sums,power));}
#endif
 _mm256_storeu_si256((__m256i *)out,sums);
}
#endif
