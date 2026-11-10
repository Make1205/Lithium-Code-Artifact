#include "copa_irs64.h"
#include "copa_round.h"
#ifdef COPA_VEC_POLY_HEADER
#include COPA_VEC_POLY_HEADER
#endif
#include <valgrind/memcheck.h>
#include <stdint.h>
static volatile uint64_t sink;
int main(void){
 for(unsigned which=0;which<4;++which){uint64_t norm=COPA_MIN_NORM+(COPA_BK_SQ-COPA_MIN_NORM)*which/3;
  copa_irs64_context ctx;VALGRIND_MAKE_MEM_UNDEFINED(&norm,sizeof(norm));copa_irs64_prepare(&ctx,norm);VALGRIND_MAKE_MEM_DEFINED(&norm,sizeof(norm));int64_t dot=(int64_t)norm-1;uint64_t random[6]={1,2,3,4,5,6};unsigned minus;
  VALGRIND_MAKE_MEM_UNDEFINED(&ctx,sizeof(ctx));VALGRIND_MAKE_MEM_UNDEFINED(&dot,sizeof(dot));VALGRIND_MAKE_MEM_UNDEFINED(random,sizeof(random));
  unsigned status=copa_irs64_step(&minus,dot,&ctx,random);VALGRIND_MAKE_MEM_DEFINED(&status,sizeof(status));VALGRIND_MAKE_MEM_DEFINED(&minus,sizeof(minus));sink=status+minus;
 }
#ifdef COPA_VEC_POLY_HEADER
 for(unsigned first=0;first<COPA_Q;first+=COPAV_W){int32_t x[COPAV_W],h[COPAV_W];for(unsigned j=0;j<COPAV_W;++j)x[j]=(first+j)%COPA_Q;
  VALGRIND_MAKE_MEM_UNDEFINED(x,sizeof(x));cv_store(h,copa_poly_high(cv_load(x)));VALGRIND_MAKE_MEM_DEFINED(h,sizeof(h));sink=(uint32_t)h[0];
 }
#else
 for(int32_t v=0;v<COPA_Q;++v){int32_t x=v;VALGRIND_MAKE_MEM_UNDEFINED(&x,sizeof(x));int32_t y=copa_highbits(x);VALGRIND_MAKE_MEM_DEFINED(&y,sizeof(y));sink=(uint32_t)y;}
#endif
 return 0;
}
