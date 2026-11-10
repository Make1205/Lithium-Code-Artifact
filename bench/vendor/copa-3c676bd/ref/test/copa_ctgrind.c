/* Dynamic undefined-value tracing: no secret-dependent branches/addresses in
 * fixed-schedule kernels. Rejection progress and returned statuses are outside
 * this test and are declassified before harness assertions/output. */
#include "copa_sign.h"
#include "copa_sampler.h"
#include "copa_scalar64.h"
#include "ntt_crt.h"
#include "poly_inv.h"
#include "quant.h"
#include <valgrind/memcheck.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
int main(int argc,char **argv){
 const char *which=argc>1?argv[1]:"inverse";int32_t f[256],g[256],out[256];int64_t product[256];unsigned ok=0;
 for(unsigned i=0;i<256;++i){f[i]=(int32_t)(i%7)-3;g[i]=(int32_t)(i%5)-2;}
 if(!strcmp(which,"inverse")){VALGRIND_MAKE_MEM_UNDEFINED(f,sizeof(f));ok=(unsigned)copa_poly_inv(out,f,COPA_Q);VALGRIND_MAKE_MEM_DEFINED(out,sizeof(out));}
 else if(!strcmp(which,"ntt")){VALGRIND_MAKE_MEM_UNDEFINED(f,sizeof(f));VALGRIND_MAKE_MEM_UNDEFINED(g,sizeof(g));copa_ntt_crt_mul(product,f,g,256);uint32_t a[1][2][256],b[1][2][256],v[2][256];copa_ntt_forward(a[0],f,256);copa_ntt_forward(b[0],g,256);copa_ntt_dot(v,a,b,1,256);copa_ntt_inverse_modq(out,v,256,COPA_Q);VALGRIND_MAKE_MEM_DEFINED(product,sizeof(product));VALGRIND_MAKE_MEM_DEFINED(out,sizeof(out));}
 else if(!strcmp(which,"irs")){copa_irs64_context ctx;copa_irs64_prepare(&ctx,63);int64_t dot=31;uint64_t random[6]={1,2,3,4,5,6};unsigned minus;
  VALGRIND_MAKE_MEM_UNDEFINED(&ctx,sizeof(ctx));VALGRIND_MAKE_MEM_UNDEFINED(&dot,sizeof(dot));VALGRIND_MAKE_MEM_UNDEFINED(random,sizeof(random));ok=copa_irs64_step(&minus,dot,&ctx,random);VALGRIND_MAKE_MEM_DEFINED(&minus,sizeof(minus));}
 else if(!strcmp(which,"irs_setup")){copa_irs64_context ctx;uint64_t norm=63;VALGRIND_MAKE_MEM_UNDEFINED(&norm,sizeof(norm));copa_irs64_prepare(&ctx,norm);VALGRIND_MAKE_MEM_DEFINED(&ctx,sizeof(ctx));}
 else if(!strcmp(which,"quant")){VALGRIND_MAKE_MEM_UNDEFINED(f,sizeof(f));copa_quant_block(out,f,COPA_MODE);VALGRIND_MAKE_MEM_DEFINED(out,sizeof(out));}
 else if(!strcmp(which,"noise")){uint8_t seed[32]={0};int32_t s[1+COPA_ELL][COPA_N];VALGRIND_MAKE_MEM_UNDEFINED(seed,sizeof(seed));copa_expand_s(s,seed);VALGRIND_MAKE_MEM_DEFINED(s,sizeof(s));}
 else if(!strcmp(which,"exp")){uint64_t num[4]={10000,100000,1000000,10000000};copa_fp63 result[4];
  VALGRIND_MAKE_MEM_UNDEFINED(num,sizeof(num));
  /* The context is public r, rather than the secret norm ratio. */
  extern void copa_ctgrind_exp4(copa_fp63[4],const uint64_t[4]);copa_ctgrind_exp4(result,num);VALGRIND_MAKE_MEM_DEFINED(result,sizeof(result));}
 else return 2;
 VALGRIND_MAKE_MEM_DEFINED(&ok,sizeof(ok));VALGRIND_MAKE_MEM_DEFINED(f,sizeof(f));VALGRIND_MAKE_MEM_DEFINED(g,sizeof(g));printf("mode=%d component=%s trace complete status=%u\n",COPA_MODE,which,ok);return 0;
}
