#include "copa_irs64.h"
#include "copa_scalar64_tables.h"
typedef struct {uint64_t bm,gm,reciprocal,fraction_reciprocal;int16_t be,ge;uint8_t shift;} copa64_norm_row;
#include "copa_irs64_tables.h"
#ifdef COPA_VEC_FP_HEADER
#include "../common/copa_irs_setup_vectors.h"
#include "../common/copa_irs_horner_vectors.h"
#endif
static copa_fp63 choose(copa_fp63 a,copa_fp63 b,unsigned bit){copa_fp63 r={copa64_select(a.m,b.m,bit),(int32_t)copa64_select((uint32_t)a.e,(uint32_t)b.e,bit)};return r;}
static uint64_t polynomial(const uint64_t *coeff,unsigned degree,uint64_t x){uint64_t v=coeff[degree];for(int k=(int)degree-1;k>=0;--k)v=coeff[k]-copa64_mulq63(v,x);return v;}
static copa_fp63 fraction(uint64_t num,const copa64_ratio_ctx *c){
 uint64_t lo,whole=copa64_product(num,c->reciprocal,&lo),rem=num-whole*c->den;
 unsigned b=rem>=c->den;whole+=b;rem-=(uint64_t)b*c->den;
 uint64_t f=copa64_fraction(rem,c);copa_fp63 small=copa64_normalize(f,-1);
 copa_fp63 big={ (UINT64_C(1)<<62)+(f>>2),1};
 copa_fp63 two={UINT64_C(1)<<62,2};big=choose(two,big,whole==2);return choose(big,small,whole!=0);
}
void copa_irs64_prepare(copa_irs64_context *ctx,uint64_t norm){
 ctx->norm=norm;ctx->norm_ratio.den=(uint32_t)norm;ctx->norm_ratio.reciprocal=ctx->norm_ratio.fraction_reciprocal=0;ctx->norm_ratio.shift=0;ctx->norm_ratio.log_reciprocal=0;
 ctx->bump.m=ctx->gap.m=0;uint32_t be=0,ge=0;
 /* Sequential scan keeps the secret norm out of table addresses. */
#ifdef COPA_VEC_FP_HEADER
 uint64_t selected[5];copa_backend_irs_constants(selected,norm);ctx->bump.m=selected[0];ctx->gap.m=selected[1];ctx->norm_ratio.reciprocal=selected[2];ctx->norm_ratio.fraction_reciprocal=selected[3];be=(uint32_t)(int32_t)(int16_t)selected[4];ge=(uint32_t)(int32_t)(int16_t)(selected[4]>>16);ctx->norm_ratio.shift=(unsigned)(selected[4]>>32);
#else
 for(unsigned k=1;k<=COPA_BK_SQ;++k){uint64_t mask=copa64_mask(k==norm);const copa64_norm_row *r=&irs64_rows[k-1];
  ctx->bump.m|=r->bm&mask;ctx->gap.m|=r->gm&mask;be|=(uint32_t)(int32_t)r->be&(uint32_t)mask;ge|=(uint32_t)(int32_t)r->ge&(uint32_t)mask;
  ctx->norm_ratio.reciprocal|=r->reciprocal&mask;ctx->norm_ratio.fraction_reciprocal|=r->fraction_reciprocal&mask;ctx->norm_ratio.shift|=r->shift&(uint32_t)mask;
 }
#endif
 ctx->bump.e=(int32_t)be;ctx->gap.e=(int32_t)ge;
 ctx->gamma=copa64_normalize(copa64_fraction(norm*10000,&copa64_r2),-1);
 for(unsigned k=0;k<COPA64_IRS_TERMS;k+=4){uint64_t num[4];copa_fp63 value[4];for(unsigned i=0;i<4;++i)num[i]=norm*10000*(k+i)*(k+i);copa64_exp_ratio4(value,num,&copa64_r2);for(unsigned i=0;i<4;++i)ctx->c[k+i]=copa64_to_q63(value[i]);}
 ctx->exp_gamma=ctx->c[1];
#ifdef COPA_VEC_FP_HEADER
 uint64_t ordered[COPA64_IRS_TERMS];
 for(unsigned j=0;j<COPA64_IRS_TERMS/2;++j){ordered[j]=ctx->c[2*j];ordered[COPA64_IRS_TERMS/2+j]=ctx->c[2*j+1];}
 for(unsigned j=0;j<COPA64_IRS_TERMS;++j)ctx->c[j]=ordered[j];
#endif
}
COPA_FLATTEN void copa_irs64_probabilities(copa_fp63 *abort,copa_fp63 *conditional_small,int64_t dot,const copa_irs64_context *c){
 uint64_t neg=copa64_mask(dot<0),a=((uint64_t)dot^neg)-neg;
 uint64_t clamped=copa64_select(a,c->norm,a<c->norm),hnum=c->norm-clamped;
 copa_fp63 h=fraction(hnum,&c->norm_ratio),oneplus=fraction(2*c->norm-hnum,&c->norm_ratio);
 copa_fp63 t=copa64_multiply(c->gamma,copa64_multiply(h,oneplus));uint64_t tq=copa64_to_q63(t);
 copa_fp63 pi={UINT64_C(7244019458077122842),2}; /* pi = (pi/4)*2^2 */
 copa_fp63 eta=copa64_multiply(pi,h),eta2=copa64_multiply(eta,eta),arg=eta2;arg.e-=4;
 copa_fp63 factor,cosc,phie,exp_t;
#ifdef COPA_VEC_FP_HEADER
 uint64_t poly[3];
 copa64_irs_polynomials(&factor,poly,(2*a+c->norm)*10000,tq,copa64_to_q63(arg),&copa64_r2);
 phie=copa64_normalize(poly[0],0);exp_t=copa64_normalize(poly[1],0);cosc=copa64_normalize(poly[2],0);
#else
 cosc=copa64_normalize(polynomial(copa64_cosc_coeff,16,copa64_to_q63(arg)),0);
 phie=copa64_normalize(polynomial(copa64_expm1_coeff,20,tq),0);
 exp_t=copa64_normalize(polynomial(copa64_exp_coeff,20,tq),0);
 factor=copa64_exp_ratio((2*a+c->norm)*10000,&copa64_r2);
#endif
 copa_fp63 one_minus_cos=copa64_multiply(eta2,cosc);
 copa_fp63 one_minus_exp=copa64_multiply(t,phie);
 copa_fp63 shape=copa64_add(one_minus_exp,copa64_multiply(exp_t,one_minus_cos));
 copa_fp63 central=copa64_multiply(copa64_add(c->gap,copa64_multiply(c->bump,shape)),irs64_inverse_m);
 unsigned inside=(dot>=-(int64_t)c->norm)&(dot<(int64_t)c->norm);
 *abort=choose(central,irs64_outside_abort,inside);
 /* 1-S(a) = exp(-a-gamma) S(a+2 gamma); this avoids subtracting
  * two Q63 numbers near one when the selected branch is tiny. */
 copa_fp63 gamma_factor={c->exp_gamma,0}; /* exp(-gamma) is in [1/2,1). */
 copa_fp63 base=copa64_multiply(factor,gamma_factor);
 uint64_t w=copa64_to_q63(base),w2=copa64_mulq63(w,w),w4=copa64_mulq63(w2,w2),w8=copa64_mulq63(w4,w4),sums[4]={0};
 /* Four independent positive-pair Horner chains, recombined in w^2.
  * All partials stay <=1; independent multiplies hide hardware latency. */
#ifdef COPA_VEC_FP_HEADER
 copa_backend_positive_horner(sums,c->c,COPA64_IRS_TERMS,w,w8);
#else
 for(int k=COPA64_IRS_TERMS/8-1;k>=0;--k)for(unsigned j=0;j<4;++j){uint64_t coefficient=c->c[8*k+2*j]-copa64_mulq63(c->c[8*k+2*j+1],w);sums[j]=coefficient+copa64_mulq63(sums[j],w8);}
#endif
 uint64_t s=sums[3];for(int j=2;j>=0;--j)s=sums[j]+copa64_mulq63(s,w2);
 copa_fp63 small=copa64_multiply(factor,copa64_normalize(s,0));small=copa64_multiply(small,irs64_inverse_m);
 /* 1/(1-abort) by a fixed geometric Horner polynomial, abort<0.04. */
 uint64_t aq=copa64_to_q63(*abort),inv=COPA64_ONE;
 for(unsigned k=0;k<16;++k)inv=COPA64_ONE+copa64_mulq63(inv,aq);
 *conditional_small=copa64_multiply(small,copa64_normalize(inv,0));
}
unsigned copa_irs64_step(unsigned *minus,int64_t dot,const copa_irs64_context *c,const uint64_t random[6]){
 copa_fp63 abort,small;copa_irs64_probabilities(&abort,&small,dot,c);
 unsigned stop=copa64_bernoulli(abort,random),sm=copa64_bernoulli(small,random+3);
 *minus=sm^(unsigned)(dot>=0);return 1-stop;
}
