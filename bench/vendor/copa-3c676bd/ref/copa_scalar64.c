#include "copa_scalar64.h"
#include "copa_scalar64_tables.h"
uint64_t copa64_fraction(uint64_t rem,const copa64_ratio_ctx *c){
 uint64_t lo,hi=copa64_product(rem,c->fraction_reciprocal,&lo);
 uint64_t q=(hi<<(64-c->shift))|(lo>>c->shift),err=0-q*c->den;
 unsigned b=err>=c->den;q+=b;err-=(uint64_t)b*c->den;
 b=err>=c->den;q+=b;err-=(uint64_t)b*c->den;return q+(2*err>=c->den);
}
copa_fp63 copa64_normalize(uint64_t m,int32_t e){
 unsigned s=copa64_clz_nonzero(m|1);unsigned high=(s==0);unsigned left=(s-1)&63;
 copa_fp63 r={copa64_select(m>>1,m<<left,high),e+(int32_t)high-(int32_t)((s-1)&(unsigned)copa64_ct_mask(!high))};
 uint64_t nonzero=copa64_ct_mask(m!=0);r.m&=nonzero;r.e=(int32_t)((uint32_t)r.e&(uint32_t)nonzero);return r;
}
copa_fp63 copa64_multiply(copa_fp63 a,copa_fp63 b){return copa64_normalize(copa64_mulq63(a.m,b.m),a.e+b.e);}
uint64_t copa64_to_q63(copa_fp63 a){
 /* Callers use values <=1; positive e occurs only for exact one. */
 unsigned right=(unsigned)(-a.e)&63;uint64_t rounding=((a.m>>((right-1)&63))&1)&copa64_mask((a.e<0)&(a.e>-64));
 uint64_t v=((a.m>>right)+rounding)&copa64_mask(a.e>-64);
 return copa64_select(a.m<<((unsigned)a.e&63),v,a.e>0);
}
copa_fp63 copa64_add(copa_fp63 a,copa_fp63 b){
 unsigned take=(unsigned)((copa64_ct_mask(a.m==0)|(copa64_ct_mask(b.m!=0)&copa64_ct_mask(a.e<b.e)))&1);copa_fp63 x={copa64_select(b.m,a.m,take),(int32_t)copa64_select((uint32_t)b.e,(uint32_t)a.e,take)};
 copa_fp63 y={copa64_select(a.m,b.m,take),(int32_t)copa64_select((uint32_t)a.e,(uint32_t)b.e,take)};
 unsigned sh=(unsigned)(x.e-y.e);uint64_t ym=(y.m>>(sh&63))&copa64_mask(sh<64);
 return copa64_normalize(x.m+ym,x.e);
}
typedef struct {uint64_t x,mask,table;int32_t n;unsigned zero,over;} exp64_state;
static exp64_state exp_prepare(uint64_t num,const copa64_ratio_ctx *c){
 uint64_t limit=(uint64_t)c->den*88,over=num>=limit;num=copa64_select(limit,num,(unsigned)over);
 uint64_t lo,whole=copa64_product(num,c->reciprocal,&lo),rem=num-whole*c->den;
 unsigned b=rem>=c->den;whole+=b;rem-=(uint64_t)b*c->den;
 uint64_t fraction=copa64_fraction(rem,c);
 uint64_t m=(num*c->log_reciprocal)>>48;m+=8-(num==0);
 uint64_t lp,lh=copa64_product(m,COPA64_LN2,&lp),fl=fraction<<3;
 uint64_t rl=lp-fl,ih=lh-whole*8-(fraction>>61)-(lp<fl);
 uint64_t r=(ih<<61)|(rl>>3);
 unsigned correction=r<COPA64_RMIN;m+=correction;r+=COPA64_STEP*(uint64_t)correction;
 uint64_t n=m>>3,j=m&7,diff=r-COPA64_CENTER,sign=diff>>63;
 uint64_t mask=copa64_mask((unsigned)sign),x=((diff^mask)-mask)>>1;
 uint64_t table=0;for(unsigned k=0;k<8;++k)table|=copa64_power_eighth[k]&copa64_mask(j==k);
 exp64_state st={x,mask,table,(int32_t)n,num==0,(unsigned)over};return st;
}
static copa_fp63 exp_finish(uint64_t p,const exp64_state *st){
 copa_fp63 out=copa64_normalize(copa64_mulq63(p,st->table),-st->n);
 copa_fp63 unit={UINT64_C(1)<<62,1};out.m=copa64_select(unit.m,out.m,st->zero);out.e=(int32_t)copa64_select((uint32_t)unit.e,(uint32_t)out.e,st->zero);
 out.m&=~copa64_mask(st->over);return out;
}
copa_fp63 copa64_exp_ratio(uint64_t num,const copa64_ratio_ctx *c){
 exp64_state st=exp_prepare(num,c);uint64_t p=copa64_center_coeff[9];
 for(int k=8;k>=0;--k){uint64_t v=copa64_mulq63(p,st.x);p=copa64_center_coeff[k]+((v^st.mask)-st.mask);}
 return exp_finish(p,&st);
}
COPA_FLATTEN void copa64_exp_ratio4(copa_fp63 out[4],const uint64_t num[4],const copa64_ratio_ctx *c){
 exp64_state st[4];uint64_t p[4];for(unsigned i=0;i<4;++i){st[i]=exp_prepare(num[i],c);p[i]=copa64_center_coeff[9];}
 /* Independent lane chains expose multiply instruction-level parallelism. */
#ifdef COPA_VEC_FP_HEADER
 uint64_t x[4];for(unsigned i=0;i<4;++i)x[i]=st[i].x;
 for(int k=8;k>=0;--k){uint64_t v[4];copa_backend_mulq63x4(v,p,x);for(unsigned i=0;i<4;++i)p[i]=copa64_center_coeff[k]+((v[i]^st[i].mask)-st[i].mask);}
#else
 for(int k=8;k>=0;--k)for(unsigned i=0;i<4;++i){uint64_t v=copa64_mulq63(p[i],st[i].x);p[i]=copa64_center_coeff[k]+((v^st[i].mask)-st[i].mask);}
#endif
 for(unsigned i=0;i<4;++i)out[i]=exp_finish(p[i],&st[i]);
}
#ifdef COPA_VEC_FP_HEADER
#include "../common/copa_irs_polynomial_tables.h"
/* Four independent per-step chains. Leading zero coefficients preserve each
 * original Horner order, including all intermediate Q63 rounding. */
COPA_FLATTEN void copa64_irs_polynomials(copa_fp63 *factor,uint64_t out[3],uint64_t num,uint64_t tq,uint64_t cq,const copa64_ratio_ctx *ctx){
 exp64_state st=exp_prepare(num,ctx);
 uint64_t p[4];
#if COPA_SIMD_BITS == 256
 __m256i x=_mm256_set_epi64x((int64_t)st.x,(int64_t)cq,(int64_t)tq,(int64_t)tq);
 __m256i mask=_mm256_set_epi64x((int64_t)st.mask,-1,-1,-1);
 __m256i v=_mm256_loadu_si256((const __m256i *)copa64_irs_poly_coeff[20]);
 for(int k=19;k>=0;--k){__m256i q=copa_fp_mul4(v,x);v=_mm256_add_epi64(_mm256_loadu_si256((const __m256i *)copa64_irs_poly_coeff[k]),_mm256_sub_epi64(_mm256_xor_si256(q,mask),mask));}
 _mm256_storeu_si256((__m256i *)p,v);
#else
 __m256i x4=_mm256_set_epi64x((int64_t)st.x,(int64_t)cq,(int64_t)tq,(int64_t)tq),m4=_mm256_set_epi64x((int64_t)st.mask,-1,-1,-1);
 __m512i x=_mm512_broadcast_i64x4(x4),mask=_mm512_broadcast_i64x4(m4);
 __m512i v=_mm512_broadcast_i64x4(_mm256_loadu_si256((const __m256i *)copa64_irs_poly_coeff[20]));
 for(int k=19;k>=0;--k){__m512i q=copa_fp_mul8(v,x);v=_mm512_add_epi64(_mm512_broadcast_i64x4(_mm256_loadu_si256((const __m256i *)copa64_irs_poly_coeff[k])),_mm512_sub_epi64(_mm512_xor_si512(q,mask),mask));}
 _mm256_storeu_si256((__m256i *)p,_mm512_castsi512_si256(v));
#endif
 for(unsigned j=0;j<3;++j)out[j]=p[j];
 *factor=exp_finish(p[3],&st);
}
#endif
unsigned copa64_bernoulli(copa_fp63 p,const uint64_t u[3]){
 /* The supported exponents require at most 192 uniform bits. No secret
  * indexed loads: scan/select all three words and all preceding words. */
 unsigned n=(unsigned)(-p.e),block=n>>6,s=n&63,sh=(64-s)&63;
 uint64_t a=0,b=0,prefix=0;
 for(unsigned j=0;j<3;++j){a|=u[j]&copa64_ct_mask(block==j);b|=u[j]&copa64_ct_mask(block+1==j);prefix|=u[j]&copa64_ct_mask(j<block);}
 uint64_t top=(a<<s)|((b>>sh)&copa64_ct_mask(s!=0));
 prefix|=(a>>sh)&copa64_ct_mask(s!=0);
 unsigned accept=(prefix==0)&((top>>1)<p.m)&(n<=129)&(p.m!=0);
 return accept|(unsigned)((p.e>0)&(p.m!=0));
}
