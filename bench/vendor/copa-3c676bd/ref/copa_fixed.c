#include "copa_fixed.h"
static copa_u128 umul(copa_u128 a,copa_u128 b)
{
    uint64_t a0=(uint64_t)a,b0=(uint64_t)b,a1=(uint64_t)(a>>64),b1=(uint64_t)(b>>64);
    copa_u128 p00=(copa_u128)a0*b0,p01=(copa_u128)a0*b1,p10=(copa_u128)a1*b0,p11=(copa_u128)a1*b1;
    copa_u128 mid=(p00>>64)+(uint64_t)p01+(copa_u128)(uint64_t)p10;
    copa_u128 hi=p11+(p01>>64)+(p10>>64)+(mid>>64);
    return (hi<<8)|((uint64_t)mid>>56);
}
copa_i128 copa_mul_fixed(copa_i128 a,copa_i128 b)
{
    copa_u128 ma=0-(copa_u128)(a<0),mb=0-(copa_u128)(b<0);
    copa_u128 aa=((copa_u128)a^ma)-ma,bb=((copa_u128)b^mb)-mb;
    copa_u128 v=umul(aa,bb),m=ma^mb;
    return (copa_i128)((v^m)-m); /* truncate towards zero */
}
copa_u128 copa_ratio_fixed(uint64_t num,uint32_t den)
{
    /* All division operands are public. Quotient and fractional expansion
     * of secret num use only fixed-schedule multiply, shifts and masks. */
    uint64_t reciprocal=(uint64_t)(((copa_u128)1<<64)/den);
    uint64_t whole=(uint64_t)(((copa_u128)num*reciprocal)>>64);
    uint64_t rem=num-whole*den;whole+=(rem>=den);rem=num-whole*den;
    uint64_t capped=(whole>=88),mask=0-capped;
    whole=(whole&~mask)|(88&mask);rem&=~mask;
    copa_u128 frac=0;
    for(unsigned i=0;i<120;++i){rem<<=1;unsigned bit=rem>=den;rem-=(uint64_t)den*bit;frac=(frac<<1)|bit;}
    return ((copa_u128)whole<<120)|frac;
}
static copa_u128 polynomial(copa_i128 t)
{
    copa_i128 v=(copa_i128)copa_exp_coeff[48];
    for(int k=47;k>=0;--k)v=(copa_i128)copa_exp_coeff[k]+copa_mul_fixed(v,t);
    return (copa_u128)v;
}
copa_u128 copa_exp_negative(copa_u128 a)
{
    unsigned n=0;
    for(unsigned k=1;k<=127;++k)n+=(a>=k*COPA_LN2);
    copa_u128 rem=a-n*COPA_LN2;
    return polynomial(-(copa_i128)rem)>>n;
}
copa_u128 copa_exp_positive(copa_u128 a){return polynomial((copa_i128)a);}
