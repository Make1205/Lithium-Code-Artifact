#include "copa_gaussian_fp.h"
#include "copa_scalar64.h"
#include "copa_scalar64_tables.h"
void exp64(uint64_t out[2],uint64_t num){copa_fp63 v=copa64_exp_ratio(num,&copa64_r2);out[0]=v.m;out[1]=(uint32_t)v.e;}
void fp63_op(uint64_t out[2],unsigned op,uint64_t am,int32_t ae,uint64_t bm,int32_t be){
 copa_fp63 a={am,ae},b={bm,be},v;
 if(op==0)v=copa64_normalize(am,ae);else if(op==1)v=copa64_multiply(a,b);else v=copa64_add(a,b);
 out[0]=v.m;out[1]=(uint32_t)v.e;
}
unsigned coin63(uint64_t m,int32_t e,uint64_t a,uint64_t b,uint64_t c){copa_fp63 p={m,e};uint64_t u[3]={a,b,c};return copa64_bernoulli(p,u);}
uint64_t sampled_probability64(uint64_t num){copa_gaussian_env old=copa_gaussian_enter();uint64_t p=copa_gaussian_threshold(num);copa_gaussian_leave(old);return p;}
void exp64_batch(uint64_t out[8],const uint64_t num[4]){copa_fp63 p[4];copa64_exp_ratio4(p,num,&copa64_r2);for(unsigned j=0;j<4;++j){out[2*j]=p[j].m;out[2*j+1]=(uint32_t)p[j].e;}}
