/* Independent lattice blocks occupy SIMD lanes. All cosets, lexicographic
 * ties and parity corrections follow quant.c exactly. No secret indexing. */
#ifndef COPA_LATTICE_VECTORS_H
#define COPA_LATTICE_VECTORS_H
#include "copa_simd.h"
static inline copav_vec qv_div_positive(copav_vec x,unsigned d){
 copav_vec q=cv_high32(x,cv_set((uint32_t)(UINT64_C(0x100000000)/d))),r=cv_sub(x,cv_mul32(q,cv_set(d)));
 return cv_add(q,cv_and(cv_gt(r,cv_set(d-1)),cv_set(1)));
}
static void copa_quant_all(int32_t *out,int32_t *error,const int32_t *in){
 const unsigned mode=COPA_MODE,n=COPA_MODE==5?16:8,scale=COPA_MODE==5?2:COPA_MODE==2?5:COPA_MODE==6?2:1,spacing=2*scale,count=COPA_MODE==5?32:COPA_MODE==6?16:2;
 const int32_t *code=mode==5?&rm14[0][0]:mode==6?&rm13[0][0]:&rm03[0][0];
 for(unsigned block=0;block<64;block+=COPAV_W){copav_vec x[16],y[16],best[16],be[16];
  for(unsigned j=0;j<n;++j)x[j]=cv_load(in+(j/4)*256+(j%4)*64+block);
  if(mode==5){for(unsigned j=0;j<16;j+=4){y[pi_w[j]]=cv_sub(x[j],x[j+2]);y[pi_w[j+1]]=cv_sub(x[j+1],x[j+3]);y[pi_w[j+2]]=cv_add(x[j],x[j+2]);y[pi_w[j+3]]=cv_add(x[j+1],x[j+3]);}}
  else for(unsigned j=0;j<n;++j)y[pi(mode,j)]=x[j];
  copav_vec bd=cv_set(INT32_MAX);for(unsigned j=0;j<n;++j)best[j]=be[j]=cv_set(0);
  for(unsigned c=0;c<count;++c){copav_vec f[16],err[16],p[16],e[16],odd=cv_set(0),max=cv_set(-1),index=cv_set(0),dist=cv_set(0),less=cv_set(0),equal=cv_set(-1);
   for(unsigned j=0;j<n;++j){copav_vec v=cv_sub(y[j],cv_set(scale*code[c*n+j])),t=cv_add(v,cv_set(scale));
    if(mode==2) f[j]=qv_div_positive(t,10); else if(spacing==4)f[j]=cv_sra32(t,2);else f[j]=cv_sra32(t,1);
    err[j]=cv_sub(v,cv_mul32(f[j],cv_set(spacing)));odd=cv_add(odd,f[j]);copav_vec sign=cv_sra32(err[j],31),mag=cv_sub(cv_xor(err[j],sign),sign),take=cv_gt(mag,max);max=cv_select(mag,max,take);index=cv_select(cv_set(j),index,take);
   }
   for(unsigned j=0;j<n;++j){copav_vec correction=cv_or(cv_sra32(err[j],31),cv_set(1)),mask=mode==6?cv_set(0):cv_and(cv_sub(cv_set(0),cv_and(odd,cv_set(1))),cv_eq(index,cv_set(j)));
    f[j]=cv_add(f[j],cv_and(correction,mask));p[j]=cv_add(cv_mul32(f[j],cv_set(spacing)),cv_set(scale*code[c*n+j]));e[j]=cv_sub(p[j],y[j]);dist=cv_add(dist,cv_mul32(e[j],e[j]));less=cv_or(less,cv_and(equal,cv_gt(be[j],e[j])));equal=cv_and(equal,cv_eq(e[j],be[j]));}
   copav_vec take=cv_or(cv_gt(bd,dist),cv_and(cv_eq(bd,dist),less));bd=cv_select(dist,bd,take);
   for(unsigned j=0;j<n;++j){best[j]=cv_select(p[j],best[j],take);be[j]=cv_select(e[j],be[j],take);}
  }
  copav_vec result[16];if(mode==5){for(unsigned j=0;j<16;j+=4){copav_vec a=best[pi_w[j]],b=best[pi_w[j+1]],c=best[pi_w[j+2]],d=best[pi_w[j+3]];result[j]=cv_sra32(cv_add(a,c),1);result[j+1]=cv_sra32(cv_add(b,d),1);result[j+2]=cv_sra32(cv_sub(c,a),1);result[j+3]=cv_sra32(cv_sub(d,b),1);}}
  else for(unsigned j=0;j<n;++j)result[j]=best[pi(mode,j)];
  for(unsigned j=0;j<n;++j){unsigned offset=(j/4)*256+(j%4)*64+block;cv_store(out+offset,result[j]);cv_store(error+offset,cv_sub(result[j],x[j]));}
  if(mode==6)for(unsigned j=0;j<4;++j){unsigned offset=512+j*64+block;copav_vec v=cv_load(in+offset),p=cv_mul32(qv_div_positive(cv_add(v,cv_set(1)),3),cv_set(3));cv_store(out+offset,p);cv_store(error+offset,cv_sub(p,v));}
 }
}
#endif
