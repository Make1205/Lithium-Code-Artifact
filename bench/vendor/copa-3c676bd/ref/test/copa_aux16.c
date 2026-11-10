#include <stdint.h>
typedef struct {uint16_t p,qinv;const uint16_t *twist,*post,*forward,*inverse;} aux16_prime;
#include "copa_aux16_tables.h"
static inline uint16_t mul16(uint16_t a,uint16_t b,const aux16_prime *p){
 uint32_t z=(uint32_t)a*b;uint16_t lo=(uint16_t)z,t=(uint16_t)((uint32_t)lo*p->qinv);uint32_t tp=(uint32_t)t*p->p;
 uint16_t zh=(uint16_t)(z>>16),h=(uint16_t)(zh+(uint16_t)(tp>>16));unsigned carry=h<zh;
 uint16_t u=(uint16_t)(h+(lo!=0));carry|=u<h;return(uint16_t)(u-p->p*(carry|(u>=p->p)));
}
static inline uint16_t add16(uint16_t a,uint16_t b,uint16_t p){uint16_t sum=(uint16_t)(a+b);return(uint16_t)(sum-p*((sum<a)|(sum>=p)));}
static inline uint16_t sub16(uint16_t a,uint16_t b,uint16_t p){return(uint16_t)((uint16_t)(a-b)+p*(a<b));}
static void transform16(uint16_t a[256],const aux16_prime *p,unsigned inverse){
 for(unsigned j=0;j<256;j+=2){uint16_t x=a[j],y=a[j+1];a[j]=add16(x,y,p->p);a[j+1]=sub16(x,y,p->p);}
 const uint16_t *w=inverse?p->inverse:p->forward;
 for(unsigned len=4;len<=256;len<<=1){unsigned h=len/2;for(unsigned base=0;base<256;base+=len)for(unsigned j=0;j<h;++j){uint16_t x=a[base+j],y=mul16(a[base+j+h],w[h-1+j],p);a[base+j]=add16(x,y,p->p);a[base+j+h]=sub16(x,y,p->p);}}
}
static void forward16(uint16_t out[256],const int32_t in[256],const aux16_prime *p){
 uint16_t temp[256];for(unsigned j=0;j<256;++j){uint32_t mask=0u-(in[j]<0),a=((uint32_t)in[j]^mask)-mask;for(unsigned k=0;k<5;++k)a-=p->p*(a>=p->p);uint32_t neg=(p->p-a)*(a!=0);temp[j]=(uint16_t)((neg&mask)|(a&~mask));}
 for(unsigned j=0;j<256;++j)out[j]=mul16(temp[j],p->twist[j],p);
 for(unsigned j=0;j<256;++j){unsigned r=aux16_rev[j];if(j<r){uint16_t t=out[j];out[j]=out[r];out[r]=t;}}
 transform16(out,p,0);
}
void copa_aux16_mul(int64_t out[256],const int32_t a[256],const int32_t b[256]){
 uint16_t v[3][256];
 for(unsigned k=0;k<3;++k){const aux16_prime *p=&aux16_primes[k];uint16_t x[256],y[256];forward16(x,a,p);forward16(y,b,p);
  for(unsigned j=0;j<256;++j)v[k][aux16_rev[j]]=mul16(x[j],y[j],p);
  transform16(v[k],p,1);
  for(unsigned j=0;j<256;++j)v[k][j]=mul16(v[k][j],p->post[j],p);
 }
 const uint32_t p01=UINT32_C(3963743233);const uint64_t product=UINT64_C(235418601837569);
 for(unsigned j=0;j<256;++j){int32_t d=(int32_t)v[1][j]-v[0][j];uint32_t u=(uint32_t)(d+(d<0)*122882);u-=61441u*(u>=61441);
  uint32_t t1=(u*AUX16_INV01)%61441u,x01=v[0][j]+64513u*t1;
  int32_t dd=(int32_t)v[2][j]-(int32_t)(x01%59393u);uint32_t uu=(uint32_t)(dd+(dd<0)*59393);
  uint32_t t2=(uu*AUX16_INV012)%59393u;uint64_t x=x01+(uint64_t)p01*t2;out[j]=(int64_t)x-(int64_t)(product*(x>product/2));
 }
}
