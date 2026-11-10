#include "copa_prof.h"
#include "pack_pk.h"
#include "quant.h"
#include <string.h>
#include <assert.h>
size_t copa_pk_bytes(unsigned m){return m==2?648:m==3?840:m==5?1480:1160;}
static unsigned lowbits(unsigned p){unsigned k=0;while(!(p&1)){++k;p>>=1;}return k;}
typedef unsigned __int128 copa_pk_wide;
static unsigned bitlength(copa_pk_wide v){unsigned k=0;while(v){++k;v>>=1;}return k;}
/* The original little-endian low-part/mixed-radix wire layout is retained.
 * v2 needs up to 121 mixed bits; bounded 56-bit windows avoid wide shifts
 * or reads past the eight-byte padding. Loop counts depend only on mode. */
static uint64_t load64(const uint8_t *p){uint64_t v=0;for(unsigned j=0;j<8;++j)v|=(uint64_t)p[j]<<(8*j);return v;}
static void store64(uint8_t *p,uint64_t v){for(unsigned j=0;j<8;++j)p[j]=(uint8_t)(v>>(8*j));}
static void put(uint8_t *buf,size_t *pos,copa_pk_wide v,unsigned bits){
 while(bits){unsigned n=bits>56?56:bits;uint8_t *p=buf+*pos/8;uint64_t mask=(UINT64_C(1)<<n)-1;
  store64(p,load64(p)|(((uint64_t)v&mask)<<(*pos%8)));*pos+=n;v>>=n;bits-=n;}
}
static copa_pk_wide get(const uint8_t *buf,size_t *pos,unsigned bits){
 copa_pk_wide v=0;unsigned shift=0;
 while(bits){unsigned n=bits>56?56:bits;uint64_t x=(load64(buf+*pos/8)>>(*pos%8))&((UINT64_C(1)<<n)-1);
  v|=(copa_pk_wide)x<<shift;*pos+=n;shift+=n;bits-=n;}return v;
}
/* Only public key labels are divided; these constant radices also avoid
 * a general wide divide for the four installed parameter sets. */
/* Split the public mixed integer in groups of four radix digits.
 * Only n/4 wide quotients remain; each group fits in uint64_t. */
/* floor(v/d) from floor(2^128/d): the estimate is at most one low.
 * The 128x128 high product is four exact 64x64 products; r<2d. */
static copa_pk_wide pk_high_product(copa_pk_wide x,copa_pk_wide y){
 uint64_t xl=(uint64_t)x,xh=(uint64_t)(x>>64),yl=(uint64_t)y,yh=(uint64_t)(y>>64);
 copa_pk_wide a=(copa_pk_wide)xl*yl,b=(copa_pk_wide)xh*yl,c=(copa_pk_wide)xl*yh;
 copa_pk_wide mid=(a>>64)+(uint64_t)b+(uint64_t)c;
 return (copa_pk_wide)xh*yh+(b>>64)+(c>>64)+(mid>>64);
}
static copa_pk_wide divide_four(copa_pk_wide v,unsigned odd){
 uint64_t d,hi,lo;switch(odd){
 case 3:d=UINT64_C(81);hi=UINT64_C(227737581156908044);lo=UINT64_C(11842354220159218321);break;
 case 1573:d=UINT64_C(6122304000241);hi=UINT64_C(3013039);lo=UINT64_C(10097079412225518701);break;
 case 189:d=UINT64_C(1275989841);hi=UINT64_C(14456811081);lo=UINT64_C(15863798462901322420);break;
 case 525:d=UINT64_C(75969140625);hi=UINT64_C(242818912);lo=UINT64_C(381723863037492715);break;
 default:assert(odd==1);return v;}
 copa_pk_wide q=pk_high_product(v,((copa_pk_wide)hi<<64)|lo);uint64_t rem=(uint64_t)(v-q*d);return q+(rem>=d);
}
static uint64_t divide_digit(uint64_t v,unsigned odd){
 switch(odd){case 3:return v/3;case 1573:return v/1573;case 189:return v/189;case 525:return v/525;default:assert(odd==1);return v;}
}
static void unpack_digits(uint32_t label[16],const unsigned shift[16],unsigned n,unsigned odd,copa_pk_wide mixed){
 uint64_t radix4=(uint64_t)odd*odd*odd*odd;
 for(unsigned base=0;base<n;base+=4){
  copa_pk_wide quotient=divide_four(mixed,odd);uint64_t group=(uint64_t)(mixed-quotient*radix4);mixed=quotient;
  for(unsigned j=base;j<base+4;++j){uint64_t q=divide_digit(group,odd);label[j]|=(uint32_t)(group-q*odd)<<shift[j];group=q;}
 }
}
int copa_pack_pk(uint8_t *pk,size_t len,const uint8_t rho[32],const int32_t *b,unsigned mode)
{
    if(len!=copa_pk_bytes(mode))return 0;
    unsigned n=copa_label_width(mode),q=copa_modulus(mode);size_t pos=256;
    uint8_t padded[1480+8]={0};uint8_t *wire=pk;pk=padded;memcpy(pk,rho,32);
    unsigned shift[16],low_total=0;copa_pk_wide powers[16],power=1;
    for(unsigned i=0;i<n;++i){unsigned p=copa_label_radix(mode,i);shift[i]=lowbits(p);low_total+=shift[i];powers[i]=power;power*=p>>shift[i];}
    unsigned mixed_bits=bitlength(power-1);assert(low_total<=64);
    for(unsigned k=0;k<64;++k){
        int32_t x[16];uint32_t label[16];copa_pk_wide mixed=0;uint64_t low=0;unsigned offset=0;
        for(unsigned r=0;r<n;++r)x[r]=b[(r/4)*256+k+(r%4)*64];
        if(!copa_lab(label,x,mode))return 0;
        for(unsigned i=0;i<n;++i){low|=(uint64_t)(label[i]&((1u<<shift[i])-1))<<offset;offset+=shift[i];mixed+=powers[i]*(label[i]>>shift[i]);}
        put(pk,&pos,low,low_total);put(pk,&pos,mixed,mixed_bits);
    }
    if(mode==6)for(unsigned i=0;i<256;++i){
        int32_t v=b[512+i];if(v<0||(unsigned)v>=q||v%3)return 0;
        put(pk,&pos,(uint32_t)v/3,bitlength(q/3-1));
    }
    assert(pos==8*len);memcpy(wire,pk,len);return 1;
}
#ifndef COPA_VEC_LATTICE
int copa_unpack_pk(uint8_t rho[32],int32_t *b,const uint8_t *pk,size_t len,unsigned mode)
{
    if(len!=copa_pk_bytes(mode))return 0;
    uint8_t padded[1480+8];memcpy(padded,pk,len);memset(padded+len,0,8);pk=padded;
    unsigned n=copa_label_width(mode),shift[16],odd[16],low_total=0;copa_pk_wide power=1;size_t pos=256;memcpy(rho,pk,32);
    for(unsigned i=0;i<n;++i){unsigned p=copa_label_radix(mode,i);shift[i]=lowbits(p);odd[i]=p>>shift[i];power*=odd[i];low_total+=shift[i];}
    unsigned mixed_bits=bitlength(power-1);
    for(unsigned k=0;k<64;++k){
        uint32_t label[16];
        uint64_t low=(uint64_t)get(pk,&pos,low_total);
        for(unsigned i=0;i<n;++i){label[i]=(uint32_t)low&((1u<<shift[i])-1);low>>=shift[i];}
        copa_pk_wide mixed=get(pk,&pos,mixed_bits);if(mixed>=power)return 0;
        unpack_digits(label,shift,n,odd[0],mixed);
        int32_t x[16];if(!copa_unlab(x,label,mode))return 0;
        for(unsigned r=0;r<n;++r)b[(r/4)*256+k+(r%4)*64]=x[r];
    }
    if(mode==6)for(unsigned i=0;i<256;++i){unsigned v=(unsigned)get(pk,&pos,bitlength(copa_modulus(mode)/3-1));if(v>=copa_modulus(mode)/3)return 0;b[512+i]=3*(int32_t)v;}
    assert(pos==8*len);return 1;
}

#else
#include "../common/copa_unpack_vectors.h"
#endif
