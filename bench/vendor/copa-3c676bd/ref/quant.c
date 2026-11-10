#include "copa_prof.h"
#include "quant.h"
#include "ntt_crt.h"
#include "copa_scalar64.h"
#include "fips202.h"
#include "copa_xofbuf.h"
#ifdef COPA_PARALLEL_XOF
#include "copa_xof8.h"
#endif
#include <assert.h>
#include <string.h>
#include "copa_lattice_tables.h"
static const unsigned pi_e[8]={1,2,4,7,0,3,5,6};
static const unsigned pi_w[16]={1,2,4,7,8,11,13,14,0,3,5,6,9,10,12,15};
unsigned copa_modulus(unsigned m){assert(m==2||m==3||m==5||m==6);return m==2?7680:m==3?12584:m==5?6048:8400;}
unsigned copa_rows(unsigned m){return m==5?4:m==6?3:2;}
unsigned copa_quant_width(unsigned m){return m==5?16:m==6?12:8;}
unsigned copa_label_width(unsigned m){return m==5?16:8;}
static unsigned sigma(unsigned m){return m==2?5:m==6?2:1;}
static unsigned pi(unsigned m,unsigned i){return m==3?i:m==5?pi_w[i]:m==2?i:pi_e[i];}
static unsigned diagonal(unsigned m,unsigned i){return m==5?diag_phiw[i]:m==6?diag_e[i]:diag_phie[i];}
static unsigned basis(unsigned m,unsigned i,unsigned j){return m==5?basis_phiw[i][j]:m==6?basis_e[i][j]:basis_phie[i][j];}
unsigned copa_label_radix(unsigned m,unsigned i){return copa_modulus(m)/(sigma(m)*diagonal(m,i));}
/* The only division is 2^64 / public d; the dividend of that instruction
 * is fixed, not the secret x. Remaining quotient arithmetic has no division. */
static int64_t floor_div(int64_t x, unsigned d)
{
    if(d==1)return x; /* public divisor */
    uint64_t neg=(uint64_t)0-((uint64_t)x>>63), a=((uint64_t)x^neg)-neg;
    uint64_t rec=UINT64_MAX/d+(UINT64_MAX%d==d-1),lo;
    uint64_t v=copa64_product(a,rec,&lo), rem=a-v*d;
    v+=(rem>=d); rem=a-v*d;
    uint64_t minus=(v+(rem!=0));
    return (int64_t)((v&~neg)|((0-minus)&neg));
}
static int32_t select32(int32_t a,int32_t b,unsigned choose_a)
{uint32_t mask=(uint32_t)copa64_ct_mask(choose_a);return (int32_t)(((uint32_t)a&mask)|((uint32_t)b&~mask));}
static void inner(int32_t *out,const int32_t *y,unsigned n,unsigned spacing,unsigned parity)
{
    int32_t f[16],err[16];uint32_t odd=0;uint64_t maxerr=0;unsigned index=0;
    for(unsigned j=0;j<n;++j){
        f[j]=(int32_t)floor_div((int64_t)y[j]+spacing/2,spacing);
        err[j]=y[j]-(int32_t)spacing*f[j];odd+=(uint32_t)f[j];
        int64_t e=err[j];uint64_t sign=0u-(uint64_t)(e<0), mag=((uint64_t)e^sign)-sign;
        unsigned take=mag>maxerr;uint64_t mask=0u-(uint64_t)take;
        maxerr=(mag&mask)|(maxerr&~mask);index=(unsigned)select32((int32_t)j,(int32_t)index,take);
    }
    for(unsigned j=0;j<n;++j){
        int32_t correction=1-2*(int32_t)(err[j]<0);
        f[j]+=correction*(int32_t)(parity&(odd&1u)&(j==index));
        out[j]=(int32_t)spacing*f[j];
    }
}
static void cosets(int32_t *out,const int32_t *x,unsigned n,unsigned scale,
                   unsigned parity,const int32_t *code,unsigned count)
{
    int32_t best[16]={0},be[16]={0};int64_t bd=INT64_MAX;
    for(unsigned c=0;c<count;++c){
        int32_t y[16],p[16],e[16];int64_t dist=0;
        for(unsigned j=0;j<n;++j)y[j]=x[j]-(int32_t)scale*code[c*n+j];
        inner(p,y,n,2*scale,parity);
        unsigned less=0,equal=1;
        for(unsigned j=0;j<n;++j){
            p[j]+=(int32_t)scale*code[c*n+j];e[j]=p[j]-x[j];dist+=(int64_t)e[j]*e[j];
            less|=equal&(e[j]<be[j]);equal&=(e[j]==be[j]);
        }
        unsigned take=(dist<bd)|((dist==bd)&less);uint64_t mask=0u-(uint64_t)take;
        bd=(int64_t)(((uint64_t)dist&mask)|((uint64_t)bd&~mask));
        for(unsigned j=0;j<n;++j){best[j]=select32(p[j],best[j],take);be[j]=select32(e[j],be[j],take);}
    }
    memcpy(out,best,n*sizeof(*out));
}
void copa_quant_block(int32_t *out,const int32_t *x,unsigned mode)
{
    unsigned n=copa_label_width(mode),s=sigma(mode);int32_t y[16],p[16],ring[16];
    if(mode==5){
        for(unsigned k=0;k<16;k+=4){ring[k]=x[k]-x[k+2];ring[k+1]=x[k+1]-x[k+3];ring[k+2]=x[k]+x[k+2];ring[k+3]=x[k+1]+x[k+3];}
        for(unsigned j=0;j<n;++j)y[pi_w[j]]=ring[j];
        cosets(p,y,16,2,1,&rm14[0][0],32);
        for(unsigned j=0;j<n;++j)ring[j]=p[pi_w[j]];
        for(unsigned k=0;k<16;k+=4){out[k]=(ring[k]+ring[k+2])/2;out[k+1]=(ring[k+1]+ring[k+3])/2;out[k+2]=(-ring[k]+ring[k+2])/2;out[k+3]=(-ring[k+1]+ring[k+3])/2;}
    }else{
        for(unsigned j=0;j<n;++j)y[pi(mode,j)]=x[j];
        cosets(p,y,n,s,mode!=6,mode==6?&rm13[0][0]:&rm03[0][0],mode==6?16:2);
        for(unsigned j=0;j<n;++j)out[j]=p[pi(mode,j)];
    }
    if(mode==6)for(unsigned j=8;j<12;++j)out[j]=(int32_t)(3*floor_div((int64_t)x[j]+1,3));
}
int copa_lab(uint32_t *label,const int32_t *x,unsigned mode)
{
    unsigned n=copa_label_width(mode),s=sigma(mode);int64_t rm[16]={0};uint32_t bad=0;
    for(unsigned j=0;j<n;++j)rm[pi(mode,j)]=x[j];
    for(unsigned i=0;i<n;++i){
        int64_t v=rm[i];for(unsigned j=0;j<i;++j)v-=(int64_t)s*basis(mode,i,j)*label[j];
        unsigned d=s*diagonal(mode,i);bad|=copa_mod_i64(v,d==1?2:d)*(uint32_t)(d!=1);
        label[i]=copa_mod_i64(floor_div(v,d),copa_label_radix(mode,i));
    }
    return bad==0;
}
int copa_unlab(int32_t *x,const uint32_t *label,unsigned mode)
{
    CP_BEGIN(punlab);unsigned n=copa_label_width(mode),s=sigma(mode),q=copa_modulus(mode);uint32_t bad=0,reciprocal=(uint32_t)(UINT64_C(0x100000000)/q);int32_t rm[16];
    for(unsigned i=0;i<n;++i){
        int64_t v=0;bad|=(label[i]>=copa_label_radix(mode,i));
        for(unsigned j=0;j<=i;++j)v+=(int64_t)s*basis(mode,i,j)*label[j];
        /* Valid labels give 0 <= v < n*q. Invalid labels are rejected
         * by bad; unsigned narrowing remains defined for those inputs. */
        uint32_t value=(uint32_t)v,quotient=(uint32_t)(((uint64_t)value*reciprocal)>>32);
        value-=quotient*q;rm[i]=(int32_t)(value-q*(value>=q));
    }
    for(unsigned i=0;i<n;++i)x[i]=rm[pi(mode,i)];
    CP_END(CP_UNLAB,punlab);return bad==0;
}
#ifdef COPA_VEC_LATTICE
#include "../common/copa_expand_d_vectors.h"
#endif
void copa_expand_d(int32_t *out,const uint8_t rho[32],unsigned mode)
{
#ifdef COPA_VEC_LATTICE
    copa_expand_d_vectors(out,rho,mode);
#else
    unsigned n=copa_label_width(mode),total=copa_rows(mode)*256,per=total/8,s=sigma(mode);
    unsigned inverse[16];for(unsigned r=0;r<n;++r)inverse[pi(mode,r)]=r;
#ifdef COPA_PARALLEL_XOF
    copa_xof8 engine;uint8_t blocks[8][136];unsigned counts[8]={0},done=0;
    copa_xof8_seed_init(&engine,rho,0x14);
    while(done<total){copa_xof8_block(&engine,blocks);for(unsigned lane=0;lane<8;++lane)for(unsigned j=0;j<136&&counts[lane]<per;j+=2){
     unsigned idx=lane*per+counts[lane],range=idx<64*n?s*diagonal(mode,idx%n):3,v=blocks[lane][j]+256u*blocks[lane][j+1];
     if(v>=65536/range*range)continue;
     v%=range;++counts[lane];++done;
     if(idx<64*n){unsigned block=idx/n,r=inverse[idx%n];out[(r/4)*256+block+(r%4)*64]=(int32_t)v;}else out[512+idx-64*n]=(int32_t)v;
    }}
#else
    for(unsigned lane=0;lane<8;++lane){
        uint8_t ext[35];copa_xofbuf st;memcpy(ext,rho,32);ext[32]=0x14;ext[33]=(uint8_t)lane;ext[34]=0;
        copa_xof_init(&st,ext,sizeof(ext));
        for(unsigned j=0;j<per;++j){
            unsigned idx=lane*per+j,range=idx<64*n?s*diagonal(mode,idx%n):3;
            unsigned threshold=65536/range*range,v;
            do{v=copa_xof_u16(&st);}while(v>=threshold);
            v%=range; /* public seed-derived dithering */
            if(idx<64*n){unsigned block=idx/n,rmidx=idx%n;unsigned r=inverse[rmidx];out[(r/4)*256+block+(r%4)*64]=(int32_t)v;}
            else out[2*256+idx-64*n]=(int32_t)v;
        }
    }
#endif
#endif
}

#ifdef COPA_VEC_LATTICE
#include "../common/copa_lattice_vectors.h"
#endif
