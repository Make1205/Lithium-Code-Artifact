#include "ntt_crt.h"
#include "poly_inv.h"
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
static uint64_t state = 1;
static uint32_t rng(void) { state ^= state << 13; state ^= state >> 7; state ^= state << 17; return (uint32_t)state; }
static void school(int64_t *out, const int32_t *a, const int32_t *b, size_t n)
{
    memset(out, 0, n*sizeof(*out));
    for (size_t i=0;i<n;++i) for(size_t j=0;j<n;++j) {
        int64_t v=(int64_t)a[i]*b[j];
        if(i+j<n)out[i+j]+=v;else out[i+j-n]-=v;
    }
}
/* Independent Euclidean polynomial gcd over each prime factor of q. */
static int gcd_unit(const int32_t *f, unsigned p)
{
    unsigned a[COPA_N+1]={0}, b[COPA_N+1]={0};
    a[0]=a[COPA_N]=1;
    for(size_t i=0;i<COPA_N;++i) b[i]=(unsigned)(((int64_t)f[i]%p+p)%p);
    int da=COPA_N,db=COPA_N-1;
    while(db>=0 && !b[db])--db;
    while(db>=0){
        unsigned inv=1;
        while((inv*b[db])%p!=1)++inv; /* public test only */
        while(da>=db){
            unsigned c=a[da]*inv%p;
            for(int i=0;i<=db;++i) a[i+da-db]=(a[i+da-db]+p-c*b[i]%p)%p;
            while(da>=0 && !a[da])--da;
        }
        for(size_t i=0;i<=COPA_N;++i){unsigned t=a[i];a[i]=b[i];b[i]=t;}
        int t=da;da=db;db=t;
    }
    return da==0;
}
int main(void)
{
    int32_t a[COPA_N],b[COPA_N],inv[COPA_N];
    int64_t got[COPA_N],want[COPA_N];
    unsigned qlist[]={7680,12584,6048,8400};
    size_t products=0,inversions=0,units=0;
    for(size_t n=1;n<=COPA_N;n*=2) for(int k=0;k<40;++k){
        for(size_t i=0;i<n;++i){a[i]=(int32_t)(rng()%524289)-262144;b[i]=(int32_t)(rng()%524289)-262144;}
        if(k==0)for(size_t i=0;i<n;++i)a[i]=b[i]=262144;
        school(want,a,b,n);copa_ntt_crt_mul(got,a,b,n);
        assert(!memcmp(want,got,n*sizeof(*got)));++products;
    }
    for(size_t j=0;j<4;++j){
        unsigned q=qlist[j];
        int64_t edge[]={INT64_MIN,INT64_MAX,-12584,-1,0,1,12584};
        for(size_t i=0;i<sizeof(edge)/sizeof(*edge);++i)
            assert(copa_mod_i64(edge[i],q)==(uint32_t)((edge[i]%q+q)%q));
        for(int k=0;k<100;++k){
            for(size_t i=0;i<COPA_N;++i)a[i]=(int32_t)(rng()%7)-3;
            if(k<3){memset(a,0,sizeof(a));a[k==1?1:0]=k==2?2:1;}
            int expected=gcd_unit(a,2);
            if(q==7680)expected &=gcd_unit(a,3)&gcd_unit(a,5);
            if(q==12584)expected &=gcd_unit(a,11)&gcd_unit(a,13);
            if(q==6048)expected &=gcd_unit(a,3)&gcd_unit(a,7);
            if(q==8400)expected &=gcd_unit(a,3)&gcd_unit(a,5)&gcd_unit(a,7);
            int ok=copa_poly_inv(inv,a,q);assert(ok==expected);
            school(want,a,inv,COPA_N);
            if(ok){++units;for(size_t i=0;i<COPA_N;++i) assert(copa_mod_i64(want[i],q)==(i==0));}
            else for(size_t i=0;i<COPA_N;++i) assert(inv[i]==0);
            memcpy(b,a,sizeof(b));assert(copa_poly_inv(b,b,q)==ok);assert(!memcmp(b,inv,sizeof(b)));
            ++inversions;
        }
    }
    printf("PASS Copa arithmetic: exact_products=%zu inversions=%zu units=%zu (independent prime-factor gcd)\n",products,inversions,units);
    return 0;
}
