#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "adapter.h"
#include "cpucycles.h"
#include "speed_print.h"
int bench_testing=1;
static uint64_t state;
void bench_seed(uint64_t s) { state=s+0x9e3779b97f4a7c15ULL; }
void bench_fill(void *out,size_t len) {
 uint8_t *p=out;
 while(len--) { state^=state>>12;state^=state<<25;state^=state>>27;*p++=(uint8_t)((state*2685821657736338717ULL)>>56); }
}
static void must(int ok,const char *s) {if(!ok){fprintf(stderr,"FAIL %s\n",s);exit(1);}}
static uint8_t pk[B_PK],sk[B_SK],sig[B_SIG],msg[32];
static size_t slen;
static void put(FILE*f,const void*p,size_t n){must(fwrite(p,1,n,f)==n,"write");}
static void get(FILE*f,void*p,size_t n){must(fread(p,1,n,f)==n,"read");}
static void sign_one(void){slen=B_SIG; must(b_sign(sig,&slen,msg,sizeof msg,sk)==0,"sign");must(slen<=B_SIG,"sig size");}
int main(int argc,char**argv) {
 must(argc>=3,"mode count [file]"); unsigned count=(unsigned)strtoul(argv[2],0,10);must(count>0,"count");
 printf("META %d %d %d %d %d\n",B_N,B_Q,B_PK,B_SK,B_SIG);
 bench_seed(1234567);
 if(!strcmp(argv[1],"generate")) {
  must(argc==4,"file");FILE*f=fopen(argv[3],"wb");must(f!=NULL,"open");
  must(!b_keygen(pk,sk),"keygen");put(f,pk,sizeof pk);put(f,sk,sizeof sk);
  unsigned long long sum=0;size_t lo=B_SIG,hi=0;
  for(unsigned i=0;i<count;i++){
   for(unsigned j=0;j<32;j++)msg[j]=(uint8_t)(i+j);
   sign_one();must(!b_verify(sig,slen,msg,32,pk),"self verify");
   size_t n=b_unpadded(sig,slen);must(n>0&&n<=slen,"unpadded bound");sum+=n;if(n<lo)lo=n;if(n>hi)hi=n;
   uint64_t l=slen;put(f,&l,8);put(f,sig,slen);
  }
  must(!fclose(f),"close");printf("LENGTH %llu %u %zu %zu\n",sum,count,lo,hi);puts("PASS generate");
 }else if(!strcmp(argv[1],"check")){
  must(argc==4,"file");FILE*f=fopen(argv[3],"rb");must(f!=NULL,"open");get(f,pk,sizeof pk);get(f,sk,sizeof sk);
  for(unsigned i=0;i<count;i++){
   for(unsigned j=0;j<32;j++)msg[j]=(uint8_t)(i+j);
   uint64_t l;get(f,&l,8);must(l<=B_SIG,"length");get(f,sig,(size_t)l);
   must(!b_verify(sig,(size_t)l,msg,32,pk),"cross verify");msg[0]^=1;must(b_verify(sig,(size_t)l,msg,32,pk)!=0,"tamper rejected");
  }must(fgetc(f)==EOF,"vector count");must(!fclose(f),"close");printf("PASS cross %u\n",count);
 }else if(!strcmp(argv[1],"time")){
  bench_testing=0;
  const unsigned warm=20;
  uint64_t *kg=calloc(count+1,sizeof *kg),*sg=calloc(count+1,sizeof *sg),*vf=calloc(count+1,sizeof *vf);
  must(kg&&sg&&vf,"alloc");
  /* Fresh byte keys and signatures on every round, matching local_check's
     sample population. Only the byte-API call is inside each interval. */
  for(unsigned i=0;i<count+warm;i++){
   uint64_t a=cpucycles();int rc=b_keygen(pk,sk);uint64_t dk=cpucycles()-a;must(!rc,"timed keygen");
   slen=B_SIG;a=cpucycles();rc=b_sign(sig,&slen,msg,32,sk);uint64_t ds=cpucycles()-a;must(!rc,"timed sign");
   a=cpucycles();rc=b_verify(sig,slen,msg,32,pk);uint64_t dv=cpucycles()-a;must(!rc,"timed verify");
   if(i>=warm){unsigned j=i-warm;kg[j+1]=kg[j]+dk;sg[j+1]=sg[j]+ds;vf[j+1]=vf[j]+dv;}
  }
  print_results("KeyGen",kg,count+1);print_results("Sign",sg,count+1);print_results("Verify",vf,count+1);
#ifdef B_EXPANDED
  must(!b_expand(sk),"expanded key preparation");
  uint64_t *ex=calloc(count+1,sizeof *ex);must(ex!=NULL,"expanded allocation");
  for(unsigned i=0;i<count+warm;i++){
   slen=B_SIG;uint64_t a=cpucycles();int rc=b_expanded_sign(sig,&slen,msg,32);uint64_t dt=cpucycles()-a;must(!rc,"expanded signing");
   if(i>=warm){unsigned j=i-warm;ex[j+1]=ex[j]+dt;}
  }
  must(!b_verify(sig,slen,msg,32,pk),"expanded verification");print_results("ExpandedSign",ex,count+1);free(ex);
#endif
  free(kg);free(sg);free(vf);

 }else must(0,"mode");return 0;
}
