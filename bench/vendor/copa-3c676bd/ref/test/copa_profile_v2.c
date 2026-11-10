#ifdef NDEBUG
#error Benchmarks require active correctness assertions
#endif
#include "copa_api.h"
#include "copa_prof.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
uint64_t copa_profile_cycles[CP_COUNT],copa_profile_calls[CP_COUNT];
static const char *names[]={"parse_sk_cache","ExpandY_including_XOF","Ay_NTT_CRT","HighBits_challenge","IRS","Bs","MakeHint_z2_Bv","rANS","pk_parse_Unlab","ExpandA","ExpandD","matrix_prepare_NTT","sig_decode","verify_NTT_CRT","challenge_hash","verify_reconstruct","ExpandS","inverse","quantization","packing","NTT_forward","NTT_pointwise","NTT_inverse","sk_decode_norm","sk_hash","secret_NTT_correlation","IRS_prepare","mu_rhoprime_setup","bprime_rebuild","SHAKE_Y","IRS_projection","IRS_transition","IRS_rotation","Unlab","bprime_NTT","A_NTT","verify_norm","KG_hash","KG_norm","KG_compute_b","KG_matrix"};
int main(int argc,char **argv){assert(argc==3);unsigned count=(unsigned)strtoul(argv[2],0,10);uint8_t pk[COPA_PK_BYTES],skbytes[COPA_SK_BYTES],sig[COPA_SIG_BYTES],m[59]={0};size_t len;copa_secret_key sk;
 assert(!copa_keygen(pk,&sk,0));assert(!copa_pack_sk(skbytes,&sk));assert(!copa_sign(sig,m,sizeof(m),&sk,0));memset(copa_profile_cycles,0,sizeof(copa_profile_cycles));memset(copa_profile_calls,0,sizeof(copa_profile_calls));uint64_t tries=0,steps=0,kgtries=0,noninv=0,norm=0;uint64_t start=copa_profile_tick();
 for(unsigned i=0;i<count;++i){if(!strcmp(argv[1],"keygen")){copa_keygen_stats st;assert(!copa_keypair_bytes(pk,skbytes,&st));kgtries+=st.attempts;noninv+=st.noninvertible;norm+=st.norm_reject;}else if(!strcmp(argv[1],"sign")){copa_sign_stats st;assert(!copa_sign_signature(sig,&len,m,sizeof(m),skbytes,sizeof(skbytes),&st));tries+=st.attempts;steps+=st.irs_steps;}else{assert(!strcmp(argv[1],"verify"));assert(!copa_verify(sig,sizeof(sig),m,sizeof(m),pk,sizeof(pk)));}}
 uint64_t elapsed=copa_profile_tick()-start;
 uint64_t fixed=copa_profile_cycles[CP_PARSE_SK]+copa_profile_cycles[CP_SIGN_SETUP],iteration=0;
 const unsigned ids[]={CP_Y,CP_AY,CP_HIGH_CH,CP_IRS,CP_BS,CP_HINT_BV,CP_RANS};
 for(unsigned j=0;j<sizeof(ids)/sizeof(*ids);++j)iteration+=copa_profile_cycles[ids[j]];
 printf("{\"fixed_cycles\":%llu,\"iteration_cycles\":%llu,\"fixed_fraction\":%.9f,\"actual_mean_iteration_cycles\":%.9f,",(unsigned long long)fixed,(unsigned long long)iteration,(double)fixed/elapsed,tries?(double)iteration/tries:0.0);
 printf("\"mode\":%d,\"op\":\"%s\",\"samples\":%u,\"total_cycles\":%llu,\"iterations\":%llu,\"irs_steps\":%llu,\"keygen_attempts\":%llu,\"noninvertible\":%llu,\"norm_reject\":%llu,\"groups\":{",COPA_MODE,argv[1],count,(unsigned long long)elapsed,(unsigned long long)tries,(unsigned long long)steps,(unsigned long long)kgtries,(unsigned long long)noninv,(unsigned long long)norm);
 for(unsigned j=0;j<CP_COUNT;++j)printf("%s\"%s\":{\"cycles\":%llu,\"calls\":%llu}",j?",":"",names[j],(unsigned long long)copa_profile_cycles[j],(unsigned long long)copa_profile_calls[j]);puts("}}");return 0;}
