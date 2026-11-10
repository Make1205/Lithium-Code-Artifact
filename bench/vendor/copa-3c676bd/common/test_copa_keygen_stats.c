#include "copa_api.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
int main(int argc,char **argv){unsigned n=argc>1?(unsigned)strtoul(argv[1],0,10):10000;uint8_t pk[COPA_PK_BYTES],sk[COPA_SK_BYTES];copa_keygen_stats sum={0};
 for(unsigned j=0;j<n;++j){copa_keygen_stats s;assert(!copa_keypair_bytes(pk,sk,&s));sum.attempts+=s.attempts;sum.noninvertible+=s.noninvertible;sum.norm_reject+=s.norm_reject;sum.lower_norm_reject+=s.lower_norm_reject;sum.packing_failure+=s.packing_failure;}
 assert(sum.attempts==n+sum.noninvertible+sum.norm_reject+sum.lower_norm_reject);
 printf("{\"mode\":%d,\"keys\":%u,\"attempts\":%llu,\"noninvertible\":%llu,\"upper_norm\":%llu,\"lower_norm\":%llu,\"packing_failure\":%llu,\"mean_attempts\":%.9f,\"unit_probability\":%.9f,\"upper_gate_acceptance_given_unit\":%.9f}\n",COPA_MODE,n,(unsigned long long)sum.attempts,(unsigned long long)sum.noninvertible,(unsigned long long)sum.norm_reject,(unsigned long long)sum.lower_norm_reject,(unsigned long long)sum.packing_failure,(double)sum.attempts/n,1-(double)sum.noninvertible/sum.attempts,(double)n/(sum.attempts-sum.noninvertible));return 0;}
