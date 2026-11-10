#define _GNU_SOURCE
#include <x86intrin.h>
static inline void copa_test_fence(void){_mm_mfence();_mm_lfence();}
#define _mm_lfence() copa_test_fence()
#include "../ref/test/copa_dudect.c"
