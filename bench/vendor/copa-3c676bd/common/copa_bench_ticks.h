#ifndef COPA_BENCH_TICKS_H
#define COPA_BENCH_TICKS_H
#include <stdint.h>
#include <x86intrin.h>
static inline uint64_t copa_bench_ticks(void){unsigned aux;_mm_lfence();uint64_t t=__rdtscp(&aux);_mm_lfence();return t;}
#endif
