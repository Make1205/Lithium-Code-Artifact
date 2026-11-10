#include <stddef.h>
#include <stdint.h>
extern int bench_testing;
void bench_fill(void *,size_t);
void bench_native_randombytes(uint8_t *,size_t);
void randombytes(uint8_t *p,size_t n){if(bench_testing)bench_fill(p,n);else bench_native_randombytes(p,n);}
