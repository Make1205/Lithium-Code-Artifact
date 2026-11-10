#define _GNU_SOURCE
/* SDE traces only the fixed primitive between SSC marks; print public
 * address anchors so runs with different ASLR stack bases can be compared.
 */
#include <assert.h>
#include <dlfcn.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "copa_sampler.h"
#include "copa_scalar64_tables.h"
#include "copa_sign.h"
#include "ntt_crt.h"
static copa_irs64_context ctx;
static int32_t a[256], b[256], sample[(1 + COPA_ELL) * COPA_N];
static int64_t product[256];
static uint8_t seed[32];
static uint64_t random_bits[6] = {1, 2, 3, 4, 5, 6}, nums[4];
static copa_fp63 out[4];
static volatile unsigned sink;
static volatile int64_t trace_dot;
#define START                                                      \
    __asm__ volatile("movl $0x111, %%ebx; .byte 0x64,0x67,0x90" :: \
                         : "ebx", "memory")
#define STOP                                                       \
    __asm__ volatile("movl $0x222, %%ebx; .byte 0x64,0x67,0x90" :: \
                         : "ebx", "memory")
int main(int argc, char **argv)
{
    assert(argc == 3);
    unsigned cls = (unsigned)strtoul(argv[2], 0, 10);
    assert(cls < 2);
    for (unsigned i = 0; i < 256; ++i) {
        a[i] = cls ? (int32_t)(i % 7) - 3 : (int32_t)(i == 0);
        b[i] = cls ? (int32_t)(i % 5) - 2 : (int32_t)(i == 0);
    }
    for (unsigned i = 0; i < 32; ++i)
        seed[i] = cls ? (uint8_t)(i * 37 + 11) : 0;
    for (unsigned i = 0; i < 4; ++i)
        nums[i] = 10000 * (1 + cls * (100 + i));
    copa_irs64_prepare(&ctx, cls ? 64 : 1);
    trace_dot = cls ? -64 : 0;
    volatile unsigned char aligned_stack[64] __attribute__((aligned(64)));
    aligned_stack[0] = 0;
    __asm__ volatile(""::"m"(aligned_stack):"memory");
    Dl_info info;
    assert(dladdr(dlsym(RTLD_DEFAULT,"memcpy"),&info));
    printf("libc_anchor=%" PRIxPTR "\n",(uintptr_t)info.dli_fbase);
    uintptr_t sp;
    __asm__ volatile("mov %%rsp,%0" : "=r"(sp));
    printf("stack_anchor=%" PRIxPTR "\n", sp);
    fflush(stdout);
    if (!strcmp(argv[1], "ntt")) {
        START;
        copa_ntt_crt_mul(product, a, b, 256);
        STOP;
    } else if (!strcmp(argv[1], "noise")) {
        START;
        copa_expand_s((int32_t(*)[COPA_N])sample, seed);
        STOP;
    } else if (!strcmp(argv[1], "irs")) {
        unsigned minus;
        START;
        sink = copa_irs64_step(&minus, trace_dot, &ctx, random_bits);
        STOP;
    } else if (!strcmp(argv[1], "exp")) {
        START;
        copa64_exp_ratio4(out, nums, &copa64_r2);
        STOP;
    } else
        return 2;
    return 0;
}
