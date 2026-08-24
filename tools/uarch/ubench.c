/*
 * ubench.c - on-device instruction latency/throughput microbench for the
 * one AVX-512-capable benchmark system.
 *
 * the test environment had no usable PMU (perf/rdpmc) and no exposed
 * P-state, so we time with rdtsc (invariant TSC @ 3.6 GHz reference) and
 * NORMALISE every result to a vpaddq-latency anchor (= exactly 1 core
 * cycle, measured in the same 512-bit AVX-512-license frequency domain).
 * The normalised column is thus frequency-independent *core cycles*, which
 * is what determines IPC and the Turbo-off / Zen5 picture the design
 * targets.
 *
 *   latency  : a dependent chain (op_k reads op_{k-1}) -> core-cycle
 * latency throughput: 16 independent accumulators (zmm0..15) -> the loop
 * is port/issue bound, so cyc/op = 1/sustained-IPC.
 *
 * Build:  gcc -O2 -march=native tools/uarch/ubench.c -o /tmp/ubench
 * Run  :  taskset -c 3 /tmp/ubench
 */
#include <stdint.h>
#include <stdio.h>

/* ---- unrolled-body repetition ---- */
#define R2(x) x x
#define R4(x) R2(x) R2(x)
#define R8(x) R4(x) R4(x)
#define R16(x) R8(x) R8(x)
#define R32(x) R16(x) R16(x)

/* 16 copies of instruction-prefix I, one per dest zmm0..15 (throughput).
 */
#define A16(I)                                                      \
    I "%%zmm0\n\t" I "%%zmm1\n\t" I "%%zmm2\n\t" I "%%zmm3\n\t" I   \
      "%%zmm4\n\t" I "%%zmm5\n\t" I "%%zmm6\n\t" I "%%zmm7\n\t" I   \
      "%%zmm8\n\t" I "%%zmm9\n\t" I "%%zmm10\n\t" I "%%zmm11\n\t" I \
      "%%zmm12\n\t" I "%%zmm13\n\t" I "%%zmm14\n\t" I "%%zmm15\n\t"

/* interleave two instruction-prefixes A,B over 8 dest pairs (16 ops):
 * tells us whether A and B share an execution port (mix tput -> 1/cyc) or
 * use distinct ports (mix tput -> 2/cyc). */
#define MIX16(A, B)                                                 \
    A "%%zmm0\n\t" B "%%zmm8\n\t" A "%%zmm1\n\t" B "%%zmm9\n\t" A   \
      "%%zmm2\n\t" B "%%zmm10\n\t" A "%%zmm3\n\t" B "%%zmm11\n\t" A \
      "%%zmm4\n\t" B "%%zmm12\n\t" A "%%zmm5\n\t" B "%%zmm13\n\t" A \
      "%%zmm6\n\t" B "%%zmm14\n\t" A "%%zmm7\n\t" B "%%zmm15\n\t"

/* rdtscp + lfence timestamp into %[c0]/%[c1]. */
#define TS0                                                          \
    "lfence\n\trdtscp\n\tshlq $32,%%rdx\n\torq %%rdx,%%rax\n\tmovq " \
    "%%rax,%[c0]\n\tlfence\n\t"
#define TS1                                                          \
    "lfence\n\trdtscp\n\tshlq $32,%%rdx\n\torq %%rdx,%%rax\n\tmovq " \
    "%%rax,%[c1]\n\tlfence\n\t"

/* zmm16/17 = 52-bit all-ones, zmm18 = 0 (perm idx / shift cnt), zmm0 seed.
 */
#define SET512                                                      \
    "vpternlogq $0xff,%%zmm16,%%zmm16,%%zmm16\n\t vpsrlq "          \
    "$12,%%zmm16,%%zmm16\n\t vmovdqa64 %%zmm16,%%zmm17\n\t vpxorq " \
    "%%zmm18,%%zmm18,%%zmm18\n\t vmovdqa64 %%zmm16,%%zmm0\n\t"

#define CL512                                                            \
    "rax", "rdx", "rcx", "zmm0", "zmm1", "zmm2", "zmm3", "zmm4", "zmm5", \
        "zmm6", "zmm7", "zmm8", "zmm9", "zmm10", "zmm11", "zmm12",       \
        "zmm13", "zmm14", "zmm15", "zmm16", "zmm17", "zmm18", "cc",      \
        "memory"

#define LAT512(NAME, INSN)                                               \
    static uint64_t NAME(uint64_t it)                                    \
    {                                                                    \
        uint64_t c0 = 0, c1 = 0;                                         \
        __asm__ volatile(                                                \
            SET512 TS0 "1:\n\t" R32(INSN) "decq %[it]\n\tjnz 1b\n\t" TS1 \
            : [c0] "=&r"(c0), [c1] "=&r"(c1), [it] "+r"(it)              \
            :                                                            \
            : CL512);                                                    \
        return c1 - c0;                                                  \
    }
#define TPUT512(NAME, PREFIX)                                            \
    static uint64_t NAME(uint64_t it)                                    \
    {                                                                    \
        uint64_t c0 = 0, c1 = 0;                                         \
        __asm__ volatile(SET512 TS0 "1:\n\t" R4(                         \
                             A16(PREFIX)) "decq %[it]\n\tjnz 1b\n\t" TS1 \
                         : [c0] "=&r"(c0), [c1] "=&r"(c1), [it] "+r"(it) \
                         :                                               \
                         : CL512);                                       \
        return c1 - c0;                                                  \
    }

/* ---- 512-bit ops: latency (chain on zmm0) + throughput (16 indep) ----
 */
LAT512(lat_vpaddq, "vpaddq %%zmm16,%%zmm0,%%zmm0\n\t") /* ANCHOR = 1cyc */
TPUT512(tput_vpaddq, "vpaddq %%zmm16,%%zmm17,")

LAT512(lat_madd52lo, "vpmadd52luq %%zmm16,%%zmm17,%%zmm0\n\t")
TPUT512(tput_madd52lo, "vpmadd52luq %%zmm16,%%zmm17,")
LAT512(lat_madd52hi, "vpmadd52huq %%zmm16,%%zmm17,%%zmm0\n\t")
TPUT512(tput_madd52hi, "vpmadd52huq %%zmm16,%%zmm17,")

LAT512(lat_permt2q, "vpermt2q %%zmm17,%%zmm18,%%zmm0\n\t")
TPUT512(tput_permt2q, "vpermt2q %%zmm17,%%zmm18,")
TPUT512(tput_permq, "vpermq $0x1b,%%zmm17,") /* p5 shuffle proxy */

LAT512(lat_sllvq, "vpsllvq %%zmm18,%%zmm0,%%zmm0\n\t")
TPUT512(tput_sllvq, "vpsllvq %%zmm16,%%zmm17,")
LAT512(lat_sllq, "vpsllq $1,%%zmm0,%%zmm0\n\t")
TPUT512(tput_sllq, "vpsllq $1,%%zmm17,")

LAT512(lat_rolq, "vprolq $7,%%zmm0,%%zmm0\n\t")
TPUT512(tput_rolq, "vprolq $7,%%zmm17,")
LAT512(lat_ternlog, "vpternlogq $0x96,%%zmm16,%%zmm17,%%zmm0\n\t")
TPUT512(tput_ternlog, "vpternlogq $0x96,%%zmm16,%%zmm17,")

LAT512(lat_mullw, "vpmullw %%zmm16,%%zmm0,%%zmm0\n\t")
TPUT512(tput_mullw, "vpmullw %%zmm16,%%zmm17,")
LAT512(lat_mulhw, "vpmulhw %%zmm16,%%zmm0,%%zmm0\n\t")
TPUT512(tput_mulhw, "vpmulhw %%zmm16,%%zmm17,")

LAT512(lat_subd, "vpsubd %%zmm16,%%zmm0,%%zmm0\n\t")
TPUT512(tput_subd, "vpsubd %%zmm16,%%zmm17,")
LAT512(lat_packssdw, "vpackssdw %%zmm16,%%zmm0,%%zmm0\n\t")
TPUT512(tput_packssdw, "vpackssdw %%zmm16,%%zmm17,")

/* ---- mixed-port probes: does IFMA share a port with shifts? ---- */
#define MIXT(NAME, A, B)                                                 \
    static uint64_t NAME(uint64_t it)                                    \
    {                                                                    \
        uint64_t c0 = 0, c1 = 0;                                         \
        __asm__ volatile(SET512 TS0 "1:\n\t" R4(                         \
                             MIX16(A, B)) "decq %[it]\n\tjnz 1b\n\t" TS1 \
                         : [c0] "=&r"(c0), [c1] "=&r"(c1), [it] "+r"(it) \
                         :                                               \
                         : CL512);                                       \
        return c1 - c0;                                                  \
    }
/* IFMA + immediate-shift: both p0 -> expect ~1/cyc total (they contend).
 */
MIXT(mix_madd_sll, "vpmadd52luq %%zmm16,%%zmm17,", "vpsllq $1,%%zmm17,")
/* IFMA + shuffle (p5): distinct ports -> expect ~2/cyc total. */
MIXT(mix_madd_perm, "vpmadd52luq %%zmm16,%%zmm17,",
     "vpermq $0x1b,%%zmm17,")
/* IFMA + add (p0/p5): expect ~2/cyc (add steals p5, IFMA keeps p0). */
MIXT(mix_madd_add, "vpmadd52luq %%zmm16,%%zmm17,",
     "vpaddq %%zmm16,%%zmm17,")
/* IFMA-lo + IFMA-hi (same op pair as mulsh52): both p0 -> ~1/cyc. */
MIXT(mix_lo_hi, "vpmadd52luq %%zmm16,%%zmm17,",
     "vpmadd52huq %%zmm16,%%zmm17,")

/* ---- scalar x64 baseline (for the integer-exp comparison) ---- */
static uint64_t lat_add64(uint64_t it) /* ANCHOR = 1 core cycle (scalar) */
{
    uint64_t c0 = 0, c1 = 0;
    __asm__ volatile(
        "movq $0x9E3779B97F4A7C15,%%rbx\n\t" TS0
        "movq $1,%%rax\n\t"
        "1:\n\t" R32("addq %%rbx,%%rax\n\t") "decq %[it]\n\tjnz 1b\n\t" TS1
        : [c0] "=&r"(c0), [c1] "=&r"(c1), [it] "+r"(it)
        :
        : "rax", "rbx", "rdx", "rcx", "cc", "memory");
    return c1 - c0;
}
static uint64_t lat_mulx(uint64_t it) /* dependent chain on the HIGH 64b */
{
    uint64_t c0 = 0, c1 = 0;
    __asm__ volatile(
        "movq $0x9E3779B97F4A7C15,%%rsi\n\t" TS0
        "movq $0xD1B54A32D192ED03,%%rdx\n\t"
        "1:\n\t" R32(
            "mulx %%rsi,%%r8,%%rdx\n\t") "decq %[it]\n\tjnz 1b\n\t" TS1
        : [c0] "=&r"(c0), [c1] "=&r"(c1), [it] "+r"(it)
        :
        : "rax", "rdx", "rcx", "rsi", "r8", "cc", "memory");
    return c1 - c0;
}
static uint64_t tput_mulx(uint64_t it) /* 4 independent mulx per group */
{
    uint64_t c0 = 0, c1 = 0;
    __asm__ volatile(
        "movq $0x9E3779B97F4A7C15,%%rsi\n\t" TS0
        "movq $0xD1B54A32D192ED03,%%rdx\n\t"
        "1:\n\t" R8("mulx %%rsi,%%r8,%%r9\n\t mulx %%rsi,%%r10,%%r11\n\t "
                    "mulx %%rsi,%%r12,%%r13\n\t mulx "
                    "%%rsi,%%r14,%%r15\n\t") "decq %[it]\n\tjnz 1b\n\t" TS1
        : [c0] "=&r"(c0), [c1] "=&r"(c1), [it] "+r"(it)
        :
        : "rax", "rdx", "rcx", "rsi", "r8", "r9", "r10", "r11", "r12",
          "r13", "r14", "r15", "cc", "memory");
    return c1 - c0;
}

/* ---- driver ---- */
typedef uint64_t (*fn_t)(uint64_t);
typedef struct {
    const char *name;
    fn_t f;
    uint64_t iters;
    int ops_per_iter;
    int scalar; /* 1 -> normalise by scalar anchor, else 512 anchor */
    int lat;    /* 1 -> latency test (else throughput) */
} test_t;

static double run_min(fn_t f, uint64_t iters, int ops)
{
    f(iters); /* warm-up: AVX-512 power-up + frequency transition */
    uint64_t best = ~0ULL;
    for (int r = 0; r < 11; r++) {
        uint64_t c = f(iters);
        if (c < best)
            best = c;
    }
    return (double)best / (double)(iters * (uint64_t)ops);
}

int main(void)
{
    const uint64_t L = 400000, T = 150000;
    test_t tests[] = {
        {"vpaddq", lat_vpaddq, L, 32, 0, 1},
        {"vpaddq", tput_vpaddq, T, 64, 0, 0},
        {"vpmadd52luq", lat_madd52lo, L, 32, 0, 1},
        {"vpmadd52luq", tput_madd52lo, T, 64, 0, 0},
        {"vpmadd52huq", lat_madd52hi, L, 32, 0, 1},
        {"vpmadd52huq", tput_madd52hi, T, 64, 0, 0},
        {"vpermt2q", lat_permt2q, L, 32, 0, 1},
        {"vpermt2q", tput_permt2q, T, 64, 0, 0},
        {"vpermq", tput_permq, T, 64, 0, 0},
        {"vpsllvq", lat_sllvq, L, 32, 0, 1},
        {"vpsllvq", tput_sllvq, T, 64, 0, 0},
        {"vpsllq(imm)", lat_sllq, L, 32, 0, 1},
        {"vpsllq(imm)", tput_sllq, T, 64, 0, 0},
        {"vprolq", lat_rolq, L, 32, 0, 1},
        {"vprolq", tput_rolq, T, 64, 0, 0},
        {"vpternlogq", lat_ternlog, L, 32, 0, 1},
        {"vpternlogq", tput_ternlog, T, 64, 0, 0},
        {"vpmullw", lat_mullw, L, 32, 0, 1},
        {"vpmullw", tput_mullw, T, 64, 0, 0},
        {"vpmulhw", lat_mulhw, L, 32, 0, 1},
        {"vpmulhw", tput_mulhw, T, 64, 0, 0},
        {"vpsubd", lat_subd, L, 32, 0, 1},
        {"vpsubd", tput_subd, T, 64, 0, 0},
        {"vpackssdw", lat_packssdw, L, 32, 0, 1},
        {"vpackssdw", tput_packssdw, T, 64, 0, 0},
        {"mix lo+sllq", mix_madd_sll, T, 64, 0, 0},
        {"mix lo+permq", mix_madd_perm, T, 64, 0, 0},
        {"mix lo+addq", mix_madd_add, T, 64, 0, 0},
        {"mix lo+hi", mix_lo_hi, T, 64, 0, 0},
        {"add r64", lat_add64, L, 32, 1, 1},
        {"mul r64->hi", lat_mulx, L, 32, 1, 1},
        {"mulx r64", tput_mulx, T, 32, 1, 0},
    };
    int N = (int)(sizeof(tests) / sizeof(tests[0]));

    /* anchors first (1 core cycle each, in their freq domain). */
    double anc512 = run_min(lat_vpaddq, L, 32);
    double anc64 = run_min(lat_add64, L, 32);
    printf("# anchors (rdtsc ref-cyc/op): vpaddq-lat=%.4f  add-lat=%.4f\n",
           anc512, anc64);
    printf(
        "# implied core freq: 512-domain=%.2f GHz  scalar=%.2f GHz "
        "(ref TSC 3.6 GHz)\n",
        3.6 / anc512, 3.6 / anc64);
    printf("# %-13s %-6s  %10s  %12s  %10s\n", "insn", "kind",
           "ref-cyc/op", "core-cyc/op", "ops/cyc");
    for (int i = 0; i < N; i++) {
        double raw =
            run_min(tests[i].f, tests[i].iters, tests[i].ops_per_iter);
        double anc = tests[i].scalar ? anc64 : anc512;
        double core = raw / anc;
        printf("  %-13s %-6s  %10.4f  %12.3f  %10.3f\n", tests[i].name,
               tests[i].lat ? "lat" : "tput", raw, core,
               tests[i].lat ? 0.0 : 1.0 / core);
    }
    return 0;
}
