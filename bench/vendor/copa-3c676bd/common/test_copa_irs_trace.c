/* Compare deterministic fixed-step event counts independently of wall
 * time. */
#include <stdint.h>
#include <stdlib.h>
#include <valgrind/callgrind.h>

#include "copa_irs64.h"
static volatile uint64_t sink;
int main(int argc, char **argv)
{
    unsigned cls = argc > 1 ? (unsigned)strtoul(argv[1], 0, 10) : 0;
    copa_irs64_context ctx;
    copa_irs64_prepare(&ctx, cls ? 64 : 1);
    uint64_t random[6] = {1, 2, 3, 4, 5, 6};
    unsigned minus = 0;
    CALLGRIND_START_INSTRUMENTATION;
    CALLGRIND_TOGGLE_COLLECT;
    for (unsigned i = 0; i < 100; ++i) {
        unsigned ok = copa_irs64_step(&minus, cls ? -64 : 0, &ctx, random);
        sink += (uint64_t)ok + minus;
    }
    CALLGRIND_TOGGLE_COLLECT;
    CALLGRIND_STOP_INSTRUMENTATION;
    return 0;
}
