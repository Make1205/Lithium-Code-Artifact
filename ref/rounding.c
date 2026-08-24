/*
 * rounding.c - IRS mod-2q HighBits / LSB / hint helpers.
 *
 * tau is a power of two (512/512/1024) and divides 2(q-1); the bucket
 * count W1_MAX = 2(q-1)/tau (225/121/126) is NOT a power of two, so the
 * bucket index uses %.  "0 and 2(q-1) are the same bucket" (Lithium.md
 * A.1) falls out of the (mod W1_MAX) reduction.
 */
/* All reductions here are division-free / constant-time (no `% const`):
 * the mod-2q goes through reduce_mod_2q (Barrett, reduce.c); the `%
 * W1_MAX` operands are already nearly reduced (bounded ranges proven
 * below), so a single masked add/subtract suffices.  This runs on the
 * SECRET commitment w during signing. */
#include "rounding.h"

#include <stdint.h>

#include "params.h"
#include "reduce.h"

#if TAU == 1024
#    define TAU_LOG2 10
#elif TAU == 512
#    define TAU_LOG2 9
#else
#    error "TAU must be a power of two (512/1024)"
#endif

int32_t reduce2q(int32_t a)
{
    return reduce_mod_2q(a); /* Barrett mod-2q in [0, 2q) */
}

/* HighBits when the input is ALREADY in [0, 2q): skips the redundant
 * Barrett mod-2q that highbits() applies up front.  The avx2/avx512
 * highbits8/16 already omit that reduce on pre-bounded inputs (byte-exact
 * KAT), which proves it is a no-op here -- so callers feeding the matmul
 * mod-2q lift output / usehint's wtilde can use this directly. */
int32_t highbits_reduced(int32_t x)
{
    int32_t b = (x + (TAU >> 1)) >>
                TAU_LOG2; /* explicit shift; never idiv under -Os */
    /* b %= W1_MAX, b in [0, W1_MAX]: only b==W1_MAX needs the single
     * subtract. */
    b -= (-(int32_t)(b >= W1_MAX)) & W1_MAX;
    return b; /* [0, W1_MAX) */
}

int32_t highbits(int32_t w)
{
    return highbits_reduced(
        reduce2q(w)); /* reduce arbitrary w to [0,2q) */
}

int32_t hbvalue(int32_t bucket)
{
    return bucket * TAU;
}

int32_t lsb(int32_t a)
{
    return reduce2q(a) & 1;
}

/* Centered representative in (-m/2, m/2].  Precondition (holds for the
 * only caller, makehint with m = W1_MAX): |x| < m, so x % m is the
 * identity and we just normalise branchlessly. */
static int32_t centermod(int32_t x, int32_t m)
{
    x += (x >> 31) & m;               /* x<0 -> +m: now [0, m) */
    x -= (-(int32_t)(x > m / 2)) & m; /* x>m/2 -> -m  (m/2 = m>>1) */
    return x;
}

int32_t makehint(int32_t w, int32_t wtilde)
{
    return centermod(highbits(w) - highbits(wtilde), W1_MAX);
}

/* makehint when HighBits(w) is already in hand (the signing buckets[][]):
 * avoids recomputing highbits(w).  wtilde is the per-coeff reduce2q
 * output, always in [0,2q), so highbits_reduced applies. */
int32_t makehint_hi(int32_t hi_w, int32_t wtilde)
{
    return centermod(hi_w - highbits_reduced(wtilde), W1_MAX);
}

int32_t usehint(int32_t wtilde, int32_t hint)
{
    /* b0 = highbits(wtilde) + hint in (-W1_MAX, 2*W1_MAX); reduce to [0,
     * W1_MAX) with one masked add (b0<0) + one masked subtract
     * (b0>=W1_MAX).  wtilde is in [0,2q) at every call site (verifier wt,
     * and the test/dudect callers all pass a reduce2q output), so use the
     * reduced HighBits. */
    int32_t b = highbits_reduced(wtilde) + hint;
    b += (b >> 31) & W1_MAX; /* -> [0, 2*W1_MAX) */
    b -= (-(int32_t)(b >= W1_MAX)) & W1_MAX;
    return b; /* [0, W1_MAX) */
}
