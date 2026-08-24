/*
 * prof.h - macro-gated profiling probes for the IRS sign path.
 *
 * Two independent, opt-in instrumentation systems (both no-ops unless
 * their macro is defined, so normal/KAT/CT builds are byte-identical and
 * zero-cost):
 *
 *   -DPROF_TIME : per-component cycle accounting.  Disjoint
 * regions are timed with PROF_START/PROF_STOP; top-level buckets sum to
 * the sign total, sub-buckets (gauss internals, NTT stages, sign_irs exp)
 * are reported as a breakdown of their parent.  (Probes use rdtscp; their
 * ~50-cycle overhead is reported and matters only for the very hot inner
 * buckets.)
 *
 *   -DPROF_RAND : randomness accounting.  PROF_SQ in the SHAKE
 * squeeze functions tallies bytes squeezed against the current PROF_CTX
 * context (the "who is pulling randomness" label set around each sign
 * region); PROF_USE tallies the fine per-field consumption inside the
 * Gaussian sampler (signs / sigma2 / y / rejection-tail).
 *
 * `make profile` (avx2/ or ref/) defines both and runs
 * test/speed_profile.c.
 */
#ifndef LITHIUM_PROF_H
#define LITHIUM_PROF_H

#include <stddef.h>
#include <stdint.h>

/* ---- timing buckets (order = report order; parents precede their
 * children) ---- */
enum {
    /* tr/mu/rhoprime derivation (SHAKE). */
    PT_SETUP = 0,
    /* expand_A0 matrix expansion. */
    PT_EXPAND_A0,
    /* Parent bucket for Gaussian y-vector sampling. */
    PT_SAMPLE_Y,
    /* Gaussian sampler: Keccak absorb and squeeze. */
    PT_G_SHAKE,
    /* Gaussian sampler: sigma2 via RCDT/BaseSampler. */
    PT_G_BASESAMP,
    /* Gaussian sampler: ApproxExp batch. */
    PT_G_APPROXEXP,
    /* Gaussian sampler: a_q60, finalize, and y-unpack. */
    PT_G_FINAL,
    /* Parent bucket for commitment w = A y. */
    PT_COMMIT,
    /* Commitment NTT forward pass. */
    PT_NTT_FWD,
    /* Commitment pointwise Montgomery products. */
    PT_NTT_PW,
    /* Commitment inverse NTT pass. */
    PT_NTT_INV,
    /* HighBits/LSB bucketing. */
    PT_HIGHBITS,
    /* Challenge hashing and Fisher-Yates sampling. */
    PT_CHALLENGE,
    /* Parent bucket for iterative rejection sampling. */
    PT_SIGN,
    /* IRS ApproxExp work. */
    PT_EXP,
    /* ||z|| <= B_s check. */
    PT_NORMCHECK,
    /* Hint generation. */
    PT_MAKEHINT,
    /* Signature packing, including rANS. */
    PT_PACK,
    /* ---- keygen buckets (crypto_sign_keypair) ---- */
    PT_KG_EXPAND_A0,
    /* Keygen noise sampling: f0, s0, e. */
    PT_KG_NOISE,
    /* Keygen norm window. */
    PT_KG_NORM,
    /* compute_b: f^-1 + (A0 s0+e) f^-1. */
    PT_KG_COMPUTEB,
    /* pack_pk, tr, and pack_sk. */
    PT_KG_PACK,
    /* ---- verify buckets (crypto_sign_verify) ---- */
    PT_VF_UNPACK,
    /* Verify tr and mu derivation. */
    PT_VF_SETUP,
    /* Verify expand_A0_ntt and b-hat NTT. */
    PT_VF_A0,
    /* Verify commitment reconstruction. */
    PT_VF_MATMUL,
    /* Verify usehint / lsb bucketing. */
    PT_VF_HINT,
    /* Verify challenge hashing and sampling. */
    PT_VF_CHALLENGE,
    /* Verify ||z|| reconstruction and check. */
    PT_VF_NORM,
    PT_NBUCKETS
};

/* ---- randomness contexts (SHAKE squeeze attribution) ---- */
enum {
    /* tr/mu/rhoprime. */
    PC_SETUP = 0,
    /* Matrix expansion via SHAKE128. */
    PC_A0,
    /* Gaussian y sampling via SHAKE256. */
    PC_GAUSS,
    /* Challenge hash and Fisher-Yates sampling. */
    PC_CHALLENGE,
    /* sign_irs randomness. */
    PC_IRS,
    /* Keygen noise via SHAKE256. */
    PC_NOISE,
    /* Anything else. */
    PC_OTHER,
    PC_NCTX
};

/* ---- fine randomness consumption (within the samplers) ---- */
enum {
    /* Gaussian sign bits. */
    PU_SIGNS = 0,
    /* Gaussian RCDT draws. */
    PU_SIGMA2,
    /* Gaussian y bits. */
    PU_Y,
    /* Gaussian rejection-tail bytes. */
    PU_GTAIL,
    /* Keygen noise bytes. */
    PU_NOISE,
    PU_NCONS
};

#if defined(PROF_TIME) || defined(PROF_RAND)
extern uint64_t prof_cyc[PT_NBUCKETS];
extern uint64_t prof_cnt[PT_NBUCKETS];
extern uint64_t prof_sq[PC_NCTX];
extern uint64_t prof_use[PU_NCONS];
extern int prof_ctx; /* current squeeze attribution context */
void prof_reset(void);
void prof_report(const char *title, uint64_t nsig);
static inline uint64_t prof_rdtsc(void)
{
    unsigned lo, hi, aux;
    __asm__ __volatile__("rdtscp" : "=a"(lo), "=d"(hi), "=c"(aux));
    return ((uint64_t)hi << 32) | lo;
}
#endif

#ifdef PROF_TIME
#    define PROF_START(v) uint64_t v = prof_rdtsc()
#    define PROF_STOP(b, v)                      \
        do {                                     \
            prof_cyc[(b)] += prof_rdtsc() - (v); \
            prof_cnt[(b)]++;                     \
        } while (0)
#else
#    define PROF_START(v)
#    define PROF_STOP(b, v)
#endif

#ifdef PROF_RAND
#    define PROF_CTX(c) (prof_ctx = (c))
#    define PROF_CTX_GET() (prof_ctx)
#    define PROF_CTX_SET(c) (prof_ctx = (c))
#    define PROF_SQ(nbytes) (prof_sq[prof_ctx] += (uint64_t)(nbytes))
#    define PROF_USE(t, n) (prof_use[(t)] += (uint64_t)(n))
#else
#    define PROF_CTX(c) ((void)0)
#    define PROF_CTX_GET() (0)
#    define PROF_CTX_SET(c) ((void)0)
#    define PROF_SQ(nbytes) ((void)0)
#    define PROF_USE(t, n) ((void)0)
#endif

#endif /* LITHIUM_PROF_H */
