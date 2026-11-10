#ifndef COPA_IRS_LIMITS_H
#define COPA_IRS_LIMITS_H
/* Generated v3.1 limits. Public fixed depths; no per-key shortening.
 * Exact square-convolution certificate: test/copa/v31/irs_lower_norm.json.
 * Keep legacy FULL/TIGHT research controls and the original 5a depth. */
#if !defined(COPA_FULL_IRS) && !defined(COPA_TIGHT_IRS)
#if COPA_MODE==2
#define COPA_PARTIAL_MIN_NORM 447
#define COPA_MIN_NORM 447
#define COPA64_IRS_TERMS 48
#elif COPA_MODE==3
#define COPA_PARTIAL_MIN_NORM 253
#define COPA_MIN_NORM 253
#define COPA64_IRS_TERMS 24
#elif COPA_MODE==5
#define COPA_PARTIAL_MIN_NORM 0
#define COPA_MIN_NORM 158
#define COPA64_IRS_TERMS 40
#elif COPA_MODE==6
#define COPA_PARTIAL_MIN_NORM 344
#define COPA_MIN_NORM 344
#define COPA64_IRS_TERMS 24
#endif
#else
#define COPA_PARTIAL_MIN_NORM 0
#if defined(COPA_TIGHT_IRS) || (!defined(COPA_FULL_IRS) && (COPA_MODE==2 || COPA_MODE==5))
#if COPA_MODE==2
#define COPA_MIN_NORM 28
#define COPA64_IRS_TERMS 168
#elif COPA_MODE==3
#define COPA_MIN_NORM 98
#define COPA64_IRS_TERMS 40
#elif COPA_MODE==5
#define COPA_MIN_NORM 158
#define COPA64_IRS_TERMS 40
#elif COPA_MODE==6
#define COPA_MIN_NORM 160
#define COPA64_IRS_TERMS 40
#endif
#else
#define COPA_MIN_NORM 1
#if COPA_MODE==2
#define COPA64_IRS_TERMS 800
#elif COPA_MODE==3
#define COPA64_IRS_TERMS 312
#elif COPA_MODE==5
#define COPA64_IRS_TERMS 392
#elif COPA_MODE==6
#define COPA64_IRS_TERMS 384
#endif
#endif
#endif
#endif
