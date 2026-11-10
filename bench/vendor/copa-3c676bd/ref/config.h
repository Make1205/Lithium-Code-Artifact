#ifndef LITHIUM_CONFIG_H
#define LITHIUM_CONFIG_H

/*
 * Compile-time configuration for IRS.
 *
 * LITHIUM_MODE selects the parameter set:
 *   120 -> Lithium-120  (n=256, q=57601, incomplete NTT base_deg=2)
 *   180 -> Lithium-180  (n=256, q=30977, incomplete NTT base_deg=2)
 *   260 -> Lithium-260  (n=512, q=64513, full NTT base_deg=1)
 *
 * Override at build time with: -DLITHIUM_MODE=120|180|260
 */

#ifndef LITHIUM_MODE
#    define LITHIUM_MODE 120
#endif

#if LITHIUM_MODE == 120
#    define CRYPTO_ALGNAME "Lithium-120"
#    define LITHIUM_NAMESPACETOP lithium_nist120_ref
#    define LITHIUM_NAMESPACE(s) lithium_nist120_ref_##s
#elif LITHIUM_MODE == 180
#    define CRYPTO_ALGNAME "Lithium-180"
#    define LITHIUM_NAMESPACETOP lithium_nist180_ref
#    define LITHIUM_NAMESPACE(s) lithium_nist180_ref_##s
#elif LITHIUM_MODE == 260
#    define CRYPTO_ALGNAME "Lithium-260"
#    define LITHIUM_NAMESPACETOP lithium_nist260_ref
#    define LITHIUM_NAMESPACE(s) lithium_nist260_ref_##s
#else
#    error "Unsupported LITHIUM_MODE (expected 120, 180, or 260)"
#endif

/* Symbol renames (lithium_nist<mode>_ref_*) all live in one place; see
 * namespace.h. Included here, after LITHIUM_NAMESPACE is defined, so every
 * translation unit that reaches config.h (the whole tree, via params.h)
 * gets the renames.  Disable with -DDISABLE_NAMESPACE=1. */
#include "namespace.h"

#endif /* LITHIUM_CONFIG_H */
