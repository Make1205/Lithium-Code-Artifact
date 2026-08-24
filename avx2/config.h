#ifndef LITHIUM_CONFIG_H
#define LITHIUM_CONFIG_H

/*
 * Compile-time configuration for Lithium (AVX2 build).
 *
 * Identical parameter selection to ref/config.h, but the symbol namespace
 * carries an `_avx2` suffix so the AVX2 library can coexist with the
 * reference library without symbol clashes.  The KAT (pk/sk/sig bytes) is
 * independent of symbol names, so the AVX2 default build (integer exp)
 * stays byte-identical to the reference KAT.
 *
 * Override at build time with: -DLITHIUM_MODE=120|180|260
 */

#ifndef LITHIUM_MODE
#    define LITHIUM_MODE 120
#endif

#if LITHIUM_MODE == 120
#    define CRYPTO_ALGNAME "Lithium-120-AVX2"
#    define LITHIUM_NAMESPACETOP lithium_nist120_avx2
#    define LITHIUM_NAMESPACE(s) lithium_nist120_avx2_##s
#elif LITHIUM_MODE == 180
#    define CRYPTO_ALGNAME "Lithium-180-AVX2"
#    define LITHIUM_NAMESPACETOP lithium_nist180_avx2
#    define LITHIUM_NAMESPACE(s) lithium_nist180_avx2_##s
#elif LITHIUM_MODE == 260
#    define CRYPTO_ALGNAME "Lithium-260-AVX2"
#    define LITHIUM_NAMESPACETOP lithium_nist260_avx2
#    define LITHIUM_NAMESPACE(s) lithium_nist260_avx2_##s
#else
#    error "Unsupported LITHIUM_MODE (expected 120, 180, or 260)"
#endif

/* Symbol renames (lithium_nist<mode>_avx2_*) all live in one place; see
 * namespace.h (symlinked to ../ref/namespace.h -- it uses
 * LITHIUM_NAMESPACE, which carries the _avx2 suffix here).  Included after
 * LITHIUM_NAMESPACE is defined.  Disable with -DDISABLE_NAMESPACE=1. */
#include "namespace.h"

#endif /* LITHIUM_CONFIG_H */
