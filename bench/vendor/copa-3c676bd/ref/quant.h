#ifndef COPA_QUANT_H
#define COPA_QUANT_H
#include <stdint.h>
#include <stddef.h>
/* Compile-time mode aliases will use this shared byte-level implementation.
 * These functions accept public mode 2,3,5,6 and row-major ring blocks.
 * Coordinate inputs satisfy |x_i| <= 2^20.
 * Mode 6's ring block has 12 coordinates: eight for 2E plus four 3Z.
 */
unsigned copa_quant_width(unsigned mode);
unsigned copa_modulus(unsigned mode);
unsigned copa_rows(unsigned mode);
void copa_quant_block(int32_t *out, const int32_t *in, unsigned mode);
/* Lab/Unlab operate on the E/BW block only (first 8 coordinates for mode 6).
 * Return 0 for a nonmember / noncanonical label; otherwise return 1.
 */
unsigned copa_label_width(unsigned mode);
unsigned copa_label_radix(unsigned mode, unsigned index);
int copa_lab(uint32_t *label, const int32_t *x, unsigned mode);
int copa_unlab(int32_t *x, const uint32_t *label, unsigned mode);
/* Dithering: SHAKE256(rho || 0x14 || LE16(lane)), eight contiguous chunks
 * of equal length in the canonical block/RM stream.  Rejection sampling
 * draws LE16 values with threshold floor(65536/range)*range.
 * Output: row-major coefficient array of length rows(mode)*256.
 */
void copa_expand_d(int32_t *out, const uint8_t rho[32], unsigned mode);
#endif
