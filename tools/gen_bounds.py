#!/usr/bin/env python3
"""Generate EXACT integer squared norm thresholds for IRS and patch ref/params.h.

The norm checks compare an EXACT integer squared-norm `norm_sq` against a real
threshold B^2. We must NOT hard-truncate the full-precision bounds from
param.md (e.g. B_k = 46.00893646...). Instead we emit the unique integer
threshold that makes the integer test bit-for-bit equivalent to the real test:

  KeyGen : accept iff ||s||  <  B_k   (strict)   -> norm_sq <= BK_SQ,  BK_SQ = ceil(B_k^2) - 1
  Sign   : accept iff ||z||  <= B_s   (B_s int)  -> norm_sq <= BS_SQ,  BS_SQ = floor(B_s^2)
  Verify : accept iff ||z'|| <= B_v              -> norm_sq <= BV_SQ,  BV_SQ = floor(B_v^2)

For a non-integer B, `norm_sq < B^2  <=>  norm_sq <= floor(B^2)` (norm_sq integer);
ceil(B^2)-1 == floor(B^2) when B^2 is not an integer and == B^2-1 when it is,
giving the correct strict-`<` threshold in both cases.

Full-precision bounds are taken verbatim from docs/design-notes/param.md.
Worst-case norm_sq ~ B_v^2 ~ 4.3e8 fits int64; no scaling needed.

Run:  python3 gen_bounds.py
"""

import os
import autogen
from decimal import Decimal, getcontext

getcontext().prec = 80

# Full-precision (verbatim from param.md): (B_k, B_k', B_s, B_v).
# B_k' = bk_prime is the keygen LOWER gate: accepted keys satisfy
#   B_k' <= ||s|| < B_k.  The lower gate bounds gamma = ||s||^2/(2r^2) from
# below (gamma_min at ||s||=B_k'), which is what fixes the IRS truncation depth
# T (see tools/irs_ref.py).  Both gates -> ~50% keygen acceptance.
BOUNDS = {
    "120": (Decimal("46.00893646172169"),  Decimal("44.5007315301564"),
            Decimal("8090"),  Decimal("10106.329462945841")),
    "180": (Decimal("51.3038533806998"),   Decimal("49.889031744688346"),
            Decimal("10396"), Decimal("12793.701585698811")),
    "260": (Decimal("64.75410246578265"),  Decimal("63.245897534217356"),
            Decimal("15210"), Decimal("20646.021465119095")),
}


def floor_dec(x):
    return int(x.to_integral_value(rounding="ROUND_FLOOR"))


def ceil_dec(x):
    return int(x.to_integral_value(rounding="ROUND_CEILING"))


BLOCK = '''\
#if LITHIUM_MODE == {sec}
/* B_k = {bk}  B_k' = {bkp}  B_s = {bs}  B_v = {bv} */
#define BK_SQ      INT64_C({bk_sq})
#define BK_LOW_SQ  INT64_C({bkl_sq})
#define BS_SQ      INT64_C({bs_sq})
#define BV_SQ      INT64_C({bv_sq})
#endif
'''


def main():
    blocks, summary = [], []
    print("IRS gen_bounds.py -- exact integer squared thresholds\n")
    for sec, (bk, bkp, bs, bv) in BOUNDS.items():
        bk2, bkp2, bs2, bv2 = bk * bk, bkp * bkp, bs * bs, bv * bv
        bk_sq = ceil_dec(bk2) - 1     # ||s|| < B_k     -> norm_sq <= bk_sq
        bkl_sq = ceil_dec(bkp2)       # ||s|| >= B_k'   -> norm_sq >= bkl_sq
        bs_sq = floor_dec(bs2)        # ||z|| <= B_s
        bv_sq = floor_dec(bv2)        # ||z'|| <= B_v
        print(f"--- IRS-{sec} ---")
        print(f"  B_k ={bk}  B_k^2 ={bk2}  -> BK_SQ=ceil(B_k^2)-1={bk_sq}  (||s|| < B_k)")
        print(f"  B_k'={bkp}  B_k'^2={bkp2}  -> BK_LOW_SQ=ceil(B_k'^2)={bkl_sq}  (||s|| >= B_k')")
        print(f"  B_s ={bs}  B_s^2 ={bs2}  -> BS_SQ=floor(B_s^2)={bs_sq}")
        print(f"  B_v ={bv}  B_v^2 ={bv2}  -> BV_SQ=floor(B_v^2)={bv_sq}")
        summary.append(f" *   IRS-{sec}: BK_LOW_SQ={bkl_sq} <= ||s||^2 <= BK_SQ={bk_sq}, "
                       f"BS_SQ={bs_sq}, BV_SQ={bv_sq}")
        blocks.append(BLOCK.format(sec=sec, bk=bk, bkp=bkp, bs=bs, bv=bv,
                                   bk_sq=bk_sq, bkl_sq=bkl_sq, bs_sq=bs_sq, bv_sq=bv_sq))
    body = ("/* EXACT integer squared norm thresholds (see tools/gen_bounds.py).\n"
            " *   KeyGen accept iff BK_LOW_SQ <= norm_sq <= BK_SQ (B_k' <= ||s|| < B_k)\n"
            " *   Sign   accept iff norm_sq <= BS_SQ;  Verify accept iff norm_sq <= BV_SQ\n"
            + "\n".join(summary) + "\n */\n"
            + "\n".join(blocks))
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.normpath(os.path.join(here, "..", "ref", "params.h"))
    autogen.patch_region(path, "bounds", body)
    print(f"\nPatched @@AUTOGEN:bounds@@ region in {path}")


if __name__ == "__main__":
    main()
