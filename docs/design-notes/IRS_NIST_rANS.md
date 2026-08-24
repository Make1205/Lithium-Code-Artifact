# IRS rANS Notes

The compact signature uses one merged two-state rANS stream for rounded z-head and MakeHint symbols, plus a packed z-low bitstream. The reserve `LITHIUM_RANS_RESERVED_BYTES` is sized by `tools/SigSize.py` for a strict `2^-100` overflow target under the true emitted PMFs.

`tools/gen_rans_tables.py` emits the canonical encode/decode tables into `ref/irs_rans.h`. `tools/check_rans.py` and the C tests validate byte-exact agreement between the Python model and the C implementation. The verifier authenticates zero padding in the reserved rANS tail, preventing length-padding malleability.
