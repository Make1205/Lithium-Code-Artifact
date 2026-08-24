# RCDT Table Notes

`tools/gen_rcdt.py` generates the reverse-CDT tables embedded in the `@@AUTOGEN:rcdt@@` regions of `ref/sampler.c` and `avx2/sampler.c`. The generated tables cover the BLISS-style base sigma2 distribution and the sigma=1 key-noise distribution used by the Lithium sampler.

The generator verifies that every table is monotone in the expected direction and emits the canonical limb representation consumed by the scalar and SIMD borrow-chain comparisons. Regenerate through `make tables` from `ref/`, or through `make check-consts` in optimized builds when validating that committed tables are unchanged.
