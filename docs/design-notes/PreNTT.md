# Pre-NTT Matrix Expansion Notes

`A0` is a public uniform matrix used in signing and verification products. The implementation expands it directly into each backend's NTT-domain layout so repeated products avoid a separate forward transform of `A0`.

Uniformity is preserved because each backend NTT is a bijection on the corresponding product ring. The AVX layouts are checked by `tools/gen_ntt_perm.py`, which proves that the AVX NTT output order is transpose-equivalent to the scalar `ntt16` order and therefore can be reached by the emitted `nttunpack` shuffle ladder.
