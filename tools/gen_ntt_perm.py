#!/usr/bin/env python3
"""gen_ntt_perm.py - transpose-equivalence GUARD for the avx NTT vs the scalar
ntt16 and verifies that no separate scramble table is needed.

Background (Pre-NTT / seed-to-NTT, docs/design-notes/PreNTT.md): A0 is generated
DIRECTLY in the NTT domain (expand_A0_ntt) -- the uniform [0,q) samples ARE the
forward-NTT coefficients in the *scalar* ntt16 (standard) order.  The avx
expand_A0_ntt lays them into the avx ntt.S slot order with `nttunpack_avx`
(Kyber gen_matrix / Dilithium ExpandA style: the forward shuffle network with no
butterflies).  That shuffle ladder is only correct if the avx NTT and the scalar
ntt16 are TRANSPOSE-EQUIVALENT, i.e. their outputs differ by a pure
bit-permutation of the coefficient index:

    avx_ntt(p)[s] == ref_ntt(p)[P[s]]   with P a bit-permutation (P(a^b)=P(a)^P(b)).

This holds because the avx generators (tools/ntt/<qset>/gen_ntt.py) now use the
SAME primitive root as ntt16 (gen_zetas.py's smallest primitive root).  A
multiplicative-scramble root (the old PSI=rho^k choice) would make P NON-affine,
which no shuffle ladder can realize -- so this guard is what pins the alignment.

This script DERIVES P purely in Python from the two reference NTT models
(tools/ntt_ref.py mirrors ntt16; tools/ntt/<qset>/gen_ntt.py is the
instruction-accurate avx model) by value-matching forward outputs, VALIDATES it
on fresh inputs, and ASSERTS it is a bit-permutation.  It emits no files; it is a
pass/fail gate (wired into `make check-consts`).  The runtime correctness of the
emitted nttunpack_avx (that it realizes exactly this permutation) is pinned
separately by tools/ntt/<qset>/test_ntt.c ("nttunpack(ntt_ref(p)) == ntt(p)") and
end-to-end by the avx2-vs-ref KAT (`make check-kat`).

Run:  python3 gen_ntt_perm.py
"""
import importlib.util
import os
import random
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ntt_ref as R  # forward ntt16 model + montgomery reduce

# mode -> (in-tree NTT generator qset, q, n, base_deg)
MODES = {
    "120": ("q57601n256", 57601, 256, 2),
    "180": ("q30977n256", 30977, 256, 2),
    "260": ("q64513n512", 64513, 512, 1),
}


def parse_ntt_h(sec):
    """Pull QINV + the zetas[] table for one mode out of ref/ntt.h (the same
    generated source of truth ntt16 uses).  Handles IndentPPDirectives
    (`#    define`)."""
    txt = open(os.path.join(HERE, "..", "ref", "ntt.h")).read()
    blk = re.search(
        r"#if LITHIUM_MODE == %s\b(.*?)#endif /\* LITHIUM_MODE == %s"
        % (sec, sec), txt, re.S).group(1)
    qinv = int(re.search(r"#\s*define QINV\s+(\d+)", blk).group(1))
    zm = re.search(r"zetas\[NTT_ZETAS_SIZE\] = \{(.*?)\};", blk, re.S)
    zetas = [int(x) for x in re.findall(r"-?\d+", zm.group(1))]
    return qinv, zetas


def load_avx_ntt(qset, sec):
    path = os.path.join(HERE, "ntt", qset, "gen_ntt.py")
    spec = importlib.util.spec_from_file_location(
        "irs_ntt_" + sec, os.path.abspath(path))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)  # import-only; main() is __main__-guarded
    return m


def derive_perm(sec, qset, q, n, bd):
    qinv, zetas = parse_ntt_h(sec)
    avx = load_avx_ntt(qset, sec)
    assert getattr(avx, "Q", q) == q, "in-tree gen_ntt q mismatch (mode %s)" % sec
    mr = R.make_mreduce(q, qinv)
    rng = random.Random(0x9E3779B9 ^ q)

    def ref_fwd(a):
        return R.ntt(a, q, zetas, bd, mr)

    def avx_fwd(a):
        return [x % q for x in avx.ntt(a[:])]

    # Intersect candidate ref-slots for each avx-slot over many inputs.
    cand = [set(range(n)) for _ in range(n)]
    for _ in range(200):
        a = [rng.randrange(q) for _ in range(n)]
        ro, ao = ref_fwd(a), avx_fwd(a)
        pos = {}
        for r, v in enumerate(ro):
            pos.setdefault(v, []).append(r)
        for s in range(n):
            cand[s] &= set(pos.get(ao[s], ()))
        if all(len(c) == 1 for c in cand):
            break
    if any(len(c) != 1 for c in cand):
        raise SystemExit("perm underdetermined (mode %s)" % sec)
    P = [next(iter(c)) for c in cand]
    if sorted(P) != list(range(n)):
        raise SystemExit("perm not a bijection (mode %s)" % sec)

    # Validate: avx[s] == ref[P[s]] EXACTLY on fresh inputs.
    vr = random.Random(0xC0FFEE ^ q)
    for _ in range(2000):
        a = [vr.randrange(q) for _ in range(n)]
        ro, ao = ref_fwd(a), avx_fwd(a)
        for s in range(n):
            if ao[s] != ro[P[s]]:
                raise SystemExit("perm validate FAIL mode %s s=%d" % (sec, s))
    return P


def is_bit_permutation(P, n):
    """True iff P is a GF(2)-linear bijection with P(0)=0 -- i.e. it permutes the
    binary digits of the index.  Exactly the permutations a shuffle ladder
    (vperm2i128 / vpunpck / vpblend) can realize; equivalently, avx and ref are
    transpose-equivalent.  A multiplicative-scramble (root mismatch) breaks this."""
    if P[0] != 0:
        return False
    return all(P[a ^ b] == (P[a] ^ P[b]) for a in range(n) for b in range(n))


def main():
    print("gen_ntt_perm.py - avx-vs-ref NTT transpose-equivalence guard\n")
    ok = True
    for sec, (qset, q, n, bd) in MODES.items():
        P = derive_perm(sec, qset, q, n, bd)
        bp = is_bit_permutation(P, n)
        ok &= bp
        print(f"  IRS-{sec} (n={n}): P validated; bit-permutation (nttunpack-realizable, "
              f"avx==transpose(ref)): {'YES' if bp else 'NO -- ROOT MISALIGNED'}")
    if not ok:
        print("\nFAIL: avx NTT is NOT transpose-equivalent to ntt16 -- the avx "
              "generator's primitive root must match gen_zetas.py (smallest "
              "primitive root) so expand_A0_ntt's nttunpack ladder is valid.")
        return 1
    print("\nOK: all modes transpose-equivalent; expand_A0_ntt's nttunpack_avx "
          "ladder is valid (runtime-pinned by test_ntt.c [nttunpack] + the KAT).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
