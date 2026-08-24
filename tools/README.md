# IRS Python Tooling

Offline helpers for the IRS reference implementation. Every design-level constant the C code uses is **generated** by a script here **and** the same script doubles as a **bit-exact / numerically-exact golden reference** the C implementation is checked against (the project's "Python full-reproduction" principle).

The rANS design rationale and all numeric targets live in `docs/design-notes/IRS_NIST_rANS.md`; the parameter notes that fix n/q/ell/m/kappa/tau/r/bounds live in `docs/design-notes/param.md`.

## Files

### Generators (patch an `@@AUTOGEN:<tag>@@` region inside a hand-written file)

Each generator rewrites only the text between its `@@AUTOGEN:<tag>@@ BEGIN/END` markers (via `autogen.py`), leaving the hand-written parts intact.  `make tables` runs them all; the avx2 `make check-consts` re-runs and diffs (drift $\Rightarrow$ failure).

| file | patches region in | pins down |
|---|---|---|
| `gen_zetas.py`      | `ntt.h` (`@@ntt@@`)            | NTT twiddles + Montgomery constants; full NTT (base_deg=1) for Lithium-260, incomplete deg-2 NTT for Lithium-120/180. |
| `gen_bounds.py`     | `params.h` (`@@bounds@@`)      | exact integer squared norm thresholds (no truncation; keygen window B_k' <= ||s|| < B_k). |
| `gen_approx_exp.py` | `approx_exp.h` (`@@approx_exp@@`) | sampler minimax exp(-a) constants (B=8, degree-6, target $2^{-43}$) plus wide/small Taylor helpers. |
| `gen_rcdt.py`       | `sampler.c` (`@@rcdt@@`, ref + avx2) | 93-bit RCDT tables (BLISS base sigma2=r/128 and sigma=1 key noise) from `docs/design-notes/RCDT.md`. |
| `gen_ntt_perm.py`   | (guard, emits no file) | per-mode transpose-equivalence guard: derives the avx-to-ref NTT permutation `P` from both Python NTT models (`ntt_ref.py`, `tools/ntt/<qset>/gen_ntt.py`) and asserts it is a **bit-permutation** — i.e. the avx `ntt.S` and the scalar `ntt16` differ only by a shuffle network, so the Pre-NTT `expand_A0_ntt` can lay A0 into the avx layout with the `nttunpack_avx` ladder (no scramble table). |
| `gen_rans_tables.py`| `irs_rans.h` (`@@rans_tables@@`) | the two rANS frequency tables/mode (true rounded z-head, true MakeHint bucket crossing), prob_bits=10; plus the per-slot multiply-by-reciprocal (RCP/RSH/BIAS) replacing the encoder's `x/freq` divide (ryg rANS; self-checked exact). |
| `SigSize.py`        | `irs_rans.h` (`@@rans_meta@@`) | `LITHIUM_RANS_RESERVED_BYTES` for a $2^{-100}$ overflow target (true rounded-z-head/MakeHint PMFs + Bernstein). |
| `irs_ref.py`        | `irs_rans.h` (`@@irs_consts@@`) | exp-recurrence repetition rate M, truncation depth `LITHIUM_REJ_T = 33`, and the conditional-output error gate (`max_rel < $2^{-40}$`); golden ref for `decide`. |

### Pure golden references (no header output)

| file | reference for |
|---|---|
| `ntt_ref.py`    | forward/inverse NTT + basemul (both base degrees), checked element-wise vs C. |
| `reduce_ref.py` | uint16 Montgomery (R=2^16) reduction model + numeric no-overflow / `[0,q)` range proof per modulus. |
| `rans.py`       | byte-exact rANS engine (two-table single stream); C `irs_rans` matches byte-for-byte. |
| `check_rans.py` | reproduces the C `test_rans` deterministic dump from `rans.py`. |
| `rans_model.py` | shared true PMF model for `gen_rans_tables.py` and `SigSize.py`. |


### In-tree support directories

| directory | purpose |
|---|---|
| `ntt/` | AVX2/AVX-512 NTT assembly and constant generators used by optimized-build `check-consts`. |
| `avxexp/` | AVX2 floating-point and AVX-512 IFMA ApproxExp constant generators and local audit helpers. |
| `uarch/` | Instruction-port audit notes and benchmark helper used to explain the recorded optimization campaign. |
| `LogPolyApprox/` | Production scalar ApproxExp polynomial log consumed by `gen_approx_exp.py` and `ref/test/explore/sweep.sh`. |

## Dependencies

Python 3.8+ (no external libraries).

## Quick start

```bash
cd Lithium-Code/ref && make tables       # regenerate every header from the generators above
# or run one generator directly:
cd Lithium-Code/tools
python3 gen_rans_tables.py           # patches irs_rans.h @@rans_tables@@
python3 SigSize.py                   # prints the size model, patches irs_rans.h @@rans_meta@@
python3 rans.py                      # Python reference self-test
```

## Where the frequency tables come from

We do *not* sample any empirical histogram. Both PMFs are evaluated analytically (`gen_rans_tables.py`):

z-head context: the symbol is $\left\lfloor(z+2^{t-1})/2^t\right\rfloor$ for $z \leftarrow D_{\mathbb Z,r}$ truncated by $|z| \le B_s$. Hint context: the symbol is the MakeHint bucket-crossing count for `HighBits(w) - HighBits(w - 2z)` with a uniform bucket phase. For $|2z| = a\tau + b$, it emits magnitude $a$ with probability $1-b/\tau$ and $a+1$ with probability $b/\tau$, carrying the sign of $z$. This is why the PMF includes $|h| \ge 2$ tails.

The vocabulary radius comes from **norm-gating**, not a sigma-multiple truncation: a signature is only emitted when `||z||_2 <= Bs`, so every coefficient obeys `|z_i| <= Bs`, giving the hard bounds `M1 = ceil(Bs / 2^t)` (z-head) and `M2 = ceil(2*Bs / tau) + 1` (hint). The alphabet `[-M, M]` therefore covers every possible coefficient, so out-of-vocabulary failure is mathematically impossible (`p_OOV ≡ 0`).

**All three modes (Lithium-120/180/260) keep two separate tables.** The z-head and hint scales differ by ~2:1 in every mode (sigma_z/sigma_h $\approx$ 1.42/0.71, 1.59/0.80, 1.00/0.50) and are never equal — there is no shared table. The two tables are nonetheless encoded into a single merged rANS byte stream (one two-state flush, one length prefix).

Quantization is the standard "round each `p(s) * 2^prob_bits` to an integer, nudge the largest bucket to absorb the slack" recipe, with **prob_bits = 10** (total frequency 2^10 = 1024). Do not confuse `prob_bits` with the rANS head/tail split `t = floor(log2 r) = 7/7/8`.

## Reservation budget

The single-stream byte reserve in `ref/irs_rans.h` (`LITHIUM_RANS_RESERVED_BYTES` = 481/688/746 for Lithium-120/180/260) comes from `tools/SigSize.py`. It models the merged-stream encoded length over the true emitted symbols: rounded z-head and the full MakeHint bucket-crossing PMF, including $|h| \ge 2$ tails. It sets `R = ceil(mu_b + max(Delta_CLT, Delta_Bern))` for a $2^{-100}$ overflow probability and includes the current two-state rANS flush cost of 8 bytes. The **Bernstein** bound is the binding constraint in all three modes because each symbol's cost is bounded by prob_bits=10 bits and that bounded tail dominates at the $2^{-100}$ quantile. Re-run `SigSize.py` (or `make tables`) after changing any spec parameter.

## rANS reference

`rans.py` is a literal transliteration of the C engine. State is 32-bit, renormalization is byte-wise, and the decoder verifies `x == RANS_L` after the last symbol — any byte-level corruption that survives the renorm loop will overwhelmingly land here. `check_rans.py` reproduces the C `test_rans` deterministic byte dump from it.
