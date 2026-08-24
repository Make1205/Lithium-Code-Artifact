#!/usr/bin/env python3
"""SigSize.py - rANS single-stream length model for IRS.

Emits RANS_RESERVED_BYTES sized for a 2^-100 overflow probability and
prints the compact signature size.

Merged stream of n_b = n(ell+1+m) symbols.  The COST of symbol s is
ell(s) = -log2(g(s)/2^prob_bits), with g the QUANTIZED rANS table.  The
DISTRIBUTION over which that cost is averaged is the TRUE emission PMF p(s):
rounded z-head for z1, and the full MakeHint bucket-crossing distribution for
hint, including |h| >= 2 tails.

    mu_b    = (sum n_i * E_{p_i}[ell_i]) / 8 + 8        (8 = two-state flush)
    sigma_b = sqrt(sum n_i * Var_{p_i}[ell_i]) / 8
    R       = mu_b + max( Phi^{-1}(1-2^-100)*sigma_b , Bernstein(2^-100) )

Bernstein uses ell(s) in [0, prob_bits].  Tables read from the generated rANS block in ref/irs_rans.h
(single source of truth).

Run:  python3 SigSize.py
"""
import os, re, math
import autogen
import rans_model

PROB_BITS = rans_model.PROB_BITS
PROB_SCALE = rans_model.PROB_SCALE
PHI_INV = 11.489  # Phi^{-1}(1 - 2^-100)


def parse_tables(path):
    txt = autogen.read_region(path, "rans_tables")  # only the rANS-tables region
    out = {}
    for sec in ("120", "180", "260"):
        blk = re.search(r"#if LITHIUM_MODE == %s\b(.*?)#endif" % sec, txt, re.S).group(1)
        def get(name):
            fr = [int(x) for x in re.search(r"%s_FREQ\[\d+\] = \{([^}]*)\}" % name, blk).group(1).split(",")]
            lo = int(re.search(r"%s_LO\s+\((-?\d+)\)" % name, blk).group(1))
            return fr, lo
        zf, zlo = get("RANS_ZH")
        hf, hlo = get("RANS_HINT")
        out[sec] = dict(zf=zf, zlo=zlo, hf=hf, hlo=hlo)
    return out


def bernstein_delta(V, p=2.0 ** -100, ell_max=PROB_BITS):
    tgt = -math.log(p)
    b = tgt * (2.0 / 3.0) * ell_max
    c = tgt * 2.0 * V
    x = (b + math.sqrt(b * b + 4 * c)) / 2.0    # x = 8*Delta (bits)
    return x / 8.0


BODY = '''\
/* Single merged rANS stream reservation sized for 2^-100 overflow.  Symbol
 * costs come from the quantized rANS tables; expectations and variances are
 * over the TRUE rounded z-head and MakeHint bucket-crossing PMFs.
{summary}
 */
{blocks}'''


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    tabs = parse_tables(os.path.normpath(os.path.join(here, "..", "ref", "irs_rans.h")))
    print("SigSize.py - IRS rANS reservation (2^-100 overflow, true PMFs)\n")
    blocks, summary = [], []
    for sec, gg in rans_model.LEVELS.items():
        n, ell, m, r, Bs, tau, t = (gg["n"], gg["ell"], gg["m"], gg["r"],
                                    gg["Bs"], gg["tau"], gg["t"])
        T = tabs[sec]
        zlo, zhi, zdist = rans_model.zhead_pmf(r, Bs, t)
        hlo, hhi, hdist = rans_model.hint_pmf(r, Bs, tau)
        assert zlo == T["zlo"] and len(T["zf"]) == zhi - zlo + 1
        assert hlo == T["hlo"] and len(T["hf"]) == hhi - hlo + 1
        zpmf = [zdist[zlo + i] for i in range(len(T["zf"]))]
        hpmf = [hdist[hlo + i] for i in range(len(T["hf"]))]
        mu_z, var_z = rans_model.cost_moments(T["zf"], zpmf)
        mu_h, var_h = rans_model.cost_moments(T["hf"], hpmf)
        n_z = n * (ell + 1); n_h = n * m; nb = n_z + n_h
        mean_bits = n_z * mu_z + n_h * mu_h
        var_bits = n_z * var_z + n_h * var_h
        mu_b = mean_bits / 8.0 + rans_model.RANS_FLUSH_BYTES
        sigma_b = math.sqrt(var_bits) / 8.0
        d_clt = PHI_INV * sigma_b
        d_bern = bernstein_delta(var_bits)
        tail_bits = max(d_clt, d_bern)
        R = math.ceil(mu_b + max(d_clt, d_bern))
        bound = math.exp(-(tail_bits * 8) ** 2 /
                         (2 * var_bits + (2.0 / 3.0) * PROB_BITS * tail_bits * 8))
        challenge = (n + 7) // 8
        zlo = (ell + 1) * ((n * t + 7) // 8)
        sig = challenge + (2 + R) + zlo
        print(f"--- IRS-{sec} ---")
        print(f"  z-head E[l]={mu_z:.3f} sd={math.sqrt(var_z):.3f} bit | "
              f"hint E[l]={mu_h:.3f} sd={math.sqrt(var_h):.3f} bit")
        print(f"  n_b={nb}  mu_b={mu_b:.1f}B sigma_b={sigma_b:.2f}B  "
              f"Delta_CLT={d_clt:.1f} Delta_Bern={d_bern:.1f}  -> R={R}B")
        print(f"  Bernstein bound at R: Pr[len > R] <= 2^{math.log2(bound):.1f}")
        print(f"  est compact sig = challenge({challenge}) + rANS_block({2+R}) + z_lo({zlo}) = {sig}B")
        summary.append(f" *   IRS-{sec}: R={R} (mu_b={mu_b:.0f}, sigma_b={sigma_b:.1f}, "
                       f"Delta_Bern={d_bern:.0f}); est sig={sig}B")
        blocks.append(f"#if LITHIUM_MODE == {sec}\n#define RANS_RESERVED_BYTES {R}\n#endif")
    body = BODY.format(summary="\n".join(summary), blocks="\n".join(blocks))
    autogen.patch_region(os.path.normpath(os.path.join(here, "..", "ref", "irs_rans.h")),
                         "rans_meta", body)
    print("\nPatched @@AUTOGEN:rans_meta@@ region in ref/irs_rans.h")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
