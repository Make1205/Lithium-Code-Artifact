#!/usr/bin/env python3
"""gen_rans_tables.py - generate the rANS frequency tables for IRS.

z-head table ZH: the actual high part emitted by packing.c,
    high = floor((z + 2^(t-1)) / 2^t), z ~ D_{Z,r}, |z| <= B_s.

hint table HINT: the actual MakeHint bucket-crossing distribution induced by
    h = HighBits(w) - HighBits(w - 2z), with uniform bucket phase and
    z ~ D_{Z,r}, |z| <= B_s.  This includes |h| >= 2 tails.

Quantization: f(s) = max(1, round(PMF(s)*2^PROB_BITS)); the largest bucket is
adjusted so each table sums to exactly 2^PROB_BITS = 1024.

Run:  python3 gen_rans_tables.py
"""
import os, math, random
import autogen
import rans_model

PROB_BITS = rans_model.PROB_BITS
PROB_SCALE = rans_model.PROB_SCALE
RANS_L = 1 << 23          # encoder lower bound (matches irs_rans.c)


# ---- multiply-by-reciprocal for the encoder's x/freq (ryg rANS) ----
# The encoder replaces the per-symbol  x_new = (x/freq)<<PB + x%freq + start
# divide with  q = (x*RCP >> 32) >> RSH ;  x_new = x + BIAS + q*(PSCALE-freq).
# This is the Granlund-Montgomery "round-up reciprocal": with
#   shift = ceil(log2 freq),  RCP = ceil(2^(shift+31) / freq),  RSH = shift-1,
# q == floor(x/freq) EXACTLY for all x in [0, 2^31).  The encoder's x is always
# < x_max = (RANS_L>>PB<<8)*freq <= 2^31, so the substitution is bit-exact (the
# rANS byte stream -- hence the KAT -- is unchanged).  freq==1 is the ryg special
# case (RCP=~0, RSH=0, BIAS=start+PSCALE-1).  Verified per slot below.
def rans_sym_init(start, freq):
    assert 1 <= freq <= PROB_SCALE
    if freq < 2:
        return 0xFFFFFFFF, 0, start + PROB_SCALE - 1
    shift = 0
    while freq > (1 << shift):
        shift += 1                                  # ceil(log2 freq)
    rcp = ((1 << (shift + 31)) + freq - 1) // freq  # ceil(2^(shift+31)/freq)
    assert 0 <= rcp <= 0xFFFFFFFF
    return rcp, shift - 1, start


def _enc_ref(x, start, freq):                       # plain-division reference
    return ((x // freq) << PROB_BITS) + (x % freq) + start


def _enc_rcp(x, start, freq, rcp, rsh, bias):       # reciprocal (mirrors C, u32)
    q = (((x * rcp) >> 32) & 0xFFFFFFFF) >> rsh
    return (x + bias + q * (PROB_SCALE - freq)) & 0xFFFFFFFF


def verify_sym(start, freq, rcp, rsh, bias):
    """Double-check: the reciprocal state-update equals the plain divide over a
    thorough sample of the encoder's x range [1, x_max) (small x in full, both
    edges of each floor step, and random x).  GM proves exactness; this guards
    against any transcription/formula error."""
    x_max = ((RANS_L >> PROB_BITS) << 8) * freq     # == 2^21 * freq
    kmax = x_max // freq                            # == 2^21
    rng = random.Random(0x9E3779B9 ^ (freq << 8) ^ (start & 0xFFFF))
    ks = set(range(0, min(kmax, 2048) + 1)) | {kmax}
    ks |= {rng.randrange(0, kmax + 1) for _ in range(4000)}
    xs = {1, 2, x_max - 1}
    for k in ks:
        for x in (k * freq - 1, k * freq, k * freq + 1):
            if 1 <= x < x_max:
                xs.add(x)
    xs |= {rng.randrange(1, x_max) for _ in range(4000)}
    for x in xs:
        if _enc_rcp(x, start, freq, rcp, rsh, bias) != _enc_ref(x, start, freq):
            raise SystemExit("rcp verify FAIL: freq=%d start=%d x=%d" %
                             (freq, start, x))

def zhead_table(r, Bs, t):
    sigma = r / (1 << t)
    lo, hi, pmf = rans_model.zhead_pmf(r, Bs, t)
    syms, freqs = rans_model.quantize(pmf, lo, hi)
    return sigma, max(abs(lo), abs(hi)), pmf, syms, freqs


def hint_table(r, Bs, tau):
    sigma_h = 2 * r / tau
    lo, hi, pmf = rans_model.hint_pmf(r, Bs, tau)
    syms, freqs = rans_model.quantize(pmf, lo, hi)
    return sigma_h, max(abs(lo), abs(hi)), pmf, syms, freqs


BODY = '''\
/* Two rANS frequency tables per mode (prob_bits=10, sum=1024; see tools/gen_rans_tables.py):
 *   RANS_ZH_*  : true rounded z-head PMF, z ~ D_Z,r, |z| <= B_s
 *   RANS_HINT_*: true bucket-crossing PMF for HighBits(w)-HighBits(w-2z)
 * Alphabets are contiguous [SYM_LO, SYM_LO + N - 1]; CDF has N+1 prefix sums.
{summary}
 */
#define RANS_PROB_BITS 10

{blocks}'''


def _fmt_rows(vals, per=32):
    lines = []
    for i in range(0, len(vals), per):
        lines.append("    " + ", ".join(str(v) for v in vals[i:i + per]) + ",")
    return "\n".join(lines)


def emit_table(name, syms, freqs):
    lo = syms[0]; n = len(syms)
    cdf = [0]
    for f in freqs:
        cdf.append(cdf[-1] + f)
    rcp, rsh, bias = [], [], []
    for start, f in zip(cdf, freqs):                 # start[slot] = cdf[slot]
        rc, sh, bi = rans_sym_init(start, f)
        verify_sym(start, f, rc, sh, bi)
        rcp.append(rc); rsh.append(sh); bias.append(bi)
    # Decoder slot lookup: SLOT[v] = the slot s with CDF[s] <= v < CDF[s+1],
    # for every v in [0, 2^PROB_BITS).  Replaces the O(n) linear scan with one
    # load.  n <= 256 so uint8 suffices; KAT-neutral (identical decode result).
    assert n <= 256, f"{name}: alphabet {n} > 256, SLOT needs wider type"
    slot = []
    s = 0
    for v in range(PROB_SCALE):
        while not (cdf[s] <= v < cdf[s + 1]):
            s += 1
        slot.append(s)
    assert slot[0] == 0 and slot[-1] == n - 1
    fb = ", ".join(str(f) for f in freqs)
    cb = ", ".join(str(c) for c in cdf)
    rb = ", ".join("%uu" % v for v in rcp)
    sb = ", ".join(str(v) for v in rsh)
    bb = ", ".join("%uu" % v for v in bias)
    return (f"#define {name}_LO   ({lo})\n"
            f"#define {name}_N    ({n})\n"
            f"static const uint16_t {name}_FREQ[{n}] = {{ {fb} }};\n"
            f"static const uint16_t {name}_CDF[{n + 1}] = {{ {cb} }};\n"
            f"/* encoder x/freq via multiply-by-reciprocal (ryg rANS; exact, "
            f"see gen_rans_tables.py) */\n"
            f"static const uint32_t {name}_RCP[{n}] = {{ {rb} }};\n"
            f"static const uint8_t {name}_RSH[{n}] = {{ {sb} }};\n"
            f"static const uint32_t {name}_BIAS[{n}] = {{ {bb} }};\n"
            f"/* decoder val->slot direct lookup ({PROB_SCALE} entries; "
            f"O(1), replaces linear CDF scan) */\n"
            f"static const uint8_t {name}_SLOT[{PROB_SCALE}] = {{\n"
            f"{_fmt_rows(slot)}\n}};\n")


def main():
    print("gen_rans_tables.py - IRS two-table rANS frequency tables\n")
    blocks, summary = [], []
    for sec, p in rans_model.LEVELS.items():
        sig_zh, M1, zpmf, zsyms, zfreq = zhead_table(p["r"], p["Bs"], p["t"])
        sig_h, M2, hpmf, hsyms, hfreq = hint_table(p["r"], p["Bs"], p["tau"])
        H_zh = rans_model.entropy(zpmf)
        H_hint = rans_model.entropy(hpmf)
        tail2 = sum(v for k, v in hpmf.items() if abs(k) >= 2)
        tail3 = sum(v for k, v in hpmf.items() if abs(k) >= 3)
        print(f"--- IRS-{sec} ---")
        print(f"  z-head: sigma={sig_zh:.4f} M1={M1} alphabet=[{-M1},{M1}] "
              f"({2*M1+1} syms) H={H_zh:.3f} bit/coef")
        print(f"  hint:   sigma={sig_h:.4f} M2={M2} alphabet=[{-M2},{M2}] "
              f"({2*M2+1} syms) H={H_hint:.3f} bit/coef "
              f"P(|h|>=2)={tail2:.4e} P(|h|>=3)={tail3:.4e}")
        summary.append(f" *   IRS-{sec}: z-head sigma={sig_zh:.3f} M1={M1} H={H_zh:.3f}; "
                       f"hint sigma={sig_h:.3f} M2={M2} H={H_hint:.3f}, "
                       f"P(|h|>=2)={tail2:.3e}")
        block = (f"#if LITHIUM_MODE == {sec}\n"
                 + emit_table("RANS_ZH", zsyms, zfreq)
                 + emit_table("RANS_HINT", hsyms, hfreq)
                 + "#endif\n")
        blocks.append(block)
    body = BODY.format(summary="\n".join(summary), blocks="\n".join(blocks))
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.normpath(os.path.join(here, "..", "ref", "irs_rans.h"))
    autogen.patch_region(path, "rans_tables", body)
    print(f"\nPatched @@AUTOGEN:rans_tables@@ region in {path}")
    print("  reciprocal tables (RCP/RSH/BIAS) self-check: PASS (x/freq exact)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
