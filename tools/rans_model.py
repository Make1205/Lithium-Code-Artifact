#!/usr/bin/env python3
"""Shared probability model for the IRS rANS generators.

The model uses the actual symbols emitted by packing.c:
  * z-head is floor((z + 2^(t-1)) / 2^t) for z <- D_{Z,r};
  * hint is the bucket crossing count of HighBits(w) vs HighBits(w - 2z),
    with a uniform bucket phase.

Both PMFs are truncated to |z| <= B_s, matching the signer's norm gate.  The
omitted one-dimensional Gaussian tail is far below the rANS overflow budget, and
no accepted signature can contain a coordinate outside this bound.
"""
import math

PROB_BITS = 10
PROB_SCALE = 1 << PROB_BITS
RANS_N = 2
RANS_FLUSH_BYTES = 4 * RANS_N

# n, ell, m, r, B_s, tau, t per mode (docs/design-notes/param.md).
LEVELS = {
    "120": dict(n=256, ell=2, m=2, r=182, Bs=8090, tau=512, t=7),
    "180": dict(n=256, ell=3, m=3, r=204, Bs=10396, tau=512, t=7),
    "260": dict(n=512, ell=2, m=2, r=257, Bs=15210, tau=1024, t=8),
}


def gaussian_z_pmf(r, Bs):
    """Discrete Gaussian D_{Z,r}, truncated and normalised on [-Bs, Bs]."""
    weights = [math.exp(-(z * z) / (2 * r * r)) for z in range(-Bs, Bs + 1)]
    Z = sum(weights)
    return [(z, w / Z) for z, w in zip(range(-Bs, Bs + 1), weights)]


def zhead_symbol(z, t):
    """Mirror packing.c: head = (z + 2^(t-1)) >> t."""
    return (z + (1 << (t - 1))) // (1 << t)


def zhead_pmf(r, Bs, t):
    M = math.ceil(Bs / (1 << t))
    pmf = {s: 0.0 for s in range(-M, M + 1)}
    for z, p in gaussian_z_pmf(r, Bs):
        pmf[zhead_symbol(z, t)] += p
    return -M, M, pmf


def _add_hint_for_z(out, z, p, tau):
    """Bucket-crossing law for displacement 2z with a uniform bucket phase.

    For |2z| = a*tau + rem, the crossing count has magnitude a with
    probability 1-rem/tau and a+1 with probability rem/tau.  The sign is the
    sign of z.  This is the exact circular bucket-crossing distribution for
    tau-sized buckets; it is the distribution implemented by MakeHint up to the
    two extra residues in Z/(2q), a <2/q effect.
    """
    d = 2 * z
    if d == 0:
        out[0] = out.get(0, 0.0) + p
        return
    sign = 1 if d > 0 else -1
    a, rem = divmod(abs(d), tau)
    out[sign * a] = out.get(sign * a, 0.0) + p * (1.0 - rem / tau)
    if rem:
        out[sign * (a + 1)] = out.get(sign * (a + 1), 0.0) + p * (rem / tau)


def hint_pmf(r, Bs, tau):
    M = math.ceil(2 * Bs / tau) + 1
    pmf = {s: 0.0 for s in range(-M, M + 1)}
    for z, p in gaussian_z_pmf(r, Bs):
        _add_hint_for_z(pmf, z, p, tau)
    s = sum(pmf.values())
    return -M, M, {k: v / s for k, v in pmf.items()}


def quantize(pmf, lo, hi):
    syms = list(range(lo, hi + 1))
    raw = [max(1, round(pmf.get(s, 0.0) * PROB_SCALE)) for s in syms]
    imax = max(range(len(raw)), key=lambda i: raw[i])
    raw[imax] += PROB_SCALE - sum(raw)
    assert raw[imax] >= 1 and sum(raw) == PROB_SCALE
    return syms, raw


def entropy(pmf):
    return -sum(p * math.log2(p) for p in pmf.values() if p > 0)


def cost_moments(freqs, pmf_vals):
    mu = m2 = 0.0
    for g, p in zip(freqs, pmf_vals):
        ell = -math.log2(g / PROB_SCALE)
        mu += p * ell
        m2 += p * ell * ell
    return mu, m2 - mu * mu
