#!/usr/bin/env python3
"""
gen_ntt.py -- verified model + .S/constant generator for an INCOMPLETE n=256
negacyclic NTT over Z_q[x]/(x^256+1), q=30977, 16-bit AVX2, Kyber-style level
merging + 3-tier shuffle network, UNSIGNED [0,q) + LAZY REDUCTION.

PARAM B (q=30977, n=256):
  q is prime, q-1 = 2^8 * 11^2, so only n | q-1 (2n = 512 does NOT divide q-1).
  => the cyclotomic x^256+1 splits only into degree-2 factors:
        x^256 + 1 = prod_{i<128} (x^2 - gamma_i),   gamma_i = psi^(2*brv7(i)+1)
  => INCOMPLETE NTT: 7 butterfly levels (len 128,64,...,2; STOP at len=2), 128
     zetas, bit-reversal width 7, and the pointwise step is a BASEMUL of 128
     coefficient-pairs (each a product mod x^2 - gamma_i), NOT a coefficient-wise
     multiply.  Final inverse scale = (2^7)^{-1} = 128^{-1} mod q.

WHY UNSIGNED + LAZY (this is the rewrite from the old SIGNED+Barrett path):
  q = 30977 < 2^15, BUT 2q = 61954 fits uint16 (< 2^16 = 65536) while NOT int16
  (> 32767).  So the UNSIGNED [0,q) representation has ~2x the dynamic range of
  the signed-centered one for THIS q.  Consequences:
    * The OLD signed montmul output was (-q,q); a+t reached +-2q = 61952 > 32767
      and OVERFLOWED int16, forcing a Barrett red16 (3 mults) on EVERY level.
    * The UNSIGNED montmul (q57601-style: vpmullw + vpmulhuw + vpmulhuw + borrow
      correction) outputs [0,q) DIRECTLY -- reduction folded in, NO separate
      Barrett.  montmul = 3 multiplies total (vs signed montmul 3 + red16 3 = 6).
    * a+t in [0,2q) FITS uint16, so a "lazy" level emits a PLAIN vpaddw (a value
      in [0,q) plus t in [0,q) -> [0,2q)).  When the running bound would exceed
      uint16 the next level does a MULTIPLY-FREE conditional subtract (csub:
      s-q via vpsubw + sign-mask add-back), 0 multiplies.

LAZY SCHEDULE: derived by an EXACT per-register bound tracker (run_*_half tracks
  each register's worst-case max value through every butterfly and shuffle).  A
  butterfly is emitted LAZY unless a lazy op would (a) overflow uint16 on the add
  or (b) feed a montmul an input >= 2^16 or with input*zeta >= 2^16*q; in that
  case it is emitted CSUB.  The tracker DERIVES the forward (CT) and inverse (GS)
  schedules INDEPENDENTLY -- they differ because GS reduces one output for free
  (its b-output passes through montmul).

Pipeline (mirrors the proven q=57601 UNSIGNED template):
  1. Build the forward transform as an explicit op-list (cross-register vertical
     butterflies for len>=16, shuffle+butterfly for len<16, STOPPING at len=2),
     baking the per-lane Montgomery zeta vectors by coefficient-index tracking.
  2. The inverse is the LITERAL REVERSE of that program (CT->GS, shuffles are
     involutions and replay unchanged), using z^{-1} twiddles.
  3. The basemul (mod x^2 - gamma_i) reads the adjacent even/odd lanes the len=2
     butterfly leaves (q57601 layout: even lane = coef0, odd lane = coef1).
  4. Verify ALL of it with a pure-Python AVX2 lane-permutation interpreter that
     ALSO tracks exact integer bounds, confirming the lazy schedule never leaves
     [0,2^16) and keeps every montmul input valid.  Then emit ntt.S + ntt_consts.c.
"""
import sys

Q = 30977
R = 1 << 16
MONT = R % Q                       # 3582   = 2^16 mod q
QINV = pow(Q, -1, R)               # 34561  = q^-1 mod 2^16
RMQ  = R - Q                       # 34559  = 2^16 - q  (borrow-correction const)
def _factorize(m):
    f = []; d = 2
    while d * d <= m:
        if m % d == 0:
            f.append(d)
            while m % d == 0: m //= d
        d += 1
    if m > 1: f.append(m)
    return f
def _mod_order(g, q):                # multiplicative order of g mod q
    n = q - 1
    for p in _factorize(n):
        while n % p == 0 and pow(g, n // p, q) == 1: n //= p
    return n
def _primitive_root(q, order):       # SMALLEST primitive order-th root (matches gen_zetas.py)
    for g in range(2, q):
        if _mod_order(g, q) == order: return g
    raise ValueError("no primitive %d-th root mod %d" % (order, q))
# Use ref's root convention (gen_zetas.py: smallest primitive root) so the avx
# NTT is a true transpose of the scalar ntt16 -- this lets the A0 layout
# permutation be a clean shuffle ladder (nttunpack) instead of a scramble.
PSI = _primitive_root(Q, 256)      # 430: smallest primitive 256th root, psi^128 == -1
RINV = pow(R, -1, Q)
LEVELS = 7                         # incomplete: stop forward at len=2
SCALE = pow(1 << LEVELS, -1, Q)    # (2^7)^-1 = 128^-1 mod q = 30735
NINV  = pow(256, -1, Q)            # 256^-1 mod q = 30856 (informational; NOT used)
R2_MOD_Q  = MONT * MONT % Q        # 6246
SCALE_TOMONT = SCALE * R2_MOD_Q % Q # = 6341  (final-scale to-Montgomery multiplier)

BRV = 7                            # bit-reversal width = #levels (incomplete)
def brv(x, w=BRV):
    r = 0
    for i in range(w):
        r |= ((x >> i) & 1) << (w - 1 - i)
    return r

# Forward zeta[k] (k=1..127) in the standard CT loop order is psi^brv7(k); store
# Montgomery.  k indexes the (len,start) pairs exactly like the scalar loop.
ZETA_MONT = [pow(PSI, brv(k), Q) * MONT % Q for k in range(128)]

def fqmul(a, b):
    """plain modular multiply with one Montgomery factor removed (unsigned)."""
    return (a * b * RINV) % Q

# (len,start) -> global zeta index, standard CT loop order, len 128..2.
FWD_K = {}
_k = 1; _len = 128
while _len >= 2:
    _s = 0
    while _s < 256:
        FWD_K[(_len, _s)] = _k; _k += 1; _s += 2 * _len
    _len >>= 1

def zeta_of(length, coeff_low):
    start = (coeff_low // (2 * length)) * (2 * length)
    return ZETA_MONT[FWD_K[(length, start)]]

# basemul twiddle: pair i (coeffs 2i,2i+1) lives mod (x^2 - gamma_i),
# gamma_i = psi^(2*brv7(i)+1) mod q.
def gamma_of(i):
    return pow(PSI, 2 * brv(i) + 1, Q)

# --------------------------------------------------------------------------
# AVX2 lane permutations, modelled INSTRUCTION-ACCURATELY (16 int16 lanes,
# lanes 0-7 = low 128-bit, 8-15 = high 128-bit).  Each shuffleK returns exactly
# what the corresponding asm macro computes.  Lanes carry (value, bound) cells
# in the interpreter; the permutations move whole cells, so bounds follow.
# --------------------------------------------------------------------------
def _dwords(x):  return [x[2*i:2*i+2] for i in range(8)]
def _undw(dws):  return [v for dw in dws for v in dw]
def vperm2i128(a, b, imm):
    lo = {0: a[0:8], 1: a[8:16], 2: b[0:8], 3: b[8:16]}
    return lo[imm & 3] + lo[(imm >> 4) & 3]
def vpunpcklqdq(a, b): return a[0:4] + b[0:4] + a[8:12] + b[8:12]
def vpunpckhqdq(a, b): return a[4:8] + b[4:8] + a[12:16] + b[12:16]
def vmovsldup(x):
    d = _dwords(x); return _undw([d[0],d[0],d[2],d[2],d[4],d[4],d[6],d[6]])
def vpsrlq(x, n):
    assert n == 32
    d = _dwords(x); z = [(0, 0), (0, 0)]
    return _undw([d[1],z,d[3],z,d[5],z,d[7],z])
def vpblendd(a, b, mask):
    da, db = _dwords(a), _dwords(b)
    return _undw([db[i] if (mask>>i)&1 else da[i] for i in range(8)])

def shuffle8(a, b):
    return vperm2i128(a, b, 0x20), vperm2i128(a, b, 0x31)
def shuffle4(a, b):
    return vpunpcklqdq(a, b), vpunpckhqdq(a, b)
def shuffle2(a, b):
    c = vpblendd(a, vmovsldup(b), 0xAA)
    d = vpblendd(vpsrlq(a, 32), b, 0xAA)
    return c, d
SHUF = {8: shuffle8, 4: shuffle4, 2: shuffle2}

# --------------------------------------------------------------------------
# Build the forward op-list for one half (base = 0 or 128), baking zeta vectors
# by tracking coefficient indices through the exact shuffles.  STOPS at len=2
# (incomplete), so the deepest in-register tier is shuf2 (no shuf1).
# Ops:  ('xb', ra, rb, zvec, lvl)   cross-register vertical CT butterfly
#       ('sh', k, r2m, r2m1, lvl)   in-place shuffle of a register pair
#       ('ib', r2m, r2m1, zvec, lvl) vertical CT butterfly after a shuffle
# --------------------------------------------------------------------------
LVL = {64: 1, 32: 2, 16: 3, 8: 4, 4: 5, 2: 6}   # len -> level (level 0 = len128 cross-half)

def build_half(base):
    idx = [[base + 16*r + l for l in range(16)] for r in range(8)]
    prog = []
    for length in (64, 32, 16):                  # cross-register
        step = length // 16
        for r in range(8):
            if (r // step) % 2 == 0:
                p = r + step
                z = zeta_of(length, idx[r][0])
                prog.append(('xb', r, p, [z]*16, LVL[length]))
    for length in (8, 4, 2):                      # in-register, stop at len=2
        sh = SHUF[length]
        for m in range(4):
            a, b = idx[2*m], idx[2*m+1]
            c, d = sh(a, b)
            zvec = [zeta_of(length, c[l]) for l in range(16)]
            prog.append(('sh', length, 2*m, 2*m+1, LVL[length]))
            prog.append(('ib', 2*m, 2*m+1, zvec, LVL[length]))
            idx[2*m], idx[2*m+1] = c, d
    return prog, idx                              # idx = final coeff layout per reg

FWD0, IDX0 = build_half(0)
FWD1, IDX1 = build_half(128)
Z0   = zeta_of(128, 0)                            # level-0 (cross-half) zeta
FINAL_IDX = [IDX0, IDX1]

# --------------------------------------------------------------------------
# BIT-EXACT 16-bit AVX2 UNSIGNED arithmetic, mirroring ntt.S instruction by
# instruction (q57601 template), so the Python proof is identical to what the
# assembly computes.  A "cell" is (uint16_value, max_bound); we track an exact
# worst-case bound alongside each value so the lazy schedule is PROVEN to keep
# every value < 2^16 and every montmul input valid (< 2^16 and *zeta < 2^16*q).
# --------------------------------------------------------------------------
def u16(v):  return v & 0xFFFF
def vpmullw(a, b):  return u16(a * b)
def vpmulhuw(a, b): return (u16(a) * u16(b)) >> 16        # UNSIGNED high word

def montmul_u(inp, zeta):
    """ASM montmul (UNSIGNED): out in [0,q).  zl = zeta*qinv mod 2^16.
      m=in*zl mod 2^16; out=mulhi_u(in,zeta)-mulhi_u(m,q); +RMQ on borrow."""
    inp = u16(inp); zeta = u16(zeta)
    zl  = u16(zeta * QINV)
    m   = vpmullw(inp, zl)
    out = vpmulhuw(inp, zeta)
    mq  = vpmulhuw(m, Q)
    # bit-exact borrow correction (mirrors ntt.S montmul):
    #   out = out - mq (wraps mod 2^16); if mq > out the subtract underflowed, so
    #   subtract RMQ (= 2^16 - q): the wrap already added 2^16, net effect is +q.
    res = u16(out - mq)
    if mq > out:                       # underflow -> -RMQ keeps result in [0,q)
        res = u16(res - RMQ)
    return res

# montmul input validity caps (PROVEN bounds the tracker enforces).
MM_CAP   = R                            # input must be < 2^16
MMZ_CAP  = R * Q                        # input*zeta must be < 2^16 * q

def montmul_valid(bound, zmax=Q-1):
    return bound < MM_CAP and bound * zmax < MMZ_CAP

def csub_u(s):
    """conditional subtract: s in [0,2q) -> [0,q).  s-q; if borrow add q back."""
    assert 0 <= s < 2*Q
    t = s - Q
    return t if t >= 0 else t + Q

# CT butterfly (UNSIGNED).  out a' = a+t, b' = a+(q-t), t = montmul(b,zeta).
#   LAZY  : inputs in [0,q); outputs in [0,2q) (no csub).
#   CSUB  : inputs in [0,2q); csub a first, then csub both outputs -> [0,q).
def ct_lazy(a, t):
    va, ba = a; vt, bt = t
    assert ba + Q <= R, "lazy CT add overflow: a-bound %d + q" % ba
    qt = Q - vt
    return (u16(va + vt), ba + (Q - 1)), (u16(va + qt), ba + (Q - 1))

def ct_csub(a, t):
    va, ba = a; vt, bt = t
    var = csub_u(va) if ba >= Q else va     # reduce input a to [0,q)
    qt = Q - vt
    return (csub_u(var + vt), Q - 1), (csub_u(var + qt), Q - 1)

def mul_t(z, b):
    """t = montmul(b, z), bit-exact.  Asserts the input is a valid montmul input.
    Output is in [0,q) so bound = Q-1."""
    vb, bb = b
    assert montmul_valid(bb), "montmul input out of range: bound %d" % bb
    return (montmul_u(vb, u16(z)), Q - 1)

def run_fwd_half(regs, prog=None, fixed=None):
    """Run a FWD0/FWD1 program on 8 (value,bound) registers.  If `fixed` is None,
    DERIVE the csub-vs-lazy choice per butterfly from the live bounds and return
    the derived schedule; if `fixed` is a schedule list, REPLAY it exactly (so the
    model matches the emitted asm bit-for-bit).  Both halves are symmetric, so the
    schedule derived from FWD0 is the one emitted for both."""
    if prog is None: prog = FWD0
    sched = []; fi = 0
    for op in prog:
        if op[0] in ('xb', 'ib'):
            _, ra, rb, z, lvl = op
            if fixed is not None:
                do_csub = fixed[fi][1]; fi += 1
            else:
                # decide from the WORST (max) lane bound of the operand register
                # (lanes can diverge after shuffles).
                ba = max(c[1] for c in regs[ra])
                bb = max(c[1] for c in regs[rb])
                # lazy OK iff a+t stays < 2^16 AND the montmul input b is valid.
                do_csub = not ((ba + Q <= R) and montmul_valid(bb))
            sched.append((lvl, do_csub))
            for l in range(16):
                t = mul_t(z[l], regs[rb][l])
                if do_csub:
                    regs[ra][l], regs[rb][l] = ct_csub(regs[ra][l], t)
                else:
                    regs[ra][l], regs[rb][l] = ct_lazy(regs[ra][l], t)
        else:  # 'sh'
            _, k, r0, r1, _ = op
            regs[r0], regs[r1] = SHUF[k](regs[r0], regs[r1])
    return regs, sched

def zinv_mont(zm):
    plain = zm * RINV % Q
    return pow(plain, -1, Q) * MONT % Q

# GS butterfly (UNSIGNED).  s = u+v ; d = u+(q-v) ; b = montmul(d, zinv).
#   The b-output ALWAYS passes through montmul -> always [0,q) (reduced free).
#   Only s = u+v can be lazy.
#   LAZY : inputs [0,q); s in [0,2q), b in [0,q).
#   CSUB : inputs may be [0,2q); csub u,v first, then csub s -> [0,q).
def gs_lazy(u, v, zi):
    vu, bu = u; vv, bv = v
    assert bu + Q <= R, "lazy GS add overflow: u-bound %d + q" % bu
    s = (u16(vu + vv), bu + (Q - 1))
    d = (u16(vu + (Q - vv)), bu + (Q - 1))
    assert montmul_valid(d[1]), "GS montmul input out of range: %d" % d[1]
    b = (montmul_u(d[0], u16(zi)), Q - 1)
    return s, b

def gs_csub(u, v, zi):
    vu, bu = u; vv, bv = v
    vur = csub_u(vu) if bu >= Q else vu
    vvr = csub_u(vv) if bv >= Q else vv
    s = (csub_u(vur + vvr), Q - 1)
    d = (u16(vur + (Q - vvr)), Q - 1)
    b = (montmul_u(d[0], u16(zi)), Q - 1)
    return s, b

def run_inv_half(regs, prog=None, fixed=None):
    if prog is None: prog = FWD0
    sched = []; fi = 0
    for op in reversed(prog):
        if op[0] in ('xb', 'ib'):
            _, r0, r1, z, lvl = op
            if fixed is not None:
                do_csub = fixed[fi][1]; fi += 1
            else:
                # WORST (max) lane bounds of both operand registers (lanes diverge
                # after shuffles + GS reducing only the b-output).
                bu = max(c[1] for c in regs[r0]); bv = max(c[1] for c in regs[r1])
                # lazy GS needs s = u+v < 2^16 (bu + bv <= 2^16) AND the montmul
                # input d = u+(q-v) be valid (bound <= bu+q, d<2^16, d*zeta<2^16*q).
                d_bound = bu + Q                    # u + (q - v) <= bu + q
                do_csub = not ((bu + bv <= R) and (d_bound <= R) and montmul_valid(d_bound))
            sched.append((lvl, do_csub))
            for l in range(16):
                zi = zinv_mont(z[l])
                if do_csub:
                    regs[r0][l], regs[r1][l] = gs_csub(regs[r0][l], regs[r1][l], zi)
                else:
                    regs[r0][l], regs[r1][l] = gs_lazy(regs[r0][l], regs[r1][l], zi)
        else:
            _, k, r0, r1, _ = op
            regs[r0], regs[r1] = SHUF[k](regs[r0], regs[r1])
    return regs, sched

# Derive the schedules ONCE from a fully-reduced [0,q) input (the worst case for
# bound growth is the bound, not the value, so a symbolic max-bound run suffices).
def _bound_regs(base, init_bound):
    return [[(0, init_bound) for _ in range(16)] for _ in range(8)]

# --- forward schedule derivation (level 0 happens outside the half program) ---
# level 0 (cross-half): inputs [0,q), lazy add -> [0,2q).  No montmul-input issue
# (q-1 < 2^16).  Always lazy (first level, inputs reduced).
FWD_L0_CSUB = False    # level 0 lazy: out [0,2q)
_fwd_regs, FWD_SCHED = run_fwd_half(_bound_regs(0, 2*Q - 1))   # half sees [0,2q) from L0
# regs after a full forward half -> per-register output bound (feeds basemul)
FWD_OUT_BOUND = max(c[1] for reg in _fwd_regs for c in reg)

# --- inverse schedule derivation ---
# The inverse half is fed the FORWARD output, whose worst-case per-register bound
# is FWD_OUT_BOUND = 2q-1 (the level-6 lazy output).  (When called after basemul
# the input is [0,q), a strict sub-case.)  We therefore DERIVE the inverse csub
# schedule from the 2q-1 start bound, NOT q-1 -- the GS csub butterfly handles a
# [0,2q) input by csub'ing u,v first.  After the half, level 0 (cross-half GS) and
# final scaling run.
INV_IN_BOUND = FWD_OUT_BOUND                 # = 2q-1
_inv_regs, INV_SCHED = run_inv_half(_bound_regs(0, INV_IN_BOUND))
INV_OUT_BOUND = max(c[1] for reg in _inv_regs for c in reg)

# convenience: which forward / inverse LEVELS are csub vs lazy (for docs).
def levels_csub(sched):
    d = {}
    for lvl, c in sched:
        d.setdefault(lvl, c)
        d[lvl] = d[lvl] or c
    return d
FWD_LEVEL_CSUB = levels_csub(FWD_SCHED); FWD_LEVEL_CSUB[0] = FWD_L0_CSUB
INV_LEVEL_CSUB = levels_csub(INV_SCHED)

# --------------------------------------------------------------------------
# Full transform model (value + bound), used by verify().
# --------------------------------------------------------------------------
def ntt(poly, raw=False):
    """Forward NTT, bit-exact to ntt_avx.  Returns the RAW lane uint16 values (in
    [0,2q), as the level-6 lazy asm leaves them) when raw=True, else reduced %Q.
    The raw form is what the asm leaves; intt()/basemul() accept it directly."""
    a = [(x % Q, Q - 1) for x in poly]
    # level 0 (cross-half): lazy CT (inputs [0,q) -> [0,2q))
    for j in range(128):
        t = mul_t(Z0, a[j+128])
        if FWD_L0_CSUB:
            a[j], a[j+128] = ct_csub(a[j], t)
        else:
            a[j], a[j+128] = ct_lazy(a[j], t)
    out = [c[0] for c in a]
    for h, prog in ((0, FWD0), (1, FWD1)):
        regs = [[a[128*h+16*r+l] for l in range(16)] for r in range(8)]
        regs, _ = run_fwd_half(regs, prog, fixed=FWD_SCHED)
        for r in range(8):
            for l in range(16):
                out[128*h+16*r+l] = regs[r][l][0]
    return out if raw else [x % Q for x in out]

def intt(poly):
    # mirror the asm: the inverse does NOT pre-reduce -- it receives the forward
    # output in [0,2q) (level-6 lazy) directly.  Track each input's true bound so
    # the derived INV_SCHED applies.  (basemul output [0,q) is a strict sub-case.)
    a = [(x % R, (Q - 1) if (x % R) < Q else (2*Q - 1)) for x in poly]
    for h, prog in ((0, FWD0), (1, FWD1)):
        regs = [[a[128*h+16*r+l] for l in range(16)] for r in range(8)]
        regs, _ = run_inv_half(regs, prog, fixed=INV_SCHED)
        for r in range(8):
            for l in range(16):
                a[128*h+16*r+l] = regs[r][l]
    zi0 = zinv_mont(Z0)
    # level 0 cross-half GS.  inputs may be [0,2q) (from a lazy GS s-output);
    # use csub form to be safe (derived: the cross-half level reduces).
    for j in range(128):
        u, v = a[j], a[j+128]
        bu = u[1]
        if bu + v[1] <= R and bu + Q <= R:
            a[j], a[j+128] = gs_lazy(u, v, zi0)
        else:
            a[j], a[j+128] = gs_csub(u, v, zi0)
    # final *1/2^7 and to-Montgomery: montmul by SCALE_TOMONT; this makes
    # intt_tomont(ntt(a)) = a*R, while product pipelines return normal values.
    out = []
    for x in a:
        xr = csub_u(x[0]) if x[1] >= Q else x[0]
        out.append(montmul_u(xr, SCALE_TOMONT) % Q)
    return out

# level-0 cross-half GS schedule (for asm emission): decide from a [0,q) start
# (inverse halves leave their cross-half operands... derived live below).
def _inv_l0_csub():
    # after a full inverse half the cross-half operands have INV_OUT_BOUND each.
    bu = INV_OUT_BOUND
    return not (bu + INV_OUT_BOUND <= R and bu + Q <= R)
INV_L0_CSUB = _inv_l0_csub()

# --------------------------------------------------------------------------
# BASEMUL.  After 7 levels the storage layout (FINAL_IDX) places, in each data
# register, 8 quadratic pairs: even lane l = coef0 of a pair, odd lane l+1 =
# coef1 (q57601 layout, no shuf1).  For pair i with (a0,a1),(b0,b1) over
# (x^2 - gamma_i):
#     c0 = a0*b0 + gamma_i * (a1*b1)
#     c1 = a0*b1 + a1*b0
# Data-data products use one Montgomery multiply directly, so basemul outputs
# plain products with one missing R factor.  montmul(P, gamma_mont)=P*gamma
# preserves that R^-1 factor on the gamma term.
# --------------------------------------------------------------------------
def pair_index_of(coef): return coef // 2

def check_basemul_adjacency():
    for h in range(2):
        for r in range(8):
            lc = FINAL_IDX[h][r]
            for l in range(0, 16, 2):
                assert lc[l] % 2 == 0 and lc[l+1] == lc[l] + 1, \
                    "basemul adjacency broken h=%d r=%d l=%d: %d,%d" % (h,r,l,lc[l],lc[l+1])
check_basemul_adjacency()

def gamma_vec(h, r):
    """16-lane Montgomery gamma vector: both lanes of a pair carry gamma_i."""
    lc = FINAL_IDX[h][r]
    v = [0]*16
    for l in range(0, 16, 2):
        gi = gamma_of(pair_index_of(lc[l])) * MONT % Q
        v[l] = gi; v[l+1] = gi
    return v

def dmul_u(x, y):
    """GENERAL unsigned Montgomery mul of two data vals: (x*y)*R^-1 mod q in [0,q)."""
    x = u16(x); y = u16(y)
    lo = vpmullw(x, y)
    m  = vpmullw(lo, u16(QINV))
    out = vpmulhuw(x, y)
    mq  = vpmulhuw(m, Q)
    res = u16(out - mq)
    if mq > out:                       # bit-exact borrow correction (see montmul_u)
        res = u16(res - RMQ)
    return res

def addmod_u(a, t):
    """(a+t) mod q for a,t in [0,q)."""
    s = a + t
    return s - Q if s >= Q else s

def basemul_model(na, nb):
    """Pure-Python basemul on NTT-domain arrays na,nb (storage order), bit-exact
    to basemul_avx: data-data dmul gives R^-1-deficient products, gamma is a
    Montgomery constant, then the pair terms are combined modulo q."""
    nc = [0]*256
    for h in range(2):
        for r in range(8):
            lc = FINAL_IDX[h][r]
            base = 128*h + 16*r
            for l in range(0, 16, 2):
                a0 = na[base+l] % Q; a1 = na[base+l+1] % Q
                b0 = nb[base+l] % Q; b1 = nb[base+l+1] % Q
                a0b0 = dmul_u(a0, b0)       # a0*b0*R^-1
                a1b1 = dmul_u(a1, b1)       # a1*b1*R^-1
                a0b1 = dmul_u(a0, b1)       # a0*b1*R^-1
                a1b0 = dmul_u(a1, b0)       # a1*b0*R^-1
                gi = gamma_of(pair_index_of(lc[l])) * MONT % Q
                gab = montmul_u(a1b1, gi)   # gamma*(a1*b1*R^-1)
                nc[base+l]   = addmod_u(a0b0, gab)
                nc[base+l+1] = addmod_u(a0b1, a1b0)
    return nc

def schoolbook(a, b):
    c = [0]*256
    for i in range(256):
        for j in range(256):
            p = a[i]*b[j] % Q; k = i+j
            if k < 256: c[k] = (c[k]+p) % Q
            else:       c[k-256] = (c[k-256]-p) % Q
    return c

# --------------------------------------------------------------------------
def verify():
    import random
    random.seed(2026)
    out = []
    ok = True
    for _ in range(400):
        a = [random.randrange(Q) for _ in range(256)]
        # raw=True exercises the REAL asm path: intt receives the [0,2q) lazy
        # forward output, not a pre-reduced [0,q) copy.
        want = [(x * MONT) % Q for x in a]
        if intt(ntt(a, raw=True)) != want:
            ok = False; break
    out.append("roundtrip-to-mont(400): " + ("PASS" if ok else "FAIL"))
    ok2 = True
    for _ in range(120):
        a = [random.randrange(Q) for _ in range(256)]
        b = [random.randrange(Q) for _ in range(256)]
        nc = basemul_model(ntt(a, raw=True), ntt(b, raw=True))
        if intt(nc) != schoolbook(a, b):
            ok2 = False; break
    out.append("basemul/convolution(120): " + ("PASS" if ok2 else "FAIL"))
    fl = sorted(k for k, v in FWD_LEVEL_CSUB.items() if v)
    fz = sorted(k for k, v in FWD_LEVEL_CSUB.items() if not v)
    il = sorted(k for k, v in INV_LEVEL_CSUB.items() if v)
    iz = sorted(k for k, v in INV_LEVEL_CSUB.items() if not v)
    out.append("forward lazy schedule: CSUB levels %s, LAZY levels %s (L0 lazy=%s)"
               % (fl, fz, not FWD_L0_CSUB))
    out.append("inverse lazy schedule: CSUB levels %s, LAZY levels %s (cross-half L0 csub=%s)"
               % (il, iz, INV_L0_CSUB))
    out.append("bound check: fwd-out max %d < 2^16 (%s); inv-out max %d < 2^16 (%s)"
               % (FWD_OUT_BOUND, FWD_OUT_BOUND < R, INV_OUT_BOUND, INV_OUT_BOUND < R))
    return ok and ok2, out

# --------------------------------------------------------------------------
def s16(v):
    v %= 65536
    return v - 65536 if v >= 32768 else v

def emit_consts(path):
    fwd = []
    for prog in (FWD0, FWD1):
        for op in prog:
            if op[0] in ('xb', 'ib'):
                fwd.append(op[3])
    inv = []
    for prog in (FWD0, FWD1):
        for op in reversed(prog):
            if op[0] in ('xb', 'ib'):
                inv.append([zinv_mont(z) for z in op[3]])

    def zl16(z): return (z * QINV) % R
    def pair_table(table):
        o = []
        for zvec in table:
            o.append([zl16(z) for z in zvec]); o.append(list(zvec))
        return o
    def vecs(name, table):
        lines = [f"const int16_t {name}[{len(table)}][16] = {{"]
        for v in table:
            lines.append("  {" + ",".join(str(s16(x)) for x in v) + "},")
        lines.append("};")
        return "\n".join(lines)
    def pair32(name, z):
        zl = ",".join([str(s16(zl16(z)))]*16)
        zh = ",".join([str(s16(z))]*16)
        return ("const int16_t %s[32] __attribute__((aligned(32))) = {\n  %s,  /* zl */\n  %s   /* zh */\n};\n"
                % (name, zl, zh))
    fwd_p, inv_p = pair_table(fwd), pair_table(inv)

    # basemul gamma table: 16 vectors (h0 r0..7, h1 r0..7) as Seiler pairs.
    gam = []
    for h in range(2):
        for r in range(8):
            gam.append(gamma_vec(h, r))
    gam_p = pair_table(gam)

    with open(path, "w") as f:
        f.write("#include <stdint.h>\n")
        f.write("/* Auto-generated by gen_ntt.py -- do not edit. q=%d, INCOMPLETE n=256,\n"
                "   UNSIGNED [0,q) + LAZY reduction.  Twiddles use the Seiler precompute:\n"
                "   each butterfly stores zl=z*qinv mod 2^16 then zh=z, as two 16-lane vecs. */\n\n" % Q)
        # qdata: 3 vectors consumed by the macros: [0]=q [1]=RMQ(2^16-q) [2]=zero
        f.write("const int16_t ntt_qdata[48] __attribute__((aligned(32))) = {\n")
        f.write("  " + ",".join([str(s16(Q))]*16) + ",  /* _16XQ    (off 0)  */\n")
        f.write("  " + ",".join([str(s16(RMQ))]*16) + ",  /* _16XRMQ  (off 32, 2^16-q borrow const) */\n")
        f.write("  " + ",".join(["0"]*16) + "   /* _16XZERO (off 64) */\n};\n\n")
        # qinv vector (for data*data basemul, where the Seiler precompute can't apply)
        f.write("const int16_t ntt_qinv[16] __attribute__((aligned(32))) = {\n  "
                + ",".join([str(s16(QINV))]*16) + "\n};\n\n")
        # R2 is kept for diagnostics/experiments; Kyber-style basemul no longer
        # pushes an operand through R2.
        f.write(pair32("ntt_r2", R2_MOD_Q))
        f.write(pair32("ntt_z0", Z0))
        f.write(pair32("ntt_z0inv", zinv_mont(Z0)))
        f.write(pair32("ntt_scale", SCALE_TOMONT) + "\n")
        f.write("__attribute__((aligned(32)))\n" + vecs("ntt_zetas_fwd", fwd_p) + "\n\n")
        f.write("__attribute__((aligned(32)))\n" + vecs("ntt_zetas_inv", inv_p) + "\n\n")
        f.write("/* basemul gamma table (Seiler pair per data register, h0 r0..7, h1 r0..7).\n"
                "   Both lanes of a pair carry gamma_i (Montgomery); the odd copy is unused. */\n")
        f.write("__attribute__((aligned(32)))\n" + vecs("ntt_basemul_gamma", gam_p) + "\n")
    return len(fwd_p), len(inv_p), len(gam_p)

# --------------------------------------------------------------------------
# Emit ntt.S from the proven op-list (UNSIGNED + LAZY).  Register map:
#   ymm0 = q   ymm1 = zh   ymm2 = RMQ(2^16-q)   ymm3 = 0
#   ymm4..ymm11 = 8 data regs of one half   ymm12,13 = scratch   ymm14 = t
#   ymm15 = zl (precomputed Montgomery-factor twiddle)
# csub uses ymm12/ymm13 as scratch.  A butterfly is emitted lazy or csub per the
# DERIVED schedule (FWD_SCHED / INV_SCHED).
# --------------------------------------------------------------------------
ASM_PREAMBLE = r"""/* ==========================================================================
 * Auto-generated by gen_ntt.py -- DO NOT EDIT (edit gen_ntt.py, run `make gen`).
 *
 * INCOMPLETE n=256 negacyclic NTT / INTT over Z_q[x]/(x^256+1), q = 30977,
 * 16-bit AVX2, UNSIGNED [0,q) + LAZY REDUCTION.  7 butterfly levels (len 128..2),
 * then a 128-pair BASEMUL (mod x^2 - gamma_i).
 *
 * Entry points (System V AMD64 ABI):
 *   ntt_avx     (rdi=poly[256] i16, rsi=qdata, rdx=zetas_fwd, rcx=z0,    r8=scale)
 *   invntt_tomont_avx(rdi=poly[256] i16, rsi=qdata, rdx=zetas_inv, rcx=z0inv, r8=scale)
 *   basemul_avx (rdi=c, rsi=a, rdx=b, rcx=qdata, r8=gammas)
 *
 * Register map (held live across a half-transform):
 *   ymm0 = q (16x)   ymm1 = zh = zeta vec   ymm2 = RMQ = 2^16-q   ymm3 = 0
 *   ymm4..ymm11 = the 8 data registers of one 128-coeff half
 *   ymm12,ymm13 = scratch   ymm14 = butterfly temp t
 *   ymm15 = zl = zeta*q^-1 mod 2^16 (precomputed Montgomery-factor twiddle)
 *
 * DOMAIN (the REWRITE from the old signed+Barrett path): q = 30977 < 2^15 but
 * 2q = 61954 fits uint16 (< 2^16) while NOT int16 (> 32767).  We therefore use
 * the UNSIGNED [0,q) representation (vpmulhuw + borrow correction), which has ~2x
 * the dynamic range of the signed-centered one for this q.  The montmul folds the
 * reduction in (output [0,q) DIRECTLY -- NO separate Barrett), so it is 3 mults
 * (vs the old signed montmul 3 + Barrett red16 3 = 6).
 *
 * LAZY REDUCTION (the point): a value in [0,q) plus t in [0,q) lands in [0,2q),
 * which fits uint16, so a LAZY butterfly emits a PLAIN vpaddw (no reduce).  When
 * the running bound would exceed uint16 (or violate a montmul input cap) the next
 * level does a MULTIPLY-FREE conditional subtract (csub: s-q via vpsubw + a sign-
 * mask add-back).  The exact per-register bound tracker in gen_ntt.py DERIVES the
 * csub schedule for BOTH forward (CT) and inverse (GS) -- they differ because GS
 * reduces one output for free (its b-output passes through montmul).
 *
 * PERFORMANCE (Seiler precompute, eprint 2018/039): per twiddle store
 *   zl = zeta*qinv mod 2^16   and   zh = zeta
 * so montmul is 1 vpmullw + 2 vpmulhuw (4->3 mul, dep chain 3->2).
 * ========================================================================== */

/* montmul: out = (in * zeta) * R^-1 mod q, result in [0,q).
 *   zl in ymm15 (= zeta*qinv mod 2^16),  zh in ymm1 (= zeta).
 *   m    = in*zl mod 2^16 (vpmullw); out = mulhi_u(in,zh) (vpmulhuw, UNSIGNED);
 *   mq   = mulhi_u(m,q) (vpmulhuw); out -= mq; borrow-correct +RMQ if underflow.
 * Clobbers ymm12,ymm13. */
.macro montmul out,in
  vpmullw  %ymm15,\in,%ymm12
  vpmulhuw %ymm1,\in,\out
  vpmulhuw %ymm0,%ymm12,%ymm12
  vpsubusw \out,%ymm12,%ymm13
  vpsubw   %ymm12,\out,\out
  vpcmpeqw %ymm3,%ymm13,%ymm13
  vpandn   %ymm2,%ymm13,%ymm13
  vpsubw   %ymm13,\out,\out
.endm

/* csub: conditional subtract.  in s in [0,2q) -> [0,q).  t=s-q; mask=t>>15
 * (all-ones where t<0, i.e. s<q); s = t + (mask & q).  Clobbers ymm12,ymm13. */
.macro csub a
  vpsubw   %ymm0,\a,\a
  vpsraw   $15,\a,%ymm13
  vpand    %ymm0,%ymm13,%ymm13
  vpaddw   %ymm13,\a,\a
.endm

/* ctbf_lazy: forward CT, LAZY (inputs [0,q), outputs [0,2q), NO reduce).
 *   t = montmul(b,zeta);  b = a + (q - t);  a = a + t.   (q-t keeps b >= 0) */
.macro ctbf_lazy a,b
  montmul %ymm14,\b
  vpsubw  %ymm14,%ymm0,%ymm12
  vpaddw  %ymm12,\a,\b
  vpaddw  %ymm14,\a,\a
.endm

/* ctbf_csub: forward CT, CSUB (inputs [0,2q); csub a, then csub both outputs ->
 * [0,q)).  t = montmul(b,zeta); ar = csub(a); b = csub(ar+(q-t)); a = csub(ar+t). */
.macro ctbf_csub a,b
  montmul %ymm14,\b
  csub    \a
  vpsubw  %ymm14,%ymm0,%ymm12
  vpaddw  %ymm12,\a,\b
  csub    \b
  vpaddw  %ymm14,\a,\a
  csub    \a
.endm

/* gsbf_lazy: inverse GS, LAZY.  s = u+v (lazy [0,2q)); d = u+(q-v); b = montmul(d).
 *   the b-output is reduced [0,q) for free by montmul; only s stays lazy.
 *   d is built in ymm14 (montmul clobbers ymm12/ymm13, so the input must NOT live
 *   in either -- using ymm12 would self-clobber inside montmul). */
.macro gsbf_lazy a,b
  vpsubw  \b,%ymm0,%ymm14
  vpaddw  \a,%ymm14,%ymm14
  vpaddw  \b,\a,\a
  montmul \b,%ymm14
.endm

/* gsbf_csub: inverse GS, CSUB.  csub u,v; s = csub(u+v); d = u+(q-v); b = montmul(d).
 *   d in ymm14 (see gsbf_lazy note). */
.macro gsbf_csub a,b
  csub    \a
  csub    \b
  vpsubw  \b,%ymm0,%ymm14
  vpaddw  \a,%ymm14,%ymm14
  vpaddw  \b,\a,\a
  csub    \a
  montmul \b,%ymm14
.endm

/* dmul: GENERAL Montgomery multiply of two DATA vectors (no Seiler precompute,
 * both operands vary).  out = (x*y)*R^-1 mod q in [0,q).
 *   m = (x*y mod 2^16)*qinv mod 2^16; out = mulhi_u(x,y) - mulhi_u(m,q) + borrow.
 * qinv vector in ymm9 (basemul reserves it).  Clobbers ymm12,ymm13. */
.macro dmul out,x,y
  vpmullw  \x,\y,%ymm12
  vpmullw  %ymm9,%ymm12,%ymm12
  vpmulhuw \x,\y,\out
  vpmulhuw %ymm0,%ymm12,%ymm12
  vpsubusw \out,%ymm12,%ymm13
  vpsubw   %ymm12,\out,\out
  vpcmpeqw %ymm3,%ymm13,%ymm13
  vpandn   %ymm2,%ymm13,%ymm13
  vpsubw   %ymm13,\out,\out
.endm

/* ---- 3-tier shuffle network (incomplete NTT stops at len=2, no shuf1) ----
 * shuf8 (granularity 8, the only cross-128-bit tier) / shuf4 (64-bit) /
 * shuf2 (32-bit dword).  All involutions, so the inverse replays them. */
.macro shuf8 r0,r1
  vperm2i128 $0x20,\r1,\r0,%ymm12
  vperm2i128 $0x31,\r1,\r0,\r1
  vmovdqa  %ymm12,\r0
.endm
.macro shuf4 r0,r1
  vpunpcklqdq \r1,\r0,%ymm12
  vpunpckhqdq \r1,\r0,\r1
  vmovdqa  %ymm12,\r0
.endm
.macro shuf2 r0,r1
  vmovsldup \r1,%ymm12
  vpblendd $0xAA,%ymm12,\r0,%ymm12
  vpsrlq   $32,\r0,\r0
  vpblendd $0xAA,\r1,\r0,\r1
  vmovdqa  %ymm12,\r0
.endm
"""

def emit_input_conversion_asm(L):
    rec = (1 << 47) // Q
    rec_lo = rec & 0xffffffff
    rec_hi = rec >> 32
    loops = globals().get("N", 256) // 8
    L += [".text",
          "",
          "/* ============ NTT input conversion helpers (standard order) ============",
          " * C ABI:",
          " *   poly_freeze_to_u16_avx(dst:int16_t*, src:int32_t*)       // freeze(src) in [0,q)",
          " *   poly_freeze_small_to_u16_avx(dst:int16_t*, src:int32_t*) // src in [-(q-1),q-1]",
          " *   poly_cast_modq_avx(dst:int16_t*, src:int32_t*)           // src already [0,q)",
          " * These helpers feed ntt_avx's standard-order input buffer; unlike",
          " * nttunpack_avx, they do not apply the AVX NTT slot shuffle ladder. */",
          ".p2align 5",
          ".global poly_freeze_to_u16_avx",
          "poly_freeze_to_u16_avx:",
          "  vpbroadcastd .Lfreeze_q_avx(%rip),%ymm7",
          "  vpbroadcastd .Lfreeze_rec_lo_avx(%rip),%ymm5",
          "  vpbroadcastd .Lfreeze_rec_hi_avx(%rip),%ymm8",
          f"  mov ${loops},%ecx",
          ".Lfreeze_full_loop_avx:",
          "  vmovdqu  (%rsi),%ymm0",
          "  vpsrad   $31,%ymm0,%ymm1",
          "  vpxor    %ymm1,%ymm0,%ymm2",
          "  vpsubd   %ymm1,%ymm2,%ymm2",
          "  vpmuludq %ymm5,%ymm2,%ymm3",
          "  vpsrlq   $47,%ymm3,%ymm3",
          "  vpsrlq   $32,%ymm2,%ymm4",
          "  vpmuludq %ymm5,%ymm4,%ymm4",
          "  vpsrlq   $47,%ymm4,%ymm4",
          "  vpsllq   $32,%ymm4,%ymm4",
          "  vpor     %ymm4,%ymm3,%ymm3",
          "  vpsrld   $15,%ymm2,%ymm4",
          "  vpmulld  %ymm8,%ymm4,%ymm4",
          "  vpaddd   %ymm4,%ymm3,%ymm3",
          "  vpmulld  %ymm7,%ymm3,%ymm3",
          "  vpsubd   %ymm3,%ymm2,%ymm2",
          "  vpsubd   %ymm7,%ymm2,%ymm2",
          "  vpsrad   $31,%ymm2,%ymm4",
          "  vpand    %ymm7,%ymm4,%ymm4",
          "  vpaddd   %ymm4,%ymm2,%ymm2",
          "  vpxor    %ymm1,%ymm2,%ymm2",
          "  vpsubd   %ymm1,%ymm2,%ymm2",
          "  vpsrad   $31,%ymm2,%ymm4",
          "  vpand    %ymm7,%ymm4,%ymm4",
          "  vpaddd   %ymm4,%ymm2,%ymm2",
          "  vextracti128 $1,%ymm2,%xmm3",
          "  vpackusdw %xmm3,%xmm2,%xmm2",
          "  vmovdqu  %xmm2,(%rdi)",
          "  add      $32,%rsi",
          "  add      $16,%rdi",
          "  dec      %ecx",
          "  jnz      .Lfreeze_full_loop_avx",
          "  ret",
          "",
          ".p2align 5",
          ".global poly_freeze_small_to_u16_avx",
          "poly_freeze_small_to_u16_avx:",
          "  vpbroadcastd .Lfreeze_q_avx(%rip),%ymm7",
          f"  mov ${loops},%ecx",
          ".Lfreeze_small_loop_avx:",
          "  vmovdqu  (%rsi),%ymm0",
          "  vpsrad   $31,%ymm0,%ymm1",
          "  vpand    %ymm7,%ymm1,%ymm1",
          "  vpaddd   %ymm1,%ymm0,%ymm0",
          "  vextracti128 $1,%ymm0,%xmm3",
          "  vpackusdw %xmm3,%xmm0,%xmm0",
          "  vmovdqu  %xmm0,(%rdi)",
          "  add      $32,%rsi",
          "  add      $16,%rdi",
          "  dec      %ecx",
          "  jnz      .Lfreeze_small_loop_avx",
          "  ret",
          "",
          ".p2align 5",
          ".global poly_cast_modq_avx",
          "poly_cast_modq_avx:",
          f"  mov ${loops},%ecx",
          ".Lfreeze_cast_loop_avx:",
          "  vmovdqu  (%rsi),%ymm0",
          "  vextracti128 $1,%ymm0,%xmm3",
          "  vpackusdw %xmm3,%xmm0,%xmm0",
          "  vmovdqu  %xmm0,(%rdi)",
          "  add      $32,%rsi",
          "  add      $16,%rdi",
          "  dec      %ecx",
          "  jnz      .Lfreeze_cast_loop_avx",
          "  ret",
          "",
          ".section .rodata",
          ".p2align 4",
          f".Lfreeze_q_avx:       .long {Q}",
          f".Lfreeze_rec_lo_avx:  .long {rec_lo}",
          f".Lfreeze_rec_hi_avx:  .long {rec_hi}",
          ".text"]

def emit_asm(path):
    ymm = lambda r: "%%ymm%d" % (4 + r)
    def load_consts():
        return ["  vmovdqa   0(%rsi),%ymm0",      # q
                "  vmovdqa  32(%rsi),%ymm2",      # RMQ
                "  vmovdqa  64(%rsi),%ymm3"]      # zero
    def load_tw(reg, disp=0):
        return ["  vmovdqa  %d(%s),%%ymm15" % (disp, reg),
                "  vmovdqa  %d(%s),%%ymm1"  % (disp + 32, reg)]
    def load_half(h):
        return ["  vmovdqu  %d(%%rdi),%s" % (256*h + 32*r, ymm(r)) for r in range(8)]
    def store_half(h):
        return ["  vmovdqu  %s,%d(%%rdi)" % (ymm(r), 256*h + 32*r) for r in range(8)]
    XB_LENS = [64]*4 + [32]*4 + [16]*4

    # The forward/inverse butterfly schedules (per butterfly, in program order).
    fsched = iter(FWD_SCHED)
    # ---- forward ----
    L = list(ASM_PREAMBLE.splitlines())
    L += [".text", ".global ntt_avx", ".global invntt_tomont_avx", ".global basemul_avx", ".global reduce_avx",
          "", "/* ============ forward NTT (7 levels, len 128..2) ============ */",
          "ntt_avx:"]
    L += load_consts()
    L += ["", "/* --- level 0 (len=128): cross-HALF CT, partners 128 apart.",
          "       inputs [0,q) -> LAZY add, outputs [0,2q) (no reduce). --- */"]
    L += load_tw("%rcx")
    bf0 = "ctbf_csub" if FWD_L0_CSUB else "ctbf_lazy"
    for r in range(8):
        L += ["  vmovdqu  %d(%%rdi),%%ymm4" % (32*r),
              "  vmovdqu  %d(%%rdi),%%ymm5" % (256 + 32*r),
              "  %s %%ymm4,%%ymm5" % bf0,
              "  vmovdqu  %%ymm4,%d(%%rdi)" % (32*r),
              "  vmovdqu  %%ymm5,%d(%%rdi)" % (256 + 32*r)]
    ki = 0
    # forward butterfly schedule repeats identically for both halves; rebuild it.
    half_sched = list(FWD_SCHED)
    for h, prog in ((0, FWD0), (1, FWD1)):
        lo, hi = 128*h, 128*h + 127
        L += ["", "/* ===== half %d (coeffs %d..%d): levels 1-6 merged ===== */" % (h, lo, hi)]
        L += load_half(h)
        xi = 0; cur = None; si = 0
        for op in prog:
            if op[0] == 'xb':
                ln = XB_LENS[xi]; xi += 1; lvl, do_csub = half_sched[si]; si += 1
                bf = "ctbf_csub" if do_csub else "ctbf_lazy"
                tag = "" if do_csub else " [LAZY]"
                if ln != cur:
                    cur = ln
                    L += ["  /* level %d (len=%d): cross-REGISTER CT%s */" % (lvl, ln, tag)]
                _, r, p, _, _ = op
                L += load_tw("%rdx", ki*64) + ["  %s %s,%s" % (bf, ymm(r), ymm(p))]
                ki += 1
            elif op[0] == 'sh':
                _, k, r0, r1, lvl = op
                if k != cur:
                    cur = k
                    L += ["  /* level %d (len=%d): IN-REGISTER shuf%d then CT */" % (lvl, k, k)]
                L += ["  shuf%d %s,%s" % (k, ymm(r0), ymm(r1))]
            else:  # ib
                _, r0, r1, _, lvl = op
                _, do_csub = half_sched[si]; si += 1
                bf = "ctbf_csub" if do_csub else "ctbf_lazy"
                tag = "" if do_csub else "  /* [LAZY] */"
                L += load_tw("%rdx", ki*64) + ["  %s %s,%s%s" % (bf, ymm(r0), ymm(r1), tag)]
                ki += 1
        L += store_half(h)
    L += ["  ret", ""]

    # ---- inverse ----  (rdx=zinv table, rcx=z0inv, r8=scale).
    L += ["/* ============ inverse NTT (GS, levels 6..0) ============ */",
          "invntt_tomont_avx:"]
    L += load_consts()
    ki = 0
    inv_half_sched = list(INV_SCHED)
    for h, prog in ((0, FWD0), (1, FWD1)):
        lo, hi = 128*h, 128*h + 127
        L += ["", "/* ===== half %d (coeffs %d..%d): levels 6-1 reversed ===== */" % (h, lo, hi)]
        L += load_half(h)
        rev = list(reversed(prog))
        fl = []; xi = 0; last_sh = None
        for op in prog:
            if op[0] == 'xb': fl.append(XB_LENS[xi]); xi += 1
            elif op[0] == 'sh': last_sh = op[1]; fl.append(op[1])
            else: fl.append(last_sh)
        rl = list(reversed(fl)); cur = None; si = 0
        for op, ln in zip(rev, rl):
            if op[0] in ('xb', 'ib'):
                lvl, do_csub = inv_half_sched[si]; si += 1
                bf = "gsbf_csub" if do_csub else "gsbf_lazy"
                tag = "" if do_csub else " [LAZY]"
                kind = "cross-REGISTER GS" if op[0] == 'xb' else "IN-REGISTER GS (after shuf%d)" % ln
                if ln != cur:
                    cur = ln
                    L += ["  /* level %d (len=%d): %s%s */" % (lvl, ln, kind, tag)]
                _, r0, r1, _, _ = op
                L += load_tw("%rdx", ki*64) + ["  %s %s,%s" % (bf, ymm(r0), ymm(r1))]
                ki += 1
            else:
                _, k, r0, r1, _ = op
                L += ["  shuf%d %s,%s" % (k, ymm(r0), ymm(r1))]
        L += store_half(h)
    L += ["", "/* --- level 0 (len=128): cross-HALF GS (z0inv) --- */"]
    L += load_tw("%rcx")
    bfl0 = "gsbf_csub" if INV_L0_CSUB else "gsbf_lazy"
    for r in range(8):
        L += ["  vmovdqu  %d(%%rdi),%%ymm4" % (32*r),
              "  vmovdqu  %d(%%rdi),%%ymm5" % (256 + 32*r),
              "  %s %%ymm4,%%ymm5" % bfl0,
              "  vmovdqu  %%ymm4,%d(%%rdi)" % (32*r),
              "  vmovdqu  %%ymm5,%d(%%rdi)" % (256 + 32*r)]
    L += ["", "/* --- final scaling: multiply every coeff by 1/2^7 and R.",
          "       csub first (lazy GS s-outputs may be [0,2q)); ntt_scale stores",
          "       128^-1*R^2, so this is invntt_tomont. --- */"]
    L += load_tw("%r8")
    for r in range(16):
        L += ["  vmovdqu  %d(%%rdi),%%ymm4" % (32*r),
              "  csub %ymm4",
              "  montmul %ymm5,%ymm4",
              "  vmovdqu  %%ymm5,%d(%%rdi)" % (32*r)]
    L += ["  ret", ""]

    # ---- basemul ----  c = a o b mod (x^2 - gamma_i), q57601 adjacent-lane layout.
    # ABI: rdi=c, rsi=a, rdx=b, rcx=qdata, r8=gammas (one Seiler pair per reg).
    L += ["/* ============ basemul: c = a o b mod (x^2 - gamma_i) ============",
          "   ABI: rdi=c rsi=a rdx=b rcx=qdata r8=gamma (Seiler pairs, 1/reg).",
          "   even lane=coef0, odd lane=coef1.  csub lazy inputs, then direct dmul(a,b)=a*b*R^-1.",
          "   P=a0b0|a1b1 ; X=a0b1|a1b0 ; G=montmul(P,gamma).",
          "   c0=a0b0+gamma*a1b1 (even) ; c1=a0b1+a1b0 (odd). swap16=vpshuflw/hw $0xB1. */",
          "basemul_avx:",
          "  vmovdqa   0(%rcx),%ymm0",          # q
          "  vmovdqa  32(%rcx),%ymm2",          # RMQ
          "  vmovdqa  64(%rcx),%ymm3",          # zero
          "  vmovdqa  ntt_qinv(%rip),%ymm9"]    # qinv for dmul, kept live
    for ri in range(16):
        off  = 32*ri
        goff = ri*64
        L += [
          "  /* --- register %d (pairs %d..%d) --- */" % (ri, ri*8, ri*8+7),
          "  vmovdqu  %d(%%rsi),%%ymm4" % off,                 # a
          "  vmovdqu  %d(%%rdx),%%ymm5" % off,                 # b
          "  csub     %ymm4",
          "  csub     %ymm5",
          "  vpshuflw $0xB1,%ymm5,%ymm7",
          "  vpshufhw $0xB1,%ymm7,%ymm7",                      # ymm7 = swap16(b)
          "  dmul     %ymm8,%ymm4,%ymm5",                      # ymm8  = P = a0b0|a1b1, R^-1
          "  dmul     %ymm10,%ymm4,%ymm7",                     # ymm10 = X = a0b1|a1b0
          "  vmovdqa  %d(%%r8),%%ymm15" % goff,
          "  vmovdqa  %d(%%r8),%%ymm1"  % (goff+32),
          "  montmul  %ymm11,%ymm8",                           # ymm11 = G = a0b0*g | a1b1*g
          "  vpshuflw $0xB1,%ymm11,%ymm13",
          "  vpshufhw $0xB1,%ymm13,%ymm13",                    # swap16(G): even=a1b1*g
          # q30977 has 2q<2^16, so these basemul sums fit in uint16.  Use
          # lazy add followed by the cheaper csub instead of full addmod.
          "  vpaddw   %ymm13,%ymm8,%ymm11",                    # ymm11 = c0 (even valid, <2q)
          "  csub     %ymm11",
          # c1 = X_even + (swap X)_even = a0b1 + a1b0, then csub.
          "  vpshuflw $0xB1,%ymm10,%ymm13",
          "  vpshufhw $0xB1,%ymm13,%ymm13",                    # swap16(X)
          "  vpaddw   %ymm13,%ymm10,%ymm10",                   # ymm10 = c1 (odd valid, <2q)
          "  csub     %ymm10",
          "  vpblendw $0xAA,%ymm10,%ymm11,%ymm11",             # even<-c0, odd<-c1
          "  vmovdqu  %%ymm11,%d(%%rdi)" % off,
        ]
    L += ["  ret", ""]

    L += ["/* ============ reduce_avx: canonicalize 256 uint16 coeffs mod q ============ */",
          "/* q=30977 has 65535 < 3q, so two unsigned conditional subtracts are enough. */",
          "reduce_avx:",
          "  vmovdqa  ntt_qdata(%rip),%ymm0"]
    for ri in range(16):
        off = 32 * ri
        L += ["  vmovdqu  %d(%%rdi),%%ymm4" % off,
              "  vpsubw   %ymm0,%ymm4,%ymm5",
              "  vpminuw  %ymm5,%ymm4,%ymm4",
              "  vpsubw   %ymm0,%ymm4,%ymm5",
              "  vpminuw  %ymm5,%ymm4,%ymm4",
              "  vmovdqu  %%ymm4,%d(%%rdi)" % off]
    L += ["  ret", ""]

    # ---- nttunpack: standard-order int32 samples -> avx NTT layout ----
    # Kyber/Dilithium-style (poly_nttunpack / nttunpack_avx): replay the forward
    # shuffle network with NO butterflies, so out[s] = src[FINAL_IDX[s]].  The A0
    # seed-to-NTT samples arrive int32 in [0,q) and standard order; this narrows
    # them to uint16 and lays them out in the avx NTT slot order.
    L += ["/* ============ nttunpack: standard int32 samples -> avx NTT layout ============",
          " * Kyber-style: forward shuffle ladder, no butterflies, so",
          " * out[s] = src[FINAL_IDX[s]].  Narrows int32 (in [0,q)) -> uint16 on load.",
          " *   rdi = dst (int16[256]), rsi = src (int32[256]). */",
          ".global nttunpack_avx", "nttunpack_avx:"]
    for h, prog in ((0, FWD0), (1, FWD1)):
        for r in range(8):
            so = (128 * h + 16 * r) * 4
            L += ["  vmovdqu  %d(%%rsi),%%ymm12" % so,
                  "  vmovdqu  %d(%%rsi),%%ymm13" % (so + 32),
                  "  vpackusdw %%ymm13,%%ymm12,%s" % ymm(r),
                  "  vpermq   $0xD8,%s,%s" % (ymm(r), ymm(r))]
        for op in prog:
            if op[0] == 'sh':
                k, r0, r1 = op[1], op[2], op[3]
                L += ["  shuf%d %s,%s" % (k, ymm(r0), ymm(r1))]
        for r in range(8):
            L += ["  vmovdqu  %s,%d(%%rdi)" % (ymm(r), (128 * h + 16 * r) * 2)]
    L += ["  ret", ""]

    # ---- nttunpack16: standard-order uint16 samples -> avx NTT layout ----
    # Same shuffle ladder as nttunpack_avx, but the source is ALREADY int16 (in
    # [0,q)), so load it directly -- no int32 widen / vpackusdw / vpermq.  Lets
    # poly_ntt_import unpack a poly16 in place without an int32 staging copy.
    # In-place safe: within each 128-coeff half all 8 regs are loaded before any
    # store, and the two halves touch disjoint [0,256)/[256,512) byte ranges.
    L += ["/* ====== nttunpack16: standard uint16 samples -> avx NTT layout ======",
          " * Same shuffle ladder as nttunpack_avx; src is ALREADY int16 in [0,q),",
          " * so load directly (no int32 widen / vpackusdw / vpermq).  In-place safe.",
          " *   rdi = dst (int16[256]), rsi = src (int16[256]). */",
          ".global nttunpack16_avx", "nttunpack16_avx:"]
    for h, prog in ((0, FWD0), (1, FWD1)):
        for r in range(8):
            so = (128 * h + 16 * r) * 2
            L += ["  vmovdqu  %d(%%rsi),%s" % (so, ymm(r))]
        for op in prog:
            if op[0] == 'sh':
                k, r0, r1 = op[1], op[2], op[3]
                L += ["  shuf%d %s,%s" % (k, ymm(r0), ymm(r1))]
        for r in range(8):
            L += ["  vmovdqu  %s,%d(%%rdi)" % (ymm(r), (128 * h + 16 * r) * 2)]
    L += ["  ret", ""]

    emit_input_conversion_asm(L)

    L += ['.section .note.GNU-stack,"",@progbits']

    body = "\n".join(L) + "\n"
    with open(path, "w") as f:
        f.write(body)
    return len(body.splitlines())

if __name__ == "__main__":
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from avx512_codegen import Avx512Config, emit_avx512
    from namespace_generated import namespace_avx2, namespace_avx512
    ok, lines = verify()
    print("q=%d psi=%d (psi^128==-1:%s) MONT=%d QINV=%d RMQ=%d SCALE_TOMONT=%d R2=%d SCALE=%d"
          % (Q, PSI, pow(PSI,128,Q)==Q-1, MONT, QINV, RMQ, SCALE_TOMONT, R2_MOD_Q, SCALE))
    for l in lines:
        print(l)
    nf, ni, ng = emit_consts("ntt_consts.c")
    print("emitted ntt_consts.c  fwd_vecs=%d inv_vecs=%d gamma_vecs=%d" % (nf, ni, ng))
    nl = emit_asm("ntt.S")
    namespace_avx2()
    print("emitted ntt.S (%d lines, forward+inverse+basemul)" % nl)
    print("FWD0_ops=%d FWD1_ops=%d" % (len(FWD0), len(FWD1)))
    print("FWD csub butterflies=%d/%d, INV csub butterflies=%d/%d"
          % (sum(1 for _,c in FWD_SCHED if c), len(FWD_SCHED),
             sum(1 for _,c in INV_SCHED if c), len(INV_SCHED)))
    cfg512 = Avx512Config(q=Q, n=256, root=PSI, brv_width=BRV,
                          levels=(128, 64, 32, 16, 8, 4, 2),
                          scale_mont=SCALE_TOMONT, incomplete=True, lazy=True,
                          basemul_abi="qdata_fourth")
    ok512, lines512 = emit_avx512(cfg512)
    namespace_avx512()
    for l in lines512:
        print(l)
    sys.exit(0 if (ok and ok512) else 1)
