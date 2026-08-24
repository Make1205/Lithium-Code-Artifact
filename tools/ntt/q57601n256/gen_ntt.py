#!/usr/bin/env python3
"""
gen_ntt.py -- verified model + .S/constant generator for an INCOMPLETE n=256
negacyclic NTT over Z_q[x]/(x^256+1), q=57601, 16-bit AVX2, Kyber-style level
merging + 3-tier shuffle network, stopping at len=2 (basemul).

WHY INCOMPLETE: q-1 = 2^8*3^2*5^2 = 57600.  n=256 | q-1 but 2n=512 does NOT
divide q-1, so there is NO primitive 512-th root of unity; x^256+1 splits only
into 128 quadratic factors (x^2 - gamma_i), not 256 linear ones.  Hence the
forward transform performs 7 butterfly levels (len 128,64,...,2) and STOPS at
len=2; the pointwise step is a BASEMUL (per-pair product mod x^2 - gamma_i)
instead of a coefficient-wise product.

VALLEY (q>2^15): coefficients live in PLAIN [0,q); since 2q > 2^16 there is ZERO
lazy-reduction headroom -- every add/sub is an overflow-aware modular reduction
and the Montgomery high word must be UNSIGNED (vpmulhuw + borrow correction).

Pipeline (mirrors the q=45569 COMPLETE template):
  1. Build the forward transform as an explicit op-list "program" (cross-register
     vertical butterflies for len>=16, shuffle+butterfly for len<16), baking the
     per-lane Montgomery zeta vectors via coefficient-index tracking.  STOPS at
     len=2 (no len=1 / shuf1 tier).
  2. The inverse is the LITERAL REVERSE of that program (CT->GS, shuffle stays),
     so a forward/inverse roundtrip is correct by construction.
  3. Verify with a pure-Python AVX2 lane-permutation interpreter: roundtrip AND
     full negacyclic convolution a*b via ntt -> basemul -> intt == schoolbook.
  4. Emit ntt.S (forward+inverse+basemul) and ntt_consts.c.
"""
import sys

Q = 57601
N = 256
LEVELS = 7                         # incomplete: stop at len=2 -> 7 CT levels
BRV_W = 7                          # brv width = #levels
R = 1 << 16
MONT = R % Q                       # 7935   = 2^16 mod q
QINV = pow(Q, -1, R)               # 7937   = q^-1 mod 2^16
RMQ  = R - Q                       # 7935   = 2^16 - q  (borrow-correction const)
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
# NTT is a true transpose of the scalar ntt16 -- this is what lets the A0
# layout permutation be a clean shuffle ladder (nttunpack) instead of a scramble.
PSI = _primitive_root(Q, 256)      # 68: smallest primitive 256th root, psi^128 == -1
RINV = pow(R, -1, Q)
# scale = (2^LEVELS)^{-1} = 128^{-1} mod q  (undo the GS doublings across 7 levels)
SCALE = pow(1 << LEVELS, -1, Q)
NINV  = pow(N, -1, Q)              # 256^-1 mod q (informational)

def brv(x, w):
    r = 0
    for i in range(w):
        r |= ((x >> i) & 1) << (w - 1 - i)
    return r

# zeta[k] = psi^brv(k, BRV_W) * MONT mod q   (standard CT loop order, k=1..127)
ZETA_MONT = [pow(PSI, brv(k, BRV_W), Q) * MONT % Q for k in range(128)]
R2_MOD_Q   = MONT * MONT % Q
SCALE_TOMONT = SCALE * R2_MOD_Q % Q

def fqmul(a, b):
    return (a * b * RINV) % Q

# (len,start) -> global zeta index, standard CT loop order, len 128..2
FWD_K = {}
_k = 1; _len = 128
while _len >= 2:
    _s = 0
    while _s < N:
        FWD_K[(_len, _s)] = _k; _k += 1; _s += 2 * _len
    _len >>= 1

def zeta_of(length, coeff_low):
    start = (coeff_low // (2 * length)) * (2 * length)
    return ZETA_MONT[FWD_K[(length, start)]]

# basemul twiddle: pair i (coeffs 2i,2i+1) lives mod (x^2 - gamma_i),
# gamma_i = psi^(2*brv(i,BRV_W)+1) mod q.
def gamma_of(i):
    return pow(PSI, 2 * brv(i, BRV_W) + 1, Q)

# --------------------------------------------------------------------------
# AVX2 lane permutations, modelled INSTRUCTION-ACCURATELY (16 int16 lanes,
# lanes 0-7 = low 128-bit, lanes 8-15 = high 128-bit).
# --------------------------------------------------------------------------
def _dwords(x):  return [x[2*i:2*i+2] for i in range(8)]
def _undw(dws):  return [v for dw in dws for v in dw]

def vperm2i128(a, b, imm):
    lo = {0: a[0:8], 1: a[8:16], 2: b[0:8], 3: b[8:16]}
    return lo[imm & 3] + lo[(imm >> 4) & 3]
def vpunpcklqdq(a, b):
    return a[0:4] + b[0:4] + a[8:12] + b[8:12]
def vpunpckhqdq(a, b):
    return a[4:8] + b[4:8] + a[12:16] + b[12:16]
def vmovsldup(x):
    d = _dwords(x); return _undw([d[0],d[0],d[2],d[2],d[4],d[4],d[6],d[6]])
def vpsrlq(x, n):
    assert n == 32
    d = _dwords(x); z = [0,0]
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
# Build the forward op-list for one half (base = 0 or 128).
#   ('xb', ra, rb, zvec)    cross-register vertical CT butterfly  (len 64,32,16)
#   ('sh', k, r2m, r2m1)    in-place shuffle of a register pair   (len 8,4,2)
#   ('ib', r2m, r2m1, zvec) vertical CT butterfly after a shuffle (len 8,4,2)
# Stops at len=2 (shuf2 then butterfly), leaving each register holding pairs.
# --------------------------------------------------------------------------
def build_half(base):
    idx = [[base + 16*r + l for l in range(16)] for r in range(8)]
    prog = []
    for length in (64, 32, 16):                 # cross-register, no shuffle
        step = length // 16
        for r in range(8):
            if (r // step) % 2 == 0:
                p = r + step
                z = zeta_of(length, idx[r][0])
                prog.append(('xb', r, p, [z]*16))
    for length in (8, 4, 2):                     # in-register, shuffle first
        sh = SHUF[length]
        for m in range(4):
            a, b = idx[2*m], idx[2*m+1]
            c, d = sh(a, b)
            zvec = [zeta_of(length, c[l]) for l in range(16)]
            prog.append(('sh', length, 2*m, 2*m+1))
            prog.append(('ib', 2*m, 2*m+1, zvec))
            idx[2*m], idx[2*m+1] = c, d           # chained: stay shuffled
    return prog, idx

FWD0, IDX0 = build_half(0)
FWD1, IDX1 = build_half(128)
Z0   = zeta_of(128, 0)                            # level-0 (cross-half) zeta
FINAL_IDX = [IDX0, IDX1]                          # final storage layout (coeff per lane)

# --------------------------------------------------------------------------
# Interpreter: run a forward program on 8 value-registers of a half.
# --------------------------------------------------------------------------
def run_fwd_half(regs, prog):
    for op in prog:
        if op[0] in ('xb', 'ib'):
            _, ra, rb, z = op
            a, b = regs[ra], regs[rb]
            for l in range(16):
                t = fqmul(z[l], b[l])
                a[l], b[l] = (a[l]+t) % Q, (a[l]-t) % Q
        else:  # 'sh'
            _, k, r0, r1 = op
            regs[r0], regs[r1] = SHUF[k](regs[r0], regs[r1])
    return regs

def zinv_mont(zm):
    plain = zm * RINV % Q
    return pow(plain, -1, Q) * MONT % Q

def run_inv_half(regs, prog):
    for op in reversed(prog):
        if op[0] in ('xb', 'ib'):
            _, r0, r1, z = op
            a, b = regs[r0], regs[r1]
            for l in range(16):
                zi = zinv_mont(z[l])
                u, v = a[l], b[l]
                a[l] = (u + v) % Q
                b[l] = fqmul(zi, (u - v) % Q)
        else:  # 'sh'
            _, k, r0, r1 = op
            regs[r0], regs[r1] = SHUF[k](regs[r0], regs[r1])
    return regs

def to_regs(poly, base):
    return [[poly[base+16*r+l] for l in range(16)] for r in range(8)]
def from_regs(regs, poly, base):
    for r in range(8):
        for l in range(16):
            poly[base+16*r+l] = regs[r][l]

def ntt(poly):
    a = poly[:]
    for j in range(128):                         # level 0 (cross-half)
        t = fqmul(Z0, a[j+128])
        a[j+128] = (a[j]-t) % Q; a[j] = (a[j]+t) % Q
    for h, prog in ((0, FWD0), (1, FWD1)):
        regs = run_fwd_half(to_regs(a, 128*h), prog)
        from_regs(regs, a, 128*h)
    return a

def intt(poly):
    a = poly[:]
    for h, prog in ((0, FWD0), (1, FWD1)):
        regs = run_inv_half(to_regs(a, 128*h), prog)
        from_regs(regs, a, 128*h)
    zi0 = zinv_mont(Z0)
    for j in range(128):                          # reversed CT level 0
        A, B = a[j], a[j+128]
        a[j]     = (A + B) % Q
        a[j+128] = fqmul(zi0, (A - B) % Q)
    return [fqmul(x, SCALE_TOMONT) for x in a]    # final 1/2^7 and to-Montgomery

# --------------------------------------------------------------------------
# BASEMUL.  After 7 levels the storage layout (FINAL_IDX) places, in each data
# register, 8 quadratic pairs: even lane l carries coef0 of a pair, the odd lane
# l+1 carries coef1.  For pair i with (a0,a1),(b0,b1) over (x^2 - gamma_i):
#     c0 = a0*b0 + gamma_i * (a1*b1)
#     c1 = a0*b1 + a1*b0
# Data-data products use one Montgomery multiply directly, so basemul outputs
# plain products with one missing R factor.  gamma is stored *R, so
# montmul(t11,gamma)=t11*gamma and preserves that R^-1 factor.
# --------------------------------------------------------------------------
def pair_index_of(coef):
    return coef // 2

def check_basemul_adjacency():
    for h in range(2):
        for r in range(8):
            lc = FINAL_IDX[h][r]
            for l in range(0, 16, 2):
                assert lc[l] % 2 == 0 and lc[l+1] == lc[l] + 1, \
                    "basemul adjacency broken h=%d r=%d l=%d: %d,%d" % (h,r,l,lc[l],lc[l+1])
check_basemul_adjacency()

def gamma_vec(h, r):
    """16-lane Montgomery gamma vector: both lanes of a pair carry gamma_i (the
    odd lane copy is harmless -- c1 never multiplies by gamma)."""
    lc = FINAL_IDX[h][r]
    v = [0]*16
    for l in range(0, 16, 2):
        gi = gamma_of(pair_index_of(lc[l])) * MONT % Q
        v[l] = gi; v[l+1] = gi
    return v

def basemul_model(na, nb):
    """Pure-Python basemul on NTT-domain arrays na,nb (storage order)."""
    nc = [0]*N
    for h in range(2):
        for r in range(8):
            lc = FINAL_IDX[h][r]
            base = 128*h + 16*r
            for l in range(0, 16, 2):
                a0 = na[base+l]; a1 = na[base+l+1]
                b0 = nb[base+l]; b1 = nb[base+l+1]
                gi = gamma_of(pair_index_of(lc[l]))
                p0 = fqmul(a0, b0)
                p1 = fqmul(a1, b1)
                x0 = fqmul(a0, b1)
                x1 = fqmul(a1, b0)
                c0 = (p0 + fqmul(p1, gi * MONT % Q)) % Q
                c1 = (x0 + x1) % Q
                nc[base+l] = c0; nc[base+l+1] = c1
    return nc

def schoolbook(a, b):
    c = [0]*N
    for i in range(N):
        for j in range(N):
            p = a[i]*b[j] % Q; k = i+j
            if k < N: c[k] = (c[k]+p) % Q
            else:     c[k-N] = (c[k-N]-p) % Q
    return c

# --------------------------------------------------------------------------
def verify():
    import random
    random.seed(2026)
    out = []
    ok = all(intt(ntt(a := [random.randrange(Q) for _ in range(N)])) == [(x * MONT) % Q for x in a]
             for _ in range(300))
    out.append("roundtrip-to-mont(300): " + ("PASS" if ok else "FAIL"))
    ok2 = True
    for _ in range(80):
        a = [random.randrange(Q) for _ in range(N)]
        b = [random.randrange(Q) for _ in range(N)]
        na, nb = ntt(a), ntt(b)
        nc = basemul_model(na, nb)
        if intt(nc) != schoolbook(a, b):
            ok2 = False; break
    out.append("basemul/convolution(80): " + ("PASS" if ok2 else "FAIL"))
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
                fwd.append(op[-1])
    inv = []
    for prog in (FWD0, FWD1):
        for op in reversed(prog):
            if op[0] in ('xb', 'ib'):
                inv.append([zinv_mont(z) for z in op[-1]])

    def zl16(z): return (z * QINV) % R
    def pair_table(table):
        out = []
        for zvec in table:
            out.append([zl16(z) for z in zvec]); out.append(list(zvec))
        return out
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
        f.write("/* Auto-generated by gen_ntt.py -- do not edit. q=%d (INCOMPLETE NTT).\n"
                "   Twiddles use the Seiler precompute: each butterfly stores zl=z*qinv\n"
                "   mod 2^16 then zh=z, as two consecutive 16-lane vectors. */\n\n" % Q)
        f.write("const int16_t ntt_qdata[48] __attribute__((aligned(32))) = {\n")
        f.write("  " + ",".join([str(Q)]*16) + ",  /* _16XQ    */\n")
        f.write("  " + ",".join([str(RMQ)]*16) + ",  /* _16XRMQ  */\n")
        f.write("  " + ",".join(["0"]*16) + "   /* _16XZERO */\n};\n\n")
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
# Emit ntt.S from the proven op-list.
# C ABI (System V): rdi=poly, rsi=qdata[48], rdx=ztab[][16], rcx=z0 pair, r8=scale.
# --------------------------------------------------------------------------
ASM_PREAMBLE = r"""/* ==========================================================================
 * Auto-generated by gen_ntt.py -- DO NOT EDIT (edit gen_ntt.py, run `make gen`).
 *
 * INCOMPLETE n=256 negacyclic NTT / INTT over Z_q[x]/(x^256+1), q = 57601,
 * 16-bit AVX2, Kyber-style level merging + 3-tier shuffle network, stopping at
 * len=2 (basemul).
 *
 * WHY INCOMPLETE: q-1 = 2^8*3^2*5^2; n=256 | q-1 but 2n=512 does NOT, so there
 * is no primitive 512-th root: x^256+1 splits into 128 quadratics (x^2-gamma_i),
 * not 256 linears.  7 CT levels (len 128..2); the pointwise step is a BASEMUL.
 *
 * Entry points (System V AMD64 ABI):
 *   ntt_avx     (rdi=poly[256] i16, rsi=qdata[48], rdx=zetas_fwd, rcx=z0,    r8=scale)
 *   invntt_tomont_avx(rdi=poly[256] i16, rsi=qdata[48], rdx=zetas_inv, rcx=z0inv, r8=scale)
 *   basemul_avx (rdi=c[256] i16,    rsi=qdata[48], rdx=a[256], rcx=b[256], r8=gamma[16][16])
 *
 * Register map (held live across a whole half-transform):
 *   ymm0 = q (16x)   ymm1 = zh = current zeta vec   ymm2 = 2^16-q   ymm3 = 0
 *   ymm4..ymm11 = the 8 data registers of one 128-coeff half
 *   ymm12,ymm13 = scratch   ymm14 = butterfly temp t
 *   ymm15 = zl = (zeta * q^-1) mod 2^16 = the PRECOMPUTED Montgomery-factor twiddle
 *
 * DOMAIN: coefficients stay in PLAIN [0,q); twiddles are *R mod q (Montgomery),
 * so montmul(zeta,x) = zeta*x mod q keeps x's domain.  Because q > 2^15 there is
 * ZERO lazy-reduction headroom: every add/sub is a full overflow-aware modular
 * reduction (2q > 2^16, so a+b does NOT fit 16 bits) and the Montgomery high
 * word must be UNSIGNED (vpmulhuw + a borrow correction).
 *
 * PERFORMANCE (Seiler precompute, eprint 2018/039): per twiddle store
 *   zl = zeta*qinv mod 2^16   and   zh = zeta
 * so montmul is 1 vpmullw + 2 vpmulhuw (4->3 mul, dep chain 3->2).
 * ========================================================================== */

/* montmul: out = (in * zeta) * R^-1 mod q, result in [0,q).
 *   zl in ymm15 (= zeta*qinv mod 2^16),  zh in ymm1 (= zeta).
 *   m    = in*zl mod 2^16   (vpmullw); out = floor(in*zh/2^16) (vpmulhuw, UNSIGNED);
 *   mqhi = floor(m*q/2^16)  (vpmulhuw); out -= mqhi; borrow-correct +RMQ if underflow.
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

/* addmod: dst = (a + t) mod q via a-(q-t) with borrow-driven +RMQ. Clobbers ymm12,ymm13. */
.macro addmod dst,a,t
  vpsubw   \t,%ymm0,%ymm12
  vpsubusw \a,%ymm12,%ymm13
  vpsubw   %ymm12,\a,%ymm12
  vpcmpeqw %ymm3,%ymm13,%ymm13
  vpandn   %ymm2,%ymm13,%ymm13
  vpsubw   %ymm13,%ymm12,\dst
.endm

/* submod: dst = (a - t) mod q with borrow-driven +RMQ. Clobbers ymm12,ymm13. */
.macro submod dst,a,t
  vpsubusw \a,\t,%ymm13
  vpsubw   \t,\a,%ymm12
  vpcmpeqw %ymm3,%ymm13,%ymm13
  vpandn   %ymm2,%ymm13,%ymm13
  vpsubw   %ymm13,%ymm12,\dst
.endm

/* ctbf: Cooley-Tukey butterfly (forward).  t=b*zeta; b=a-t; a=a+t. */
.macro ctbf a,b
  montmul %ymm14,\b
  submod  \b,\a,%ymm14
  addmod  \a,\a,%ymm14
.endm

/* gsbf: Gentleman-Sande butterfly (inverse).  t=a-b; a=a+b; b=t*zeta. */
.macro gsbf a,b
  submod  %ymm14,\a,\b
  addmod  \a,\a,\b
  montmul \b,%ymm14
.endm

/* dmul: GENERAL Montgomery multiply of two DATA vectors (no Seiler precompute
 * possible, both operands vary).  out = (x*y)*R^-1 mod q in [0,q).
 *   m    = (x*y mod 2^16)*qinv mod 2^16   (vpmullw x,y ; vpmullw *qinv)
 *   out  = mulhi_u(x,y) - mulhi_u(m,q)    (two vpmulhuw) + borrow-correct +RMQ.
 * qinv vector must be in ymm9 (basemul reserves it there).  Clobbers ymm12,ymm13. */
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
    loops = N // 8
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
        return ["  vmovdqa   0(%rsi),%ymm0", "  vmovdqa  32(%rsi),%ymm2",
                "  vmovdqa  64(%rsi),%ymm3"]
    def load_tw(reg, disp=0):
        return ["  vmovdqa  %d(%s),%%ymm15" % (disp, reg),
                "  vmovdqa  %d(%s),%%ymm1"  % (disp + 32, reg)]
    def load_half(h):
        return ["  vmovdqu  %d(%%rdi),%s" % (256*h + 32*r, ymm(r)) for r in range(8)]
    def store_half(h):
        return ["  vmovdqu  %s,%d(%%rdi)" % (ymm(r), 256*h + 32*r) for r in range(8)]

    LVL = {64: 1, 32: 2, 16: 3, 8: 4, 4: 5, 2: 6}
    XB_LENS = [64]*4 + [32]*4 + [16]*4

    L = list(ASM_PREAMBLE.splitlines())
    L += [".text", ".global ntt_avx", ".global invntt_tomont_avx", ".global basemul_avx",
          "", "/* ============ forward NTT (7 levels, len 128..2) ============ */",
          "ntt_avx:"]
    L += load_consts()
    L += ["",
          "/* --- level 0 (len=128): cross-HALF CT butterflies, partners 128 apart. --- */"]
    L += load_tw("%rcx")
    for r in range(8):
        L += ["  vmovdqu  %d(%%rdi),%%ymm4" % (32*r),
              "  vmovdqu  %d(%%rdi),%%ymm5" % (256 + 32*r),
              "  ctbf %ymm4,%ymm5",
              "  vmovdqu  %%ymm4,%d(%%rdi)" % (32*r),
              "  vmovdqu  %%ymm5,%d(%%rdi)" % (256 + 32*r)]
    ki = 0
    for h, prog in ((0, FWD0), (1, FWD1)):
        lo, hi = 128*h, 128*h + 127
        L += ["",
              "/* ===== half %d (coeffs %d..%d): levels 1-6 merged, no stores until done ===== */" % (h, lo, hi)]
        L += load_half(h)
        xi = 0; cur = None
        for op in prog:
            if op[0] == 'xb':
                ln = XB_LENS[xi]; xi += 1
                if ln != cur:
                    cur = ln
                    L += ["  /* level %d (len=%d): cross-REGISTER CT, partners %d lanes apart */"
                          % (LVL[ln], ln, ln)]
                _, r, p, _ = op
                L += load_tw("%rdx", ki*64) + ["  ctbf %s,%s" % (ymm(r), ymm(p))]
                ki += 1
            elif op[0] == 'sh':
                _, k, r0, r1 = op
                if k != cur:
                    cur = k
                    L += ["  /* level %d (len=%d): IN-REGISTER, shuf%d interleave then vertical CT */"
                          % (LVL[k], k, k)]
                L += ["  shuf%d %s,%s" % (k, ymm(r0), ymm(r1))]
            else:  # ib
                _, r0, r1, _ = op
                L += load_tw("%rdx", ki*64) + ["  ctbf %s,%s" % (ymm(r0), ymm(r1))]
                ki += 1
        L += store_half(h)
    L += ["  ret", ""]

    # ---- inverse ----
    L += ["/* ============ inverse NTT (GS butterflies, levels 6..0) ============ */",
          "invntt_tomont_avx:"]
    L += load_consts()
    ki = 0
    for h, prog in ((0, FWD0), (1, FWD1)):
        lo, hi = 128*h, 128*h + 127
        L += ["",
              "/* ===== half %d (coeffs %d..%d): levels 6-1 merged (forward order reversed) ===== */" % (h, lo, hi)]
        L += load_half(h)
        rev = list(reversed(prog))
        fl = []; xi = 0; last_sh = None
        for op in prog:
            if op[0] == 'xb': fl.append(XB_LENS[xi]); xi += 1
            elif op[0] == 'sh': last_sh = op[1]; fl.append(op[1])
            else: fl.append(last_sh)
        rl = list(reversed(fl)); cur = None
        for op, ln in zip(rev, rl):
            if op[0] in ('xb', 'ib'):
                kind = "cross-REGISTER GS" if op[0] == 'xb' else "IN-REGISTER GS (after shuf%d)" % ln
                if ln != cur:
                    cur = ln
                    L += ["  /* level %d (len=%d): %s */" % (LVL[ln], ln, kind)]
                _, r0, r1, _ = op
                L += load_tw("%rdx", ki*64) + ["  gsbf %s,%s" % (ymm(r0), ymm(r1))]
                ki += 1
            else:
                _, k, r0, r1 = op
                L += ["  shuf%d %s,%s" % (k, ymm(r0), ymm(r1))]
        L += store_half(h)
    L += ["",
          "/* --- level 0 (len=128): cross-HALF GS butterflies (z0inv) --- */"]
    L += load_tw("%rcx")
    for r in range(8):
        L += ["  vmovdqu  %d(%%rdi),%%ymm4" % (32*r),
              "  vmovdqu  %d(%%rdi),%%ymm5" % (256 + 32*r),
              "  gsbf %ymm4,%ymm5",
              "  vmovdqu  %%ymm4,%d(%%rdi)" % (32*r),
              "  vmovdqu  %%ymm5,%d(%%rdi)" % (256 + 32*r)]
    L += ["",
          "/* --- final scaling: multiply every coefficient by 1/2^7 and R.",
          "       ntt_scale stores 128^-1*R^2, so this is invntt_tomont. --- */"]
    L += load_tw("%r8")
    for r in range(16):
        L += ["  vmovdqu  %d(%%rdi),%%ymm4" % (32*r),
              "  montmul %ymm5,%ymm4",
              "  vmovdqu  %%ymm5,%d(%%rdi)" % (32*r)]
    L += ["  ret", ""]

    # ---- basemul ----
    # rdi=c, rsi=qdata, rdx=a, rcx=b, r8=gamma table (16 Seiler pairs, [zl,zh] per reg)
    # Per data register (16 lanes = 8 pairs); even lane = coef0, odd lane = coef1.
    #   P   = dmul(a, b)                -> even=a0*b0, odd=a1*b1, with R^-1
    #   bs  = swap16(b) (even<->odd per pair); Xc = dmul(a, bs)
    #                                   -> even=a0*b1, odd=a1*b0, with R^-1
    #   c1  = Xc + swap16(Xc)  on odd? simpler: c1 = a0*b1 + a1*b0 = even(Xc)+odd(Xc).
    #         We build a vector with even=a0*b1, odd=a1*b0 (= Xc) and its swap, add.
    #   gP  = even lanes of P times gamma (montmul P_odd-broadcast?)...
    # We need t11 (=odd of P) multiplied by gamma and added to t00 (=even of P).
    # Approach: gmul = montmul(P, gamma) where gamma vector has gamma_i in BOTH
    # lanes of a pair (so odd lane of gmul = a1b1*gamma; we use that).  Then
    #   c0 = even(P) + odd(gmul)         -> but even/odd live in different lanes.
    # Align by swapping: let Ps=swap16(P) (even=a1b1, odd=a0b0). gmul=montmul on
    # P gives (even=a0b0*gamma, odd=a1b1*gamma). We want a0b0 + gamma*a1b1 in the
    # EVEN lane, and a0b1 + a1b0 in the ODD lane.
    #   c0_even = even(P) + odd(montmul(P,gamma_in_odd))  -- align via swap.
    # Concretely (matches basemul_model):
    #   P     = dmul(a,b)                      P[2i]=a0b0, P[2i+1]=a1b1, R^-1
    #   G     = montmul(P, GAMMA)              G[2i]=a0b0*g, G[2i+1]=a1b1*g  (GAMMA=g in both lanes)
    #   Gs    = swap16(G)                      Gs[2i]=a1b1*g, Gs[2i+1]=a0b0*g
    #   c0    = addmod(P_even-take, Gs_even)   -> need a0b0 (P even) + a1b1*g (Gs even)
    #           => c0_even = addmod(P, Gs) then keep EVEN lane
    #   Xc    = dmul(a, swap16(b))             Xc[2i]=a0b1, Xc[2i+1]=a1b0, R^-1
    #   Xcs   = swap16(Xc)                     Xcs[2i]=a1b0, Xcs[2i+1]=a0b1
    #   c1    = addmod(Xc, Xcs)                c1[2i]=c1[2i+1]=a0b1+a1b0 (both lanes equal)
    #   result: even lane <- c0 (a0b0+g*a1b1), odd lane <- c1 (a0b1+a1b0)
    #   blend even from c0vec, odd from c1vec (vpblendw 0xAA picks odd from 2nd).
    L += ["/* ============ basemul: c = a*b mod (x^2 - gamma_i) per pair ============ */",
          "/* rdi=c rsi=qdata rdx=a rcx=b r8=gamma[16] (Seiler pairs).",
          "   even lane=coef0, odd lane=coef1.  Direct dmul(a,b)=a*b*R^-1.",
          "   P=a0b0|a1b1 ; X=a0b1|a1b0 ; G=montmul(P,gamma).",
          "   c0=a0b0+gamma*a1b1 (even) ; c1=a0b1+a1b0 (odd). swap16 = vpshuflw/hw $0xB1. */",
          "basemul_avx:"]
    L += load_consts()
    L += ["  vmovdqa  ntt_qinv(%rip),%ymm9"]   # qinv vector for dmul, kept live
    for ri in range(16):
        off  = 32*ri
        goff = ri*64
        L += [
          "  /* --- register %d (pairs %d..%d) --- */" % (ri, ri*8, ri*8+7),
          "  vmovdqu  %d(%%rdx),%%ymm4" % off,                 # a
          "  vmovdqu  %d(%%rcx),%%ymm5" % off,                 # b
          "  vpshuflw $0xB1,%ymm5,%ymm7",
          "  vpshufhw $0xB1,%ymm7,%ymm7",                      # ymm7 = swap16(b)
          "  dmul     %ymm8,%ymm4,%ymm5",                      # ymm8 = P  = a0b0|a1b1, R^-1
          "  dmul     %ymm10,%ymm4,%ymm7",                     # ymm10 = X = a0b1|a1b0
          # G = montmul(P, gamma)
          "  vmovdqa  %d(%%r8),%%ymm15" % goff,
          "  vmovdqa  %d(%%r8),%%ymm1"  % (goff+32),
          "  montmul  %ymm11,%ymm8",                           # ymm11 = G = a0b0*g | a1b1*g
          "  vpshuflw $0xB1,%ymm11,%ymm13",
          "  vpshufhw $0xB1,%ymm13,%ymm13",                    # ymm13 = swap16(G): even=a1b1*g
          "  addmod   %ymm11,%ymm8,%ymm13",                    # c0 (even lane = a0b0 + g*a1b1)
          # c1 = X + swap16(X)
          "  vpshuflw $0xB1,%ymm10,%ymm13",
          "  vpshufhw $0xB1,%ymm13,%ymm13",                    # swap16(X): even=a1b0, odd=a0b1
          "  addmod   %ymm10,%ymm10,%ymm13",                   # c1 (both lanes = a0b1 + a1b0)
          # blend: even <- c0 (ymm11), odd <- c1 (ymm10).  vpblendw 0xAA picks odd from 2nd src.
          "  vpblendw $0xAA,%ymm10,%ymm11,%ymm11",
          "  vmovdqu  %%ymm11,%d(%%rdi)" % off,
        ]
    L += ["  ret", ""]

    # ---- nttunpack: standard-order int32 samples -> avx NTT layout ----
    # Kyber/Dilithium-style (poly_nttunpack / nttunpack_avx): replay the forward
    # shuffle network with NO butterflies, so out[s] = src[FINAL_IDX[s]].  The
    # A0 seed-to-NTT samples arrive as int32 in [0,q) and standard order; this
    # narrows them to uint16 and lays them out in the avx NTT slot order.
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
                _, k, r0, r1 = op
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
                _, k, r0, r1 = op
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
    cfg512 = Avx512Config(q=Q, n=N, root=PSI, brv_width=BRV_W,
                          levels=(128, 64, 32, 16, 8, 4, 2),
                          scale_mont=SCALE_TOMONT,
                          incomplete=True, lazy=False,
                          basemul_abi="qdata_second")
    ok512, lines512 = emit_avx512(cfg512)
    namespace_avx512()
    for l in lines512:
        print(l)
    sys.exit(0 if (ok and ok512) else 1)
