#!/usr/bin/env python3
"""
gen_ntt.py -- verified model + .S/constant generator for a COMPLETE n=256
negacyclic NTT over Z_q[x]/(x^256+1), q=45569, 16-bit AVX2, Kyber-style level
merging + the 4-tier shuffle network.

Pipeline:
  1. Build the forward transform as an explicit op-list "program" (cross-register
     vertical butterflies for len>=16, shuffle+butterfly for len<16), baking the
     per-lane Montgomery zeta vectors via coefficient-index tracking.
  2. The inverse is the LITERAL REVERSE of that program (CT->GS, shuffle stays),
     so a forward/inverse roundtrip is correct by construction.
  3. Verify both with a pure-Python AVX2 lane-permutation interpreter: roundtrip
     and full negacyclic convolution vs schoolbook.
  4. Emit ntt.S, invntt.S (generated bodies) and ntt_consts.c ({q, 2^16-q, 0}
     + the Seiler-precomputed (zl,zh) twiddle tables in exact consumption order).
     The assembly only replays a proven sequence, so any transcription error is
     caught by the C harness.
"""
import sys

Q = 45569
R = 1 << 16
MONT = R % Q                       # 19967
QINV = pow(Q, -1, R)               # 19969
RMQ  = R - Q                       # 19967  (2^16 - q, the reduction constant)
NINV = pow(256, Q - 2, Q)
G = 3
OMEGA = pow(G, (Q - 1) // 512, Q)  # primitive 512-th root, omega^256 = -1
RINV = pow(R, -1, Q)

def brv8(x):
    r = 0
    for i in range(8):
        r |= ((x >> i) & 1) << (7 - i)
    return r

ZETA_MONT = [pow(OMEGA, brv8(k), Q) * MONT % Q for k in range(256)]
R2_MOD_Q  = MONT * MONT % Q
NINV_TOMONT = NINV * R2_MOD_Q % Q

def fqmul(a, b):
    return (a * b * RINV) % Q

# (len,start) -> global zeta index, standard CT loop order
FWD_K = {}
_k = 1; _len = 128
while _len >= 1:
    _s = 0
    while _s < 256:
        FWD_K[(_len, _s)] = _k; _k += 1; _s += 2 * _len
    _len >>= 1

def zeta_of(length, coeff_low):
    start = (coeff_low // (2 * length)) * (2 * length)
    return ZETA_MONT[FWD_K[(length, start)]]

# --------------------------------------------------------------------------
# AVX2 lane permutations, modelled INSTRUCTION-ACCURATELY (16 int16 lanes,
# lanes 0-7 = low 128-bit, lanes 8-15 = high 128-bit).  Each shuffleK returns
# (c, d) = exactly what the corresponding asm macro (below) computes, so that
# the generated zeta tables and the emitted assembly are guaranteed consistent.
# --------------------------------------------------------------------------
def _dwords(x):  return [x[2*i:2*i+2] for i in range(8)]      # eight 32-bit dwords
def _undw(dws):  return [v for dw in dws for v in dw]

def vperm2i128(a, b, imm):                                    # AT&T: src2=b? we pass (a,b)
    sel = {0: a[0:4], 1: a[4:8], 2: b[0:4], 3: b[4:8]}        # 128-bit halves (8 int16)
    lo = {0: a[0:8], 1: a[8:16], 2: b[0:8], 3: b[8:16]}
    return lo[imm & 3] + lo[(imm >> 4) & 3]
def vpunpcklqdq(a, b):                                        # dst=punpcklqdq(a,b)
    return a[0:4] + b[0:4] + a[8:12] + b[8:12]
def vpunpckhqdq(a, b):
    return a[4:8] + b[4:8] + a[12:16] + b[12:16]
def vmovsldup(x):                                            # dup even dwords
    d = _dwords(x); return _undw([d[0],d[0],d[2],d[2],d[4],d[4],d[6],d[6]])
def vpsrlq(x, n):                                            # n==32: bring odd dword down
    assert n == 32
    d = _dwords(x); z = [0,0]
    return _undw([d[1],z,d[3],z,d[5],z,d[7],z])
def vpblendd(a, b, mask):                                    # dword blend: bit set -> b
    da, db = _dwords(a), _dwords(b)
    return _undw([db[i] if (mask>>i)&1 else da[i] for i in range(8)])
def vpslld(x, n):                                            # n==16: word i -> hi half of its dword
    assert n == 16
    return _undw([[0, x[2*i]] for i in range(8)])
def vpsrld(x, n):
    assert n == 16
    return _undw([[x[2*i+1], 0] for i in range(8)])
def vpblendw(a, b, mask):                                    # word blend (mask repeats per 128)
    return [b[i] if (mask >> (i % 8)) & 1 else a[i] for i in range(16)]

def shuffle8(a, b):                                          # vperm2i128 $0x20 / $0x31
    return vperm2i128(a, b, 0x20), vperm2i128(a, b, 0x31)
def shuffle4(a, b):                                          # vpunpcklqdq / vpunpckhqdq
    return vpunpcklqdq(a, b), vpunpckhqdq(a, b)
def shuffle2(a, b):                                          # Kyber shuffle2 sequence
    # asm: vmovsldup b->t; vpblendd $0xAA,t,a->c; vpsrlq $32,a->a; vpblendd $0xAA,b,a->d
    c = vpblendd(a, vmovsldup(b), 0xAA)
    d = vpblendd(vpsrlq(a, 32), b, 0xAA)
    return c, d
def shuffle1(a, b):                                          # Kyber shuffle1 sequence
    c = vpblendw(a, vpslld(b, 16), 0xAA)                     # word i: odd -> (b<<16), even -> a
    d = vpblendw(vpsrld(a, 16), b, 0xAA)
    return c, d
SHUF = {8: shuffle8, 4: shuffle4, 2: shuffle2, 1: shuffle1}

# --------------------------------------------------------------------------
# Build the forward op-list for one half (base = 0 or 128), baking zeta vectors
# by tracking coefficient indices through the exact shuffles.
# Ops:  ('xb', ra, rb, zvec)              cross-register vertical CT butterfly
#       ('sh', k, r2m, r2m1)              in-place shuffle of a register pair
#       ('ib', r2m, r2m1, zvec)           vertical CT butterfly after a shuffle
# --------------------------------------------------------------------------
def build_half(base):
    idx = [[base + 16*r + l for l in range(16)] for r in range(8)]
    prog = []
    for length in (64, 32, 16):                 # cross-register
        step = length // 16
        for r in range(8):
            if (r // step) % 2 == 0:
                p = r + step
                z = zeta_of(length, idx[r][0])
                prog.append(('xb', r, p, [z]*16))
                # vertical butterfly does not move indices
    for length in (8, 4, 2, 1):                  # in-register
        sh = SHUF[length]
        for m in range(4):
            a, b = idx[2*m], idx[2*m+1]
            c, d = sh(a, b)
            zvec = [zeta_of(length, c[l]) for l in range(16)]
            prog.append(('sh', length, 2*m, 2*m+1))
            prog.append(('ib', 2*m, 2*m+1, zvec))
            idx[2*m], idx[2*m+1] = c, d           # chained: stay shuffled
    return prog

FWD0 = build_half(0)
FWD1 = build_half(128)
Z0   = zeta_of(128, 0)                            # level-0 (cross-half) zeta

# --------------------------------------------------------------------------
# Interpreter: run a forward program on 8 value-registers of a half.
# --------------------------------------------------------------------------
def run_fwd_half(regs, prog):
    for op in prog:
        if op[0] == 'xb':
            _, ra, rb, z = op
            a, b = regs[ra], regs[rb]
            for l in range(16):
                t = fqmul(z[l], b[l])
                a[l], b[l] = (a[l]+t) % Q, (a[l]-t) % Q
        elif op[0] == 'sh':
            _, k, r0, r1 = op
            regs[r0], regs[r1] = SHUF[k](regs[r0], regs[r1])
        else:  # 'ib'
            _, r0, r1, z = op
            a, b = regs[r0], regs[r1]
            for l in range(16):
                t = fqmul(z[l], b[l])
                a[l], b[l] = (a[l]+t) % Q, (a[l]-t) % Q
    return regs

def zinv_mont(zm):
    """Montgomery zeta z*R -> (z^{-1})*R, the twiddle that inverts CT_z."""
    plain = zm * RINV % Q
    return pow(plain, -1, Q) * MONT % Q

def run_inv_half(regs, prog):
    # literal reverse of forward: each CT_z butterfly is undone by
    #   out_a = a'+b',  out_b = z^{-1} (a'-b')   (2x factors folded into 1/n);
    # shuffles are involutions so they replay unchanged.
    for op in reversed(prog):
        if op[0] in ('xb', 'ib'):
            tag, r0, r1, z = op
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
    for j in range(128):                         # level 0
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
    zi0 = zinv_mont(Z0)                            # level 0 inverse twiddle
    for j in range(128):                          # reversed CT level 0
        A, B = a[j], a[j+128]
        a[j]     = (A + B) % Q
        a[j+128] = fqmul(zi0, (A - B) % Q)
    # final 1/n and to-Montgomery: invntt_tomont(ntt(a)) = a*R.
    return [fqmul(x, NINV_TOMONT) for x in a]

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
    ok = all(intt(ntt(a := [random.randrange(Q) for _ in range(256)])) ==
             [(x * MONT) % Q for x in a]
             for _ in range(300))
    out.append("roundtrip_to_mont(300): " + ("PASS" if ok else "FAIL"))
    ok2 = True
    for _ in range(80):
        a = [random.randrange(Q) for _ in range(256)]
        b = [random.randrange(Q) for _ in range(256)]
        na, nb = ntt(a), ntt(b)
        # Kyber-style pointwise: direct data*data Montgomery multiply, producing
        # a product that is deficient by one R; invntt_tomont repairs it.
        nc = [fqmul(na[i], nb[i]) for i in range(256)]
        if intt(nc) != schoolbook(a, b):
            ok2 = False; break
    out.append("convolution(80): " + ("PASS" if ok2 else "FAIL"))
    return ok and ok2, out

# --------------------------------------------------------------------------
# Emit constant table (consumption order) and the generated .S bodies.
# --------------------------------------------------------------------------
def s16(v):                                       # signed 16-bit literal
    v %= 65536
    return v - 65536 if v >= 32768 else v

def emit_consts(path):
    # forward zeta vectors in execution order: half0 then half1 (xb + ib carry zvec)
    fwd = []
    for prog in (FWD0, FWD1):
        for op in prog:
            if op[0] in ('xb', 'ib'):
                fwd.append(op[-1])
    # inverse consumes the reversed program per half, with INVERTED twiddles
    # (out_b = z^{-1}(a'-b')), matching run_inv_half / the proven model.
    inv = []
    for prog in (FWD0, FWD1):
        for op in reversed(prog):
            if op[0] in ('xb', 'ib'):
                inv.append([zinv_mont(z) for z in op[-1]])
    # --- Seiler precompute: each montgomery twiddle z -> (zl, zh) where
    #     zl = z*qinv mod 2^16 (gives the Montgomery factor in one multiply),
    #     zh = z (gives the product high word).  Stored interleaved zl,zh so the
    #     butterfly loads two consecutive 32-byte vectors. ---
    def zl16(z): return (z * QINV) % R          # 16-lane "low" operand value
    def vec_zl(zvec): return [zl16(z) for z in zvec]
    def vec_zh(zvec): return list(zvec)
    def pair_table(table):                       # [zvec,...] -> [zl,zh,zl,zh,...]
        out = []
        for zvec in table:
            out.append(vec_zl(zvec)); out.append(vec_zh(zvec))
        return out
    def vecs(name, table):
        lines = [f"const int16_t {name}[{len(table)}][16] = {{"]
        for v in table:
            lines.append("  {" + ",".join(str(s16(x)) for x in v) + "},")
        lines.append("};")
        return "\n".join(lines)
    def pair32(name, z):                          # 32-int16: zl[16] then zh[16]
        zl = ",".join([str(s16(zl16(z)))]*16)
        zh = ",".join([str(s16(z))]*16)
        return ("const int16_t %s[32] __attribute__((aligned(32))) = {\n  %s,  /* zl */\n  %s   /* zh */\n};\n"
                % (name, zl, zh))
    fwd_p, inv_p = pair_table(fwd), pair_table(inv)
    with open(path, "w") as f:
        f.write("#include <stdint.h>\n")
        f.write("/* Auto-generated by gen_ntt.py -- do not edit. q=%d.\n"
                "   Twiddles use the Seiler precompute: each butterfly stores zl=z*qinv\n"
                "   mod 2^16 then zh=z, as two consecutive 16-lane vectors. */\n\n" % Q)
        f.write("const int16_t ntt_qdata[48] __attribute__((aligned(32))) = {\n")
        f.write("  " + ",".join([str(Q)]*16) + ",  /* _16XQ    */\n")
        f.write("  " + ",".join([str(RMQ)]*16) + ",  /* _16XRMQ  */\n")
        f.write("  " + ",".join(["0"]*16) + "   /* _16XZERO */\n};\n\n")
        f.write(pair32("ntt_z0", Z0))
        f.write(pair32("ntt_z0inv", zinv_mont(Z0)))
        f.write(pair32("ntt_ninv", NINV_TOMONT) + "\n")
        f.write("const int16_t ntt_qinv[16] __attribute__((aligned(32))) = {\n")
        f.write("  " + ",".join([str(s16(QINV))]*16) + "\n};\n\n")
        f.write("__attribute__((aligned(32)))\n" + vecs("ntt_zetas_fwd", fwd_p) + "\n\n")
        f.write("__attribute__((aligned(32)))\n" + vecs("ntt_zetas_inv", inv_p) + "\n")
    return len(fwd_p), len(inv_p)

# --------------------------------------------------------------------------
# Emit ntt.S / invntt.S from the proven op-list (correct by construction).
# Register map: ymm0=q ymm2=2^16-q ymm3=zero ; ymm1=zh ymm15=zl (the precomputed
# twiddle pair) ; ymm4..11 = 8 data regs of a half (logical r -> ymm(4+r)) ;
# ymm12,13 = scratch ; ymm14 = butterfly temp t.
# C ABI (System V): rdi=poly, rsi=qdata[48], rdx=ztab[][16] (zl,zh pairs),
# rcx=z0 pair, r8=ninv pair.
# --------------------------------------------------------------------------
ASM_PREAMBLE = r"""/* ==========================================================================
 * Auto-generated by gen_ntt.py -- DO NOT EDIT (edit gen_ntt.py, run `make gen`).
 *
 * Complete n=256 negacyclic NTT / INTT over Z_q[x]/(x^256+1), q = 45569,
 * 16-bit AVX2, Kyber-style level merging + 4-tier shuffle network.
 *
 * Entry points (System V AMD64 ABI):
 *   ntt_avx   (rdi=poly[256] int16, rsi=qdata, rdx=zetas_fwd, rcx=z0,    r8=ninv)
 *   invntt_tomont_avx(rdi=poly[256] int16, rsi=qdata, rdx=zetas_inv, rcx=z0inv, r8=ninv)
 *   pointwise_avx(rdi=c[256], rsi=a[256], rdx=b[256], rcx=qdata)
 *
 * Register map (held live across a whole half-transform):
 *   ymm0 = q (16x)       ymm1 = zh = current zeta vec   ymm2 = 2^16 - q   ymm3 = 0
 *   ymm4..ymm11 = the 8 data registers of one 128-coeff half (logical reg r
 *                 lives in ymm(4+r); 8 regs x 16 lanes = 128 coefficients)
 *   ymm12,ymm13 = scratch   ymm14 = butterfly temp t
 *   ymm15 = zl = (zeta * q^-1) mod 2^16 = the PRECOMPUTED Montgomery-factor twiddle
 *
 * DOMAIN: coefficients stay in PLAIN [0,q); twiddles are stored *R mod q
 * (Montgomery domain), so montmul(zeta,x) = zeta*x mod q keeps x's domain.
 * Because q > 2^15 there is ZERO lazy-reduction headroom: every add/sub is a
 * full overflow-aware modular reduction (you cannot even form a+b in 16 bits,
 * since 2q > 2^16) and the Montgomery high word must be UNSIGNED (vpmulhuw).
 *
 * PERFORMANCE (Kyber/Seiler trick, eprint 2018/039 "fqmulprecomp"): instead of
 * deriving the Montgomery factor at run time (lo=in*zeta, then m=lo*qinv -- a
 * THREE-multiply dependency chain), we PRECOMPUTE per twiddle two vectors:
 *   zl = zeta*qinv mod 2^16   (the "low" operand, gives the Montgomery factor
 *                              m directly in ONE multiply)
 *   zh = zeta                 (the "high" operand, gives the product high word)
 * Both come from ntt_zetas_* (interleaved zl,zh per butterfly).  This cuts the
 * per-butterfly multiplies 4 -> 3 and the multiply dependency chain 3 -> 2,
 * relieving the p0/p1 multiply-port bottleneck (~25% fewer multiplies).
 * ========================================================================== */

/* montmul: out = (in * zeta) * R^-1 mod q, result in [0,q).
 *   zl in ymm15 (= zeta*qinv mod 2^16),  zh in ymm1 (= zeta).
 *   m    = in*zl mod 2^16                (vpmullw  -- Montgomery factor, 1 mul)
 *   out  = floor(in*zh / 2^16)           (vpmulhuw -- UNSIGNED high, q>2^15)
 *   mqhi = floor(m*q / 2^16)             (vpmulhuw)
 *   out  = out - mqhi                    (high halves; result in (-q,q))
 *   then a borrow correction: if out underflowed (Thi<mqhi) add (2^16-q) so the
 *   result lands in [0,q).  vpsubusw+vpcmpeqw+vpandn emulate the unsigned test.
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

/* dmul: GENERAL Montgomery multiply of two DATA vectors (no Seiler precompute).
 * qinv is kept in ymm9.  This is used for Kyber-style pointwise multiplication:
 * normal-domain NTT outputs multiply to an R^-1-deficient product, and
 * invntt_tomont repairs the missing factor. */
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

/* addmod: dst = (a + t) mod q, a,t in [0,q).  a+t can reach 2q-2 > 2^16, so it
 * is computed as a-(q-t)=a+t-q with a borrow-driven +(2^16-q) correction.
 * Clobbers ymm12,ymm13. */
.macro addmod dst,a,t
  vpsubw   \t,%ymm0,%ymm12
  vpsubusw \a,%ymm12,%ymm13
  vpsubw   %ymm12,\a,%ymm12
  vpcmpeqw %ymm3,%ymm13,%ymm13
  vpandn   %ymm2,%ymm13,%ymm13
  vpsubw   %ymm13,%ymm12,\dst
.endm

/* submod: dst = (a - t) mod q, a,t in [0,q), with a borrow-driven correction.
 * Clobbers ymm12,ymm13. */
.macro submod dst,a,t
  vpsubusw \a,\t,%ymm13
  vpsubw   \t,\a,%ymm12
  vpcmpeqw %ymm3,%ymm13,%ymm13
  vpandn   %ymm2,%ymm13,%ymm13
  vpsubw   %ymm13,%ymm12,\dst
.endm

/* ctbf: Cooley-Tukey butterfly (forward NTT).  t = b*zeta; b = a-t; a = a+t. */
.macro ctbf a,b
  montmul %ymm14,\b
  submod  \b,\a,%ymm14
  addmod  \a,\a,%ymm14
.endm

/* gsbf: Gentleman-Sande butterfly (inverse NTT).  t = a-b; a = a+b; b = t*zeta. */
.macro gsbf a,b
  submod  %ymm14,\a,\b
  addmod  \a,\a,\b
  montmul \b,%ymm14
.endm

/* ---- the 4-tier shuffle network (granularity halves each tier) ----
 * Each shufK takes two registers whose butterfly partners are <16 lanes apart
 * (i.e. inside one register) and interleaves them at granularity K coefficients,
 * so the partners become same-numbered lanes of TWO registers -> a plain
 * vertical butterfly works again.  Only shuf8 must cross the 128-bit lane
 * barrier (vperm2i128); the rest stay within 128 bits.  All are involutions
 * (applying twice = identity), which is what lets the inverse replay them. */
.macro shuf8 r0,r1            /* granularity 8 coeffs; the ONLY cross-128-bit tier */
  vperm2i128 $0x20,\r1,\r0,%ymm12
  vperm2i128 $0x31,\r1,\r0,\r1
  vmovdqa  %ymm12,\r0
.endm
.macro shuf4 r0,r1            /* granularity 4 coeffs (64-bit); in-128 unpack */
  vpunpcklqdq \r1,\r0,%ymm12
  vpunpckhqdq \r1,\r0,\r1
  vmovdqa  %ymm12,\r0
.endm
.macro shuf2 r0,r1            /* granularity 2 coeffs (32-bit dword) */
  vmovsldup \r1,%ymm12
  vpblendd $0xAA,%ymm12,\r0,%ymm12
  vpsrlq   $32,\r0,\r0
  vpblendd $0xAA,\r1,\r0,\r1
  vmovdqa  %ymm12,\r0
.endm
.macro shuf1 r0,r1            /* granularity 1 coeff (16-bit word) */
  vpslld   $16,\r1,%ymm12
  vpblendw $0xAA,%ymm12,\r0,%ymm12
  vpsrld   $16,\r0,\r0
  vpblendw $0xAA,\r1,\r0,\r1
  vmovdqa  %ymm12,\r0
.endm
"""

def emit_asm(fwd_path, inv_path):
    ymm = lambda r: "%%ymm%d" % (4 + r)
    def load_consts():       # ymm0=q, ymm2=2^16-q, ymm3=0 (qinv no longer needed)
        return ["  vmovdqa   0(%rsi),%ymm0", "  vmovdqa  32(%rsi),%ymm2",
                "  vmovdqa  64(%rsi),%ymm3"]
    def load_tw(reg, disp=0):   # load precomputed pair zl->ymm15, zh->ymm1 (zh 32B after zl)
        return ["  vmovdqa  %d(%s),%%ymm15" % (disp, reg),
                "  vmovdqa  %d(%s),%%ymm1"  % (disp + 32, reg)]
    def load_half(h):
        return ["  vmovdqu  %d(%%rdi),%s" % (256*h + 32*r, ymm(r)) for r in range(8)]
    def store_half(h):
        return ["  vmovdqu  %s,%d(%%rdi)" % (ymm(r), 256*h + 32*r) for r in range(8)]

    # length -> NTT level number (level 0 is the cross-half len=128 stage)
    LVL = {64: 1, 32: 2, 16: 3, 8: 4, 4: 5, 2: 6, 1: 7}
    XB_LENS = [64]*4 + [32]*4 + [16]*4   # build_half emits 4 xb per cross-reg level

    # ---- forward ----
    L = list(ASM_PREAMBLE.splitlines())
    L += [".text", ".global ntt_avx", ".global invntt_tomont_avx", ".global pointwise_avx",
          "", "/* ================= forward NTT (8 levels, len 128..1) ================= */",
          "ntt_avx:"]
    L += load_consts()
    L += ["",
          "/* --- level 0 (len=128): cross-HALF CT butterflies, partners 128 apart.",
          "       This is the ONLY level that touches memory twice; the remaining 7",
          "       levels of each 128-coeff half are then done register-resident. --- */"]
    L += load_tw("%rcx")
    for r in range(8):
        L += ["  vmovdqu  %d(%%rdi),%%ymm4" % (32*r),
              "  vmovdqu  %d(%%rdi),%%ymm5" % (256 + 32*r),
              "  ctbf %ymm4,%ymm5",
              "  vmovdqu  %%ymm4,%d(%%rdi)" % (32*r),
              "  vmovdqu  %%ymm5,%d(%%rdi)" % (256 + 32*r)]
    # per-half merged levels 1..7  (3 cross-register + 4 in-register-with-shuffle)
    ki = 0
    for h, prog in ((0, FWD0), (1, FWD1)):
        lo, hi = 128*h, 128*h + 127
        L += ["",
              "/* ===== half %d (coeffs %d..%d): levels 1-7 merged, NO stores until done ===== */" % (h, lo, hi)]
        L += load_half(h)
        xi = 0; cur = None
        for op in prog:
            if op[0] == 'xb':
                ln = XB_LENS[xi]; xi += 1
                if ln != cur:
                    cur = ln
                    L += ["  /* level %d (len=%d): cross-REGISTER CT, partners %d lanes apart, no shuffle */"
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

    # ---- inverse ----  (rdx=zinv table, rcx=z0inv, r8=ninv)
    L += ["/* ================= inverse NTT (GS butterflies, levels 7..0) ================= */",
          "invntt_tomont_avx:"]
    L += load_consts()
    ki = 0
    for h, prog in ((0, FWD0), (1, FWD1)):
        lo, hi = 128*h, 128*h + 127
        L += ["",
              "/* ===== half %d (coeffs %d..%d): levels 7-1 merged (forward order reversed) ===== */" % (h, lo, hi)]
        L += load_half(h)
        # walk the forward program backwards; track length to annotate
        rev = list(reversed(prog))
        # precompute the length of each forward op so we can label in reverse
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
            else:  # sh
                _, k, r0, r1 = op
                L += ["  shuf%d %s,%s" % (k, ymm(r0), ymm(r1))]
        L += store_half(h)
    # level 0 GS (cross-half) with z0inv
    L += ["",
          "/* --- level 0 (len=128): cross-HALF GS butterflies (z0inv) --- */"]
    L += load_tw("%rcx")
    for r in range(8):
        L += ["  vmovdqu  %d(%%rdi),%%ymm4" % (32*r),
              "  vmovdqu  %d(%%rdi),%%ymm5" % (256 + 32*r),
              "  gsbf %ymm4,%ymm5",
              "  vmovdqu  %%ymm4,%d(%%rdi)" % (32*r),
              "  vmovdqu  %%ymm5,%d(%%rdi)" % (256 + 32*r)]
    # final *1/n and to-Montgomery (Montgomery multiply by ninv_tomont)
    L += ["",
          "/* --- final scaling: multiply every coefficient by 1/n and R --- */"]
    L += load_tw("%r8")
    for r in range(16):
        L += ["  vmovdqu  %d(%%rdi),%%ymm4" % (32*r),
              "  montmul %ymm5,%ymm4",
              "  vmovdqu  %%ymm5,%d(%%rdi)" % (32*r)]
    L += ["  ret", ""]

    # pointwise multiply for the complete NTT: c_i = a_i*b_i*R^-1.
    L += ["/* ================= pointwise: c = a*b*R^-1 coefficient-wise ================= */",
          "pointwise_avx:"]
    L += ["  vmovdqa   0(%rcx),%ymm0",
          "  vmovdqa  32(%rcx),%ymm2",
          "  vmovdqa  64(%rcx),%ymm3",
          "  vmovdqa  ntt_qinv(%rip),%ymm9"]
    for r in range(16):
        off = 32*r
        L += ["  vmovdqu  %d(%%rsi),%%ymm4" % off,
              "  vmovdqu  %d(%%rdx),%%ymm5" % off,
              "  dmul     %ymm7,%ymm4,%ymm5",
              "  vmovdqu  %%ymm7,%d(%%rdi)" % off]
    L += ["  ret", ""]

    # ---- nttunpack: standard-order int32 samples -> avx NTT layout ----
    # Kyber/Dilithium-style (poly_nttunpack / nttunpack_avx): replay the forward
    # shuffle network with NO butterflies, so out[s] = src[FINAL_IDX[s]].  The
    # source samples arrive int32 in [0,q) and standard order; this narrows them
    # to uint16 and lays them out in the avx NTT slot order.
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

    # noexec stack note
    L += ['.section .note.GNU-stack,"",@progbits']

    body = "\n".join(L) + "\n"
    with open(fwd_path, "w") as f:
        f.write(body)
    return len(body.splitlines())

if __name__ == "__main__":
    ok, lines = verify()
    print("q=%d omega=%d (omega^256==-1:%s) MONT=%d QINV=%d RMQ=%d NINV_TOMONT=%d R2=%d"
          % (Q, OMEGA, pow(OMEGA,256,Q)==Q-1, MONT, QINV, RMQ, NINV_TOMONT, R2_MOD_Q))
    for l in lines:
        print(l)
    nf, ni = emit_consts("ntt_consts.c")
    print("emitted ntt_consts.c  fwd_vecs=%d inv_vecs=%d" % (nf, ni))
    nl = emit_asm("ntt.S", None)
    print("emitted ntt.S (%d lines, forward+inverse)" % nl)
    print("FWD0_ops=%d FWD1_ops=%d" % (len(FWD0), len(FWD1)))
    sys.exit(0 if ok else 1)
