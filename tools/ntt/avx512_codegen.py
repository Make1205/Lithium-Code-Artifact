#!/usr/bin/env python3
"""Shared full-ZMM AVX-512 code generator for the NTT parameter directories.

The AVX2 generators in the subdirectories keep their original output.  This
module adds an opt-in AVX-512BW path that uses the 32 architectural ZMM
registers to keep a whole n=256 or n=512 polynomial register-resident.
"""

from dataclasses import dataclass
from pathlib import Path


R = 1 << 16
VL = 32


@dataclass(frozen=True)
class Avx512Config:
    q: int
    n: int
    root: int
    brv_width: int
    levels: tuple
    scale_mont: int
    incomplete: bool
    lazy: bool
    basemul_abi: str = ""
    prefix: str = "ntt512"


def brv(x, w):
    r = 0
    for i in range(w):
        r |= ((x >> i) & 1) << (w - 1 - i)
    return r


def s16(v):
    v %= R
    return v - R if v >= (R >> 1) else v


def vec_literal(v):
    return "{" + ",".join(str(s16(x)) for x in v) + "}"


def fqmul_factory(q):
    rinv = pow(R, -1, q)

    def fqmul(a, b):
        return (a * b * rinv) % q

    return fqmul


def zeta_table(cfg):
    q = cfg.q
    mont = R % q
    return [pow(cfg.root, brv(k, cfg.brv_width), q) * mont % q
            for k in range(cfg.n)]


def fwd_k_map(cfg):
    out = {}
    k = 1
    for length in cfg.levels:
        s = 0
        while s < cfg.n:
            out[(length, s)] = k
            k += 1
            s += 2 * length
    return out


def make_zeta_of(cfg):
    zetas = zeta_table(cfg)
    fmap = fwd_k_map(cfg)

    def zeta_of(length, coeff_low):
        start = (coeff_low // (2 * length)) * (2 * length)
        return zetas[fmap[(length, start)]]

    return zeta_of


def shuffle_pair(a, b, k):
    c = []
    d = []
    for base in range(0, VL, 2 * k):
        c.extend(a[base:base + k])
        c.extend(b[base:base + k])
        d.extend(a[base + k:base + 2 * k])
        d.extend(b[base + k:base + 2 * k])
    assert len(c) == VL and len(d) == VL
    return c, d


def shuffle_indices(k, which):
    a = list(range(VL))
    b = list(range(VL, 2 * VL))
    c, d = shuffle_pair(a, b, k)
    return c if which == "c" else d


def build_prog(cfg):
    regs = cfg.n // VL
    idx = [[VL * r + l for l in range(VL)] for r in range(regs)]
    zeta_of = make_zeta_of(cfg)
    prog = []
    for length in cfg.levels:
        if length >= VL:
            step = length // VL
            for r in range(regs):
                if (r // step) % 2 == 0:
                    p = r + step
                    z = zeta_of(length, idx[r][0])
                    prog.append(("xb", length, r, p, [z] * VL))
        else:
            for m in range(regs // 2):
                r0, r1 = 2 * m, 2 * m + 1
                c, d = shuffle_pair(idx[r0], idx[r1], length)
                z = [zeta_of(length, c[l]) for l in range(VL)]
                prog.append(("sh", length, r0, r1))
                prog.append(("ib", length, r0, r1, z))
                idx[r0], idx[r1] = c, d
    return prog, idx


def zinv_mont(cfg, z):
    q = cfg.q
    rinv = pow(R, -1, q)
    plain = z * rinv % q
    return pow(plain, -1, q) * (R % q) % q


def make_pair_table(cfg, zvecs):
    qinv = pow(cfg.q, -1, R)
    out = []
    for zvec in zvecs:
        out.append([(z * qinv) % R for z in zvec])
        out.append(list(zvec))
    return out


def collect_zvecs_forward(prog):
    return [op[-1] for op in prog if op[0] in ("xb", "ib")]


def collect_zvecs_inverse(cfg, prog):
    out = []
    for op in reversed(prog):
        if op[0] in ("xb", "ib"):
            out.append([zinv_mont(cfg, z) for z in op[-1]])
    return out


def model_run(cfg, regs, prog, inverse=False, sched=None):
    q = cfg.q
    fqmul = fqmul_factory(q)
    seq = list(reversed(prog)) if inverse else list(prog)
    si = 0
    for op in seq:
        if op[0] == "sh":
            _, k, r0, r1 = op
            regs[r0], regs[r1] = shuffle_pair(regs[r0], regs[r1], k)
            continue
        _, length, r0, r1, z = op
        if inverse:
            z = [zinv_mont(cfg, x) for x in z]
        use_csub = False
        if cfg.lazy:
            use_csub = sched[si][1]
            si += 1
        a = regs[r0]
        b = regs[r1]
        for l in range(VL):
            if not inverse:
                t = fqmul(b[l] % q, z[l])
                if cfg.lazy and not use_csub:
                    a[l], b[l] = a[l] + t, a[l] + (q - t)
                else:
                    ar = a[l] % q
                    a[l], b[l] = (ar + t) % q, (ar - t) % q
            else:
                if cfg.lazy and not use_csub:
                    s = a[l] + b[l]
                    d = (a[l] % q - b[l] % q) % q
                    a[l], b[l] = s, fqmul(d, z[l])
                else:
                    ar = a[l] % q
                    br = b[l] % q
                    a[l], b[l] = (ar + br) % q, fqmul((ar - br) % q, z[l])
    return regs


def model_ntt(cfg, poly, raw=False, fwd_sched=None):
    prog, _ = build_prog(cfg)
    regs = [poly[VL * r:VL * (r + 1)] for r in range(cfg.n // VL)]
    model_run(cfg, regs, prog, inverse=False, sched=fwd_sched)
    out = [x for reg in regs for x in reg]
    return out if raw else [x % cfg.q for x in out]


def model_intt(cfg, poly, inv_sched=None):
    q = cfg.q
    fqmul = fqmul_factory(q)
    prog, _ = build_prog(cfg)
    regs = [poly[VL * r:VL * (r + 1)] for r in range(cfg.n // VL)]
    model_run(cfg, regs, prog, inverse=True, sched=inv_sched)
    out = [x for reg in regs for x in reg]
    return [fqmul(x % q, cfg.scale_mont) for x in out]


def schoolbook(cfg, a, b):
    q = cfg.q
    n = cfg.n
    c = [0] * n
    for i in range(n):
        ai = a[i]
        for j in range(n):
            k = i + j
            p = ai * b[j] % q
            if k < n:
                c[k] = (c[k] + p) % q
            else:
                c[k - n] = (c[k - n] - p) % q
    return c


def gamma_of(cfg, i):
    assert cfg.incomplete
    q = cfg.q
    return pow(cfg.root, 2 * brv(i, cfg.brv_width) + 1, q)


def basemul_model(cfg, na, nb):
    q = cfg.q
    fqmul = fqmul_factory(q)
    _, idx_after = build_prog(cfg)
    idx_flat = [x for reg in idx_after for x in reg]
    out = [0] * cfg.n
    for i in range(cfg.n // 2):
        a0, a1 = na[2 * i] % q, na[2 * i + 1] % q
        b0, b1 = nb[2 * i] % q, nb[2 * i + 1] % q
        p0 = fqmul(a0, b0)
        p1 = fqmul(a1, b1)
        x0 = fqmul(a0, b1)
        x1 = fqmul(a1, b0)
        assert idx_flat[2 * i] // 2 == idx_flat[2 * i + 1] // 2
        g = gamma_of(cfg, idx_flat[2 * i] // 2)
        out[2 * i] = (p0 + fqmul(p1, g * (R % q) % q)) % q
        out[2 * i + 1] = (x0 + x1) % q
    return out


def pointwise_model(cfg, na, nb):
    q = cfg.q
    fqmul = fqmul_factory(q)
    return [fqmul(na[i] % q, nb[i] % q) for i in range(cfg.n)]


def derive_lazy_schedules(cfg):
    if not cfg.lazy:
        return None, None

    q = cfg.q
    prog, _ = build_prog(cfg)
    regs = [[(0, q - 1) for _ in range(VL)] for _ in range(cfg.n // VL)]
    fwd = []
    for op in prog:
        if op[0] == "sh":
            _, k, r0, r1 = op
            regs[r0], regs[r1] = shuffle_pair(regs[r0], regs[r1], k)
            continue
        _, length, r0, r1, _ = op
        maxa = max(x[1] for x in regs[r0])
        maxb = max(x[1] for x in regs[r1])
        # Lazy CT is legal if a+t fits uint16 and b is a valid montmul input.
        lazy_ok = maxa + q <= R and maxb < R and maxb * (q - 1) < R * q
        csub = not lazy_ok
        fwd.append((length, csub))
        for l in range(VL):
            ba = regs[r0][l][1]
            if csub:
                regs[r0][l] = (0, q - 1)
                regs[r1][l] = (0, q - 1)
            else:
                regs[r0][l] = (0, ba + q)
                regs[r1][l] = (0, ba + q)

    out_bound = max(x[1] for reg in regs for x in reg)
    regs = [[(0, out_bound) for _ in range(VL)] for _ in range(cfg.n // VL)]
    inv = []
    for op in reversed(prog):
        if op[0] == "sh":
            _, k, r0, r1 = op
            regs[r0], regs[r1] = shuffle_pair(regs[r0], regs[r1], k)
            continue
        _, length, r0, r1, _ = op
        maxu = max(x[1] for x in regs[r0])
        maxv = max(x[1] for x in regs[r1])
        d_bound = maxu + q
        lazy_ok = maxu + maxv <= R and d_bound <= R and d_bound * (q - 1) < R * q
        csub = not lazy_ok
        inv.append((length, csub))
        for l in range(VL):
            if csub:
                regs[r0][l] = (0, q - 1)
                regs[r1][l] = (0, q - 1)
            else:
                regs[r0][l] = (0, maxu + maxv)
                regs[r1][l] = (0, q - 1)
    return fwd, inv


ASM_PREAMBLE = r"""/* Auto-generated AVX-512BW full-ZMM NTT.  Edit gen_ntt.py. */
.text
.global ntt_avx512
.global invntt_tomont_avx512
"""


def asm_macros(lazy=False, include_dmul=False):
    csub = r"""
.macro csub512 a
  vpcmpuw $1,%zmm16,\a,%k1
  vpsubw  %zmm16,\a,\a
  vpaddw  %zmm16,\a,\a{%k1}
.endm
"""
    common = r"""
.macro load_consts512 qptr
  vmovdqa64 0(\qptr),%zmm16
  vmovdqa64 64(\qptr),%zmm17
.endm
.macro load_tw512 base,disp
  vmovdqa64 \disp(\base),%zmm18
  vmovdqa64 \disp+64(\base),%zmm19
.endm
.macro montmul512 out,in
  vpmullw  %zmm18,\in,%zmm20
  vpmulhuw %zmm19,\in,\out
  vpmulhuw %zmm16,%zmm20,%zmm20
  vpcmpuw  $1,%zmm20,\out,%k1
  vpsubw   %zmm20,\out,\out
  vpsubw   %zmm17,\out,\out{%k1}
.endm
.macro montmul512_pair out,in,zl,zh
  vpmullw  \zl,\in,%zmm20
  vpmulhuw \zh,\in,\out
  vpmulhuw %zmm16,%zmm20,%zmm20
  vpcmpuw  $1,%zmm20,\out,%k1
  vpsubw   %zmm20,\out,\out
  vpsubw   %zmm17,\out,\out{%k1}
.endm
.macro addmod512 dst,a,t
  vpsubw  \t,%zmm16,%zmm21
  vpcmpuw $1,%zmm21,\a,%k1
  vpsubw  %zmm21,\a,\dst
  vpsubw  %zmm17,\dst,\dst{%k1}
.endm
.macro submod512 dst,a,t
  vpcmpuw $1,\t,\a,%k1
  vpsubw  \t,\a,\dst
  vpsubw  %zmm17,\dst,\dst{%k1}
.endm
.macro ctbf512 a,b
  montmul512 %zmm22,\b
  submod512  \b,\a,%zmm22
  addmod512  \a,\a,%zmm22
.endm
.macro gsbf512 a,b
  submod512  %zmm22,\a,\b
  addmod512  \a,\a,\b
  montmul512 \b,%zmm22
.endm
.macro ctbf512_lazy a,b
  montmul512 %zmm22,\b
  vpsubw  %zmm22,%zmm16,%zmm21
  vpaddw  %zmm21,\a,\b
  vpaddw  %zmm22,\a,\a
.endm
.macro ctbf512_csub a,b
  montmul512 %zmm22,\b
  csub512 \a
  vpsubw  %zmm22,%zmm16,%zmm21
  vpaddw  %zmm21,\a,\b
  csub512 \b
  vpaddw  %zmm22,\a,\a
  csub512 \a
.endm
.macro gsbf512_lazy a,b
  vpsubw  \b,%zmm16,%zmm22
  vpaddw  \a,%zmm22,%zmm22
  vpaddw  \b,\a,\a
  montmul512 \b,%zmm22
.endm
.macro gsbf512_csub a,b
  csub512 \a
  csub512 \b
  vpsubw  \b,%zmm16,%zmm22
  vpaddw  \a,%zmm22,%zmm22
  vpaddw  \b,\a,\a
  csub512 \a
  montmul512 \b,%zmm22
.endm
"""
    dmul = r"""
.macro dmul512 out,x,y
  vpmullw  \x,\y,%zmm20
  vpmullw  %zmm28,%zmm20,%zmm20
  vpmulhuw \x,\y,\out
  vpmulhuw %zmm16,%zmm20,%zmm20
  vpcmpuw  $1,%zmm20,\out,%k1
  vpsubw   %zmm20,\out,\out
  vpsubw   %zmm17,\out,\out{%k1}
.endm
"""
    shuf = ""
    for k in (16, 8, 4, 2, 1):
        shuf += f""".macro shuf{k}_512 r0,r1
  vmovdqa64 \\r0,%zmm31
  vmovdqa64 shuf{k}_c(%rip),%zmm30
  vpermt2w \\r1,%zmm30,%zmm31
  vmovdqa64 shuf{k}_d(%rip),%zmm30
  vpermt2w \\r1,%zmm30,\\r0
  vmovdqa64 \\r0,\\r1
  vmovdqa64 %zmm31,\\r0
.endm
"""
    return common + (csub if lazy else "") + (dmul if include_dmul else "") + shuf


def emit_consts(cfg, prog, idx_after, path):
    q = cfg.q
    qinv = pow(q, -1, R)
    rmq = R - q
    fwd = make_pair_table(cfg, collect_zvecs_forward(prog))
    inv = make_pair_table(cfg, collect_zvecs_inverse(cfg, prog))
    qdata = [[q] * VL, [rmq] * VL]
    scale = make_pair_table(cfg, [[cfg.scale_mont] * VL])
    lines = ["#include <stdint.h>",
             f"/* Auto-generated AVX-512 constants. q={q} n={cfg.n}. */", ""]
    lines.append(f"const int16_t {cfg.prefix}_qdata[64] __attribute__((aligned(64))) = {{")
    lines.append("  " + ",".join(str(s16(x)) for x in qdata[0]) + ",")
    lines.append("  " + ",".join(str(s16(x)) for x in qdata[1]))
    lines.append("};\n")
    lines.append(f"const int16_t {cfg.prefix}_scale[64] __attribute__((aligned(64))) = {{")
    lines.append("  " + ",".join(str(s16(x)) for x in scale[0]) + ",")
    lines.append("  " + ",".join(str(s16(x)) for x in scale[1]))
    lines.append("};\n")
    for k in (16, 8, 4, 2, 1):
        for which in ("c", "d"):
            lines.append(f"const uint16_t shuf{k}_{which}[32] __attribute__((aligned(64))) = "
                         + vec_literal(shuffle_indices(k, which)) + ";\n")
    lines.append(f"__attribute__((aligned(64)))\nconst int16_t {cfg.prefix}_zetas_fwd[{len(fwd)}][32] = {{")
    lines.extend("  " + vec_literal(v) + "," for v in fwd)
    lines.append("};\n")
    lines.append(f"__attribute__((aligned(64)))\nconst int16_t {cfg.prefix}_zetas_inv[{len(inv)}][32] = {{")
    lines.extend("  " + vec_literal(v) + "," for v in inv)
    lines.append("};\n")

    if cfg.incomplete:
        r2 = (R % q) * (R % q) % q
        r2pair = make_pair_table(cfg, [[r2] * VL])
        lines.append(f"const int16_t {cfg.prefix}_qinv[32] __attribute__((aligned(64))) = "
                     + vec_literal([qinv] * VL) + ";\n")
        lines.append(f"const int16_t {cfg.prefix}_r2[64] __attribute__((aligned(64))) = {{")
        lines.append("  " + ",".join(str(s16(x)) for x in r2pair[0]) + ",")
        lines.append("  " + ",".join(str(s16(x)) for x in r2pair[1]))
        lines.append("};\n")
        gammas = []
        for r, reg in enumerate(idx_after):
            gvec = []
            for l, coef in enumerate(reg):
                if l % 2 == 0:
                    assert coef // 2 == reg[l + 1] // 2
                gvec.append(gamma_of(cfg, coef // 2) * (R % q) % q)
            gammas.extend(make_pair_table(cfg, [gvec]))
        lines.append(f"__attribute__((aligned(64)))\nconst int16_t {cfg.prefix}_basemul_gamma[{len(gammas)}][32] = {{")
        lines.extend("  " + vec_literal(v) + "," for v in gammas)
        lines.append("};\n")
    else:
        r2 = (R % q) * (R % q) % q
        r2pair = make_pair_table(cfg, [[r2] * VL])
        lines.append(f"const int16_t {cfg.prefix}_qinv[32] __attribute__((aligned(64))) = "
                     + vec_literal([qinv] * VL) + ";\n")
        lines.append(f"const int16_t {cfg.prefix}_r2[64] __attribute__((aligned(64))) = {{")
        lines.append("  " + ",".join(str(s16(x)) for x in r2pair[0]) + ",")
        lines.append("  " + ",".join(str(s16(x)) for x in r2pair[1]))
        lines.append("};\n")

    Path(path).write_text("\n".join(lines), encoding="ascii")
    return len(fwd), len(inv)


def zmm(r):
    return f"%zmm{r}"


def emit_asm(cfg, prog, path, fwd_sched=None, inv_sched=None):
    regs = cfg.n // VL
    lines = [ASM_PREAMBLE.rstrip()]
    if cfg.incomplete:
        lines.append(".global basemul_avx512")
    else:
        lines.append(".global pointwise_avx512")
    lines.append(asm_macros(cfg.lazy, include_dmul=True))

    def load_all(base="%rdi"):
        return [f"  vmovdqu16 {64*r}({base}),{zmm(r)}" for r in range(regs)]

    def store_all(base="%rdi"):
        return [f"  vmovdqu16 {zmm(r)},{64*r}({base})" for r in range(regs)]

    def load_tw(base, idx):
        disp = idx * 128
        return [f"  vmovdqa64 {disp}({base}),%zmm18",
                f"  vmovdqa64 {disp + 64}({base}),%zmm19"]

    # Forward.
    lines += ["", "ntt_avx512:", "  load_consts512 %rsi"]
    lines += load_all("%rdi")
    ki = 0
    si = 0
    for op in prog:
        if op[0] == "sh":
            _, k, r0, r1 = op
            lines.append(f"  shuf{k}_512 {zmm(r0)},{zmm(r1)}")
        else:
            _, _length, r0, r1, _z = op
            lines += load_tw("%rdx", ki)
            if cfg.lazy:
                csub = fwd_sched[si][1]
                si += 1
                macro = "ctbf512_csub" if csub else "ctbf512_lazy"
            else:
                macro = "ctbf512"
            lines.append(f"  {macro} {zmm(r0)},{zmm(r1)}")
            ki += 1
    lines += store_all("%rdi")
    lines += ["  vzeroupper", "  ret", ""]

    # Inverse.
    lines += ["invntt_tomont_avx512:", "  load_consts512 %rsi"]
    lines += load_all("%rdi")
    ki = 0
    si = 0
    for op in reversed(prog):
        if op[0] == "sh":
            _, k, r0, r1 = op
            lines.append(f"  shuf{k}_512 {zmm(r0)},{zmm(r1)}")
        else:
            _, _length, r0, r1, _z = op
            lines += load_tw("%rdx", ki)
            if cfg.lazy:
                csub = inv_sched[si][1]
                si += 1
                macro = "gsbf512_csub" if csub else "gsbf512_lazy"
            else:
                macro = "gsbf512"
            lines.append(f"  {macro} {zmm(r0)},{zmm(r1)}")
            ki += 1
    lines += ["  vmovdqa64 0(%rcx),%zmm18", "  vmovdqa64 64(%rcx),%zmm19"]
    for r in range(regs):
        if cfg.lazy:
            lines.append(f"  csub512 {zmm(r)}")
        lines.append(f"  montmul512 {zmm(r)},{zmm(r)}")
    lines += store_all("%rdi")
    lines += ["  vzeroupper", "  ret", ""]

    if cfg.incomplete:
        emit_basemul_asm(cfg, lines)
    else:
        emit_pointwise_asm(cfg, lines)

    emit_nttunpack_asm(cfg, prog, lines)
    emit_nttunpack16_asm(cfg, prog, lines)
    emit_input_conversion_asm512(cfg, lines)

    lines.append('.section .note.GNU-stack,"",@progbits')
    Path(path).write_text("\n".join(lines) + "\n", encoding="ascii")
    return len(lines)


def emit_input_conversion_asm512(cfg, lines):
    rec = (1 << 47) // cfg.q
    rec_lo = rec & 0xffffffff
    rec_hi = rec >> 32
    loops = cfg.n // 16
    lines += [".text",
              "",
              "/* ============ NTT input conversion helpers (standard order) ============",
              " * C ABI:",
              " *   poly_freeze_to_u16_avx512(dst:int16_t*, src:int32_t*)",
              " *   poly_freeze_small_to_u16_avx512(dst:int16_t*, src:int32_t*)",
              " *   poly_cast_modq_avx512(dst:int16_t*, src:int32_t*)",
              " * These helpers feed ntt_avx512's standard-order input buffer; unlike",
              " * nttunpack_avx512, they do not apply the AVX-512 NTT slot shuffle ladder. */",
              ".p2align 5",
              ".global poly_freeze_to_u16_avx512",
              "poly_freeze_to_u16_avx512:",
              "  vpbroadcastd .Lfreeze_q_avx512(%rip),%zmm7",
              "  vpbroadcastd .Lfreeze_rec_lo_avx512(%rip),%zmm5",
              "  vpbroadcastd .Lfreeze_rec_hi_avx512(%rip),%zmm8",
              f"  mov ${loops},%ecx",
              ".Lfreeze_full_loop_avx512:",
              "  vmovdqu32 (%rsi),%zmm0",
              "  vpsrad    $31,%zmm0,%zmm1",
              "  vpxord    %zmm1,%zmm0,%zmm2",
              "  vpsubd    %zmm1,%zmm2,%zmm2",
              "  vpmuludq  %zmm5,%zmm2,%zmm3",
              "  vpsrlq    $47,%zmm3,%zmm3",
              "  vpsrlq    $32,%zmm2,%zmm4",
              "  vpmuludq  %zmm5,%zmm4,%zmm4",
              "  vpsrlq    $47,%zmm4,%zmm4",
              "  vpsllq    $32,%zmm4,%zmm4",
              "  vpord     %zmm4,%zmm3,%zmm3",
              "  vpsrld    $15,%zmm2,%zmm4",
              "  vpmulld   %zmm8,%zmm4,%zmm4",
              "  vpaddd    %zmm4,%zmm3,%zmm3",
              "  vpmulld   %zmm7,%zmm3,%zmm3",
              "  vpsubd    %zmm3,%zmm2,%zmm2",
              "  vpsubd    %zmm7,%zmm2,%zmm2",
              "  vpsrad    $31,%zmm2,%zmm4",
              "  vpandd    %zmm7,%zmm4,%zmm4",
              "  vpaddd    %zmm4,%zmm2,%zmm2",
              "  vpxord    %zmm1,%zmm2,%zmm2",
              "  vpsubd    %zmm1,%zmm2,%zmm2",
              "  vpsrad    $31,%zmm2,%zmm4",
              "  vpandd    %zmm7,%zmm4,%zmm4",
              "  vpaddd    %zmm4,%zmm2,%zmm2",
              "  vpmovdw   %zmm2,%ymm0",
              "  vmovdqu   %ymm0,(%rdi)",
              "  add       $64,%rsi",
              "  add       $32,%rdi",
              "  dec       %ecx",
              "  jnz       .Lfreeze_full_loop_avx512",
              "  ret",
              "",
              ".p2align 5",
              ".global poly_freeze_small_to_u16_avx512",
              "poly_freeze_small_to_u16_avx512:",
              "  vpbroadcastd .Lfreeze_q_avx512(%rip),%zmm7",
              f"  mov ${loops},%ecx",
              ".Lfreeze_small_loop_avx512:",
              "  vmovdqu32 (%rsi),%zmm0",
              "  vpsrad    $31,%zmm0,%zmm1",
              "  vpandd    %zmm7,%zmm1,%zmm1",
              "  vpaddd    %zmm1,%zmm0,%zmm0",
              "  vpmovdw   %zmm0,%ymm0",
              "  vmovdqu   %ymm0,(%rdi)",
              "  add       $64,%rsi",
              "  add       $32,%rdi",
              "  dec       %ecx",
              "  jnz       .Lfreeze_small_loop_avx512",
              "  ret",
              "",
              ".p2align 5",
              ".global poly_cast_modq_avx512",
              "poly_cast_modq_avx512:",
              f"  mov ${loops},%ecx",
              ".Lfreeze_cast_loop_avx512:",
              "  vmovdqu32 (%rsi),%zmm0",
              "  vpmovdw   %zmm0,%ymm0",
              "  vmovdqu   %ymm0,(%rdi)",
              "  add       $64,%rsi",
              "  add       $32,%rdi",
              "  dec       %ecx",
              "  jnz       .Lfreeze_cast_loop_avx512",
              "  ret",
              "",
              ".section .rodata",
              ".p2align 4",
              f".Lfreeze_q_avx512:       .long {cfg.q}",
              f".Lfreeze_rec_lo_avx512:  .long {rec_lo}",
              f".Lfreeze_rec_hi_avx512:  .long {rec_hi}",
              ".text"]

def emit_basemul_asm(cfg, lines):
    regs = cfg.n // VL
    if cfg.basemul_abi == "qdata_second":
        c_reg, q_reg, a_reg, b_reg, g_reg = "%rdi", "%rsi", "%rdx", "%rcx", "%r8"
    else:
        c_reg, a_reg, b_reg, q_reg, g_reg = "%rdi", "%rsi", "%rdx", "%rcx", "%r8"
    lines += ["basemul_avx512:", f"  load_consts512 {q_reg}",
              f"  vmovdqa64 {cfg.prefix}_qinv(%rip),%zmm28",
              "  mov $0xaaaaaaaa,%eax", "  kmovd %eax,%k3"]
    for r in range(regs):
        lines.append(f"  vmovdqu16 {64*r}({a_reg}),{zmm(r)}")
        lines.append(f"  vmovdqu16 {64*r}({b_reg}),{zmm(8+r)}")
    if cfg.lazy:
        for r in range(regs):
            lines.append(f"  csub512 {zmm(r)}")
            lines.append(f"  csub512 {zmm(8+r)}")
    for r in range(regs):
        lines += ["  vprold $16," + zmm(8 + r) + ",%zmm21",
                  f"  dmul512 %zmm22,{zmm(r)},{zmm(8+r)}",
                  f"  dmul512 %zmm23,{zmm(r)},%zmm21",
                  f"  vmovdqa64 {128*r}({g_reg}),%zmm18",
                  f"  vmovdqa64 {128*r+64}({g_reg}),%zmm19",
                  "  montmul512 %zmm24,%zmm22",
                  "  vprold $16,%zmm24,%zmm25",
                  "  addmod512 %zmm24,%zmm22,%zmm25",
                  "  vprold $16,%zmm23,%zmm25",
                  "  addmod512 %zmm23,%zmm23,%zmm25",
                  "  vmovdqa64 %zmm24,%zmm26",
                  "  vmovdqu16 %zmm23,%zmm26{%k3}",
                  f"  vmovdqu16 %zmm26,{64*r}({c_reg})"]
    lines += ["  vzeroupper", "  ret", ""]


def emit_nttunpack_asm(cfg, prog, lines):
    # nttunpack: standard-order int32 samples -> avx512 NTT layout.
    # Kyber/Dilithium-style (poly_nttunpack): replay the forward shuffle network
    # with NO butterflies, so out[s] = src[FINAL_IDX[s]].  Narrows int32 [0,q) ->
    # int16 on load (vpmovdw).  rdi = dst (int16[n]), rsi = src (int32[n]).
    regs = cfg.n // VL
    lines += ["",
              "/* nttunpack: standard-order int32 samples -> avx512 NTT layout.",
              " * Kyber-style forward shuffle ladder, no butterflies; narrows int32 [0,q) -> int16.",
              " *   rdi = dst (int16[n]), rsi = src (int32[n]). */",
              ".global nttunpack_avx512", "nttunpack_avx512:"]
    for r in range(regs):
        so = 128 * r  # 32 int32 per reg (VL coeffs) = 128 bytes; two 16-lane halves
        lines += [f"  vmovdqu32 {so}(%rsi),%zmm20",
                  f"  vpmovdw %zmm20,%ymm{r}",
                  f"  vmovdqu32 {so + 64}(%rsi),%zmm20",
                  f"  vpmovdw %zmm20,%ymm21",
                  f"  vinserti64x4 $1,%ymm21,%zmm{r},%zmm{r}"]
    for op in prog:
        if op[0] == "sh":
            _, k, r0, r1 = op
            lines.append(f"  shuf{k}_512 {zmm(r0)},{zmm(r1)}")
    for r in range(regs):
        lines.append(f"  vmovdqu16 {zmm(r)},{64 * r}(%rdi)")
    lines += ["  vzeroupper", "  ret", ""]


def emit_nttunpack16_asm(cfg, prog, lines):
    # nttunpack16: standard-order UINT16 samples -> avx512 NTT layout.  Same
    # shuffle ladder as nttunpack_avx512, but the source is ALREADY int16 (in
    # [0,q)), so load it directly -- no int32 widen / vpmovdw narrow.  Lets
    # poly_ntt_import unpack a poly16 in place without an int32 staging copy.
    # In-place safe: the whole poly is loaded into zmm regs before any store.
    regs = cfg.n // VL
    lines += ["",
              "/* nttunpack16: standard-order uint16 samples -> avx512 NTT layout.",
              " * Same shuffle ladder as nttunpack_avx512 but src is ALREADY int16",
              " * (in [0,q)): load directly, no int32 widen / vpmovdw.  In-place safe",
              " * (whole poly loaded before any store).",
              " *   rdi = dst (int16[n]), rsi = src (int16[n]). */",
              ".global nttunpack16_avx512", "nttunpack16_avx512:"]
    for r in range(regs):
        lines.append(f"  vmovdqu16 {64 * r}(%rsi),{zmm(r)}")
    for op in prog:
        if op[0] == "sh":
            _, k, r0, r1 = op
            lines.append(f"  shuf{k}_512 {zmm(r0)},{zmm(r1)}")
    for r in range(regs):
        lines.append(f"  vmovdqu16 {zmm(r)},{64 * r}(%rdi)")
    lines += ["  vzeroupper", "  ret", ""]


def emit_pointwise_asm(cfg, lines):
    regs = cfg.n // VL
    lines += ["pointwise_avx512:", "  load_consts512 %rcx",
              f"  vmovdqa64 {cfg.prefix}_qinv(%rip),%zmm28"]
    for r in range(regs):
        lines += [f"  vmovdqu16 {64*r}(%rsi),%zmm0",
                  f"  vmovdqu16 {64*r}(%rdx),%zmm1",
                  "  dmul512 %zmm21,%zmm0,%zmm1",
                  f"  vmovdqu16 %zmm21,{64*r}(%rdi)"]
    lines += ["  vzeroupper", "  ret", ""]


def verify(cfg, fwd_sched=None, inv_sched=None):
    import random

    random.seed(20260601 + cfg.q + cfg.n)
    ok = True
    for _ in range(80):
        a = [random.randrange(cfg.q) for _ in range(cfg.n)]
        mont = R % cfg.q
        want = [(x * mont) % cfg.q for x in a]
        if model_intt(cfg, model_ntt(cfg, a, raw=True, fwd_sched=fwd_sched), inv_sched=inv_sched) != want:
            ok = False
            break
    ok2 = True
    for _ in range(24 if cfg.n == 512 else 50):
        a = [random.randrange(cfg.q) for _ in range(cfg.n)]
        b = [random.randrange(cfg.q) for _ in range(cfg.n)]
        na = model_ntt(cfg, a, raw=True, fwd_sched=fwd_sched)
        nb = model_ntt(cfg, b, raw=True, fwd_sched=fwd_sched)
        nc = basemul_model(cfg, na, nb) if cfg.incomplete else pointwise_model(cfg, na, nb)
        if model_intt(cfg, nc, inv_sched=inv_sched) != schoolbook(cfg, a, b):
            ok2 = False
            break
    return ok and ok2, [f"avx512 roundtrip-to-mont: {'PASS' if ok else 'FAIL'}",
                       f"avx512 convolution: {'PASS' if ok2 else 'FAIL'}"]


def emit_avx512(cfg, out_dir="."):
    prog, idx_after = build_prog(cfg)
    fwd_sched, inv_sched = derive_lazy_schedules(cfg)
    ok, lines = verify(cfg, fwd_sched=fwd_sched, inv_sched=inv_sched)
    nf, ni = emit_consts(cfg, prog, idx_after, Path(out_dir) / "ntt_consts_avx512.c")
    nl = emit_asm(cfg, prog, Path(out_dir) / "ntt_avx512.S", fwd_sched, inv_sched)
    sched_lines = []
    if cfg.lazy:
        f = sorted(set(length for length, csub in fwd_sched if csub))
        i = sorted(set(length for length, csub in inv_sched if csub))
        sched_lines.append(f"avx512 lazy fwd CSUB lengths: {f}")
        sched_lines.append(f"avx512 lazy inv CSUB lengths: {i}")
    return ok, lines + sched_lines + [
        f"emitted ntt_consts_avx512.c fwd_vecs={nf} inv_vecs={ni}",
        f"emitted ntt_avx512.S ({nl} lines)",
    ]
