#!/usr/bin/env python3
"""Compile and objdump Lithium secret-handling objects for forbidden CT opcodes."""
from __future__ import annotations

import argparse
import datetime as _dt
import os
import re
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

QSETS = {120: "q57601n256", 180: "q30977n256", 260: "q64513n512"}
BACKENDS = ("ref", "avx2", "avx512")
MODES = (120, 180, 260)
CCS = ("gcc", "clang", "cc")
OPTS = ("-O3", "-Os")

BAD_MNEMONIC_RE = re.compile(
    r"^(?:"
    r"idiv.*|div.*|sdiv.*|udiv.*|vdiv.*|"
    r"sqrt.*|vsqrt.*|"
    r".*gather.*|.*scatter.*|"
    r"v?cvtsi2sd.*|v?cvttsd2si.*|v?cvtsd2si.*|v?comisd.*|ucomisd.*"
    r")$",
    re.IGNORECASE,
)

ALLOW_GATHER_SOURCES = {
    ("avx2", "fips202x4.c"): "public SHAKE lane pointers/offsets only",
    ("avx512", "fips202x8.c"): "public SHAKE lane pointers/offsets only",
}

ALLOW_PUBLIC_DIVS = [
    (
        ("ref", "avx2", "avx512"),
        {"sampler.c"},
        re.compile(r"^(?:sample_gauss_N|sample_gauss_N\.part\.\d+)$"),
        "public/fixed Gaussian length block-count arithmetic",
    ),
    (
        ("ref", "avx2", "avx512"),
        {"fips202.c"},
        re.compile(r".*_shake(?:128|256)$"),
        "public one-shot SHAKE output-length block count",
    ),
    (
        ("avx2",),
        {"fips202x4.c"},
        re.compile(r".*_shake(?:128|256)x4$"),
        "public four-lane SHAKE output-length block count",
    ),
    (
        ("avx512",),
        {"fips202x8.c"},
        re.compile(r".*_shake(?:128|256)x8$"),
        "public eight-lane SHAKE output-length block count",
    ),
]

FUNC_LABEL_RE = re.compile(r"^[0-9a-fA-F]+ <([^>]+)>:$")

BASE_CFLAGS = ["-fstack-protector-strong", "-Wall", "-Wextra", "-std=c99"]
REF_SRCS = [
    "sampler.c",
    "polyvec.c",
    "approx_exp.c",
    "irs.c",
    "sign.c",
    "packing.c",
    "reduce.c",
    "rounding.c",
    "fips202.c",
    "symmetric-shake.c",
    "ntt.c",
    "irs_rans.c",
]
AVX2_SRCS = REF_SRCS + ["fips202x4.c", "f1600x4.S"]
AVX512_SRCS = REF_SRCS + ["fips202x8.c", "keccakf1600x8.c"]


@dataclass(frozen=True)
class Variant:
    name: str
    extra_cflags: tuple[str, ...]
    extra_sources: tuple[str, ...]


def run(cmd: list[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, cwd=str(cwd), text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def cc_version(cc: str) -> str:
    p = subprocess.run([cc, "--version"], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    first = (p.stdout or p.stderr).splitlines()
    return first[0] if first else cc


def backend_cflags(backend: str, opt: str, variant: Variant) -> list[str]:
    if backend == "ref":
        return [opt, *BASE_CFLAGS, "-I.", "-Itest", *variant.extra_cflags]
    if backend == "avx2":
        return [
            opt,
            *BASE_CFLAGS,
            "-I.",
            "-I../ref",
            "-I../ref/test",
            "-mavx2",
            "-mbmi2",
            "-mpopcnt",
            "-march=native",
            "-mtune=native",
            "-mno-avx512f",
            "-DUSE_AVX2_NTT",
            "-DUSE_AVX2_SAMPLER",
            "-DUSE_AVX2_SHAKE4X",
            *variant.extra_cflags,
        ]
    if backend == "avx512":
        return [
            opt,
            *BASE_CFLAGS,
            "-I.",
            "-I../ref",
            "-I../ref/test",
            "-mavx512f",
            "-mavx512bw",
            "-mavx512vl",
            "-mavx512dq",
            "-mavx512ifma",
            "-mbmi2",
            "-mpopcnt",
            "-march=native",
            "-mtune=native",
            "-DUSE_AVX512_NTT",
            "-DUSE_AVX512_SAMPLER",
            "-DUSE_AVX512_SHAKE8X",
            *variant.extra_cflags,
        ]
    raise ValueError(backend)


def backend_sources(backend: str, mode: int, variant: Variant) -> list[str]:
    if backend == "ref":
        base = REF_SRCS
    elif backend == "avx2":
        base = AVX2_SRCS + [f"{QSETS[mode]}/ntt.S"]
    elif backend == "avx512":
        base = AVX512_SRCS + [f"{QSETS[mode]}/ntt_avx512.S"]
    else:
        raise ValueError(backend)
    out: list[str] = []
    for src in [*base, *variant.extra_sources]:
        if src not in out:
            out.append(src)
    return out


def variants_for(backend: str, include_optional: bool) -> list[Variant]:
    variants = [Variant("default", (), ())]
    if include_optional and backend == "avx2":
        variants.append(Variant("avx2exp", ("-DAVX2_EXP", "-ffp-contract=off", "-fno-fast-math"), ("approx_exp_avx2.c",)))
    if include_optional and backend == "avx512":
        variants.append(Variant("ifmaexp", ("-DAVX512IFMA_EXP",), ("approx_exp_avx512ifma.c",)))
    return variants


def extract_mnemonic(line: str) -> str | None:
    if "\t" not in line:
        return None
    fields = [f for f in line.split("\t") if f.strip()]
    if len(fields) < 2:
        return None
    asm = fields[-1].strip()
    if not asm or asm.endswith(":"):
        return None
    return asm.split(None, 1)[0].lower()


def public_div_reason(backend: str, src: str, func: str | None, mnemonic: str) -> str | None:
    if not func or not re.match(r"^(?:i?div|sdiv|udiv)", mnemonic, re.IGNORECASE):
        return None
    src_name = Path(src).name
    func_base = func.rsplit("/", 1)[-1]
    for backends, sources, func_re, reason in ALLOW_PUBLIC_DIVS:
        if backend in backends and src_name in sources and func_re.match(func_base):
            return reason
    return None


def scan_objdump(text: str, backend: str, src: str) -> tuple[list[str], list[str]]:
    forbidden: list[str] = []
    allowed: list[str] = []
    current_func: str | None = None
    for line in text.splitlines():
        label = FUNC_LABEL_RE.match(line.strip())
        if label:
            current_func = label.group(1)
            continue
        mnemonic = extract_mnemonic(line)
        if not mnemonic or not BAD_MNEMONIC_RE.match(mnemonic):
            continue
        if "gather" in mnemonic.lower() and (backend, Path(src).name) in ALLOW_GATHER_SOURCES:
            allowed.append(line.rstrip())
            continue
        div_reason = public_div_reason(backend, src, current_func, mnemonic)
        if div_reason:
            allowed.append(f"{line.rstrip()}  # {div_reason}")
            continue
        forbidden.append(line.rstrip())
    return forbidden, allowed


def compile_and_scan(backend: str, mode: int, cc: str, opt: str, variant: Variant, objdir: Path) -> tuple[bool, list[str]]:
    cwd = ROOT / backend
    cflags = backend_cflags(backend, opt, variant)
    lines = [f"## backend={backend} mode={mode} cc={cc} opt={opt} variant={variant.name}"]
    lines.append(f"cc_version: {cc_version(cc)}")
    fail = False
    for src in backend_sources(backend, mode, variant):
        src_path = cwd / src
        if not src_path.exists():
            lines.append(f"FAIL missing source: {src}")
            fail = True
            continue
        stem = src.replace("/", "_").replace(".", "_")
        obj = objdir / f"{backend}_{variant.name}_{cc}_{opt.replace('-', '')}_{mode}_{stem}.o"
        cmd = [cc, *cflags, f"-DLITHIUM_MODE={mode}", "-c", src, "-o", str(obj)]
        p = run(cmd, cwd)
        if p.returncode != 0:
            lines.append(f"FAIL compile {src}: {' '.join(shlex.quote(x) for x in cmd)}")
            lines.extend((p.stderr or p.stdout).splitlines()[:40])
            fail = True
            continue
        d = run(["objdump", "-d", str(obj)], cwd)
        if d.returncode != 0:
            lines.append(f"FAIL objdump {src}: {d.stderr.strip()}")
            fail = True
            continue
        forbidden, allowed = scan_objdump(d.stdout, backend, src)
        if forbidden:
            lines.append(f"FAIL {src}: forbidden CT mnemonics")
            lines.extend(f"  {x}" for x in forbidden[:60])
            if len(forbidden) > 60:
                lines.append(f"  ... {len(forbidden) - 60} more")
            fail = True
        else:
            allow_note = ""
            if allowed:
                gather_count = sum(1 for x in allowed if "gather" in x.lower())
                public_div_count = len(allowed) - gather_count
                notes = []
                if gather_count:
                    reason = ALLOW_GATHER_SOURCES[(backend, Path(src).name)]
                    notes.append(f"allowlisted_gather={gather_count} ({reason})")
                if public_div_count:
                    notes.append(f"allowlisted_public_div={public_div_count}")
                allow_note = "; " + "; ".join(notes)
            lines.append(f"PASS {src}{allow_note}")
    lines.append("RESULT: " + ("FAIL" if fail else "PASS"))
    return not fail, lines


def parse_list(raw: str, allowed: tuple) -> list:
    vals = raw.replace(",", " ").split()
    if vals == ["all"]:
        return list(allowed)
    cast = type(allowed[0])
    out = [cast(v) for v in vals]
    bad = [v for v in out if v not in allowed]
    if bad:
        raise SystemExit(f"unsupported values: {bad}; allowed={allowed}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--backends", default="all")
    ap.add_argument("--modes", default="all")
    ap.add_argument("--ccs", default="gcc")
    ap.add_argument("--opts", default="-O3")
    ap.add_argument("--include-optional", action="store_true")
    ap.add_argument("--report", default=None)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        sample = "0000000000000000 <sample_gauss_N.part.0>:\n   0:\t48 f7 f1             \tdiv    %rcx\n0000000000000003 <x>:\n   3:\tc5 fd 92 00          \tvgatherdps (%rax,%ymm0,4),%ymm1\n   8:\tc5 fb 2a c0          \tvcvtsi2sd %rax,%xmm0,%xmm0\n"
        bad, allowed = scan_objdump(sample, "ref", "x.c")
        assert len(bad) == 3 and not allowed, (bad, allowed)
        bad, allowed = scan_objdump(sample, "ref", "sampler.c")
        assert len(bad) == 2 and len(allowed) == 1, (bad, allowed)
        bad, allowed = scan_objdump(sample, "avx2", "fips202x4.c")
        assert len(bad) == 2 and len(allowed) == 1, (bad, allowed)
        print("ct_scan.py self-test: PASS")
        return 0

    backends = parse_list(args.backends, BACKENDS)
    modes = parse_list(args.modes, MODES)
    ccs = parse_list(args.ccs, CCS)
    opts = parse_list(args.opts, OPTS)
    objdir = ROOT / "test" / ".ct_scan_objs"
    objdir.mkdir(parents=True, exist_ok=True)

    all_lines = [
        "# Lithium-Code CT Scan Matrix Report",
        "",
        f"date_utc: {_dt.datetime.now(_dt.timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}",
        "root: <artifact-root>",
        f"backends: {' '.join(backends)}",
        f"modes: {' '.join(map(str, modes))}",
        f"compilers: {' '.join(ccs)}",
        f"opts: {' '.join(opts)}",
        f"include_optional: {int(args.include_optional)}",
        "allowlist:",
        "  - gather: avx2/fips202x4.c and avx512/fips202x8.c public SHAKE lane pointers/offsets only",
        "  - public_div: fips202 one-shot SHAKE output-length block counts and sampler Gaussian initial refill block counts",
        "",
    ]
    ok = True
    for backend in backends:
        for variant in variants_for(backend, args.include_optional):
            for mode in modes:
                for cc in ccs:
                    for opt in opts:
                        one_ok, lines = compile_and_scan(backend, mode, cc, opt, variant, objdir)
                        ok = ok and one_ok
                        all_lines.extend(lines)
                        all_lines.append("")

    all_lines.append("OVERALL: " + ("PASS" if ok else "FAIL"))
    text = "\n".join(all_lines) + "\n"
    if args.report:
        report = ROOT / args.report
        report.parent.mkdir(parents=True, exist_ok=True)
        report.write_text(text)
    print(text)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
