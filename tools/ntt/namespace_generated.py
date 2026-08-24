#!/usr/bin/env python3
"""Namespace generated NTT assembly and constants for Lithium builds."""
from pathlib import Path

AVX2_SYMBOLS = [
    "ntt_avx", "invntt_tomont_avx", "basemul_avx", "pointwise_avx",
    "reduce_avx", "nttunpack_avx", "nttunpack16_avx",
    "poly_freeze_to_u16_avx",
    "poly_freeze_small_to_u16_avx", "poly_cast_modq_avx",
    "ntt_qdata", "ntt_ninv", "ntt_cross_fwd", "ntt_cross_inv",
    "ntt_zetas_fwd", "ntt_zetas_inv", "ntt_z0", "ntt_z0inv",
    "ntt_scale", "ntt_basemul_gamma", "ntt_qinv", "ntt_r2",
]

AVX512_SYMBOLS = [
    "ntt_avx512", "invntt_tomont_avx512", "basemul_avx512",
    "pointwise_avx512", "nttunpack_avx512", "nttunpack16_avx512",
    "poly_freeze_to_u16_avx512",
    "poly_freeze_small_to_u16_avx512", "poly_cast_modq_avx512",
    "ntt512_qdata", "ntt512_scale", "ntt512_qinv", "ntt512_r2",
    "ntt512_zetas_fwd", "ntt512_zetas_inv", "ntt512_basemul_gamma",
]
for k in (16, 8, 4, 2, 1):
    AVX512_SYMBOLS.extend([f"shuf{k}_c", f"shuf{k}_d"])


def _add_config_include(text):
    marker = '#include <stdint.h>\n'
    include = '#include "../config.h"\n'
    if include in text:
        return text
    if marker in text:
        return text.replace(marker, marker + include, 1)
    return include + text


def _namespace_c_symbol_defs(text, symbols):
    text = _add_config_include(text)
    for sym in sorted(symbols, key=len, reverse=True):
        text = text.replace(f"const int16_t {sym}",
                            f"const int16_t LITHIUM_NAMESPACE({sym})")
        text = text.replace(f"const uint16_t {sym}",
                            f"const uint16_t LITHIUM_NAMESPACE({sym})")
    return text


def _namespace_asm(text, symbols):
    include = '#include "../config.h"\n'
    if include not in text:
        # Keep the include before .text so the assembler sees only expanded labels.
        pos = text.find(".text")
        if pos >= 0:
            text = text[:pos] + include + text[pos:]
        else:
            text = include + text
    for sym in sorted(symbols, key=len, reverse=True):
        ns = f"LITHIUM_NAMESPACE({sym})"
        text = text.replace(f".global {sym}", f".global {ns}")
        text = text.replace(f"\n{sym}:", f"\n{ns}:")
        text = text.replace(f" {sym}(%rip)", f" {ns}(%rip)")
        text = text.replace(f"\t{sym}(%rip)", f"\t{ns}(%rip)")
    return text


def namespace_avx2(out_dir="."):
    d = Path(out_dir)
    p = d / "ntt_consts.c"
    if p.exists():
        p.write_text(_namespace_c_symbol_defs(p.read_text(), AVX2_SYMBOLS), encoding="ascii")
    p = d / "ntt.S"
    if p.exists():
        p.write_text(_namespace_asm(p.read_text(), AVX2_SYMBOLS), encoding="ascii")


def namespace_avx512(out_dir="."):
    d = Path(out_dir)
    p = d / "ntt_consts_avx512.c"
    if p.exists():
        p.write_text(_namespace_c_symbol_defs(p.read_text(), AVX512_SYMBOLS), encoding="ascii")
    p = d / "ntt_avx512.S"
    if p.exists():
        p.write_text(_namespace_asm(p.read_text(), AVX512_SYMBOLS), encoding="ascii")
