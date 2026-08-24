#!/usr/bin/env bash
#
# Methodology:
#   1. Each external code repository is pinned to a fixed commit hash, using the
#      version available as of 2026/06/05. The URLs and exact hashes are listed
#      in `third_party/upstream-commits.txt` and mirrored below where needed so
#      every reproduced benchmark is tied to an immutable source state.
#   2. Benchmark output follows the same style as the `make speed_txt` reports in
#      our Lithium implementation: each mode prints public-key/signature sizes
#      and their sum, e.g., `mode=120, pk=1056B, sig=1187B, pk+sig=2243B, mlen=59B`.
#      All signing benchmarks use 59-byte short messages, matching the
#      eBACS/SUPERCOP benchmarking convention.
#      Keygen and signing benchmarks also print detailed statistical metrics, e.g.,
#      `sign (10000 runs) in cycles: p50: 353950, p90: 360546, p99: 369714, avg: 354802`
#      and
#      `sign stability: stddev: 6348, cv: 1.79%, outlier_count(>p50*1.5): 1/10000, p90/p50: 1.02x, p99/p50: 1.04x`,
#      so the distribution, variance, and outliers can be analyzed rather than
#      relying on a single aggregate timing value.
#   3. After each benchmark binary is compiled, the script sleeps for 2 seconds
#      before running it, so short build bursts do not immediately precede
#      timing collection.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXTERNAL_DIR="$ROOT_DIR/external"
THIRD_PARTY_DIR="$ROOT_DIR/third_party"
TAIL_SRC_DIR="$ROOT_DIR/ref/test"

DILITHIUM_URL="https://github.com/pq-crystals/dilithium.git"
DILITHIUM_COMMIT="bba9534a3c5cbef66481cd6518d786e23ae601b5"
FALCON_REF_URL="https://github.com/pornin/c-fn-dsa.git"
FALCON_REF_COMMIT="3aa50485d5f992adbca8ce32a307401d1ed37dc8"
FALCON_OPT_URL="$(awk -F'[()]' '/commits.falcon_opt/{url=$2; sub(/,.*/, "", url); print url; exit}' "$THIRD_PARTY_DIR/upstream-commits.txt")"
FALCON_OPT_COMMIT="61b326e9bb6afe2ea4927e7a25832e7a64485367"
HAETAE_URL="https://github.com/kpqc-cryptocraft/KpqClean_ver2.git"
HAETAE_COMMIT="0896ebec2441804f13ce2fca9c942f32fc62805c"
HAETAE_TAG="v2.0.0"
NTRUPLUS_URL="https://github.com/KU-Cryptographic-Protocol-Lab/ntruplus_sign.git"
NTRUPLUS_COMMIT="4a540aa184214f9de1d4b130f5268da00b1468b4"
RACCOON_URL="https://github.com/masksign/raccoon.git"
RACCOON_COMMIT="e789b4b72a2b7e8a2205df49c487736985fc8417"
GALACTICS_URL="https://github.com/espitau/GALACTICS.git"
GALACTICS_COMMIT="79141a6b4032dbeb79fa8cab0d3c715c0c3eb442"
HAWK_URL="https://github.com/hawk-sign/dev.git"
HAWK_COMMIT="1b9fef52559273fe7b40fe3e22968eaedd3a4c2a"
DILITHIUM_XU_URL="https://github.com/xurq5953/Optimizing-Dilithium-Implementation-with-AVX2-512.git"
DILITHIUM_XU_COMMIT="100eb09973bce5e41fd14373aa279ddca06ef1d2"
BENCH_SIGN_RUNS="${REPRO_BENCH_SIGN_RUNS:-10000}"


usage() {
    cat <<'USAGE'
Usage: ./reproduce.sh [--all] [--dilithium] [--falcon] [--haetae] [--ntruplus] [--raccoon] [--bliss] [--hawk] [--Dilithium-Xu]

Options:
  --all        Reproduce all supported third-party benchmarks in order.
  --dilithium  Reproduce only Dilithium ref/ and avx2/ results.
  --falcon     Reproduce only Falcon-Ref and Falcon-Opt results.
  --haetae     Reproduce only HAETAE clean/ and avx2/ results.
  --ntruplus   Reproduce only NTRU+Sign ref results.
  --raccoon    Reproduce only Raccoon ref results.
  --bliss      Reproduce only BLISS-I from GALACTICS ref results.
  --hawk       Reproduce only HAWK-512/1024 ref and AVX2 results.
  --Dilithium-Xu
               Reproduce Xu et al. Dilithium AVX2/AVX-512 results.
USAGE
}

log() {
    printf '[reproduce] %s\n' "$*"
}

check_required_tools() {
    local missing=()
    local tool
    for tool in git gcc clang make awk cp date uname grep head tr find python3; do
        if ! command -v "$tool" >/dev/null 2>&1; then
            missing+=("$tool")
        fi
    done

    if [ "${#missing[@]}" -ne 0 ]; then
        printf 'ERROR: missing required tool(s): %s\n' "${missing[*]}" >&2
        printf 'Please install the missing tool(s) and rerun this script.\n' >&2
        exit 1
    fi
}

run_logged() {
    local logfile="$1"
    shift
    printf '$' >> "$logfile"
    printf ' %q' "$@" >> "$logfile"
    printf '\n' >> "$logfile"
    "$@" >> "$logfile" 2>&1
}

append_run_with_blank() {
    local outfile="$1"
    shift
    sleep 2
    "$@" >> "$outfile" 2>&1
    printf '\n' >> "$outfile"
}

cpu_model() {
    if command -v lscpu >/dev/null 2>&1; then
        lscpu | awk -F: '/^Model name:/ { sub(/^[[:space:]]+/, "", $2); print $2; exit }'
    elif [ -r /proc/cpuinfo ]; then
        awk -F: '/model name/ { sub(/^[[:space:]]+/, "", $2); print $2; exit }' /proc/cpuinfo
    else
        uname -m
    fi
}

compiler_line() {
    local cc="${CC:-cc}"
    if command -v "$cc" >/dev/null 2>&1; then
        "$cc" --version | head -1
    elif command -v gcc >/dev/null 2>&1; then
        gcc --version | head -1
    else
        printf 'unknown\n'
    fi
}

has_cpu_flag() {
    # Match a CPU flag ignoring '_' differences between the kernel's
    # /proc/cpuinfo / lscpu spelling (e.g. avx512_vbmi2) and the compiler's
    # intrinsic name (e.g. avx512vbmi2).
    local flag
    flag="$(printf '%s' "$1" | tr -d '_')"
    if command -v lscpu >/dev/null 2>&1; then
        lscpu | awk -F: '/^Flags:/ { print $2 }' | tr ' ' '\n' | tr -d '_' | grep -qx "$flag"
    elif [ -r /proc/cpuinfo ]; then
        awk -F: '/flags/ { print $2; exit }' /proc/cpuinfo | tr ' ' '\n' | tr -d '_' | grep -qx "$flag"
    else
        return 1
    fi
}

write_machine_info() {
    local outfile="$1"
    {
        printf 'date    : '; date -u '+%Y-%m-%dT%H:%M:%SZ'
        printf 'kernel  : '; uname -sr
        printf 'cpu     : '; cpu_model
        printf 'cc      : '; compiler_line
        printf 'cflags  : %s\n' "${CFLAGS:-upstream Makefile defaults}"
        printf '\n'
    } > "$outfile"
}

append_commit_info() {
    local outfile="$1"
    {
        printf 'commits.dilithium = %s (%s, HEAD of master on 2026-06-03, no exact tag)\n' "$DILITHIUM_COMMIT" "$DILITHIUM_URL"
        printf 'commits.falcon_ref = %s (%s, HEAD of main on 2026-06-03, no exact tag)\n' "$FALCON_REF_COMMIT" "$FALCON_REF_URL"
        printf 'commits.falcon_opt = %s (%s, HEAD of main on 2026-06-03, no exact tag)\n' "$FALCON_OPT_COMMIT" "$FALCON_OPT_URL"
        printf 'commits.haetae = %s (%s, exact tag %s)\n' "$HAETAE_COMMIT" "$HAETAE_URL" "$HAETAE_TAG"
        printf 'commits.ntruplus = %s (%s, HEAD of main on 2026-06-05, no exact tag)\n' "$NTRUPLUS_COMMIT" "$NTRUPLUS_URL"
        printf 'commits.raccoon = %s (%s, HEAD of main on 2026-06-05, no exact tag)\n' "$RACCOON_COMMIT" "$RACCOON_URL"
        printf 'commits.galactics = %s (%s, HEAD of master on 2026-06-05, no exact tag)\n' "$GALACTICS_COMMIT" "$GALACTICS_URL"
        printf 'commits.hawk = %s (%s, HEAD of main on 2026-06-06, no exact tag)\n' "$HAWK_COMMIT" "$HAWK_URL"
        printf 'commits.dilithium_xu = %s (%s, HEAD of main on 2026-06-06, no exact tag)\n' "$DILITHIUM_XU_COMMIT" "$DILITHIUM_XU_URL"
        printf '\n'
    } >> "$outfile"
}


append_stability_info() {
    local outfile="$1"
    {
        printf 'stability.p25_p50_p75_p90_p99 = nearest-rank percentiles of keygen/sign cycles\n'
        printf 'stability.StQ1_StQ2_StQ3 = stabilized quartiles: [1/8,3/8], [3/8,5/8], [5/8,7/8] means; see Daniel J. Bernstein, Robust statistics for rejection-sampling timings, https://cr.yp.to/papers/rsrst-20250727.pdf\n'
        printf 'stability.stddev = population standard deviation of keygen/sign cycles; lower means less run-to-run spread\n'
        printf 'stability.cv = stddev/avg; rough guide: <1%% very stable, 1-5%% stable, 5-15%% noisy, >15%% unstable or long-tailed\n'
        printf 'stability.outlier_count = number of keygen/sign runs slower than 1.5x p50\n'
        printf 'stability.p90_p50 = p90/p50 tail-inflation ratio; values near 1.00x mean a short tail\n'
        printf 'stability.p99_p50 = p99/p50 tail-inflation ratio; values near 1.00x mean a short tail\n'
        printf '\n'
    } >> "$outfile"
}

write_commit_manifest() {
    local outfile="$THIRD_PARTY_DIR/upstream-commits.txt"
    write_machine_info "$outfile"
    append_commit_info "$outfile"
}

setup_repo() {
    local name="$1"
    local url="$2"
    local commit="$3"
    local dir="$EXTERNAL_DIR/$name"

    mkdir -p "$EXTERNAL_DIR"
    if [ ! -d "$dir/.git" ]; then
        log "cloning $name"
        git clone "$url" "$dir"
    fi
    log "checking out $name@$commit"
    git -C "$dir" fetch --tags origin
    git -C "$dir" checkout --detach "$commit"
    git -C "$dir" reset --hard "$commit"
    git -C "$dir" clean -fdx
}

copy_tailbench_test_dir() {
    local test_dir="$1"
    cp "$TAIL_SRC_DIR/cpucycles.h" "$test_dir/cpucycles.h"
    cp "$TAIL_SRC_DIR/cpucycles.c" "$test_dir/cpucycles.c"
}

copy_tailbench_named_dir() {
    local bench_dir="$1"
    mkdir -p "$bench_dir"
    cp "$TAIL_SRC_DIR/cpucycles.h" "$bench_dir/cpucycles.h"
    cp "$TAIL_SRC_DIR/cpucycles.c" "$bench_dir/cpucycles.c"
}

allow_bench_sign_runs_override() {
    local header="$1"
    python3 - "$header" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()
s = s.replace("#define BENCH_SIGN_RUNS 10000",
              "#ifndef BENCH_SIGN_RUNS\n#define BENCH_SIGN_RUNS 10000\n#endif")
p.write_text(s)
PY
}

set_bench_sign_runs_default() {
    local header="$1"
    local runs="$2"
    python3 - "$header" "$runs" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
runs = sys.argv[2]
s = p.read_text()
s = s.replace("#define BENCH_SIGN_RUNS 10000",
              f"#define BENCH_SIGN_RUNS {runs}")
p.write_text(s)
PY
}

write_speed_nist_c() {
    local dst="$1"
    cat > "$dst" <<'EOF'
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "api.h"
#if __has_include("crypto_sign.h")
#include "crypto_sign.h"
#endif
#include "cpucycles.h"

#define VERIFY_RUNS 10000
#define MLEN 59
#define FEEDBACK_MLEN (MLEN + CRYPTO_BYTES)

static uint64_t verify_cycles[VERIFY_RUNS];

typedef struct keygen_bench_ctx {
    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char sk[CRYPTO_SECRETKEYBYTES];
} keygen_bench_ctx;

typedef struct sign_bench_ctx {
    unsigned char sk[CRYPTO_SECRETKEYBYTES];
    unsigned char sm[FEEDBACK_MLEN + CRYPTO_BYTES];
    unsigned long long smlen;
} sign_bench_ctx;

static int keygen_once(void *vctx)
{
    keygen_bench_ctx *ctx = (keygen_bench_ctx *)vctx;

    return crypto_sign_keypair(ctx->pk, ctx->sk);
}

static int sign_once(void *vctx)
{
    sign_bench_ctx *ctx = (sign_bench_ctx *)vctx;
    /* Feed each signed message back as the next message. */
    return crypto_sign(ctx->sm, &ctx->smlen, ctx->sm, MLEN, ctx->sk);
}

int main(void)
{
    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char sk[CRYPTO_SECRETKEYBYTES];
    unsigned char msg[MLEN] = {0};
    unsigned char opened[MLEN + CRYPTO_BYTES];
    unsigned char sm[MLEN + CRYPTO_BYTES];
    unsigned long long opened_len = 0;
    unsigned long long smlen = 0;
    keygen_bench_ctx keygen_ctx;
    sign_bench_ctx sign_ctx;
    tail_sign_result keygen_result;
    tail_sign_result sign_result;

    if (bench_tail_sign("keygen", &keygen_result, keygen_once,
                        &keygen_ctx) != 0)
        return 1;

    memcpy(pk, keygen_ctx.pk, sizeof pk);
    memcpy(sk, keygen_ctx.sk, sizeof sk);
    memset(&sign_ctx, 0, sizeof sign_ctx);
    memcpy(sign_ctx.sk, sk, sizeof sign_ctx.sk);

    if (bench_tail_sign("sign", &sign_result, sign_once, &sign_ctx) != 0)
        return 1;

    if (crypto_sign(sm, &smlen, msg, sizeof msg, sk) != 0)
        return 1;

    for (int i = 0; i < VERIFY_RUNS; i++) {
        uint64_t begin = cpucycles();
        if (crypto_sign_open(opened, &opened_len, sm, smlen, pk) != 0)
            return 1;
        verify_cycles[i] = cpucycles() - begin;
    }

    printf("mode=%s, pk=%dB, sig=%dB, pk+sig=%dB, mlen=%dB\n", CRYPTO_ALGNAME,
           CRYPTO_PUBLICKEYBYTES, CRYPTO_BYTES,
           CRYPTO_PUBLICKEYBYTES + CRYPTO_BYTES, MLEN);
    bench_print_tail(&keygen_result);
    bench_print_tail(&sign_result);
    printf("verify median (%d runs): %llu cycles\n", VERIFY_RUNS,
           (unsigned long long)cpucycles_median(verify_cycles, VERIFY_RUNS));
    return 0;
}
EOF
}

write_speed_ntruplus_c() {
    local dst="$1"
    cat > "$dst" <<'EOF'
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "api.h"
#include "sign.h"
#include "cpucycles.h"

#define crypto_sign_signature crypto_sign

#define VERIFY_RUNS 10000
#define MLEN 59

static uint64_t verify_cycles[VERIFY_RUNS];

typedef struct keygen_bench_ctx {
    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char sk[CRYPTO_SECRETKEYBYTES];
} keygen_bench_ctx;

typedef struct sign_bench_ctx {
    unsigned char sk[CRYPTO_SECRETKEYBYTES];
    unsigned char sig[CRYPTO_BYTES];
    unsigned long long siglen;
} sign_bench_ctx;

static int keygen_once(void *vctx)
{
    keygen_bench_ctx *ctx = (keygen_bench_ctx *)vctx;

    return crypto_sign_keypair(ctx->pk, ctx->sk);
}

static int sign_once(void *vctx)
{
    sign_bench_ctx *ctx = (sign_bench_ctx *)vctx;
    /* Feed each detached signature back as the next message. */
    return crypto_sign_signature(ctx->sig, &ctx->siglen, ctx->sig, MLEN, ctx->sk);
}

int main(void)
{
    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char sk[CRYPTO_SECRETKEYBYTES];
    unsigned char msg[MLEN] = {0};
    unsigned char sig[CRYPTO_BYTES];
    unsigned long long siglen = 0;
    keygen_bench_ctx keygen_ctx;
    sign_bench_ctx sign_ctx;
    tail_sign_result keygen_result;
    tail_sign_result sign_result;

    if (bench_tail_sign("keygen", &keygen_result, keygen_once,
                        &keygen_ctx) != 0)
        return 1;

    memcpy(pk, keygen_ctx.pk, sizeof pk);
    memcpy(sk, keygen_ctx.sk, sizeof sk);
    memset(&sign_ctx, 0, sizeof sign_ctx);
    memcpy(sign_ctx.sk, sk, sizeof sign_ctx.sk);

    if (bench_tail_sign("sign", &sign_result, sign_once, &sign_ctx) != 0)
        return 1;

    if (crypto_sign(sig, &siglen, msg, sizeof msg, sk) != 0)
        return 1;

    for (int i = 0; i < VERIFY_RUNS; i++) {
        uint64_t begin = cpucycles();
        if (crypto_sign_verify(msg, sizeof msg, sig, siglen, pk) != 0)
            return 1;
        verify_cycles[i] = cpucycles() - begin;
    }

    printf("mode=%s, pk=%dB, sig=%dB, pk+sig=%dB, mlen=%dB\n", CRYPTO_ALGNAME,
           CRYPTO_PUBLICKEYBYTES, CRYPTO_BYTES,
           CRYPTO_PUBLICKEYBYTES + CRYPTO_BYTES, MLEN);
    bench_print_tail(&keygen_result);
    bench_print_tail(&sign_result);
    printf("verify median (%d runs): %llu cycles\n", VERIFY_RUNS,
           (unsigned long long)cpucycles_median(verify_cycles, VERIFY_RUNS));
    return 0;
}
EOF
}

write_speed_haetae_c() {
    local dst="$1"
    cat > "$dst" <<'EOF'
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "api.h"
#include "sign.h"
#include "cpucycles.h"

#define VERIFY_RUNS 10000
#define MLEN 59

static uint64_t verify_cycles[VERIFY_RUNS];

typedef struct keygen_bench_ctx {
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
} keygen_bench_ctx;

typedef struct sign_bench_ctx {
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t sig[CRYPTO_BYTES];
    size_t siglen;
} sign_bench_ctx;

static int keygen_once(void *vctx)
{
    keygen_bench_ctx *ctx = (keygen_bench_ctx *)vctx;

    return crypto_sign_keypair(ctx->pk, ctx->sk);
}

static int sign_once(void *vctx)
{
    sign_bench_ctx *ctx = (sign_bench_ctx *)vctx;
    /* Feed each detached signature back as the next message. */
    return crypto_sign_signature(ctx->sig, &ctx->siglen, ctx->sig, MLEN, ctx->sk);
}

int main(void)
{
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t msg[MLEN] = {0};
    uint8_t sig[CRYPTO_BYTES];
    size_t siglen = 0;
    keygen_bench_ctx keygen_ctx;
    sign_bench_ctx sign_ctx;
    tail_sign_result keygen_result;
    tail_sign_result sign_result;

    if (bench_tail_sign("keygen", &keygen_result, keygen_once,
                        &keygen_ctx) != 0)
        return 1;

    memcpy(pk, keygen_ctx.pk, sizeof pk);
    memcpy(sk, keygen_ctx.sk, sizeof sk);
    memset(&sign_ctx, 0, sizeof sign_ctx);
    memcpy(sign_ctx.sk, sk, sizeof sign_ctx.sk);

    if (bench_tail_sign("sign", &sign_result, sign_once, &sign_ctx) != 0)
        return 1;

    if (crypto_sign_signature(sig, &siglen, msg, sizeof msg, sk) != 0)
        return 1;

    for (int i = 0; i < VERIFY_RUNS; i++) {
        uint64_t begin = cpucycles();
        if (crypto_sign_verify(sig, siglen, msg, sizeof msg, pk) != 0)
            return 1;
        verify_cycles[i] = cpucycles() - begin;
    }

    printf("mode=%s, pk=%dB, sig=%dB, pk+sig=%dB, mlen=%dB\n", CRYPTO_ALGNAME,
           CRYPTO_PUBLICKEYBYTES, CRYPTO_BYTES,
           CRYPTO_PUBLICKEYBYTES + CRYPTO_BYTES, MLEN);
    bench_print_tail(&keygen_result);
    bench_print_tail(&sign_result);
    printf("verify median (%d runs): %llu cycles\n", VERIFY_RUNS,
           (unsigned long long)cpucycles_median(verify_cycles, VERIFY_RUNS));
    return 0;
}
EOF
}



write_speed_dilithium_xu_c() {
    local dst="$1"
    local mode_name="$2"
    local impl_name="$3"
    local pk_macro="$4"
    local sk_macro="$5"
    local sig_macro="$6"
    local keypair_fn="$7"
    local sign_fn="$8"
    local verify_fn="$9"
    cat > "$dst" <<EOF
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "../sign.h"
#include "cpucycles.h"

#define MODE_NAME "$mode_name"
#define IMPL_NAME "$impl_name"
#define CRYPTO_PUBLICKEYBYTES $pk_macro
#define CRYPTO_SECRETKEYBYTES $sk_macro
#define CRYPTO_BYTES $sig_macro
#define KEYPAIR $keypair_fn
#define SIGN_SIGNATURE $sign_fn
#define VERIFY_SIGNATURE $verify_fn

#define VERIFY_RUNS 10000
#define MLEN 59

static uint64_t verify_cycles[VERIFY_RUNS];
static volatile int verify_acc;

typedef struct keygen_bench_ctx {
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
} keygen_bench_ctx;

typedef struct sign_bench_ctx {
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t sig[CRYPTO_BYTES];
    size_t siglen;
} sign_bench_ctx;

static int keygen_once(void *vctx)
{
    keygen_bench_ctx *ctx = (keygen_bench_ctx *)vctx;

    return KEYPAIR(ctx->pk, ctx->sk);
}

static int sign_once(void *vctx)
{
    sign_bench_ctx *ctx = (sign_bench_ctx *)vctx;
    /* Feed each detached signature back as the next message. */
    return SIGN_SIGNATURE(ctx->sig, &ctx->siglen, ctx->sig, MLEN, ctx->sk);
}

int main(void)
{
    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];
    uint8_t msg[MLEN] = {0};
    uint8_t sig[CRYPTO_BYTES];
    size_t siglen = 0;
    keygen_bench_ctx keygen_ctx;
    sign_bench_ctx sign_ctx;
    tail_sign_result keygen_result;
    tail_sign_result sign_result;

    if (bench_tail_sign("keygen", &keygen_result, keygen_once,
                        &keygen_ctx) != 0)
        return 1;

    memcpy(pk, keygen_ctx.pk, sizeof pk);
    memcpy(sk, keygen_ctx.sk, sizeof sk);
    memset(&sign_ctx, 0, sizeof sign_ctx);
    memcpy(sign_ctx.sk, sk, sizeof sign_ctx.sk);

    if (bench_tail_sign("sign", &sign_result, sign_once, &sign_ctx) != 0)
        return 1;

    if (SIGN_SIGNATURE(sig, &siglen, msg, sizeof msg, sk) != 0)
        return 1;

    for (int i = 0; i < VERIFY_RUNS; i++) {
        uint64_t begin = cpucycles();
        verify_acc |= VERIFY_SIGNATURE(sig, siglen, msg, sizeof msg, pk);
        verify_cycles[i] = cpucycles() - begin;
    }

    printf("mode=%s, impl=%s, pk=%dB, sig=%dB, pk+sig=%dB, mlen=%dB\n",
           MODE_NAME, IMPL_NAME, CRYPTO_PUBLICKEYBYTES, CRYPTO_BYTES,
           CRYPTO_PUBLICKEYBYTES + CRYPTO_BYTES, MLEN);
    bench_print_tail(&keygen_result);
    bench_print_tail(&sign_result);
    printf("verify median (%d runs): %llu cycles\n", VERIFY_RUNS,
           (unsigned long long)cpucycles_median(verify_cycles, VERIFY_RUNS));
    return 0;
}
EOF
}

write_speed_hawk_c() {
    local dst="$1"
    cat > "$dst" <<'EOF'
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "api.h"
#include "rng.h"
#include "cpucycles.h"

#ifndef SECURITY_BITS
#define SECURITY_BITS 0
#endif
#ifndef HAWK_N
#define HAWK_N 0
#endif

#define VERIFY_RUNS 10000
#define MLEN 59
#define FEEDBACK_MLEN (MLEN + CRYPTO_BYTES)

static uint64_t verify_cycles[VERIFY_RUNS];

typedef struct keygen_bench_ctx {
    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char sk[CRYPTO_SECRETKEYBYTES];
} keygen_bench_ctx;

typedef struct sign_bench_ctx {
    unsigned char sk[CRYPTO_SECRETKEYBYTES];
    unsigned char sm[FEEDBACK_MLEN + CRYPTO_BYTES];
    unsigned long long smlen;
} sign_bench_ctx;

static int keygen_once(void *vctx)
{
    keygen_bench_ctx *ctx = (keygen_bench_ctx *)vctx;

    return crypto_sign_keypair(ctx->pk, ctx->sk);
}

static int sign_once(void *vctx)
{
    sign_bench_ctx *ctx = (sign_bench_ctx *)vctx;
    /* Feed each signed message back as the next message. */
    return crypto_sign(ctx->sm, &ctx->smlen, ctx->sm, MLEN, ctx->sk);
}

int main(void)
{
    unsigned char entropy[48] = {0};
    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char sk[CRYPTO_SECRETKEYBYTES];
    unsigned char msg[MLEN] = {0};
    unsigned char opened[MLEN + CRYPTO_BYTES];
    unsigned char sm[MLEN + CRYPTO_BYTES];
    unsigned long long opened_len = 0;
    unsigned long long smlen = 0;
    keygen_bench_ctx keygen_ctx;
    sign_bench_ctx sign_ctx;
    tail_sign_result keygen_result;
    tail_sign_result sign_result;

    randombytes_init(entropy, NULL, SECURITY_BITS);

    if (bench_tail_sign("keygen", &keygen_result, keygen_once,
                        &keygen_ctx) != 0)
        return 1;

    memcpy(pk, keygen_ctx.pk, sizeof pk);
    memcpy(sk, keygen_ctx.sk, sizeof sk);
    memset(&sign_ctx, 0, sizeof sign_ctx);
    memcpy(sign_ctx.sk, sk, sizeof sign_ctx.sk);

    if (bench_tail_sign("sign", &sign_result, sign_once, &sign_ctx) != 0)
        return 1;

    if (crypto_sign(sm, &smlen, msg, sizeof msg, sk) != 0)
        return 1;

    for (int i = 0; i < VERIFY_RUNS; i++) {
        uint64_t begin = cpucycles();
        if (crypto_sign_open(opened, &opened_len, sm, smlen, pk) != 0)
            return 1;
        verify_cycles[i] = cpucycles() - begin;
    }

    printf("mode=%s, n=%d, security=%d, pk=%dB, sig=%dB, pk+sig=%dB, mlen=%dB\n",
           CRYPTO_ALGNAME, HAWK_N, SECURITY_BITS, CRYPTO_PUBLICKEYBYTES,
           CRYPTO_BYTES, CRYPTO_PUBLICKEYBYTES + CRYPTO_BYTES, MLEN);
    bench_print_tail(&keygen_result);
    bench_print_tail(&sign_result);
    printf("verify median (%d runs): %llu cycles\n", VERIFY_RUNS,
           (unsigned long long)cpucycles_median(verify_cycles, VERIFY_RUNS));
    return 0;
}
EOF
}

append_skip() {
    local outfile="$1"
    shift
    printf '%s\n\n' "$*" >> "$outfile"
}

run_dilithium() {
    local repo="$EXTERNAL_DIR/dilithium"
    local outfile="$THIRD_PARTY_DIR/dilithium.txt"
    local buildlog="$THIRD_PARTY_DIR/dilithium.build.log"

    mkdir -p "$THIRD_PARTY_DIR"
    : > "$buildlog"
    write_machine_info "$outfile"
    append_commit_info "$outfile"
    append_stability_info "$outfile"
    printf '=== Dilithium Performance ===\n\n' >> "$outfile"

    setup_repo "dilithium" "$DILITHIUM_URL" "$DILITHIUM_COMMIT"
    git -C "$repo" apply "$THIRD_PARTY_DIR/dilithium-tailbench.patch"
    copy_tailbench_test_dir "$repo/ref/test"
    copy_tailbench_test_dir "$repo/avx2/test"

    printf -- '-- ref/ --\n' >> "$outfile"
    run_logged "$buildlog" make -j1 -s -C "$repo/ref" clean
    run_logged "$buildlog" make -j1 -s -C "$repo/ref" speed
    append_run_with_blank "$outfile" "$repo/ref/test/test_speed2"
    append_run_with_blank "$outfile" "$repo/ref/test/test_speed3"
    append_run_with_blank "$outfile" "$repo/ref/test/test_speed5"

    printf -- '-- avx2/ --\n' >> "$outfile"
    if has_cpu_flag avx2; then
        run_logged "$buildlog" make -j1 -s -C "$repo/avx2" clean
        run_logged "$buildlog" make -j1 -s -C "$repo/avx2" speed
        append_run_with_blank "$outfile" "$repo/avx2/test/test_speed2"
        append_run_with_blank "$outfile" "$repo/avx2/test/test_speed3"
        append_run_with_blank "$outfile" "$repo/avx2/test/test_speed5"
    else
        append_skip "$outfile" 'SKIP: CPU does not report avx2 support.'
    fi
    log "wrote $outfile"
}

run_falcon_ref() {
    local outfile="$1"
    local buildlog="$THIRD_PARTY_DIR/falcon-ref.build.log"
    local repo="$EXTERNAL_DIR/falcon-ref"

    : > "$buildlog"
    setup_repo "falcon-ref" "$FALCON_REF_URL" "$FALCON_REF_COMMIT"
    git -C "$repo" apply "$THIRD_PARTY_DIR/falcon-ref-tailbench.patch"
    copy_tailbench_named_dir "$repo/irs_tailbench"

    printf '=== Falcon-Ref Performance ===\n\n' >> "$outfile"
    run_logged "$buildlog" make -j1 -s -C "$repo" clean
    run_logged "$buildlog" make -j1 -s -C "$repo" CFLAGS="-W -Wextra -Wundef -Wshadow -O2" speed_fndsa
    append_run_with_blank "$outfile" "$repo/speed_fndsa"
}

run_falcon_opt() {
    local outfile="$1"
    local buildlog="$THIRD_PARTY_DIR/falcon-opt.build.log"
    local repo="$EXTERNAL_DIR/falcon-opt"
    local opt="$repo/opt"

    : > "$buildlog"
    setup_repo "falcon-opt" "$FALCON_OPT_URL" "$FALCON_OPT_COMMIT"
    git -C "$repo" apply "$THIRD_PARTY_DIR/falcon-opt-tailbench.patch"
    copy_tailbench_named_dir "$opt/irs_tailbench"

    printf '=== Falcon-Opt Performance ===\n\n' >> "$outfile"
    run_logged "$buildlog" make -j1 -s -C "$opt" clean

    printf -- '-- opt/ speed_fndsa_avx2 --\n' >> "$outfile"
    if has_cpu_flag avx2; then
        run_logged "$buildlog" make -j1 -s -C "$opt" CFLAGS_COMM="-W -Wextra -Wundef -Wshadow -O2 \$(PG) \$(X4) -DWARMUP_N=2000 -DTESTS_N=10000 -g" out/speed_fndsa_avx2
        append_run_with_blank "$outfile" "$opt/out/speed_fndsa_avx2"
    else
        append_skip "$outfile" 'SKIP: CPU does not report avx2 support.'
    fi

    printf -- '-- opt/ speed_fndsa_avx512f --\n' >> "$outfile"
    if has_cpu_flag avx512f; then
        run_logged "$buildlog" make -j1 -s -C "$opt" CFLAGS_COMM="-W -Wextra -Wundef -Wshadow -O2 \$(PG) \$(X4) -DWARMUP_N=2000 -DTESTS_N=10000 -g" out/speed_fndsa_avx512f
        append_run_with_blank "$outfile" "$opt/out/speed_fndsa_avx512f"
    else
        append_skip "$outfile" 'SKIP: CPU does not report avx512f support.'
    fi
}

run_falcon() {
    local outfile="$THIRD_PARTY_DIR/falcon.txt"

    mkdir -p "$THIRD_PARTY_DIR"
    write_machine_info "$outfile"
    append_commit_info "$outfile"
    append_stability_info "$outfile"
    run_falcon_ref "$outfile"
    run_falcon_opt "$outfile"
    log "wrote $outfile"
}

run_haetae_impl() {
    local label="$1"
    local dir="$2"
    local impl="$3"
    local outfile="$4"
    local buildlog="$5"
    local bench_dir="$dir/irs_tailbench"
    local common_dir="../../../common"

    printf -- '-- %s --\n' "$label" >> "$outfile"
    copy_tailbench_named_dir "$bench_dir"
    allow_bench_sign_runs_override "$bench_dir/cpucycles.h"
    write_speed_haetae_c "$bench_dir/speed_tailbench.c"
    run_logged "$buildlog" make -j1 -s -C "$dir" clean
    if [ "$impl" = "clean" ]; then
        run_logged "$buildlog" bash -c "cd '$dir' && \${CC:-cc} -Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls -Wshadow -Wpointer-arith -O3 -fomit-frame-pointer -I$common_dir -I. -Iirs_tailbench -DBENCH_SIGN_RUNS=$BENCH_SIGN_RUNS -o irs_tailbench/speed_tailbench irs_tailbench/speed_tailbench.c irs_tailbench/cpucycles.c decompose.c encoding.c fft.c fixpoint.c ntt.c packing.c poly.c polyfix.c polymat.c polyvec.c reduce.c sampler.c sign.c symmetric-shake.c $common_dir/randombytes.c $common_dir/fips202.c -lm"
    else
        run_logged "$buildlog" bash -c "cd '$dir' && \${CC:-cc} -Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls -Wshadow -Wpointer-arith -mavx2 -mbmi2 -mpopcnt -maes -march=native -mtune=native -O3 -fomit-frame-pointer -I$common_dir -I. -Iirs_tailbench -I$common_dir/keccak4x -I$common_dir/Keccak_avx2 -DBENCH_SIGN_RUNS=$BENCH_SIGN_RUNS -o irs_tailbench/speed_tailbench irs_tailbench/speed_tailbench.c irs_tailbench/cpucycles.c consts.c decompose.c encoding.c fft.c fixpoint.c packing.c poly.c polyfix.c polymat.c polyvec.c reduce.c sampler.c samplerx4.c sign.c symmetric-shake.c f1600x4.S invntt.S ntt.S pointwise.S shuffle.S keccakx4_vec.c $common_dir/randombytes.c $common_dir/aes256ctr.c $common_dir/keccak4x/fips202x4.c $common_dir/keccak4x/KeccakP-1600-times4-SIMD256.o $common_dir/Keccak_avx2/fips202.c $common_dir/Keccak_avx2/KeccakP-1600-AVX2.o -lcrypto -lm"
    fi
    append_run_with_blank "$outfile" "$bench_dir/speed_tailbench"
}

run_haetae() {
    local repo="$EXTERNAL_DIR/haetae-kpqclean"
    local outfile="$THIRD_PARTY_DIR/haetae.txt"
    local buildlog="$THIRD_PARTY_DIR/haetae.build.log"

    mkdir -p "$THIRD_PARTY_DIR"
    : > "$buildlog"
    write_machine_info "$outfile"
    append_commit_info "$outfile"
    append_stability_info "$outfile"
    printf '=== HAETAE Performance ===\n\n' >> "$outfile"

    setup_repo "haetae-kpqclean" "$HAETAE_URL" "$HAETAE_COMMIT"

    run_haetae_impl "clean/ HAETAE2" "$repo/crypto_sign/HAETAE2/clean" "clean" "$outfile" "$buildlog"
    run_haetae_impl "clean/ HAETAE3" "$repo/crypto_sign/HAETAE3/clean" "clean" "$outfile" "$buildlog"
    run_haetae_impl "clean/ HAETAE5" "$repo/crypto_sign/HAETAE5/clean" "clean" "$outfile" "$buildlog"
    if has_cpu_flag avx2; then
        run_haetae_impl "avx2/ HAETAE2" "$repo/crypto_sign/HAETAE2/avx2" "avx2" "$outfile" "$buildlog"
        run_haetae_impl "avx2/ HAETAE3" "$repo/crypto_sign/HAETAE3/avx2" "avx2" "$outfile" "$buildlog"
        run_haetae_impl "avx2/ HAETAE5" "$repo/crypto_sign/HAETAE5/avx2" "avx2" "$outfile" "$buildlog"
    else
        append_skip "$outfile" '-- avx2/ --'
        append_skip "$outfile" 'SKIP: CPU does not report avx2 support.'
    fi
    log "wrote $outfile"
}

run_ntruplus_impl() {
    local label="$1"
    local dir="$2"
    local outfile="$3"
    local buildlog="$4"
    local bench_dir="$dir/irs_tailbench"

    printf -- '-- %s --\n' "$label" >> "$outfile"
    copy_tailbench_named_dir "$bench_dir"
    allow_bench_sign_runs_override "$bench_dir/cpucycles.h"
    write_speed_ntruplus_c "$bench_dir/speed_tailbench.c"
    run_logged "$buildlog" make -j1 -s -C "$dir" clean
    run_logged "$buildlog" bash -c "cd '$dir' && \${CC:-cc} -Wall -Wextra -O3 -fomit-frame-pointer -I. -Iirs_tailbench -DBENCH_SIGN_RUNS=$BENCH_SIGN_RUNS -o irs_tailbench/speed_tailbench irs_tailbench/speed_tailbench.c irs_tailbench/cpucycles.c poly.c ntt.c fft.c fips202.c packing.c sampler.c sign.c symmetric-shake.c randombytes.c decompose.c gaussian.c rejection.c encoding.c -lcrypto -lssl -lm"
    append_run_with_blank "$outfile" "$bench_dir/speed_tailbench"
}

run_ntruplus() {
    local repo="$EXTERNAL_DIR/ntruplus_sign"
    local outfile="$THIRD_PARTY_DIR/ntruplus.txt"
    local buildlog="$THIRD_PARTY_DIR/ntruplus.build.log"

    mkdir -p "$THIRD_PARTY_DIR"
    : > "$buildlog"
    write_machine_info "$outfile"
    append_commit_info "$outfile"
    append_stability_info "$outfile"
    printf '=== NTRU+Sign Performance ===\n\n' >> "$outfile"

    setup_repo "ntruplus_sign" "$NTRUPLUS_URL" "$NTRUPLUS_COMMIT"
    run_ntruplus_impl "ref/ mode-512" "$repo/Reference_Implementation/NTRU+Sign512" "$outfile" "$buildlog"
    run_ntruplus_impl "ref/ mode-1024" "$repo/Reference_Implementation/NTRU+Sign1024" "$outfile" "$buildlog"
    log "wrote $outfile"
}

run_raccoon_impl() {
    local mode="$1"
    local dir="$2"
    local outfile="$3"
    local buildlog="$4"
    local bench_dir="$dir/irs_tailbench"

    printf -- '-- ref/ %s --\n' "$mode" >> "$outfile"
    printf '#define %s\n' "$mode" > "$dir/param_select.h"
    copy_tailbench_named_dir "$bench_dir"
    allow_bench_sign_runs_override "$bench_dir/cpucycles.h"
    write_speed_nist_c "$bench_dir/speed_tailbench.c"
    run_logged "$buildlog" make -j1 -s -C "$dir" clean
    run_logged "$buildlog" bash -c "cd '$dir' && \${CC:-cc} -Wall -Wextra -Ofast -march=native -Iinc -I. -Iirs_tailbench -DBENCH_SIGN_RUNS=$BENCH_SIGN_RUNS -o irs_tailbench/speed_tailbench irs_tailbench/speed_tailbench.c irs_tailbench/cpucycles.c racc_api.c racc_core.c racc_serial.c polyr.c ntt32.c ntt64.c xof_sample.c mask_random.c util/ct_util.c util/keccakf1600.c util/nist_random.c util/sha3_t.c util/test_aes1kt.c"
    append_run_with_blank "$outfile" "$bench_dir/speed_tailbench"
}

run_raccoon() {
    local repo="$EXTERNAL_DIR/raccoon"
    local dir="$repo/ref-c"
    local outfile="$THIRD_PARTY_DIR/raccoon.txt"
    local buildlog="$THIRD_PARTY_DIR/raccoon.build.log"

    mkdir -p "$THIRD_PARTY_DIR"
    : > "$buildlog"
    write_machine_info "$outfile"
    append_commit_info "$outfile"
    append_stability_info "$outfile"
    printf '=== Raccoon Performance ===\n\n' >> "$outfile"

    setup_repo "raccoon" "$RACCOON_URL" "$RACCOON_COMMIT"
    run_raccoon_impl "RACCOON_128_1" "$dir" "$outfile" "$buildlog"
    run_raccoon_impl "RACCOON_192_1" "$dir" "$outfile" "$buildlog"
    run_raccoon_impl "RACCOON_256_1" "$dir" "$outfile" "$buildlog"
    log "wrote $outfile"
}

run_bliss() {
    local repo="$EXTERNAL_DIR/galactics"
    local dir="$repo/Artifacts/CT-Implementation/bliss1_galactics/ref"
    local outfile="$THIRD_PARTY_DIR/bliss.txt"
    local buildlog="$THIRD_PARTY_DIR/bliss.build.log"
    local bench_dir="$dir/irs_tailbench"

    mkdir -p "$THIRD_PARTY_DIR"
    : > "$buildlog"
    write_machine_info "$outfile"
    append_commit_info "$outfile"
    append_stability_info "$outfile"
    printf '=== BLISS-I from GALACTICS Performance ===\n\n' >> "$outfile"

    setup_repo "galactics" "$GALACTICS_URL" "$GALACTICS_COMMIT"
    git -C "$repo" apply "$THIRD_PARTY_DIR/galactics-bliss-tailbench.patch"
    copy_tailbench_named_dir "$bench_dir"
    allow_bench_sign_runs_override "$bench_dir/cpucycles.h"
    write_speed_nist_c "$bench_dir/speed_tailbench.c"

    printf -- '-- ref/ BLISS-I --\n' >> "$outfile"
    run_logged "$buildlog" bash -c "cd '$dir' && \${CC:-cc} -Wall -Wextra -O3 -fomit-frame-pointer -I. -Iirs_tailbench -DBENCH_SIGN_RUNS=$BENCH_SIGN_RUNS -o irs_tailbench/speed_tailbench irs_tailbench/speed_tailbench.c irs_tailbench/cpucycles.c sign.c poly.c ntt.c reduce.c fips202.c packing.c gaussian_sampling.c rej_sampling.c randombytes.c -lm"
    append_run_with_blank "$outfile" "$bench_dir/speed_tailbench"
    log "wrote $outfile"
}


run_hawk_impl() {
    local label="$1"
    local dir="$2"
    local mode_n="$3"
    local security="$4"
    local outfile="$5"
    local buildlog="$6"
    local bench_dir="$dir/irs_tailbench"

    printf -- '-- %s --\n' "$label" >> "$outfile"
    copy_tailbench_named_dir "$bench_dir"
    allow_bench_sign_runs_override "$bench_dir/cpucycles.h"
    write_speed_hawk_c "$bench_dir/speed_tailbench.c"
    run_logged "$buildlog" make -j1 -s -C "$dir" clean
    run_logged "$buildlog" make -j1 -s -C "$dir"
    run_logged "$buildlog" bash -c "cd '$dir' && \${CC:-c99} -Wall -Wextra -Wshadow -Wundef -O2 -I. -Iirs_tailbench -DBENCH_SIGN_RUNS=$BENCH_SIGN_RUNS -DHAWK_N=$mode_n -DSECURITY_BITS=$security -o irs_tailbench/speed_tailbench irs_tailbench/speed_tailbench.c irs_tailbench/cpucycles.c api.o rng.o hawk_kgen.o hawk_sign.o hawk_vrfy.o ng_fxp.o ng_hawk.o ng_mp31.o ng_ntru.o ng_poly.o ng_zint31.o sha3.o"
    append_run_with_blank "$outfile" "$bench_dir/speed_tailbench"
}

run_hawk() {
    local repo="$EXTERNAL_DIR/hawk-dev"
    local outfile="$THIRD_PARTY_DIR/hawk.txt"
    local buildlog="$THIRD_PARTY_DIR/hawk.build.log"

    mkdir -p "$THIRD_PARTY_DIR"
    : > "$buildlog"
    write_machine_info "$outfile"
    append_commit_info "$outfile"
    append_stability_info "$outfile"
    printf '=== HAWK Performance ===\n\n' >> "$outfile"

    setup_repo "hawk-dev" "$HAWK_URL" "$HAWK_COMMIT"
    run_logged "$buildlog" python3 "$repo/build.py" ref avx2

    run_hawk_impl "ref/ HAWK-512" "$repo/NIST/Reference_Implementation/hawk512" 512 128 "$outfile" "$buildlog"
    run_hawk_impl "ref/ HAWK-1024" "$repo/NIST/Reference_Implementation/hawk1024" 1024 256 "$outfile" "$buildlog"
    if has_cpu_flag avx2; then
        run_hawk_impl "avx2/ HAWK-512" "$repo/NIST/Optimized_Implementation/avx2/hawk512" 512 128 "$outfile" "$buildlog"
        run_hawk_impl "avx2/ HAWK-1024" "$repo/NIST/Optimized_Implementation/avx2/hawk1024" 1024 256 "$outfile" "$buildlog"
    else
        append_skip "$outfile" '-- avx2/ --'
        append_skip "$outfile" 'SKIP: CPU does not report avx2 support.'
    fi
    log "wrote $outfile"
}


dilithium_xu_cpu_ok() {
    local impl="$1"

    if [ "$impl" = "avx2" ]; then
        has_cpu_flag avx2
        return $?
    fi
    has_cpu_flag avx512f && has_cpu_flag avx512vl && \
        has_cpu_flag avx512cd && has_cpu_flag avx512bw && \
        has_cpu_flag avx512dq && has_cpu_flag avx512vbmi2 && \
        has_cpu_flag avx512vbmi
}

run_dilithium_xu_impl() {
    local label="$1"
    local dir="$2"
    local impl="$3"
    local mode_name="$4"
    local pk_macro="$5"
    local sk_macro="$6"
    local sig_macro="$7"
    local keypair_fn="$8"
    local sign_fn="$9"
    local verify_fn="${10}"
    local outfile="${11}"
    local buildlog="${12}"

    printf -- '-- %s --\n' "$label" >> "$outfile"
    if ! dilithium_xu_cpu_ok "$impl"; then
        append_skip "$outfile" "SKIP: CPU does not report required $impl support."
        return 0
    fi

    copy_tailbench_test_dir "$dir/test"
    set_bench_sign_runs_default "$dir/test/cpucycles.h" "$BENCH_SIGN_RUNS"
    write_speed_dilithium_xu_c "$dir/test/test_speed.c" "$mode_name" "$impl" \
        "$pk_macro" "$sk_macro" "$sig_macro" "$keypair_fn" "$sign_fn" "$verify_fn"
    run_logged "$buildlog" make -j1 -s -C "$dir" clean
    run_logged "$buildlog" make -j1 -s -C "$dir"
    append_run_with_blank "$outfile" "$dir/test_speed"
}

run_dilithium_xu() {
    local repo="$EXTERNAL_DIR/dilithium-xu"
    local outfile="$THIRD_PARTY_DIR/dilithium-xu.txt"
    local buildlog="$THIRD_PARTY_DIR/dilithium-xu.build.log"

    mkdir -p "$THIRD_PARTY_DIR"
    : > "$buildlog"
    write_machine_info "$outfile"
    append_commit_info "$outfile"
    append_stability_info "$outfile"
    printf '=== Xu et al. Dilithium Performance ===\n\n' >> "$outfile"

    setup_repo "dilithium-xu" "$DILITHIUM_XU_URL" "$DILITHIUM_XU_COMMIT"

    run_dilithium_xu_impl "avx2/ Dilithium2" "$repo/dilithium2_avx2" "avx2" "Dilithium2-Xu" \
        "DILITHIUM2_CRYPTO_PUBLICKEYBYTES" "DILITHIUM2_AVX2_CRYPTO_SECRETKEYBYTES" \
        "DILITHIUM2_AVX2_CRYPTO_BYTES" "dilithium2_keypair_avx2" \
        "dilithium2_sign_avx2" "dilithium2_verify_avx2" "$outfile" "$buildlog"
    run_dilithium_xu_impl "avx512/ Dilithium2" "$repo/dilithium2_avx512" "avx512" "Dilithium2-Xu" \
        "PQCLEAN_DILITHIUM2_AVX2_CRYPTO_PUBLICKEYBYTES" "PQCLEAN_DILITHIUM2_AVX2_CRYPTO_SECRETKEYBYTES" \
        "PQCLEAN_DILITHIUM2_AVX2_CRYPTO_BYTES" "XURQ_AVX512_crypto_sign_keypair" \
        "XURQ_AVX512_crypto_sign_signature" "XURQ_AVX512_crypto_sign_verify" "$outfile" "$buildlog"
    run_dilithium_xu_impl "avx2/ Dilithium3" "$repo/dilithium3_avx2" "avx2" "Dilithium3-Xu" \
        "DILITHIUM3_AVX2_CRYPTO_PUBLICKEYBYTES" "DILITHIUM3_AVX2_CRYPTO_SECRETKEYBYTES" \
        "DILITHIUM3_AVX2_CRYPTO_BYTES" "dilithium3_keypair_avx2" \
        "dilithium3_sign_avx2" "dilithium3_verify_avx2" "$outfile" "$buildlog"
    run_dilithium_xu_impl "avx512/ Dilithium3" "$repo/dilithium3_avx512" "avx512" "Dilithium3-Xu" \
        "PQCLEAN_DILITHIUM3_AVX2_CRYPTO_PUBLICKEYBYTES" "PQCLEAN_DILITHIUM3_AVX2_CRYPTO_SECRETKEYBYTES" \
        "PQCLEAN_DILITHIUM3_AVX2_CRYPTO_BYTES" "XURQ_AVX512_crypto_sign_keypair" \
        "XURQ_AVX512_crypto_sign_signature" "XURQ_AVX512_crypto_sign_verify" "$outfile" "$buildlog"
    run_dilithium_xu_impl "avx2/ Dilithium5" "$repo/dilithium5_avx2" "avx2" "Dilithium5-Xu" \
        "PQCLEAN_DILITHIUM5_AVX2_CRYPTO_PUBLICKEYBYTES" "PQCLEAN_DILITHIUM5_AVX2_CRYPTO_SECRETKEYBYTES" \
        "PQCLEAN_DILITHIUM5_AVX2_CRYPTO_BYTES" "PQCLEAN_DILITHIUM5_AVX2_crypto_sign_keypair" \
        "PQCLEAN_DILITHIUM5_AVX2_crypto_sign_signature" "PQCLEAN_DILITHIUM5_AVX2_crypto_sign_verify" "$outfile" "$buildlog"
    run_dilithium_xu_impl "avx512/ Dilithium5" "$repo/dilithium5_avx512" "avx512" "Dilithium5-Xu" \
        "DILITHIUM5_AVX2_CRYPTO_PUBLICKEYBYTES" "DILITHIUM5_AVX2_CRYPTO_SECRETKEYBYTES" \
        "DILITHIUM5_AVX2_CRYPTO_BYTES" "dilithium5_sign_keypair_avx512" \
        "dilithium5_sign_signature_avx512" "dilithium5_sign_verify_avx512" "$outfile" "$buildlog"
    log "wrote $outfile"
}

if [ "$#" -eq 0 ]; then
    usage
    exit 1
fi

run_dil=0
run_fal=0
run_hae=0
run_ntr=0
run_rac=0
run_bli=0
run_haw=0
run_dil_xu=0
for arg in "$@"; do
    case "$arg" in
        --all)
            run_dil=1
            run_fal=1
            run_hae=1
            run_ntr=1
            run_rac=1
            run_bli=1
            run_haw=1
            run_dil_xu=1
            ;;
        --dilithium)
            run_dil=1
            ;;
        --falcon)
            run_fal=1
            ;;
        --haetae)
            run_hae=1
            ;;
        --ntruplus)
            run_ntr=1
            ;;
        --raccoon)
            run_rac=1
            ;;
        --bliss)
            run_bli=1
            ;;
        --hawk)
            run_haw=1
            ;;
        --Dilithium-Xu|--dilithium-xu)
            run_dil_xu=1
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            usage
            exit 1
            ;;
    esac
done

check_required_tools

mkdir -p "$EXTERNAL_DIR" "$THIRD_PARTY_DIR"

# Keep all benchmark execution serial; -j1 also neutralizes inherited MAKEFLAGS.
if [ "$run_dil" -eq 1 ]; then
    run_dilithium
fi
if [ "$run_fal" -eq 1 ]; then
    run_falcon
fi
if [ "$run_hae" -eq 1 ]; then
    run_haetae
fi
if [ "$run_ntr" -eq 1 ]; then
    run_ntruplus
fi
if [ "$run_rac" -eq 1 ]; then
    run_raccoon
fi
if [ "$run_bli" -eq 1 ]; then
    run_bliss
fi
if [ "$run_haw" -eq 1 ]; then
    run_hawk
fi
if [ "$run_dil_xu" -eq 1 ]; then
    run_dilithium_xu
fi
