#!/bin/sh
# (B, degree) x k-way throughput/latency sweep for the sampler ApproxExp.
# Replays the production ApproxExp point from docs/design-notes/LITHIUM_x64_ApproxExp.md.
#   - precision comes from the committed production log
#     tools/LogPolyApprox/expx_7div8mulln2_ln2_43_64.txt;
#   - throughput/latency are measured here (machine-specific; -O3 -march=native).
# Usage:  sh Lithium-Code/ref/test/explore/sweep.sh
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(python3 -c 'import pathlib,sys; print(pathlib.Path(sys.argv[1]).resolve().parents[2])' "$HERE")
REF=$ROOT/ref
echo "============ ApproxExp production B=8 k-way sweep ============"
for spec in "3:7div8"; do
  beta=${spec%%:*}; frac=${spec##*:}
  log=$ROOT/tools/LogPolyApprox/expx_${frac}mulln2_ln2_43_64.txt
  python3 "$HERE/gen_explore.py" "$beta" "$log" "$HERE/cfg.h" >/dev/null
  for K in 1 2 3 4 5 6 8; do
    gcc -O3 -march=native -std=c99 -I"$HERE" -I"$REF/test" -DKWAY=$K \
        "$HERE/explore_exp.c" "$REF/test/cpucycles.c" \
        -lm -o "$HERE/run" || { echo "build fail beta=$beta K=$K"; continue; }
    "$HERE/run"
  done
  echo "---"
done
rm -f "$HERE/cfg.h" "$HERE/run"
