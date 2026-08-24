# ApproxExp design-space sweep

Replays the production ApproxExp point from `docs/design-notes/x64_ApproxExp.md`: **B=8, degree-6**, with a local k-way throughput/latency sweep.

- `gen_explore.py BETA LOG OUT` — emits a config header (`RCP`, `LB`, table, Q63 poly parsed from a GALACTICS log) for table size `B = 2^BETA`.
- `explore_exp.c` — scalar + `KWAY`-way ApproxExp; measures latency-bound and k-way throughput, with a bit-identity + libm-accuracy correctness gate.
- `sweep.sh` — runs the production B=8 point over k in {1,2,3,4,5,6,8}.

Precision for the production point is reproducible from `tools/LogPolyApprox/expx_7div8mulln2_ln2_43_64.txt`; cycle counts are machine specific (measured on one i7-11700K Rocket Lake native-Ubuntu test system with Turbo Boost and Hyper-Threading disabled, `-O3 -march=native`, rdtsc).

Run: `sh sweep.sh`
