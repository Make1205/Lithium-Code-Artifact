# Scalar ApproxExp Notes

The scalar ApproxExp path uses a lookup-table plus polynomial approximation for `exp(-a)` in the Gaussian rejection sampler. `tools/gen_approx_exp.py` regenerates the default integer constants embedded in `ref/approx_exp.h`; optimized builds may opt into the AVX-specific kernels documented in `docs/design-notes/AVXExp.md`.

The production parameters use table size `B = 8`, degree 6, and a multi-way evaluation schedule selected by local accuracy and throughput sweeps. The committed generators check the relevant relative-error target and keep all constants reproducible through `make tables` and `make check-consts`.

Benchmark numbers in the design documents are measurements on one i7-11700K (Rocket Lake, 2021) native-Ubuntu test system with Turbo Boost and Hyper-Threading disabled, and are reported only as reproducibility context, not as a design restriction.
