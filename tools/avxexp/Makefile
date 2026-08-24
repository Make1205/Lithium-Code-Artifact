# tools/avxexp -- exp(-a) kernels: scalar-int, AVX2-double(4), AVX-512 IFMA(8)
CC     ?= gcc
CFLAGS ?= -O2 -mavx2 -mfma -mavx512f -mavx512dq -mavx512ifma
LDLIBS ?= -lm

all: bench ifma_ref

bench: bench.c          # correctness + throughput + dudect-style constant-time
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

ifma_ref: ifma_ref.c    # minimal scalar reference (no AVX-512 needed)
	$(CC) -O2 -o $@ $< $(LDLIBS)

run: bench
	./bench | tee LogAVXApprox/bench_$(shell uname -n).txt

clean:
	rm -f bench ifma_ref

.PHONY: all run clean
