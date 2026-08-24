# Complete n=512 negacyclic NTT for q=64513, 16-bit AVX2, Kyber-style.
#   make            -> build + run correctness test, then the speed benchmark
#   make test       -> build + run correctness test only
#   make speed      -> build + run the scalar-vs-vector speed benchmark
#   make gen        -> regenerate ntt.S + ntt_consts.c from the verified model
#   make clean
CC      = gcc
CFLAGS  = -O3 -mavx2 -Wall
AVX512FLAGS = -mavx512f -mavx512bw

ASM     = ntt.S ntt_consts.c
ASM512  = ntt_avx512.S ntt_consts_avx512.c
REF     = ntt_ref.c
SPEEDINF= cpucycles.c speed_print.c

all: test speed

test: test_ntt
	./test_ntt

speed: speed_ntt
	./speed_ntt

test512:
	@if [ "$(AVX512)" != "1" ]; then echo "Use: make AVX512=1 test512"; exit 1; fi
	$(MAKE) test_ntt512 AVX512=1
	./test_ntt512

speed512:
	@if [ "$(AVX512)" != "1" ]; then echo "Use: make AVX512=1 speed512"; exit 1; fi
	$(MAKE) speed_ntt512 AVX512=1
	./speed_ntt512

test_ntt: test_ntt.c $(ASM) $(REF) ntt_ref.h
	$(CC) $(CFLAGS) -o $@ test_ntt.c $(ASM) $(REF)

speed_ntt: speed_ntt.c $(ASM) $(REF) $(SPEEDINF) ntt_ref.h cpucycles.h speed_print.h
	$(CC) $(CFLAGS) -o $@ speed_ntt.c $(ASM) $(REF) $(SPEEDINF)

test_ntt512: test_ntt.c $(ASM) $(ASM512) $(REF) ntt_ref.h
	$(CC) $(CFLAGS) $(AVX512FLAGS) -DTEST_AVX512 -o $@ test_ntt.c $(ASM) $(ASM512) $(REF)

speed_ntt512: speed_ntt.c $(ASM) $(ASM512) $(REF) $(SPEEDINF) ntt_ref.h cpucycles.h speed_print.h
	$(CC) $(CFLAGS) $(AVX512FLAGS) -DTEST_AVX512 -o $@ speed_ntt.c $(ASM) $(ASM512) $(REF) $(SPEEDINF)

# regenerate the assembly + constant tables (also self-verifies the model)
gen:
	python3 gen_ntt.py

clean:
	rm -f test_ntt speed_ntt test_ntt512 speed_ntt512

.PHONY: all test speed test512 speed512 gen clean
