# Complete n=256 negacyclic NTT for q=45569, 16-bit AVX2, Kyber-style.
#   make            -> build + run correctness test, then the speed benchmark
#   make test       -> build + run correctness test only
#   make speed      -> build + run the scalar-vs-vector speed benchmark
#   make gen        -> regenerate ntt.S + ntt_consts.c from the verified model
#   make clean
CC      = gcc
CFLAGS  = -O3 -mavx2 -Wall

ASM     = ntt.S ntt_consts.c
REF     = ntt_ref.c
SPEEDINF= cpucycles.c speed_print.c

all: test speed

test: test_ntt
	./test_ntt

speed: speed_ntt
	./speed_ntt

test_ntt: test_ntt.c $(ASM) $(REF) ntt_ref.h
	$(CC) $(CFLAGS) -o $@ test_ntt.c $(ASM) $(REF)

speed_ntt: speed_ntt.c $(ASM) $(REF) $(SPEEDINF) ntt_ref.h cpucycles.h speed_print.h
	$(CC) $(CFLAGS) -o $@ speed_ntt.c $(ASM) $(REF) $(SPEEDINF)

# regenerate the assembly + constant tables (also self-verifies the model)
gen:
	python3 gen_ntt.py

clean:
	rm -f test_ntt speed_ntt

.PHONY: all test speed gen clean
