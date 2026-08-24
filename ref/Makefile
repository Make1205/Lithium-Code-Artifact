# Makefile for the Lithium signature reference implementation.
#
#   make                 # build all test binaries for all modes -> out/
#   make check           # build + run every test for every mode
#   make tables          # regenerate all generated headers from ../tools
#   make clean
#
# MODES selects the parameter sets (default: all three).
# Generated constants live in @@AUTOGEN@@ regions inside hand-written files
# (params.h, ntt.h, approx_exp.h, irs_rans.h, sampler.c); they are committed, so a
# plain build needs no Python.  `make tables` regenerates every region via ../tools.

MODES  ?= 120 180 260
TESTS ?= test_ntt t_poly t_exp speed_approx_exp test_sampler test_irs \
         test_mod2q test_rans test_keys test_sign test_kat
CC     ?= gcc
CFLAGS ?= -O3 -fstack-protector-strong -Wall -Wextra -std=c99 -I. -Itest
LDLIBS ?= -lm
BINDIR  = out
CT_OPT ?= -O3

# Per-test source lists (POSIX-sh strings; test_kat omits randombytes.c -- it
# provides its own deterministic DRBG randombytes()).
SH_SRCS = \
'test_ntt=ntt.c reduce.c' \
't_poly=polyvec.c ntt.c reduce.c fips202.c' \
't_exp=approx_exp.c' \
'speed_approx_exp=approx_exp.c test/cpucycles.c' \
'test_sampler=sampler.c approx_exp.c fips202.c symmetric-shake.c' \
'test_irs=irs.c polyvec.c ntt.c reduce.c approx_exp.c fips202.c symmetric-shake.c' \
'test_mod2q=rounding.c reduce.c' \
'test_rans=irs_rans.c' \
'test_keys=polyvec.c ntt.c reduce.c sampler.c approx_exp.c fips202.c symmetric-shake.c randombytes.c' \
'test_sign=ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c randombytes.c' \
'test_kat=ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c' \
'parser_negative=ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c' \
'security_fault_injection=ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c' \
'trace_kat=ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c' \
'diff_fuzz=ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c' \
'dudect_components=polyvec.c ntt.c reduce.c rounding.c fips202.c test/cpucycles.c'

.PHONY: all check tables check-consts clean
all:
	@mkdir -p $(BINDIR)
	@for m in $(MODES); do \
	  for entry in $(SH_SRCS); do \
	    t=$${entry%%=*}; srcs=$${entry#*=}; \
	    for want in $(TESTS); do [ "$$want" = "$$t" ] || continue; \
	      stats=""; case "$$t" in trace_kat) stats="-DTEST_RANS_STATS";; diff_fuzz) stats="-DDIFF_FUZZ_N=$(DIFF_FUZZ_N)";; dudect_components) stats="-DDUDECT_COMPONENT_N=$(DUDECT_COMPONENT_N)";; esac; \
	      $(CC) $(CFLAGS) $$stats -DLITHIUM_MODE=$$m test/$$t.c $$srcs $(LDLIBS) \
	        -o $(BINDIR)/$${t}_$${m} || exit 1; \
	    done; \
	  done; \
	done
	@echo "built: [$(TESTS)] x [$(MODES)]"

check: all
	@for m in $(MODES); do \
	  echo "=== Lithium-$$m ==="; \
	  for want in $(TESTS); do \
	    [ -x $(BINDIR)/$${want}_$${m} ] \
	      || { echo "missing test binary: $(BINDIR)/$${want}_$${m}"; exit 1; }; \
	    $(BINDIR)/$${want}_$${m} || exit 1; \
	  done; \
	done

# Per-component time + randomness profile of crypto_sign (scalar reference).
# Compare against the AVX2 build's `make profile` in ../avx2.
PROF_SRC = ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c \
           packing.c irs_rans.c sign.c fips202.c symmetric-shake.c randombytes.c test/prof.c

RANS_STATS_N ?= 512
RANS_FAILRATE_SRC = ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c \
                    irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c randombytes.c
.PHONY: check-rans-failrate
check-rans-failrate:
	@mkdir -p $(BINDIR)
	@for m in $(MODES); do \
	  $(CC) $(CFLAGS) -DTEST_RANS_STATS -DRANS_STATS_N=$(RANS_STATS_N) \
	    -DLITHIUM_MODE=$$m test/test_rans_failrate.c $(RANS_FAILRATE_SRC) \
	    $(LDLIBS) -o $(BINDIR)/test_rans_failrate_$$m || exit 1; \
	  $(BINDIR)/test_rans_failrate_$$m || exit 1; \
	done

DIFF_FUZZ_N ?= 32
DUDECT_COMPONENT_N ?= 512

DUDECT_N ?= 2000
DUDECT_SIGN_SRC = ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c \
                  irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c test/cpucycles.c
.PHONY: dudect-sign
dudect-sign:
	@mkdir -p $(BINDIR)
	@for m in $(MODES); do \
	  $(CC) $(CFLAGS) -DDUDECT_N=$(DUDECT_N) -DLITHIUM_MODE=$$m \
	    test/dudect_sign.c $(DUDECT_SIGN_SRC) $(LDLIBS) \
	    -o $(BINDIR)/dudect_sign_$$m || exit 1; \
	  $(BINDIR)/dudect_sign_$$m || exit 1; \
	done

.PHONY: profile
profile:
	@mkdir -p $(BINDIR)
	@for m in $(MODES); do \
	  $(CC) $(CFLAGS) -DPROF_TIME -DPROF_RAND -DLITHIUM_MODE=$$m \
	    test/speed_profile.c $(PROF_SRC) $(LDLIBS) -o $(BINDIR)/speed_profile_$$m || exit 1; \
	  $(BINDIR)/speed_profile_$$m; \
	  echo; \
	done

# End-to-end keygen/sign/verify cycle benchmark (scalar reference).
.PHONY: speed
speed:
	@mkdir -p $(BINDIR)
	@for m in $(MODES); do \
	  $(CC) $(CFLAGS) -DLITHIUM_MODE=$$m test/speed_sign.c $(PROF_SRC:test/prof.c=) test/cpucycles.c \
	    $(LDLIBS) -o $(BINDIR)/speed_sign_$$m || exit 1; \
	  $(BINDIR)/speed_sign_$$m; \
	  echo; \
	done

# Human-readable machine-info preamble shared by speed_txt / profile_txt.
MACHINFO = printf ' date    : '; date -u '+%Y-%m-%dT%H:%M:%SZ'; \
           printf ' kernel  : '; uname -sr; \
           printf ' cpu     : '; ( lscpu 2>/dev/null | sed -n 's/^Model name:[[:space:]]*//p' | head -1 || grep -m1 'model name' /proc/cpuinfo | sed 's/.*: //' ); \
           printf ' cc      : '; $(CC) --version | head -1; \
           printf ' cflags  : %s\n' '$(CFLAGS)'

# Dump end-to-end speed (tail/median rdtsc cycles) for every mode to test/speed.txt,
# prefixed with a machine-info header.  Pure-scalar reference baseline.
.PHONY: speed_txt
speed_txt:
	@{ echo "=== Lithium-Code  ref/  end-to-end speed (tail/median rdtsc cycles, speed_sign) ==="; \
	   $(MACHINFO); \
	   echo " stability metrics:"; \
	   echo "   p25, p50, p75, p90, p99: nearest-rank percentiles of keygen/sign cycles."; \
	   echo "   StQ1, StQ2, StQ3: stabilized quartiles, i.e., [1/8,3/8], [3/8,5/8], [5/8,7/8] means."; \
	   echo "     Reference: Daniel J. Bernstein, Robust statistics for rejection-sampling timings, https://cr.yp.to/papers/rsrst-20250727.pdf."; \
	   echo "   stddev: population standard deviation of keygen/sign cycles; lower means less run-to-run spread."; \
	   echo "   cv: stddev/avg; rough guide: <1% very stable, 1-5% stable, 5-15% noisy, >15% unstable or long-tailed."; \
	   echo "   outlier_count(>p50*1.5): keygen/sign runs slower than 1.5x the median."; \
	   echo "   p90/p50 and p99/p50: tail-inflation ratios; values near 1.00x mean a short tail."; \
	   echo; \
	   echo "-- pure-scalar reference --"; } > test/speed.txt
	@$(MAKE) --no-print-directory speed >> test/speed.txt
	@echo "wrote test/speed.txt"; cat test/speed.txt

# Dump the scalar approx_exp microbenchmark to
# test/speed_exp.txt, prefixed with a machine-info header.
.PHONY: exp_txt
exp_txt:
	@$(MAKE) --no-print-directory MODES="$(firstword $(MODES))" TESTS=speed_approx_exp all >/dev/null
	@{ echo "=== Lithium-Code  ref/  approx_exp speed (median rdtsc cycles, speed_approx_exp) ==="; \
	   $(MACHINFO); \
	   echo; \
	   echo "exp performance is mode-independent; using mode $(firstword $(MODES))."; \
	   echo "-- pure-scalar reference --"; \
	   m="$(firstword $(MODES))"; echo "=== Lithium-$$m ==="; $(BINDIR)/speed_approx_exp_$$m; echo; } > test/speed_exp.txt
	@echo "wrote test/speed_exp.txt"; cat test/speed_exp.txt

# Dump the per-component profile (cycles/op for keygen/sign/verify) for every
# mode to test/profiling.txt, prefixed with the same machine-info header.
.PHONY: profile_txt
profile_txt:
	@{ echo "=== Lithium-Code  ref/  per-component profile (cycles/op, make profile) ==="; \
	   $(MACHINFO); echo; \
	   echo "-- pure-scalar reference --"; } > test/profiling.txt
	@$(MAKE) --no-print-directory profile >> test/profiling.txt
	@echo "wrote test/profiling.txt"

tables:
	cd ../tools && python3 gen_zetas.py && python3 gen_bounds.py \
	  && python3 gen_approx_exp.py && python3 gen_rcdt.py \
	  && python3 gen_rans_tables.py && python3 SigSize.py && python3 irs_ref.py

# Reproducible-constants gate for the reference generator output that is shared
# by ref/avx2/avx512.  Keep this separate from `tables`: `tables` regenerates,
# while `check-consts` proves the committed bytes do not drift.
check-consts:
	@mkdir -p $(BINDIR)
	@echo "== approx_exp.h (@@AUTOGEN:approx_exp@@; tools/gen_approx_exp.py) =="
	@cp approx_exp.h $(BINDIR)/approx_exp.h.bak
	@( cd ../tools && python3 gen_approx_exp.py >/dev/null 2>&1 ) \
	  || { echo "  FAIL: gen_approx_exp.py"; rm -f $(BINDIR)/approx_exp.h.bak; exit 1; }
	@cmp -s approx_exp.h $(BINDIR)/approx_exp.h.bak \
	  && echo "  PASS: approx_exp.h reproducible" \
	  || { echo "  FAIL: approx_exp.h drift"; rm -f $(BINDIR)/approx_exp.h.bak; exit 1; }
	@rm -f $(BINDIR)/approx_exp.h.bak

.PHONY: ct-scan
ct-scan:
	@python3 ../tools/ct_scan.py --backends ref --modes "$(MODES)" --ccs "$(CC)" --opts="$(CT_OPT)" --report ref/test/ct_scan_ref.txt

clean:
	rm -rf $(BINDIR)
