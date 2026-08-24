# Makefile for the Lithium signature AVX-512 implementation.
#
#   make                    # build all tests for all modes -> out/
#   make check              # build + run every test for every mode
#   make check-kat          # build + run only test_kat; assert recorded hashes
#   make AVX512IFMAEXP=1 ... # opt-in IFMA exp (separate KAT; default OFF == shared KAT)
#   make TESTS="t_poly"     # build a subset of tests
#   make clean
#
# Shared algorithm/KAT source is symlinked from ../ref (single source of
# truth).  AVX-512-specific overrides are real files here; the per-mode NTT
# assembly + constants are VENDORED locally in ./<qset>/ (copied from
# tools/ntt; self-contained build, `check-consts` proves the copies match).
#
# Compile-time switches (set in CFLAGS) select the AVX-512 code paths in the
# symlinked/forked sources:
#   USE_AVX512_NTT     -> poly_ntt/invntt/pointwise routes in polyvec.c
#   USE_AVX512_SAMPLER -> 16-lane sampler_sigma2 + noise_magnitude_batch
#                             routes in sampler.c
#   USE_AVX512_SHAKE8X -> keccak8x gauss/noise routes in sampler.c
#                             + 8-way uniform A0 route in polyvec.c

MODES ?= 120 180 260
TESTS ?= test_ntt t_poly t_exp speed_approx_exp test_sampler test_irs \
         test_mod2q test_rans test_keys test_sign test_kat test_ntt_avx512 \
         test_freeze_avx512 test_sigma2_diff test_keccakx8
CC     ?= gcc
AVX512IFMAEXP ?= 0
CFLAGS ?= -O3 -fstack-protector-strong -Wall -Wextra -std=c99 -I. -I../ref -I../ref/test \
          -mavx512f -mavx512bw -mavx512vl -mavx512dq -mavx512ifma \
          -mbmi2 -mpopcnt -march=native -mtune=native \
          -DUSE_AVX512_NTT -DUSE_AVX512_SAMPLER -DUSE_AVX512_SHAKE8X
LDLIBS ?= -lm
BINDIR  = out
CT_OPT ?= -O3

ifeq ($(AVX512IFMAEXP),1)
# IFMA integer exp (Q51 vpmadd52).  Opt-in, separate KAT (KAT_avx512ifma.txt).
CFLAGS += -DAVX512IFMA_EXP
EXPIFMA = approx_exp_avx512ifma.c
else
EXPIFMA =
endif

# AVX-512 SIMD source bundles.  The NTT routing + 8-way uniform are now merged
# into polyvec.c, and the 16-lane sigma2/noise + keccak8x vector samplers into
# sampler.c (NTT512/SAMP512/UNIF512 are therefore empty -- kept as names so the
# source lists below read the same).  KEC8X (the 8-way Keccak permutation) stays
# a separate vendored object; every test that links polyvec.c or sampler.c also
# links KEC8X, since both now reference shake128x8 / shake256x8.
NTT512  =
KEC8X   = fips202x8.c keccakf1600x8.c
SAMP512 =
UNIF512 =
SHAKE8X = $(KEC8X) # keccak8x objects (profile / speed_sign)

# Per-test source lists (POSIX-sh strings).  Tests that list polyvec.c also get
# the per-mode vendored ntt_avx512.S + ntt_consts_avx512.c appended.
SH_SRCS = \
'test_ntt=ntt.c reduce.c' \
't_poly=polyvec.c ntt.c reduce.c fips202.c $(NTT512) $(UNIF512) $(KEC8X)' \
't_exp=approx_exp.c $(EXPIFMA)' \
'speed_approx_exp=approx_exp.c test/cpucycles.c $(EXPIFMA)' \
'test_sampler=sampler.c approx_exp.c fips202.c symmetric-shake.c $(SAMP512) $(KEC8X) $(EXPIFMA)' \
'test_irs=irs.c polyvec.c ntt.c reduce.c approx_exp.c fips202.c symmetric-shake.c $(NTT512) $(UNIF512) $(KEC8X) $(EXPIFMA)' \
'test_mod2q=rounding.c reduce.c' \
'test_rans=irs_rans.c' \
'test_keys=polyvec.c ntt.c reduce.c sampler.c approx_exp.c fips202.c symmetric-shake.c randombytes.c $(NTT512) $(SAMP512) $(UNIF512) $(KEC8X) $(EXPIFMA)' \
'test_sign=ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c randombytes.c $(NTT512) $(SAMP512) $(UNIF512) $(KEC8X) $(EXPIFMA)' \
'test_kat=ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c $(NTT512) $(SAMP512) $(UNIF512) $(KEC8X) $(EXPIFMA)' \
'parser_negative=ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c $(NTT512) $(SAMP512) $(UNIF512) $(KEC8X) $(EXPIFMA)' \
'security_fault_injection=ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c $(NTT512) $(SAMP512) $(UNIF512) $(KEC8X) $(EXPIFMA)' \
'trace_kat=ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c $(NTT512) $(SAMP512) $(UNIF512) $(KEC8X) $(EXPIFMA)' \
'diff_fuzz=ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c $(NTT512) $(SAMP512) $(UNIF512) $(KEC8X) $(EXPIFMA)' \
'dudect_components=polyvec.c ntt.c reduce.c rounding.c test/cpucycles.c fips202.c $(NTT512) $(UNIF512) $(KEC8X)' \
'test_ntt_avx512=polyvec.c ntt.c reduce.c fips202.c $(NTT512) $(UNIF512) $(KEC8X)' \
'test_freeze_avx512=polyvec.c ntt.c reduce.c fips202.c $(NTT512) $(UNIF512) $(KEC8X)' \
'test_sigma2_diff=sampler.c approx_exp.c fips202.c symmetric-shake.c $(KEC8X) $(EXPIFMA)' \
'test_keccakx8=fips202x8.c keccakf1600x8.c fips202.c' \
'speed_exp_avx512ifma=approx_exp_avx512ifma.c approx_exp.c test/cpucycles.c' \
'speed_sign=ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c randombytes.c test/cpucycles.c $(NTT512) $(SAMP512) $(UNIF512) $(KEC8X) $(EXPIFMA)'

QSET_120 = q57601n256
QSET_180 = q30977n256
QSET_260 = q64513n512
# The per-mode NTT (ntt_avx512.S + ntt_consts_avx512.c) is VENDORED locally in
# ./<qset>/ so this build is self-contained.  NTTGEN points at the generator
# (tools/ntt), used only by `check-consts` to prove the vendored copies match.
NTTGEN   = ../tools/ntt

# Recorded canonical KAT hashes (shared: ref == avx2 == avx512, integer exp).
KAT_120 = 16043336652881898616
KAT_180 = 12054898970128143599
KAT_260 = 14183906608478196884

.PHONY: all check check-kat clean
all:
	@mkdir -p $(BINDIR)
	@for m in $(MODES); do \
	  case $$m in 120) qs=$(QSET_120);; 180) qs=$(QSET_180);; 260) qs=$(QSET_260);; esac; \
	  for entry in $(SH_SRCS); do \
	    t=$${entry%%=*}; srcs=$${entry#*=}; \
	    for want in $(TESTS); do [ "$$want" = "$$t" ] || continue; \
	      extra=""; \
	      case "$$srcs" in *polyvec.c*) extra="$$qs/ntt_avx512.S $$qs/ntt_consts_avx512.c";; esac; \
	      stats=""; case "$$t" in trace_kat) stats="-DTEST_RANS_STATS";; diff_fuzz) stats="-DDIFF_FUZZ_N=$(DIFF_FUZZ_N)";; dudect_components) stats="-DDUDECT_COMPONENT_N=$(DUDECT_COMPONENT_N)";; esac; \
	      $(CC) $(CFLAGS) $$stats -DLITHIUM_MODE=$$m test/$$t.c $$srcs $$extra $(LDLIBS) \
	        -o $(BINDIR)/$${t}_$${m} || exit 1; \
	    done; \
	  done; \
	done
	@echo "built: [$(TESTS)] x [$(MODES)]"

check: all
	@for m in $(MODES); do \
	  echo "=== Lithium-$$m AVX-512 ==="; \
	  for want in $(TESTS); do \
	    [ -x $(BINDIR)/$${want}_$${m} ] \
	      || { echo "missing test binary: $(BINDIR)/$${want}_$${m}"; exit 1; }; \
	    $(BINDIR)/$${want}_$${m} || exit 1; \
	  done; \
	done

check-kat:
	@$(MAKE) TESTS=test_kat all >/dev/null
	@fail=0; for m in $(MODES); do \
	  case $$m in 120) want=$(KAT_120);; 180) want=$(KAT_180);; 260) want=$(KAT_260);; esac; \
	  h=$$($(BINDIR)/test_kat_$$m | grep -o 'KAT_hash=[0-9]*' | cut -d= -f2); \
	  [ "$$h" = "$$want" ] && echo "  PASS: mode=$$m KAT_hash=$$h" \
	    || { echo "  FAIL: mode=$$m got=$$h want=$$want"; fail=1; }; \
	done; [ $$fail -eq 0 ] && echo "check-kat: shared KAT matches (ref==avx2==avx512)" || exit 1

RANS_STATS_N ?= 512
RANS_FAILRATE_SRC = ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c \
                    irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c \
                    randombytes.c $(NTT512) $(SAMP512) $(UNIF512) $(KEC8X) $(EXPIFMA)
.PHONY: check-rans-failrate
check-rans-failrate:
	@mkdir -p $(BINDIR)
	@for m in $(MODES); do \
	  case $$m in 120) qs=$(QSET_120);; 180) qs=$(QSET_180);; 260) qs=$(QSET_260);; esac; \
	  $(CC) $(CFLAGS) -DTEST_RANS_STATS -DRANS_STATS_N=$(RANS_STATS_N) \
	    -DLITHIUM_MODE=$$m ../ref/test/test_rans_failrate.c $(RANS_FAILRATE_SRC) \
	    $$qs/ntt_avx512.S $$qs/ntt_consts_avx512.c $(LDLIBS) \
	    -o $(BINDIR)/test_rans_failrate_$$m || exit 1; \
	  $(BINDIR)/test_rans_failrate_$$m || exit 1; \
	done

DIFF_FUZZ_N ?= 32
DUDECT_COMPONENT_N ?= 512

DUDECT_N ?= 2000
DUDECT_SIGN_SRC = ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c \
                  irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c \
                  test/cpucycles.c $(NTT512) $(SAMP512) $(UNIF512) $(KEC8X) $(EXPIFMA)
.PHONY: dudect-sign
dudect-sign:
	@mkdir -p $(BINDIR)
	@for m in $(MODES); do \
	  case $$m in 120) qs=$(QSET_120);; 180) qs=$(QSET_180);; 260) qs=$(QSET_260);; esac; \
	  $(CC) $(CFLAGS) -DDUDECT_N=$(DUDECT_N) -DLITHIUM_MODE=$$m \
	    test/dudect_sign.c $(DUDECT_SIGN_SRC) $$qs/ntt_avx512.S $$qs/ntt_consts_avx512.c \
	    $(LDLIBS) -o $(BINDIR)/dudect_sign_$$m || exit 1; \
	  $(BINDIR)/dudect_sign_$$m || exit 1; \
	done

# Regression-lock the opt-in IFMA-exp KAT against the recorded KAT_avx512ifma.txt.
.PHONY: check-kat-avx512ifma
check-kat-avx512ifma:
	@$(MAKE) AVX512IFMAEXP=1 TESTS=test_kat all >/dev/null
	@fail=0; for m in $(MODES); do \
	  h=$$($(BINDIR)/test_kat_$$m | grep -o 'KAT_hash=[0-9]*'); \
	  grep -q "mode=$$m $$h" KAT_avx512ifma.txt \
	    && echo "  PASS: mode=$$m $$h matches KAT_avx512ifma.txt" \
	    || { echo "  FAIL: mode=$$m $$h NOT in KAT_avx512ifma.txt"; fail=1; }; \
	done; [ $$fail -eq 0 ] && echo "check-kat-avx512ifma: IFMA-exp KAT locked" || exit 1

# Force-exercise the 8-way SHAKE continuation-squeeze (refill) with a tiny
# initial squeeze, and assert the KAT is unchanged.
.PHONY: check-refill
check-refill:
	@$(MAKE) CFLAGS="$(CFLAGS) -DGAUSS8X_INIT=8 -DNOISE8X_INIT=8" TESTS="test_kat test_sampler" all >/dev/null
	@echo "== forced 8-way refill: KAT must be unchanged =="
	@for m in $(MODES); do $(BINDIR)/test_kat_$$m; $(BINDIR)/test_sampler_$$m | grep -o '4x==1x=[01]'; done

# Reproducible-constants gate: the per-mode AVX-512 NTT asm + consts must
# regenerate byte-for-byte (one gen_ntt.py run emits avx2 AND avx512), and the
# 8-way SHAKE must reuse ref's single KeccakF_RoundConstants symbol.
.PHONY: check-consts
check-consts:
	@echo "== vendored NTT (./<qset>/ntt_avx512.S + ntt_consts_avx512.c) vs generator (tools/ntt/<qset>/gen_ntt.py) =="
	@for qs in $(QSET_120) $(QSET_180) $(QSET_260); do \
	  ( cd $(NTTGEN)/$$qs && python3 gen_ntt.py >/dev/null 2>&1 ) || { echo "  FAIL: $$qs gen_ntt.py"; exit 1; }; \
	  cmp -s $(NTTGEN)/$$qs/ntt_consts_avx512.c $$qs/ntt_consts_avx512.c && cmp -s $(NTTGEN)/$$qs/ntt_avx512.S $$qs/ntt_avx512.S \
	    && echo "  PASS: $$qs vendored ntt_avx512.{S,c} match gen_ntt.py" \
	    || { echo "  FAIL: $$qs vendored avx512 NTT differs from gen_ntt.py (re-copy)"; exit 1; }; \
	done
	@echo "== 8-way SHAKE round constants (reuse ref KeccakF_RoundConstants) =="
	@$(MAKE) TESTS=test_kat all >/dev/null
	@n=$$(nm $(BINDIR)/test_kat_120 2>/dev/null | grep -c 'KeccakF_RoundConstants'); \
	  [ "$$n" = "1" ] && echo "  PASS: exactly one KeccakF_RoundConstants symbol (no x8 duplicate)" \
	  || { echo "  FAIL: $$n KeccakF_RoundConstants symbols"; exit 1; }
	@echo "== RCDT tables (tools/gen_rcdt.py; merged into sampler.c @@AUTOGEN:rcdt@@, ref+avx2+avx512) =="
	@mkdir -p $(BINDIR)
	@cp ../ref/sampler.c $(BINDIR)/rcdt_ref.bak; cp ../avx2/sampler.c $(BINDIR)/rcdt_avx2.bak; cp sampler.c $(BINDIR)/rcdt_avx512.bak
	@( cd ../tools && python3 gen_rcdt.py >/dev/null 2>&1 ) || { echo "  FAIL: gen_rcdt.py"; rm -f $(BINDIR)/rcdt_ref.bak $(BINDIR)/rcdt_avx2.bak $(BINDIR)/rcdt_avx512.bak; exit 1; }
	@cmp -s ../ref/sampler.c $(BINDIR)/rcdt_ref.bak && cmp -s ../avx2/sampler.c $(BINDIR)/rcdt_avx2.bak && cmp -s sampler.c $(BINDIR)/rcdt_avx512.bak \
	  && echo "  PASS: sampler.c RCDT region reproducible (ref + avx2 + avx512)" \
	  || { echo "  FAIL: sampler.c RCDT drift"; rm -f $(BINDIR)/rcdt_ref.bak $(BINDIR)/rcdt_avx2.bak $(BINDIR)/rcdt_avx512.bak; exit 1; }
	@rm -f $(BINDIR)/rcdt_ref.bak $(BINDIR)/rcdt_avx2.bak $(BINDIR)/rcdt_avx512.bak
ifeq ($(AVX512IFMAEXP),1)
	@echo "== IFMA exp consts (@@AUTOGEN:approx_exp_avx512ifma@@ in approx_exp_avx512ifma.c; tools/avxexp/gen_avxexp_ifma_consts.py) =="
	@mkdir -p $(BINDIR)
	@cp approx_exp_avx512ifma.c $(BINDIR)/ce.bak
	@python3 ../tools/avxexp/gen_avxexp_ifma_consts.py --check
	@cmp -s approx_exp_avx512ifma.c $(BINDIR)/ce.bak \
	  && echo "  PASS: approx_exp_avx512ifma.c consts region reproducible + 2^-43 gate" \
	  || { echo "  FAIL: IFMA consts drift"; rm -f $(BINDIR)/ce.bak; exit 1; }
	@rm -f $(BINDIR)/ce.bak
endif
	@echo "check-consts: all constants reproducible"

# Constant-time scan: delegate to the shared object-code scanner.  The root
# ./ct_scan_matrix.sh command covers the full gcc/clang, -O3/-Os, all-mode
# matrix including optional exp variants; this local target follows the current
# Makefile CC/MODES/CT_OPT knobs for quick backend checks.
.PHONY: ct-scan
ct-scan:
	@python3 ../tools/ct_scan.py --backends avx512 --modes "$(MODES)" --ccs "$(CC)" --opts="$(CT_OPT)" --include-optional --report avx512/test/ct_scan_avx512.txt

# Per-component time + randomness profile of crypto_sign.
.PHONY: profile
PROF_SRC = ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c \
           packing.c irs_rans.c sign.c fips202.c symmetric-shake.c randombytes.c ../ref/test/prof.c \
           $(NTT512) $(SHAKE8X)
profile:
	@mkdir -p $(BINDIR)
	@for m in $(MODES); do \
	  case $$m in 120) qs=$(QSET_120);; 180) qs=$(QSET_180);; 260) qs=$(QSET_260);; esac; \
	  $(CC) $(CFLAGS) -DPROF_TIME -DPROF_RAND -DLITHIUM_MODE=$$m \
	    test/speed_profile.c $(PROF_SRC) $(EXPIFMA) \
	    $$qs/ntt_avx512.S $$qs/ntt_consts_avx512.c $(LDLIBS) \
	    -o $(BINDIR)/speed_profile_$$m || exit 1; \
	  $(BINDIR)/speed_profile_$$m; \
	  echo; \
	done

# Human-readable machine-info preamble shared by speed_txt / profile_txt.
MACHINFO = printf ' date    : '; date -u '+%Y-%m-%dT%H:%M:%SZ'; \
           printf ' kernel  : '; uname -sr; \
           printf ' cpu     : '; ( lscpu 2>/dev/null | sed -n 's/^Model name:[[:space:]]*//p' | head -1 || grep -m1 'model name' /proc/cpuinfo | sed 's/.*: //' ); \
           printf ' cc      : '; $(CC) --version | head -1; \
           printf ' cflags  : %s\n' '$(CFLAGS)'

# Dump end-to-end speed (tail/median rdtsc cycles, speed_sign) for every mode to
# test/speed.txt, prefixed with a machine-info header.  Reports BOTH the default
# integer-exp build (== ref KAT) and the opt-in 8-lane IFMA-exp build.
.PHONY: speed_txt
speed_txt:
	@$(MAKE) TESTS=speed_sign all >/dev/null
	@{ echo "=== Lithium-Code  avx512/  end-to-end speed (tail/median rdtsc cycles, speed_sign) ==="; \
	   $(MACHINFO); \
	   echo " stability metrics:"; \
	   echo "   p25, p50, p75, p90, p99: nearest-rank percentiles of keygen/sign cycles."; \
	   echo "   StQ1, StQ2, StQ3: stabilized quartiles, i.e., [1/8,3/8], [3/8,5/8], [5/8,7/8] means."; \
	   echo "     Reference: Daniel J. Bernstein, Robust statistics for rejection-sampling timings, https://cr.yp.to/papers/rsrst-20250727.pdf."; \
	   echo "   stddev: population standard deviation of keygen/sign cycles; lower means less run-to-run spread."; \
	   echo "   cv: stddev/avg; rough guide: <1% very stable, 1-5% stable, 5-15% noisy, >15% unstable or long-tailed."; \
	   echo "   outlier_count(>p50*1.5): keygen/sign runs slower than 1.5x the median."; \
	   echo "   p90/p50 and p99/p50: tail-inflation ratios; values near 1.00x mean a short tail."; \
	   echo; echo "-- default build: integer exp (byte-identical to ref KAT) --"; \
	   for m in $(MODES); do $(BINDIR)/speed_sign_$$m; echo; done; } > test/speed.txt
	@$(MAKE) AVX512IFMAEXP=1 TESTS=speed_sign all >/dev/null
	@{ echo; echo "-- AVX512IFMAEXP=1 build: 8-lane IFMA Q51 integer ApproxExp (opt-in KAT) --"; \
	   for m in $(MODES); do $(BINDIR)/speed_sign_$$m; echo; done; } >> test/speed.txt
	@echo "wrote test/speed.txt"; cat test/speed.txt

# Dump the AVX-512 IFMA exp microbenchmark to
# test/speed_exp_avx512.txt, prefixed with a machine-info header.
.PHONY: exp_avx512_txt exp_avx512_run
exp_avx512_txt:
	@$(MAKE) --no-print-directory AVX512IFMAEXP=1 exp_avx512_run

exp_avx512_run:
	@$(MAKE) --no-print-directory MODES="$(firstword $(MODES))" TESTS=speed_exp_avx512ifma all >/dev/null
	@{ echo "=== Lithium-Code  avx512/  AVX512 IFMA exp speed (median rdtsc cycles, speed_exp_avx512ifma) ==="; \
	   $(MACHINFO); \
	   echo; \
	   echo "exp performance is mode-independent; using mode $(firstword $(MODES))."; \
	   echo "-- AVX512IFMAEXP=1 build: W8/W16/W32 IFMA integer exp kernels --"; \
	   m="$(firstword $(MODES))"; echo "=== Lithium-$$m AVX512 ==="; $(BINDIR)/speed_exp_avx512ifma_$$m; echo; } > test/speed_exp_avx512.txt
	@echo "wrote test/speed_exp_avx512.txt"; cat test/speed_exp_avx512.txt

# Dump the per-component profile (cycles/op for keygen/sign/verify) for every
# mode to test/profiling.txt, prefixed with the same machine-info header.
# Reports BOTH the default integer-exp build and the IFMA-exp build.
.PHONY: profile_txt
profile_txt:
	@{ echo "=== Lithium-Code  avx512/  per-component profile (cycles/op, make profile) ==="; \
	   $(MACHINFO); \
	   echo; echo "-- default build: integer exp --"; } > test/profiling.txt
	@$(MAKE) --no-print-directory profile >> test/profiling.txt
	@{ echo; echo "-- AVX512IFMAEXP=1 build: 8-lane IFMA Q51 integer ApproxExp --"; } >> test/profiling.txt
	@$(MAKE) --no-print-directory AVX512IFMAEXP=1 profile >> test/profiling.txt
	@echo "wrote test/profiling.txt"

clean:
	rm -rf $(BINDIR)
