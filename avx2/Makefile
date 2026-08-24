# Makefile for the Lithium signature AVX2 implementation.
#
#   make                 # build all tests for all modes -> out/
#   make check           # build + run every test for every mode
#   make check-kat       # build + run only test_kat; assert the recorded hashes
#   make TESTS="t_poly"  # build a subset of tests
#   make AVX2EXP=1 ...    # opt-in float exp (separate KAT; default OFF == ref KAT)
#   make clean
#
# Shared algorithm/KAT source is symlinked from ../ref (single source of truth).
# AVX2-specific overrides are real files here; the per-mode NTT assembly +
# constants are VENDORED locally in ./<qset>/ (copied from tools/ntt; the build
# is self-contained, `check-consts` proves the copies match the generator).
#
# Compile-time switches (set in CFLAGS below) select the AVX2 code paths in the
# symlinked ref sources:
#   USE_AVX2_NTT     -> poly_ntt/invntt/pointwise routes in polyvec.c
#   USE_AVX2_SAMPLER -> 8-way sigma2/noise sampler routes in sampler.c
#   USE_AVX2_SHAKE4X -> keccak4x sampler/uniform routes in sampler.c/polyvec.c

MODES ?= 120 180 260
TESTS ?= test_ntt t_poly t_exp speed_approx_exp test_sampler test_irs \
         test_mod2q test_rans test_keys test_sign test_kat test_ntt_avx \
         test_freeze_avx test_sigma2_diff
CC     ?= gcc
AVX2EXP ?= 0
# -mno-avx512f forces a PURE AVX2 path: -march=native would otherwise enable
# AVX-512 on AVX-512-capable systems and let GCC auto-vectorize the symlinked scalar C
# into zmm ops.  -mno-avx512f transitively disables every AVX-512 subset (they
# all depend on F) while keeping AVX2/FMA/BMI2.  Enforced by `make check-no-avx512`.
CFLAGS ?= -O3 -fstack-protector-strong -Wall -Wextra -std=c99 -I. -I../ref -I../ref/test \
          -mavx2 -mbmi2 -mpopcnt -march=native -mtune=native -mno-avx512f \
          -DUSE_AVX2_NTT -DUSE_AVX2_SAMPLER -DUSE_AVX2_SHAKE4X
LDLIBS ?= -lm
BINDIR  = out
CT_OPT ?= -O3

ifeq ($(AVX2EXP),1)
# Float exp: deterministic across compilers (explicit FMA intrinsics + no
# contraction/fast-math).  Separate, opt-in KAT.
CFLAGS += -DAVX2_EXP -ffp-contract=off -fno-fast-math
EXPAVX  = approx_exp_avx2.c
else
EXPAVX  =
endif

# Per-test source lists (POSIX-sh strings).  Tests that list polyvec.c also get
# the per-mode vendored ntt.S + ntt_consts.c appended in the build loop.  The
# AVX2 NTT routing + 4-way uniform now live in polyvec.c, and the 8-way sigma2 /
# noise + keccak4x samplers in sampler.c (so both pull fips202x4.c + f1600x4.S).
SH_SRCS = \
'test_ntt=ntt.c reduce.c' \
't_poly=polyvec.c ntt.c reduce.c fips202.c fips202x4.c f1600x4.S' \
't_exp=approx_exp.c $(EXPAVX)' \
'speed_approx_exp=approx_exp.c test/cpucycles.c $(EXPAVX)' \
'test_sampler=sampler.c approx_exp.c fips202.c symmetric-shake.c fips202x4.c f1600x4.S $(EXPAVX)' \
'test_irs=irs.c polyvec.c ntt.c reduce.c approx_exp.c fips202.c symmetric-shake.c fips202x4.c f1600x4.S $(EXPAVX)' \
'test_mod2q=rounding.c reduce.c' \
'test_rans=irs_rans.c' \
'test_keys=polyvec.c ntt.c reduce.c sampler.c approx_exp.c fips202.c symmetric-shake.c randombytes.c fips202x4.c f1600x4.S $(EXPAVX)' \
'test_sign=ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c randombytes.c fips202x4.c f1600x4.S $(EXPAVX)' \
'test_kat=ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c fips202x4.c f1600x4.S $(EXPAVX)' \
'parser_negative=ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c fips202x4.c f1600x4.S $(EXPAVX)' \
'security_fault_injection=ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c fips202x4.c f1600x4.S $(EXPAVX)' \
'trace_kat=ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c fips202x4.c f1600x4.S $(EXPAVX)' \
'diff_fuzz=ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c fips202x4.c f1600x4.S $(EXPAVX)' \
'dudect_components=polyvec.c ntt.c reduce.c rounding.c test/cpucycles.c fips202.c fips202x4.c f1600x4.S' \
'test_ntt_avx=polyvec.c ntt.c reduce.c fips202.c fips202x4.c f1600x4.S' \
'test_freeze_avx=polyvec.c ntt.c reduce.c fips202.c fips202x4.c f1600x4.S' \
'test_sigma2_diff=sampler.c approx_exp.c fips202.c symmetric-shake.c fips202x4.c f1600x4.S $(EXPAVX)' \
'speed_exp_avx2=approx_exp_avx2.c test/cpucycles.c' \
'speed_sign=ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c randombytes.c test/cpucycles.c fips202x4.c f1600x4.S $(EXPAVX)'

QSET_120 = q57601n256
QSET_180 = q30977n256
QSET_260 = q64513n512
# The per-mode NTT (ntt.S + ntt_consts.c) is VENDORED locally in ./<qset>/ so
# this build is self-contained.  NTTGEN points at the generator (tools/ntt),
# used only by `check-consts` to prove the vendored copies are reproducible.
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
	      case "$$srcs" in *polyvec.c*) extra="$$qs/ntt.S $$qs/ntt_consts.c";; esac; \
	      stats=""; case "$$t" in trace_kat) stats="-DTEST_RANS_STATS";; diff_fuzz) stats="-DDIFF_FUZZ_N=$(DIFF_FUZZ_N)";; dudect_components) stats="-DDUDECT_COMPONENT_N=$(DUDECT_COMPONENT_N)";; esac; \
	      $(CC) $(CFLAGS) $$stats -DLITHIUM_MODE=$$m test/$$t.c $$srcs $$extra $(LDLIBS) \
	        -o $(BINDIR)/$${t}_$${m} || exit 1; \
	    done; \
	  done; \
	done
	@echo "built: [$(TESTS)] x [$(MODES)]"

check: all
	@for m in $(MODES); do \
	  echo "=== Lithium-$$m AVX2 ==="; \
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
                    randombytes.c fips202x4.c f1600x4.S $(EXPAVX)
.PHONY: check-rans-failrate
check-rans-failrate:
	@mkdir -p $(BINDIR)
	@for m in $(MODES); do \
	  case $$m in 120) qs=$(QSET_120);; 180) qs=$(QSET_180);; 260) qs=$(QSET_260);; esac; \
	  $(CC) $(CFLAGS) -DTEST_RANS_STATS -DRANS_STATS_N=$(RANS_STATS_N) \
	    -DLITHIUM_MODE=$$m ../ref/test/test_rans_failrate.c $(RANS_FAILRATE_SRC) \
	    $$qs/ntt.S $$qs/ntt_consts.c $(LDLIBS) \
	    -o $(BINDIR)/test_rans_failrate_$$m || exit 1; \
	  $(BINDIR)/test_rans_failrate_$$m || exit 1; \
	done

DIFF_FUZZ_N ?= 32
DUDECT_COMPONENT_N ?= 512

DUDECT_N ?= 2000
DUDECT_SIGN_SRC = ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c \
                  irs.c packing.c irs_rans.c sign.c fips202.c symmetric-shake.c \
                  fips202x4.c f1600x4.S test/cpucycles.c $(EXPAVX)
.PHONY: dudect-sign
dudect-sign:
	@mkdir -p $(BINDIR)
	@for m in $(MODES); do \
	  case $$m in 120) qs=$(QSET_120);; 180) qs=$(QSET_180);; 260) qs=$(QSET_260);; esac; \
	  $(CC) $(CFLAGS) -DDUDECT_N=$(DUDECT_N) -DLITHIUM_MODE=$$m \
	    test/dudect_sign.c $(DUDECT_SIGN_SRC) $$qs/ntt.S $$qs/ntt_consts.c \
	    $(LDLIBS) -o $(BINDIR)/dudect_sign_$$m || exit 1; \
	  $(BINDIR)/dudect_sign_$$m || exit 1; \
	done

# Regression-lock the opt-in float-exp KAT against the recorded KAT_avx2exp.txt.
.PHONY: check-kat-avx2exp
check-kat-avx2exp:
	@$(MAKE) AVX2EXP=1 TESTS=test_kat all >/dev/null
	@fail=0; for m in $(MODES); do \
	  h=$$($(BINDIR)/test_kat_$$m | grep -o 'KAT_hash=[0-9]*'); \
	  grep -q "mode=$$m $$h" KAT_avx2exp.txt \
	    && echo "  PASS: mode=$$m $$h matches KAT_avx2exp.txt" \
	    || { echo "  FAIL: mode=$$m $$h NOT in KAT_avx2exp.txt"; fail=1; }; \
	done; [ $$fail -eq 0 ] && echo "check-kat-avx2exp: float-exp KAT locked" || exit 1

# Force-exercise the 4-way SHAKE continuation-squeeze (refill) path with a tiny
# initial squeeze, and assert the KAT is unchanged -- proves the rarely-hit
# refill code is correct (the production INIT sizes usually avoid it).
.PHONY: check-refill
check-refill:
	@$(MAKE) CFLAGS="$(CFLAGS) -DGAUSS4X_INIT=4 -DNOISE4X_INIT=4" TESTS="test_kat test_sampler" all >/dev/null
	@echo "== forced 4-way refill (GAUSS4X_INIT=4 NOISE4X_INIT=4): KAT must be unchanged =="
	@for m in $(MODES); do $(BINDIR)/test_kat_$$m; $(BINDIR)/test_sampler_$$m | grep -o '4x==1x=[01]'; done

# Reproducible-constants gate: every generated/derived constant must regenerate
# byte-for-byte and (for the float exp) meet its accuracy bound.  Run before
# `make check`.  Uses cmp against a pre-regeneration copy (robust for both
# tracked and freshly-generated headers).
.PHONY: check-consts
check-consts:
	@echo "== float exp consts (@@AUTOGEN:approx_exp_avx2@@ in approx_exp_avx2.h; tools/avxexp/gen_avxexp_consts.py) =="
	@mkdir -p $(BINDIR)
	@cp approx_exp_avx2.h $(BINDIR)/consts.bak
	@python3 ../tools/avxexp/gen_avxexp_consts.py --check
	@cmp -s approx_exp_avx2.h $(BINDIR)/consts.bak \
	  && echo "  PASS: approx_exp_avx2.h consts region reproducible" \
	  || { echo "  FAIL: approx_exp_avx2.h drift"; rm -f $(BINDIR)/consts.bak; exit 1; }
	@rm -f $(BINDIR)/consts.bak
	@echo "== vendored NTT (./<qset>/ntt.S + ntt_consts.c) vs generator (tools/ntt/<qset>/gen_ntt.py) =="
	@for qs in $(QSET_120) $(QSET_180) $(QSET_260); do \
	  ( cd $(NTTGEN)/$$qs && python3 gen_ntt.py >/dev/null 2>&1 ) || { echo "  FAIL: $$qs gen_ntt.py"; exit 1; }; \
	  cmp -s $(NTTGEN)/$$qs/ntt_consts.c $$qs/ntt_consts.c && cmp -s $(NTTGEN)/$$qs/ntt.S $$qs/ntt.S \
	    && echo "  PASS: $$qs vendored ntt_consts.c + ntt.S match gen_ntt.py" \
	    || { echo "  FAIL: $$qs vendored NTT differs from gen_ntt.py (re-copy)"; exit 1; }; \
	done
	@echo "== avx-vs-ref NTT transpose-equivalence (tools/gen_ntt_perm.py guard) =="
	@( cd ../tools && python3 gen_ntt_perm.py >/dev/null 2>&1 ) \
	  && echo "  PASS: avx NTT transpose-equivalent to ntt16 (expand_A0_ntt nttunpack ladder valid)" \
	  || { echo "  FAIL: avx NTT not transpose-equivalent to ntt16 (avx generator root must match gen_zetas.py)"; exit 1; }
	@echo "== RCDT tables (tools/gen_rcdt.py; merged into sampler.c @@AUTOGEN:rcdt@@, ref+avx2) =="
	@mkdir -p $(BINDIR)
	@cp ../ref/sampler.c $(BINDIR)/rcdt_ref.bak; cp sampler.c $(BINDIR)/rcdt_avx2.bak
	@( cd ../tools && python3 gen_rcdt.py >/dev/null 2>&1 ) || { echo "  FAIL: gen_rcdt.py"; rm -f $(BINDIR)/rcdt_ref.bak $(BINDIR)/rcdt_avx2.bak; exit 1; }
	@cmp -s ../ref/sampler.c $(BINDIR)/rcdt_ref.bak && cmp -s sampler.c $(BINDIR)/rcdt_avx2.bak \
	  && echo "  PASS: sampler.c RCDT region reproducible (ref + avx2)" \
	  || { echo "  FAIL: sampler.c RCDT drift"; rm -f $(BINDIR)/rcdt_ref.bak $(BINDIR)/rcdt_avx2.bak; exit 1; }
	@rm -f $(BINDIR)/rcdt_ref.bak $(BINDIR)/rcdt_avx2.bak
	@echo "== 4-way SHAKE round constants (reuse ref KeccakF_RoundConstants) =="
	@$(MAKE) TESTS=test_kat all >/dev/null
	@n=$$(nm $(BINDIR)/test_kat_120 2>/dev/null | grep -c 'KeccakF_RoundConstants'); \
	  [ "$$n" = "1" ] && echo "  PASS: exactly one KeccakF_RoundConstants symbol (no x4 duplicate)" \
	  || { echo "  FAIL: $$n KeccakF_RoundConstants symbols"; exit 1; }
	@echo "check-consts: all constants reproducible + accuracy gate met"

# Constant-time scan: delegate to the shared object-code scanner.  The root
# ./ct_scan_matrix.sh command covers the full gcc/clang, -O3/-Os, all-mode
# matrix including optional exp variants; this local target follows the current
# Makefile CC/MODES/CT_OPT knobs for quick backend checks.
.PHONY: ct-scan
ct-scan:
	@python3 ../tools/ct_scan.py --backends avx2 --modes "$(MODES)" --ccs "$(CC)" --opts="$(CT_OPT)" --include-optional --report avx2/test/ct_scan_avx2.txt

# Pure-AVX2 gate: assert NO AVX-512 (zmm / EVEX mask / AVX-512-only mnemonic)
# leaked into any built object.  -mno-avx512f should make this impossible; the
# gate proves the flag is effective (and would catch its accidental removal,
# which -march=native on the recorded benchmark system would otherwise let auto-vectorize to zmm).
.PHONY: check-no-avx512
check-no-avx512:
	@rm -rf $(BINDIR)            # avoid scanning stale (e.g. -march=native) binaries
	@$(MAKE) all >/dev/null
	@fail=0; \
	for b in $(BINDIR)/*; do \
	  [ -f "$$b" ] || continue; \
	  hits=$$(objdump -d "$$b" 2>/dev/null | grep -iE '%zmm[0-9]|\{%k[0-7]\}|vpternlog|vpcompress|vpexpand|vpmadd52|vperm[ti]2[a-z]|kmov|kshift|kortest|kunpck|valign[dq]' || true); \
	  if [ -n "$$hits" ]; then echo "  FAIL $$b:"; echo "$$hits" | head -3; fail=1; fi; \
	done; \
	[ $$fail -eq 0 ] && echo "check-no-avx512: pure AVX2 (no zmm/EVEX/AVX-512 op in any build object)" || { echo "check-no-avx512: AVX-512 LEAKED"; exit 1; }

# Per-component time + randomness profile of crypto_sign (test/speed_profile.c).
# `make profile`            -> AVX2 default (integer exp)
# `make AVX2EXP=1 profile`  -> AVX2 float exp
.PHONY: profile
PROF_SRC = ntt.c reduce.c polyvec.c rounding.c sampler.c approx_exp.c irs.c \
           packing.c irs_rans.c sign.c fips202.c symmetric-shake.c randombytes.c ../ref/test/prof.c \
           fips202x4.c f1600x4.S
profile:
	@mkdir -p $(BINDIR)
	@for m in $(MODES); do \
	  case $$m in 120) qs=$(QSET_120);; 180) qs=$(QSET_180);; 260) qs=$(QSET_260);; esac; \
	  $(CC) $(CFLAGS) -DPROF_TIME -DPROF_RAND -DLITHIUM_MODE=$$m \
	    test/speed_profile.c $(PROF_SRC) $(EXPAVX) \
	    $$qs/ntt.S $$qs/ntt_consts.c $(LDLIBS) \
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
# integer-exp build (== ref KAT) and the opt-in AVX2Exp 4-wide binary64 build.
.PHONY: speed_txt
speed_txt:
	@$(MAKE) TESTS=speed_sign all >/dev/null
	@{ echo "=== Lithium-Code  avx2/  end-to-end speed (tail/median rdtsc cycles, speed_sign) ==="; \
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
	@$(MAKE) AVX2EXP=1 TESTS=speed_sign all >/dev/null
	@{ echo; echo "-- AVX2EXP=1 build: 4-wide binary64 ApproxExp (opt-in KAT) --"; \
	   for m in $(MODES); do $(BINDIR)/speed_sign_$$m; echo; done; } >> test/speed.txt
	@echo "wrote test/speed.txt"; cat test/speed.txt

# Dump the AVX2 float exp microbenchmark to
# test/speed_exp_avx2.txt, prefixed with a machine-info header.
.PHONY: exp_avx2_txt exp_avx2_run
exp_avx2_txt:
	@$(MAKE) --no-print-directory AVX2EXP=1 exp_avx2_run

exp_avx2_run:
	@$(MAKE) --no-print-directory MODES="$(firstword $(MODES))" TESTS=speed_exp_avx2 all >/dev/null
	@{ echo "=== Lithium-Code  avx2/  AVX2 exp speed (median rdtsc cycles, speed_exp_avx2) ==="; \
	   $(MACHINFO); \
	   echo; \
	   echo "exp performance is mode-independent; using mode $(firstword $(MODES))."; \
	   echo "-- AVX2EXP=1 build: W4/W8/W16 binary64 exp kernels --"; \
	   m="$(firstword $(MODES))"; echo "=== Lithium-$$m AVX2 ==="; $(BINDIR)/speed_exp_avx2_$$m; echo; } > test/speed_exp_avx2.txt
	@echo "wrote test/speed_exp_avx2.txt"; cat test/speed_exp_avx2.txt

# Dump the per-component profile (cycles/op for keygen/sign/verify) for every
# mode to test/profiling.txt, prefixed with the same machine-info header.
# Reports BOTH the default integer-exp build and the AVX2Exp float build.
.PHONY: profile_txt
profile_txt:
	@{ echo "=== Lithium-Code  avx2/  per-component profile (cycles/op, make profile) ==="; \
	   $(MACHINFO); \
	   echo; echo "-- default build: integer exp --"; } > test/profiling.txt
	@$(MAKE) --no-print-directory profile >> test/profiling.txt
	@{ echo; echo "-- AVX2EXP=1 build: 4-wide binary64 ApproxExp --"; } >> test/profiling.txt
	@$(MAKE) --no-print-directory AVX2EXP=1 profile >> test/profiling.txt
	@echo "wrote test/profiling.txt"

clean:
	rm -rf $(BINDIR)
