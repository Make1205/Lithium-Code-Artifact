CC ?= gcc

SRC_DIR = ..

CFLAGS_BASE = -Wall -Wextra -Wpedantic -O3 -std=c99 \
              -I. -I$(SRC_DIR) -I$(SRC_DIR)/test
LDLIBS = -lm
OUTDIR = $(SRC_DIR)/out/vector

# Same source set as the ref/ test-sign target, minus randombytes.c
# (we link against our deterministic rng.c instead).
SOURCES = $(SRC_DIR)/fips202.c $(SRC_DIR)/symmetric-shake.c \
          $(SRC_DIR)/approx_exp.c \
          $(SRC_DIR)/sampler.c \
          $(SRC_DIR)/irs.c \
          $(SRC_DIR)/reduce.c $(SRC_DIR)/ntt.c \
          $(SRC_DIR)/polyvec.c \
          $(SRC_DIR)/rounding.c $(SRC_DIR)/irs_rans.c \
          $(SRC_DIR)/packing.c $(SRC_DIR)/sign.c

MODES ?= 120 180 260
NROUNDS ?= 10
GENERATORS = $(addprefix $(OUTDIR)/gen_vectors_mode,$(MODES))

.PHONY: all build vectors clean clean-bin clean-vectors

# Default: build all modes and produce the vector files.
all: vectors

build: $(GENERATORS)

$(OUTDIR)/gen_vectors_mode%: gen_vectors.c rng.c rng.h $(SOURCES) | $(OUTDIR)
	$(CC) $(CFLAGS_BASE) -DLITHIUM_MODE=$* \
	    -o $@ gen_vectors.c rng.c $(SOURCES) $(LDLIBS)

$(OUTDIR):
	mkdir -p $(OUTDIR)

# Run each mode's generator.  Output files land directly in this
# directory (ref/vector/) so they're easy to ship alongside the source.
vectors: build
	@for m in $(MODES); do \
	  $(OUTDIR)/gen_vectors_mode$$m . $(NROUNDS) || exit 1; \
	done

clean: clean-bin

clean-bin:
	rm -rf $(OUTDIR)

clean-vectors:
	rm -f $(foreach m,$(MODES),keygen_mode$(m).kat sign_mode$(m).kat verify_mode$(m).kat)
