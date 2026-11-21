# Build everything with one command:  make
# -ffp-contract=off is required by the project (no fused multiply-adds),
# so results are bit-for-bit reproducible across machines.
MPICXX   ?= mpicxx
# OMPI_SKIP_MPICXX: skip Open MPI's deprecated C++ bindings (we use the C API).
CXXFLAGS ?= -O3 -std=c++17 -ffp-contract=off -Wall -Wextra -DOMPI_SKIP_MPICXX
BUILD    := build

PROGRAMS := $(BUILD)/hello $(BUILD)/pingpong

all: $(PROGRAMS)

$(BUILD)/%: src/%.cpp | $(BUILD)
	$(MPICXX) $(CXXFLAGS) -o $@ $<

$(BUILD):
	mkdir -p $@

clean:
	rm -rf $(BUILD)

.PHONY: all clean
