# Build everything with one command:  make
# -ffp-contract=off is required by the project (no fused multiply-adds),
# so results are bit-for-bit reproducible across machines.
MPICXX   ?= mpicxx
# OMPI_SKIP_MPICXX: skip Open MPI's deprecated C++ bindings (we use the C API).
CXXFLAGS ?= -O3 -std=c++17 -ffp-contract=off -Wall -Wextra -DOMPI_SKIP_MPICXX
CPPFLAGS += -Isrc
BUILD    := build

# The solver is every .cpp under src/ except the two stand-alone programs.
HEAT_SRC := $(filter-out src/hello.cpp src/pingpong.cpp,$(shell find src -name '*.cpp'))
HEAT_OBJ := $(HEAT_SRC:src/%.cpp=$(BUILD)/obj/%.o)

all: $(BUILD)/heat $(BUILD)/hello $(BUILD)/pingpong

$(BUILD)/heat: $(HEAT_OBJ)
	$(MPICXX) $(CXXFLAGS) -o $@ $^

$(BUILD)/obj/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(MPICXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(BUILD)/hello $(BUILD)/pingpong: $(BUILD)/%: src/%.cpp
	@mkdir -p $(BUILD)
	$(MPICXX) $(CXXFLAGS) -o $@ $<

# Rebuild objects when a header they include changes.
-include $(HEAT_OBJ:.o=.d)

clean:
	rm -rf $(BUILD)

.PHONY: all clean
