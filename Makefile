CXX       = mpicxx
CXXFLAGS ?= -O2 -Wall -Wextra -std=c++17 -Isrc
SRC      := $(shell find src -name '*.cpp')
OBJ      := $(SRC:src/%.cpp=build/%.o)

heat: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

build/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf build heat
.PHONY: clean