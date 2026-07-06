# gensep — self-contained C++ port of LDAK SumHer (sum-hers/sum-cors) + block-jackknife
# No external links: Eigen is header-only (vendored under third_party/eigen).
# Use the system GCC, which ships static libc/libm/libstdc++ (the conda toolchain
# does not, so -static fails there). Override on the command line: make CXX=...
CXX      := /usr/bin/g++
EIGEN    ?= third_party/eigen
# Use ':=' (not '?=') so these win over any CXXFLAGS exported by the environment
# (e.g. a conda activate script) — otherwise -I$(EIGEN) would be dropped.
# -MMD -MP emit per-object .d files listing header deps, so editing a .hpp
# (e.g. gensep.hpp) forces every dependent .cpp to recompile — avoids stale
# objects with a mismatched struct layout. -MP adds phony header targets so a
# deleted/renamed header doesn't break the build.
CXXFLAGS := -O3 -std=c++17 -Wall -Wextra -I$(EIGEN) -Isrc -MMD -MP
# Fully static link (no runtime libc/libstdc++/libgcc/libm dependency) so the
# binary runs on any compatible Linux host; -s strips symbols to shrink it.
LDFLAGS  := -static -static-libstdc++ -static-libgcc -s

# Optional OpenMP: make OMP=1
ifdef OMP
CXXFLAGS += -fopenmp
LDFLAGS  += -fopenmp
endif

SRC := $(wildcard src/*.cpp)
OBJ := $(SRC:.cpp=.o)
DEP := $(OBJ:.o=.d)
BIN := gensep

all: $(BIN)

$(BIN): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LDFLAGS)

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

test: $(BIN)
	bash tests/run_tests.sh ./$(BIN)

clean:
	rm -f $(OBJ) $(DEP) $(BIN)

# Pull in the auto-generated header dependencies (silent if not built yet).
-include $(DEP)

.PHONY: all clean test
