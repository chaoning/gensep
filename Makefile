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

# OpenMP is ON by default: the block-jackknife loop is parallel, and the thread count is
# chosen at run time with --max-threads (default 1 = single-threaded). The static binary
# stays fully static ("not a dynamic executable"). Disable entirely with `make OMP=0`.
OMP ?= 1
ifeq ($(OMP),1)
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

# Refresh the committed prebuilt binary (app/linux/gensep) after a source change,
# gated on the regression suite passing. Then `git add app/linux/gensep && git commit`.
release: test
	mkdir -p app/linux
	cp -f $(BIN) app/linux/$(BIN)
	@echo "updated app/linux/$(BIN)  ->  commit it:  git add app/linux/$(BIN) && git commit"

clean:
	rm -f $(OBJ) $(DEP) $(BIN)

# Pull in the auto-generated header dependencies (silent if not built yet).
-include $(DEP)

.PHONY: all clean test release
