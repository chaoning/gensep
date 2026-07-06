---
layout: page
title: Installation
description: ~
order: 2
---

`gensep` is a single self-contained C++17 command-line tool. There is no external library
to install: [Eigen](https://eigen.tuxfamily.org) is header-only and vendored under
`third_party/`, and the program links statically. Building is a plain `make`.

## 1. Prerequisites

- A C++17 compiler that can produce a static binary. The Makefile pins the **system**
  GCC (`/usr/bin/g++`) because it ships static `libc` / `libm` / `libstdc++`; the conda
  toolchain usually does **not**, so `-static` fails there.
- GNU `make`.
- (Optional) `bash`, to run the regression test suite.

No MKL, no BLAS, no Boost — the numerics use vendored header-only Eigen.

## 2. Build from source

```bash
cd code/gensep
make            # -> ./gensep
```

`make` compiles everything under `src/` and links a single executable, `./gensep`.
The binary is **fully static** (`ldd ./gensep` → "not a dynamic executable", ~1.2 MB
stripped), so it runs on any compatible Linux host with no runtime dependencies.

### Overriding the compiler

If the pinned `/usr/bin/g++` is not what you want, override `CXX` on the command line:

```bash
make CXX=/path/to/your/g++
```

Note that the compiler you choose must be able to satisfy `-static -static-libstdc++
-static-libgcc`; otherwise the final link step will fail.

### Optional OpenMP

An OpenMP build is available if you want threaded numerics:

```bash
make OMP=1
```

### Eigen location

Eigen is expected at `third_party/eigen`. In this repository that is a symlink to the
Eigen copy vendored with `fastgxe`. If you relocate the project, point `EIGEN` at your
Eigen headers:

```bash
make EIGEN=/path/to/eigen
```

The Makefile deliberately uses `CXXFLAGS :=` (not `?=`) so that a `CXXFLAGS` exported by
a conda `activate` script cannot silently drop `-Ithird_party/eigen`.

## 3. Quick validation

Confirm the executable starts and prints its usage:

```bash
./gensep --help
```

You should see the `--se-method <jackknife|mc|delta|none>` synopsis.

## 4. Run the regression tests

A minimal, data-free regression suite is included. It exercises the three point-mode SE
methods, point-value correctness, Monte-Carlo determinism, and every input-validation
guard (bad numbers, negative SEs, `--num-draws 0`, malformed tagfile, …):

```bash
make test
```

All checks should report `PASS` and the run ends with `N passed, 0 failed`.

## 5. Troubleshooting

- **`-static` link errors** (e.g. "cannot find `-lc`"): your compiler lacks the static
  runtime libraries. Use the system GCC (`make CXX=/usr/bin/g++`, the default) rather than
  a conda toolchain, or install the static `glibc`/`libstdc++` packages.
- **`-Ithird_party/eigen` seems ignored / Eigen headers not found**: check that
  `third_party/eigen` exists (or pass `make EIGEN=...`), and make sure no environment
  `CXXFLAGS` is interfering.
- **Rebuild from scratch**: `make clean && make`.
