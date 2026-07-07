# gensep: **Genetic separation between disease subtypes from GWAS summary statistics**

## Overview

`gensep` quantifies the **genetic separation** between two disease subtypes (or two related
traits) from GWAS summary statistics. Given a tagging file and the two traits' summaries —
or simply their heritability / genetic-correlation point estimates — it reports, on one
common SNP set:

- per-subtype SNP heritability, observed and liability scale (`hsq1` / `hsq2`);
- the genetic correlation `rg`;
- the genetic-separation variance `VS` and the case–case heritability `h2cc = VS/(VS+4)`;
- the upper-limit (ceiling) case–case `auc`, and — given per-subtype PRS AUCs — the
  achievable **finite-PRS** case–case AUC.

Every quantity comes with a **standard error**. `gensep` is a self-contained C++ port of
LDAK SumHer (`--sum-hers` / `--sum-cors`), extended with the joint uncertainty on the
derived separation quantities that LDAK does not expose — via a fused block-jackknife (from
summary statistics) or Monte-Carlo / delta propagation (from point estimates).

## Installation

`gensep` is a single self-contained C++17 program. There is no external library to
install — [Eigen](https://eigen.tuxfamily.org) is header-only and vendored — and it links
statically.

### Linux Executable

A statically compiled executable for 64-bit Linux systems is available:
[**gensep Linux Executable**](https://github.com/chaoning/gensep/raw/refs/heads/main/app/linux/gensep).
It has no runtime dependencies (`ldd` → "not a dynamic executable") and can be used
directly on compatible systems:

```bash
chmod +x gensep
./gensep --help
```

To build from source instead, follow the steps below.

### Prerequisites

- **C++ compiler** with C++17 support. The system GCC is recommended, because it ships the
  static `libc` / `libm` / `libstdc++` needed for a fully static binary (the conda
  toolchain usually does not).
- **GNU make**.

### Installation Steps

#### 1. Clone the Repository

```bash
git clone https://github.com/chaoning/gensep.git
cd gensep
```

#### 2. Provide Eigen

`gensep` expects Eigen headers at `third_party/eigen`. Point it at your Eigen copy:

```bash
ln -s /path/to/eigen third_party/eigen      # or build with: make EIGEN=/path/to/eigen
```

#### 3. Build the Project

```bash
make            # -> ./gensep
```

The binary is fully static (`ldd ./gensep` → "not a dynamic executable", ~1.2 MB
stripped), so it runs on any compatible Linux host with no runtime dependencies.

#### 4. Run gensep

```bash
./gensep --help
```

### Notes

- The Makefile pins `CXX := /usr/bin/g++` for static linking; override with `make CXX=...`.
- Optional OpenMP build: `make OMP=1`.
- Run the regression suite with `make test`.
- Rebuild from scratch with `make clean && make`.

## How to use `gensep`

See the full tutorial in our [documentation](https://chaoning.github.io/gensep).

## Citing the work

`gensep` builds on the SumHer method (Speed & Balding, *Nature Genetics* 2019) implemented
in [LDAK](https://dougspeed.com/ldak/), and on the liability / AUC theory of Wray et al. and
Lee et al. A `gensep` manuscript is in preparation.

## License

`gensep` is free software released under the [GNU General Public License v3.0 (GPL-3.0)](LICENSE).
