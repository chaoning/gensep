# gensep: **Quantifying the genetic separability of disease subtypes**

## Overview

`gensep` implements a liability-threshold framework for **case–case subtype discrimination**
— how well genetics can tell two disease subtypes apart. The central quantity is the
**genetic separation variance**

```
V_S = λ1² h1² + λ2² h2² − 2 λ1 λ2 rg h1 h2,
```

built from the two subtypes' liability-scale SNP heritabilities (`h1²`, `h2²`), their
genetic correlation (`rg`), and the prevalence-driven selection intensities (`λ_i`). `V_S`
is large when the subtypes have strong subtype-specific genetic components and small when
their genetic effects are largely shared. From `V_S`, `gensep` derives:

- the **oracle case–case AUC** — the maximum achievable AUC for predicting subtype if the
  true genetic values were known — and its leading-order approximation;
- the balanced observed-scale **case–case heritability** `h²_cc = V_S / (V_S + 4)`, a
  bounded, AUC-linked summary of genetic separation.

`gensep` computes these three interchangeable ways, all writing the same output:

1. **From GWAS summary statistics** — supply a tagging file and the two subtypes' summaries;
   `gensep` estimates the subtype heritabilities and genetic correlation (SumHer) on one
   common SNP set and derives everything with a fused block-jackknife standard error.
2. **From point estimates** — if you already have the subtypes' observed-scale
   heritabilities and genetic correlation with their SEs (from **LDAK / SumHer, LDSC, or any
   other method**), `gensep` computes the separation quantities directly, propagating the SE
   by Monte-Carlo or the delta method.
3. **Plus finite-PRS AUC** — additionally supply each subtype's PRS case/control AUC and
   `gensep` reports the AUC achievable with those finite-accuracy polygenic scores, together
   with **PRS recovery** (`V_PRS / V_S`, the fraction of the oracle genetic separation the
   PRS captures).

Every quantity comes with a standard error. `gensep` is a self-contained C++ tool
accompanying *Quantifying the Genetic Separability of Disease Subtypes* (Ning, Hof & Speed);
the heritability and genetic-correlation solvers are a port of LDAK SumHer (`--sum-hers` /
`--sum-cors`).

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
