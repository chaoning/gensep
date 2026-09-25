---
title: Home
layout: home
nav_order: 1
description: Quantifying the genetic separability of disease subtypes
permalink: /
---
# GenSep
{: .fs-9 }

Quantifying the genetic separability of disease subtypes.
{: .fs-6 .fw-300 }

[Download GenSep (Linux)](https://github.com/chaoning/gensep/raw/refs/heads/main/app/linux/gensep){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 }
[Installation](documentation/02_installation.html){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }
[Tutorial](documentation/03_Tutorial.html){: .btn .fs-5 .mb-4 .mb-md-0 .mr-2 }
[View on GitHub](https://github.com/chaoning/gensep){: .btn .fs-5 .mb-4 .mb-md-0 }

---

## Overview

**GenSep** implements a liability-threshold framework for **case–case subtype
discrimination** — how well genetics can tell two disease subtypes apart. The central
quantity is the **genetic separation variance**

```
V_S = λ1² h1² + λ2² h2² − 2 λ1 λ2 rg h1 h2,
```

built from the two subtypes' liability-scale SNP heritabilities (`h1²`, `h2²`), their
genetic correlation (`rg`), and the prevalence-driven selection intensities (`λ_i`). `V_S`
is large when the subtypes have strong subtype-specific genetic components and small when
their genetic effects are largely shared. From `V_S`, GenSep derives the **oracle case–case
AUC** (the maximum achievable AUC if the true genetic values were known) and its
leading-order approximation, and the balanced observed-scale **case–case heritability**
`h²_cc = V_S / (V_S + 4)` — a bounded, AUC-linked summary of genetic separation. Every
quantity comes with a **standard error**.

GenSep computes these three interchangeable ways:

- **From GWAS summary statistics** — a tagging file and the two subtypes' summaries; GenSep
  estimates the heritabilities and genetic correlation (SumHer) on one common SNP set and
  derives everything with a fused **block-jackknife** SE.
- **From point estimates** — observed-scale heritabilities and genetic correlation with
  their SEs, from **LDAK / SumHer, LDSC, or any other method**; the SE is propagated by
  **Monte-Carlo** or the **delta** method.
- **Plus finite-PRS AUC** — given each subtype's PRS case/control AUC, GenSep also reports
  the AUC achievable with those polygenic scores and the **PRS recovery** `V_PRS / V_S`.

GenSep is an open-source, statically linked C++17 program with no runtime dependencies
(Eigen is header-only and vendored). GenSep is developed and tested on 64-bit Linux, for
which a prebuilt executable is provided; on Windows it runs under WSL, and on macOS a Linux
container or VM is the supported route.

## Citation
Chao Ning, Jasper Hof and Doug Speed. Quantifying the genetic separability of disease
subtypes. *medRxiv* (2026).
[doi:10.64898/2026.09.18.26363395](https://doi.org/10.64898/2026.09.18.26363395)

The heritability and genetic-correlation solvers are a C++ re-implementation of
SumHer (Speed & Balding, *Nature Genetics* 2019), originally implemented in
[LDAK](https://dougspeed.com/ldak/).

## Contact
For questions, open an issue on [GitHub](https://github.com/chaoning/gensep/issues) or
email me at chao.ning@qgg.au.dk

For other tools, see [chaoning.github.io/software.html](https://chaoning.github.io/software.html).
