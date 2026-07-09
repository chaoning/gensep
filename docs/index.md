---
layout: full
homepage: true
disable_anchors: true
description: Quantifying the genetic separability of disease subtypes
---
## gensep Overview

**gensep** implements a liability-threshold framework for **case–case subtype
discrimination** — how well genetics can tell two disease subtypes apart. The central
quantity is the **genetic separation variance**

```
V_S = λ1² h1² + λ2² h2² − 2 λ1 λ2 rg h1 h2,
```

built from the two subtypes' liability-scale SNP heritabilities (`h1²`, `h2²`), their
genetic correlation (`rg`), and the prevalence-driven selection intensities (`λ_i`). `V_S`
is large when the subtypes have strong subtype-specific genetic components and small when
their genetic effects are largely shared. From `V_S`, gensep derives the **oracle case–case
AUC** (the maximum achievable AUC if the true genetic values were known) and its
leading-order approximation, and the balanced observed-scale **case–case heritability**
`h²_cc = V_S / (V_S + 4)` — a bounded, AUC-linked summary of genetic separation. Every
quantity comes with a **standard error**.

gensep computes these three interchangeable ways:

- **From GWAS summary statistics** — a tagging file and the two subtypes' summaries; gensep
  estimates the heritabilities and genetic correlation (SumHer) on one common SNP set and
  derives everything with a fused **block-jackknife** SE.
- **From point estimates** — observed-scale heritabilities and genetic correlation with
  their SEs, from **LDAK / SumHer, LDSC, or any other method**; the SE is propagated by
  **Monte-Carlo** or the **delta** method.
- **Plus finite-PRS AUC** — given each subtype's PRS case/control AUC, gensep also reports
  the AUC achievable with those polygenic scores and the **PRS recovery** `V_PRS / V_S`.

gensep is an open-source, statically linked C++17 program with no runtime dependencies
(Eigen is header-only and vendored).

## User's Guide: [Installation](./documentation/02_installation.html) · [Tutorial](./documentation/03_Tutorial.html)

## Citation
Chao Ning, Jasper Hof and Doug Speed. *Quantifying the Genetic Separability of Disease
Subtypes* (in preparation). The heritability and genetic-correlation solvers are a port of
SumHer (Speed & Balding, *Nature Genetics* 2019), implemented in
[LDAK](https://dougspeed.com/ldak/).

## Contact
For questions, open an issue on [GitHub](https://github.com/chaoning/gensep/issues) or
email me at ningchao91@gmail.com
