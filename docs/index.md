---
layout: full
homepage: true
disable_anchors: true
description: Genetic separation between two disease subtypes from GWAS summary statistics
---
## gensep Overview

**gensep** quantifies the **genetic separation** between two disease subtypes (or two
related traits) directly from GWAS summary statistics. Given a tagging file and two
traits' summary statistics — or simply their heritability/genetic-correlation point
estimates — it reports, on one common SNP set:

- **hsq1 / hsq2** — per-trait SNP heritability, observed and liability scale;
- **rg** — genetic correlation between the two subtypes;
- **VS** — the genetic-separation variance built from the liability-scale heritabilities and rg;
- **h2cc** — the case–case heritability, `VS / (VS + 4)`;
- **auc** — the upper-limit (ceiling) AUC of a genetic classifier separating the two subtypes;
- **auc_lo** — the leading-order AUC approximation.

Every quantity comes with a **standard error**. gensep is a self-contained C++ port of
LDAK SumHer (`--sum-hers` / `--sum-cors`), extended with the joint uncertainty on the
derived separation quantities that LDAK does not expose. You choose how the SE is
computed with a single switch, `--se-method`:

- **`jackknife`** — from summary statistics, via a fused block-jackknife over one common SNP set;
- **`mc`** / **`delta`** — from given point estimates, via Monte-Carlo or first-order (delta) propagation;
- **`none`** — point estimates only, no SE.

gensep is implemented as an open-source, statically linked C++17 program with no runtime
dependencies (Eigen is header-only and vendored).

## User's Guide: [Installation](./documentation/02_installation.html) · [Tutorial](./documentation/03_Tutorial.html) · [Method](./documentation/04_Method.html)

## Citation
gensep builds on the SumHer method (Speed & Balding, *Nature Genetics* 2019) implemented
in [LDAK](https://dougspeed.com/ldak/), and on the liability/AUC theory of Wray et al. and
Lee et al. Manuscript in preparation.

## Contact
For questions, open an issue on [GitHub](https://github.com/chaoning/gensep/issues) or
email me at ningchao91@gmail.com
