---
layout: page
title: Tutorial
description:
order: 3
---

This tutorial follows the current command-line workflow implemented in the codebase. The
executable is assumed to be on your `PATH` as `gensep`; if you built from source and it is
not, replace `gensep` with the full path, e.g. `./gensep`.

GenSep is one command. A single required switch, **`--se-method`**, selects both the
input you provide and how the standard error is computed:

| `--se-method` | Input | SE for VS / h2cc / auc |
| --- | --- | --- |
| `jackknife` | summary statistics + tagging file | fused block-jackknife |
| `mc` | point estimates + their SEs | Monte-Carlo propagation |
| `delta` | point estimates + their SEs | first-order (delta) propagation |
| `none` | point estimates only | none (SE column is `NA`) |

`--se-method` has **no default** — you must pass one of the four values.

All four routes write the same `PREFIX.gensep` output (see [Output](#3-output-file)).

## 1. Run from summary statistics (`jackknife`)

This is the full pipeline: estimate per-trait heritability (SumHer), the genetic
correlation (sum-cors), and the derived separation quantities on one common SNP set, with
the SE from a fused block-jackknife.

```bash
gensep --se-method jackknife \
       --tagfile HumDef.tagging \
       --summary trait1.summaries \
       --summary2 trait2.summaries \
       --K1 0.01 --K2 0.02 \
       --P1 0.3  --P2 0.45 \
       --num-blocks 200 \
       --out output/pair
```

### Inputs

- **`--tagfile`** — an LDAK tagging file. Use a **ready-made** one (see
  [Ready-made tagging files](#ready-made-tagging-files) below), or build your own with
  `ldak --calc-tagging` ([Calculate Taggings](https://dougspeed.com/calculate-taggings/)).
- **`--summary` / `--summary2`** — the two traits' summary statistics, in LDAK
  `.summaries` format: a header line followed by rows of

  ```
  Predictor  A1  A2  Z  n
  ```

  where `Predictor` is the SNP ID (matched against the tagging file), `A1`/`A2` are the
  alleles, `Z` is the signed Z-score, and `n` is the per-SNP sample size. Non-numeric or
  malformed `Z`/`n` fields are rejected rather than silently coerced.
- **`--K1` / `--K2`** — population **prevalences** of subtype 1 and subtype 2 (each in
  `(0, 1)`): the fraction of the general population affected by that subtype. These drive
  the selection intensity λ and the Lee observed→liability transform. `K1` and `K2` are
  independent quantities and need not sum to 1 (or to the prevalence of the parent
  disease).
- **`--P1` / `--P2`** — the **case fractions in each GWAS sample** (each in `(0, 1)`):
  `P1 = N_case / (N_case + N_control)` for the GWAS of subtype 1 (`--summary`), and
  likewise `P2` for the GWAS of subtype 2 (`--summary2`). Each GWAS has its own case and
  control set, so `P1` and `P2` are unrelated and need not sum to 1; e.g. two case–control
  GWAS with 3,000 and 12,000 cases against the same 300,000 controls have `P1 ≈ 0.0099`
  and `P2 ≈ 0.0385`. These are not stored in the summaries, so you must supply them;
  they enter the Lee factor.
- **`--num-blocks`** — number of jackknife blocks (default `200`, must be ≥ 2).
- **`--cutoff`** — exclude strong-effect loci: drop any SNP explaining ≥ `cutoff` of
  phenotypic variance (`rho² = chi/(chi+n)`) in **either** trait, since such loci can bias
  SumHer h²/rg (as LDAK `--sum-cors --cutoff`). Off by default; must be in `(0, 0.5)`,
  e.g. `--cutoff 0.01`. If not set and some SNP exceeds 1%, GenSep prints a reminder. Note
  that with well-powered GWAS (large `n`) per-SNP variance explained is small, so `0.01`
  often excludes nothing.
- **`--intercept`** — `YES` or `NO` (default `NO`). `NO`: each SNP heritability is fit with
  a **fixed intercept of 1** (standard SumHer, assumes no confounding inflation). `YES`: fit
  a **free (LDSC-style) intercept** to absorb inflation from stratification / cryptic
  relatedness — the heritabilities then come from the sum-cors per-trait fit, so `h1`, `h2`
  and `rg` share one model, scale and block set. Use `YES` if your summaries may be inflated
  and not genomic-control corrected. (The genetic correlation always uses a free intercept
  and a sample-overlap term regardless of this flag.)
- **`--max-threads`** — threads for the block-jackknife loop (default `1` = single-threaded).
  e.g. `--max-threads 8` parallelizes the 200 block solves across 8 cores; results are
  identical to single-threaded.
- **`--out`** — output prefix; results are written to `PREFIX.gensep`.

### Ready-made tagging files

You do not need to build a tagging file yourself. We provide pre-computed **HapMap3
non-ambiguous SNP** tagging files, one per ancestry — courtesy of Doug Speed. They are
single-category, built with the **LDAK heritability model** (`ldak --calc-tagging --power
-0.25`) on an ancestry-matched HapMap3 reference panel. Pick the one matching your GWAS
ancestry; the method is robust to SNP overlap, so an exact SNP match is not needed.

| Ancestry | Download (gzip) | md5 (`.gz`) |
| --- | --- | --- |
| European / UK | [tag.HAPMAP.UK.tagging.gz](https://github.com/chaoning/gensep/releases/download/tagging-v1/tag.HAPMAP.UK.tagging.gz) | `88b22d43b17735114619a7dfbf2a7816` |
| Finnish | [tag.HAPMAP.FIN.tagging.gz](https://github.com/chaoning/gensep/releases/download/tagging-v1/tag.HAPMAP.FIN.tagging.gz) | `86e8c9aac5d764911f07ac9d85036e70` |
| East Asian | [tag.HAPMAP.EAS.tagging.gz](https://github.com/chaoning/gensep/releases/download/tagging-v1/tag.HAPMAP.EAS.tagging.gz) | `526f6fff81b1d779c62f976c08373cf8` |
| African / Caribbean | [tag.HAPMAP.CARAFR.tagging.gz](https://github.com/chaoning/gensep/releases/download/tagging-v1/tag.HAPMAP.CARAFR.tagging.gz) | `bc2016c0e5c98a375280a0b3b04472df` |
| Admixed American | [tag.HAPMAP.AMR.tagging.gz](https://github.com/chaoning/gensep/releases/download/tagging-v1/tag.HAPMAP.AMR.tagging.gz) | `c91abda8dacca7d81c1821efd54e2a83` |

Each covers ~1.1M HapMap3 SNPs. Download, decompress, and pass it to `--tagfile`:

```bash
wget https://github.com/chaoning/gensep/releases/download/tagging-v1/tag.HAPMAP.UK.tagging.gz
gunzip tag.HAPMAP.UK.tagging.gz          # -> tag.HAPMAP.UK.tagging (~56 MB)
gensep --se-method jackknife --tagfile tag.HAPMAP.UK.tagging \
       --summary trait1.summaries --summary2 trait2.summaries \
       --K1 <prev1> --K2 <prev2> --P1 <casefrac1> --P2 <casefrac2> --out output/pair
```

To build a tagging file for a different reference or SNP set instead, use LDAK
`--calc-tagging` (see [Calculate Taggings](https://dougspeed.com/calculate-taggings/)).

The jackknife shows live progress on `stderr`, with each stage reporting its elapsed time
(reading + QC, heritability/rg estimation, then a `block N / B` counter for the jackknife):

```
Read tagging + summaries and QC: 1.8 s
Trait 1: 527122 SNPs, weighted GIF 1.009, max variance explained 0.0016, weighted mean N 325246
Trait 2: 527122 SNPs, weighted GIF 1.136, max variance explained 0.0021, weighted mean N 348968
Estimating heritabilities and genetic correlation (sum-cors)... 0.7 s
Running 200-block jackknife...
  block 200 / 200  8.1 s
sum-cors diagnostics: intercept1=0.9357 intercept2=0.9003 overlap=0.0062
```

The `Trait 1/2` and `sum-cors diagnostics` lines are printed to `stderr` only (never written
to `PREFIX.gensep`), matching LDAK's SumHer diagnostics:

- **weighted GIF** — tagging-weighted genomic inflation factor (`> 1` flags inflation from
  stratification / cryptic relatedness / strong polygenicity → consider `--intercept YES` or
  genomic-control-corrected summaries).
- **max variance explained** — the largest single-SNP `rho² = chi/(chi+n)`; informs `--cutoff`.
- **weighted mean N** — effective sample size.
- **intercept1/2, overlap** — the per-trait LDSC-style intercepts and cross-trait
  sample-overlap term from sum-cors.

Each leave-one-block heritability solve is warm-started from the full-data fit, so even
single-threaded it is fast, and `--max-threads` scales the block loop across cores.
Progress goes to `stderr` only (the `PREFIX.gensep` result and the stdout summary are
unaffected); redirect it with `2>/dev/null` if you want silence.

## 2. Run from point estimates (`mc` / `delta` / `none`)

When you already have observed-scale subtype heritabilities and their genetic correlation
(for example from LDSC or a published table), skip the solvers and compute the derived
quantities directly.

### With standard errors (`mc` or `delta`)

```bash
gensep --se-method mc \
       --h1 0.10 --h2 0.15 --rg 0.30 \
       --K1 0.01 --K2 0.02 --P1 0.3 --P2 0.45 \
       --se-h1 0.02 --se-h2 0.03 --se-rg 0.10 \
       --num-draws 100000 --seed 1 \
       --out output/pair_point
```

- **`--h1` / `--h2`** — **observed-scale** SNP heritability of subtype 1 / 2.
- **`--rg`** — genetic correlation between the two subtypes.
- **`--se-h1` / `--se-h2` / `--se-rg`** — the standard errors of the three inputs. For
  `mc` and `delta` **all three are required**; each must be ≥ 0. `(h1, h2, rg)` are treated
  as **independent** — with only marginal SEs their estimation covariance is unknown, so it
  is assumed zero (this can mis-estimate `VS_se` via the `−2 λ1 λ2 rg √(h1 h2)` cross term;
  the `jackknife` route does not have this limitation).
- **`--num-draws`** — Monte-Carlo draws for `mc` (default `100000`, must be ≥ 2); ignored
  by `delta`.
- **`--seed`** — RNG seed for `mc` (default `1`), so runs are reproducible.

`mc` and `delta` agree to well under 1% away from boundaries; near `VS ≈ 0` they diverge
(there `delta` under-spreads or returns `NA`) and `mc` is the one to trust.

### Without standard errors (`none`)

If you only have point estimates and no SEs, use `none`. It computes point values only;
the SE column is `NA`. Passing any `--se-*` with `none` is an error (and, conversely,
`mc`/`delta` without all three `--se-*` is an error).

```bash
gensep --se-method none \
       --h1 0.10 --h2 0.15 --rg 0.30 \
       --K1 0.01 --K2 0.02 --P1 0.3 --P2 0.45 \
       --out output/pair_point_only
```

## 3. Output file

Every mode writes a whitespace-delimited `PREFIX.gensep`:

```
Quantity   Value      SE
hsq1_obs   0.100000   0.020000
hsq1_liab  0.055191   0.011038
hsq2_obs   0.150000   0.030000
hsq2_liab  0.098321   0.019664
rg         0.300000   0.100000
VS         0.683101   0.137592
h2cc       0.145865   0.024935
auc        0.723200   0.020546
auc_lo     0.720532   0.019819
# rg_used 0.300000  lam1 2.665214  lam2 2.420907  n_used(VS,h2cc) 100000  n_used(auc) 100000
```

| Row | Meaning |
| --- | --- |
| `hsq1_obs` / `hsq2_obs` | observed-scale SNP heritability of each subtype |
| `hsq1_liab` / `hsq2_liab` | liability-scale heritability, `= obs × Lee(K, P)` |
| `rg` | genetic correlation between the two subtypes |
| `VS` | genetic-separation variance (from the liability heritabilities and rg) |
| `h2cc` | case–case heritability, `VS / (VS + 4)` |
| `auc` | upper-limit (ceiling) AUC |
| `auc_lo` | leading-order AUC approximation |

The `SE` column is `NA` when a standard error is undefined (e.g. `--se-method none`, or a
quantity that is out of domain at the point). The trailing `#` comment records the clipped
rg used, the two selection intensities λ, and `n_used` — the number of jackknife blocks
(`jackknife`) or Monte-Carlo draws (`mc`) that were retained after dropping degenerate ones
(`VS ≤ 0`, or `denom ≤ 0` for `auc`).

A one-line summary of the same numbers is also printed to standard output.

## 4. Finite-PRS case-case AUC (`--auc1` / `--auc2`)

`auc` above is the genetic **ceiling** — the AUC if the total genetic value were known
exactly. If you also pass the per-subtype **PRS case/control AUC** measured on a test set,
GenSep reports the AUC achievable with those finite-accuracy PRS. These options work in
**every** `--se-method` mode (they use only the point `hsq*_liab`/`rg`):

```bash
gensep --se-method jackknife \
       --tagfile HumDef.tagging --summary trait1.summaries --summary2 trait2.summaries \
       --K1 0.01 --K2 0.02 --P1 0.3 --P2 0.45 \
       --auc1 0.75 --auc2 0.68 \
       --out output/pair
```

- **`--auc1` / `--auc2`** — PRS case/control AUC of subtype 1 / 2 on a held-out test set.
  Both-or-neither; each must lie in `(0.5, 0.9999)` (the upper bound is the
  `auc_to_corr_liab` domain — an AUC `≥ 0.9999` is rejected rather than silently `NA`).

GenSep converts each AUC to a PRS accuracy internally,
`Rsq_i = auc_to_corr_liab(auc_i, K_i)² / hsq_i_liab` (clipped — the same chain the
real-data pipeline uses), then evaluates the finite-PRS case-case AUC. Six **point-only**
rows are appended (`SE` is always `NA`), and the footer gains `Rsq1 Rsq2`:

| Row | Meaning |
| --- | --- |
| `prs_auc` | finite-PRS moment-corrected case-case AUC (optimal weights `w_B`) |
| `prs_auc_lo` | finite-PRS leading-order AUC (weights `w_LO`) |
| `h2cc_prs` | PRS case-case heritability, `V_PRS / (V_PRS + 4)` |
| `prs_eff` | PRS efficiency `V_PRS / VS_tgv ∈ [0, 1]` — fraction of the genetic separation the PRS captures |
| `prs_w1` / `prs_w2` | the optimal `w_B` weights for combining the two subtype PRS into the case-case discriminant: `D = prs_w1·PRS1 + prs_w2·PRS2` (higher `D` favours subtype 1). Scale is arbitrary — only the ratio matters |

`prs_auc ≤ auc` always (a finite PRS cannot beat the genetic ceiling). No SE is propagated
for the PRS-based quantities. The computation is a port of
`case_case_auc.compute_case_case_auc_prs`; at `Rsq1 = Rsq2 = 1` it recovers the ceiling `auc`.

## 5. Input rules and validation

- `--se-method` is required; `--K1/--K2/--P1/--P2` must all lie in `(0, 1)`, on **every**
  route. Neither pair is required to sum to 1 (see [Inputs](#inputs)).
- All numeric options are parsed strictly: a bad value such as `--h1 foo` or
  `--num-draws 3.5` produces a clear error instead of a crash or a silent `0`.
- SEs must be ≥ 0; `--num-draws` (for `mc`) must be ≥ 2; `--num-blocks` must be ≥ 2.
- The tagging file must be single-category; summary-statistic numeric fields are validated.

## 6. Worked example: type 1 vs type 2 diabetes in FinnGen

A complete run on public data, from download to result, takes about ten minutes (most of
it downloading). It uses two [FinnGen](https://www.finngen.fi/en) release-13 GWAS
(GRCh38; [Kurki *et al.*, *Nature* 2023](https://doi.org/10.1038/s41586-022-05473-8)),
which are freely downloadable from FinnGen's public bucket, and the ready-made Finnish
tagging file from the table above. Case and control counts are from the
[FinnGen R13 manifest](https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv).

| Subtype | FinnGen endpoint | Cases / controls | Download (~0.8 GB each) |
| --- | --- | --- | --- |
| Type 1 diabetes | `T1D_WIDE` (wide definition) | 11,197 / 396,409 | [finngen_R13_T1D_WIDE.gz](https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_T1D_WIDE.gz) |
| Type 2 diabetes | `T2D` (definitions combined) | 89,727 / 396,292 | [finngen_R13_T2D.gz](https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_T2D.gz) |

The whole example is scripted in
[`examples/run_finngen_t1d_t2d.sh`](https://github.com/chaoning/gensep/blob/main/examples/run_finngen_t1d_t2d.sh);
the steps it runs are:

**1. Download.**

```bash
FG=https://storage.googleapis.com/finngen-public-data-r13/summary_stats
wget $FG/finngen_R13_T1D_WIDE.gz
wget $FG/finngen_R13_T2D.gz
wget https://github.com/chaoning/gensep/releases/download/tagging-v1/tag.HAPMAP.FIN.tagging.gz
gunzip tag.HAPMAP.FIN.tagging.gz
```

**2. Convert to `.summaries`.** FinnGen files report `beta`/`sebeta` for the `alt`
allele and carry no per-SNP sample size, so
[`examples/finngen_to_summaries.py`](https://github.com/chaoning/gensep/blob/main/examples/finngen_to_summaries.py)
writes `Predictor A1 A2 Z n` with `A1 = alt`, `A2 = ref`, `Z = beta/sebeta` and
`n = cases + controls`, dropping variants without an rsID and rsIDs that occur more than
once. `--tagfile` keeps only the ~1.1M tagging SNPs (GenSep uses no others), which shrinks
each output to ~34 MB.

```bash
python3 finngen_to_summaries.py --in finngen_R13_T1D_WIDE.gz --cases 11197 --controls 396409 \
        --tagfile tag.HAPMAP.FIN.tagging --out t1d.summaries
python3 finngen_to_summaries.py --in finngen_R13_T2D.gz --cases 89727 --controls 396292 \
        --tagfile tag.HAPMAP.FIN.tagging --out t2d.summaries
```

**3. Run GenSep.** `K1`/`K2` are the population prevalences used for this pair in the
paper (0.5% for type 1, 7% for type 2; Supplementary Table 1); `P1 = 11,197 / 407,606 =
0.02747` and `P2 = 89,727 / 486,019 = 0.18462` are the case fractions of the two GWAS.

```bash
gensep --se-method jackknife \
       --tagfile tag.HAPMAP.FIN.tagging \
       --summary t1d.summaries --summary2 t2d.summaries \
       --K1 0.005 --K2 0.07 --P1 0.02747 --P2 0.18462 \
       --num-blocks 200 --out finngen_t1d_t2d
```

This takes about 30 s single-threaded. The log reports ~1,095,000 SNPs matched for both
traits, and `finngen_t1d_t2d.gensep` should read (also in
[`examples/finngen_t1d_t2d.expected.gensep`](https://github.com/chaoning/gensep/blob/main/examples/finngen_t1d_t2d.expected.gensep)):

```
Quantity Value SE
hsq1_obs 0.073107 0.040882
hsq1_liab 0.323938 0.181148
hsq2_obs 0.126523 0.005181
hsq2_liab 0.197582 0.008091
rg 0.297508 0.068786
VS 2.601128 1.360858
h2cc 0.394043 0.156602
auc 0.887240 0.096638
auc_lo 0.872945 0.087810
```

The two subtypes are genetically well separated: `rg` ≈ 0.30 and an oracle case–case AUC
of 0.89 (SE 0.10). The large SE on `hsq1` reflects the modest number of type 1 cases.

## 7. Summary workflow

1. **From GWAS**: build a tagging file with LDAK
   ([Calculate Taggings](https://dougspeed.com/calculate-taggings/)), prepare the two
   `.summaries` files, then run `gensep --se-method jackknife …`.
2. **From existing estimates**: run `gensep --se-method mc …` (or `delta`) with your h²/rg
   point estimates and their SEs, or `none` if you have no SEs.
3. Read the derived separation quantities and their SEs from `PREFIX.gensep`.
