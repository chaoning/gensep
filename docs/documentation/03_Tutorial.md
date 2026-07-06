---
layout: page
title: Tutorial
description:
order: 3
---

This tutorial follows the current command-line workflow implemented in the codebase. The
executable is assumed to be on your `PATH` as `gensep`; if you built from source and it is
not, replace `gensep` with the full path, e.g. `./gensep`.

`gensep` is one command. A single required switch, **`--se-method`**, selects both the
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
       --P1 0.5  --P2 0.5 \
       --num-blocks 200 \
       --out output/pair
```

### Inputs

- **`--tagfile`** — an LDAK tagging file, produced by `ldak --calc-tagging`. gensep
  supports **single-category** tagging only (e.g. `HumDef` / an `ldak-thin` model); a
  multi-category file is rejected with an error.
- **`--summary` / `--summary2`** — the two traits' summary statistics, in LDAK
  `.summaries` format: a header line followed by rows of

  ```
  Predictor  A1  A2  Z  n
  ```

  where `Predictor` is the SNP ID (matched against the tagging file), `A1`/`A2` are the
  alleles, `Z` is the signed Z-score, and `n` is the per-SNP sample size. Non-numeric or
  malformed `Z`/`n` fields are rejected rather than silently coerced.
- **`--K1` / `--K2`** — population **prevalences** of the two subtypes (each in `(0, 1)`).
  These drive the selection intensity λ and the Lee observed→liability transform.
- **`--P1` / `--P2`** — the **sample case fractions** of the two subtypes (each in
  `(0, 1)`). These are not stored in the summaries, so you must supply them; they enter the
  Lee factor.
- **`--num-blocks`** — number of jackknife blocks (default `200`, must be ≥ 2).
- **`--out`** — output prefix; results are written to `PREFIX.gensep`.

## 2. Run from point estimates (`mc` / `delta` / `none`)

When you already have observed-scale subtype heritabilities and their genetic correlation
(for example from LDSC or a published table), skip the solvers and compute the derived
quantities directly.

### With standard errors (`mc` or `delta`)

```bash
gensep --se-method mc \
       --h1 0.10 --h2 0.15 --rg 0.30 \
       --K1 0.01 --K2 0.02 --P1 0.5 --P2 0.5 \
       --se-h1 0.02 --se-h2 0.03 --se-rg 0.10 \
       --num-draws 100000 --seed 1 \
       --out output/pair_point
```

- **`--h1` / `--h2`** — **observed-scale** SNP heritability of subtype 1 / 2.
- **`--rg`** — genetic correlation between the two subtypes.
- **`--se-h1` / `--se-h2` / `--se-rg`** — the standard errors of the three inputs. For
  `mc` and `delta` **all three are required**; each must be ≥ 0. `(h1, h2, rg)` are treated
  as independent (see [Method](./04_Method.html) for the caveat).
- **`--num-draws`** — Monte-Carlo draws for `mc` (default `100000`, must be ≥ 2); ignored
  by `delta`.
- **`--seed`** — RNG seed for `mc` (default `1`), so runs are reproducible.

`mc` and `delta` agree to well under 1% away from boundaries; near `VS ≈ 0` they diverge
and `mc` is the one to trust. See [Method](./04_Method.html) for the full comparison.

### Without standard errors (`none`)

If you only have point estimates and no SEs, use `none`. It computes point values only;
the SE column is `NA`. Passing any `--se-*` with `none` is an error (and, conversely,
`mc`/`delta` without all three `--se-*` is an error).

```bash
gensep --se-method none \
       --h1 0.10 --h2 0.15 --rg 0.30 \
       --K1 0.01 --K2 0.02 --P1 0.5 --P2 0.5 \
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
gensep reports the AUC achievable with those finite-accuracy PRS. These options work in
**every** `--se-method` mode (they use only the point `hsq*_liab`/`rg`):

```bash
gensep --se-method jackknife \
       --tagfile HumDef.tagging --summary trait1.summaries --summary2 trait2.summaries \
       --K1 0.01 --K2 0.02 --P1 0.5 --P2 0.5 \
       --auc1 0.75 --auc2 0.68 \
       --out output/pair
```

- **`--auc1` / `--auc2`** — PRS case/control AUC of subtype 1 / 2 on a held-out test set.
  Both-or-neither; each must lie in `(0.5, 1)`.

gensep converts each AUC to a PRS accuracy internally,
`Rsq_i = auc_to_corr_liab(auc_i, K_i)² / hsq_i_liab` (clipped — the same chain the
real-data pipeline uses), then evaluates the finite-PRS case-case AUC. Four **point-only**
rows are appended (`SE` is always `NA`), and the footer gains `Rsq1 Rsq2`:

| Row | Meaning |
| --- | --- |
| `prs_auc` | finite-PRS moment-corrected case-case AUC (optimal weights `w_B`) |
| `prs_auc_lo` | finite-PRS leading-order AUC (weights `w_LO`) |
| `h2cc_prs` | PRS case-case heritability, `V_PRS / (V_PRS + 4)` |
| `prs_eff` | PRS efficiency `V_PRS / VS_tgv ∈ [0, 1]` — fraction of the genetic separation the PRS captures |

`prs_auc ≤ auc` always (a finite PRS cannot beat the genetic ceiling). No SE is propagated
for the PRS-based quantities. See the [Method](./04_Method.html) page for the formula.

## 5. Input rules and validation

- `--se-method` is required; `--K1/--K2/--P1/--P2` must all lie in `(0, 1)`, on **every**
  route.
- All numeric options are parsed strictly: a bad value such as `--h1 foo` or
  `--num-draws 3.5` produces a clear error instead of a crash or a silent `0`.
- SEs must be ≥ 0; `--num-draws` (for `mc`) must be ≥ 2; `--num-blocks` must be ≥ 2.
- The tagging file must be single-category; summary-statistic numeric fields are validated.

## 6. Summary workflow

1. **From GWAS**: build a tagging file with LDAK, prepare the two `.summaries` files, then
   run `gensep --se-method jackknife …`.
2. **From existing estimates**: run `gensep --se-method mc …` (or `delta`) with your h²/rg
   point estimates and their SEs, or `none` if you have no SEs.
3. Read the derived separation quantities and their SEs from `PREFIX.gensep`.
