---
layout: page
title: Method
description:
order: 4
---

This page documents the quantities gensep computes and the three standard-error
estimators. The point-estimate formulas are shared by every `--se-method`; only how the SE
is obtained changes.

## 1. From heritability and rg to the separation quantities

For each subtype `i`, let `K_i` be the population prevalence and `P_i` the sample case
fraction. Define the standard-normal threshold, selection intensity, and Lee factor:

```
t_i    = Φ⁻¹(1 − K_i)
λ_i    = φ(t_i) / K_i
d_i    = λ_i · (λ_i − t_i)
Lee_i  = ( K_i (1 − K_i) / φ(t_i) )² / ( P_i (1 − P_i) )
```

The observed-scale heritability is put on the liability scale with the Lee transform,
`h_i(liab) = h_i(obs) · Lee_i`. From the two liability heritabilities `h1, h2` and the
genetic correlation `rg` (clipped to ±0.999), gensep forms:

```
VS       = λ1² h1 + λ2² h2 − 2 λ1 λ2 rg √(h1 h2)     (genetic-separation variance)
h2cc     = VS / (VS + 4)                              (case–case heritability)
auc_lo   = Φ( √(VS / 2) )                             (leading-order AUC)
auc      = Φ( VS / √(2 VS − a1² d1 − a2² d2) )        (upper-limit / ceiling AUC)
          with  a1 = λ1 h1 − λ2 rg √(h1 h2),  a2 = λ2 h2 − λ1 rg √(h1 h2)
```

`VS`, `h2cc` and the two AUCs are **nonlinear** functions of `(h1, h2, rg)`, with domain
boundaries (`VS ≤ 0`, or `denom ≤ 0` for `auc`). This is what makes their standard error
non-trivial, and it is why gensep offers a joint estimator rather than propagating the
LDAK marginal SEs independently.

## 2. Finite-PRS case-case AUC (`--auc1` / `--auc2`)

`auc` above is the genetic **ceiling** — it assumes the total genetic value (TGV) is known
exactly. A real classifier uses a **polygenic score** of finite accuracy, so its case-case
AUC is lower. When you supply the per-subtype PRS case/control AUC (`--auc1`, `--auc2`),
gensep computes that finite-PRS AUC (a port of `case_case_auc.compute_case_case_auc_prs`).

**Step 1 — AUC → PRS accuracy.** Each PRS case/control AUC is converted to a PRS-to-TGV
squared correlation via the liability threshold model (the same chain the real-data
pipeline uses, and validated on simulation to `r_loglog ≈ 0.998` against the oracle):

```
ρ_i   = auc_to_corr_liab(auc_i, K_i)         # Corr(PRS_i, liability_i), Wray 2010
Rsq_i = clip( ρ_i² / hsq_i_liab , 1e-6, 0.999 )   # = Corr(PRS_i, TGV_i)²
```

**Step 2 — finite-PRS AUC.** With `h_i = √hsq_i`, `r_i = √Rsq_i`, the two PRS estimates
`ĝ1, ĝ2` are correlated (`ρ = r1 r2 rg`), so the separation runs through a 2×2 system:

```
a_PRS = [ r1(λ1 h1 − rg λ2 h2),  r2(λ2 h2 − rg λ1 h1) ]
Σ_PRS = [[1,−ρ],[−ρ,1]]
V_PRS = a_PRSᵀ Σ_PRS⁻¹ a_PRS        →  h2cc_prs = V_PRS/(V_PRS+4),  prs_eff = V_PRS/VS_tgv
B_PRS = 2Σ_PRS − δ1 d1d1ᵀ − δ2 d2d2ᵀ,   d1 = h1[r1,−r2 rg],  d2 = h2[r1 rg,−r2]
w_LO  = Σ_PRS⁻¹ a_PRS ,   w_B = B_PRS⁻¹ a_PRS         (optimal discriminant weights)
prs_auc    = Φ( wᵀa / √(2 wᵀΣw − δ1(wᵀd1)² − δ2(wᵀd2)²) )   using w_B  (moment-corrected)
prs_auc_lo = same with w_LO                                  (leading-order weights)
```

At `Rsq1 = Rsq2 = 1` this recovers the ceiling `auc`. All 2×2 inverses are closed form (no
Eigen). Outputs `prs_auc`, `prs_auc_lo`, `h2cc_prs`, `prs_eff` are **point-only** — no SE
is propagated for them — and are computed identically in every `--se-method` mode, since
they depend only on the point `hsq*_liab`/`rg` plus the supplied AUCs.

## 3. Standard errors

### `jackknife` — fused block-jackknife (summary-statistics route)

One common SNP set (SNPs with summary statistics for both traits, strand-ambiguous
dropped) is partitioned into `--num-blocks` blocks **once**. For each leave-one-block-out
replicate, `h1`/`h2` are re-estimated with the SumHer model and `rg` with the cross-trait
model on the **same** retained SNPs, then each derived quantity is recomputed. The SE is

```
SE = √[ (B − 1)/B · Σ_b (θ_b − θ̄)² ]        (LDAK convention)
```

Because all three of `h1, h2, rg` are recomputed together on each block, their joint
sampling covariance — including the strong coupling through the `−2 λ1 λ2 rg √(h1 h2)`
cross term — is captured automatically, with no delta method or covariance matrix.
Degenerate blocks (`VS_b ≤ 0`; `denom_b ≤ 0` for `auc`) are dropped from that quantity's
SE, and `n_used` reports how many remained.

### `mc` — Monte-Carlo propagation (point route)

Treat each input as a random variable with the supplied SE:

```
h1 ~ N(ĥ1, se_h1²),  h2 ~ N(ĥ2, se_h2²),  rg ~ N(r̂g, se_rg²)
```

Draw `--num-draws` samples (reproducibly, via `--seed`), push each through the exact same
point-estimate formulas, and take the **sample standard deviation** of the resulting
`VS / h2cc / auc / auc_lo` values as the SE. Degenerate draws are dropped, mirroring the
jackknife's block deletion. Monte-Carlo reconstructs the true (possibly skewed, bounded)
sampling distribution, so it is accurate through the nonlinearity and at the boundaries;
its only approximation here is the independence assumption below.

### `delta` — first-order propagation (point route)

Propagate with a first-order Taylor expansion,

```
Var(Q) = Σ_i ( ∂Q/∂θ_i · se_i )²,   θ = (h1_obs, h2_obs, rg)
```

where the gradient is obtained by **central finite differences on the same point-estimate
function** (no hand-derived analytic gradients). This is exact and deterministic when the
function is locally near-linear, and it needs no random draws — but it is a symmetric
linear approximation and becomes **unreliable near a boundary** (`VS` or `denom → 0`),
where it under-spreads or returns `NA`. `n_used` is not applicable and prints `0`.

### Which to use

- The `hsq*_obs/liab` and `rg` SEs are always exact (`liab_se = Lee · obs_se`; `rg_se` as
  given), independent of the method.
- Away from boundaries, `mc` and `delta` agree to well under 1% — a useful mutual check.
- Near `VS ≈ 0` they diverge; prefer `mc`, whose spread reflects the truncated,
  skewed distribution that `delta` cannot represent.

## 4. The independence caveat (point route)

Both `mc` and `delta` propagate whatever covariance you give them, and with only **marginal**
SEs (`se_h1, se_h2, se_rg`) the estimation covariance among `h1, h2, rg` is unknown — so it
is assumed **zero**. Because `VS` contains the `−2 λ1 λ2 rg √(h1 h2)` cross term, a real
(usually positive) correlation among the inputs would change `VS_se`. Treat the point-mode
`VS_se` / `auc_se` as an **independence-assumption approximation**. The `jackknife` route
does not have this limitation, because it recomputes `h1, h2, rg` jointly on each block and
therefore carries their full covariance.

## 5. Limitations

- **Single-category tagging only.** The sum-cors solver supports one heritability category
  (`num_parts == 1`); a multi-category tagging file is rejected.
- **Out-of-domain points report `NA`.** A non-positive `VS` (or `denom` for `auc`) is out
  of the AUC/`h2cc` domain; the point value is still reported where defined, but the
  corresponding SE (and undefined AUC) is `NA`.
- **Point-route SEs assume input independence** (Section 4).

## 6. Reference

gensep is a self-contained C++ port of LDAK SumHer (`--sum-hers` / `--sum-cors`), validated
to six decimals against LDAK for the underlying `.hers` / `.cors` solvers (including the
jackknife SEs), and against a reference `case_case_auc` implementation for the VS / h2cc /
AUC point estimates. See [LDAK](https://dougspeed.com/ldak/) for the SumHer method
(Speed & Balding, *Nature Genetics* 2019).
