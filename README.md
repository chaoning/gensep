# gensep

Genetic separation from GWAS summary statistics. One command: given a tagging file,
two traits' summaries, and prevalences/sample fractions, it outputs

    hsq1 (obs+liab)  hsq2 (obs+liab)  rg  VS  h2cc  auc_exact  auc_approx

each with a **block-jackknife SE**.

A self-contained C++ port of LDAK SumHer (`--sum-hers` / `--sum-cors`), extended with
the block-jackknife SE on the derived genetic-separation quantities that LDAK does not
expose (VS / case-case h² / upper-limit AUC are nonlinear functions of jointly-estimated
h1/h2/rg, so their SE needs a joint jackknife, not LDAK's marginal SEs).

Design/math in [plan.md](plan.md).

## Build

```sh
make            # -> ./gensep   (C++17, Eigen header-only, statically linked)
```

The binary is **fully static** (`ldd` → "not a dynamic executable", ~1.2 MB stripped),
so it runs on any compatible Linux host with no libc/libstdc++ dependency. This needs
static libc/libm/libstdc++, which the system GCC provides but the conda toolchain does
not — so the Makefile pins `CXX := /usr/bin/g++` (override with `make CXX=...`).

Eigen is vendored at `third_party/eigen` (symlink to
`/faststorage/project/gensubtypes/chao/code/fastgxe/external/eigen-5.0.0`).
The Makefile uses `CXXFLAGS :=` (not `?=`) so a conda-exported `CXXFLAGS` cannot
drop `-Ithird_party/eigen`.

## Run

One command; `--se-method` (**required, no default**) selects both the input and the SE
estimator. All three write the same `PREFIX.gensep` layout.

### `--se-method jackknife` — from summary statistics

```sh
gensep --se-method jackknife \
       --tagfile HumDef.tagging --summary g1.summaries --summary2 g2.summaries \
       --K1 <prev1> --K2 <prev2> --P1 <casefrac1> --P2 <casefrac2> \
       [--num-blocks 200] --out PREFIX
```

Estimates h1/h2 (SumHer), rg (sum-cors) and the derived quantities on one common SNP
set, SE via the fused block-jackknife (see Method).

> **Single-category tagging only.** The sum-cors solver supports one heritability
> category (`num_parts == 1`, e.g. `HumDef`/`ldak-thin` GBAT); a multi-category tagging
> file is rejected with an error. `--num-blocks` must be ≥ 2.

### `--se-method mc|delta|none` — from given point estimates

When you already have observed-scale subtype h² and rg (e.g. from LDSC / a published
table), skip the solvers and compute the derived quantities directly:

```sh
gensep --se-method mc --h1 <h2obs1> --h2 <h2obs2> --rg <rg> \
       --K1 <prev1> --K2 <prev2> --P1 <casefrac1> --P2 <casefrac2> \
       --se-h1 s --se-h2 s --se-rg s [--num-draws 100000] [--seed 1] --out PREFIX
```

`--h1/--h2` are **observed-scale** h² for the two subtypes. The point estimates use the
exact same `derive()` as the jackknife mode.

The `--se-*` inputs are gated by the method:

- **`mc` / `delta`** require **all three** of `--se-h1`, `--se-h2`, `--se-rg` (SEs are
  emitted; how they propagate is below).
- **`none`** is for when you have **no SEs**: it computes point values only, the SE column
  is `NA`, and passing any `--se-*` is an error. (This is the answer to "I only have point
  estimates" — use `--se-method none`, not a dummy `mc`.)

- `hsq*_obs/liab` and `rg` SEs are analytic (`liab_se = Lee(K,P) · obs_se`, `rg_se` as given).
- `VS`, `h2cc`, `auc`, `auc_lo` SEs are propagated per `--se-method`:
  - **`mc`** — Monte-Carlo: draw `(h1,h2,rg)` as **independent** `Normal(point, se)`
    (`--num-draws`, default 1e5; `--seed`, default 1), push each draw through `derive()`,
    take the sample SD across draws. Degenerate draws (`VS≤0`; `denom≤0` for `auc`) are
    dropped, as the jackknife pair-deletion; `n_used(...)` reports how many draws remained.
  - **`delta`** — first-order propagation `Var(Q) = Σᵢ (∂Q/∂θᵢ · seᵢ)²`, gradient by
    central finite-difference on the same `derive()`. Deterministic and exact when
    `derive()` is locally near-linear; **unreliable near a boundary** (`VS`/`denom→0`,
    where it under-spreads or returns `NA`). `n_used` prints 0 (not applicable).
  - The two agree to <1% away from boundaries (good mutual check); near `VS≈0` they
    diverge and `mc` is the one to trust.
- **Independence caveat (both methods):** with only marginal SEs the h1/h2/rg estimation
  covariance is unknown, so it is assumed zero. Because VS contains the
  `−2·λ1·λ2·rg·√(h1·h2)` cross term, a real (usually positive) h1/h2/rg correlation
  would change `VS_se` — treat the point-mode `VS_se`/`auc_se` as an independence
  approximation.

Both write `PREFIX.gensep`:

```
Quantity      Value     SE
hsq1_obs      ...       ...     # observed-scale h2, trait 1 (sep: SumHer)
hsq1_liab     ...       ...     # liability-scale h2, trait 1 (= obs * Lee(K1,P1))
hsq2_obs      ...       ...     # observed-scale h2, trait 2
hsq2_liab     ...       ...     # liability-scale h2, trait 2
rg            ...       ...     # genetic correlation (sep: sum-cors)
VS            ...       ...     # genetic separation variance (built from liability h2)
h2cc          ...       ...     # case-case h2 = VS/(VS+4)  (observed 50/50 scale)
auc           ...       ...     # upper-limit (ceiling) AUC
auc_lo        ...       ...     # leading-order AUC
# rg_used .. lam1 .. lam2 .. n_used(VS,h2cc) .. n_used(auc) ..   # n_used = blocks (sep) / draws (point)
```

`hsq*_liab` SE is just `Lee(K,P) * hsq*_obs SE` (Lee is a constant scaling).

`K1,K2` (population prevalences) drive the selection intensity λ and the Lee
observed→liability transform; `P1,P2` (sample case fractions, not in the summaries)
are needed by the Lee factor.

## Method (fused jackknife)

One common SNP set (SNPs with summary stats for both traits, strand-ambiguous dropped),
partitioned into 200 blocks **once**. Per leave-one-block-out:
- h1, h2 via the SumHer model (Newton, no intercept),
- rg via the cross-trait model (with intercept),

all on the same retained SNPs, so h1_b/h2_b/rg_b are aligned → VS_b, h2cc_b, AUC_b →
`SE = sqrt[(B-1)/B · Σ(θ_b - θ̄)²]` (LDAK convention). Point estimates use the same
common-set estimators, so point and SE are consistent. The covariance among h1,h2,rg is
captured by recomputing all three on the same blocks (no delta method / covariance matrix).

Degenerate blocks (`VS_b<=0`; `denom<=0` for auc_exact) are dropped from that quantity's
SE; `B_used(...)` reports how many blocks remained. The point estimate still reports the
real VS even if negative. rg is clipped to ±0.999 per block (and at the point).

## Validation

The underlying solvers were validated to 6 decimals against LDAK (`.hers`, `.cors`
including all jackknife SEs and the `-nan` propagation for degenerate Cor blocks) and
against the Python `case_case_auc` module (VS/h2cc/AUC point estimates) before the tool
was reduced to this single command. Defaults matched from LDAK source: hers
`amb=1→here amb=0 (common set), gcon=0, cept=0, chisol=1, tol=1e-3`; cors `amb=0, gcon=0,
cept=1, oversamp=1, num_blocks=200, tol=1e-4`.

Note: because h1/h2 are estimated on the common (both-trait, amb=0) SNP set, they differ
slightly from a standalone LDAK `--sum-hers` run (which uses the full single-trait,
amb=1 set).
