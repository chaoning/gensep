// gensep — derived genetic-separation quantities (the `sep` command).
// Reuses code/common/case_case_auc.py formulas (validated there); adds the
// block-jackknife SE that step5 currently lacks.
//
//   lam_i = phi(Phi^-1(1-K_i))/K_i ;  Lee_i = (K_i(1-K_i)/phi(t_i))^2 / (P_i(1-P_i))
//   h_i(liab) = h_i(obs) * Lee_i
//   VS            = lam1^2 h1 + lam2^2 h2 - 2 lam1 lam2 rg sqrt(h1 h2)
//   h2cc_formula  = VS/(VS+4)                          (case-case h2, formula route only)
//   auc_exact     = Phi( VS / sqrt(2VS - a1^2 d1 - a2^2 d2) )
//   auc_approx    = Phi( sqrt(VS/2) )
//
// Point estimates keep the current pipeline sources (h1/h2 from sum-hers, rg from
// sum-cors); the SE is the delete-one-block jackknife over the sum-cors per-block
// (Her1,Her2,Cor) — recomputing each quantity per block with the same Lee/lam
// constants, then SE = sqrt[(B-1)/B * sum(theta_b - mean)^2] (LDAK convention).
#pragma once
#include <string>
#include <vector>
#include "qc.hpp"

namespace gs {

struct SepResult {
    double VS = 0,    VS_se = 0;
    double h2cc = 0,  h2cc_se = 0;     // formula route: VS/(VS+4)
    double auc_exact = 0, auc_exact_se = 0;
    double auc_approx = 0, auc_approx_se = 0;
    double rg_used = 0, lam1 = 0, lam2 = 0;
    int B_used = 0;       // blocks with VS_b>0 (used for VS_se / h2cc_se)
    int B_used_auc = 0;   // blocks with VS_b>0 AND denom>0 (used for auc_exact_se)
    // heritabilities, observed + liability scale, each with jackknife SE.
    // liability = observed * Lee(K,P); since Lee is a constant, liab_se = Lee*obs_se.
    double hsq1_obs = 0, hsq1_obs_se = 0, hsq1_liab = 0, hsq1_liab_se = 0;
    double hsq2_obs = 0, hsq2_obs_se = 0, hsq2_liab = 0, hsq2_liab_se = 0;
    double rg = 0,       rg_se = 0;
};

// Point estimate inputs: h1obs/h2obs (observed-scale, from sum-hers), rg (sum-cors).
// her1_b/her2_b/cor_b: the per-block observed-scale values from the .cors.jackknife.
// Fused jackknife (option A): one common SNP set (PairData), 200 shared blocks.
// Per leave-one-block-out: h1/h2 via the sum-hers solver (cept=0), rg via sum-cors;
// point estimates of h1/h2/rg also on the common set -> point and SE fully consistent.
SepResult gene_sep_fused(const PairData& D, double K1, double K2, double P1, double P2,
                         int num_blocks = 200);

// SE method for the nonlinear derived quantities (VS/h2cc/auc/auc_lo) in point mode.
//   SE_MC    — Monte-Carlo: draw (h1,h2,rg) ~ independent Normal(point,se), push each
//              through derive(), take the sample SD across draws (degenerate draws
//              VS<=0 / denom<=0 dropped, as the jackknife pair-deletion). Captures the
//              nonlinearity/boundaries; needs a seed; n_used = draws retained.
//   SE_DELTA — first-order propagation Var = sum_i (dQ/dtheta_i * se_i)^2 with the
//              gradient by central finite-difference on the SAME derive(). Exact &
//              deterministic when derive() is locally near-linear; unreliable near a
//              boundary (VS/denom->0), where it yields NA. n_used is not applicable.
//   SE_NONE  — point-input, point estimates only (caller passes have_se=false); no SE.
enum SeMethod { SE_MC = 0, SE_DELTA = 1, SE_NONE = 2 };

// Point-input mode ("point" command): the caller supplies the observed-scale
// subtype h2 (h1obs/h2obs), rg and K/P directly — no summary stats / solvers. The
// derived quantities use the SAME derive() as gene_sep_fused. If have_se, the derived
// SEs use `method` (above); the h1/h2/rg obs/liab SEs are always analytic (liab = Lee *
// obs, a linear map). Independence is assumed for BOTH methods (no h1/h2/rg covariance
// is available from marginal SEs) — this can mis-estimate VS_se via the
// -2 lam1 lam2 rg sqrt(h1 h2) cross term. Without se, only point values are filled.
SepResult gene_sep_point(double h1obs, double h2obs, double rg,
                         double K1, double K2, double P1, double P2,
                         bool have_se, double se_h1, double se_h2, double se_rg,
                         SeMethod method = SE_MC,
                         long num_draws = 100000, unsigned long long seed = 1);

void write_gensep(const std::string& prefix, const SepResult& r);

}  // namespace gs
