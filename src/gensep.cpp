#include "gensep.hpp"
#include "common.hpp"
#include "sumhers.hpp"
#include "sumcors.hpp"
#include <cmath>
#include <cstdio>
#include <random>
#include <atomic>
#include <chrono>

namespace gs {

static inline double lee_factor(double K, double P) {
    if (!(K > 0 && K < 1 && P > 0 && P < 1)) return NAN;
    double t = normal_inv(1 - K);
    double v = K * (1 - K) / normal_pdf(t);
    return v * v / (P * (1 - P));
}
static inline double clip_rg(double rg) {
    if (rg < -0.999) return -0.999;
    if (rg >  0.999) return  0.999;
    return rg;
}

// The four derived quantities from liability h1,h2 + rg (case_case_auc.py).
// Returns false-ish via NaN where out of domain (VS<0 / denom<0).
struct Derived { double VS, h2cc, auc_exact, auc_approx; };
static Derived derive(double h1, double h2, double rg,
                      double lam1, double lam2, double d1, double d2) {
    Derived o{NAN, NAN, NAN, NAN};
    if (!(h1 > 0 && h2 > 0)) return o;
    double sq = std::sqrt(h1 * h2);
    double VS = lam1 * lam1 * h1 + lam2 * lam2 * h2 - 2 * lam1 * lam2 * rg * sq;
    o.VS = VS;
    o.h2cc = VS / (VS + 4.0);
    if (VS > 0) {
        o.auc_approx = normal_cdf(std::sqrt(VS / 2.0));
        double a1 = lam1 * h1 - lam2 * rg * sq;
        double a2 = lam2 * h2 - lam1 * rg * sq;
        double denom = 2.0 * VS - a1 * a1 * d1 - a2 * a2 * d2;
        if (denom > 0) o.auc_exact = normal_cdf(VS / std::sqrt(denom));
    }
    return o;
}

// ---- finite-PRS case-case AUC (port of case_case_auc.compute_case_case_auc_prs) ----

// PRS case/control AUC -> Corr(PRS, liability) via the liability threshold model
// (Wray 2010); matches step5_metrics.auc_to_corr_liab. NaN outside the valid domain.
static double auc_to_corr_liab(double auc, double K) {
    if (!(auc > 0.5 && auc < 0.9999) || !(K > 0 && K < 1)) return NAN;
    double t = normal_inv(1 - K), phi = normal_pdf(t);
    double ic = phi / K, ictrl = -phi / (1 - K);
    double S = ic * (ic - t) + ictrl * (ictrl - t);
    double m = phi / (K * (1 - K));
    double z = normal_inv(auc);
    double rho2 = 2.0 * z * z / (m * m + z * z * S);
    if (!(rho2 > 0) || !std::isfinite(rho2)) return NAN;
    return std::sqrt(rho2 < 0.999 ? rho2 : 0.999);
}

// hsq1/hsq2: liability heritabilities. Rsq_i = Corr(PRS_i, TGV_i)^2. lam/del as in
// derive(). Two 2x2 inverses in closed form (no Eigen). auc_b uses moment-corrected
// weights w_B, auc_lo the leading-order weights w_LO (both with the full moment
// denominator). NaN where out of domain (Sigma/B singular, variance<=0).
struct PrsDerived { double auc_b, auc_lo, h2cc_prs, prs_eff; };
static PrsDerived derive_prs(double hsq1, double hsq2, double rg, double Rsq1, double Rsq2,
                             double lam1, double lam2, double del1, double del2) {
    PrsDerived o{NAN, NAN, NAN, NAN};
    if (!(hsq1 > 0 && hsq2 > 0 && Rsq1 > 0 && Rsq2 > 0)) return o;
    double h1 = std::sqrt(hsq1), h2 = std::sqrt(hsq2);
    double r1 = std::sqrt(Rsq1), r2 = std::sqrt(Rsq2);
    double VS_tgv = lam1 * lam1 * hsq1 + lam2 * lam2 * hsq2 - 2 * lam1 * lam2 * rg * h1 * h2;
    double rho = r1 * r2 * rg;
    if (!(std::fabs(rho) < 1.0)) return o;                // Sigma_PRS singular
    double a1 = r1 * (lam1 * h1 - rg * lam2 * h2);
    double a2 = r2 * (lam2 * h2 - rg * lam1 * h1);
    // Sigma = [[1,-rho],[-rho,1]] ; Sigma^-1 = 1/(1-rho^2) [[1,rho],[rho,1]]
    double inv = 1.0 / (1.0 - rho * rho);
    double wlo1 = inv * (a1 + rho * a2), wlo2 = inv * (rho * a1 + a2);
    double V_PRS = a1 * wlo1 + a2 * wlo2;                 // a^T Sigma^-1 a
    o.h2cc_prs = V_PRS / (V_PRS + 4.0);
    o.prs_eff  = VS_tgv > 0 ? V_PRS / VS_tgv : 0.0;
    // conditional-variance vectors: d1 = h1[r1,-r2 rg], d2 = h2[r1 rg,-r2]
    double d1x = h1 * r1, d1y = -h1 * r2 * rg;
    double d2x = h2 * r1 * rg, d2y = -h2 * r2;
    // B = 2 Sigma - del1 d1 d1^T - del2 d2 d2^T (2x2 symmetric)
    double B11 = 2.0        - del1 * d1x * d1x - del2 * d2x * d2x;
    double B12 = -2.0 * rho - del1 * d1x * d1y - del2 * d2x * d2y;
    double B22 = 2.0        - del1 * d1y * d1y - del2 * d2y * d2y;
    double detB = B11 * B22 - B12 * B12;
    auto auc_for_w = [&](double w1, double w2) -> double {
        double dw = w1 * a1 + w2 * a2;
        if (dw < 0) { w1 = -w1; w2 = -w2; dw = -dw; }
        double wSw = w1 * w1 - 2.0 * rho * w1 * w2 + w2 * w2;   // w^T Sigma w
        double wd1 = w1 * d1x + w2 * d1y, wd2 = w1 * d2x + w2 * d2y;
        double var = 2.0 * wSw - wd1 * wd1 * del1 - wd2 * wd2 * del2;
        return var > 0 ? normal_cdf(dw / std::sqrt(var)) : NAN;
    };
    o.auc_lo = auc_for_w(wlo1, wlo2);
    if (detB != 0.0 && std::isfinite(detB)) {
        double wb1 = (B22 * a1 - B12 * a2) / detB, wb2 = (-B12 * a1 + B11 * a2) / detB;
        o.auc_b = auc_for_w(wb1, wb2);
    }
    return o;
}

// Fill the finite-PRS outputs from the point estimates (no SE). auc_i -> rho_i
// (liability inversion) -> Rsq_i = rho_i^2 / hsq_i_liab (clipped, as step5), then
// derive_prs(). Identical across every --se-method mode.
static void fill_prs_auc(SepResult& r, double auc1, double auc2, double K1, double K2,
                         double lam1, double lam2, double del1, double del2) {
    r.have_prs = true;
    auto rsq = [](double auc, double K, double hl) -> double {
        double rho = auc_to_corr_liab(auc, K);
        if (!(std::isfinite(rho) && hl > 0)) return NAN;
        double v = rho * rho / hl;
        return v < 1e-6 ? 1e-6 : (v > 0.999 ? 0.999 : v);
    };
    r.Rsq1 = rsq(auc1, K1, r.hsq1_liab);
    r.Rsq2 = rsq(auc2, K2, r.hsq2_liab);
    PrsDerived p = derive_prs(r.hsq1_liab, r.hsq2_liab, r.rg_used, r.Rsq1, r.Rsq2,
                              lam1, lam2, del1, del2);
    r.prs_auc = p.auc_b; r.prs_auc_lo = p.auc_lo;
    r.h2cc_prs = p.h2cc_prs; r.prs_eff = p.prs_eff;
}

// Shared core: given liability-scale point + per-block (h1,h2,rg), fill all derived
// quantities + their block-jackknife SE. A block with VS_b<=0 is degenerate (VS is a
// variance; h2cc=VS/(VS+4) near its pole; AUC undefined) -> dropped from ALL derived
// SEs; auc_exact additionally drops denom<=0 blocks. The point estimate is untouched
// (reports real VS even if <0). (sep is new functionality, not an LDAK output, so
// pair-deletion + B_used is by design.)
static void fill_derived(SepResult& r,
                         double h1l, double h2l, double rg_pt,
                         const std::vector<double>& h1l_b,
                         const std::vector<double>& h2l_b,
                         const std::vector<double>& rg_b,
                         double lam1, double lam2, double d1, double d2) {
    r.lam1 = lam1; r.lam2 = lam2;
    r.hsq1_liab = h1l; r.hsq2_liab = h2l; r.rg_used = clip_rg(rg_pt);

    Derived pe = derive(h1l, h2l, r.rg_used, lam1, lam2, d1, d2);
    r.VS = pe.VS; r.h2cc = pe.h2cc; r.auc_exact = pe.auc_exact; r.auc_approx = pe.auc_approx;

    const int B = (int)h1l_b.size();
    std::vector<double> vs(B), hc(B), ae(B), aa(B);
    for (int b = 0; b < B; ++b) {
        Derived d = derive(h1l_b[b], h2l_b[b], clip_rg(rg_b[b]), lam1, lam2, d1, d2);
        if (!(d.VS > 0)) { vs[b] = hc[b] = ae[b] = aa[b] = NAN; }
        else { vs[b] = d.VS; hc[b] = d.h2cc; ae[b] = d.auc_exact; aa[b] = d.auc_approx; }
    }
    int nb;
    r.VS_se         = jackknife_se(vs, &nb); r.B_used = nb;
    r.h2cc_se       = jackknife_se(hc, &nb);
    r.auc_approx_se = jackknife_se(aa, &nb);
    r.auc_exact_se  = jackknife_se(ae, &nb); r.B_used_auc = nb;
}

SepResult gene_sep_fused(const PairData& D, double K1, double K2, double P1, double P2,
                         int num_blocks, bool have_auc, double auc1, double auc2,
                         bool intercept) {
    const int L = D.n;
    int B = num_blocks; if (B > L) B = L;
    double t1 = normal_inv(1 - K1), t2 = normal_inv(1 - K2);
    double lam1 = normal_pdf(t1) / K1, lam2 = normal_pdf(t2) / K2;
    double d1 = lam1 * (lam1 - t1), d2 = lam2 * (lam2 - t2);
    double c1 = lee_factor(K1, P1), c2 = lee_factor(K2, P2);

    const double* stg = D.stags.data();
    const double* sv  = D.svars[0].data();
    const double ss00 = D.ssums[0][0];

    // one fixed scale per trait on the full common set (mirrors sum-cors)
    double s1 = 0, s2 = 0, s3 = 0;
    for (int j = 0; j < L; ++j) { s1 += D.snss[j]/stg[j]; s2 += D.snss2[j]/stg[j]; s3 += 1.0/stg[j]; }
    double scale1 = s1 / s3, scale2 = s2 / s3;

    using clk = std::chrono::steady_clock;
    auto secs = [](clk::time_point a) { return std::chrono::duration<double>(clk::now() - a).count(); };

    // rg (sum-cors), which also yields per-trait heritabilities fit WITH a free intercept
    // (cept=1) plus their block-jackknife — used directly when --intercept YES.
    std::fprintf(stderr, "Estimating heritabilities and genetic correlation (sum-cors)...");
    std::fflush(stderr);
    auto t_pt = clk::now();
    CorsResult cors = sum_cors(D, B);
    std::vector<double> h1_b(B), h2_b(B), h1l_b(B), h2l_b(B);
    double h1obs, h2obs;

    if (intercept) {
        // --intercept YES: heritabilities from the sum-cors per-trait fit (free intercept),
        // already jackknifed on the same B blocks -> h1/h2/rg fully consistent (same
        // solver / scale / blocks), and no extra SumHer solves needed.
        h1obs = cors.her1; h2obs = cors.her2;
        for (int p = 0; p < B; ++p) {
            h1_b[p] = cors.her1_b[p]; h2_b[p] = cors.her2_b[p];
            h1l_b[p] = h1_b[p] * c1;  h2l_b[p] = h2_b[p] * c2;
        }
        std::fprintf(stderr, " %.1f s\n", secs(t_pt));
    } else {
        // --intercept NO (default): standalone SumHer h2 with a fixed intercept of 1.
        h1obs = sumhers_h2_block(L, stg, D.schis.data(),  D.snss.data(),  sv, ss00, scale1, 0, 0);
        h2obs = sumhers_h2_block(L, stg, D.schis2.data(), D.snss2.data(), sv, ss00, scale2, 0, 0);
        std::fprintf(stderr, " %.1f s\n", secs(t_pt));

        // warm-start seeds for the leave-one-block solves (h2 = theta*ss00/scale)
        double th1 = h1obs * scale1 / ss00, th2 = h2obs * scale2 / ss00;
        std::fprintf(stderr, "Running %d-block jackknife...\n", B);
        auto t_jk = clk::now();
        std::atomic<int> done{0};
        const int step = B / 20 > 0 ? B / 20 : 1;             // progress every ~5%
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int p = 0; p < B; ++p) {
            int s = (int)((double)p * L / B), e = (int)((double)(p + 1) * L / B);
            h1_b[p] = sumhers_h2_block(L, stg, D.schis.data(),  D.snss.data(),  sv, ss00, scale1, s, e, 0.001, 100, th1);
            h2_b[p] = sumhers_h2_block(L, stg, D.schis2.data(), D.snss2.data(), sv, ss00, scale2, s, e, 0.001, 100, th2);
            h1l_b[p] = h1_b[p] * c1;
            h2l_b[p] = h2_b[p] * c2;
            int n = ++done;
            if (n % step == 0 || n == B) std::fprintf(stderr, "\r  block %d / %d", n, B);
        }
        std::fprintf(stderr, "  %.1f s\n", secs(t_jk));
    }

    // sum-cors diagnostics to the log only (NOT written to PREFIX.gensep): the per-trait
    // LDSC-style intercepts (>1 = inflation) and the cross-trait sample-overlap term.
    std::fprintf(stderr, "sum-cors diagnostics: intercept1=%.4f intercept2=%.4f overlap=%.4f\n",
                 cors.intercept1, cors.intercept2, cors.overlap);

    SepResult r;
    fill_derived(r, h1obs * c1, h2obs * c2, cors.cor, h1l_b, h2l_b, cors.cor_b, lam1, lam2, d1, d2);
    r.hsq1_obs = h1obs; r.hsq1_obs_se = jackknife_se(h1_b);
    r.hsq2_obs = h2obs; r.hsq2_obs_se = jackknife_se(h2_b);
    r.hsq1_liab_se = c1 * r.hsq1_obs_se;   // Lee is a constant multiplier
    r.hsq2_liab_se = c2 * r.hsq2_obs_se;
    r.rg = clip_rg(cors.cor); r.rg_se = cors.cor_se;   // rg SE faithful to sum-cors (may be NaN)
    if (have_auc) fill_prs_auc(r, auc1, auc2, K1, K2, lam1, lam2, d1, d2);
    return r;
}

// Streaming mean/SD accumulator with pair-deletion (skips non-finite samples).
namespace {
struct Acc {
    double sum = 0, sumsq = 0; long n = 0;
    void add(double x) { if (std::isfinite(x)) { sum += x; sumsq += x * x; ++n; } }
    double sd() const {                      // sample SD = Monte-Carlo SE of the quantity
        if (n < 2) return NAN;
        double m = sum / n, v = (sumsq - sum * m) / (n - 1);
        return v > 0 ? std::sqrt(v) : 0.0;
    }
};
}  // namespace

// Point-input mode. See gensep.hpp. Reuses derive() / lee_factor() / clip_rg() so the
// point estimate is byte-identical to the fused mode given the same liability h/rg.
SepResult gene_sep_point(double h1obs, double h2obs, double rg,
                         double K1, double K2, double P1, double P2,
                         bool have_se, double se_h1, double se_h2, double se_rg,
                         SeMethod method, long num_draws, unsigned long long seed,
                         bool have_auc, double auc1, double auc2) {
    double t1 = normal_inv(1 - K1), t2 = normal_inv(1 - K2);
    double lam1 = normal_pdf(t1) / K1, lam2 = normal_pdf(t2) / K2;
    double d1 = lam1 * (lam1 - t1), d2 = lam2 * (lam2 - t2);
    double c1 = lee_factor(K1, P1), c2 = lee_factor(K2, P2);

    SepResult r;
    r.lam1 = lam1; r.lam2 = lam2;
    r.hsq1_obs = h1obs; r.hsq2_obs = h2obs;
    r.hsq1_liab = h1obs * c1; r.hsq2_liab = h2obs * c2;
    r.rg = clip_rg(rg); r.rg_used = r.rg;

    // Point estimate of the derived quantities (same derive() as the fused path).
    Derived pe = derive(r.hsq1_liab, r.hsq2_liab, r.rg_used, lam1, lam2, d1, d2);
    r.VS = pe.VS; r.h2cc = pe.h2cc; r.auc_exact = pe.auc_exact; r.auc_approx = pe.auc_approx;

    // Finite-PRS AUC (point-only) — set before any SE branch so every return carries it.
    if (have_auc) fill_prs_auc(r, auc1, auc2, K1, K2, lam1, lam2, d1, d2);

    if (!have_se) {                          // point-only: mark every SE as NA
        r.hsq1_obs_se = r.hsq1_liab_se = r.hsq2_obs_se = r.hsq2_liab_se = NAN;
        r.rg_se = r.VS_se = r.h2cc_se = r.auc_exact_se = r.auc_approx_se = NAN;
        return r;
    }

    // Analytic SEs: liability = Lee(K,P) * observed, Lee a constant multiplier.
    r.hsq1_obs_se = se_h1; r.hsq1_liab_se = c1 * se_h1;
    r.hsq2_obs_se = se_h2; r.hsq2_liab_se = c2 * se_h2;
    r.rg_se = se_rg;

    // --- Delta method: Var(Q) = sum_i (dQ/dtheta_i * se_i)^2, gradient by central
    // finite-difference on the SAME derive() (theta = h1obs, h2obs, rg; independent).
    // A boundary crossing within a step makes a perturbed Q non-finite -> SE NA. ---
    if (method == SE_DELTA) {
        auto eval = [&](double h1o, double h2o, double rgv) {
            return derive(h1o * c1, h2o * c2, clip_rg(rgv), lam1, lam2, d1, d2);
        };
        const double e1 = 1e-5, e2 = 1e-5, er = 1e-5;
        Derived p1 = eval(h1obs + e1, h2obs, rg), m1 = eval(h1obs - e1, h2obs, rg);
        Derived p2 = eval(h1obs, h2obs + e2, rg), m2 = eval(h1obs, h2obs - e2, rg);
        Derived pr = eval(h1obs, h2obs, rg + er), mr = eval(h1obs, h2obs, rg - er);
        auto se_of = [&](double a1, double b1, double a2, double b2, double ar, double br) {
            double g1 = (a1 - b1) / (2 * e1), g2 = (a2 - b2) / (2 * e2), gr = (ar - br) / (2 * er);
            double v = g1 * g1 * se_h1 * se_h1 + g2 * g2 * se_h2 * se_h2 + gr * gr * se_rg * se_rg;
            return std::sqrt(v);
        };
        r.VS_se         = se_of(p1.VS, m1.VS, p2.VS, m2.VS, pr.VS, mr.VS);
        r.h2cc_se       = se_of(p1.h2cc, m1.h2cc, p2.h2cc, m2.h2cc, pr.h2cc, mr.h2cc);
        r.auc_approx_se = se_of(p1.auc_approx, m1.auc_approx, p2.auc_approx, m2.auc_approx, pr.auc_approx, mr.auc_approx);
        r.auc_exact_se  = se_of(p1.auc_exact, m1.auc_exact, p2.auc_exact, m2.auc_exact, pr.auc_exact, mr.auc_exact);
        r.B_used = r.B_used_auc = 0;         // n_used is Monte-Carlo-only
        return r;
    }

    // Monte-Carlo SE for the nonlinear derived quantities. Draw (h1,h2,rg) as
    // INDEPENDENT Normal(point, se), convert to liability, push through derive();
    // the SD across draws is the SE. Degenerate draws (VS<=0, and denom<=0 for
    // auc_exact) are dropped, matching fill_derived's jackknife pair-deletion.
    Acc aVS, aHC, aAE, aAA;
    std::mt19937_64 gen(seed);
    std::normal_distribution<double> Z(0.0, 1.0);
    for (long i = 0; i < num_draws; ++i) {
        double h1l = (h1obs + se_h1 * Z(gen)) * c1;
        double h2l = (h2obs + se_h2 * Z(gen)) * c2;
        double rgd = clip_rg(rg + se_rg * Z(gen));
        Derived d = derive(h1l, h2l, rgd, lam1, lam2, d1, d2);
        if (!(d.VS > 0)) continue;           // degenerate: drop from ALL derived SEs
        aVS.add(d.VS); aHC.add(d.h2cc); aAA.add(d.auc_approx); aAE.add(d.auc_exact);
    }
    r.VS_se = aVS.sd();           r.B_used = (int)aVS.n;
    r.h2cc_se = aHC.sd();
    r.auc_approx_se = aAA.sd();
    r.auc_exact_se = aAE.sd();    r.B_used_auc = (int)aAE.n;
    return r;
}

void write_gensep(const std::string& prefix, const SepResult& r) {
    std::string path = prefix + ".gensep";
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) die("cannot write " + path);
    std::fprintf(f, "Quantity Value SE\n");
    auto row = [&](const char* n, double v, double se) {
        if (std::isfinite(v) && std::isfinite(se)) std::fprintf(f, "%s %.6f %.6f\n", n, v, se);
        else if (std::isfinite(v))                 std::fprintf(f, "%s %.6f NA\n", n, v);
        else                                       std::fprintf(f, "%s NA NA\n", n);
    };
    row("hsq1_obs",   r.hsq1_obs,  r.hsq1_obs_se);
    row("hsq1_liab",  r.hsq1_liab, r.hsq1_liab_se);
    row("hsq2_obs",   r.hsq2_obs,  r.hsq2_obs_se);
    row("hsq2_liab",  r.hsq2_liab, r.hsq2_liab_se);
    row("rg",         r.rg,        r.rg_se);
    row("VS",         r.VS,        r.VS_se);
    row("h2cc",       r.h2cc,      r.h2cc_se);
    row("auc",        r.auc_exact, r.auc_exact_se);
    row("auc_lo",     r.auc_approx, r.auc_approx_se);
    if (r.have_prs) {                         // finite-PRS case-case AUC, point-only (SE=NA)
        row("prs_auc",    r.prs_auc,    NAN);
        row("prs_auc_lo", r.prs_auc_lo, NAN);
        row("h2cc_prs",   r.h2cc_prs,   NAN);
        row("prs_eff",    r.prs_eff,    NAN);
    }
    std::fprintf(f, "# rg_used %.6f  lam1 %.6f  lam2 %.6f"
                 "  n_used(VS,h2cc) %d  n_used(auc) %d",
                 r.rg_used, r.lam1, r.lam2, r.B_used, r.B_used_auc);
    if (r.have_prs) std::fprintf(f, "  Rsq1 %.6f  Rsq2 %.6f", r.Rsq1, r.Rsq2);
    std::fprintf(f, "\n");
    std::fclose(f);
}

}  // namespace gs
