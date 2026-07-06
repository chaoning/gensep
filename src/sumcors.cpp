#include "sumcors.hpp"
#include "common.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <cstdio>

namespace gs {

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Solve symmetric A x = b, replicating LDAK eigen_invert: diagonal-scale so the
// diagonal is 1, eigendecompose, zero scaled eigenvalues with |e|<1e-6, pseudo-invert.
static VectorXd ei_solve(MatrixXd A, const VectorXd& b) {
    const int n = (int)A.rows();
    VectorXd sc(n);
    for (int j = 0; j < n; ++j) sc[j] = A(j, j) > 0 ? 1.0 / std::sqrt(A(j, j)) : 1.0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) A(i, j) *= sc[i] * sc[j];
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(A);
    const VectorXd& ev = es.eigenvalues();
    const MatrixXd& V  = es.eigenvectors();
    VectorXd inv(n);
    for (int j = 0; j < n; ++j) inv[j] = (std::fabs(ev[j]) >= 1e-6) ? 1.0 / ev[j] : 0.0;
    VectorXd rb(n);
    for (int j = 0; j < n; ++j) rb[j] = b[j] * sc[j];
    VectorXd y = V * (inv.asDiagonal() * (V.transpose() * rb));
    for (int j = 0; j < n; ++j) y[j] *= sc[j];
    return y;
}

// One IRLS fit. cols holds the k weighted design columns (col-major), sY the
// weighted response, after convergence (the last build used to get `theta`).
struct Fit {
    VectorXd theta;
    std::vector<double> exps;            // final updated per-SNP expectation
    std::vector<std::vector<double>> cols;
    std::vector<double> sY;
};

// LDAK jackknife block bounds: start=(int)(p*L/B), end=(int)((p+1)*L/B).
static inline void block_bounds(int p, int L, int B, int& s, int& e) {
    s = (int)((double)p * L / B);
    e = (int)((double)(p + 1) * L / B);
}

// Full X'X and X'y from weighted columns, then per-block leave-one-out theta.
static std::vector<VectorXd> jack_thetas(const Fit& f, int L, int B) {
    const int k = (int)f.cols.size();
    MatrixXd XtX = MatrixXd::Zero(k, k);
    VectorXd Xty = VectorXd::Zero(k);
    for (int j = 0; j < L; ++j) {
        for (int a = 0; a < k; ++a) {
            Xty[a] += f.cols[a][j] * f.sY[j];
            for (int b = 0; b < k; ++b) XtX(a, b) += f.cols[a][j] * f.cols[b][j];
        }
    }
    std::vector<VectorXd> out(B);
    for (int p = 0; p < B; ++p) {
        int s, e; block_bounds(p, L, B, s, e);
        MatrixXd bX = MatrixXd::Zero(k, k);
        VectorXd bY = VectorXd::Zero(k);
        for (int j = s; j < e; ++j) {
            for (int a = 0; a < k; ++a) {
                bY[a] += f.cols[a][j] * f.sY[j];
                for (int b = 0; b < k; ++b) bX(a, b) += f.cols[a][j] * f.cols[b][j];
            }
        }
        out[p] = ei_solve(XtX - bX, Xty - bY);
    }
    return out;
}

// IRLS for one trait: design = [ n/scale*svar , n/scale (cept) ], response = chi-1,
// weights sW = stags*exps^2. (num_parts=1, gcon=0, cept=1 => k=2)
static Fit irls_trait(int L, const double* stg, const double* chi, const double* nss,
                      const double* sv, double scale, double ss00,
                      double tol, int maxiter) {
    const int k = 2;
    Fit f;
    f.cols.assign(k, std::vector<double>(L));
    f.sY.assign(L, 0.0);
    f.exps.assign(L, 1.0);
    std::vector<double> sW(L);
    VectorXd theta = VectorXd::Zero(k);
    double sumold = 0;
    int count = 0;
    while (true) {
        for (int j = 0; j < L; ++j) {
            sW[j] = stg[j] * f.exps[j] * f.exps[j];
            double v = 1.0 / std::sqrt(sW[j]);
            f.cols[0][j] = nss[j] / scale * sv[j] * v;
            f.cols[1][j] = nss[j] / scale * v;          // cept column
            f.sY[j] = (chi[j] - 1) * v;
        }
        MatrixXd XtX = MatrixXd::Zero(k, k);
        VectorXd Xty = VectorXd::Zero(k);
        for (int j = 0; j < L; ++j)
            for (int a = 0; a < k; ++a) {
                Xty[a] += f.cols[a][j] * f.sY[j];
                for (int b = 0; b < k; ++b) XtX(a, b) += f.cols[a][j] * f.cols[b][j];
            }
        theta = ei_solve(XtX, Xty);
        double sumhers = theta[0] * ss00 / scale;        // gc=1
        double diff = sumhers - sumold;
        sumold = sumhers;
        // update exps = 1 + (X*theta)*sqrt(sW)
        for (int j = 0; j < L; ++j) {
            double xt = f.cols[0][j] * theta[0] + f.cols[1][j] * theta[1];
            double e = 1 + xt * std::sqrt(sW[j]);
            f.exps[j] = e <= 0 ? 1e-6 : e;
        }
        if (std::fabs(diff) < tol) break;
        if (count == maxiter) break;
        ++count;
    }
    f.theta = theta;
    return f;
}

CorsResult sum_cors(const PairData& D, int num_blocks, double tol, int maxiter) {
    if (D.num_parts != 1) die("sum_cors: only single-category tagging supported");
    const int L = D.n;
    int B = num_blocks; if (B > L) B = L;
    const double* stg = D.stags.data();
    const double* sv  = D.svars[0].data();
    const double ss00 = D.ssums[0][0];

    // scales (sumfuns.c:1105)
    double s1 = 0, s2 = 0, s3 = 0;
    for (int j = 0; j < L; ++j) { s1 += D.snss[j]/stg[j]; s2 += D.snss2[j]/stg[j]; s3 += 1.0/stg[j]; }
    double scale = s1 / s3, scale2 = s2 / s3, scale3 = std::sqrt(scale) * std::sqrt(scale2);

    // snss3 = sqrt(n1 n2);  schis3 = +/- sqrt(chi1 chi2) with sign(rho1*rho2)
    std::vector<double> snss3(L), schis3(L);
    for (int j = 0; j < L; ++j) {
        snss3[j] = std::sqrt(D.snss[j]) * std::sqrt(D.snss2[j]);
        double m = std::sqrt(D.schis[j]) * std::sqrt(D.schis2[j]);
        schis3[j] = (D.srhos[j] * D.srhos2[j] > 0) ? m : -m;
    }

    // ---- Trait 1 & Trait 2 ----
    Fit f1 = irls_trait(L, stg, D.schis.data(),  D.snss.data(),  sv, scale,  ss00, tol, maxiter);
    Fit f2 = irls_trait(L, stg, D.schis2.data(), D.snss2.data(), sv, scale2, ss00, tol, maxiter);

    CorsResult r;
    r.her1 = f1.theta[0] * ss00 / scale;
    r.her2 = f2.theta[0] * ss00 / scale2;
    r.intercept1 = 1 + f1.theta[1] / 1.0;   // gc=1
    r.intercept2 = 1 + f2.theta[1] / 1.0;

    // ---- Coheritability (uses f1.exps, f2.exps) ----
    const int k = 2;
    Fit fc;
    fc.cols.assign(k, std::vector<double>(L));
    fc.sY.assign(L, 0.0);
    fc.exps.assign(L, 0.0);                  // exps3 starts at 0
    std::vector<double> sW(L);
    VectorXd thc = VectorXd::Zero(k);
    double sumold = 0; int count = 0;
    while (true) {
        for (int j = 0; j < L; ++j) {
            double w = stg[j] * (f1.exps[j] * f2.exps[j] + fc.exps[j] * fc.exps[j]);
            if (w <= 0) w = 1e-6;
            sW[j] = w;
            double v = 1.0 / std::sqrt(w);
            fc.cols[0][j] = snss3[j] / scale3 * sv[j] * v;
            fc.cols[1][j] = v;                // oversamp overlap column
            fc.sY[j] = schis3[j] * v;
        }
        MatrixXd XtX = MatrixXd::Zero(k, k);
        VectorXd Xty = VectorXd::Zero(k);
        for (int j = 0; j < L; ++j)
            for (int a = 0; a < k; ++a) {
                Xty[a] += fc.cols[a][j] * fc.sY[j];
                for (int b = 0; b < k; ++b) XtX(a, b) += fc.cols[a][j] * fc.cols[b][j];
            }
        thc = ei_solve(XtX, Xty);
        double sumh3 = thc[0] * ss00 / scale3;   // gc3=1
        double diff = sumh3 - sumold;
        sumold = sumh3;
        for (int j = 0; j < L; ++j) {
            double xt = fc.cols[0][j] * thc[0] + fc.cols[1][j] * thc[1];
            fc.exps[j] = xt * std::sqrt(sW[j]);  // no +1
        }
        if (std::fabs(diff) < tol) break;
        if (count == maxiter) break;
        ++count;
    }
    fc.theta = thc;
    r.coher   = thc[0] * ss00 / scale3;
    r.overlap = thc[1];                          // gc3=1
    r.cor     = r.coher * std::pow(r.her1 * r.her2, -0.5);

    // ---- jackknife (delete-one-block) ----
    auto j1 = jack_thetas(f1, L, B);
    auto j2 = jack_thetas(f2, L, B);
    auto jc = jack_thetas(fc, L, B);
    r.B = B;
    r.her1_b.resize(B); r.her2_b.resize(B); r.coher_b.resize(B); r.cor_b.resize(B);
    for (int p = 0; p < B; ++p) {
        double h1 = j1[p][0] * ss00 / scale;
        double h2 = j2[p][0] * ss00 / scale2;
        double ch = jc[p][0] * ss00 / scale3;
        r.her1_b[p] = h1;
        r.her2_b[p] = h2;
        r.coher_b[p] = ch;
        r.cor_b[p] = ch * std::pow(h1 * h2, -0.5);
    }
    // LDAK-faithful: NaN in any block propagates to SE (e.g. Cor when Her1*Her2<0
    // gives a NaN block -> SE=NaN, matching LDAK's "-nan"). Drop-in for --sum-cors.
    r.her1_se  = jackknife_se_ldak(r.her1_b);
    r.her2_se  = jackknife_se_ldak(r.her2_b);
    r.coher_se = jackknife_se_ldak(r.coher_b);
    r.cor_se   = jackknife_se_ldak(r.cor_b);
    return r;
}

}  // namespace gs
