#include "sumhers.hpp"
#include "common.hpp"
#include <vector>
#include <cmath>

namespace gs {

double sumhers_h2_block(int L, const double* stg, const double* chi,
                        const double* nss, const double* sv, double ss00,
                        double scale, int skip_s, int skip_e,
                        double tol, int maxiter) {
    std::vector<double> sT(L);
    for (int j = 0; j < L; ++j) sT[j] = nss[j] / scale * sv[j];
    auto inblk = [&](int j) { return j >= skip_s && j < skip_e; };
    auto loglik = [&](const std::vector<double>& e) {
        double sum = 0, s2 = 0, val = 0;
        for (int j = 0; j < L; ++j) {
            if (inblk(j)) continue;
            sum += chi[j] / e[j] / stg[j];
            s2  += 1.0 / stg[j];
            val += std::log(chi[j] * e[j]) / stg[j];
        }
        return -0.5 * sum - 0.5 * val - 0.5 * s2 * std::log(2 * M_PI);
    };
    auto make_exps = [&](double th, std::vector<double>& e) {
        for (int j = 0; j < L; ++j) { double v = 1 + sT[j] * th; e[j] = v <= 0 ? 1e-6 : v; }
    };
    double theta = 0, like, likeold, diff = 0;
    std::vector<double> exps(L), exps2(L);
    make_exps(theta, exps); like = loglik(exps);
    int count = 0, rflag = 0;
    while (true) {
        if (count > 0) diff = like - likeold;
        likeold = like;
        if (count > 0 && std::fabs(diff) < tol && rflag == 0) break;
        if (count == maxiter) break;
        double AI = 0, BI = 0;
        for (int j = 0; j < L; ++j) {
            if (inblk(j)) continue;
            double e = exps[j];
            AI += (chi[j] - 0.5 * e) / stg[j] * std::pow(e, -3) * sT[j] * sT[j];
            BI += 0.5 * (chi[j] - e) / stg[j] * std::pow(e, -2) * sT[j];
        }
        double td = BI / AI, relax = 1; bool moved = false;
        while (relax > 0.001) {
            double tth = theta + relax * td;
            make_exps(tth, exps2);
            double l2 = loglik(exps2);
            if (l2 > like - tol) { theta = tth; like = l2; exps = exps2; rflag = 0; moved = true; break; }
            relax *= 0.5;
        }
        if (!moved) { rflag = 1; break; }
        ++count;
    }
    return theta * ss00 / scale;
}

}  // namespace gs
