// gensep — sum-cors (bivariate SumHer). Port of LDAK solve_cors (sumfuns.c:1047)
// for the default --sum-cors config with single-category tagging:
//   num_parts=1, gcon=0, cept=1, oversamp=1, num_blocks=200.
// Three IRLS least-squares fits (Trait1 h2, Trait2 h2, coheritability), each with
// delete-one-block jackknife; correlation rg = Coher / sqrt(Her1*Her2).
// The per-block [Her1,Her2,Coher,Cor] are kept for the derived-quantity SE (`sep`).
#pragma once
#include "qc.hpp"

namespace gs {

struct CorsResult {
    double her1 = 0, her2 = 0, coher = 0, cor = 0;
    double her1_se = 0, her2_se = 0, coher_se = 0, cor_se = 0;
    double intercept1 = 0, intercept2 = 0, overlap = 0;
    int B = 0;
    std::vector<double> her1_b, her2_b, coher_b, cor_b;  // [B] leave-one-block-out
};

CorsResult sum_cors(const PairData& D, int num_blocks = 200,
                    double tol = 0.0001, int maxiter = 100);

}  // namespace gs
