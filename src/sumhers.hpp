// gensep — SumHer single-trait heritability point solve, block-aware.
// Port of LDAK solve_sums (sumfuns.c:245) Newton path for the default --sum-hers
// config with a single-category tagging: gcon=0, cept=0, chisol=1.
//   E[chi_j] = 1 + (n_j/scale) * svar_j * theta ;   h2 = theta * ss00 / scale
#pragma once

namespace gs {

// Newton solve on raw per-SNP arrays, optionally EXCLUDING SNPs in [skip_s, skip_e)
// (set them equal to use all). `scale` is supplied (computed once on the full set)
// so jackknife replicates share one scale, mirroring the sum-cors jackknife.
// Returns h2 (observed scale).
double sumhers_h2_block(int L, const double* stags, const double* schis,
                        const double* snss, const double* svar, double ss00,
                        double scale, int skip_s, int skip_e,
                        double tol = 0.001, int maxiter = 100);

}  // namespace gs
