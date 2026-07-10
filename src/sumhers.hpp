// gensep — SumHer single-trait heritability point solve, block-aware.
// Port of LDAK solve_sums (sumfuns.c:245) Newton path for the default --sum-hers
// config with a single-category tagging: gcon=0, cept=0, chisol=1.
//   E[chi_j] = 1 + (n_j/scale) * svar_j * theta ;   h2 = theta * ss00 / scale
#pragma once

namespace gs {

// Newton solve on raw per-SNP arrays, optionally EXCLUDING SNPs in [skip_s, skip_e)
// (set them equal to use all). `scale` is supplied (computed once on the full set)
// so jackknife replicates share one scale, mirroring the sum-cors jackknife.
// `theta_init` seeds the Newton iteration: for a leave-one-block replicate, passing the
// full-data theta (= h2_full * scale / ss00) is a near-exact warm start (the block drops
// ~1/num_blocks of SNPs), so the solve converges in a couple of iterations instead of
// many. It only affects the starting point, not the optimum reached — the result is
// identical to starting from 0.
// Returns h2 (observed scale).
double sumhers_h2_block(int L, const double* stags, const double* schis,
                        const double* snss, const double* svar, double ss00,
                        double scale, int skip_s, int skip_e,
                        double tol = 0.001, int maxiter = 100, double theta_init = 0.0);

}  // namespace gs
