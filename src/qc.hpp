// gensep — QC + dataset assembly. Port of the sumsa.c QC block (L116-158):
//   * keep only SNPs with a summary statistic (snss>0), squeeze down
//   * stags<1  -> 1
//   * schis==0 -> 1e-6
// (--cutoff / --cutoff2 truncation paths are LDAK defaults -9999 = off.)
#pragma once
#include <vector>
#include <string>
#include "tagging.hpp"
#include "summary.hpp"

namespace gs {

// Two traits' per-SNP data for the fused estimation.
struct PairData {
    int n = 0;
    int num_parts = 0;
    std::vector<double> stags;
    std::vector<std::vector<double>> svars;   // [num_parts][n]
    std::vector<std::vector<double>> ssums;    // [num_parts][num_parts+2]
    std::vector<double> snss,  schis,  srhos;   // trait 1
    std::vector<double> snss2, schis2, srhos2;  // trait 2
};

// Paired QC (for `cors`): keep SNPs with snss>0 AND snss2>0 (sumsb.c:135).
// cutoff > 0 drops predictors explaining >= cutoff of phenotypic variance in either trait
// (strong-effect loci that can bias SumHer h2/rg). cutoff <= 0 disables it (default).
PairData qc_pair(const Tagging& T, const SummaryAligned& S1, const SummaryAligned& S2,
                 double cutoff = 0.0);

}  // namespace gs
