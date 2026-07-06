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
PairData qc_pair(const Tagging& T, const SummaryAligned& S1, const SummaryAligned& S2);

}  // namespace gs
