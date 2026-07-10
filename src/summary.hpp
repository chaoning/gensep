// gensep — summary-statistics reader. Port of LDAK read_sumsfile (filemain.c:1204)
// for the Z-column format produced by step2 (--linear): "Predictor A1 A2 Z n A1Freq".
// Aligns to the tagging predictor order; SNPs absent / allele-inconsistent /
// strand-ambiguous (when amb=0) stay snss=0 and are dropped later in QC.
#pragma once
#include <vector>
#include <string>
#include "tagging.hpp"

namespace gs {

struct SummaryAligned {
    std::vector<double> snss;   // [T.n] sample size (0 if unused)
    std::vector<double> schis;  // [T.n] chi-squared (= Z^2 * scaling)
    std::vector<double> srhos;  // [T.n] signed sqrt(chi/(chi+n)), allele-oriented
    int n_matched = 0;          // predictors actually used (LDAK count3)
};

// amb=0 drops strand-ambiguous (A/T, C/G) SNPs (LDAK default); scaling default 1.
SummaryAligned read_sumsfile(const std::string& path, const Tagging& T,
                             int amb = 0, double scaling = 1.0, bool verbose = true);

}  // namespace gs
