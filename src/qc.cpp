#include "qc.hpp"
#include "common.hpp"
#include <cstdio>

namespace gs {

PairData qc_pair(const Tagging& T, const SummaryAligned& S1, const SummaryAligned& S2,
                 double cutoff) {
    const int P = T.num_parts;
    std::vector<int> keep;
    keep.reserve(T.n);
    int overlap = 0, excl = 0;
    double maxr2 = 0;                              // max variance explained (either trait)
    for (int j = 0; j < T.n; ++j) {
        if (!(S1.snss[j] > 0 && S2.snss[j] > 0)) continue;      // both traits present
        ++overlap;
        double r1 = S1.srhos[j] * S1.srhos[j], r2 = S2.srhos[j] * S2.srhos[j];
        if (r1 > maxr2) maxr2 = r1;
        if (r2 > maxr2) maxr2 = r2;
        // --cutoff: drop strong-effect loci (rho^2 = variance explained >= cutoff in either
        // trait). Equivalent to LDAK sum-cors (per-trait n=0 then intersection).
        if (cutoff > 0 && (r1 >= cutoff || r2 >= cutoff)) { ++excl; continue; }
        keep.push_back(j);
    }
    if (keep.empty()) die("no predictors with summary statistics for both traits");
    std::fprintf(stderr, "Overlap: %d predictors with summary statistics for both traits\n", overlap);
    if (cutoff > 0)
        std::fprintf(stderr, "--cutoff %.4g: excluded %d strong-effect predictor(s) "
                     "(variance explained >= %.4g in either trait); %d remain\n",
                     cutoff, excl, cutoff, (int)keep.size());
    else if (maxr2 > 0.01)
        std::fprintf(stderr, "Warning: a predictor explains up to %.4f of phenotypic variance; "
                     "strong-effect loci can bias SumHer h2/rg -- consider --cutoff (e.g. --cutoff 0.01)\n",
                     maxr2);

    PairData D;
    D.n = (int)keep.size();
    D.num_parts = P;
    D.stags.resize(D.n);
    D.svars.assign(P, std::vector<double>(D.n));
    D.ssums = T.ssums;
    D.snss.resize(D.n);  D.schis.resize(D.n);  D.srhos.resize(D.n);
    D.snss2.resize(D.n); D.schis2.resize(D.n); D.srhos2.resize(D.n);

    for (int i = 0; i < D.n; ++i) {
        int j = keep[i];
        double tag = T.stags[j];
        if (tag < 1) tag = 1;
        D.stags[i] = tag;
        for (int q = 0; q < P; ++q) D.svars[q][i] = T.svars[q][j];
        D.snss[i]   = S1.snss[j];
        D.schis[i]  = (S1.schis[j] == 0) ? 1e-6 : S1.schis[j];
        D.srhos[i]  = S1.srhos[j];
        D.snss2[i]  = S2.snss[j];
        D.schis2[i] = (S2.schis[j] == 0) ? 1e-6 : S2.schis[j];
        D.srhos2[i] = S2.srhos[j];
    }
    return D;
}

}  // namespace gs
