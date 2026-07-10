#include "summary.hpp"
#include "common.hpp"
#include <unordered_map>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cmath>

namespace gs {

static int find_head(const std::vector<char*>& hdr, const char* name) {
    for (int i = 0; i < (int)hdr.size(); ++i)
        if (std::strcmp(hdr[i], name) == 0) return i;
    return -1;
}

SummaryAligned read_sumsfile(const std::string& path, const Tagging& T,
                             int amb, double scaling, bool verbose) {
    std::FILE* f = std::fopen(path.c_str(), "r");
    if (!f) die("cannot open summary file " + path);

    // tagging predictor name -> index (reference order)
    std::unordered_map<std::string, int> idx;
    idx.reserve(T.n * 2);
    for (int j = 0; j < T.n; ++j) idx.emplace(T.preds[j], j);

    SummaryAligned S;
    S.snss.assign(T.n, 0.0);
    S.schis.assign(T.n, 0.0);
    S.srhos.assign(T.n, 0.0);

    char* buf = nullptr; size_t cap = 0; ssize_t len;
    std::vector<char*> tok;

    // ---- header: locate columns ----
    if ((len = getline(&buf, &cap, f)) == -1) die("empty summary file " + path);
    std::string hdr_line(buf, (size_t)std::max<ssize_t>(len, 0));
    int nc = split_ws(&hdr_line[0], tok);
    int cPred = find_head(tok, "Predictor"); if (cPred < 0) cPred = find_head(tok, "SNP");
    int cA1   = find_head(tok, "A1");
    int cA2   = find_head(tok, "A2");
    int cZ    = find_head(tok, "Z");
    int cN    = find_head(tok, "n");  if (cN < 0) cN = find_head(tok, "N");
    if (cPred < 0 || cA1 < 0 || cA2 < 0 || cZ < 0 || cN < 0)
        die("summary file " + path + " must have columns Predictor/SNP A1 A2 Z n");
    int need = std::max({cPred, cA1, cA2, cZ, cN});

    int ccount = 0, acount = 0, lcount = 0;  // inconsistent / ambiguous / long-allele
    std::vector<char> seen(T.n, 0);          // dedup: first match wins (find_strings)

    while ((len = getline(&buf, &cap, f)) != -1) {
        if (len <= 0) continue;
        nc = split_ws(buf, tok);
        if (nc <= need) continue;

        const char* a1s = tok[cA1];
        const char* a2s = tok[cA2];
        if (std::strlen(a1s) > 1 || std::strlen(a2s) > 1) { ++lcount; continue; }  // multi-char

        auto it = idx.find(tok[cPred]);
        if (it == idx.end()) continue;
        int j2 = it->second;
        if (seen[j2]) continue;

        char g1 = a1s[0], g2 = a2s[0];
        char t1 = T.al1[j2], t2 = T.al2[j2];

        // allele consistency (LDAK L1431)
        bool inconsistent = ((t1 != g1 && t1 != g2) || (t2 != g1 && t2 != g2));
        if (inconsistent) { ++ccount; continue; }

        // strand-ambiguous: A/T (65+84=149) or C/G (67+71=138)
        int s = (int)g1 + (int)g2;
        if (s == 138 || s == 149) {
            ++acount;
            if (amb == 0) continue;
        }

        // strtod with endptr: reject non-numeric / trailing garbage (LDAK uses
        // sscanf "%lf%c"!=1). Without this a bad Z would become 0.0 -> chi=0 ->
        // QC silently rewrites it to 1e-6.
        char* end;
        double z = std::strtod(tok[cZ], &end);
        if (end == tok[cZ] || *end != '\0' || !std::isfinite(z))
            die(std::string("non-numeric Z for predictor ") + tok[cPred] + " (" + tok[cZ] + ")");
        double n = std::strtod(tok[cN], &end);
        if (end == tok[cN] || *end != '\0' || !std::isfinite(n) || n < 0.5)
            die(std::string("bad sample size for predictor ") + tok[cPred] + " (" + tok[cN] + ")");

        double chi = z * z;
        int sign = (z >= 0) ? 1 : -1;
        double rho = sign * std::sqrt(chi / (chi + n));
        if (t1 != g1) rho = -rho;   // orient to tagging A1

        S.snss[j2]  = n;
        S.schis[j2] = chi * scaling;
        S.srhos[j2] = rho;
        seen[j2] = 1;
        ++S.n_matched;
    }
    std::free(buf);
    std::fclose(f);

    if (S.n_matched == 0) die("summary file " + path + " matches none of the tagging predictors");
    if (verbose)
        std::fprintf(stderr, "[%s] matched %d / %d tagging predictors "
                     "(inconsistent=%d ambiguous=%d long=%d)\n",
                     path.c_str(), S.n_matched, T.n, ccount, acount, lcount);
    return S;
}

}  // namespace gs
