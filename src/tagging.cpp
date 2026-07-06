#include "tagging.hpp"
#include "common.hpp"
#include <cstdio>
#include <cstring>

namespace gs {

Tagging read_tagfile(const std::string& path) {
    std::FILE* f = std::fopen(path.c_str(), "r");
    if (!f) die("cannot open tagfile " + path);

    // Read every non-empty line into memory (~33 MB for UKB chip — fine).
    std::vector<std::string> lines;
    {
        char* buf = nullptr; size_t cap = 0; ssize_t len;
        while ((len = getline(&buf, &cap, f)) != -1) {
            while (len > 0 && (buf[len-1]=='\n' || buf[len-1]=='\r')) buf[--len] = '\0';
            if (len > 0) lines.emplace_back(buf, (size_t)len);
        }
        std::free(buf);
    }
    std::fclose(f);
    if (lines.size() < 3) die("tagfile " + path + " too short");

    std::vector<char*> tok;
    Tagging T;

    // ---- header: num_parts = ncols - 9 ----
    {
        std::string h = lines[0];
        int nc = split_ws(&h[0], tok);
        if (nc < 10 || std::strcmp(tok[0], "Predictor") != 0)
            die("tagfile must begin with header line starting \"Predictor\"");
        T.num_parts = nc - 9;
        for (int q = 0; q < T.num_parts; ++q) T.catlabels.emplace_back(tok[9 + q]);
    }

    const int P = T.num_parts;
    const int total = (int)lines.size();
    T.n = total - 2 - P;                 // data rows
    if (T.n <= 0) die("tagfile has no predictor rows");

    T.preds.resize(T.n);
    T.al1.resize(T.n); T.al2.resize(T.n);
    T.stags.resize(T.n);
    T.svars.assign(P, std::vector<double>(T.n));
    T.ssums.assign(P, std::vector<double>(P + 2, 0.0));

    // ---- predictor rows ----
    for (int j = 0; j < T.n; ++j) {
        std::string& ln = lines[1 + j];
        int nc = split_ws(&ln[0], tok);
        if (nc < 9 + P)
            die("tagfile row " + std::to_string(j + 2) + " has too few columns");
        T.preds[j] = tok[0];
        T.al1[j]   = tok[1][0];
        T.al2[j]   = tok[2][0];
        T.stags[j] = parse_double(tok[4], "tagfile Tagging column");   // reject bad numeric
        for (int q = 0; q < P; ++q)
            T.svars[q][j] = parse_double(tok[9 + q], "tagfile category svar");
    }

    // ---- P footer rows: "...to each category" <P values> ----
    for (int q2 = 0; q2 < P; ++q2) {
        std::string& ln = lines[1 + T.n + q2];
        int nc = split_ws(&ln[0], tok);
        if (nc < P) die("tagfile footer row malformed");
        // values are the last P tokens
        for (int q = 0; q < P; ++q)
            T.ssums[q][q2] = parse_double(tok[nc - P + q], "tagfile footer value");
    }

    // ---- last row: "There are <m> ... predictors" <P values> ----
    {
        std::string& ln = lines[total - 1];
        int nc = split_ws(&ln[0], tok);
        if (nc < 3 || std::strcmp(tok[0], "There") != 0)
            die("tagfile final row should start with \"There\" (remake with current LDAK)");
        T.tag_count = parse_double(tok[2], "tagfile predictor count");   // m
        for (int q = 0; q < P; ++q)
            T.ssums[q][P] = parse_double(tok[nc - P + q], "tagfile final-row value");
    }

    // proportion of predictors in each category (matches read_tagfile L1588)
    for (int q = 0; q < P; ++q)
        T.ssums[q][P + 1] = T.tag_count > 0 ? T.ssums[q][P] / T.tag_count : 0.0;

    return T;
}

}  // namespace gs
