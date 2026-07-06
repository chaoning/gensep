// gensep — common math, data structures, and small helpers.
// Math (normal pdf/cdf/inv) uses the C++ standard library (exact erf/erfc) to
// MATCH scipy.stats.norm used in step5_metrics.py — the derived-quantity point
// estimates stay defined by step5, gensep only adds the block-jackknife SE,
// so both must use the same normal functions.
#pragma once
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cerrno>

namespace gs {

constexpr double INV_SQRT_2PI = 0.3989422804014327;  // 1/sqrt(2*pi)
constexpr double INV_SQRT2    = 0.7071067811865476;  // 1/sqrt(2)
constexpr double MISSING      = -9999.0;

inline double normal_pdf(double x) { return INV_SQRT_2PI * std::exp(-0.5 * x * x); }
inline double normal_cdf(double x) { return 0.5 * std::erfc(-x * INV_SQRT2); }

// Standard-normal inverse CDF: Acklam initial guess + one Halley refinement
// using std::erfc (matches scipy.stats.norm.ppf to ~1e-15).
inline double normal_inv(double p) {
    if (p <= 0.0) return -INFINITY;
    if (p >= 1.0) return  INFINITY;
    static const double a[6] = {-3.969683028665376e+01, 2.209460984245205e+02,
                                -2.759285104469687e+02, 1.383577518672690e+02,
                                -3.066479806614716e+01, 2.506628277459239e+00};
    static const double b[5] = {-5.447609879822406e+01, 1.615858368580409e+02,
                                -1.556989798598866e+02, 6.680131188771972e+01,
                                -1.328068155288572e+01};
    static const double c[6] = {-7.784894002430293e-03, -3.223964580411365e-01,
                                -2.400758277161838e+00, -2.549732539343734e+00,
                                 4.374664141464968e+00,  2.938163982698783e+00};
    static const double d[4] = { 7.784695709041462e-03,  3.224671290700398e-01,
                                 2.445134137142996e+00,  3.754408661907416e+00};
    const double plow = 0.02425, phigh = 1.0 - plow;
    double x, q, r;
    if (p < plow) {
        q = std::sqrt(-2.0 * std::log(p));
        x = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
    } else if (p <= phigh) {
        q = p - 0.5; r = q * q;
        x = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
            (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0);
    } else {
        q = std::sqrt(-2.0 * std::log(1.0 - p));
        x = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
             ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
    }
    // one Halley step
    double e = normal_cdf(x) - p;
    double u = e / normal_pdf(x);
    x = x - u / (1.0 + 0.5 * x * u);
    return x;
}

// Delete-one-block jackknife SE EXACTLY as LDAK (sumfuns.c:947 / 1570):
//   var = (B-1) * ( mean(theta^2) - mean(theta)^2 );  SE = sqrt(var)
// Computed over ALL B blocks with no filtering, so a single NaN block (e.g. a
// Cor block where Her1*Her2<0) propagates to SE=NaN, matching LDAK's "-nan".
// Use this for the sum-hers / sum-cors outputs that must be a drop-in for LDAK.
inline double jackknife_se_ldak(const std::vector<double>& theta) {
    const int B = (int)theta.size();
    if (B < 2) return NAN;
    double sum = 0.0, sumsq = 0.0;
    for (double t : theta) { sum += t; sumsq += t * t; }
    double mean = sum / B;
    double var = (B - 1) * (sumsq / B - mean * mean);
    return std::sqrt(var);          // NaN/neg propagate to NaN, as in LDAK
}

// Pair-deleting variant: skips non-finite blocks and reports how many remained.
// Used ONLY for gensep's derived quantities (VS / h2cc / AUC), which are new
// (no LDAK equivalent) — degenerate blocks are dropped on purpose, see gensep.cpp.
inline double jackknife_se(const std::vector<double>& theta, int* n_used = nullptr) {
    double sum = 0.0, sumsq = 0.0; int B = 0;
    for (double t : theta) {
        if (std::isfinite(t)) { sum += t; sumsq += t * t; ++B; }
    }
    if (n_used) *n_used = B;
    if (B < 2) return NAN;
    double mean = sum / B;
    double var = (B - 1) * (sumsq / B - mean * mean);
    return var > 0 ? std::sqrt(var) : 0.0;
}

// ---- fatal error helper ----
[[noreturn]] inline void die(const std::string& msg) {
    std::fprintf(stderr, "Error: %s\n", msg.c_str());
    std::exit(1);
}

// ---- strict numeric parsing (shared by the CLI and all file readers) ----
// Reject empty / non-numeric / trailing-garbage / non-finite input instead of
// silently yielding 0 (std::strtod) or aborting on an exception (std::stod).
// `what` names the field for the error; kept as const char* so the hot tagfile
// loop builds no std::string on the success path.
inline double parse_double(const char* s, const char* what) {
    char* end;
    double v = std::strtod(s, &end);
    if (end == s || *end != '\0' || !std::isfinite(v))
        die(std::string("bad numeric value for ") + what + " ('" + (s ? s : "") + "')");
    return v;
}
inline double parse_double(const std::string& s, const char* what) {
    return parse_double(s.c_str(), what);
}
inline long parse_long(const std::string& s, const char* what) {
    const char* c = s.c_str();
    char* end;
    errno = 0;
    long v = std::strtol(c, &end, 10);
    if (end == c || *end != '\0' || errno == ERANGE)
        die(std::string("bad integer value for ") + what + " ('" + s + "')");
    return v;
}

// Fast in-place whitespace tokenizer: writes NULs into `s`, fills `out` with
// pointers to each token. Returns the token count.
inline int split_ws(char* s, std::vector<char*>& out) {
    out.clear();
    char* p = s;
    while (*p) {
        while (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r') ++p;
        if (!*p) break;
        out.push_back(p);
        while (*p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r') ++p;
        if (*p) *p++ = '\0';
    }
    return (int)out.size();
}

}  // namespace gs
