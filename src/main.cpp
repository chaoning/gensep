// gensep — genetic separation from GWAS summary statistics.
// Single purpose: given a tagging file + two traits' summaries + prevalences/
// sample fractions, output h1, h2, rg, VS, case-case h2 (formula) and the
// upper-limit AUC, each with a block-jackknife SE.
//
// Method (fused, option A): one common SNP set + shared 200 blocks; per
// leave-one-block-out h1/h2 are estimated with the SumHer model (cept=0) and rg
// with the cross-trait model (cept=1), so the point estimate and the jackknife SE
// use the same estimator on the same blocks.
#include <cstdio>
#include <cstring>
#include <string>
#include <map>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "common.hpp"
#include "tagging.hpp"
#include "summary.hpp"
#include "qc.hpp"
#include "gensep.hpp"

using namespace gs;

static std::map<std::string, std::string> parse_opts(int argc, char** argv) {
    std::map<std::string, std::string> o;
    for (int i = 1; i < argc; ++i) {
        if (std::strncmp(argv[i], "--", 2) == 0 && i + 1 < argc) { o[argv[i] + 2] = argv[i + 1]; ++i; }
    }
    return o;
}
static std::string need(std::map<std::string, std::string>& o, const std::string& k) {
    auto it = o.find(k);
    if (it == o.end()) die("missing required --" + k);
    return it->second;
}
// Strict numeric option accessors (clean die on bad input, never std::stod abort).
static double opt_d(std::map<std::string, std::string>& o, const char* k) {
    return parse_double(need(o, k), k);
}
static long opt_l(std::map<std::string, std::string>& o, const char* k, long def) {
    return o.count(k) ? parse_long(o[k], k) : def;
}
// Shared K/P domain check — identical for every --se-method path.
static void check_KP(double K1, double K2, double P1, double P2) {
    if (!(K1 > 0 && K1 < 1 && K2 > 0 && K2 < 1)) die("--K1,--K2 must be in (0,1)");
    if (!(P1 > 0 && P1 < 1 && P2 > 0 && P2 < 1)) die("--P1,--P2 must be in (0,1)");
}
// Optional per-subtype PRS case/control AUC (--auc1/--auc2): both-or-neither, in
// (0.5, 0.9999). The upper bound matches auc_to_corr_liab()'s domain in gensep.cpp — an
// auc >= 0.9999 would pass here but yield NaN internally (Rsq nan -> prs_auc NA).
// When present, gensep additionally reports the finite-PRS case-case AUC (point-only).
static bool parse_auc(std::map<std::string, std::string>& o, double& auc1, double& auc2) {
    int n = (int)o.count("auc1") + (int)o.count("auc2");
    if (n == 0) return false;
    if (n != 2) die("--auc1 and --auc2 must be given together");
    auc1 = opt_d(o, "auc1"); auc2 = opt_d(o, "auc2");
    if (!(auc1 > 0.5 && auc1 < 0.9999 && auc2 > 0.5 && auc2 < 0.9999))
        die("--auc1,--auc2 (PRS case/control AUC) must be in (0.5, 0.9999)");
    return true;
}

static void usage() {
    std::fprintf(stderr,
        "gensep --se-method <jackknife|mc|delta|none> ...\n"
        "\n"
        "  --se-method jackknife    — from summary statistics (SumHer + fused block-jackknife)\n"
        "    gensep --se-method jackknife --tagfile T --summary S1 --summary2 S2 \\\n"
        "           --K1 <prev1> --K2 <prev2> --P1 <casefrac1> --P2 <casefrac2> \\\n"
        "           [--num-blocks 200] [--max-threads 1] [--cutoff 0.01] --out PREFIX\n"
        "\n"
        "  --se-method mc|delta|none — from given point estimates (h2_obs, rg, K, P)\n"
        "    gensep --se-method mc --h1 <h2obs1> --h2 <h2obs2> --rg <rg> \\\n"
        "           --K1 <prev1> --K2 <prev2> --P1 <casefrac1> --P2 <casefrac2> \\\n"
        "           --se-h1 s --se-h2 s --se-rg s [--num-draws 100000] [--seed 1] --out PREFIX\n"
        "    mc = Monte-Carlo, delta = finite-difference propagation (both need all three\n"
        "    --se-*, (h1,h2,rg) treated as independent); none = point estimates only, no SE\n"
        "    (do NOT pass --se-*).\n"
        "\n"
        "--se-method is required (no default).\n"
        "\n"
        "  Optional (any mode): --auc1 <v> --auc2 <v>  per-subtype PRS case/control AUC.\n"
        "    Adds the finite-PRS case-case AUC (prs_auc, prs_auc_lo, h2cc_prs, prs_eff),\n"
        "    point-only (no SE), computed from the point hsq*_liab/rg + K.\n"
        "\n"
        "Writes PREFIX.gensep: hsq1/2 (obs+liab), rg, VS, h2cc, auc, auc_lo, each with SE\n"
        "(+ prs_auc/prs_auc_lo/h2cc_prs/prs_eff when --auc1/--auc2 are given).\n");
}

static int run_jackknife(std::map<std::string, std::string>& o) {
    std::string tagfile = need(o, "tagfile");
    std::string sum1 = need(o, "summary");
    std::string sum2 = need(o, "summary2");
    std::string out  = need(o, "out");
    double K1 = opt_d(o, "K1"), K2 = opt_d(o, "K2");
    double P1 = opt_d(o, "P1"), P2 = opt_d(o, "P2");
    check_KP(K1, K2, P1, P2);
    long B = opt_l(o, "num-blocks", 200);
    if (B < 2) die("--num-blocks must be >= 2");
    double auc1 = 0, auc2 = 0; bool have_auc = parse_auc(o, auc1, auc2);
    double cutoff = 0.0;
    if (o.count("cutoff")) {
        cutoff = opt_d(o, "cutoff");
        if (!(cutoff > 0 && cutoff < 0.5)) die("--cutoff must be in (0, 0.5)");
    }

    auto t_read = std::chrono::steady_clock::now();
    Tagging T = read_tagfile(tagfile);
    SummaryAligned S1 = read_sumsfile(sum1, T, /*amb=*/0);
    SummaryAligned S2 = read_sumsfile(sum2, T, /*amb=*/0);
    PairData D = qc_pair(T, S1, S2, cutoff);
    std::fprintf(stderr, "Read tagging + summaries and QC: %.1f s\n",
                 std::chrono::duration<double>(std::chrono::steady_clock::now() - t_read).count());

    SepResult r = gene_sep_fused(D, K1, K2, P1, P2, (int)B, have_auc, auc1, auc2);
    write_gensep(out, r);

    std::printf("hsq1 obs %.4f(%.4f) liab %.4f(%.4f)  hsq2 obs %.4f(%.4f) liab %.4f(%.4f)  "
                "rg %.4f(%.4f)  VS %.4f(%.4f)  h2cc %.4f(%.4f)  auc %.4f(%.4f)  [n=%d/%d] -> %s.gensep\n",
                r.hsq1_obs, r.hsq1_obs_se, r.hsq1_liab, r.hsq1_liab_se,
                r.hsq2_obs, r.hsq2_obs_se, r.hsq2_liab, r.hsq2_liab_se,
                r.rg, r.rg_se, r.VS, r.VS_se, r.h2cc, r.h2cc_se,
                r.auc_exact, r.auc_exact_se, r.B_used, r.B_used_auc, out.c_str());
    return 0;
}

static int run_point(std::map<std::string, std::string>& o, SeMethod method) {
    double h1 = opt_d(o, "h1"), h2 = opt_d(o, "h2"), rg = opt_d(o, "rg");
    double K1 = opt_d(o, "K1"), K2 = opt_d(o, "K2");
    double P1 = opt_d(o, "P1"), P2 = opt_d(o, "P2");
    std::string out = need(o, "out");
    check_KP(K1, K2, P1, P2);

    // Input SEs are gated by --se-method: mc/delta require all three; none forbids them.
    int nse = (int)o.count("se-h1") + (int)o.count("se-h2") + (int)o.count("se-rg");
    bool have_se;
    if (method == SE_NONE) {
        if (nse != 0) die("--se-method none takes no --se-h1/--se-h2/--se-rg "
                          "(use --se-method mc or delta to propagate them)");
        have_se = false;
    } else {
        if (nse != 3) die("--se-method mc|delta needs all of --se-h1, --se-h2, --se-rg "
                          "(or use --se-method none for point estimates only)");
        have_se = true;
    }
    double se1 = 0, se2 = 0, ser = 0;
    if (have_se) {
        se1 = opt_d(o, "se-h1"); se2 = opt_d(o, "se-h2"); ser = opt_d(o, "se-rg");
        if (se1 < 0 || se2 < 0 || ser < 0) die("--se-h1/--se-h2/--se-rg must be >= 0");
    }
    long ndraw = opt_l(o, "num-draws", 100000);
    if (method == SE_MC && ndraw < 2) die("--num-draws must be >= 2");
    long seed = opt_l(o, "seed", 1);
    if (seed < 0) die("--seed must be >= 0");
    double auc1 = 0, auc2 = 0; bool have_auc = parse_auc(o, auc1, auc2);

    SepResult r = gene_sep_point(h1, h2, rg, K1, K2, P1, P2, have_se, se1, se2, ser,
                                 method, ndraw, (unsigned long long)seed,
                                 have_auc, auc1, auc2);
    write_gensep(out, r);

    std::printf("hsq1 obs %.4f(%.4f) liab %.4f(%.4f)  hsq2 obs %.4f(%.4f) liab %.4f(%.4f)  "
                "rg %.4f(%.4f)  VS %.4f(%.4f)  h2cc %.4f(%.4f)  auc %.4f(%.4f)  [n=%d/%d] -> %s.gensep\n",
                r.hsq1_obs, r.hsq1_obs_se, r.hsq1_liab, r.hsq1_liab_se,
                r.hsq2_obs, r.hsq2_obs_se, r.hsq2_liab, r.hsq2_liab_se,
                r.rg, r.rg_se, r.VS, r.VS_se, r.h2cc, r.h2cc_se,
                r.auc_exact, r.auc_exact_se, r.B_used, r.B_used_auc, out.c_str());
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) { usage(); return 1; }
    std::string a1 = argv[1];
    if (a1 == "-h" || a1 == "--help") { usage(); return 0; }

    auto o = parse_opts(argc, argv);

    // Thread count for the OpenMP block-jackknife (default 1 = single-threaded). No effect
    // in a non-OpenMP build (make OMP=0) or for the point modes, which are not parallel.
    long max_threads = opt_l(o, "max-threads", 1);
    if (max_threads < 1) die("--max-threads must be >= 1");
#ifdef _OPENMP
    omp_set_num_threads((int)max_threads);
#endif

    auto it = o.find("se-method");
    if (it == o.end()) die("missing required --se-method (jackknife|mc|delta|none)");
    const std::string& sm = it->second;
    if (sm == "jackknife") return run_jackknife(o);
    if (sm == "mc")        return run_point(o, SE_MC);
    if (sm == "delta")     return run_point(o, SE_DELTA);
    if (sm == "none")      return run_point(o, SE_NONE);
    die("--se-method must be 'jackknife', 'mc', 'delta', or 'none'");
}
