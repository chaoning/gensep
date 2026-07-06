// gensep — tagging-file reader. Port of LDAK read_tagfile (filemain.c:1511)
// for the common case (no --alt-tagging / --extract / --tagone).
// Tagging file layout (created by LDAK --calc-tagging):
//   row 1   : Predictor A1 A2 Neighbours Tagging Weight MAF Categories
//             Exp_Heritability <cat_label_1> ... <cat_label_P>   (P = num_parts)
//   n rows  : <pred> <A1> <A2> <Neighbours> <Tagging> <Weight> <MAF>
//             <Categories> <Exp_Heritability> <svar_1> ... <svar_P>
//   P rows  : "The relative contribution of the <cat> to each category" <P values>
//   1 row   : "There are <m> reference <m> regression <m> heritability
//             predictors" <P values>
// => num_parts = ncols(row1) - 9 ;  n = total_rows - 2 - num_parts.
#pragma once
#include <vector>
#include <string>

namespace gs {

struct Tagging {
    int num_parts = 0;
    int n = 0;                                   // number of tagging predictors
    std::vector<std::string> catlabels;          // [num_parts]
    std::vector<std::string> preds;              // [n]
    std::vector<char> al1, al2;                  // [n]
    std::vector<double> stags;                   // [n]  per-SNP LD tagging (col 5)
    std::vector<std::vector<double>> svars;      // [num_parts][n] category cols
    std::vector<std::vector<double>> ssums;      // [num_parts][num_parts+2]
    double tag_count = 0;                        // m (reference/regression count)
};

// Reads the whole tagging file. Exits on format error.
Tagging read_tagfile(const std::string& path);

}  // namespace gs
