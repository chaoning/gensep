#!/usr/bin/env bash
# End-to-end GenSep example: type 1 vs type 2 diabetes from FinnGen public summary
# statistics (release R13, GRCh38) with the ready-made Finnish tagging file.
# Runs in the current directory; needs ~2 GB of disk for the two raw FinnGen files.
#
#   bash run_finngen_t1d_t2d.sh [path/to/gensep]     (default: ./gensep)
#
# Endpoints and counts come from the FinnGen R13 manifest
# (https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv):
#   T1D_WIDE  Type 1 diabetes, wide definition           11,197 cases / 396,409 controls
#   T2D       Type 2 diabetes, definitions combined      89,727 cases / 396,292 controls
set -euo pipefail
GENSEP=${1:-./gensep}
HERE=$(cd "$(dirname "$0")" && pwd)
FG=https://storage.googleapis.com/finngen-public-data-r13/summary_stats

# 1. FinnGen summary statistics (~0.8 GB each) and the Finnish tagging file (~20 MB)
[ -s finngen_R13_T1D_WIDE.gz ] || wget -q --show-progress $FG/finngen_R13_T1D_WIDE.gz
[ -s finngen_R13_T2D.gz ]      || wget -q --show-progress $FG/finngen_R13_T2D.gz
[ -s tag.HAPMAP.FIN.tagging ]  || { wget -q --show-progress \
    https://github.com/chaoning/gensep/releases/download/tagging-v1/tag.HAPMAP.FIN.tagging.gz
    gunzip tag.HAPMAP.FIN.tagging.gz; }

# 2. Convert to GenSep .summaries (Predictor A1 A2 Z n); ~3 min per file
python3 "$HERE/finngen_to_summaries.py" --in finngen_R13_T1D_WIDE.gz \
    --cases 11197 --controls 396409 --tagfile tag.HAPMAP.FIN.tagging --out t1d.summaries
python3 "$HERE/finngen_to_summaries.py" --in finngen_R13_T2D.gz \
    --cases 89727 --controls 396292 --tagfile tag.HAPMAP.FIN.tagging --out t2d.summaries

# 3. GenSep with a 200-block jackknife (~30 s single-threaded)
#    K = population prevalences used in the paper (Supp. Table 1); P = cases / (cases + controls)
"$GENSEP" --se-method jackknife \
    --tagfile tag.HAPMAP.FIN.tagging \
    --summary t1d.summaries --summary2 t2d.summaries \
    --K1 0.005 --K2 0.07 --P1 0.02747 --P2 0.18462 \
    --num-blocks 200 --out finngen_t1d_t2d

echo; echo "=== finngen_t1d_t2d.gensep ==="; cat finngen_t1d_t2d.gensep
echo; echo "Compare with $HERE/finngen_t1d_t2d.expected.gensep"
