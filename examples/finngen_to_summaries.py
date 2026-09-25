#!/usr/bin/env python3
"""Convert a FinnGen public summary-statistics file to the GenSep .summaries format.

FinnGen files (GRCh38, tab-separated, gzipped) have the header
    #chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt ...
where beta is the effect of the *alt* allele. FinnGen files carry no per-SNP
sample size, so the case/control counts from the FinnGen manifest are supplied
on the command line and n = cases + controls is written for every SNP.

Output columns: Predictor A1 A2 Z n  (A1 = alt = effect allele, A2 = ref,
Z = beta / sebeta). Variants without an rsID and rsIDs that occur more than
once (multi-allelic sites) are dropped.

Usage:
  python3 finngen_to_summaries.py --in finngen_R13_T2D.gz --cases 89727 --controls 396292 \
      --out t2d.summaries [--tagfile tag.HAPMAP.FIN.tagging]

--tagfile is optional: when given, only predictors present in the tagging file
are written, which shrinks the output from ~1 GB to a few tens of MB. GenSep
uses only those SNPs anyway, so the results are identical either way.
"""
import argparse, gzip, sys


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="src", required=True, help="FinnGen summary-statistics .gz file")
    ap.add_argument("--cases", type=int, required=True, help="number of cases (FinnGen manifest)")
    ap.add_argument("--controls", type=int, required=True, help="number of controls (FinnGen manifest)")
    ap.add_argument("--out", required=True, help="output .summaries file")
    ap.add_argument("--tagfile", help="optional GenSep/LDAK tagging file: keep only its predictors")
    a = ap.parse_args()

    keep = None
    if a.tagfile:
        keep = set()
        with open(a.tagfile) as f:
            next(f)
            for line in f:
                keep.add(line.split(None, 1)[0])
        print(f"[tagfile] {len(keep)} predictors", file=sys.stderr)

    n_total = a.cases + a.controls
    rows, seen, dup = [], set(), set()
    n_in = n_norsid = 0
    with gzip.open(a.src, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(header)}
        for col in ("rsids", "ref", "alt", "beta", "sebeta"):
            if col not in ix:
                sys.exit(f"missing column '{col}' in {a.src}; is this a FinnGen file?")
        for line in f:
            n_in += 1
            t = line.rstrip("\n").split("\t")
            rsid = t[ix["rsids"]].split(",")[0]
            if not rsid.startswith("rs"):
                n_norsid += 1
                continue
            if keep is not None and rsid not in keep:
                continue
            try:
                z = float(t[ix["beta"]]) / float(t[ix["sebeta"]])
            except (ValueError, ZeroDivisionError):
                continue
            if rsid in seen:
                dup.add(rsid)
                continue
            seen.add(rsid)
            rows.append((rsid, t[ix["alt"]].upper(), t[ix["ref"]].upper(), z))

    with open(a.out, "w") as fo:
        fo.write("Predictor\tA1\tA2\tZ\tn\n")
        n_out = 0
        for rsid, a1, a2, z in rows:
            if rsid in dup:
                continue
            fo.write(f"{rsid}\t{a1}\t{a2}\t{z:.6f}\t{n_total}\n")
            n_out += 1
    print(f"[done] read={n_in} no_rsid={n_norsid} duplicated_rsid={len(dup)} "
          f"written={n_out} n={n_total} P={a.cases / n_total:.6f} -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
