#!/usr/bin/env bash
# Minimal regression tests for gensep. No external data needed: the point modes
# run stand-alone, and the tagfile test builds a tiny malformed fixture inline.
# Usage:  bash tests/run_tests.sh [path/to/gensep]      (default: ./gensep)
set -u

BIN=${1:-./gensep}
[ -x "$BIN" ] || { echo "gensep binary not found/executable: $BIN (run 'make' first)"; exit 2; }
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

pass=0; fail=0
ok()  { printf '  \033[32mPASS\033[0m %s\n' "$1"; pass=$((pass+1)); }
no()  { printf '  \033[31mFAIL\033[0m %s\n' "$1"; fail=$((fail+1)); }

# expect_ok  "desc"  cmd...     -> command must exit 0
expect_ok()  { local d=$1; shift; if "$@" >/dev/null 2>&1; then ok "$d"; else no "$d (expected success, got $?)"; fi; }
# expect_err "desc"  cmd...     -> command must exit non-zero (clean die, not crash)
expect_err() { local d=$1; shift; if "$@" >/dev/null 2>&1; then no "$d (expected failure, got 0)";
               else local rc=$?; if [ "$rc" -ge 128 ]; then no "$d (crashed rc=$rc, not a clean die)"; else ok "$d"; fi; fi; }

# point-mode common args
PT="--h1 0.10 --h2 0.15 --rg 0.30 --K1 0.01 --K2 0.02 --P1 0.5 --P2 0.5"
SE="--se-h1 0.02 --se-h2 0.03 --se-rg 0.10"

echo "point modes:"
expect_ok  "point none"              $BIN --se-method none  $PT --out "$TMP/none"
expect_ok  "point mc"                $BIN --se-method mc    $PT $SE --num-draws 5000 --out "$TMP/mc"
expect_ok  "point delta"             $BIN --se-method delta $PT $SE --out "$TMP/delta"

echo "point-estimate correctness (VS should be 0.683101):"
VS=$(awk '$1=="VS"{print $2}' "$TMP/none.gensep")
if awk "BEGIN{exit !(($VS-0.683101)^2 < 1e-8)}"; then ok "VS point value ($VS)"; else no "VS point value ($VS != 0.683101)"; fi

echo "mc determinism (same seed -> identical file):"
$BIN --se-method mc $PT $SE --seed 42 --out "$TMP/d1" >/dev/null 2>&1
$BIN --se-method mc $PT $SE --seed 42 --out "$TMP/d2" >/dev/null 2>&1
if diff -q "$TMP/d1.gensep" "$TMP/d2.gensep" >/dev/null; then ok "seeded runs identical"; else no "seeded runs differ"; fi

echo "input validation:"
expect_err "bad numeric --h1 foo"    $BIN --se-method none --h1 foo --h2 0.15 --rg 0.3 --K1 0.01 --K2 0.02 --P1 0.5 --P2 0.5 --out "$TMP/x"
expect_err "negative --se-h1"        $BIN --se-method mc    $PT --se-h1 -0.02 --se-h2 0.03 --se-rg 0.1 --out "$TMP/x"
expect_err "--num-draws 0"           $BIN --se-method mc    $PT $SE --num-draws 0 --out "$TMP/x"
expect_err "--num-draws non-int"     $BIN --se-method mc    $PT $SE --num-draws 3.5 --out "$TMP/x"
expect_err "K out of range"          $BIN --se-method none  --h1 0.1 --h2 0.15 --rg 0.3 --K1 1.5 --K2 0.02 --P1 0.5 --P2 0.5 --out "$TMP/x"
expect_err "missing --se-method"     $BIN $PT --out "$TMP/x"
expect_err "unknown --se-method"     $BIN --se-method bogus $PT --out "$TMP/x"
expect_err "none + stray --se-h1"    $BIN --se-method none  $PT --se-h1 0.02 --out "$TMP/x"
expect_err "mc missing some --se"    $BIN --se-method mc    $PT --se-h1 0.02 --out "$TMP/x"

echo "finite-PRS AUC (--auc1/--auc2):"
expect_ok  "prs with --auc"          $BIN --se-method none $PT --auc1 0.70 --auc2 0.65 --out "$TMP/prs"
# point value cross-checked against the Python compute_case_case_auc_prs reference
PA=$(awk '$1=="prs_auc"{print $2}' "$TMP/prs.gensep")
if awk "BEGIN{exit !(($PA-0.673704)^2 < 1e-8)}"; then ok "prs_auc value ($PA)"; else no "prs_auc value ($PA != 0.673704)"; fi
$BIN --se-method none $PT --out "$TMP/noauc" >/dev/null 2>&1
if grep -q prs_ "$TMP/noauc.gensep"; then no "no --auc must omit prs rows"; else ok "no --auc omits prs rows"; fi
if grep -q '^prs_eff' "$TMP/prs.gensep"; then ok "all four prs rows present"; else no "prs rows incomplete"; fi
expect_err "only --auc1"             $BIN --se-method none $PT --auc1 0.70 --out "$TMP/x"
expect_err "--auc out of range"      $BIN --se-method none $PT --auc1 0.4 --auc2 0.65 --out "$TMP/x"

echo "jackknife K/P check fires before file I/O (consistency with point mode):"
expect_err "jackknife bad P" $BIN --se-method jackknife --tagfile /no/such --summary a --summary2 b \
                                  --K1 0.01 --K2 0.02 --P1 2.0 --P2 0.5 --out "$TMP/x"

echo "bad tagfile rejected (strict numeric Tagging column):"
# structurally valid single-category tagfile, but Tagging value is non-numeric
printf 'Predictor A1 A2 Neighbours Tagging Info Weight Component Category Base\n' >  "$TMP/bad.tagging"
printf 'rs1 A G 5 XXX 1 1 1 1 0.5\n'                                             >> "$TMP/bad.tagging"
printf 'Categories 100\n'                                                        >> "$TMP/bad.tagging"
printf 'There are 1000 predictors here 100\n'                                    >> "$TMP/bad.tagging"
expect_err "non-numeric tagfile Tagging" $BIN --se-method jackknife --tagfile "$TMP/bad.tagging" \
                                              --summary a --summary2 b --K1 0.01 --K2 0.02 --P1 0.5 --P2 0.5 --out "$TMP/x"

echo
echo "----- $pass passed, $fail failed -----"
[ "$fail" -eq 0 ]
