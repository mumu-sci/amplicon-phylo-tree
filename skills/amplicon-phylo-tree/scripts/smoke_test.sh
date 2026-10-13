#!/usr/bin/env bash
# The end-to-end smoke test, as a script rather than a workflow body.
#
# Why: this is the only check that proves the documented baseline is still
# reproducible (20 taxa in -> 15 kept, 556 alignment columns, 46.04% gaps).  A
# workflow file needs the `workflow` token scope to edit, which not every
# contributor has; a plain script does not.  So the workflow calls this file and
# this file is edited with ordinary commits.
#
# Requires: python3, mafft, fasttree.
# Exit: 0 = the baseline reproduced, 1 = it did not, 2 = a required tool is missing.
set -euo pipefail

if ! command -v mafft >/dev/null || ! command -v fasttree >/dev/null; then
  echo "::warning::MAFFT or FastTree missing - the end-to-end smoke test did NOT run"
  exit 2
fi

# Some environments only have `python3`, some only `python`.  The scripts are
# standard-library-only Python 3, so either is fine.
PY="$(command -v python3 || command -v python)"
if [ -z "$PY" ]; then
  echo "::error::no python3 on PATH"
  exit 2
fi

# The skill folder is this script's parent; every path below is relative to it.
HERE="$(cd "$(dirname "$0")/.." && pwd)"
cd "$HERE"
P="assets/example-mini"
W="$(mktemp -d)"
trap 'rm -rf "$W"' EXIT

echo "--- probe ---"
"$PY" scripts/probe_tools.py --require mafft,tree

echo "--- taxon rule ---"
"$PY" scripts/define_taxa.py --table "$P/asv_table.tsv" --rule 'sum>=10' \
  --fasta "$P/rep_set.fasta" --out-ids "$W/kept.ids" \
  --out-fasta "$W/kept.fasta" --log "$W/taxa.log"
test "$(grep -c '^>' "$W/kept.fasta")" = "15" \
  || { echo "baseline BROKEN: expected 15 kept sequences"; exit 1; }

echo "--- input audit: an unchecked gate must FAIL without a waiver ---"
set +e
"$PY" scripts/audit_inputs.py --fasta "$W/kept.fasta" --label smoke \
  > "$W/audit_fail.txt" 2>&1
rc=$?
set -e
if [ "$rc" -eq 0 ]; then
  echo "baseline BROKEN: the audit passed with no reference database and no waiver"
  exit 1
fi
echo "    ok: exited $rc as expected (fail-closed)"

echo "--- input audit with a recorded waiver ---"
"$PY" scripts/audit_inputs.py --fasta "$W/kept.fasta" --label smoke \
  --accept-and-record --override-reason "smoke test: no reference DB available" \
  --json "$W/audit.json" > "$W/audit.txt" 2>&1

echo "--- alignment ---"
mafft --auto --quiet --thread 2 "$W/kept.fasta" > "$W/aln.fasta"

echo "--- G2 ---"
"$PY" scripts/aln_qc.py --aln "$W/aln.fasta" --fasta-raw "$W/kept.fasta" \
  --marker-class rrna --json "$W/aln.json"
"$PY" - "$W/aln.json" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
assert d["columns"] == 556, "columns: %r" % d["columns"]
assert abs(d["gap_fraction"] - 0.4604) < 5e-4, "gap fraction: %r" % d["gap_fraction"]
print("    baseline reproduced: %d columns, %.2f%% gaps"
      % (d["columns"], 100 * d["gap_fraction"]))
PY

echo "--- tree + G3 ---"
fasttree -nt -gtr -gamma "$W/aln.fasta" > "$W/tree.nwk" 2> "$W/tree.log"
"$PY" scripts/verify_tree.py --aln "$W/aln.fasta" --seqtype nt --tool fasttree \
  --tool-log "$W/tree.log" --ids "$W/kept.ids" --unrooted "$W/tree.nwk" \
  --json "$W/verify.json"

echo
echo "SMOKE TEST PASSED: the documented baseline still reproduces."
