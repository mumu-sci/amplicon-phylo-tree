#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Define the taxon set by an EXPLICIT, RECORDED rule.

Why this exists: "which sequences went into the tree" is the single most
under-documented decision in amplicon phylogenetics, and it silently changes the
sample size of every downstream metric.  This script forces the rule to be written
as a formula, applies it, and writes the resulting taxon list plus a log that
records the rule verbatim.

The rule is data, not code: nothing here is a house default, and the script never
narrows the taxon set on its own.  Changing the taxon set is a study decision.

Rule grammar (whitespace-insensitive, case-insensitive):

    sum>=10            total abundance across samples >= 10
    sum>5              total abundance > 5
    prev>=2            present (value > 0) in at least 2 samples
    nonzero>=2         same quantity as prev, under the name used in DESIGN 6
    prev_frac>=20      present in at least 20% of samples
    sum_frac>=0.5      carries at least 0.5% of the whole library's abundance

Combine with and / or; ``or`` binds loosest, so ``A and B or C`` means
``(A and B) or C``.

Examples
--------
    # keep ASVs with >=10 total reads and seen in >=2 samples
    define_taxa.py --table asv_table.tsv --rule 'sum>=10 and prev>=2' \\
                   --fasta rep_set.fasta --out-ids kept.ids --out-fasta kept.fasta

    # no table at all: every sequence in the FASTA becomes a taxon
    define_taxa.py --fasta-only rep_set.fasta --out-ids kept.ids

Exit codes: 0 ok, 1 usage/data error, 2 nothing survived the rule.
"""
import argparse
import csv
import json
import os
import re
import sys


def write_json(path, payload):
    """Write a run record as UTF-8 JSON with '\n' line endings."""
    d = os.path.dirname(os.path.abspath(path))
    if d:
        os.makedirs(d, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False)
    print("json written   : %s" % path)

ID_HEADER_HINTS = ("#otu id", "#otu_id", "otu id", "otu_id", "asv", "asv_id",
                   "feature", "featureid", "feature_id", "otu", "id", "taxon",
                   "#otuid", "seqid", "seq_id", "name")

TAXO_HEADER_HINTS = ("taxonomy", "taxon", "lineage", "classification", "tax",
                     "kingdom", "domain", "consensus_lineage")

ZERO_TOKENS = ("", "0", "0.0", "na", "nan", "n/a", "-", "none", ".")


# ----------------------------------------------------------------------------
# table reading
# ----------------------------------------------------------------------------

def looks_numeric(tok):
    try:
        float(tok)
        return True
    except (TypeError, ValueError):
        return False


def split_row(line, delim):
    # use the csv module so that quoted fields survive
    return next(csv.reader([line], delimiter=delim))


def detect_header(rows_preview, delim):
    """Decide whether the first data line is a header, and why.

    Conservative on purpose: a mis-detected header shifts every count by one row,
    which is the kind of error that reaches a paper undetected.
    """
    if not rows_preview:
        return True, "no data at all; assuming a header exists"
    first = rows_preview[0]
    if not first:
        return True, "first line is empty"

    cell0 = first[0].strip().lower()
    cell0_bare = cell0.lstrip("#").strip()
    if any(hint == cell0 or hint == cell0_bare for hint in ID_HEADER_HINTS):
        return True, "first cell %r matches a known id-column name" % first[0]

    body = first[1:]
    if body:
        non_numeric = sum(1 for c in body if not looks_numeric(c.strip()))
        if non_numeric * 2 >= len(body):
            return True, ("%d of %d cells after the first are non-numeric"
                          % (non_numeric, len(body)))
        return False, ("only %d of %d cells after the first are non-numeric, so "
                       "this looks like data" % (non_numeric, len(body)))
    return False, "single-column first line"


def find_last_annotation_col(header, rows_preview, id_col):
    """Return the index of a trailing taxonomy/annotation column, or None."""
    if header:
        for i in range(len(header) - 1, id_col, -1):
            if header[i].strip().lower() in TAXO_HEADER_HINTS:
                return i
    if rows_preview:
        width = len(rows_preview[0])
        if width > id_col + 2:
            last = width - 1
            vals = [r[last].strip() for r in rows_preview if len(r) > last]
            if vals and all(not looks_numeric(v) for v in vals):
                return last
    return None


# ----------------------------------------------------------------------------
# rule parsing
# ----------------------------------------------------------------------------

PRIM_RE = re.compile(
    r"^(sum|prev|nonzero|sum_frac|prev_frac)\s*(>=|<=|>|<|==)\s*"
    r"([0-9]*\.?[0-9]+)$")


def parse_rule(text):
    """Parse the rule into a disjunction of conjunctions of primaries."""
    clauses = []
    for or_part in re.split(r"\bor\b", text, flags=re.IGNORECASE):
        conj = []
        for and_part in re.split(r"\band\b", or_part, flags=re.IGNORECASE):
            tok = and_part.strip()
            if not tok:
                continue
            m = PRIM_RE.match(tok)
            if not m:
                raise ValueError("cannot parse rule term %r" % tok)
            conj.append((m.group(1), m.group(2), float(m.group(3))))
        if conj:
            clauses.append(conj)
    if not clauses:
        raise ValueError("empty rule")
    return clauses


def test_primary(metric, op, thr, total, prev, n_samples, lib_total):
    if metric == "sum":
        val = total
    elif metric in ("prev", "nonzero"):
        # `nonzero>=N` is the presence form documented in DESIGN 6 and
        # parameters.md D3: the number of samples where the taxon is non-zero.
        # It is the same quantity as `prev`, kept as an alias so that both the
        # documented spellings work.
        val = prev
    elif metric == "sum_frac":
        val = 100.0 * total / lib_total if lib_total else 0.0
    elif metric == "prev_frac":
        val = 100.0 * prev / n_samples if n_samples else 0.0
    else:
        raise ValueError(metric)
    return {">=": val >= thr, "<=": val <= thr, ">": val > thr,
            "<": val < thr, "==": val == thr}[op]


def passes(clauses, total, prev, n_samples, lib_total):
    for conj in clauses:
        if all(test_primary(m, o, t, total, prev, n_samples, lib_total)
               for (m, o, t) in conj):
            return True
    return False


def describe_rule(clauses):
    parts = []
    for conj in clauses:
        parts.append(" and ".join("%s %s %g" % (m, o, t) for (m, o, t) in conj))
    return " or ".join("(%s)" % p if len(clauses) > 1 else p for p in parts)


# ----------------------------------------------------------------------------
# FASTA
# ----------------------------------------------------------------------------

def read_fasta(path):
    seqs = []
    name = None
    buf = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(buf)))
                # tolerate labels with trailing description
                name = line[1:].strip().split()[0]
                buf = []
            elif line:
                buf.append(line)
    if name is not None:
        seqs.append((name, "".join(buf)))
    return seqs


def write_fasta(records, path):
    d = os.path.dirname(os.path.abspath(path))
    if d:
        os.makedirs(d, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for name, seq in records:
            fh.write(">%s\n%s\n" % (name, seq))


# ----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Define the taxon set by an explicit rule and record it.")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--table", help="abundance table (tsv or csv)")
    src.add_argument("--fasta-only", metavar="FASTA",
                     help="no table: every sequence in this FASTA is a taxon")

    ap.add_argument("--rule", default=None,
                    help="rule over abundance, e.g. 'sum>=10 and prev>=2'. "
                         "Terms: sum, prev, nonzero (alias of prev), sum_frac, "
                         "prev_frac.  Required with --table; ignored for "
                         "--fasta-only.")
    ap.add_argument("--out-ids", required=True, help="write the taxon id list here")
    ap.add_argument("--fasta", default=None,
                    help="representative sequences to subset from")
    ap.add_argument("--out-fasta", default=None,
                    help="write the subset FASTA here (requires --fasta)")
    ap.add_argument("--format", choices=["auto", "tsv", "csv"], default="auto")
    ap.add_argument("--id-col", type=int, default=0)
    ap.add_argument("--sample-cols", default=None, metavar="START:END",
                    help="0-based half-open column range holding the samples; "
                         "default = every column except the id and any trailing "
                         "annotation column")
    ap.add_argument("--no-header", action="store_true",
                    help="force 'the first line is data'")
    ap.add_argument("--drop-last-col", action="store_true",
                    help="force treating the last column as an annotation string")
    ap.add_argument("--assume-counts", action="store_true",
                    help="proceed even if the table looks like proportions")
    ap.add_argument("--json", default=None,
                    help="write a machine-readable run record here (mode, rule, "
                         "counts retained, and -- for --fasta-only -- an explicit "
                         "statement that no abundance table was used)")
    ap.add_argument("--log", default=None, help="write the processing log here")
    args = ap.parse_args()

    log_lines = []

    def log(msg=""):
        print(msg)
        log_lines.append(msg)

    log("define_taxa -- taxon set by explicit rule")
    log("date/tool : define_taxa.py (amplicon-phylo-tree)")
    log("")

    # ---------- FASTA-only ----------
    if args.fasta_only:
        recs = read_fasta(args.fasta_only)
        ids = [n for n, _ in recs]
        if len(ids) != len(set(ids)):
            log("WARNING: duplicate sequence ids in %s" % args.fasta_only)
        no_table_note = ("this run contains NO abundance table, so the tip set is "
                         "not aligned with any abundance matrix and must not be "
                         "described as if it were")
        log("mode          : fasta-only (no abundance table)")
        log("fasta         : %s" % args.fasta_only)
        log("sequences     : %d" % len(recs))
        log("taxon rule    : none applied -- every sequence with a representative "
            "is a taxon")
        log("table         : NONE")
        log("!! %s" % no_table_note)
        log("taxa retained : %d" % len(ids))
        os.makedirs(os.path.dirname(os.path.abspath(args.out_ids)) or ".",
                    exist_ok=True)
        with open(args.out_ids, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("\n".join(ids) + "\n")
        log("ids written   : %s" % args.out_ids)
        if args.log:
            with open(args.log, "w", encoding="utf-8", newline="\n") as fh:
                fh.write("\n".join(log_lines) + "\n")
        if args.json:
            write_json(args.json, {
                "mode": "fasta-only",
                "fasta": args.fasta_only,
                "abundance_table": None,
                "rule": None,
                "n_sequences": len(recs),
                "n_taxa_retained": len(ids),
                "ids_out": args.out_ids,
                "must_not_claim": no_table_note,
                "note": ("no abundance table was supplied; the tip set is defined by "
                         "sequence presence only"),
            })
        print("define_taxa -- fasta-only")
        print("  sequences      : %d" % len(recs))
        print("  taxa retained  : %d" % len(ids))
        print("  !! NO abundance table: do NOT claim the tip set is aligned with any")
        print("     abundance matrix.  (Recorded in the log and, with --json, in the")
        print("     JSON report.)")
        return 0

    if not args.rule:
        ap.error("--rule is required with --table (or use --fasta-only)")

    # ---------- read table ----------
    delim = "\t"
    if args.format == "csv":
        delim = ","
    elif args.format == "auto":
        with open(args.table, encoding="utf-8", errors="replace") as fh:
            head = fh.readline()
        delim = "," if head.count(",") > head.count("\t") else "\t"
    log("table         : %s" % args.table)
    log("delimiter     : %s" % ("TAB" if delim == "\t" else "COMMA"))

    # preview for header / annotation detection
    preview = []
    with open(args.table, encoding="utf-8", errors="replace", newline="") as fh:
        rdr = csv.reader(fh, delimiter=delim)
        for row in rdr:
            if not row or (len(row) == 1 and not row[0].strip()):
                continue
            preview.append(row)
            if len(preview) >= 12:
                break

    if args.no_header:
        has_header, why = False, "forced off by --no-header"
    else:
        has_header, why = detect_header(preview, delim)
    log("header        : %s  (%s)" % ("yes" if has_header else "no", why))

    data_preview = preview[1:] if has_header else preview
    header = preview[0] if has_header else None

    if args.drop_last_col:
        anno_col = len(preview[0]) - 1
        log("annotation col: %d (forced by --drop-last-col)" % anno_col)
    else:
        anno_col = find_last_annotation_col(header, data_preview, args.id_col)
        log("annotation col: %s" % ("none detected" if anno_col is None
                                   else "%d" % anno_col))

    width = len(preview[0])
    if args.sample_cols:
        a, b = args.sample_cols.split(":")
        s0, s1 = int(a), int(b)
    else:
        s1 = anno_col if anno_col is not None else width
        s0 = args.id_col + 1
    log("sample cols   : [%d, %d)  -> %d samples" % (s0, s1, s1 - s0))
    if s1 <= s0:
        ap.error("empty sample column range")
    log("")

    # ---------- pass 1: totals ----------
    try:
        clauses = parse_rule(args.rule)
    except ValueError as exc:
        # A rule the parser cannot read is a usage error, not a crash: the user
        # has to be told WHICH term failed, and the traceback would not say it.
        log("ERROR: could not parse --rule %r: %s" % (args.rule, exc))
        print("ERROR: could not parse --rule %r: %s" % (args.rule, exc))
        print("  grammar: sum / prev / nonzero / sum_frac / prev_frac, each with "
              ">= <= > < ==, joined by 'and' / 'or'")
        print("  example: --rule 'sum>=10 and prev>=2'")
        return 1
    log("rule (verbatim): %s" % args.rule)
    log("rule (parsed)  : %s" % describe_rule(clauses))

    stats = {}
    order = []
    lib_total = 0
    n_rows = 0
    max_cell = 0.0
    frac_hits = 0
    int_like = 0
    frac_like = 0

    with open(args.table, encoding="utf-8", errors="replace", newline="") as fh:
        rdr = csv.reader(fh, delimiter=delim)
        first = True
        for row in rdr:
            if not row or (len(row) == 1 and not row[0].strip()):
                continue
            if first and has_header:
                first = False
                continue
            first = False
            n_rows += 1
            sid = row[args.id_col].strip()
            total = 0.0
            prev = 0
            for c in range(s0, min(s1, len(row))):
                tok = row[c].strip()
                if tok.lower() in ZERO_TOKENS:
                    continue
                val = float(tok)
                if val:
                    if val == int(val):
                        int_like += 1
                    else:
                        frac_like += 1
                    total += val
                    prev += 1
                    if val > max_cell:
                        max_cell = val
            lib_total += total
            stats[sid] = (total, prev)
            order.append(sid)

    n_samples = s1 - s0
    log("")
    log("rows read      : %d" % n_rows)
    log("unique ids     : %d" % len(stats))
    log("library total  : %g" % lib_total)
    log("largest cell   : %g" % max_cell)

    # is this a count table or a proportion table?
    if n_rows:
        med = sorted(t for t, _ in stats.values())[n_rows // 2]
        log("median row sum : %g" % med)
        if not args.assume_counts and med <= 110 and max_cell <= 110:
            log("")
            log("!! VALUES LOOK LIKE PROPORTIONS, NOT COUNTS.")
            log("   A rule such as 'sum>=10' compares against this scale and may "
                "mean something very different than intended.")
            log("   Use a *_frac rule, or pass --assume-counts if these really are "
                "counts.")
            if args.log:
                os.makedirs(os.path.dirname(os.path.abspath(args.log)) or ".",
                            exist_ok=True)
                with open(args.log, "w", encoding="utf-8", newline="\n") as fh:
                    fh.write("\n".join(log_lines) + "\n")
            return 1

    # ---------- pass 2: apply ----------
    kept = [sid for sid in order
            if passes(clauses, stats[sid][0], stats[sid][1], n_samples, lib_total)]
    kept_set = set(kept)
    reads_kept = sum(stats[s][0] for s in kept)

    log("")
    log("taxa in        : %d" % len(order))
    log("taxa retained  : %d  (%.2f%% of taxa)" % (
        len(kept), 100.0 * len(kept) / len(order) if order else 0.0))
    log("abundance kept : %g  (%.2f%% of library)" % (
        reads_kept, 100.0 * reads_kept / lib_total if lib_total else 0.0))
    log("")
    log("--- report these three numbers in the methods ---")
    log("  ASVs dropped        : %d (of %d)" % (len(order) - len(kept), len(order)))
    log("  richness retained   : %.2f%%"
        % (100.0 * len(kept) / len(order) if order else 0.0))
    log("  reads retained      : %.2f%%"
        % (100.0 * reads_kept / lib_total if lib_total else 0.0))
    log("")
    log("!! IF YOU CHANGE THIS RULE, EVERYTHING DOWNSTREAM MUST BE RE-RUN:")
    log("   alignment -> tree -> rooting -> verification.")
    log("   Do NOT subset an existing tree and reuse it: that silently moves")
    log("   branch lengths and the root position, which is a conclusion-level error.")
    log("   This id list (%s) is the authoritative taxon set for the run."
        % args.out_ids)

    os.makedirs(os.path.dirname(os.path.abspath(args.out_ids)) or ".", exist_ok=True)
    with open(args.out_ids, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(kept) + "\n")
    log("ids written    : %s" % args.out_ids)

    # ---------- optional FASTA subset ----------
    if args.out_fasta or args.fasta:
        if not args.fasta:
            ap.error("--out-fasta requires --fasta")
        if not args.out_fasta:
            ap.error("--fasta requires --out-fasta")
        recs = read_fasta(args.fasta)
        have = dict(recs)
        missing = [s for s in kept if s not in have]
        extra = [n for n, _ in recs if n not in kept_set]
        subset = [(n, have[n]) for n in kept if n in have]
        write_fasta(subset, args.out_fasta)
        log("")
        log("fasta in       : %s  (%d records)" % (args.fasta, len(recs)))
        log("fasta out      : %s  (%d records)" % (args.out_fasta, len(subset)))
        log("ids with no sequence in the FASTA : %d" % len(missing))
        if missing:
            log("   e.g. %s" % ", ".join(missing[:8]))
            log("   !! those taxa cannot appear in the tree; fix the inputs or the rule")
        log("sequences not kept (dropped)      : %d" % len(extra))
        if missing:
            log("")
            log("RESULT: FAIL -- taxon list and FASTA are not consistent.")
            if args.log:
                with open(args.log, "w", encoding="utf-8", newline="\n") as fh:
                    fh.write("\n".join(log_lines) + "\n")
            return 1

    if args.log:
        os.makedirs(os.path.dirname(os.path.abspath(args.log)) or ".", exist_ok=True)
        with open(args.log, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("\n".join(log_lines) + "\n")
        print("log written    : %s" % args.log)

    if args.json:
        write_json(args.json, {
            "mode": "table",
            "table": args.table,
            "abundance_table": args.table,
            "rule": args.rule,
            "rule_parsed": describe_rule(clauses),
            "n_rows_read": n_rows,
            "n_taxa_retained": len(kept),
            "library_total": lib_total,
            "ids_out": args.out_ids,
            "must_not_claim": None,
            "note": ("the taxon rule above defines the tip set; changing it means "
                     "re-running alignment, tree, rooting and verification"),
        })

    if not kept:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
