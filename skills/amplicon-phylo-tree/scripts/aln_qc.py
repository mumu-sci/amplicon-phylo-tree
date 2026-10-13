#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Quality-check an alignment, and detect the fragmentation failure mode.

Why this exists: a "successful" run can still produce a useless alignment.  The
classic failure is that the aligner cannot relate some sequences to the others, so
it lays them out side by side and the matrix explodes to many times the sequence
length while every sequence becomes almost all gaps.  Numbers like that look fine
in a log ("wrote 128 MB") and produce a tree that is pure noise.

Two things are checked, both against the RAW sequence lengths (pass --fasta-raw):

  width ratio : columns / median raw length   -> above --max-col-ratio, suspect
  gap fraction: gaps / residues               -> above --max-gap, suspect

Note on FASTA wrapping: MAFFT and most aligners wrap output at 60 characters per
line.  Counting columns as "the length of a line" is therefore wrong and gives
absurd gap fractions.  Records are reassembled here before anything is measured.

Structural errors (hard FAIL, never overridable): not rectangular; a codon
alignment whose width is not a multiple of 3; an ILLEGAL character; and a MIXED
ALPHABET -- both nucleotide-only letters (ACGTU) and protein-only letters
(EFILPQZJO*) in the non-gap symbols.  A mixed alphabet has no single substitution
model, so it cannot be carried forward -- that is the "字母表混杂" hard stop the
workflow promises.

--ambiguous-policy decides what an ambiguous character (N ? X) MEANS for the
metrics: 'gap' counts it as a gap character (the default), 'missing' counts it as
missing data (also a gap for the gap fraction), 'drop-col' changes nothing in the
file but writes a LIST OF SUGGESTED COLUMNS whose ambiguous fraction exceeds
--ambiguous-max-frac.  The file is never modified.

Exit codes: 0 = PASS, 1 = FAIL (fragmentation suspected), 2 = input error.
"""
import argparse
import json
import math
import os
import random
import statistics
import sys

GAP_CHARS = "-."
AMBIG_CHARS = "nN?xX"
# Nucleotide-only letters (IUPAC), including the DNA/RNA forms.
NT_ONLY = set("ACGTUacgtu")
NT_IUPAC = set("RYSWKMBDHVryswkmbdhv")
# Standard amino-acid letters that are NOT IUPAC nucleotide ambiguity codes.
AA_ONLY = set("EFILPQZJOefilpqzjo*")
# Recognised gap/blank symbols across alignment formats.
GAP_LIKE = set("-.")


def pct_txt(num, den):
    """"12.34%" for a non-negative fraction; "n/a" when the denominator is 0."""
    if not den:
        return "n/a"
    return "%.2f%%" % (100.0 * num / den)


def quantiles(vals, ps=(0.0, 0.25, 0.5, 0.75, 1.0)):
    """Interpolated quantiles; [] for an empty input."""
    if not vals:
        return []
    s = sorted(vals)
    out = []
    for q in ps:
        if len(s) == 1:
            out.append(float(s[0]))
            continue
        pos = q * (len(s) - 1)
        lo = int(math.floor(pos))
        hi = int(math.ceil(pos))
        if lo == hi:
            out.append(float(s[lo]))
        else:
            out.append(float(s[lo] + (s[hi] - s[lo]) * (pos - lo)))
    return out


DIAGNOSIS = """\
Fragmentation suspected.  Ranked causes, most common first:

 1. The taxon set mixes sequences that are not homologous to each other.  The
    aligner cannot relate them, so it places them side by side.  Check whether
    the table really contains only the intended marker, and whether the taxon
    rule let non-target sequences through.
 2. Two or more different amplicon regions are mixed (e.g. 16S V3-V4 with V4,
    ITS1 with ITS2, two primer sets).  These must never be combined into one
    tree -- split them and build one tree per region.
 3. Protein-coding marker without reading-frame/strand normalisation.  Mixed
    frames and reverse-complemented sequences cannot be aligned.  Run
    frame_strand.py first.
 4. ITS or another marker with extreme length heterogeneity.  Consider splitting
    (e.g. ITSx into ITS1/ITS2) or building one tree per length class.
 5. A single very long or very short outlier dragging the matrix wide.  Look at
    the per-record length spread -- it is printed below.

Once the cause is known, either fix the input, or align in blocks
(blocks.py plan / blocks.py concat), which keeps every taxon while preventing the
unrelated part from distorting the part that is alignable.
"""


def read_clustal(path):
    """Read a Clustal alignment, preserving input order."""
    order = []
    data = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if (not line.strip() or line[0].isspace()
                    or line.upper().startswith("CLUSTAL")
                    or line.upper().startswith("MUSCLE")):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            name, chunk = parts[0], parts[1]
            if name not in data:
                data[name] = []
                order.append(name)
            data[name].append(chunk)
    for name in order:
        yield name, "".join(data[name])


def read_plain(path):
    """Yield (id, sequence) reassembling wrapped FASTA records."""
    name = None
    buf = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(buf)
                name = line[1:].strip().split()[0]
                buf = []
            elif line:
                buf.append(line)
    if name is not None:
        yield name, "".join(buf)


def pick_reader(path, fmt):
    if fmt == "clustal":
        return read_clustal
    if fmt == "fasta":
        return read_plain
    # auto: a Clustal file starts with CLUSTAL / MUSCLE, or has no '>' at all
    with open(path, encoding="utf-8", errors="replace") as fh:
        head = fh.read(2000)
    if head.upper().lstrip().startswith(("CLUSTAL", "MUSCLE")) or ">" not in head:
        return read_clustal
    return read_plain


def median_len(path, fmt="auto"):
    lens = [len(s) for _, s in pick_reader(path, fmt)(path)]
    if not lens:
        return None, 0
    return statistics.median(lens), len(lens)


def main():
    ap = argparse.ArgumentParser(
        description="Check an alignment and detect fragmentation.")
    ap.add_argument("--aln", required=True, help="alignment FASTA")
    ap.add_argument("--fasta-raw", default=None,
                    help="the UNALIGNED sequences this alignment came from; "
                         "enables the width-ratio test")
    ap.add_argument("--max-col-ratio", type=float, default=4.0,
                    help="fail if columns > this x median raw length (default 4)")
    ap.add_argument("--max-gap", type=float, default=None,
                    help="fail if gap fraction > this.  Default depends on "
                         "--marker-class: 0.80 for rrna/organelle, 0.85 for "
                         "functional-nuclear, 0.90 for its/other (ITS and other "
                         "heterogeneous amplicons are gappier by nature, so the "
                         "nucleotide number would make the exception the rule)")
    ap.add_argument("--max-gap-aa", type=float, default=None,
                    help="same, for amino-acid alignments -- they carry more gaps "
                         "than nucleotide alignments by nature (default 0.85, or "
                         "the marker-class default + 0.05)")
    ap.add_argument("--max-gap-codon", type=float, default=None,
                    help="same, for back-translated codon alignments, which are the "
                         "gappiest of all (default 0.90, or the marker-class "
                         "default + 0.10)")
    ap.add_argument("--marker-class",
                    choices=["rrna", "its", "functional-nuclear", "organelle",
                             "other"],
                    default=None,
                    help="sets the default gap thresholds for this marker family "
                         "(see references/parameters.md D2 and marker-types.md).  "
                         "An explicit --max-gap* always wins over this.")
    ap.add_argument("--codon", action="store_true",
                    help="treat the alignment as a back-translated codon alignment: "
                         "require the width to be a multiple of 3 and use "
                         "--max-gap-codon")
    ap.add_argument("--format", choices=["auto", "fasta", "clustal"], default="auto",
                    help="input format.  MAFFT --clustalout produces Clustal, and "
                         "parsing it as FASTA silently miscounts lengths")
    ap.add_argument("--ambiguous-policy", choices=["gap", "missing", "drop-col"],
                    default="gap",
                    help="what an ambiguous character (N, ?, X) MEANS for the "
                         "metrics (references/parameters.md D5).  gap = count it as "
                         "a gap character (default); missing = count it as missing "
                         "data; drop-col = do not change this file, but write the "
                         "list of columns whose ambiguous fraction exceeds "
                         "--ambiguous-max-frac")
    ap.add_argument("--ambiguous-max-frac", type=float, default=0.10,
                    help="with --ambiguous-policy drop-col, the ambiguous fraction "
                         "above which a column is suggested for removal (default "
                         "0.10)")
    ap.add_argument("--outdir", default=None,
                    help="where to write the drop-col list and other side files; "
                         "defaults to the directory of --json (or the current "
                         "directory)")
    ap.add_argument("--accept-and-record", action="store_true",
                    help="downgrade a FAIL on the RATIO metrics to a recorded "
                         "warning.  Structural errors still fail.  Requires "
                         "--override-reason: an override that is not recorded is "
                         "indistinguishable from a pass.")
    ap.add_argument("--override-reason", default=None,
                    help="why it is acceptable to override the ratio gate(s); "
                         "required with --accept-and-record and written into the "
                         "JSON report and the ops log")
    ap.add_argument("--label", default=None, help="label for reports and JSON")
    ap.add_argument("--json", default=None, help="write metrics to this JSON path")
    ap.add_argument("--occupancy", action="store_true",
                    help="also compute per-column occupancy (slower on large files)")
    ap.add_argument("--min-occupancy", type=float, default=0.05,
                    help="with --occupancy, report columns below this occupancy "
                         "(default 0.05)")
    args = ap.parse_args()

    if args.accept_and_record and not args.override_reason:
        ap.error("--accept-and-record requires --override-reason: a waived gate "
                 "must say why, or the waiver is indistinguishable from an "
                 "oversight")

    if not os.path.exists(args.aln):
        print("ERROR: no such alignment: %s" % args.aln)
        return 2

    label = args.label or os.path.basename(args.aln)

    n_seq = 0
    width = None
    residues = 0
    gaps = 0
    gap_char_count = 0
    ambig = 0
    alphabet = set()
    ungapped = []
    ragged = []
    sequences = []
    occ = None

    for name, seq in pick_reader(args.aln, args.format)(args.aln):
        n_seq += 1
        L = len(seq)
        sequences.append((name, seq))
        if width is None:
            width = L
        elif L != width:
            ragged.append((name, L))
        residues += L
        gc = seq.count("-") + seq.count(".")
        ac = sum(seq.count(c) for c in AMBIG_CHARS)
        gaps += gc
        gap_char_count += gc
        ambig += ac
        alphabet.update(seq)
        ungapped.append(L - gc)
        if args.occupancy:
            if occ is None:
                occ = [0] * (width or 0)
            for i, ch in enumerate(seq):
                if i < len(occ) and ch not in GAP_CHARS:
                    occ[i] += 1

    if n_seq == 0:
        print("ERROR: alignment has no sequences")
        return 2

    # ---- ambiguous policy: redefine what the gap metric counts ----
    # 'gap' / 'missing' both make an ambiguous character a hole in the matrix;
    # 'gap' is the historical default (so existing numbers do not move), while
    # 'missing' additionally treats it as missing DATA.  Either way the file is
    # never edited.
    if args.ambiguous_policy in ("gap", "missing"):
        effective_gaps = gaps + (ambig if args.ambiguous_policy == "missing" else 0)
    else:                                     # drop-col
        effective_gaps = gaps
    # The denominator is n_seq x width -- the SAME quantity the header prints -- not
    # the sum of record lengths.  The two agree for a rectangular alignment and
    # diverge silently for a ragged one, which is exactly when the number matters.
    denom_matrix = n_seq * (width or 0)
    denom_residues = residues
    gap_frac = effective_gaps / denom_matrix if denom_matrix else 0.0
    gap_frac_by_residues = effective_gaps / denom_residues if denom_residues else 0.0
    alpha = "".join(sorted(alphabet))
    nongap = "".join(sorted(set(alpha) - set(GAP_CHARS)))

    # crude data-type call, used to sanity-check the substitution model later
    nt_chars = set("ACGTUNacgtunRYSWKMBDHVryswkmbdhv")
    if set(nongap) <= nt_chars and nongap:
        dtype = "nucleotide"
    elif len(set(nongap) - set("ACGTUacgtu")) >= 4:
        dtype = "protein"
    else:
        dtype = "unknown"

    # ---- structural errors: alphabets and illegal symbols ----
    # These are HARD fails (never overridable).  A mixed alphabet has no single
    # substitution model, so no tree built from it can be defended.
    mixed_nt, mixed_aa, illegal = set(), set(), set()
    for _name, seq in sequences:
        syms = set(seq) - GAP_LIKE
        mixed_nt |= (syms & NT_ONLY)
        mixed_aa |= (syms & AA_ONLY)
        illegal |= {c for c in syms
                    if c not in NT_ONLY and c not in NT_IUPAC
                    and c not in AA_ONLY and c not in "Nn?Xx"}
    mixed_alphabet = bool(mixed_nt and mixed_aa)

    print("alignment QC -- %s" % label)
    print("-" * 72)
    print("  file                    : %s" % args.aln)
    print("  sequences               : %d" % n_seq)
    print("  columns                 : %s" % (width if width is not None else "?"))
    print("  residues                : %d" % residues)
    print("  gaps (- and .)          : %d" % gaps)
    print("  ambiguous policy        : %s%s"
          % (args.ambiguous_policy,
             "" if args.ambiguous_policy == "gap" else "  (N ? X counted differently)"))
    print("  gap fraction            : %.4f  (%.2f%%)"
          % (gap_frac, 100 * gap_frac))
    if abs(gap_frac - gap_frac_by_residues) > 1e-9:
        print("     !! the records are NOT all the same length, so the two possible")
        print("        denominators differ: n_seq x width gives the value above")
        print("        (the documented definition); gap chars / total residues")
        print("        would give %.4f.  Fix the alignment before using either."
              % gap_frac_by_residues)
    print("  ambiguous (N ? X ...)   : %d" % ambig)
    print("  alphabet present        : %s" % (alpha if len(alpha) <= 30
                                              else alpha[:30] + "..."))
    print("  non-gap symbols         : %s" % (nongap if len(nongap) <= 30
                                              else nongap[:30] + "..."))
    print("  inferred data type      : %s" % dtype)
    if mixed_alphabet:
        print("  !! MIXED ALPHABET: nucleotide-only %s present together with "
              "protein-only %s" % ("".join(sorted(mixed_nt)),
                                   "".join(sorted(mixed_aa))))
    if illegal:
        print("  !! ILLEGAL characters   : %s" % "".join(sorted(illegal)))
    if ragged:
        print("  !! NOT RECTANGULAR: %d record(s) differ in length, e.g. %s"
              % (len(ragged), ragged[:3]))

    if ungapped:
        s = sorted(ungapped)
        print("  ungapped length min/med/max : %d / %d / %d" % (
            s[0], s[len(s) // 2], s[-1]))

    if args.occupancy and occ:
        low = sum(1 for c in occ if c / n_seq < args.min_occupancy)
        empty = sum(1 for c in occ if c == 0)
        print("  columns below %.1f%% occupancy : %d" % (
            100 * args.min_occupancy, low))
        print("  all-gap columns              : %d" % empty)

    # ---- structural errors (never overridable) ----
    verdict = "PASS"
    reasons = []
    structural = []
    ratio_issues = []

    if ragged:
        structural.append("alignment is not rectangular (%d bad records)"
                          % len(ragged))
    if args.codon and width and width % 3 != 0:
        structural.append("codon alignment width %d is not a multiple of 3"
                          % width)
    if illegal:
        structural.append("illegal character(s) in the alignment: %s"
                          % "".join(sorted(illegal)))
    if mixed_alphabet:
        structural.append("alphabets are mixed (nucleotide-only %s + protein-only "
                          "%s), so no single substitution model exists"
                          % ("".join(sorted(mixed_nt)), "".join(sorted(mixed_aa))))

    # ---- ratio metrics (default FAIL, overridable with --accept-and-record) ----
    # Thresholds follow the marker family (M5): the historical 0.80 nucleotide
    # number is right for a conserved rRNA region and wrong for ITS or another
    # length-heterogeneous amplicon, where it would make the documented exception
    # the routine.  An explicit --max-gap* always wins, and the report says which
    # one is in force, so the number is never a hidden default.
    CLASS_GAP = {"rrna": 0.80, "organelle": 0.80, "functional-nuclear": 0.85,
                 "its": 0.90, "other": 0.90}
    base_gap = CLASS_GAP.get(args.marker_class, 0.80)
    aa_default = args.max_gap_aa if args.max_gap_aa is not None else min(0.95, base_gap + 0.05)
    codon_default = (args.max_gap_codon if args.max_gap_codon is not None
                     else min(0.95, base_gap + 0.10))
    nt_default = args.max_gap if args.max_gap is not None else base_gap
    if args.codon:
        gap_limit = codon_default
        mode = "codon"
        limit_src = ("--max-gap-codon" if args.max_gap_codon is not None
                     else "marker-class codon default")
    elif dtype == "protein":
        gap_limit = aa_default
        mode = "amino-acid"
        limit_src = ("--max-gap-aa" if args.max_gap_aa is not None
                     else "marker-class amino-acid default")
    else:
        gap_limit = nt_default
        mode = "nucleotide"
        limit_src = ("--max-gap" if args.max_gap is not None
                     else "marker-class nucleotide default")
    print("  marker class            : %s" % (args.marker_class or "(not given)"))
    print("  gap definition          : gap characters / (sequences x columns)")
    print("  gap limit for this mode : %.2f  (%s; from %s)"
          % (gap_limit, mode, limit_src))

    width_ratio = None
    overrides = []
    if args.fasta_raw:
        med, n_raw = median_len(args.fasta_raw, args.format)
        if med:
            width_ratio = (width or 0) / med
            print("  raw sequences           : %d (median length %g)" % (n_raw, med))
            print("  width / median raw len  : %.2f  (limit %g)"
                  % (width_ratio, args.max_col_ratio))
            if width_ratio > args.max_col_ratio:
                ratio_issues.append("width is %.2fx the median sequence length"
                                    % width_ratio)
                overrides.append({"gate": "G2", "metric": "width_ratio",
                                  "value": round(width_ratio, 6),
                                  "limit": args.max_col_ratio,
                                  "reason": args.override_reason})
        else:
            print("  !! --fasta-raw had no sequences; width-ratio test skipped")

    if gap_frac > gap_limit:
        ratio_issues.append("gap fraction %.3f exceeds %.3f" % (gap_frac, gap_limit))
        overrides.append({"gate": "G2", "metric": "gap_fraction",
                          "value": round(gap_frac, 6), "limit": gap_limit,
                          "reason": args.override_reason})

    if structural:
        verdict = "FAIL"
        reasons.extend(structural)
    if ratio_issues:
        if args.accept_and_record:
            print()
            print("OVERRIDE ACCEPTED (--accept-and-record): the following ratio "
                  "issues are recorded, NOT resolved:")
            for r in ratio_issues:
                print("  - %s" % r)
            print("  reason: %s" % args.override_reason)
            print("  These are written into the JSON report as gate_overrides, and "
                  "you MUST write them and their justification into the ops log.")
            print("  A gate that is waived silently is not a gate.")
            # An override records the ratio issues; it does not resurrect a
            # structural failure, and it must not be reported as a pass either.
            if not structural:
                reasons = ["OVERRIDDEN (recorded, not resolved): " + r
                           for r in ratio_issues]
        else:
            verdict = "FAIL"
            reasons.extend(ratio_issues)

    # ---- diagnostics matrix (DESIGN 5.4) ----
    safe_width = width or 0
    col_has_gap = [0] * safe_width
    col_has_amb = [0] * safe_width
    for _name, seq in sequences:
        for i, ch in enumerate(seq):
            if i >= safe_width:
                break
            if ch in GAP_CHARS:
                col_has_gap[i] += 1
            if ch in AMBIG_CHARS:
                col_has_amb[i] += 1
    gap_cols = sum(1 for c in col_has_gap if c > 0)
    all_gap_cols = sum(1 for c in col_has_gap if c == n_seq)
    amb_affected_cols = sum(1 for c in col_has_amb if c > 0)

    # terminal gap columns: columns that NO sequence fills at either end
    lead = 0
    while lead < safe_width and col_has_gap[lead] == n_seq:
        lead += 1
    trail = 0
    while trail < safe_width - lead and col_has_gap[safe_width - 1 - trail] == n_seq:
        trail += 1
    trimmed_end = lead + trail

    # pairwise gap fraction, sampled above a fixed number of pairs so that a
    # 190k-tip matrix does not become an O(n^2) job
    pair_vals = []
    npair_total = n_seq * (n_seq - 1) // 2
    MAX_PAIRS = 200000
    max_pairs = min(npair_total, MAX_PAIRS)
    if n_seq > 1 and max_pairs > 0:
        rng = random.Random(20261012)
        if npair_total <= MAX_PAIRS:
            pair_iter = ((i, j) for i in range(n_seq) for j in range(i + 1, n_seq))
            sampled = True
        else:
            def gen():
                for _ in range(MAX_PAIRS):
                    i = rng.randrange(n_seq)
                    j = rng.randrange(n_seq - 1)
                    if j >= i:
                        j += 1
                    yield (min(i, j), max(i, j))
            pair_iter = gen()
            sampled = False
        seen_pairs = set()
        for i, j in pair_iter:
            if (i, j) in seen_pairs:
                continue
            seen_pairs.add((i, j))
            si, sj = sequences[i][1], sequences[j][1]
            ident_len = min(len(si), len(sj))
            if ident_len == 0:
                continue
            either = sum(1 for k in range(ident_len)
                         if si[k] in GAP_CHARS or sj[k] in GAP_CHARS)
            pair_vals.append(either / ident_len)
    q = quantiles(pair_vals)

    # ATGC fraction: of the residues that are neither gap nor ambiguous
    atgc_n = 0
    atgc_d = 0
    for _name, seq in sequences:
        for ch in seq:
            if ch in GAP_CHARS or ch in AMBIG_CHARS:
                continue
            atgc_d += 1
            if ch in "ACGTUacgtu":
                atgc_n += 1
    atgc_fraction = (atgc_n / atgc_d) if atgc_d else None

    third_pos_rate = None
    if args.codon and safe_width:
        # fraction of column-pairs at codon position 3 that differ, sampled over
        # a bounded number of sequence pairs
        third = [i for i in range(2, safe_width, 3)]
        diff = 0
        tot = 0
        rng2 = random.Random(20261013)
        caps = min(max_pairs, 2000)
        for _ in range(caps):
            if n_seq < 2 or not third:
                break
            i = rng2.randrange(n_seq)
            j = rng2.randrange(n_seq - 1)
            if j >= i:
                j += 1
            si, sj = sequences[i][1], sequences[j][1]
            for k in third:
                if k < len(si) and k < len(sj):
                    if si[k] in GAP_CHARS or sj[k] in GAP_CHARS:
                        continue
                    tot += 1
                    if si[k].upper() != sj[k].upper():
                        diff += 1
        if tot:
            third_pos_rate = diff / tot

    # drop-col: SUGGEST columns, never edit the file
    drop_list = []
    if args.ambiguous_policy == "drop-col":
        thresh = args.ambiguous_max_frac
        for i in range(safe_width):
            frac = col_has_amb[i] / n_seq
            if frac > thresh:
                drop_list.append({"column_1based": i + 1,
                                  "ambiguous_fraction": round(frac, 6)})
        base = args.outdir or (os.path.dirname(os.path.abspath(args.json))
                               if args.json else os.getcwd())
        os.makedirs(base, exist_ok=True)
        dpath = os.path.join(base, "%s.drop_cols.txt" % os.path.basename(args.aln))
        with open(dpath, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("# columns whose ambiguous (N ? X) fraction exceeds %.3f\n"
                     % thresh)
            fh.write("# column<TAB>ambiguous_fraction   (SUGGESTION ONLY -- "
                     "this file is not modified)\n")
            for d in drop_list:
                fh.write("%d\t%.6f\n" % (d["column_1based"],
                                         d["ambiguous_fraction"]))
        print("  drop-col suggestion     : %d column(s) above %.1f%% ambiguous"
              % (len(drop_list), 100 * thresh))
        print("     list written to      : %s (the alignment file is NOT modified)"
              % dpath)

    diagnostics = {
        "gap_column_fraction": (gap_cols / safe_width) if safe_width else None,
        "all_gap_columns": all_gap_cols,
        "pairwise_gap_fraction": {
            "n_pairs_used": len(pair_vals),
            "n_pairs_total": npair_total,
            "sampled": npair_total > MAX_PAIRS,
            "min": q[0] if q else None,
            "q1": q[1] if q else None,
            "median": q[2] if q else None,
            "q3": q[3] if q else None,
            "max": q[4] if q else None,
        },
        "terminal_gap_columns": {
            "trimmed": trimmed_end,
            "lead": lead,
            "trail": trail,
            "columns_before": safe_width,
            "columns_after": safe_width - trimmed_end,
        },
        "atgc_fraction": atgc_fraction,
        "iupac_char_count": ambig,
        "iupac_sequences": sum(1 for _n, s
                               in sequences if any(c in AMBIG_CHARS for c in s)),
        "iupac_affected_columns": amb_affected_cols,
        "third_position_difference_rate": third_pos_rate,
        "ambiguous_policy": args.ambiguous_policy,
        "ambiguous_max_frac": args.ambiguous_max_frac,
        "drop_col_suggestions": drop_list,
    }

    print()
    if verdict == "PASS":
        if reasons:
            print("RESULT: PASS (with recorded override)")
        else:
            print("RESULT: PASS -- alignment looks usable.")
    else:
        print("RESULT: FAIL")
        for r in reasons:
            print("  - %s" % r)
        if not structural and ratio_issues:
            # Give the escape hatch as a copy-pasteable command rather than as a
            # sentence: this is a RATIO failure, and the documented handling is
            # "override explicitly and record why".  (A structural failure gets no
            # such text -- it cannot be overridden.)
            print()
            print("  This is a RATIO failure, so it CAN be overridden with a recorded")
            print("  reason.  Copy, fix the reason, and keep the JSON it writes:")
            print("    python scripts/aln_qc.py --aln %s \\" % args.aln)
            print("        --accept-and-record \\")
            print("        --override-reason \"<why the gappy matrix is still usable "
                  "for this question>\" \\")
            print("        --json %s" % (args.json or "qc/aln.json"))
            print("  Do NOT lower the threshold to make it pass: a lower limit is a")
            print("  claim about this study, and it has to be argued (D5).")
        print()
        print(DIAGNOSIS)

    if args.json:
        d = os.path.dirname(os.path.abspath(args.json))
        if d:
            os.makedirs(d, exist_ok=True)
        rec = {"label": label, "file": args.aln, "sequences": n_seq,
               "columns": width, "residues": residues,
               "gap_fraction": gap_frac,
               "gap_fraction_denominator": "n_seq * width",
               "gap_fraction_if_residues": (gap_frac_by_residues
                                            if abs(gap_frac - gap_frac_by_residues)
                                            > 1e-9 else None),
               "gap_char_count": gap_char_count,
               "ambiguous": ambig, "alphabet": alpha, "non_gap": nongap,
               "inferred_data_type": dtype, "ragged": ragged[:20],
               "structural_errors": structural,
               "marker_class": args.marker_class,
               "gap_limit": gap_limit,
               "gap_limit_mode": mode,
               "gap_limit_source": limit_src,
               "ambiguous_policy": args.ambiguous_policy,
               "ambiguous_max_frac": args.ambiguous_max_frac,
               "accepted_and_recorded": bool(args.accept_and_record),
               "override_reason": (args.override_reason
                                   if args.accept_and_record else None),
               "gate_overrides": (overrides if args.accept_and_record else []),
               "diagnostics": diagnostics,
               "ungapped_min": min(ungapped), "ungapped_median":
                   sorted(ungapped)[len(ungapped) // 2], "ungapped_max":
                   max(ungapped), "verdict": verdict, "reasons": reasons}
        with open(args.json, "w", encoding="utf-8", newline="\n") as fh:
            json.dump(rec, fh, indent=2, ensure_ascii=False)
        print("metrics written: %s" % args.json)

    return 1 if verdict == "FAIL" else 0


if __name__ == "__main__":
    sys.exit(main())
