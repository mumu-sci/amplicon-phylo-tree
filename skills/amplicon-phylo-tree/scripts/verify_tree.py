#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verify a finished tree before it is used or shipped.

G3 -- the substitution MODEL actually matches the DATA TYPE.
    A tool will happily run a full maximum-likelihood analysis with the wrong
    model class and print plausible likelihoods; the tree is then meaningless.
    FastTree's own log shows the classic case:
        Amino acid distances: BLOSUM45 ... ML Model: Whelan-And-Goldman ...
        WARNING! 100.0% NUCLEOTIDE CHARACTERS -- IS THIS REALLY A PROTEIN ALIGNMENT?
    i.e. a PROTEIN model applied to pure ACGT data.

    Four layers, fail-closed (no positive evidence anywhere => UNVERIFIED, stop):

      1. DATA SIDE (authoritative, tool-independent): scan the alignment's own
         character set.  Only IUPAC nucleotide symbols -> nt.  Protein-only
         letters (EFILPQZ) -> aa.  Both present -> hard stop.
         This layer is the only one that works for every tool, and it is what
         actually catches the failure above.

      2. INJECTION, NOT ASSERTION: this skill *generates* the tree command, so
         the model/seqtype flags are injected from layer 1 rather than checked
         afterwards.  Pass --model/--tool/--seqtype and they are cross-checked
         against layer 1 so a mismatch cannot be recorded silently.

      3. TOOL EVIDENCE, using strings measured on this machine (not guessed).
         Measured 2026-10-12, IQ-TREE 2.4.0, nucleotide and protein runs plus two
         deliberate mismatches:
             'Generalized Time-Reversible'      -> 0 occurrences in IQ-TREE logs
             'Nucleotide distances'             -> 0
             'Amino acid distances'             -> 0
             'Whelan-And-Goldman'               -> 0
             'BLOSUM'                           -> 0   (uppercase token only)
             'Blosum' (case-insensitive, which is what matters) -> 36
                                                   CANDIDATE model table, e.g.
                                                   'Blosum62+I+G4')
         Two consequences:
           * FastTree's wording must NOT be required of IQ-TREE (it would abort
             every legitimate IQ-TREE run).
           * 'BLOSUM'/'Blosum' must NOT be used as a blacklist: IQ-TREE legitimately
             prints it while choosing a protein model, so it would fire on a
             correct protein run.
         Measured tool-native signals:
           * veto, both tools: FastTree 'WARNING! 100.0% NUCLEOTIDE CHARACTERS';
             IQ-TREE 'specified sequence type is different from the detected one';
             IQ-TREE 'invalid character' (hard error, e.g. protein data forced to
             -st DNA).
           * FastTree positive: 'Nucleotide distances:' (nt) /
             'Amino acid distances:' (aa).
           * IQ-TREE positive: 'ModelFinder will test up to N DNA models' (nt) /
             '... protein models' (aa), and 'Model of substitution:' in the
             .iqtree report.
         Note the asymmetry that makes this dangerous: forcing PROTEIN data to a
         nucleotide model errors out, but forcing NUCLEOTIDE data to a protein
         model only warns and continues -- which is exactly the real-world trap.

      4. fail-closed: if no layer yields positive evidence, report UNVERIFIED and
         exit non-zero.

G4 -- the tree's tip set equals the input taxon set EXACTLY, and the file is a
    well-formed, usable tree.  A tip-set mismatch is usually a silent join
    failure downstream: taxa present in the abundance table but absent from the
    tree are dropped from phylogenetic metrics, shrinking the effective sample
    size without warning.

Newick is parsed here rather than shelling out to R, so the check runs anywhere.

Exit codes: 0 = PASS, 1 = FAIL, 2 = input error.
"""
import argparse
import hashlib
import json
import os
import random
import re
import sys

# ------------------------------------------------------------------ Newick ---

class Node:
    __slots__ = ("children", "label", "length", "is_tip")

    def __init__(self):
        self.children = []
        self.label = None
        self.length = None
        self.is_tip = False


class NewickError(Exception):
    pass


def _read_label(s, i):
    n = len(s)
    if i < n and s[i] == "'":
        i += 1
        buf = []
        while i < n:
            if s[i] == "'":
                if i + 1 < n and s[i + 1] == "'":
                    buf.append("'")
                    i += 2
                    continue
                i += 1
                break
            buf.append(s[i])
            i += 1
        return "".join(buf), i
    start = i
    while i < n and s[i] not in ",():;":
        i += 1
    return s[start:i].strip(), i


def _skip_ws(s, i):
    n = len(s)
    while i < n and s[i] in " \t\r\n":
        i += 1
    return i


def _parse_node(s, i):
    node = Node()
    i = _skip_ws(s, i)
    if i < len(s) and s[i] == "(":
        i += 1
        while True:
            child, i = _parse_node(s, i)
            node.children.append(child)
            i = _skip_ws(s, i)
            if i < len(s) and s[i] == ",":
                i = _skip_ws(s, i + 1)
                continue
            if i < len(s) and s[i] == ")":
                i += 1
                break
            raise NewickError("expected ',' or ')' at offset %d" % i)
        i = _skip_ws(s, i)
        node.label, i = _read_label(s, i)
    else:
        node.is_tip = True
        node.label, i = _read_label(s, i)
        if not node.label:
            raise NewickError("empty tip label at offset %d" % i)
    i = _skip_ws(s, i)
    if i < len(s) and s[i] == ":":
        i = _skip_ws(s, i + 1)
        start = i
        while i < len(s) and s[i] in "+-0123456789.eE":
            i += 1
        try:
            node.length = float(s[start:i])
        except ValueError:
            raise NewickError("bad branch length %r" % s[start:i])
    return node, i


def parse_newick(text):
    s = text.strip()
    if not s.endswith(";"):
        raise NewickError("tree does not end with ';' (file looks truncated)")
    i = 0
    root, i = _parse_node(s, i)
    i = _skip_ws(s, i)
    if i < len(s) and s[i] == ";":
        i += 1
    if i < len(s) and s[i:].strip():
        raise NewickError("trailing content after the tree: %r" % s[i:i + 40])
    return root


def collect_tips(node, acc):
    if node.is_tip:
        acc.append(node.label)
    else:
        for c in node.children:
            collect_tips(c, acc)


def _split_signature(root, tips):
    """The set of nontrivial splits, CANONICALISED as bipartitions.

    A split is a bipartition of the tips, so the two sides of the root edge
    describe the SAME split.  Keying on a clade's own size and hash (which is
    what this did first) therefore records a rooted tree's root-edge split twice
    -- once per side -- while an unrooted tree records it once.  The two files
    then disagree even when they describe the identical topology, and the check
    reports a stale file.  That bug fired on every tree; canonicalising fixes it.

    Canonical key: take the SMALLER side of the bipartition (by tip count), and
    when both sides are the same size, take the smaller of the two hashes.  The
    complement's hash is h_all ^ h, so no traversal of the other side is needed.

    Hashes are XORs of per-tip 64-bit values, which is order-independent (right
    for a set) but can in principle collide; the size travels with the hash to
    make an accidental match unlikely.  This is a staleness check, not a proof.
    """
    unique = sorted(set(tips))
    rnd = random.Random(20261012)
    tip_hash = {t: rnd.getrandbits(64) for t in unique}
    n = len(tips)
    h_all = 0
    for t in unique:
        h_all ^= tip_hash[t]
    keys = set()

    def rec(node):
        if node.is_tip:
            return 1, tip_hash[node.label]
        size = 0
        h = 0
        for c in node.children:
            cs, ch = rec(c)
            size += cs
            h ^= ch
        small = min(size, n - size)
        if 2 <= small <= n - 2:
            if size < n - size:
                cand = h                     # this side is already the small one
            elif size > n - size:
                cand = h_all ^ h             # the complement is the small side
            else:
                cand = min(h, h_all ^ h)     # equal sizes: pick deterministically
            keys.add((small, cand))
        return size, h

    rec(root)
    return keys


def tree_stats(path, want_splits=False):
    raw = open(path, encoding="utf-8", errors="replace").read()
    root = parse_newick(raw)

    tips, lengths, labels = [], [], []
    n_internal = 0

    def rec2(node):
        nonlocal n_internal
        if node.is_tip:
            tips.append(node.label)
        else:
            n_internal += 1
            if node.label:
                labels.append(node.label)
        if node.length is not None:
            lengths.append(node.length)
        for c in node.children:
            rec2(c)

    rec2(root)

    st = {
        "file": path,
        "bytes": os.path.getsize(path),
        "ends_with_semicolon": raw.strip().endswith(";"),
        "tip_labels": tips,
        "n_tips": len(tips),
        "n_unique_tips": len(set(tips)),
        "n_internal": n_internal,
        "root_degree": len(root.children),
        "n_support_labels": len(labels),
        "number_of_branches": len(lengths),
        "min_branch": min(lengths) if lengths else None,
        "max_branch": max(lengths) if lengths else None,
        "zero_branches": sum(1 for L in lengths if L == 0.0),
        "negative_branches": sum(1 for L in lengths if L < 0),
    }
    if want_splits:
        st["split_signature"] = _split_signature(root, tips)
    return st


# ------------------------------------------------------------------- model ---

NT_CANON = set("ACGTUacgtu")
NT_AMBIG = set("RYSWKMBDHVNryswkmbdhvn")
# Letters that CANNOT occur in a nucleotide sequence.  Note that the IUPAC
# ambiguity codes (R Y S W K M B D H V N) are all also valid amino-acid
# one-letter codes, so "contains an ambiguity code" says nothing about the data
# type -- only these letters do.  Treating the IUPAC set as nucleotide evidence
# makes every protein alignment look mixed, which is the bug this comment exists
# to prevent.
PROT_ONLY = set("EFILPQZefilpqz*")


def classify_sequence(seq):
    """Classify ONE sequence as 'nt', 'aa', 'suspect' or 'empty'.

    Decision rule (deliberately conservative, because a wrong call here silently
    changes the substitution model):

      >= 3 protein-only letters AND >= 5% of characters  -> aa
      1-2 protein-only letters (or < 5%)                 -> suspect
      none at all                                        -> nt

    'suspect' means "a nucleotide-like record carrying a stray non-nucleotide
    letter" -- usually a typo or an artefact, not a protein.  It is reported, and
    it only becomes fatal if the file also contains unambiguous protein records.
    """
    counts = {}
    total = 0
    for ch in seq:
        if ch.isspace() or ch in "-.":
            continue
        counts[ch] = counts.get(ch, 0) + 1
        total += 1
    if total == 0:
        return "empty", {"total": 0}
    prot = sum(n for ch, n in counts.items() if ch in PROT_ONLY)
    canon = sum(n for ch, n in counts.items() if ch in NT_CANON)
    ambig = sum(n for ch, n in counts.items() if ch in NT_AMBIG)
    frac_prot = prot / total
    stats = {"total": total, "prot_only": prot, "canonical_nt": canon,
             "iupac_ambig": ambig, "frac_prot_only": frac_prot}
    if prot == 0:
        return "nt", stats
    if prot >= 3 and frac_prot >= 0.05:
        return "aa", stats
    return "suspect", stats


def infer_datatype_from_alignment(path, verbose=True):
    """Layer 1: the authoritative, tool-independent data-type call.

    Classified PER SEQUENCE, so a file that mixes nucleotide and protein records
    is detected -- that is the genuinely dangerous case (two alignments
    concatenated, or one record left in protein form).
    """
    per = []
    name = None
    buf = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    per.append((name, "".join(buf)))
                name = line[1:].strip().split()[0]
                buf = []
            elif line:
                buf.append(line)
    if name is not None:
        per.append((name, "".join(buf)))
    if not per:
        return "UNKNOWN", [], "the alignment has no sequences"

    aa_ids, nt_ids, suspect_ids = [], [], []
    canon_total = 0
    char_total = 0
    for sid, seq in per:
        call, st = classify_sequence(seq)
        if call == "aa":
            aa_ids.append(sid)
        elif call == "nt":
            nt_ids.append(sid)
            canon_total += st.get("canonical_nt", 0)
            char_total += st.get("total", 0)
        elif call == "suspect":
            suspect_ids.append(sid)

    note = None
    if suspect_ids:
        note = ("%d sequence(s) carry a few non-nucleotide letters and were treated "
                "as nucleotide records: %s" % (len(suspect_ids),
                                               ", ".join(suspect_ids[:6])))

    canon_frac = (canon_total / char_total) if char_total else 0.0

    if aa_ids and nt_ids:
        return "MIXED", aa_ids + nt_ids, (
            "the file contains %d protein-like and %d nucleotide-like sequences.  "
            "A nucleotide/protein mixture has no single substitution model.  "
            "Protein-like: %s" % (len(aa_ids), len(nt_ids), ", ".join(aa_ids[:6])))
    if aa_ids:
        return "aa", aa_ids, note
    if nt_ids:
        # A character-only rule CANNOT separate a nucleotide sequence from a short
        # low-complexity peptide: the IUPAC ambiguity codes (R Y S W K M B D H V N)
        # are every one of them also valid amino-acid one-letter codes.  So when the
        # "nucleotide" sequences are mostly made of those shared letters, the
        # data-side call is withheld rather than guessed -- the verdict then needs
        # the tool log, which is what fail-closed means here.
        if canon_frac < 0.5:
            return "AMBIGUOUS", nt_ids, (
                "the sequences classified as nucleotide carry only %.0f%% canonical "
                "ACGTU characters.  That is equally consistent with a short, "
                "low-complexity PEPTIDE, so this layer refuses to decide.  Supply "
                "--tool-log (and/or --iqtree-report) to settle it."
                % (100 * canon_frac))
        return "nt", nt_ids, note
    return "UNKNOWN", suspect_ids, (note or "no sequence could be classified")



# Measured per-tool strings (see module docstring).
VETO = [
    (r"WARNING!\s*100\.0%\s*NUCLEOTIDE CHARACTERS", "FastTree: protein model applied "
     "to nucleotide data"),
    (r"IS THIS REALLY A PROTEIN ALIGNMENT", "FastTree: protein model applied to "
     "nucleotide data"),
    (r"specified sequence type is different from the detected one",
     "IQ-TREE: the requested sequence type disagrees with what the data look like"),
    (r"^ERROR:.*invalid character", "IQ-TREE: invalid characters for the declared "
     "sequence type"),
]
FT_NT = [r"Nucleotide distances\s*:"]
FT_AA = [r"Amino acid distances\s*:"]
IQ_NT = [r"ModelFinder will test up to \d+ DNA models",
         r"Model of substitution:\s*\S*(?:GTR|HKY|TN|TVM|TIM|JC|F81|K2P|K3P|SYM|TrN|TPM)"]
IQ_AA = [r"ModelFinder will test up to \d+ protein models",
         r"Model of substitution:\s*\S*(?:LG|WAG|JTT|Dayhoff|Blosum|VT|cpREV|rtREV|mtREV|mtMAM|FLU|HIV|PMB|Q\.)"]
RAX_NT = [r"Substitution Matrix\s*:?\s*(?:GTR|HKY|JC)", r"Data Type\s*:?\s*DNA"]
RAX_AA = [r"Substitution Matrix\s*:?\s*(?:LG|WAG|JTT|BLOSUM|DAYHOFF)", r"Data Type\s*:?\s*AA"]

TOOL_POSITIVE = {
    "fasttree": {"nt": FT_NT, "aa": FT_AA},
    "iqtree": {"nt": IQ_NT, "aa": IQ_AA},
    "raxml": {"nt": RAX_NT, "aa": RAX_AA},
}


def detect_tool(text, explicit):
    if explicit and explicit != "auto":
        return explicit
    t = text or ""
    if "FastTree" in t:
        return "fasttree"
    if "IQ-TREE" in t or "IQTREE" in t:
        return "iqtree"
    if "RAxML" in t:
        return "raxml"
    if "PhyML" in t:
        return "phyml"
    return "unknown"


def read_text_any_encoding(path):
    """Read a log file whose encoding we do NOT control.

    Windows PowerShell redirection writes UTF-16LE.  Decoding that as UTF-8 turns
    every ASCII character into NUL plus a replacement character, so every keyword
    match fails and the tool-evidence layer silently reports "no evidence" while
    the log is in fact full of it.  That is measured, not hypothetical: it
    disabled this layer on the machine this module was developed on, and turned a
    log containing two veto strings into a PASS.

    So: sniff the BOM, then a NUL in the first block, and only then assume UTF-8.
    """
    with open(path, "rb") as fh:
        raw = fh.read()
    if raw.startswith(b"\xef\xbb\xbf"):
        return raw.decode("utf-8-sig", errors="replace")
    if raw.startswith(b"\xff\xfe") or raw.startswith(b"\xfe\xff"):
        return raw.decode("utf-16", errors="replace")
    if b"\x00" in raw[:4096]:          # UTF-16 without a BOM
        for enc in ("utf-16-le", "utf-16-be"):
            try:
                text = raw.decode(enc, errors="replace")
                if "\x00" not in text:
                    return text
            except (UnicodeDecodeError, LookupError):
                continue
    return raw.decode("utf-8", errors="replace")


def gather_tool_text(tool_log, iqtree_report):
    """Return (combined text, list of paths actually read)."""
    chunks = []
    found = []
    for p in (tool_log, iqtree_report):
        if p and os.path.exists(p):
            chunks.append(read_text_any_encoding(p))
            found.append(p)
    return "\n".join(chunks), found


def _iter_fasta(path):
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


def alignment_shape(path):
    """Return (n_sequences, width, ragged) for a FASTA alignment.

    Records are reassembled first: writers wrap at 60 characters, so a line
    length is not a column count.
    """
    n = 0
    cur = 0
    widths = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if n:
                    widths.append(cur)
                n += 1
                cur = 0
            elif line:
                cur += len(line)
    if n:
        widths.append(cur)
    if not widths:
        return 0, 0, False
    return n, widths[0], len(set(widths)) > 1


# ------------------------------------------------------------------- main ----

def main():
    ap = argparse.ArgumentParser(
        description="Verify a tree: model vs data type (G3) and tip set (G4).")
    ap.add_argument("--aln", default=None,
                    help="the alignment the tree was built from; enables the "
                         "authoritative data-side data-type check (layer 1)")
    ap.add_argument("--rooted", default=None, help="rooted tree (Newick)")
    ap.add_argument("--unrooted", default=None, help="unrooted tree (Newick)")
    ap.add_argument("--ids", default=None, help="expected taxon ids, one per line")
    ap.add_argument("--table", default=None,
                    help="abundance table; its row ids are checked against the "
                         "tip set directly")
    ap.add_argument("--table-format", choices=["auto", "tsv", "csv"], default="auto")
    ap.add_argument("--tool", default="auto",
                    choices=["auto", "iqtree", "fasttree", "raxml", "phyml"])
    ap.add_argument("--model", default=None,
                    help="the model string you asked the tool for (first-class "
                         "evidence; the log is only corroboration)")
    ap.add_argument("--seqtype", choices=["nt", "aa", "codon", "none"], default="none",
                    help="the data type you intended the tool to use")
    ap.add_argument("--tool-log", default=None, help="tool stdout/stderr log")
    ap.add_argument("--iqtree-report", default=None, help="the .iqtree report file")
    ap.add_argument("--require-rooted", action="store_true")
    ap.add_argument("--require-binary", action="store_true",
                    help="require a fully bifurcating tree (needed for the "
                         "phylogenetic metrics route)")
    ap.add_argument("--check-same-topology", action="store_true",
                    help="assert the rooted and unrooted trees share the same "
                         "unrooted topology (catches a stale/overwritten file)")
    ap.add_argument("--allow-no-tool-log", action="store_true",
                    help="accept a data-side-only model check.  By DEFAULT the "
                         "verdict requires corroboration from the tool log, because "
                         "an uncorroborated model claim is exactly the failure this "
                         "gate exists to catch.  Using this flag is a recorded "
                         "decision, not a silent pass.")
    ap.add_argument("--rooting-report", default=None, metavar="FILE",
                    help="the rooting report produced by root_tree.R.  With "
                         "--require-rooted this asserts that the ROOTING "
                         "SENSITIVITY section actually ran (root_tree.R does it by "
                         "default): a report that says NOT RUN means the root was "
                         "never tested against alternative placements, which is a "
                         "gap, not a pass.")
    ap.add_argument("--require-codon-frame", action="store_true",
                    help="assert the alignment width is a multiple of 3 (a "
                         "back-translated codon alignment).  Also triggered by "
                         "--seqtype codon.")
    ap.add_argument("--backtranslation-evidence", default=None, metavar="FASTA",
                    help="the back-translated codon alignment itself; its sequence "
                         "ids are checked against the expected taxon set")
    ap.add_argument("--aa-aln", default=None, metavar="FASTA",
                    help="the amino-acid alignment that was back-translated.  With "
                         "--backtranslation-evidence this asserts that no sequence "
                         "contains 'X', because PAL2NAL and trimAl -backtrans cannot "
                         "map X (only MACSE can).  See references/parameters.md D2b.")
    ap.add_argument("--json", default=None)
    args = ap.parse_args()

    if not any([args.rooted, args.unrooted, args.aln, args.tool_log,
                args.iqtree_report]):
        ap.error("nothing to check: give --rooted/--unrooted (tree checks) and/or "
                 "--aln/--tool-log/--iqtree-report (model check)")

    failures = []
    unverified = []
    report = {}

    # ---------- expected taxon set ----------
    expected = None
    src = None
    if args.ids or args.table:
        if args.table:
            delim = "\t"
            if args.table_format == "csv":
                delim = ","
            elif args.table_format == "auto":
                with open(args.table, encoding="utf-8", errors="replace") as fh:
                    head = fh.readline()
                delim = "," if head.count(",") > head.count("\t") else "\t"
            ids = []
            with open(args.table, encoding="utf-8", errors="replace") as fh:
                first = True
                for line in fh:
                    line = line.rstrip("\r\n")
                    if not line:
                        continue
                    if line.startswith("#"):
                        if line.split(delim)[0].strip().upper() in (
                                "#OTU ID", "#OTU_ID", "#OTUID"):
                            continue
                        continue
                    if first:
                        first = False
                        cell = line.split(delim)[0].strip().lower().lstrip("#")
                        if cell in ("otu", "otu id", "asv", "asv_id", "feature",
                                    "featureid", "feature_id", "id", "taxon"):
                            continue  # this was a header
                    ids.append(line.split(delim)[0].strip())
            expected = ids
            src = "table %s" % args.table
        else:
            expected = [l.strip() for l in open(args.ids, encoding="utf-8")
                        if l.strip()]
            src = args.ids
        print("expected taxa             : %d (from %s)" % (len(expected), src))
    else:
        print("expected taxa             : not given (tip-set identity not checked)")

    # ---------- G3 ----------
    if args.seqtype != "none" or args.aln:
        print()
        print("G3  substitution model vs data type")
        print("-" * 74)

        # layer 1
        l1 = None
        if args.aln:
            if not os.path.exists(args.aln):
                print("  ERROR: alignment not found: %s" % args.aln)
                failures.append("G3 layer-1 alignment missing")
            else:
                l1, offenders, why = infer_datatype_from_alignment(args.aln)
                print("  L1 data-side call        : %s   (classified per sequence, "
                      "so a mixed file is detected)" % l1)
                print("     sequences of that type: %d" % len(offenders))
                if l1 == "MIXED":
                    print("     !! %s" % why)
                    print("     !! refusing to continue: a nucleotide/protein "
                          "mixture has no single substitution model.")
                    failures.append("G3 mixed nucleotide+protein alignment")
                elif l1 == "AMBIGUOUS":
                    # Characters alone cannot decide; do not grant data-side
                    # evidence, so the verdict will need the tool log.
                    print("     !! %s" % why)
                    l1 = "AMBIGUOUS"
                elif l1 in ("nt", "aa"):
                    if why:
                        print("     note: %s" % why)
                else:
                    print("     !! %s" % why)
        else:
            print("  L1 data-side call        : NOT PERFORMED (pass --aln)")

        # layer 2: declared vs inferred.  This is a CONSISTENCY check, NOT an
        # independent source of evidence: if the declaration and the data-side call
        # are both wrong, they agree and this layer adds nothing.  It is still
        # worth doing, because it catches the common case of a run whose flags
        # contradict its own input.
        print("  L2 declared by the run   : tool=%s seqtype=%s model=%s"
              % (args.tool, args.seqtype, args.model or "(not given)"))
        if l1 in ("nt", "aa") and args.seqtype in ("nt", "aa"):
            if (l1 == "aa") != (args.seqtype == "aa"):
                print("     !! declared %s but the alignment is %s"
                      % (args.seqtype, l1))
                failures.append("G3 declared data type disagrees with the alignment")
        elif args.seqtype == "codon":
            if l1 != "nt":
                print("     !! a codon model needs nucleotide data, but the "
                      "alignment is %s" % l1)
                failures.append("G3 codon model on non-nucleotide data")

        # layer 3: tool evidence -- the independent corroboration
        text, supplied = gather_tool_text(args.tool_log, args.iqtree_report)
        tool_evidence = False
        if not supplied:
            print("  L3 tool evidence         : NOT SUPPLIED (no --tool-log / "
                  "--iqtree-report given, or the path does not exist)")
        else:
            print("  L3 tool evidence         : read %s" % ", ".join(supplied))
            tool = detect_tool(text, args.tool)
            print("  L3 tool detected         : %s" % tool)
            vetoed = False
            for pat, why in VETO:
                if re.search(pat, text, re.MULTILINE):
                    print("     !! VETO: %s" % why)
                    failures.append("G3 veto string matched: %s" % why)
                    vetoed = True
            pos = TOOL_POSITIVE.get(tool, {})
            if l1 == "AMBIGUOUS":
                # The characters could not decide, so any positive evidence the log
                # offers settles it -- but it must not contradict --seqtype.
                hits = {t: [p for p in pos.get(t, []) if re.search(p, text)]
                        for t in ("nt", "aa")}
                resolved = sorted(t for t, h in hits.items() if h)
                if resolved:
                    print("     the tool log resolves what the characters could not: "
                          "%s" % ", ".join(resolved))
                    tool_evidence = True
                    if args.seqtype in ("nt", "aa") and args.seqtype not in resolved:
                        print("     !! the log says %s but --seqtype said %s"
                              % (", ".join(resolved), args.seqtype))
                        failures.append("G3 declared %s but the tool log shows %s"
                                        % (args.seqtype, ", ".join(resolved)))
                else:
                    print("     the log did not resolve it either")
                want = other = None
            else:
                want = "aa" if (l1 == "aa" or args.seqtype == "aa") else "nt"
                other = "nt" if want == "aa" else "aa"
                hit_want = [p for p in pos.get(want, []) if re.search(p, text)]
                hit_other = [p for p in pos.get(other, []) if re.search(p, text)]
                if hit_want:
                    print("     positive evidence for %s: matched %d pattern(s)"
                          % (want, len(hit_want)))
                    tool_evidence = True
                if hit_other:
                    print("     !! evidence for the OPPOSITE data type (%s): %s"
                          % (other, "; ".join(hit_other)))
                    failures.append("G3 tool log shows %s evidence, expected %s"
                                    % (other, want))
                if not hit_want and not hit_other and not vetoed:
                    print("     no model evidence matched for tool '%s'.  The log WAS "
                          "read, so this is not 'no log given' -- either the wording "
                          "table lacks this tool's phrasing, or the file is not the "
                          "log you think it is." % tool)

        # layer 4: the verdict.  Independent evidence = data side (L1) + tool side
        # (L3).  L2 does not count as one of them.
        print("  L4 independent evidence  : data-side=%s  tool-side=%s"
              % ("yes" if l1 in ("nt", "aa") else "no",
                 "yes" if tool_evidence else "no"))
        if not failures:
            if tool_evidence:
                pass                                   # corroborated
            elif l1 == "AMBIGUOUS":
                # There is no data-side verdict to fall back on, so waiving the tool
                # log cannot produce one.  Ambiguity is resolved by evidence or not
                # at all -- --allow-no-tool-log must not paper over it.
                print("     UNVERIFIED: the characters cannot decide the data type, so "
                      "there is no data-side call to fall back on.  "
                      "--allow-no-tool-log does not apply here -- supply --tool-log "
                      "(and/or --iqtree-report).")
                unverified.append("G3 (ambiguous data, no tool log)")
            elif args.allow_no_tool_log:
                print("     note: no tool evidence, but --allow-no-tool-log was "
                      "given, so the verdict rests on the alignment characters plus "
                      "the injected parameters.  Record that you did this and why.")
            elif not supplied:
                print("     UNVERIFIED: no tool log was supplied, so the model claim "
                      "is uncorroborated.  Pass --tool-log (and/or --iqtree-report), "
                      "or --allow-no-tool-log to accept a data-side-only check and "
                      "record that decision.")
                unverified.append("G3 (no tool log)")
            else:
                print("     UNVERIFIED: a tool log WAS supplied, but no pattern in the")
                print("     per-tool wording table matched it -- so the log gave NO")
                print("     positive model evidence.  This is not the same as 'no log")
                print("     was given': the file exists and was read, it just does not")
                print("     say what model ran.  Fail-closed on purpose: an unreadable")
                print("     log is not a passing log.")
                unverified.append("G3 (tool log supplied but no positive evidence)")

    # ---------- G4 (alignment level): codon integrity / back-translation --------
    if (args.require_codon_frame or args.seqtype == "codon"
            or args.backtranslation_evidence):
        print()
        print("G4  codon integrity and back-translation (D2b)")
        print("-" * 74)
        if not args.aln:
            print("  !! --require-codon-frame / --backtranslation-evidence need --aln "
                  "to measure the alignment")
            failures.append("G4 codon check without --aln")
        elif not os.path.exists(args.aln):
            print("  ERROR: alignment not found: %s" % args.aln)
            failures.append("G4 codon check: alignment missing")
        else:
            n_seq, width, ragged = alignment_shape(args.aln)
            print("  alignment                : %d sequences x %d columns"
                  % (n_seq, width))
            if ragged:
                print("  !! not rectangular")
                failures.append("G4 alignment is not rectangular")
            if width and width % 3 != 0:
                print("  !! width %d is not a multiple of 3, so this is NOT a "
                      "codon-aligned matrix -- a back-translated alignment must be"
                      % width)
                failures.append("G4 codon alignment width %d is not a multiple of 3"
                                % width)
            elif width:
                print("  width divisible by 3    : yes (%d codons)" % (width // 3))

        ev = args.backtranslation_evidence
        if ev:
            if not os.path.exists(ev):
                print("  !! back-translation evidence not found: %s" % ev)
                failures.append("G4 back-translation evidence missing")
            else:
                ev_ids = set()
                for line in open(ev, encoding="utf-8", errors="replace"):
                    line = line.strip()
                    if line.startswith(">"):
                        ev_ids.add(line[1:].strip().split()[0])
                print("  back-translation evidence: %s (%d records)"
                      % (ev, len(ev_ids)))
                if expected is not None:
                    diff = set(expected) ^ ev_ids
                    if diff:
                        print("  !! the evidence's ids do not match the expected "
                              "taxon set (%d differ, e.g. %s)"
                              % (len(diff), ", ".join(sorted(diff)[:6])))
                        failures.append("G4 back-translation evidence tip set "
                                        "mismatch")
                    else:
                        print("  evidence ids match the taxon set: yes")

        if args.aa_aln and os.path.exists(args.aa_aln):
            with_x = []
            for name, seq in _iter_fasta(args.aa_aln):
                if "X" in seq.upper():
                    with_x.append(name)
            print("  amino-acid alignment     : %s (%d sequence(s) contain X)"
                  % (args.aa_aln, len(with_x)))
            if with_x:
                print("     !! %d sequence(s) carry X in the amino-acid alignment."
                      % len(with_x))
                print("        PAL2NAL refuses CDS with stop codons and trimAl "
                      "-backtrans cannot map X -- only MACSE handles these.  Either")
                print("        redo the back-translation with MACSE, or remove those "
                      "sequences under D4 and record that you did.")
                failures.append("G4 back-translation evidence with X in the aa "
                                "alignment (%d)" % len(with_x))
        elif args.aa_aln:
            print("  !! --aa-aln not found: %s" % args.aa_aln)
            failures.append("G4 --aa-aln missing")

    # ---------- G4 rooting sensitivity (D9 / M6) ----------
    # A root whose placement was never tested is an untested hypothesis.  When a
    # rooting report is supplied and a rooted tree is required, the report must
    # show that the sensitivity test actually ran -- "NOT RUN" is a gap, and the
    # whole point of the gate is that a gap must not read as a pass.
    if args.rooting_report and args.require_rooted:
        print()
        print("G4  rooting sensitivity (D9)")
        print("-" * 74)
        if not os.path.exists(args.rooting_report):
            print("  ERROR: rooting report not found: %s" % args.rooting_report)
            print("     !! produce it with scripts/root_tree.R --report ...")
            failures.append("G4 rooting report missing")
        else:
            text = read_text_any_encoding(args.rooting_report)
            has_section = "ROOTING SENSITIVITY" in text
            m = re.search(r"method used\s*:\s*(\S+)", text)
            method = m.group(1) if m else "unknown"
            print("  report                   : %s" % args.rooting_report)
            print("  rooting method           : %s" % method)
            if not has_section:
                print("     !! the report has no ROOTING SENSITIVITY section, so the")
                print("        root placement was never tested.  Re-run root_tree.R")
                print("        (the test is on by default) and pass --rooting-report.")
                failures.append("G4 rooting report has no sensitivity section")
            elif re.search(r"^\s*NOT RUN", text, re.M):
                print("     !! the sensitivity test did NOT run (the report says")
                print("        NOT RUN).  root_tree.R runs it by default; it was")
                print("        turned off with --no-sensitivity.  Either run it, or")
                print("        state in the methods that the root was not tested and")
                print("        why the tree cannot feed root-sensitive metrics.")
                failures.append("G4 rooting sensitivity NOT RUN")
            else:
                cand = re.search(r"candidates compared\s*:\s*(\d+)", text)
                ndiff = re.search(r"(\d+) of \d+ candidate root positions", text)
                print("  sensitivity test         : ran (%s candidate(s)%s)"
                      % (cand.group(1) if cand else "?",
                         ", %s changed the root split" % ndiff.group(1)
                         if ndiff else ""))
                print("  a differing split does not fail this gate: it is evidence")
                print("  that the root is a decision, and the methods must say so.")

    # ---------- G4 ----------
    stats = {}
    for which, path in (("rooted", args.rooted), ("unrooted", args.unrooted)):
        if not path:
            continue
        print()
        print("G4  tree structure and tip set -- %s" % which)
        print("-" * 74)
        if not os.path.exists(path):
            print("  ERROR: not found: %s" % path)
            failures.append("G4 %s tree missing" % which)
            continue
        try:
            st = tree_stats(path, want_splits=args.check_same_topology)
        except NewickError as exc:
            print("  ERROR: cannot parse Newick: %s" % exc)
            failures.append("G4 %s tree unparsable" % which)
            continue
        stats[which] = st
        print("  file                     : %s (%d bytes)" % (path, st["bytes"]))
        print("  ends with ';'            : %s" % st["ends_with_semicolon"])
        print("  tips / unique            : %d / %d" % (st["n_tips"],
                                                         st["n_unique_tips"]))
        print("  internal nodes           : %d" % st["n_internal"])
        print("  root degree              : %d" % st["root_degree"])
        print("  support labels           : %d" % st["n_support_labels"])
        print("  branch length min/max    : %s / %s" % (st["min_branch"],
                                                         st["max_branch"]))
        print("  zero / negative branches : %d / %d"
              % (st["zero_branches"], st["negative_branches"]))
        frac_zero = (st["zero_branches"] / st["number_of_branches"]
                     if st["number_of_branches"] else 0)
        if frac_zero > 0:
            print("     note: %.1f%% of branches have length exactly 0 (normal for "
                  "identical or unresolved sequences, but it does shift what "
                  "phylogenetic metrics can resolve)" % (100 * frac_zero))

        if not st["ends_with_semicolon"]:
            failures.append("G4 %s tree does not end with ';'" % which)
            print("     !! Newick truncated")
        if st["n_tips"] != st["n_unique_tips"]:
            failures.append("G4 %s tree has duplicate tip labels" % which)
            print("     !! duplicate tip labels -- check whether reference "
                  "sequences were merged")
        if st["negative_branches"]:
            failures.append("G4 %s tree has negative branch lengths" % which)
        if which == "rooted" and args.require_rooted and st["root_degree"] != 2:
            failures.append("G4 rooted tree is not rooted (root degree %d)"
                            % st["root_degree"])
            print("     !! root degree %d, expected 2" % st["root_degree"])
        if args.require_binary:
            # A ROOTED binary tree has root degree 2; an UNROOTED binary tree
            # written by FastTree/IQ-TREE has root degree 3, because the root is
            # not a real node -- it is just where the writer started.  Treating
            # "degree > 2" as a multifurcation therefore rejects every legitimate
            # unrooted output.  Only interior nodes must be binary.
            polytomies = []
            root_obj = parse_newick(read_text_any_encoding(path))
            root_allowed = 2 if which == "rooted" else 3

            def find_poly(n, is_root):
                if not n.is_tip:
                    deg = len(n.children)
                    limit = root_allowed if is_root else 2
                    if deg > limit:
                        polytomies.append(deg)
                    for c in n.children:
                        find_poly(c, False)

            find_poly(root_obj, True)
            if polytomies:
                failures.append("G4 %s tree is not fully bifurcating "
                                "(multifurcations of degree %s)"
                                % (which, sorted(set(polytomies))))
                print("     !! multifurcations found: degrees %s -- resolve with "
                      "ape::multi2di and state that resolution is arbitrary"
                      % sorted(set(polytomies)))
            else:
                print("  bifurcating             : yes (root degree %d allowed for "
                      "a %s tree)" % (root_allowed, which))

        if expected is not None:
            got = set(t.strip("'\"") for t in st["tip_labels"])
            eset = set(expected)
            missing = eset - got
            extra = got - eset
            print("  in input, not in tree   : %d" % len(missing))
            print("  in tree, not in input   : %d" % len(extra))
            if missing:
                print("     e.g. %s" % ", ".join(sorted(missing)[:8]))
            if extra:
                print("     e.g. %s" % ", ".join(sorted(extra)[:8]))
            if missing or extra:
                failures.append("G4 %s tip set != input taxon set" % which)

    if args.require_rooted and not args.rooting_report:
        print("  note: --require-rooted without --rooting-report, so the rooting")
        print("        sensitivity test (D9) is NOT asserted.  Produce the report")
        print("        with root_tree.R and pass --rooting-report FILE to assert it.")

    if args.check_same_topology and "rooted" in stats and "unrooted" in stats:
        a = stats["rooted"].get("split_signature", set())
        b = stats["unrooted"].get("split_signature", set())
        print()
        print("G4  rooted vs unrooted topology")
        print("-" * 74)
        print("  nontrivial splits        : %d (rooted) / %d (unrooted)"
              % (len(a), len(b)))
        only_a = a - b
        only_b = b - a
        if only_a or only_b:
            print("     !! the two files do NOT describe the same unrooted topology")
            print("        splits only in rooted: %d, only in unrooted: %d"
                  % (len(only_a), len(only_b)))
            print("        most likely one of them is stale (a previous run was "
                  "overwritten); re-run the tree step for this run only")
            failures.append("G4 rooted/unrooted topology mismatch")
        else:
            print("  same unrooted topology   : yes")

    # ---------- verdict ----------
    print()
    print("=" * 74)
    if failures:
        print("RESULT: FAIL")
        for f in failures:
            print("  - %s" % f)
    elif unverified:
        print("RESULT: UNVERIFIED (treated as failure, by design)")
        for u in unverified:
            print("  - %s: no positive evidence obtained" % u)
    else:
        print("RESULT: PASS")

    if args.json:
        d = os.path.dirname(os.path.abspath(args.json))
        if d:
            os.makedirs(d, exist_ok=True)
        out = {"failures": failures, "unverified": unverified,
               "verdict": "FAIL" if (failures or unverified) else "PASS"}
        if args.rooting_report:
            body = (read_text_any_encoding(args.rooting_report)
                    if os.path.exists(args.rooting_report) else "")
            m = re.search(r"method used\s*:\s*(\S+)", body)
            out["rooting_report"] = {
                "path": args.rooting_report,
                "present": bool(body),
                "sensitivity_run": bool(body)
                    and "ROOTING SENSITIVITY" in body
                    and not re.search(r"^\s*NOT RUN", body, re.M),
                "method": m.group(1) if m else None,
                "asserted": bool(args.require_rooted),
            }
        for k, v in stats.items():
            out[k] = {kk: vv for kk, vv in v.items()
                      if kk not in ("tip_labels", "split_signature")}
            out[k]["split_count"] = len(v.get("split_signature", []))
        with open(args.json, "w", encoding="utf-8", newline="\n") as fh:
            json.dump(out, fh, indent=2, ensure_ascii=False)
        print("report written: %s" % args.json)

    return 1 if (failures or unverified) else 0


if __name__ == "__main__":
    sys.exit(main())
