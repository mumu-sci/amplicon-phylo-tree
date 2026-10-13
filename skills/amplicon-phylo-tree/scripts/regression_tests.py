#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Regression tests: one per defect found by adversarial review.

Every test names the finding it locks down (R2's IDs, plus R1 items where they
overlap), so if a later change reintroduces the bug, the failure says which one.
The tests are black-box: they run the scripts as subprocesses and assert on exit
codes and output, because the bugs that mattered were exactly the ones where the
script *looked* fine.

Run:
    python regression_tests.py            # everything
    python regression_tests.py B1 B3 M2   # a subset

Exit: 0 if all selected tests pass, 1 otherwise.

Why these particular cases exist (they are all reproducible bugs, not theory):
  B1  the topology check keyed a split by one side's hash, so a rooted and an
      unrooted tree describing the SAME split disagreed -> false "file is stale"
      on every real tree.
  B2  --require-binary treated root degree 3 as a multifurcation, so the
      documented acceptance command failed on every legitimate unrooted tree.
  B3  PowerShell writes UTF-16 logs; read as UTF-8 every keyword missed, so a log
      containing a veto string produced a PASS.
  B8  SKILL.md listed a codon %3 check that did not exist.
  M2  a short low-complexity peptide is made only of characters that are also
      IUPAC nucleotide codes, so the data-side call cannot decide it -- and it
      used to be granted a PASS.
"""
import os
import json
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
VT = os.path.join(HERE, "verify_tree.py")
FS = os.path.join(HERE, "frame_strand.py")
AI = os.path.join(HERE, "audit_inputs.py")
RS = os.path.join(HERE, "recommend_strategy.py")
AQ = os.path.join(HERE, "aln_qc.py")
RT = os.path.join(HERE, "root_tree.R")
DT = os.path.join(HERE, "define_taxa.py")
BL = os.path.join(HERE, "blocks.py")
SKILL_MD = os.path.join(HERE, os.pardir, "SKILL.md")
REFS_MD = os.path.join(HERE, os.pardir, "references")

# A stand-in for `hmmsearch`, written in Python so the M9/M10 tests do not need
# HMMER installed.  It writes a real-shaped HMMER3 --domtblout (22 columns, the
# same layout the parser documents) that the coverage criterion reads.
#
# Two things are NOT invented here:
#   * qlen (column 6) and the alignment coordinates (columns 20-21) are derived
#     from the query lengths in the profile, so the coverage fraction follows.
#   * the sequence's own length comes from the FASTA it is handed.
# That matters, because a fake that fabricated a length difference would let the
# coverage criterion "pass" for the wrong reason.  With a profile supplied, the
# sequences can be the SAME length and the only signal left is how much of the
# read the domain alignment covers.
FAKE_HMM = '''import json, os, sys
a = sys.argv[1:]
out = [a[i + 1] for i, x in enumerate(a) if x == "--domtblout"][0]
fa = a[-1]
prof = {}
pp = os.environ.get("DSH_FAKE_HMM_PROFILE")
if pp and os.path.exists(pp):
    prof = json.load(open(pp))
mode = prof.get("mode", "split")
low = set(prof.get("low_qlen_ids", ["s5", "s6", "s7", "s8"]))
hmm_len = int(prof.get("hmm_len", 100))
seqs, name = {}, None
for line in open(fa):
    line = line.rstrip("\\n")
    if line.startswith(">"):
        name = line[1:].strip().split()[0]
        seqs[name] = ""
    elif name:
        seqs[name] += line.strip()
rows = []
for n in seqs:
    L = len(seqs[n])                        # this is a PROTEIN length
    if mode == "split" and n in low:
        qlen = int(prof.get("qlen_low", 60))
    else:
        qlen = int(prof.get("qlen_high", 200))
    ali_from, ali_to = 1, min(L, qlen)      # how much of the sequence the domain covers
    span = ali_to - ali_from + 1
    rows.append(" ".join(str(x) for x in (
        n, "-", L, "Fer4_NifH", "-", qlen, "1e-60", "200.0", "0.0",
        1, 1, "1e-60", "1e-60", "200.0", "0.0", 1, 1,
        hmm_len, 1, hmm_len, ali_from, ali_to, 1, 1, span, ali_from, ali_to,
        0.99, hmm_len, "1", "100")))
open(out, "w").write("\\n".join(rows) + "\\n")
'''


def hmm_profile(tmp, qlen_low, qlen_high, low_ids=None, hmm_len=100):
    """A profile for the stand-in hmmsearch: which sequences get the SHORT query.

    `qlen` here is the SEQUENCE (query) length in the dummy table, which is what the
    coverage fraction is measured against; the domain alignment spans
    min(sequence_length, qlen) of the read.
    """
    prof = {"mode": "split", "qlen_low": qlen_low, "qlen_high": qlen_high,
            "hmm_len": hmm_len}
    if low_ids is not None:
        prof["low_qlen_ids"] = low_ids
    p = os.path.join(tmp, "profile.json")
    with open(p, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(prof, fh)
    return p


def write_fake_hmm(path, mode, profile=None):
    """Write the stand-in hmmsearch and return the path to hand to the code.

    POSIX must be able to exec it, so the script carries a shebang and is chmod
    755.  On Windows a .py file is not directly executable, so a tiny .cmd wrapper
    is written next to it.  (Forgetting the shebang let this test pass locally and
    fail on Linux CI with "[Errno 8] Exec format error" -- precisely the kind of
    platform hole CI exists to catch.)
    """
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("#!%s\n" % PY)
        fh.write(FAKE_HMM % {"mode": mode})
    os.chmod(path, 0o755)
    if profile is not None:
        with open(path[:-3] + "profile.json", "w", encoding="utf-8",
                  newline="\n") as fh:
            json.dump(profile, fh)
    if os.name == "nt":
        cmd = path[:-3] + ".cmd"
        with open(cmd, "w", encoding="ascii", newline="\r\n") as fh:
            fh.write('@"%s" "%s" %%*\n' % (PY, path))
        return cmd
    return path


TESTS = []


def test(fid, what):
    def deco(fn):
        TESTS.append((fid, what, fn))
        return fn
    return deco


def write(path, text, encoding="utf-8"):
    with open(path, "w", encoding=encoding, newline="\n") as fh:
        fh.write(text)


def run(args):
    p = subprocess.run([PY] + args, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT)
    return p.returncode, p.stdout.decode("utf-8", errors="replace")


def run_with(args, **env):
    """Run a script with extra environment variables (used to point audit_inputs
    at the stand-in hmmsearch without needing a shell wrapper)."""
    e = dict(os.environ)
    e.update(env)
    p = subprocess.run([PY] + args, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, env=e)
    return p.returncode, p.stdout.decode("utf-8", errors="replace")


# --------------------------------------------------------------- fixtures ----

NT_ALN = ">A\nACGTACGTACGT\n>B\nACGTACGAACGT\n>C\nACGTACGTTCGT\n>D\nACGTACGTACGA\n"
NT_ALN_12 = NT_ALN
NT_ALN_10 = ">A\nACGTACGTAC\n>B\nACGTACGAAC\n>C\nACGTACGTTC\n>D\nACGTACGTAC\n"

# a log with the FastTree nucleotide model line and NO mismatch warning
CLEAN_LOG = ("FastTree Version 2.1.11 Double precision\n"
             "Nucleotide distances: Jukes-Cantor Joins: balanced\n"
             "ML Model: Generalized Time-Reversible, CAT approximation with 20 "
             "rate categories\n"
             "Gamma(20) LogLk = -123.456 alpha = 0.4\n")
# ... and the same log carrying the protein-model-on-nucleotide veto
VETO_LOG = ("FastTree Version 2.2.0 Double precision\n"
            "Amino acid distances: BLOSUM45 Joins: balanced\n"
            "ML Model: Whelan-And-Goldman, CAT approximation with 20 rate "
            "categories\n"
            "WARNING! 100.0% NUCLEOTIDE CHARACTERS -- IS THIS REALLY A PROTEIN "
            "ALIGNMENT?\n")

IDS4 = "A\nB\nC\nD\n"


# ------------------------------------------------------------------- B1 ------

@test("B1", "same topology in rooted and unrooted form is NOT a stale file")
def t_b1(tmp):
    ids = os.path.join(tmp, "ids.txt")
    cases = [
        ("((A,B),(C,D));\n", "(A,B,(C,D));\n", IDS4),
        ("(((A,B),(C,D)),((E,F),(G,H)));\n", "((A,B),(C,D),((E,F),(G,H)));\n",
         "A\nB\nC\nD\nE\nF\nG\nH\n"),
    ]
    for rooted, unrooted, id_list in cases:
        write(ids, id_list)
        r = os.path.join(tmp, "r.nwk")
        u = os.path.join(tmp, "u.nwk")
        write(r, rooted)
        write(u, unrooted)
        code, out = run([VT, "--rooted", r, "--unrooted", u, "--ids", ids,
                         "--check-same-topology"])
        if code != 0:
            return False, ("expected the topology check to PASS for\n  %s\n  %s\n"
                           "got exit %d:\n%s" % (rooted.strip(), unrooted.strip(),
                                                 code, out))
        if "same unrooted topology   : yes" not in out:
            return False, "the comparison did not report a topology match:\n" + out
    return True, "both 4-tip and 8-tip rooted/unrooted pairs matched"


@test("B1b", "genuinely different topologies ARE still reported")
def t_b1b(tmp):
    ids = os.path.join(tmp, "ids.txt")
    write(ids, "A\nB\nC\nD\nE\nF\n")
    r = os.path.join(tmp, "r.nwk")
    u = os.path.join(tmp, "u.nwk")
    write(r, "(((A,B),C),((D,E),F));\n")
    write(u, "((A,B),(C,(D,(E,F))));\n")
    code, out = run([VT, "--rooted", r, "--unrooted", u, "--ids", ids,
                     "--check-same-topology"])
    if code == 0:
        return False, "different topologies were accepted:\n" + out
    return True, "a real topology difference still fails (exit %d)" % code


# ------------------------------------------------------------------- B2 ------

@test("B2", "a legitimate UNROOTED binary tree (root degree 3) passes --require-binary")
def t_b2(tmp):
    ids = os.path.join(tmp, "ids.txt")
    write(ids, "A\nB\nC\nD\nE\nF\n")
    u = os.path.join(tmp, "u.nwk")
    write(u, "((A:1,B:1):1,(C:1,D:1):1,(E:1,F:1):1);\n")
    code, out = run([VT, "--unrooted", u, "--ids", ids, "--require-binary"])
    if code != 0:
        return False, ("an ordinary FastTree/IQ-TREE unrooted tree was rejected:\n"
                       + out)
    return True, "root degree 3 accepted; interior nodes still required to be binary"


@test("B2b", "an INTERIOR multifurcation is still caught")
def t_b2b(tmp):
    ids = os.path.join(tmp, "ids.txt")
    write(ids, "A\nB\nC\nD\nE\n")
    u = os.path.join(tmp, "u.nwk")
    write(u, "((A:1,B:1,C:1):1,D:1,E:1);\n")     # interior node has degree 3
    code, out = run([VT, "--unrooted", u, "--ids", ids, "--require-binary"])
    if code == 0:
        return False, "an interior multifurcation passed:\n" + out
    return True, "interior degree 3 still fails (exit %d)" % code


# ------------------------------------------------------------------- B3 ------

@test("B3", "a UTF-16 log is READ, not silently ignored")
def t_b3(tmp):
    aln = os.path.join(tmp, "aln.fasta")
    ids = os.path.join(tmp, "ids.txt")
    write(aln, NT_ALN)
    write(ids, IDS4)

    u8 = os.path.join(tmp, "u8.log")
    u16 = os.path.join(tmp, "u16.log")
    write(u8, VETO_LOG)
    with open(u16, "wb") as fh:
        fh.write(VETO_LOG.encode("utf-16"))          # BOM + UTF-16LE, as PowerShell does

    res = {}
    for tag, log in (("utf-8", u8), ("utf-16", u16)):
        if not os.path.exists(log):
            pass
        un = os.path.join(tmp, "u.nwk")
        write(un, "(A,B,(C,D));\n")
        code, out = run([VT, "--aln", aln, "--seqtype", "nt", "--tool", "fasttree",
                         "--tool-log", log, "--ids", ids, "--unrooted", un])
        res[tag] = (code, out)

    c8, o8 = res["utf-8"]
    c16, o16 = res["utf-16"]
    if c8 == 0:
        return False, "the UTF-8 veto log was not caught:\n" + o8
    if c16 == 0:
        return False, ("the UTF-16 log was NOT read: the veto string was missed and "
                       "the run passed.\n" + o16)
    if "VETO" not in o16:
        return False, "UTF-16 log read but no veto reported:\n" + o16
    return True, "UTF-8 and UTF-16 logs now give the same FAIL"


@test("B3b", "a supplied-but-unparsed log is UNVERIFIED, not a pass")
def t_b3b(tmp):
    aln = os.path.join(tmp, "aln.fasta")
    ids = os.path.join(tmp, "ids.txt")
    write(aln, NT_ALN)
    write(ids, IDS4)
    junk = os.path.join(tmp, "junk.log")
    write(junk, "this file contains nothing recognisable at all\n")
    un = os.path.join(tmp, "u.nwk")
    write(un, "(A,B,(C,D));\n")
    code, out = run([VT, "--aln", aln, "--seqtype", "nt", "--tool", "fasttree",
                     "--tool-log", junk, "--ids", ids, "--unrooted", un])
    if code == 0:
        return False, "an unparseable log produced a PASS:\n" + out
    if "UNVERIFIED" not in out:
        return False, "expected UNVERIFIED, got:\n" + out
    return True, "unparsed log -> UNVERIFIED (exit %d)" % code


# ------------------------------------------------------------------- B8 ------

@test("B8", "codon mode rejects a width that is not a multiple of 3")
def t_b8(tmp):
    ids = os.path.join(tmp, "ids.txt")
    log = os.path.join(tmp, "clean.log")
    un = os.path.join(tmp, "u.nwk")
    write(ids, IDS4)
    write(log, CLEAN_LOG)
    write(un, "(A,B,(C,D));\n")

    bad = os.path.join(tmp, "bad.fasta")
    good = os.path.join(tmp, "good.fasta")
    write(bad, NT_ALN_10)      # 10 columns
    write(good, NT_ALN_12)     # 12 columns

    c_bad, o_bad = run([VT, "--aln", bad, "--seqtype", "codon", "--tool", "fasttree",
                        "--tool-log", log, "--ids", ids, "--unrooted", un])
    if c_bad == 0:
        return False, "a 10-column 'codon' alignment passed:\n" + o_bad

    c_good, o_good = run([VT, "--aln", good, "--seqtype", "codon", "--tool",
                          "fasttree", "--tool-log", log, "--ids", ids,
                          "--unrooted", un])
    if c_good != 0:
        return False, "a valid 12-column codon alignment was rejected:\n" + o_good
    return True, "10 columns FAIL, 12 columns PASS"


@test("B8b", "an X in the amino-acid alignment blocks back-translation evidence")
def t_b8b(tmp):
    ids = os.path.join(tmp, "ids.txt")
    log = os.path.join(tmp, "clean.log")
    un = os.path.join(tmp, "u.nwk")
    ev = os.path.join(tmp, "codon.fasta")
    aa = os.path.join(tmp, "aa.fasta")
    write(ids, IDS4)
    write(log, CLEAN_LOG)
    write(un, "(A,B,(C,D));\n")
    write(ev, NT_ALN_12)
    write(aa, ">A\nMKXTL\n>B\nMKTLL\n>C\nMKTLL\n>D\nMKTLL\n")
    code, out = run([VT, "--aln", ev, "--seqtype", "codon", "--tool", "fasttree",
                     "--tool-log", log, "--ids", ids, "--unrooted", un,
                     "--backtranslation-evidence", ev, "--aa-aln", aa])
    if code == 0:
        return False, "an X-containing aa alignment was accepted:\n" + out
    if "only MACSE" not in out:
        return False, "expected the MACSE hint, got:\n" + out
    return True, "X + PAL2NAL/trimAl incompatibility is caught"


# ------------------------------------------------------------------- M2 ------

@test("M2", "a low-complexity peptide is not granted a data-side PASS")
def t_m2(tmp):
    aln = os.path.join(tmp, "pep.fasta")
    ids = os.path.join(tmp, "ids.txt")
    write(aln, ">s1\nMASTKDVWYRHNKT\n>s2\nMASTKDVWYRHNKT\n>s3\nMASTKDVWYRHNKT\n")
    write(ids, "s1\ns2\ns3\n")
    un = os.path.join(tmp, "u.nwk")
    write(un, "(s1,s2,s3);\n")

    code, out = run([VT, "--aln", aln, "--seqtype", "nt", "--ids", ids,
                     "--unrooted", un])
    if code == 0:
        return False, ("a sequence of pure IUPAC-ambiguity letters was accepted as "
                       "nucleotides with no tool evidence:\n" + out)
    if "AMBIGUOUS" not in out and "UNVERIFIED" not in out:
        return False, "expected AMBIGUOUS/UNVERIFIED, got:\n" + out

    code2, out2 = run([VT, "--aln", aln, "--seqtype", "nt", "--ids", ids,
                       "--unrooted", un, "--allow-no-tool-log"])
    if code2 == 0:
        return False, ("--allow-no-tool-log waived an AMBIGUOUS data-side call, so "
                       "the run passed on NO evidence at all:\n" + out2)
    if "does not apply here" not in out2:
        return False, "expected the refusal to explain itself, got:\n" + out2
    return True, "ambiguity is resolved by evidence or not at all"


# ------------------------------------------------------------------- B4 ------

def domtbl_line(seqid, score=200.0):
    """One HMMER3 --domtblout row (14 fields is all the parser needs)."""
    return ("%s - 300 Fer4_NifH - 300 1e-60 %.1f 0.0 1 1 1e-60 1e-60 %.1f\n"
            % (seqid, score, score))


def fasta_count(path):
    n = 0
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def make_pick_fixture(tmp, hit_ids):
    fa = os.path.join(tmp, "in.fasta")
    write(fa, "".join(">s%d\n%s\n" % (i, "ATG" * 20) for i in (1, 2, 3)))
    fdir = os.path.join(tmp, "frames")
    os.makedirs(fdir, exist_ok=True)
    write(os.path.join(fdir, "frame_p1.domtblout"),
          "".join(domtbl_line("s%d" % i) for i in hit_ids))
    return fa, fdir


@test("B4", "pick accepts --domtblout-dir (the form SKILL.md documents)")
def t_b4(tmp):
    fa, fdir = make_pick_fixture(tmp, [1, 2])
    nt = os.path.join(tmp, "nt.fasta")
    code, out = run([FS, "pick", "--fasta", fa, "--mode", "hmm",
                     "--domtblout-dir", fdir, "--out-nt", nt])
    if "unrecognized arguments" in out:
        return False, "--domtblout-dir is still not a recognised option:\n" + out
    if code != 0:
        return False, "expected success, got exit %d:\n%s" % (code, out)
    return True, "--domtblout-dir works (no shell globbing needed)"


@test("B4b", "a missing directory is a usage error, not a traceback")
def t_b4b(tmp):
    fa, _ = make_pick_fixture(tmp, [1])
    code, out = run([FS, "pick", "--fasta", fa, "--mode", "hmm",
                     "--domtblout-dir", os.path.join(tmp, "nope")])
    if "Traceback" in out:
        return False, "raised instead of reporting:\n" + out
    if code != 2:
        return False, "expected exit 2, got %d:\n%s" % (code, out)
    return True, "missing directory reported cleanly (exit 2)"


# ------------------------------------------------------------------- B5 ------

@test("B5", "unresolved sequences are KEPT -- the taxon set never shrinks silently")
def t_b5(tmp):
    fa, fdir = make_pick_fixture(tmp, [1, 2])          # s3 has no hit
    nt = os.path.join(tmp, "nt.fasta")
    aa = os.path.join(tmp, "aa.fasta")
    tsv = os.path.join(tmp, "table.tsv")
    code, out = run([FS, "pick", "--fasta", fa, "--mode", "hmm",
                     "--domtblout-dir", fdir, "--out-nt", nt, "--out-aa", aa,
                     "--out-tsv", tsv])
    if code != 0:
        return False, "expected success, got exit %d:\n%s" % (code, out)
    n_nt = fasta_count(nt)
    n_aa = fasta_count(aa)
    if n_nt != 3 or n_aa != 3:
        return False, ("the unresolved sequence was dropped: --out-nt has %d and "
                       "--out-aa has %d, expected 3 each.\n%s" % (n_nt, n_aa, out))
    body = open(tsv, encoding="utf-8").read()
    if "orf-fallback" not in body:
        return False, ("the fallback is not recorded in the evidence table, so it is "
                       "indistinguishable from an HMM-based call:\n" + body)
    if "still go into the tree" in out and "orf-fallback" not in out:
        return False, "the printed message claims keeping without saying how"
    return True, "3/3 sequences written, s3 marked evidence=orf-fallback"


@test("B5b", "--drop-unresolved drops them AND says so")
def t_b5b(tmp):
    fa, fdir = make_pick_fixture(tmp, [1, 2])
    nt = os.path.join(tmp, "nt.fasta")
    code, out = run([FS, "pick", "--fasta", fa, "--mode", "hmm",
                     "--domtblout-dir", fdir, "--out-nt", nt,
                     "--drop-unresolved"])
    if code != 0:
        return False, "expected success, got exit %d:\n%s" % (code, out)
    if fasta_count(nt) != 2:
        return False, "expected 2 records, got %d:\n%s" % (fasta_count(nt), out)
    if "EXCLUDED" not in out or "tip-set" not in out:
        return False, ("dropping was not announced loudly enough -- the user must be "
                       "told G4 will complain:\n" + out)
    return True, "exclusion is explicit and warns about the tip-set check"


# ------------------------------------------------------------------- B6 ------

@test("B6", "a check that could not run is NOT a pass (non-coding marker, no DB)")
def t_b6(tmp):
    fa = os.path.join(tmp, "in.fasta")
    write(fa, ">a\nACGTACGTACGT\n>b\nACGTACGAACGT\n")
    code, out = run([AI, "--fasta", fa, "--marker-class", "noncoding"])
    if code == 0:
        return False, ("an unrun orientation/homology audit passed -- the gate is "
                       "decorative for non-coding markers:\n" + out)
    if "could not be checked" not in out:
        return False, "expected the unchecked items to be counted as failures:\n" + out
    return True, "unchecked -> FAIL (exit %d)" % code


@test("B6b", "the gap can be recorded, but only with a reason")
def t_b6b(tmp):
    fa = os.path.join(tmp, "in.fasta")
    write(fa, ">a\nACGTACGTACGT\n>b\nACGTACGAACGT\n")
    code, out = run([AI, "--fasta", fa, "--marker-class", "noncoding",
                     "--accept-and-record"])
    if code != 2:
        return False, ("expected a usage error, got exit %d -- a waiver with no "
                       "reason is indistinguishable from an oversight:\n%s"
                       % (code, out))
    code2, out2 = run([AI, "--fasta", fa, "--marker-class", "noncoding",
                       "--accept-and-record", "--override-reason",
                       "orientation checked separately by the lab"])
    if code2 != 0:
        return False, "a recorded override should let the run proceed:\n" + out2
    if "OVERRIDE ACCEPTED" not in out2:
        return False, "the override was not announced:\n" + out2
    return True, "waiver requires and records a reason"


# ------------------------------------------------------------------- B7 ------

@test("B7", "ITS in a community role prints NO de novo tree command")
def t_b7(tmp):
    code, out = run([RS, "--n-seq", "8000", "--median-len", "300",
                     "--marker-class", "its", "--tree-role", "skeleton"])
    if code != 0:
        return False, "expected exit 0, got %d:\n%s" % (code, out)
    low = out.lower()
    for bad in ("iqtree2 -s", "fasttree -nt", "mafft --maxiterate",
                "mafft --localpair", "mafft --auto"):
        if bad in low:
            return False, ("it forbids an ITS tree and then prints the command to "
                           "build one ('%s'):\n%s" % (bad, out))
    if "ROUTE FOR ITS" not in out:
        return False, "expected an explicit routing section:\n" + out
    if "usearch_global" not in out or "epa-ng" not in out:
        return False, "expected the no-tree and insertion routes:\n" + out
    return True, "ITS routes to no-tree / insertion with no de novo command"


@test("B7b", "a normal rRNA marker still gets a de novo route")
def t_b7b(tmp):
    code, out = run([RS, "--n-seq", "25000", "--median-len", "253",
                     "--marker-class", "rrna", "--tree-role", "skeleton"])
    if code != 0:
        return False, "expected exit 0, got %d" % code
    if "iqtree2" not in out or "mafft" not in out:
        return False, "the normal route lost its commands:\n" + out
    return True, "16S V4 still routes to FFT-NS-i + IQ-TREE"


# ------------------------------------------------------------------- B9 ------

@test("B9", "190k tips must not be sent to exact ML with full UFBoot")
def t_b9(tmp):
    code, out = run([RS, "--n-seq", "190000", "--median-len", "100",
                     "--marker-class", "rrna", "--tree-role", "skeleton",
                     "--ram-gb", "64"])
    if code != 0:
        return False, "expected exit 0, got %d:\n%s" % (code, out)
    if "iqtree2 -s" in out and "-m MFP" in out:
        return False, ("190,000 tips were sent to IQ-TREE with model selection and "
                       "UFBoot -- scaling.md says FastTree above 100,000:\n" + out)
    if "fasttree" not in out.lower():
        return False, "expected FastTree:\n" + out
    if "!! n = 190000" not in out:
        return False, "expected a loud scale warning:\n" + out
    return True, "190k -> FastTree + explicit warnings"


@test("B9b", "a mid-size run keeps exact ML but warns above 50k")
def t_b9b(tmp):
    code, out = run([RS, "--n-seq", "60000", "--median-len", "253",
                     "--marker-class", "rrna", "--tree-role", "skeleton",
                     "--ram-gb", "8"])
    if code != 0:
        return False, "expected exit 0, got %d" % code
    if "!! n = 60000" not in out:
        return False, "expected the benchmark warning at 60k:\n" + out
    if "very likely to fail" not in out:
        return False, "expected the low-RAM warning:\n" + out
    return True, "60k warns on both size and RAM"


# ------------------------------------------------------ B10 / M3 / M3b ------
# (R2: the workflow promised things the scripts did not do; these lock in the
#  promised behaviour AND keep it from being silently removed again.)

@test("B10", "D2b documents the back-translation tools' incompatible inputs")
def t_b10(tmp):
    ref = os.path.join(HERE, os.pardir, "references", "parameters.md")
    skill = os.path.join(HERE, os.pardir, "SKILL.md")
    doc = open(ref, encoding="utf-8").read() + open(skill, encoding="utf-8").read()
    # the three tool contracts, and the rule that follows from them
    for needle in ("PAL2NAL", "trimAl", "backtrans", "MACSE", "X"):
        if needle not in doc:
            return False, "the contract text no longer mentions %r" % needle
    if "终止" not in doc and "stop codon" not in doc.lower():
        return False, "the stop-codon incompatibility is not stated"
    if "MACSE" not in doc.split("PAL2NAL")[1]:
        return False, "MACSE is not presented as the only tolerant tool"
    return True, "PAL2NAL/trimAl/MACSE input contracts and the X rule are documented"


@test("M3", "a MIXED nucleotide/protein alphabet is a hard structural FAIL")
def t_m3(tmp):
    ids = os.path.join(tmp, "ids.txt")
    write(ids, IDS4)
    un = os.path.join(tmp, "u.nwk")
    write(un, "(A,B,(C,D));\n")
    mixed = os.path.join(tmp, "mixed.fasta")
    write(mixed, ">A\nACGTEFLACGT\n>B\nACGTACGAGGGT\n>C\nACGTACGTTCGT\n"
                 ">D\nACGTACGTACGA\n")
    code, out = run([AQ, "--aln", mixed])
    if code == 0:
        return False, ("ACGTU and EFILPQZJO in one alignment was accepted, yet no "
                       "single substitution model can describe it:\n" + out)
    if "MIXED ALPHABET" not in out and "alphabets are mixed" not in out:
        return False, "expected an explicit mixed-alphabet report:\n" + out

    # reverse control: the same alignment without the protein-only letters passes
    clean = os.path.join(tmp, "clean.fasta")
    write(clean, NT_ALN)
    code2, out2 = run([AQ, "--aln", clean])
    if code2 != 0:
        return False, "a pure-nucleotide alignment was rejected:\n" + out2

    # and an ILLEGAL character is a structural error too
    bad = os.path.join(tmp, "illegal.fasta")
    write(bad, ">A\nACGT@CGT\n>B\nACGTACGT\n>C\nACGTACGT\n>D\nACGTACGT\n")
    code3, out3 = run([AQ, "--aln", bad])
    if code3 == 0:
        return False, "an illegal character '@' was accepted:\n" + out3
    if "ILLEGAL" not in out3:
        return False, "expected an explicit illegal-character report:\n" + out3
    return True, ("mixed alphabet FAILs, an illegal character FAILs, and a clean alignment still passes")


@test("M3b", "--ambiguous-policy is implemented, and drop-col never edits the file")
def t_m3b(tmp):
    amb = os.path.join(tmp, "amb.fasta")
    write(amb, ">A\nACGNACGTNNNN\n>B\nACGAACGTACGT\n>C\nACGCACGTACGT\n"
               ">D\nACGTACGTACGT\n")
    before = open(amb, "rb").read()

    for policy in ("gap", "missing", "drop-col"):
        code, out = run([AQ, "--aln", amb, "--ambiguous-policy", policy,
                         "--json", os.path.join(tmp, "j_%s.json" % policy)])
        if "unrecognized arguments" in out:
            return False, "--ambiguous-policy %s is not a recognised option:\n%s" \
                % (policy, out)
        if code not in (0, 1):
            return False, "--ambiguous-policy %s crashed (exit %d):\n%s" \
                % (policy, code, out)

    # the file must be byte-identical afterwards: drop-col suggests, it never edits
    if open(amb, "rb").read() != before:
        return False, "--ambiguous-policy drop-col MODIFIED the alignment file"
    lp = os.path.join(tmp, "amb.fasta.drop_cols.txt")
    if not os.path.exists(lp):
        return False, "drop-col did not write the suggested-columns list"
    if "SUGGESTION ONLY" not in open(lp, encoding="utf-8").read():
        return False, "the suggestion file does not say it is a suggestion only"

    # reverse control: 'missing' must count the ambiguous characters as gaps
    def gapf(policy):
        import json
        p = os.path.join(tmp, "j_%s.json" % policy)
        return json.load(open(p, encoding="utf-8"))["gap_fraction"]
    if not gapf("missing") > gapf("gap"):
        return False, ("'missing' must count N as a gap: got %g vs %g"
                       % (gapf("missing"), gapf("gap")))
    return True, ("all three policies are accepted, the file is untouched, and missing counts more gaps than gap")


@test("M4", "an override needs a reason AND lands in the JSON, not only on screen")
def t_m4(tmp):
    # 60 of the 72 matrix cells are gaps -> above the 0.80 nucleotide limit
    aln = os.path.join(tmp, "gappy.fasta")
    write(aln, ">A\nACGT--------\n>B\nACGTACGT----\n>C\n------------\n"
               ">D\n------------\n>E\n------------\n>F\n------------\n")
    j = os.path.join(tmp, "aln.json")

    code, out = run([AQ, "--aln", aln, "--accept-and-record", "--json", j])
    if code != 2:
        return False, ("--accept-and-record without --override-reason gave exit %d; "
                       "an unrecorded waiver is indistinguishable from a pass:\n%s"
                       % (code, out))

    code2, out2 = run([AQ, "--aln", aln, "--accept-and-record",
                       "--override-reason", "reference alignment is sparse",
                       "--json", j])
    if code2 != 0:
        return False, ("a recorded override should let the ratio gate through, got "
                       "exit %d:\n%s" % (code2, out2))
    if "gap fraction" not in out2:
        return False, "the ratio issue was not reported at all:\n" + out2
    import json
    rec = json.load(open(j, encoding="utf-8"))
    ov = rec.get("gate_overrides") or []
    if not ov:
        return False, "the override is still not in the machine-readable record"
    if not all(o.get("reason") and o.get("metric") and o.get("limit") is not None
               for o in ov):
        return False, "a gate_overrides entry is missing gate/metric/limit/reason"
    if rec.get("override_reason") != "reference alignment is sparse":
        return False, "the reason was not carried into the report"

    # reverse control: an override must NOT cover a structural error
    ragged = os.path.join(tmp, "ragged.fasta")
    j2 = os.path.join(tmp, "ragged.json")
    write(ragged, ">A\nACGT\n>B\nACGTACGT\n>C\nACGTACGT\n>D\nACGTACGT\n")
    code3, out3 = run([AQ, "--aln", ragged, "--accept-and-record",
                       "--override-reason", "sparse", "--json", j2])
    if code3 == 0:
        return False, "a structural error (not rectangular) was overridden:\n" + out3
    rec3 = json.load(open(j2, encoding="utf-8"))
    if rec3.get("gate_overrides"):
        return False, "a structural failure still recorded a ratio override"
    return True, ("the override needs a reason, lands in gate_overrides, and does "
                  "not cover a structural error")


# -------------------------------------------------------------- M9 / M10 -----

def audit_fixture(tmp, lengths, permissive=True):
    """Write a FASTA (ids s1..sN) plus the args needed for a self-contained run."""
    fa = os.path.join(tmp, "in.fasta")
    with open(fa, "w", encoding="utf-8", newline="\n") as fh:
        for i, L in enumerate(lengths, 1):
            fh.write(">s%d\n%s\n" % (i, "ACGT" * (L // 4)))
    return fa, os.path.join(tmp, "dummy.hmm")


def passes(out):
    """True when the audit reported a clean pass (no override recorded)."""
    return "RESULT: PASS" in out and "unchecked" not in out


@test("M9", "equal-length sequences with a split domain-coverage pattern are caught")
def t_m9(tmp):
    j = os.path.join(tmp, "a.json")
    # EVERY read is 240 bp.  The length-histogram criterion therefore cannot fire,
    # so if the mixture is caught at all it is the coverage channel that caught it.
    # The profile says the domain HMM is long for one half of the reads and short
    # for the other, i.e. the same primer pair would have to amplify two different
    # parts of the marker.
    fa, hmm = audit_fixture(tmp, [240] * 8)
    prof = hmm_profile(tmp, qlen_low=150, qlen_high=320)
    fake = write_fake_hmm(os.path.join(tmp, "fake_split.py"), "split", prof)
    code, out = run_with([AI, "--fasta", fa, "--marker-class", "coding",
                          "--hmm", hmm, "--json", j],
                         DSH_HMMSEARCH=fake, DSH_FAKE_HMM_PROFILE=prof)
    if code == 0:
        return False, ("an equal-length mixture with a split domain-coverage "
                       "pattern was passed: the coverage criterion is not wired "
                       "in:\n" + out)
    if "length modes" not in out or "only one populated bin" not in out:
        return False, ("the length criterion DID fire, so this case does not "
                       "isolate the coverage channel:\n" + out)
    if "coverage" not in out.lower():
        return False, "the failure does not name the coverage criterion:\n" + out
    rec = json.load(open(j, encoding="utf-8"))
    crit = rec["mixed_region"]["criteria"]
    if crit.get("hmm_domain_coverage") is None:
        return False, "the coverage numbers are not in the machine-readable report"
    cov = crit["hmm_domain_coverage"]
    if cov["n_low"] < 1 or cov["n_high"] < 1:
        return False, "the coverage split was not reported as two groups"
    if cov["len_difference"] != 0:
        return False, ("the two coverage groups differ in length (%s), so the "
                       "length criterion could have been the trigger instead"
                       % cov["len_difference"])

    # reverse control: a uniform domain model over the same reads must NOT be
    # called a mixture, or the criterion fires on nothing
    prof2 = hmm_profile(tmp, qlen_low=300, qlen_high=300)
    fake2 = write_fake_hmm(os.path.join(tmp, "fake_uniform.py"), "split", prof2)
    j2 = os.path.join(tmp, "b.json")
    code2, out2 = run_with([AI, "--fasta", fa, "--marker-class", "coding",
                            "--hmm", hmm, "--json", j2],
                           DSH_HMMSEARCH=fake2, DSH_FAKE_HMM_PROFILE=prof2)
    if not passes(out2):
        return False, ("a uniform coverage pattern was still rejected -- the "
                       "criterion fires on nothing:\n" + out2)
    return True, ("an equal-length split is caught by the coverage channel alone, "
                  "and a uniform pattern is not")


@test("M10", "--region-map is read, validated, and recorded with a hash")
def t_m10(tmp):
    # two length classes so the mixed-region verdict is unambiguous; the region map
    # is what this test is about
    fa, hmm = audit_fixture(tmp, [200] * 4 + [80] * 4)
    prof = hmm_profile(tmp, qlen_low=300, qlen_high=300)
    fake = write_fake_hmm(os.path.join(tmp, "fake.py"), "split", prof)
    env = {"DSH_HMMSEARCH": fake, "DSH_FAKE_HMM_PROFILE": prof}
    import json

    one_col = os.path.join(tmp, "one.map")
    write(one_col, "s1\ns2\n")
    code, out = run_with([AI, "--fasta", fa, "--marker-class", "coding",
                          "--hmm", hmm, "--mixed-region-ok",
                          "--region-map", one_col], **env)
    if code != 2:
        return False, ("a one-column region map gave exit %d, expected 2:\n%s"
                       % (code, out))

    unknown = os.path.join(tmp, "unknown.map")
    write(unknown, "".join("s%d\tR1\n" % i for i in range(1, 8)) + "zzz\tR1\n")
    code2, out2 = run_with([AI, "--fasta", fa, "--marker-class", "coding",
                            "--hmm", hmm, "--mixed-region-ok",
                            "--region-map", unknown], **env)
    if code2 != 2:
        return False, ("a region map naming an id absent from the FASTA gave exit "
                       "%d, expected 2:\n%s" % (code2, out2))

    good = os.path.join(tmp, "good.map")
    write(good, "".join("s%d\tR%d\n" % (i, 1 + (i > 4)) for i in range(1, 9)))
    j = os.path.join(tmp, "good.json")
    code3, out3 = run_with([AI, "--fasta", fa, "--marker-class", "coding",
                            "--hmm", hmm, "--mixed-region-ok",
                            "--region-map", good, "--json", j], **env)
    if code3 != 0:
        return False, "a valid two-column region map was rejected:\n" + out3
    rm = json.load(open(j, encoding="utf-8")).get("region_map")
    if not rm:
        return False, "the region map is still decorative: nothing in the report"
    if len(rm.get("sha256", "")) != 64 or rm.get("n_regions") != 2:
        return False, "the recorded region map lacks its sha256 / region count"
    if rm["regions"] != {"R1": 4, "R2": 4}:
        return False, "the per-region counts are wrong: %r" % (rm["regions"],)

    # reverse control: a waiver WITHOUT a map has no traceability and must fail
    code4, out4 = run_with([AI, "--fasta", fa, "--marker-class", "coding",
                            "--hmm", hmm, "--mixed-region-ok",
                            "--json", os.path.join(tmp, "nomap.json")], **env)
    if code4 == 0:
        return False, ("--mixed-region-ok without --region-map passed, although the "
                       "region labels are recorded nowhere:\n" + out4)
    return True, ("map read and validated, sha256 and region counts recorded, and a "
                  "map-less waiver fails")


# ------------------------------------------------------------------ M6/M7 ----

def _find_rscript():
    """Rscript via DSH_RSCRIPT or PATH only -- no developer machine path is baked
    into a published skill."""
    for cand in (os.environ.get("DSH_RSCRIPT"), shutil.which("Rscript"),
                 shutil.which("Rscript.exe")):
        if cand and os.path.exists(cand):
            return cand
    return None


@test("M7", "root_tree.R writes one file per candidate root (contract in the source)")
def t_m7(tmp):
    src = open(RT, encoding="utf-8").read()
    if "cand%d.nwk" not in src and ".cand" not in src:
        return False, "root_tree.R no longer writes <rooted>.candN.nwk files"
    if "--no-sensitivity" not in src:
        return False, "the --no-sensitivity switch is gone from root_tree.R"
    if 'sensitivity = "0"' not in src and 'sensitivity = "2"' not in src:
        return False, "the sensitivity default is no longer written explicitly"
    if "ROOTING SENSITIVITY" not in src:
        return False, "the report no longer has a ROOTING SENSITIVITY section"
    return True, "candidate files + --no-sensitivity + the default are in the source"


@test("M6", "rooting defaults to a sensitivity test, and NOT RUN is a FAIL")
def t_m6(tmp):
    aln = os.path.join(tmp, "aln.fasta")
    write(aln, NT_ALN)
    ids = os.path.join(tmp, "ids.txt")
    write(ids, IDS4)
    rooted = os.path.join(tmp, "r.nwk")
    unrooted = os.path.join(tmp, "u.nwk")
    write(rooted, "((A,B),(C,D));\n")
    write(unrooted, "(A,B,(C,D));\n")

    # a report whose sensitivity test ran: the gate is satisfied
    ran = os.path.join(tmp, "ran.txt")
    write(ran, "3  ROOTING\n  method used               : midpoint\n"
               "\n4  ROOTING SENSITIVITY\n  candidates compared       : 2 (up "
               "to 4 exist), seed 20261012\n  1 of 2 candidate root positions "
               "give a DIFFERENT split across\n")
    code, out = run([VT, "--rooted", rooted, "--unrooted", unrooted, "--ids", ids,
                     "--require-rooted", "--rooting-report", ran])
    if code != 0:
        return False, "an executed sensitivity test was rejected:\n" + out

    # ... and a report that says NOT RUN is a gap, not a pass
    notrun = os.path.join(tmp, "notrun.txt")
    write(notrun, "3  ROOTING\n  method used               : midpoint\n"
                  "\n4  ROOTING SENSITIVITY\n  NOT RUN (--no-sensitivity / "
                  "--sensitivity-candidates 0).\n")
    code2, out2 = run([VT, "--rooted", rooted, "--unrooted", unrooted, "--ids", ids,
                       "--require-rooted", "--rooting-report", notrun])
    if code2 == 0:
        return False, ("a root whose sensitivity was never tested passed:\n" + out2)
    if "NOT RUN" not in out2:
        return False, "the failure does not name what was not run:\n" + out2

    # reverse control: without --rooting-report the older behaviour is unchanged
    code3, out3 = run([VT, "--rooted", rooted, "--unrooted", unrooted, "--ids", ids,
                       "--require-rooted"])
    if code3 != 0:
        return False, ("--require-rooted without --rooting-report must stay valid "
                       "(it now only prints a note):\n" + out3)

    # and the real R script is exercised when R is available (static contract
    # otherwise): the default must produce this very report shape
    r = _find_rscript()
    if r:
        t = os.path.join(tmp, "t.nwk")
        write(t, "((((A:0.1,B:0.1):0.05,C:0.12):0.03,(D:0.2,E:0.2):0.06):0.02,"
                 "((F:0.11,G:0.13):0.04,(H:0.3,(I:0.05,J:0.05):0.07):0.05):0.01);\n")
        rep = os.path.join(tmp, "rooting.txt")
        p = subprocess.run([r, RT, "--in", t,
                            "--rooted", os.path.join(tmp, "r2.nwk"),
                            "--unrooted", os.path.join(tmp, "u2.nwk"),
                            "--report", rep, "--label", "regress"],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        text = p.stdout.decode("utf-8", errors="replace")
        if "ape" in text and "not installed" in text:
            return True, ("(R ran, but ape/phangorn are missing) the Python side is "
                          "verified; the R report was not produced")
        if p.returncode != 0 or "RESULT: PASS" not in text:
            return False, "the real root_tree.R run failed:\n" + text[-1500:]
        body = open(rep, encoding="utf-8").read()
        if "NOT RUN" in body:
            return False, ("root_tree.R still defaults to NOT RUN:\n"
                           + body[body.find("4  ROOTING"):][:400])
        if ".cand1.nwk" not in body:
            return False, "no candidate root tree was reported as written"
        if not os.path.exists(os.path.join(tmp, "r2.nwk.cand1.nwk")):
            return False, "the candidate-root file does not exist on disk"
    return True, ("the Python gate fails on NOT RUN, passes on a real run, and the real R script is exercised when R is present")


# --------------------------------------------------- phase-2 revisions ------
# (the drift tool, the two-dimensional table, and the document contracts that a
#  later edit could quietly soften)

@test("M11", "check_consistency.py passes here and FAILS on a planted drift")
def t_m11(tmp):
    CC = os.path.join(HERE, "check_consistency.py")
    code, out = run([CC])
    if code != 0:
        return False, ("the repository is not self-consistent -- run "
                       "check_consistency.py and fix what it prints:\n" + out)

    # Plant a drift in a throwaway COPY of the skill: scaling.md now promises
    # E-INS-i where the code returns L-INS-i.  The checker must catch it, which is
    # what makes it more than a script that always says "ok".
    skill = os.path.dirname(HERE)
    stage = os.path.join(tmp, "skill")
    shutil.copytree(skill, stage,
                    ignore=shutil.ignore_patterns("__pycache__", "assets", "agents"))
    scal = os.path.join(stage, "references", "scaling.md")
    text = open(scal, encoding="utf-8").read()
    # the N=200-3000 row, L<=2000 cell: L-INS-i -> E-INS-i (works on the Chinese
    # original and on the English translation)
    for old, new in (
            ("**L-INS-i**：`--localpair --maxiterate 1000`",
             "**E-INS-i**：`--genafpair --maxiterate 1000`"),
            ("**L-INS-i**: `--localpair --maxiterate 1000`",
             "**E-INS-i**: `--genafpair --maxiterate 1000`")):
        broken = text.replace(old, new, 1)
        if broken != text:
            break
    if broken == text:
        return False, "could not find the cell to break in scaling.md"
    write(scal, broken)
    code2, out2 = run([os.path.join(stage, "scripts", "check_consistency.py")])
    if code2 == 0:
        return False, ("a planted scaling.md/code mismatch was NOT caught, so the "
                       "checker cannot detect the drift it exists for:\n" + out2)
    if "scaling.md says" not in out2:
        return False, "the drift was caught but not named:\n" + out2

    # A flag documented but not implemented must also be caught (R2's B4/M3).
    sk2 = os.path.join(stage, "SKILL.md")
    text2 = open(sk2, encoding="utf-8").read()
    write(sk2, text2 + "\n\nUse `--totally-made-up-flag` to enable the thing.\n")
    code3, out3 = run([os.path.join(stage, "scripts", "check_consistency.py")])
    if code3 == 0 or "--totally-made-up-flag" not in out3:
        return False, ("a documented-but-nonexistent flag was not caught:\n" + out3)
    return True, ("the checker passes here and catches both drift classes on a "
                  "planted copy")


@test("M1", "the BLOSUM counts are the corrected ones everywhere they appear")
def t_m1(tmp):
    """M1: 'BLOSUM appears 10 times' was a false measured fact.  The real counts are
    case-sensitive: uppercase BLOSUM 0, case-insensitive Blosum 36.

    The regression that matters is the CONCLUSION (do not blacklist BLOSUM) PLUS
    the numbers, because the wrong number is how the conclusion got plausible.
    """
    files = [SKILL_MD, os.path.join(REFS_MD, "pitfalls.md")]
    for path in files:
        text = open(path, encoding="utf-8").read()
        base = os.path.basename(path)
        if "BLOSUM" not in text:
            return False, "%s no longer discusses the BLOSUM blacklist" % base
    if re.search(r"出现\s*10\s*次", text) or "实测 10 次" in text \
            or "appears 10 times" in text:
        return False, "%s still carries the wrong '10 times' count" % base
    if "0" not in text or "36" not in text:
        return False, ("%s does not give the corrected counts (0 / 36)" % base)
        if "大小写" not in text and "case-insensitive" not in text.lower() \
                and "case" not in text.lower():
            return False, "%s lost the case-sensitivity rule" % base
    return True, "wrong count gone; 0/36 and the case-sensitivity rule are present"


@test("M12", "the grey zone has a usable recipe, not just the words 'build both'")
def t_m12(tmp):
    text = open(REFS_MD + os.sep + "parameters.md", encoding="utf-8").read()
    need = {
        "p_aa definition": "p_aa",
        "gapped-site rule": "gapped",
        "median vs mean": "median",
        "RF recipe": "RF.dist",
        "same tip set": "tip",
        "same region": "region",
        "branch lengths only after back-translation": "back-translat",
        "phangorn": "phangorn",
    }
    for label, needle in need.items():
        if needle.lower() not in text.lower():
            return False, "D2b no longer states the %s (%r is missing)" % (label,
                                                                          needle)
    if "需引用出处" in text:
        return False, ("D2b still claims a single citation exists for 85/90; the "
                       "honest position is 'convention starting points'")
    if "锚点示例" not in text and "惯例起点" not in text \
            and "anchor example" not in text.lower() \
            and "convention" not in text.lower():
        return False, "D2b does not say the 85/90 numbers are starting points"
    return True, "p_aa is defined, the RF recipe is written down, and 85/90 is honest"


@test("M13", "the organelle family states its interpretive ceiling")
def t_m13(tmp):
    text = open(REFS_MD + os.sep + "marker-types.md", encoding="utf-8").read()
    if "可解释上限" not in text and "Interpretive ceiling" not in text:
        return False, "marker-types.md lost the 'interpretive ceiling' rows"
    # the organelle section must carry one, with the numbers that make it usable
    sec = text.split("`organelle`")[-1].split("## 5.")[0] if "`organelle`" in text \
        else text
    for needle in ("313 bp", "104", "numt"):
        if needle not in sec:
            return False, ("the organelle ceiling is missing %r -- it is the row "
                           "outputs.md requires at delivery time" % needle)
    if "not comparable" not in sec.lower() and "cannot be compared" not in sec.lower() \
            and "不可比" not in sec:
        return False, ("the organelle ceiling no longer says short-fragment de novo "
                       "branch lengths are not comparable with a reference tree")
    # and the COI direction must not be stated backwards any more
    for backwards in ("COI 常因第三位饱和走氨基酸",
                      "COI usually goes amino-acid because the third position is "
                      "saturated"):
        if backwards in text:
            return False, ("marker-types.md still says COI usually goes amino-acid, "
                           "which is backwards for a 313 bp within-genus comparison")
    return True, "the organelle ceiling and the corrected COI direction are present"


@test("M5", "gap thresholds follow the marker family, and the escape hatch is printed")
def t_m5(tmp):
    aln = os.path.join(tmp, "g.fasta")
    # 5 gaps in 12 columns x 3 sequences = 20/36 = 0.556 -> passes everywhere
    write(aln, ">A\nACGTACGTACGT\n>B\nACGTACTTTTTT\n>C\nACGTAC------\n")
    j = os.path.join(tmp, "its.json")
    code, out = run([AQ, "--aln", aln, "--marker-class", "its", "--json", j])
    if code != 0:
        return False, "a gappy ITS alignment was rejected at the ITS threshold:\n" + out
    import json
    rec = json.load(open(j, encoding="utf-8"))
    if abs(rec.get("gap_limit", 0) - 0.90) > 1e-9:
        return False, ("ITS did not get the 0.90 nucleotide default (got %r)"
                       % rec.get("gap_limit"))
    if rec.get("gap_limit_source") != "marker-class nucleotide default":
        return False, "the report does not say where the threshold came from"

    # reverse control: an explicit --max-gap must win over the family default
    code2, out2 = run([AQ, "--aln", aln, "--marker-class", "its",
                       "--max-gap", "0.10"])
    if code2 == 0:
        return False, "an explicit --max-gap was overridden by the family default"
    if "--accept-and-record" not in out2 or "--override-reason" not in out2:
        return False, ("the ratio FAIL does not print the copy-pasteable "
                       "recorded-override command:\n" + out2)
    if "RESULT: FAIL" not in out2:
        return False, "expected a FAIL for the explicit low threshold"
    return True, ("ITS gets 0.90, an explicit --max-gap wins, and the FAIL prints "
                  "the recorded-override command")


@test("m1", "the taxon rule grammar covers nonzero, and a bad rule is a usage error")
def t_m1b(tmp):
    table = os.path.join(tmp, "t.tsv")
    write(table, "#OTU ID\ts1\ts2\ts3\tTaxonomy\n"
                 "A\t5000\t0\t7000\tx\nB\t1000\t0\t0\ty\nC\t3000\t3000\t0\tz\n")
    ids = os.path.join(tmp, "ids.txt")
    code, out = run([DT, "--table", table, "--rule", "nonzero>=2",
                     "--out-ids", ids])
    if code != 0:
        return False, "the documented 'nonzero>=N' presence form is not accepted:\n" + out
    kept = [l.strip() for l in open(ids, encoding="utf-8") if l.strip()]
    if kept != ["A", "C"]:
        return False, ("nonzero>=2 selected %r, expected the two ASVs present in "
                       ">=2 samples" % kept)

    # reverse control: `prev>=2` must select exactly the same set
    ids2 = os.path.join(tmp, "ids2.txt")
    run([DT, "--table", table, "--rule", "prev>=2", "--out-ids", ids2])
    if [l.strip() for l in open(ids2, encoding="utf-8") if l.strip()] != kept:
        return False, "nonzero and prev disagree; they are the same quantity"

    # and an unparseable rule must not raise a traceback
    code3, out3 = run([DT, "--table", table, "--rule", "sum>>10",
                       "--out-ids", os.path.join(tmp, "ids3.txt")])
    if "Traceback" in out3:
        return False, "an unparseable rule still raises:\n" + out3
    if code3 != 1 or "could not parse" not in out3:
        return False, ("expected a clean exit 1 with a message, got exit %d:\n%s"
                       % (code3, out3))
    return True, "nonzero == prev, and a bad rule is reported, not raised"


@test("m2", "the no-table route says so in the machine-readable record")
def t_m2b(tmp):
    fa = os.path.join(tmp, "a.fasta")
    write(fa, ">s1\nACGT\n>s2\nACGTACGT\n")
    j = os.path.join(tmp, "r.json")
    code, out = run([DT, "--fasta-only", fa, "--out-ids",
                     os.path.join(tmp, "ids.txt"), "--json", j])
    if code != 0:
        return False, "the documented --fasta-only path failed:\n" + out
    import json
    rec = json.load(open(j, encoding="utf-8"))
    if rec.get("abundance_table") is not None:
        return False, "the record still claims an abundance table"
    if "must_not_claim" not in rec or "abundance" not in rec["must_not_claim"]:
        return False, ("the machine-readable record does not forbid the abundance "
                       "claim:\n" + out)

    # reverse control: with a table, the field is present but empty
    table = os.path.join(tmp, "t.tsv")
    write(table, "#OTU ID\ts1\ts2\nA\t5\t7\nB\t1\t0\n")
    j2 = os.path.join(tmp, "r2.json")
    code2, out2 = run([DT, "--table", table, "--rule", "sum>=5",
                       "--out-ids", os.path.join(tmp, "i2.txt"), "--json", j2,
                       "--assume-counts"])
    if code2 != 0:
        return False, "the table route failed:\n" + out2
    rec2 = json.load(open(j2, encoding="utf-8"))
    if rec2.get("mode") != "table" or rec2.get("must_not_claim"):
        return False, "the table route did not record itself correctly"
    return True, "fasta-only records the absence of a table; the table route does not"


@test("m3", "blocks.py concat honours --gap-char and rejects an invalid one")
def t_m3c(tmp):
    # each sequence belongs to exactly ONE block (that is the concat contract), so
    # the all-gap filler appears for the sequences of the OTHER block
    b1 = os.path.join(tmp, "block_1.fasta")
    b2 = os.path.join(tmp, "block_2.fasta")
    write(b1, ">a\nACGT\n>b\nAC-T\n")
    write(b2, ">c\nCCGG\n>d\nC-GG\n")
    out1 = os.path.join(tmp, "super.fasta")
    code, so = run([BL, "concat", "--aln-glob", os.path.join(tmp, "block_*.fasta"),
                    "--out", out1])
    if code != 0:
        return False, "plain concat failed:\n" + so
    body = open(out1, encoding="utf-8").read()
    if ">c\n----CCGG" not in body.replace("\r", ""):
        return False, ("the default '-' filler was not written for the other "
                       "block's sequences:\n" + body)

    out2 = os.path.join(tmp, "super2.fasta")
    code2, so2 = run([BL, "concat", "--aln-glob",
                      os.path.join(tmp, "block_*.fasta"),
                      "--out", out2, "--gap-char", "."])
    if code2 != 0:
        return False, "--gap-char was rejected:\n" + so2
    body2 = open(out2, encoding="utf-8").read().replace("\r", "")
    if ">c\n....CCGG" not in body2:
        return False, "--gap-char . was accepted but not used:\n" + body2
    if "gap filler" not in so2:
        return False, "the filler in force was not reported:\n" + so2

    # reverse control: an unrecognised filler is a usage error, not silent corruption
    code3, so3 = run([BL, "concat", "--aln-glob",
                      os.path.join(tmp, "block_*.fasta"),
                      "--out", os.path.join(tmp, "super3.fasta"),
                      "--gap-char", "N"])
    if code3 == 0:
        return False, ("--gap-char N was accepted, although N is read as a STATE "
                       "and would corrupt the matrix:\n" + so3)
    if "not a recognised" not in so3:
        return False, "the invalid filler was not explained:\n" + so3
    return True, "fillers are honoured and an invalid one is refused"


@test("m4", "the two routes the docs arbitrate are actually arbitrated")
def t_m4b(tmp):
    mk = open(REFS_MD + os.sep + "marker-types.md", encoding="utf-8").read()
    sk = open(SKILL_MD, encoding="utf-8").read()
    # 18S: stratified trees are allowed for a conclusion, forbidden for a skeleton
    if "Arbitration rule" not in mk and "仲裁规则" not in mk:
        return False, "the 18S arbitration rule is gone from marker-types.md"
    if "Arbitration rule" in mk:
        seg = mk.split("Arbitration rule")[1][:1600]
    else:
        seg = mk.split("仲裁规则")[1][:1600]
    if "Distance skeleton" not in seg and "距离骨架" not in seg:
        return False, "the arbitration no longer names the distance-skeleton case"
    if "FORBIDDEN to stratify" not in seg and "禁止分层" not in seg:
        return False, "the arbitration no longer forbids stratification for metrics"
    if "The conclusion itself" not in seg and "结论本身" not in seg:
        return False, "the arbitration no longer mentions the conclusion case"
    if "is allowed" not in seg and "允许" not in seg:
        return False, "the arbitration no longer permits it for a conclusion"
    # ITS: the absolute wording must stay limited to cross-genus
    if "跨属以上禁止" not in mk and "cross-genus and above" not in mk:
        return False, "the ITS hard constraint lost its cross-genus limitation"
    for needle in ("RPB2", "TEF1", "tubulin", "ITS2 database", "4SALE",
                   "RNAsalsa", "ITSx"):
        if needle not in mk:
            return False, "the ITS alternatives/tools list lost %r" % needle
    # and the skeleton rule itself must still be stated where the route is chosen
    if "必须保留分析矩阵的全部类群" not in sk and "MUST be kept" not in sk \
            and "keep EVERY taxon" not in sk:
        return False, "SKILL.md no longer states the keep-every-taxon rule"
    return True, "18S and ITS are arbitrated in writing, with the alternatives named"


# ------------------------------------------------------------------ runner ---



def main(argv):
    # Error messages quote paths and tool output, which may contain characters the
    # console codepage cannot encode (a mangled path, a UTF-8 log).  Printing must
    # never be the reason a test run dies.
    try:
        sys.stdout.reconfigure(errors="replace")
    except (AttributeError, ValueError):
        pass
    want = [a.upper() for a in argv[1:]]
    sel = [t for t in TESTS
           if not want or any(t[0].upper().startswith(w) for w in want)]
    if not sel:
        print("no tests matched %s" % argv[1:])
        return 2

    print("regression tests -- %d selected" % len(sel))
    print("=" * 76)
    failed = []
    for fid, what, fn in sel:
        tmp = tempfile.mkdtemp(prefix="regress_%s_" % fid.replace("/", "_"))
        try:
            _r = fn(tmp)
            if not (isinstance(_r, tuple) and len(_r) == 2):
                raise TypeError("test %s returned %r instead of (ok, msg)"
                                % (fid, _r))
            ok, msg = _r
        except Exception as exc:                       # noqa: BLE001
            ok, msg = False, "raised %s: %s" % (type(exc).__name__, exc)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        print("%-5s %-6s %s" % (fid, "ok" if ok else "FAIL", what))
        if not ok:
            print("\n".join("        " + l for l in msg.splitlines()))
            failed.append(fid)
    print("=" * 76)
    if failed:
        print("RESULT: FAIL -- %d of %d: %s" % (len(failed), len(sel),
                                                ", ".join(failed)))
        return 1
    print("RESULT: PASS -- %d/%d" % (len(sel), len(sel)))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
