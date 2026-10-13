#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Self-test the coordinate contract that back-translation depends on.

Why this exists: frame_strand.py promises that amino-acid position i of a written
ORF corresponds to nucleotides [offset + 3i, offset + 3i + 3) of the oriented
sequence.  Everything downstream -- back-translating a protein alignment into a
codon alignment, and therefore any nucleotide tree built from amino acids -- rests
on that mapping.  Off-by-one errors here are silent: the alignment still looks
fine, only every codon is shifted.  So it is asserted, not assumed.

Definition used throughout (1-based, inclusive on both ends):
    aa position i  (i = 1, 2, 3, ...)
        <->  nucleotide positions  offset + 3i - 2  ..  offset + 3i
    where `offset` is 1-based for the first nucleotide of the chosen frame
    (offset 1 = frame p1, 2 = p2, 3 = p3; likewise on the reverse complement).

Run:  python selftest.py
Exit: 0 = all assertions hold, 1 = something is wrong (do not use frame_strand.py).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import frame_strand as fs  # noqa: E402

CODE = 1
TABLE = fs.build_table(CODE)
FAILURES = []


def check(label, got, want):
    ok = got == want
    print("  %-58s %s" % (label, "ok" if ok else "FAIL"))
    if not ok:
        print("      got : %r" % (got,))
        print("      want: %r" % (want,))
        FAILURES.append(label)


print("selftest: frame_strand coordinate contract")
print("-" * 74)

# ---------------------------------------------------------------- 5 bp case ---
# Hand-checkable: A T G | C A  -> frame p1 gives M (ATG) then a partial codon (CA)
# is dropped, so the ORF is "M" and the mapping must be aa1 <-> nt 1..3.
seq = "ATGCAGGTTCAT"
print("reference sequence (1-based):")
for i, ch in enumerate(seq, 1):
    print("   %2d %s" % (i, ch))
print()

for offset0, label in ((0, "p1 (offset 1)"), (1, "p2 (offset 2)"),
                       (2, "p3 (offset 3)")):
    sub = seq[offset0:]
    aa = fs.translate(sub, TABLE)
    print("%s -> aa %r  (from nt %d..%d, %d nt, %d whole codons)"
          % (label, aa, offset0 + 1, len(seq), len(sub), len(sub) // 3))
    # walk the codons and assert the mapping explicitly
    for i in range(1, len(aa) + 1):
        nt_start = offset0 + 3 * i - 2          # 1-based, inclusive
        nt_end = offset0 + 3 * i
        codon = seq[nt_start - 1:nt_end]
        expected = TABLE.get(codon.upper(), "X")
        if expected == "*":
            expected = "X"
        label2 = "offset %d: aa %d <-> nt %d..%d" % (offset0, i, nt_start, nt_end)
        check(label2, aa[i - 1], expected)

print()
print("sub-assertions:")
check("5 bp frame p1 translates to 'M'", fs.translate("ATGCA", TABLE), "M")
check("'ATGCAGGTTCAT' frame p1", fs.translate(seq[0:], TABLE)[0], "M")

# ------------------------------------------------------- stop handling --------
print()
print("stop-codon handling (must become X, must NOT truncate):")
check("TAA -> X", fs.translate("TAA", TABLE), "X")
check("ATG TAA GGG -> 'MXG' (not truncated)", fs.translate("ATGTAAGGG", TABLE), "MXG")
check("ATGTAAGGG is 3 codons worth", len(fs.translate("ATGTAAGGG", TABLE)), 3)
check("stop_count of ATGTAAGGG is 1", fs.stop_count("ATGTAAGGG", TABLE), 1)
check("longest stop-free run of M X G = 1", fs.longest_stop_free("MXG", "X"), 1)

# ------------------------------------------------------- revcomp ------------
print()
print("reverse complement:")
check("revcomp(ATGC) = GCAT", fs.revcomp("ATGC"), "GCAT")
check("revcomp is an involution", fs.revcomp(fs.revcomp("ACGTN")), "ACGTN")

# ------------------------------------------------- genetic code sanity -------
# These are the codons that DISTINGUISH the tables.  They were wrong once
# (table 2 had TGA as stop, table 9 had ATA as Ile), which would have
# mistranslated every mitochondrial COI run from the very first step.
print()
print("genetic code tables -- distinguishing codons (NCBI expectations):")
EXPECT = {
    1: {"TGA": "*", "AGA": "R", "ATA": "I", "AAA": "K"},
    2: {"TGA": "W", "AGA": "*", "ATA": "M", "AAA": "K"},
    3: {"TGA": "W", "AGA": "R", "ATA": "M", "AAA": "K", "CTT": "T"},
    4: {"TGA": "W", "AGA": "R", "ATA": "I", "AAA": "K"},
    5: {"TGA": "W", "AGA": "S", "ATA": "M", "AAA": "K"},
    6: {"TGA": "*", "AGA": "R", "ATA": "I", "AAA": "K", "TAA": "Q", "TAG": "Q"},
    9: {"TGA": "W", "AGA": "S", "ATA": "M", "AAA": "N"},
    11: {"TGA": "*", "AGA": "R", "ATA": "I", "AAA": "K"},
}
for cid in sorted(fs.CODES):
    t = fs.build_table(cid)
    check("code %d: 64-character table" % cid, len(fs.CODES[cid]), 64)
    check("code %d: maps all 64 codons" % cid, len(t), 64)
    for codon, want in EXPECT.get(cid, {}).items():
        check("code %d: %s -> %s" % (cid, codon, want), t[codon], want)

print()
print("the COI case (this is why table 5 exists):")
t5 = fs.build_table(5)
check("table 5 TGA is Trp, not stop", t5["TGA"], "W")
# a COI fragment must translate through TGA without a premature stop
check("COI-like fragment has no premature stop",
      fs.stop_count("TTTTGAGCTGGTGGGACTT", t5), 0)
check("the same fragment on the STANDARD table DOES hit a stop",
      fs.stop_count("TTTTGAGCTGGTGGGACTT", fs.build_table(1)), 1)


# ------------------------------------------------- back-translation ----------
# A REAL round trip, not a revcomp involution.  Back-translation is the step the
# whole amino-acid route rests on (D2b), and its only guarantee is the coordinate
# contract: protein position i comes from nucleotides offset+3i-2 .. offset+3i.
# So: translate -> back-translate -> the codons must be the ORIGINAL nucleotides,
# and translating them again must return the same protein.  The test sequence
# contains TGA, which is a stop on table 1 and Trp on table 5 -- exactly the case
# that has to survive a round trip.
print()
print("back-translation round trip (translate -> back-translate -> translate):")
print("  1-based, inclusive: aa i <-> nt (offset + 3i - 2) .. (offset + 3i)")


def backtranslate(aln, nt, offset0=0):
    """Codons for an ungapped protein alignment, from the oriented nucleotide.

    `aln` is the amino-acid string, `nt` the oriented nucleotide the frame was
    taken from, `offset0` the 0-based frame offset.  Every position whose codon
    exists in `nt` is recovered exactly -- including a codon that translates to a
    stop, because the back-translation tools need the real codon to recognise it.
    Only a gap, or a position with no nucleotides left, becomes '---' / 'NNN'.
    """
    out = []
    for i, aa in enumerate(aln, 1):
        s = offset0 + 3 * i - 3            # 0-based start of the codon
        codon = nt[s:s + 3]
        if aa == "-":
            out.append("---")
        elif len(codon) < 3 or any(c not in "ACGTUacgtu" for c in codon):
            out.append("NNN")
        else:
            out.append(codon)
    return "".join(out)


for code, seq_nt, offset0, tag in (
        (1, "ATGCAGGTTCATTGA", 0, "table 1, ending in a stop"),
        (5, "ATGCAGGTTCATTGA", 0, "table 5, TGA is Trp so no stop at all")):
    t = fs.build_table(code)
    aa = fs.translate(seq_nt[offset0:], t)
    codons = backtranslate(aa, seq_nt, offset0)
    check("%s: back-translated length == 3 x aa length" % tag,
          len(codons), 3 * len(aa))
    check("%s: the recovered codons ARE the input nucleotides" % tag,
          codons, seq_nt[offset0:offset0 + 3 * len(aa)])
    check("%s: re-translating gives the same protein" % tag,
          fs.translate(codons, t), aa)

# ... and the same walk for each frame offset, with a sequence whose length is not
# a multiple of three (the trailing partial codon must be dropped, not fudged)
seq15 = "ATGCAGGTTCATTGA"
for offset0 in (0, 1, 2):
    t = fs.build_table(1)
    aa = fs.translate(seq15[offset0:], t)
    codons = backtranslate(aa, seq15, offset0)
    check("offset %d: recovered codons match the input" % offset0,
          codons, seq15[offset0:offset0 + 3 * len(aa)])
    check("offset %d: protein survives the round trip" % offset0,
          fs.translate(codons, t), aa)

print()
print("the back-translation limits the docs promise (B10):")
check("an internal TGA becomes X on the STANDARD table (PAL2NAL would refuse it)",
      "X" in fs.translate("ATGCCATGAGGG", fs.build_table(1))[:-1], True)
check("the same alignment has no stop on the mitochondrial table",
      "X" not in fs.translate("ATGCCATGAGGG", fs.build_table(5)), True)
check("a codon that does not exist in the input becomes NNN "
      "(trimAl -backtrans cannot map it)",
      backtranslate("MXX", "ATGTG")[3:9], "NNNNNN")
check("a gap becomes '---', never a state",
      backtranslate("M-G", "ATGCAGGGG")[3:6], "---")

print()
print("-" * 74)
if FAILURES:
    print("RESULT: FAIL -- %d assertion(s) failed:" % len(FAILURES))
    for f in FAILURES:
        print("   - %s" % f)
    print()
    print("Do NOT trust frame_strand.py until these pass: a shifted mapping would")
    print("silently corrupt every back-translated codon alignment.")
    sys.exit(1)
print("RESULT: PASS -- coordinate contract, stop handling and genetic codes hold.")
sys.exit(0)
