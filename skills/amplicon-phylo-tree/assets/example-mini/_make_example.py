#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate the synthetic smoke-test dataset (asv_table.tsv + rep_set.fasta).

This data is FABRICATED for smoke-testing the pipeline.  It is deliberately built
to contain the input quirks that cause real trouble:

  * a header line that starts with '#'  ('#OTU ID')  -- a naive "skip # lines"
    parser eats it and every count shifts by one row
  * counts written as '123.0' (float formatting for integer reads)
  * a trailing 'taxonomy' column that must not be summed as a sample
  * sequences wrapped at exactly 60 characters, like MAFFT/most writers do
  * a few all-zero rows so the taxon rule has something to remove
  * five related groups plus a couple of outliers, so a tree has some structure

Re-running this reproduces the same files (fixed seed).

    python _make_example.py
"""
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
TABLE = os.path.join(HERE, "asv_table.tsv")
FASTA = os.path.join(HERE, "rep_set.fasta")

SAMPLES = ["S%02d" % i for i in range(1, 13)]
N_ASV = 20
SEED = 20261012

BASES = "ACGT"
PLANTED_TAX = [
    "k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; g__Rhizobium",
    "k__Bacteria; p__Proteobacteria; c__Betaproteobacteria; g__Burkholderia",
    "k__Bacteria; p__Cyanobacteria; c__Nostocales; g__Nostoc",
    "k__Bacteria; p__Firmicutes; c__Bacilli; g__Paenibacillus",
    "k__Bacteria; p__Actinobacteria; c__Actinomycetia; g__Frankia",
]


def make_sequence(rng, group, length=300):
    """Five related groups plus outliers; group members share a common core."""
    core_rng = random.Random(9000 + group if group < 5 else 9999)
    core = [core_rng.choice(BASES) for _ in range(length)]
    seq = []
    for ch in core:
        if rng.random() < 0.06:          # within-group variation
            seq.append(rng.choice(BASES))
        else:
            seq.append(ch)
    return "".join(seq)


def wrap(s, width=60):
    return [s[i:i + width] for i in range(0, len(s), width)]


def main():
    rng = random.Random(SEED)

    rows = []
    seqs = []
    for i in range(1, N_ASV + 1):
        asv = "ASV%04d" % i
        group = (i - 1) % 5          # five related groups
        if i > 17:
            group = 9                # two outliers, unrelated to the rest
        seq = make_sequence(rng, group)
        seqs.append((asv, seq))

        # abundance: a few ASVs present everywhere, several rare, two all-zero
        if i in (19, 20):
            counts = [0] * len(SAMPLES)              # completely absent
        elif i <= 6:
            counts = [rng.randint(40, 300) for _ in SAMPLES]
        elif i <= 12:
            counts = [0 if rng.random() < 0.4 else rng.randint(10, 90)
                      for _ in SAMPLES]
        else:
            counts = [0 if rng.random() < 0.8 else rng.randint(1, 6)
                      for _ in SAMPLES]
        tax = PLANTED_TAX[group] if group < 5 else "Unclassified"
        rows.append((asv, counts, tax))

    with open(TABLE, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("#OTU ID\t" + "\t".join(SAMPLES) + "\ttaxonomy\n")
        for asv, counts, tax in rows:
            fh.write(asv + "\t" + "\t".join("%.1f" % c for c in counts)
                     + "\t" + tax + "\n")

    with open(FASTA, "w", encoding="utf-8", newline="\n") as fh:
        for asv, seq in seqs:
            fh.write(">%s\n" % asv)
            fh.write("\n".join(wrap(seq)) + "\n")

    total = sum(sum(c) for _a, c, _t in rows)
    print("wrote %s (%d ASVs x %d samples, %d reads)"
          % (TABLE, len(rows), len(SAMPLES), total))
    print("wrote %s (%d sequences, %d bp each, wrapped at 60)"
          % (FASTA, len(seqs), len(seqs[0][1])))
    kept = [a for a, c, _t in rows if sum(c) >= 10]
    print("ASVs passing 'sum>=10': %d of %d" % (len(kept), len(rows)))


if __name__ == "__main__":
    main()
