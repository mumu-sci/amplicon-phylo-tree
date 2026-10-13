#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Suggest a pipeline for this study, and generate the benchmark that confirms it.

Why this exists: the right alignment and tree method depend on BOTH the number of
taxa and the alignment length, and picking by N alone is how people end up
running a method designed for a million sequences on a study that has 400 -- or
running L-INS-i on 2,000 sequences for two days.  Smaller studies should buy
accuracy (iterative refinement + model selection); large ones trade accuracy for
speed.

Two things this script also enforces, because they decide whether the RESULT is
defensible and are easy to get wrong:

  * The ROLE of the tree matters more than N.
      distance skeleton (UniFrac / PD / NRI-NTI / bNTI) -> keep EVERY taxon
      conclusion in itself (novel clade, within-clade relationships) -> narrow
      first, then spend compute on accuracy
    "High accuracy first" is the wrong default: it optimises the wrong objective
    for a skeleton, where the root, region homogeneity, reference DB version and
    tip consistency matter far more than alignment refinement.

  * Resource estimates are planning heuristics, NOT promises.  The output
    REQUIRES you to run the printed benchmark before committing.

Provenance for every threshold and every measured anchor: references/scaling.md.
"""
import argparse
import sys

# --------------------------------------------------------------------------
# alignment choice by (N, L) -- mirrors references/scaling.md
# --------------------------------------------------------------------------

def alignment_plan(n, L):
    """Return (method, command_template, why, cost_note).

    TWO-DIMENSIONAL on purpose: the same N can call for a different algorithm when
    L changes magnitude (1,900 x 4,000 bp already exhausted memory on the dev
    machine).  The boundaries here are the ones written in references/scaling.md;
    scripts/check_consistency.py asserts the two stay in step, so changing one
    without the other fails a test instead of drifting silently.
    """
    if n <= 200:
        if L and L > 2000:
            return ("L-INS-i (explicit: --localpair --maxiterate 1000)",
                    "mafft --localpair --maxiterate 1000 --thread {cores} in.fasta "
                    "> aln.fasta",
                    "few enough sequences to afford the high-accuracy alignment; at "
                    "this sequence length --auto's own L-INS-i choice is fine too, "
                    "but it is worth saying explicitly because it is slow",
                    "long sequences: estimate the wall clock before committing")
        return ("--auto (MAFFT selects L-INS-i itself at this size)",
                "mafft --auto --quiet --thread {cores} in.fasta > aln.fasta",
                "at <=200 sequences MAFFT's own --auto picks the high-accuracy "
                "L-INS-i, so there is nothing to decide",
                "seconds to minutes")
    if n <= 3000:
        if L and L > 2000:
            return ("E-INS-i  (long sequences usually contain unalignable stretches)",
                    "mafft --genafpair --maxiterate 1000 --thread {cores} in.fasta "
                    "> aln.fasta",
                    "L is in the thousands of columns, so allow multiple conserved "
                    "domains and long unalignable regions; E-INS-i is also the first "
                    "thing to try when an alignment comes out fragmented",
                    "minutes to hours -- ESTIMATE before committing")
        return ("L-INS-i   (one alignable domain; typical for a single amplicon)",
                "mafft --localpair --maxiterate 1000 --thread {cores} in.fasta "
                "> aln.fasta",
                "few enough sequences that high accuracy is affordable; if the "
                "region has long unalignable stretches use E-INS-i "
                "(--genafpair --maxiterate 1000) instead",
                "minutes, or longer with long sequences -- benchmark")
    if n <= 100000:
        if L and L > 2000:
            return ("E-INS-i  (long sequences usually contain unalignable stretches)",
                    "mafft --genafpair --maxiterate 1000 --thread {cores} in.fasta "
                    "> aln.fasta",
                    "L is in the thousands of columns, so the single-domain "
                    "assumption behind FFT-NS-i stops holding; E-INS-i allows "
                    "multiple conserved regions and is also the first thing to try "
                    "when an alignment comes out fragmented",
                    "it is the slowest of the three at this size -- run the benchmark "
                    "FIRST and record the result; parttree only if the benchmark "
                    "proves it necessary")
        return ("FFT-NS-i  (or plain --auto, which selects a progressive strategy)",
                "mafft --maxiterate 1000 --thread {cores} in.fasta > aln.fasta",
                "L-INS-i's practical limit is a few hundred sequences, not a few "
                "thousand; iterative refinement is the accurate-but-affordable "
                "option here.  Only escalate if the benchmark says you must.",
                "tens of minutes to a few hours -- benchmark")
    return ("--parttree --retree 2",
            "mafft --parttree --retree 2 --quiet --thread {cores} in.fasta "
            "> aln.fasta",
            "at this N a full all-pairs guide tree is the bottleneck (that is what "
            "makes --auto slow), and this is the scale parttree is designed for",
            "parttree trades alignment quality for speed -- SAY SO in the methods")


def tree_plan(n, L):
    """Return (tool, command_template, support, note).

    Mirrors the two-dimensional tree table in references/scaling.md.  The residue
    clause (N x L <= 3e7) is the only place L changes the tool, and it must never be
    allowed to override a huge taxon count: exact ML with model selection plus 1000
    UFBoot replicates on 190,000 tips is not a "big job", it is a wrong
    recommendation.
    """
    MAX_RESIDUES_IQTREE = 3e7                      # N x L, from scaling.md
    total_residues = (n * L) if (L and n) else None
    if n <= 100000 and (n <= 20000 or (total_residues is not None
                                       and total_residues <= MAX_RESIDUES_IQTREE)):
        return ("IQ-TREE (exact ML, model selection)",
                "iqtree2 -s aln.fasta -m MFP -B 1000 -alrt 1000 -T {cores} "
                "--prefix tree -pre {runid}",
                "UFBoot 1000 + SH-aLRT 1000 (report both)",
                "exact ML is affordable at this size; benchmark to confirm")
    return ("FastTree (approximate ML)",
            "fasttree -nt -gtr -gamma aln.fasta > tree.nwk 2> tree.log",
            "SH-like 1000 -- and state in the methods that this is APPROXIMATE ML",
            "at this scale exact ML with full UFBoot is not a reasonable default; "
            "if you want IQ-TREE, try -fast and be explicit about it")


def scale_warning(n, ram_gb):
    """Loud, non-optional warnings for sizes where people get burned."""
    out = []
    if n > 50000:
        out.append("!! n = %d: run the benchmark BEFORE committing.  The nearest "
                   "measured anchor (~70,000 taxa) exceeded 13.9 GB and did not "
                   "finish on a 13.9 GB machine.  Memory grows faster than linearly "
                   "in taxa, so extrapolate UPWARD." % n)
    if n > 200000:
        out.append("!! n = %d is in the regime where MAFFT --auto builds a full "
                   "all-pairs guide tree (O(n^2)); use --parttree and say so." % n)
    if ram_gb and n >= 50000 and ram_gb < 32:
        out.append("!! %d taxa with only %.0f GB RAM is very likely to fail or "
                   "swap-thrash.  Consider a bigger machine, or reduce the taxon set "
                   "as a STUDY decision (and report the costs)." % (n, ram_gb))
    return out


ANCHORS = [
    (3874, 455, "MAFFT --auto 10 s; FastTree 66 s", "well under 0.3 GB"),
    (10687, 4258, "FastTree 1,543.7 s (25.7 min)", "~0.96 GB"),
    (30258, 1560, "FastTree 3,180.9 s (53 min)", "2-4 GB"),
    (69668, 1560, "did NOT finish on a 13.9 GB laptop", "> 13.9 GB (confounded run, "
     "see scaling.md)"),
    (190032, 1560, "MAFFT --parttree --retree 2 = 158 s; FastTree 2,609.8 s on a "
     "16-core/64 GB host", "GBs"),
]


def resource_note(n, L, ram_gb):
    lines = []
    if L:
        cells = n * L
        lo = cells * 50 / 1e9
        hi = cells * 200 / 1e9
        lines.append("planning heuristic (NOT a promise): alignment peak ~%.1f-%.1f GB "
                     "for %d x %d = %.1fM cells (alpha 50-200 B/cell)"
                     % (lo, hi, n, int(L), cells / 1e6))
        lines.append("FastTree measured ~21 B/cell at 10,687 x 4,258, but the "
                     "coefficient RISES with N (memory grows faster than linearly "
                     "in taxa) -- extrapolate upward")
    else:
        lines.append("pass --median-len for a memory heuristic")
    if ram_gb:
        if n >= 50000 and ram_gb < 32:
            lines.append("!! N=%d with %.0f GB RAM: the nearest measured anchor "
                         "(~70k taxa) exceeded 14 GB.  Benchmark first; expect to "
                         "need a bigger machine or to align in parts."
                         % (n, ram_gb))
        elif n >= 20000 and ram_gb < 8:
            lines.append("!! N=%d with %.0f GB RAM is tight -- benchmark first."
                         % (n, ram_gb))
    return lines


def main():
    ap = argparse.ArgumentParser(
        description="Suggest a pipeline for this study, and print the benchmark "
                    "that must confirm it.")
    ap.add_argument("--n-seq", type=int, required=True)
    ap.add_argument("--median-len", type=float, default=None)
    ap.add_argument("--marker-class", required=True,
                    choices=["rrna", "its", "functional-nuclear", "organelle",
                             "other"])
    ap.add_argument("--tree-role", required=True,
                    choices=["skeleton", "conclusion", "placement", "figure",
                             "root-only"],
                    help="skeleton = UniFrac/PD/NRI-NTI/bNTI (needs ALL taxa and a "
                         "rooted tree); conclusion = you will publish the tree "
                         "itself; placement = reference-based; root-only = you only "
                         "need to root/verify an existing tree")
    ap.add_argument("--seqtype", choices=["nt", "aa", "codon"], default=None)
    ap.add_argument("--p-aa-identity", type=float, default=None,
                    help="median within-domain amino-acid identity; drives D2b")
    ap.add_argument("--genetic-code", type=int, default=None)
    ap.add_argument("--cores", type=int, default=8)
    ap.add_argument("--ram-gb", type=float, default=None)
    ap.add_argument("--has-outgroup", action="store_true")
    args = ap.parse_args()

    n, L = args.n_seq, args.median_len
    cores = args.cores

    print("=" * 78)
    print("STRATEGY SUGGESTION -- a starting point, NOT a conclusion")
    print("=" * 78)
    print("  taxa            : %d" % n)
    print("  median length   : %s" % ("%g bp" % L if L else "(not given)"))
    print("  marker class    : %s" % args.marker_class)
    print("  tree role       : %s" % args.tree_role)
    if L:
        print("  total residues  : %.2f M" % (n * L / 1e6))
    print()

    # ---------- what the ROLE decides (this comes before any tool choice) -----
    print("  WHAT THE ROLE DECIDES (read this before choosing tools)")
    role_notes = {
        "skeleton": [
            "keep EVERY taxon in the analysis matrix -- the distance matrix is "
            "defined over all of them, and dropping taxa changes the definition "
            "of the metric",
            "you do NOT need L-INS-i for this; spend the effort on: the ROOT, "
            "region homogeneity, reference DB version, and tip consistency",
            "rooted tree required; prefer an outgroup or reference rooting",
        ],
        "conclusion": [
            "narrow the taxon set FIRST, and write the rule into the methods "
            "before you see the tree (not tuned afterwards to look good)",
            "report the reads % and richness % retained, and do ONE sensitivity "
            "run with a looser rule",
            "now high-accuracy alignment and exact ML are worth paying for",
        ],
        "placement": [
            "reference sequences are required; consider a reference backbone "
            "alignment rather than a free de novo alignment",
            "root with an outgroup, not midpoint",
        ],
        "figure": [
            "state the support-value type explicitly, and never compare support "
            "values across programs",
        ],
        "root-only": [
            "skip alignment and tree inference; go straight to root_tree.R",
        ],
    }
    for r in role_notes[args.tree_role]:
        print("    * %s" % r)
    print()
    print("    !! never: build a full tree, then subset it for metric "
          "calculation -- that silently moves branch lengths and the root")
    print()

    # ---------- prerequisites ----------
    pre = []
    if args.marker_class in ("functional-nuclear", "organelle"):
        code = args.genetic_code
        if args.marker_class == "organelle" and not code:
            pre.append("ORGANELLE marker: you MUST set the genetic code.  For COI use "
                       "5 (invertebrate mitochondrial) or 2 (vertebrate "
                       "mitochondrial).  The standard table treats TGA as a stop "
                       "codon, so the WRONG reading frame gets chosen.")
        if args.p_aa_identity is None:
            pre.append("D2b is not decided: pass --p-aa-identity.  "
                       "<85% -> amino acids AND back-translate to a codon alignment; "
                       "85-90% -> build both and report topological agreement; "
                       ">90% -> a nucleotide tree is acceptable, say the third "
                       "position is unsaturated.")
        else:
            p = args.p_aa_identity
            if p < 85:
                pre.append("D2b: p_aa=%.1f%% < 85%% -> align at AMINO-ACID level and "
                           "BACK-TRANSLATE to codons (PAL2NAL / trimAL -backtrans / "
                           "MACSE) so the tree's branch lengths stay consistent."
                           % p)
            elif p <= 90:
                pre.append("D2b: p_aa=%.1f%% is in the grey zone 85-90%% -> build "
                           "BOTH the nucleotide and the amino-acid tree, report "
                           "topological agreement (RF / conflicting branches), then "
                           "pick one as the main tree." % p)
            else:
                pre.append("D2b: p_aa=%.1f%% > 90%% -> a nucleotide tree is "
                           "acceptable; state that the third codon position is "
                           "unsaturated." % p)
        pre.append("normalise reading frame + strand per sequence FIRST "
                   "(scripts/frame_strand.py), then run scripts/selftest.py to "
                   "confirm the coordinate contract.")
    if args.marker_class == "its":
        pre.append("ITS: do NOT build an ITS1/ITS2 tree for community metrics.  "
                   "Preferred order: (1) no tree (UNITE dynamic-clustering match), "
                   "(2) fragment insertion (SEPP/EPA-ng) into a UNITE/SILVA "
                   "reference tree, (3) narrow to genus-level + secondary-structure-"
                   "constrained alignment, supplementary only.")
    if args.marker_class == "rrna":
        pre.append("keep one amplicon region and one primer pair per tree; pair the "
                   "region with its reference DB (SILVA vs PR2 are NOT "
                   "interchangeable).")
    if args.tree_role == "skeleton":
        pre.append("make sure the biological replicate unit matches the tree's tip "
                   "level; sequencing replicates of one treatment are "
                   "pseudo-replication for a phylogenetic metric.")
    if pre:
        print("  DO THESE FIRST")
        for p in pre:
            print("    * %s" % p)
        print()

    if args.tree_role == "root-only":
        print("  next: Rscript scripts/root_tree.R --in tree.nwk --rooted ... "
              "--sensitivity-candidates 3")
        return 0

    # ---------- ITS in a community/placement role: route, do not build ----------
    # Printing "never build an ITS tree" and then a ready-to-run MAFFT + IQ-TREE
    # command line is not a warning: an agent executes what it is shown.  So this
    # branch prints NO de novo commands at all.
    if args.marker_class == "its" and args.tree_role in ("skeleton", "placement"):
        print("  ROUTE FOR ITS IN A %s ROLE -- do NOT build a de novo ITS tree"
              % args.tree_role.upper())
        print("  " + "-" * 68)
        print("  This branch deliberately prints no alignment and no tree command.")
        print("  Pick one of these instead, in order of preference:")
        print()
        print("  (1) NO TREE  -- match against UNITE and use the taxonomy, not a tree")
        print("        vsearch --usearch_global asvs.fasta \\")
        print("                --db unite.fasta --id 0.97 --uc hits.uc \\")
        print("                --strand plus --maxaccepts 1 --maxrejects 32")
        print("      (--id 0.97 for species-level, 0.99 for a stricter match;")
        print("       record which, and the UNITE version, in params.md D11)")
        print("      Split ITS1/ITS2 first if your primers are not ITS-wide:")
        print("        ITSx -i asvs.fasta --preserve T --cpu {cores} -o itsx")
        print()
        print("  (2) FRAGMENT INSERTION -- place your sequences into a reference tree")
        print("      This keeps published branch lengths; declare in the methods that")
        print("      branch lengths come from INSERTION and are not comparable with a")
        print("      de novo tree.")
        print("        # EPA-ng, given a UNITE reference alignment + tree:")
        print("        epa-ng -t ref.tree -s ref_aln.fasta -q asvs.fasta \\")
        print("               --model GTR+G --outdir epa_out")
        print("        # or SEPP (fragment insertion into a taxonomically-constrained")
        print("        # tree); note whichever you use in run_manifest.txt")
        print()
        print("  (3) A NARROWED, ALIGNMENT-CONSTRAINED tree for SUPPLEMENTARY use only")
        print("      Cross-genus relationships are not defensible from ITS.  If you")
        print("      must show a tree, narrow to one genus, use a secondary-structure")
        print("      aware alignment, and label it supplementary:")
        print("        - ITS2 database / 4SALE / RNAsalsa for structure-aware alignment")
        print("        - rrnDB or the ITS2 model for structure; Gblocks for trimming")
        print("        - then re-run this script with --tree-role conclusion")
        print()
        print("  WHY: ITS length is heterogeneous and its copies do not fully")
        print("  co-evolve, so a multiple-sequence alignment above genus level is not")
        print("  reproducible and a tree from it does not mean what a reader assumes.")
        print("  Using one for UniFrac/betaNTI compares quantities that are not")
        print("  comparable across samples.")
        return 0

    # ---------- alignment ----------
    method, cmd, why, cost = alignment_plan(n, L)
    print("  ALIGNMENT")
    print("    method : %s" % method)
    print("    why    : %s" % why)
    print("    cost   : %s" % cost)
    print("    command:")
    for line in cmd.format(cores=cores).splitlines():
        print("      %s" % line)
    print("    then ALWAYS check it:")
    print("      python scripts/aln_qc.py --aln aln.fasta --fasta-raw in.fasta")
    print()

    # ---------- tree ----------
    tool, tcmd, support, tnote = tree_plan(n, L)
    print("  TREE INFERENCE")
    print("    tool    : %s" % tool)
    print("    support : %s" % support)
    print("    note    : %s" % tnote)
    print("    command:")
    for line in tcmd.format(cores=cores, runid="run$(date +%Y%m%d%H%M)").splitlines():
        print("      %s" % line)
    if args.seqtype == "aa" or (args.p_aa_identity is not None
                                and args.p_aa_identity < 85):
        print("    (amino-acid input: use 'fasttree -lg -gamma orfs.faa' / "
              "'iqtree2 -s aln.faa -m MFP')")
    print()

    # ---------- rooting ----------
    print("  ROOTING")
    if args.has_outgroup:
        print("    outgroup available -> USE IT (a hypothesis, not a convention):")
        print("      Rscript scripts/root_tree.R --in tree.nwk --rooted tree.rooted.nwk "
              "--unrooted tree.unrooted.nwk --report rooting.txt --label <study> "
              "--outgroup outgroup.txt")
    else:
        print("    no outgroup.  If the role is 'skeleton' or 'placement', a bare "
              "midpoint root is NOT defensible -- you must either find a reference "
              "sister lineage or report a rooting-sensitivity analysis:")
        print("      Rscript scripts/root_tree.R --in tree.nwk --rooted tree.rooted.nwk "
              "--unrooted tree.unrooted.nwk --report rooting.txt --label <study> "
              "--sensitivity-candidates 3")
    print("    keep the unrooted tree too: support values live there.")
    print()

    # ---------- resources ----------
    warn = scale_warning(n, args.ram_gb)
    if warn:
        print("  WARNINGS (read before you start)")
        for w in warn:
            print("    %s" % w)
        print()
    print("  RESOURCES")
    for r in resource_note(n, L, args.ram_gb):
        print("    %s" % r)
    print("    measured anchors (taxa x columns -> step time, peak memory):")
    for a_n, a_L, t, mem in ANCHORS:
        print("      %7d x %-5d -> %-58s %s" % (a_n, a_L, t, mem))
    print()

    # ---------- benchmark ----------
    print("=" * 78)
    print("CONFIRM WITH A BENCHMARK (do not trust the table for wall-clock)")
    print("=" * 78)
    print("""  1. Take a reproducible subsample of about 2,000 sequences and WRITE DOWN
     the rule (e.g. 'seqkit sample -n 2000 -s 11', or 'first 2000 records in
     file order').
  2. Run the REAL alignment command and the REAL tree command on it, and record
     PEAK MEMORY, not just time:
       Linux/macOS : /usr/bin/time -v <command>   (or -l on macOS)
       Windows     : run inside WSL so /usr/bin/time is available, or sample the
                     process RSS periodically.
  3. Time grows at least linearly in taxa and roughly linearly in columns;
     memory grows FASTER than linearly in taxa.  Round UP.
  4. If it does not fit: the first lever is the TAXON RULE (a study decision),
     the second is a reference backbone or blocks, and only the third is a
     bigger machine.""")
    return 0


if __name__ == "__main__":
    sys.exit(main())
