#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Partitioned alignment: align in blocks, then concatenate with gap padding.

Why this exists: when a taxon set contains sequences that are not mutually
homologous (mixed amplicon regions, off-target products, extremely divergent
lineages, extreme length heterogeneity), a single global alignment cannot relate
them.  The aligner then lays the groups side by side and the matrix explodes --
measured on a real dataset: 3,933 columns and 92% gaps globally, versus 455
columns and 29% gaps once the comparable part was aligned on its own.

Block partitioning keeps EVERY taxon: each block is aligned internally, then the
block alignments are concatenated with gap padding, so a sequence is a real
character only inside its own block and all-gap elsewhere.  Relationships inside a
block are meaningful; between blocks they are not, and that boundary must be
stated in the methods.

This is a remedy for a diagnosed problem, not the default path.  Run it after
aln_qc.py reports fragmentation, or when you know the set contains distinct groups
(e.g. you have a cluster/species assignment).

Two ways to define the blocks, in order of preference:

  --assign MAP   use an existing grouping you trust (cluster map, taxonomy, region
                 label, CD-HIT families).  Preferred: the grouping is a study
                 decision you can defend.
  de novo        single linkage on shared canonical k-mers.  This is a diagnostic
                 aid, not a classification.  Sweep k, look at the table, and pick
                 on evidence.

Subcommands
    plan    report a k sweep (or read --assign), choose a partition, write blocks
    split   re-partition ONE block that still failed aln_qc (recursive refinement)
    concat  gap-pad concatenate the aligned blocks into one supermatrix
"""
import argparse
import collections
import glob as _glob
import os
import sys

COMP = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")


def read_fasta(path):
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


def write_fasta(records, path):
    d = os.path.dirname(os.path.abspath(path))
    if d:
        os.makedirs(d, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for name, seq in records:
            fh.write(">%s\n%s\n" % (name, seq))


def looks_like_protein(seqs, sample=200):
    seen = set()
    for _, s in seqs[:sample]:
        seen.update(s.upper())
    nt_chars = set("ACGTUNRYSWKMBDHV-.")
    return not (seen <= nt_chars)


def kmers(seq, k, protein):
    out = set()
    s = seq.upper()
    for j in range(len(s) - k + 1):
        a = s[j:j + k]
        if "-" in a or "." in a:
            continue
        if protein:
            out.add(a)
        else:
            b = a.translate(COMP)[::-1]
            out.add(a if a < b else b)
    return out


def components(seqs, k, protein):
    """Single-linkage components over shared k-mers.  Returns list of index lists."""
    parent = list(range(len(seqs)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    index = collections.defaultdict(list)
    for i, (_, s) in enumerate(seqs):
        for kk in kmers(s, k, protein):
            index[kk].append(i)
    for kk, lst in index.items():
        first = lst[0]
        for j in lst[1:]:
            union(first, j)

    groups = collections.defaultdict(list)
    for i in range(len(seqs)):
        groups[find(i)].append(i)
    return sorted(groups.values(), key=len, reverse=True)


def read_assign(path, idcol=0, grpcol=1, delim=None):
    m = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split(delim) if delim else line.replace(",", "\t").split("\t")
            if len(parts) <= max(idcol, grpcol):
                continue
            m[parts[idcol].strip()] = parts[grpcol].strip()
    return m


# -------------------------------------------------------------------- plan ---

def cmd_plan(args):
    seqs = list(read_fasta(args.fasta))
    if not seqs:
        print("ERROR: no sequences in %s" % args.fasta)
        return 2
    protein = looks_like_protein(seqs)
    print("input            : %s (%d sequences)" % (args.fasta, len(seqs)))
    print("alphabet         : %s" % ("protein" if protein else "nucleotide"))
    if protein and args.assign is None:
        print("  note: protein input -- k-mers are used without reverse-complement "
              "canonicalisation.")
    if len(seqs) > 50000 and args.assign is None:
        print("  warning: de-novo partitioning is memory-hungry above ~50,000 "
              "sequences; prefer --assign.")

    os.makedirs(args.outdir, exist_ok=True)
    widths = []
    os.makedirs(args.outdir, exist_ok=True)

    # ---------- external assignment ----------
    if args.assign:
        amap = read_assign(args.assign, args.assign_id_col, args.assign_group_col,
                           args.assign_delim)
        groups = collections.defaultdict(list)
        unknown = []
        for name, seq in seqs:
            g = amap.get(name)
            if g is None:
                unknown.append(name)
            else:
                groups[g].append((name, seq))
        print()
        print("assignment       : %s" % args.assign)
        print("groups found     : %d" % len(groups))
        print("sequences mapped : %d / %d" % (len(seqs) - len(unknown), len(seqs)))
        if unknown:
            print("  !! %d sequence(s) have no group.  They cannot be placed, and "
                  "silently dropping them would change the taxon set." % len(unknown))
            print("     e.g. %s" % ", ".join(unknown[:8]))
            print("     Fix the assignment, or pass a fallback group explicitly.")
            return 1
        order = sorted(groups.items(), key=lambda kv: -len(kv[1]))
        for i, (g, recs) in enumerate(order, 1):
            p = os.path.join(args.outdir, "block_%02d_%s.fasta"
                             % (i, str(g).replace(os.sep, "_")[:40]))
            write_fasta(recs, p)
            print("  block %02d  %-24s %6d seqs -> %s" % (i, str(g)[:24], len(recs), p))
        print()
        print("manifest written : %s"
              % os.path.join(args.outdir, "blocks.tsv"))
        with open(os.path.join(args.outdir, "blocks.tsv"), "w",
                  encoding="utf-8", newline="\n") as fh:
            fh.write("block\tgroup\tfasta\tn_seqs\n")
            for i, (g, recs) in enumerate(order, 1):
                fh.write("block_%02d\t%s\t%s\t%d\n"
                         % (i, g, os.path.join(args.outdir, "block_%02d_%s.fasta"
                                               % (i, str(g).replace(os.sep, "_")[:40])),
                            len(recs)))
        print()
        print("Next: align EACH block separately, then:")
        print("  python scripts/aln_qc.py --aln <block>.aligned.fasta --fasta-raw "
              "<block>.fasta --label <block>")
        print("  python scripts/blocks.py concat --aln-glob '%s/*.aligned.fasta' "
              "--ids kept.ids --out supermatrix.fasta" % args.outdir)
        return 0

    # ---------- de novo: sweep k ----------
    ks = [int(x) for x in args.k_sweep.split(",")] if args.k_sweep else [args.k]
    print()
    print("k sweep (single linkage on shared canonical k-mers)")
    print("  %-6s %-12s %-12s %-14s %s" % ("k", "components", "singletons",
                                            "top-blocks", "largest handful"))
    sweep = {}
    for k in ks:
        comps = components(seqs, k, protein)
        sizes = [len(c) for c in comps]
        singles = sum(1 for s in sizes if s == 1)
        top = sum(sizes[:args.max_blocks])
        sweep[k] = (comps, sizes, singles, top)
        print("  %-6d %-12d %-12d %-14s %s" % (
            k, len(comps), singles, "%d (%.1f%%)" % (top, 100.0 * top / len(seqs)),
            sizes[:6]))
    print()
    print("  how to read this: a SMALL k merges everything into one block (loose); "
          "a LARGE k splits finely (strict).")
    print("  you want the LARGEST k whose block count is still manageable, because "
          "each block is then as internally similar as possible.")
    print("  !! a large top block does not prove it is one alignable family -- "
          "k-mer connectivity chains through divergence.  Align it and run "
          "aln_qc.py to find out.")

    chosen = args.k
    if args.k_sweep and args.k not in sweep:
        ok = [k for k in sorted(sweep, reverse=True)
              if len(sweep[k][0]) <= args.max_blocks]
        chosen = ok[0] if ok else min(sweep)
        print()
        print("  chosen k = %d  (largest sweep k with <= %d blocks)"
              % (chosen, args.max_blocks))

    comps, sizes, singles, top = sweep.get(chosen) or components(seqs, chosen, protein)
    print()
    print("partition at k = %d" % chosen)
    if len(comps) > args.max_blocks:
        print("  !! %d blocks exceeds --max-blocks %d.  Options: lower k, use "
              "--assign with a grouping you trust, or keep the small tail in one "
              "'misc' block and state that its internal relationships are "
              "unresolved." % (len(comps), args.max_blocks))

    rows = []
    for i, idxs in enumerate(comps, 1):
        recs = [seqs[j] for j in idxs]
        gname = "block_%02d" % i
        p = os.path.join(args.outdir, "%s.fasta" % gname)
        write_fasta(recs, p)
        rows.append((gname, p, len(recs)))
    with open(os.path.join(args.outdir, "blocks.tsv"), "w",
              encoding="utf-8", newline="\n") as fh:
        fh.write("block\tfasta\tn_seqs\n")
        for g, p, n in rows:
            fh.write("%s\t%s\t%d\n" % (g, p, n))
    shown = rows if len(rows) <= 15 else rows[:15]
    for g, p, n in shown:
        print("  %-10s %6d seqs -> %s" % (g, n, p))
    if len(rows) > len(shown):
        print("  ... and %d more (see blocks.tsv)" % (len(rows) - len(shown)))
    print()
    print("Next: align each block, QC each block, then concat:")
    print("  for f in %s/block_*.fasta; do mafft --auto --quiet $f > $f.aligned; done"
          % args.outdir)
    print("  python scripts/aln_qc.py --aln <block>.fasta.aligned --fasta-raw "
          "<block>.fasta --label <block>")
    print("  # any block that fails:  python scripts/blocks.py split --fasta "
          "<block>.fasta --outdir <block>_split")
    print("  python scripts/blocks.py concat --aln-glob '%s/block_*.fasta.aligned' "
          "--ids kept.ids --out supermatrix.fasta" % args.outdir)
    return 0


# ------------------------------------------------------------------- split ---

def cmd_split(args):
    seqs = list(read_fasta(args.fasta))
    protein = looks_like_protein(seqs)
    comps = components(seqs, args.k, protein)
    os.makedirs(args.outdir, exist_ok=True)
    print("splitting %s (%d sequences) at k=%d -> %d sub-block(s)"
          % (args.fasta, len(seqs), args.k, len(comps)))
    for i, idxs in enumerate(comps, 1):
        recs = [seqs[j] for j in idxs]
        p = os.path.join(args.outdir, os.path.basename(args.fasta).rsplit(".", 1)[0]
                         + "_sub%02d.fasta" % i)
        write_fasta(recs, p)
        print("  sub%02d  %6d seqs -> %s" % (i, len(recs), p))
    print()
    print("Align each sub-block and QC it.  Repeat until aln_qc passes.")
    return 0


# ------------------------------------------------------------------ concat ---

def cmd_concat(args):
    files = sorted(_glob.glob(args.aln_glob))
    if not files:
        print("ERROR: no files matched %s" % args.aln_glob)
        return 2
    if args.ids:
        ids = [l.strip() for l in open(args.ids, encoding="utf-8") if l.strip()]
    else:
        ids = None

    blocks = []
    for path in files:
        recs = list(read_fasta(path))
        widths = set(len(s) for _, s in recs)
        if len(widths) != 1:
            print("ERROR: %s is not rectangular (widths %s).  Each block must be "
                  "an alignment." % (path, sorted(widths)[:5]))
            return 1
        blocks.append((os.path.basename(path), recs, widths.pop()))

    # membership: which block does each sequence belong to?
    owner = {}
    dupes = []
    for bname, recs, _w in blocks:
        for name, _s in recs:
            if name in owner:
                dupes.append(name)
            owner[name] = bname
    if dupes:
        print("ERROR: %d sequence(s) appear in more than one block, e.g. %s"
              % (len(dupes), ", ".join(sorted(set(dupes))[:8])))
        print("  Each sequence must belong to exactly one block, otherwise its "
              "characters are counted twice.")
        return 1

    # The MASTER id list is authoritative.  Aligning blocks separately can drop
    # or reorder sequences (aligners do this), so:
    #   sequences missing from a block -> all-gap across that block (they belong
    #   to a different block anyway, so this is the correct padding)
    #   sequences in no block at all   -> hard error, they would vanish silently
    if ids is None:
        ids = list(owner.keys())
        print("no --ids given; using the %d sequences found in the blocks" % len(ids))
        print("  NOTE: without --ids the output is defined by what the aligners "
              "happened to emit.  Pass the taxon list for a checked result.")
    missing = [i for i in ids if i not in owner]
    if missing:
        print("ERROR: %d expected sequence(s) are in NO block, e.g. %s"
              % (len(missing), ", ".join(missing[:8])))
        print("  They would silently vanish from the tree.  Fix the partitioning.")
        return 1
    unknown = [i for i in owner if i not in set(ids)]
    if unknown:
        print("  note: %d sequence(s) are in a block but not in --ids; they are "
              "kept (the tree will have more tips than --ids).  e.g. %s"
              % (len(unknown), ", ".join(sorted(unknown)[:6])))

    total_w = sum(w for _b, _r, w in blocks)
    filler_char = args.gap_char
    if args.gap_char not in "-.?":
        print("ERROR: --gap-char %r is not a recognised missing-data symbol.  Use "
              "'-' (the default, and what every tool reads), '.' or '?' -- an "
              "unrecognised filler is read as a state and silently corrupts the "
              "matrix." % args.gap_char)
        return 2
    print("blocks           : %d" % len(blocks))
    print("output sequences : %d" % len(ids))
    print("output columns   : %d" % total_w)
    print("gap filler       : %r" % filler_char)
    for bname, recs, w in blocks:
        print("   %-34s %6d seqs x %5d cols  (%.1f%% of width)"
              % (bname, len(recs), w, 100.0 * w / total_w))

    # per-block lookup, plus the all-gap filler
    lookups = []
    for bname, recs, w in blocks:
        lookups.append((dict(recs), filler_char * w))

    d = os.path.dirname(os.path.abspath(args.out))
    if d:
        os.makedirs(d, exist_ok=True)
    coverage_rows = []
    with open(args.out, "w", encoding="utf-8", newline="\n") as fh:
        for sid in ids:
            parts = []
            present = 0
            for (lk, filler) in lookups:
                if sid in lk:
                    parts.append(lk[sid])
                    present += 1
                else:
                    parts.append(filler)
            fh.write(">%s\n%s\n" % (sid, "".join(parts)))
            coverage_rows.append((sid, present, len(lookups)))
    print()
    print("supermatrix      : %s (%d bytes)" % (args.out, os.path.getsize(args.out)))

    # block coverage per sequence: how much of the matrix actually has data for it
    cov_path = args.out + ".coverage.tsv"
    with open(cov_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("seq_id\tblocks_with_data\tblocks_total\tcoverage\n")
        for sid, present, total in sorted(coverage_rows, key=lambda r: r[1]):
            fh.write("%s\t%d\t%d\t%.3f\n" % (sid, present, total, present / total))
    low = [r for r in coverage_rows if r[1] <= 1 and len(lookups) > 1]
    print("coverage report  : %s" % cov_path)
    if low:
        print("  !! %d sequence(s) have data in only ONE block.  They are kept "
              "(nothing is deleted here), but their position in the tree rests on a "
              "single block -- consider listing them as low-information in the "
              "methods, or aligning them separately." % len(low))
        print("     e.g. %s" % ", ".join(r[0] for r in low[:8]))
    print()
    print("Now:")
    print("  python scripts/aln_qc.py --aln %s --label supermatrix" % args.out)
    print("  #  then build the tree from it")
    print()
    print("STATE THIS IN THE METHODS: sequences that belong to different blocks "
          "share no columns, so branch lengths BETWEEN blocks are not biologically "
          "interpretable.  Relationships within a block are.  For protein-coding "
          "markers, concatenating blocks also breaks codon continuity across block "
          "boundaries -- prefer a reference backbone alignment if you can get one.")
    return 0


def main():
    ap = argparse.ArgumentParser(
        description="Partitioned alignment: align in blocks, concatenate with "
                    "gap padding.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("plan", help="choose a partition and write per-block FASTAs")
    p.add_argument("--fasta", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--k", type=int, default=15,
                   help="k-mer length for de-novo families (default 15)")
    p.add_argument("--k-sweep", default="21,18,15,13,11",
                   help="comma list to sweep (default 21,18,15,13,11); the largest "
                        "k with <= --max-blocks blocks is chosen")
    p.add_argument("--max-blocks", type=int, default=8)
    p.add_argument("--assign", default=None,
                   help="two-column id<TAB>group file; preferred over de-novo")
    p.add_argument("--assign-id-col", type=int, default=0)
    p.add_argument("--assign-group-col", type=int, default=1)
    p.add_argument("--assign-delim", default=None,
                   help="column delimiter for --assign (default: tab, or comma)")
    p.set_defaults(func=cmd_plan)

    s = sub.add_parser("split", help="re-partition one block that still fragments")
    s.add_argument("--fasta", required=True)
    s.add_argument("--outdir", required=True)
    s.add_argument("--k", type=int, default=15)
    s.set_defaults(func=cmd_split)

    c = sub.add_parser("concat", help="gap-pad concatenate aligned blocks")
    c.add_argument("--aln-glob", required=True)
    c.add_argument("--ids", default=None,
                   help="expected taxon ids (strongly recommended)")
    c.add_argument("--out", required=True)
    c.add_argument("--gap-char", default="-",
                   help="filler for columns a sequence has no data in (default "
                        "'-'; '.' and '?' are also accepted).  The filler is "
                        "recorded in the report and in the output header.")
    c.set_defaults(func=cmd_concat)

    args = ap.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
