#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reading frame + strand normalisation for PROTEIN-CODING amplicon markers.

Why this exists: for a protein-coding marker (nifH, amoA, nirK, pmoA, nosZ, COI,
rbcL, ...) the amplicon sequences are fragments of a coding sequence.  If some
sequences sit in a different reading frame, or are reverse-complemented relative
to the others, an aligner cannot relate them: it lays the groups side by side and
the matrix explodes (measured: a mixed set produced 3,933 columns and 92% gaps,
while the same sequences normalised to one frame produced 455 columns and 29%
gaps).  Normalising is therefore a technical requirement, not a judgement call.

Nothing is discarded here.  Sequences whose frame cannot be determined are still
written out and only counted, so the taxon set stays under the study's control.

Two modes for `pick`:

  hmm  best domain hit across the six translated frames, using a profile HMM for
       the marker with trusted cutoffs (e.g. hmmsearch --cut_tc Fer4_NifH.hmm).
       Preferred: scores are orientation-sensitive, so the frame that wins also
       fixes the strand.
  orf  no reference needed: the frame with the longest stop-free stretch wins.
       Normalises frame reliably; strand is decided by which orientation has the
       longer ORF, which is a heuristic.  Use when no HMM is available, and say so
       in the methods.

The amino-acid position i of the written ORF corresponds exactly to nucleotides
[offset + 3i, offset + 3i + 3) of the oriented sequence, which is what
back-translation needs.  Translation is NOT truncated at the first stop; stop
codons become --stop-char (default X) so the whole frame stays usable.

Genetic code: default 1.  Mitochondrial markers need 2 (vertebrate) or 5
(invertebrate); see --genetic-code.
"""
import argparse
import os
import re
import sys

# NCBI translation tables, as the conventional 64-character amino-acid string.
# Order: TTT TTC TTA TTG TCT TCC TCA TCG TAT TAC TAA TAG TGT TGC TGA TGG
#        CTT CTC CTA CTG CCT CCC CCA CCG CAT CAC CAA CAG CGT CGC CGA CGG
#        ATT ATC ATA ATG ACT ACC ACA ACG AAT AAC AAA AAG AGT AGC AGA AGG
#        GTT GTC GTA GTG GCT GCC GCA GCG GAT GAC GAA GAG GGT GGC GGA GGG
CODES = {
    1: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    2: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    3: "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    4: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    5: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
    6: "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    9: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    11: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
}
# These tables were wrong once (table 2 had TGA as stop, table 9 had ATA as Ile),
# which would have mistranslated every mitochondrial COI run from the first step.
# selftest.py now asserts the distinguishing codons of each table; keep it green.
CODE_NOTES = {
    1: "standard",
    2: "vertebrate mitochondrial (TGA=W, ATA=M, AGA/AGG=stop)",
    3: "yeast mitochondrial (TGA=W, CTN=Thr, ATA=M)",
    4: "mold/protozoan/coelenterate mitochondrial (TGA=W)",
    5: "invertebrate mitochondrial (TGA=W, ATA=M, AGA/AGG=Ser) -- use for COI",
    6: "ciliate/dasycladacean/hexamita nuclear (TAA/TAG=Gln)",
    9: "echinoderm/flatworm mitochondrial (TGA=W, ATA=M, AGA/AGG=Ser, AAA=Asn)",
    11: "bacterial/archaeal/plant plastid (same assignments as standard)",
}
BASES = "TCAG"
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
FRAMES = [("p", 1, 0), ("p", 2, 1), ("p", 3, 2), ("m", 1, 0), ("m", 2, 1), ("m", 3, 2)]


def build_table(code_id):
    aas = CODES[code_id]
    table = {}
    i = 0
    for a in BASES:
        for b in BASES:
            for c in BASES:
                table[a + b + c] = aas[i]
                i += 1
    return table


def revcomp(s):
    return s.translate(COMP)[::-1]


def translate(s, table, stop_char="X"):
    out = []
    for j in range(0, len(s) - 2, 3):
        aa = table.get(s[j:j + 3].upper())
        out.append(stop_char if (aa is None or aa == "*") else aa)
    return "".join(out)


def stop_count(s, table):
    n = 0
    for j in range(0, len(s) - 2, 3):
        if table.get(s[j:j + 3].upper()) == "*":
            n += 1
    return n


def longest_stop_free(aa, stop_char="X"):
    best = cur = 0
    for ch in aa:
        if ch == stop_char:
            cur = 0
        else:
            cur += 1
            if cur > best:
                best = cur
    return best


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


# ---------------------------------------------------------------- translate --

def cmd_translate(args):
    table = build_table(args.genetic_code)
    os.makedirs(args.outdir, exist_ok=True)
    handles = {}
    wanted = set(args.frames.split(",")) if args.frames else None
    for st, fn, off in FRAMES:
        key = "%s%d" % (st, fn)
        if wanted and key not in wanted:
            continue
        handles[key] = open(os.path.join(args.outdir, "frame_%s.faa" % key),
                            "w", encoding="utf-8", newline="\n")

    n = 0
    rows = []
    for name, seq in read_fasta(args.fasta):
        n += 1
        rc = revcomp(seq)
        for st, fn, off in FRAMES:
            key = "%s%d" % (st, fn)
            if key not in handles:
                continue
            src = seq if st == "p" else rc
            sub = src[off:]
            aa = translate(sub, table, args.stop_char)
            handles[key].write(">%s\n%s\n" % (name, aa))
            rows.append((name, key, len(seq), len(aa),
                         stop_count(sub, table),
                         longest_stop_free(aa, args.stop_char)))

    for h in handles.values():
        h.close()

    print("translate: %d sequences x %d frames -> %s" %
          (n, len(handles), args.outdir))
    print("genetic code %d, stop char %r, translation NOT truncated"
          % (args.genetic_code, args.stop_char))
    for key in sorted(handles):
        p = os.path.join(args.outdir, "frame_%s.faa" % key)
        print("  %-4s %s (%d bytes)" % (key, p, os.path.getsize(p)))
    print()
    print("Next: score each frame against your marker profile HMM with trusted "
          "cutoffs, e.g.")
    print("  for f in %s/frame_*.faa; do b=$(basename $f .faa); "
          "hmmsearch --cut_tc --noali --domtblout %s/$b.domtblout "
          "<marker.hmm> $f > /dev/null; done" % (args.outdir, args.outdir))
    print("then:  frame_strand.py pick --mode hmm ...")
    return 0


# --------------------------------------------------------------------- pick --

FRAME_RE = re.compile(r"([pm])([123])")


def pick_orf(seq, table, stop_char):
    """Choose the best of the six frames without any reference.

    Longest stop-free stretch wins; ties break on the longer frame, then on the
    plus strand.  Used by --mode orf AND as the fallback for sequences that a
    HMM search did not resolve, so those sequences can still be written out.
    Returns (strand, frame, offset, aa).
    """
    rc = revcomp(seq)
    cand = []
    for st, fn, off in FRAMES:
        src = seq if st == "p" else rc
        aa = translate(src[off:], table, stop_char)
        lorf = longest_stop_free(aa, stop_char)
        cand.append((lorf, len(aa), st, fn, off, aa))
    cand.sort(key=lambda c: (-c[0], -c[1], 0 if c[2] == "p" else 1))
    _lorf, _aalen, st, fn, off, aa = cand[0]
    return st, fn, off, aa


def parse_domtblout(path):
    """Yield (seq_id, frame_key, full_evalue, full_score, dom_ievalue, dom_score).

    HMMER3 --domtblout columns (1-based): 1 target, 3 tlen, 4 query, 6 qlen,
    7 full E-value, 8 full score, 13 i-Evalue, 14 domain score, ...
    """
    out = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 14:
                continue
            try:
                out.append((f[0], float(f[6]), float(f[7]),
                            float(f[12]), float(f[13])))
            except ValueError:
                continue
    return out


def cmd_pick(args):
    table = build_table(args.genetic_code)
    seqs = list(read_fasta(args.fasta))
    if not seqs:
        print("ERROR: no sequences in %s" % args.fasta)
        return 2

    best = {}       # id -> (score, strand, frame, offset, evidence...)
    n_rows = 0

    if args.mode == "hmm":
        import glob as _glob
        # --domtblout-dir wins over --domtblout-glob: a directory needs no shell
        # globbing, so it is the form that works on every platform.
        if args.domtblout_dir:
            if not os.path.isdir(args.domtblout_dir):
                print("ERROR: --domtblout-dir is not a directory: %s"
                      % args.domtblout_dir)
                return 2
            files = sorted(_glob.glob(os.path.join(args.domtblout_dir,
                                                   "*.domtblout")))
            if not files:
                files = sorted(p for p in _glob.glob(os.path.join(args.domtblout_dir,
                                                                  "*"))
                               if os.path.isfile(p))
            if not files:
                print("ERROR: no domtblout files found in %s" % args.domtblout_dir)
                return 2
            print("domtblout directory        : %s (%d file(s))"
                  % (args.domtblout_dir, len(files)))
        elif args.domtblout_glob:
            files = sorted(_glob.glob(args.domtblout_glob))
            if not files:
                print("ERROR: no files matched %s" % args.domtblout_glob)
                return 2
        else:
            print("ERROR: --mode hmm requires --domtblout-dir (preferred) or "
                  "--domtblout-glob")
            return 2
        seen_frames = set()
        for path in files:
            m = FRAME_RE.search(os.path.basename(path))
            if not m:
                print("  !! cannot tell the frame from %s -- expected a filename "
                      "containing p1..m3 (e.g. frame_p2.domtblout); skipped"
                      % os.path.basename(path))
                continue
            st, fn = m.group(1), int(m.group(2))
            key = "%s%d" % (st, fn)
            seen_frames.add(key)
            hits = parse_domtblout(path)
            n_rows += len(hits)
            for (sid, fev, fsc, iev, dsc) in hits:
                if dsc < args.min_score:
                    continue
                prev = best.get(sid)
                if prev is None or dsc > prev[0]:
                    best[sid] = (dsc, st, fn, 0, fsc, fev, iev)
        print("frames seen in domtblout filenames : %s"
              % ", ".join(sorted(seen_frames)))
        print("domain rows parsed                 : %d" % n_rows)
        if seen_frames and len(seen_frames) < 6:
            print("  !! only %d of 6 frames present -- did you score every frame?"
                  % len(seen_frames))
    else:
        # ORF mode: no reference needed
        for name, seq in seqs:
            st, fn, off, aa = pick_orf(seq, table, args.stop_char)
            lorf = longest_stop_free(aa, args.stop_char)
            best[name] = (lorf, st, fn, off, None, None, None, aa)

    # ---- one pass builds the table AND the outputs ----
    # They must come from the same computation, otherwise the table can describe
    # sequences the FASTA does not contain (which is what used to happen).
    rows = []
    resolved = 0
    fallback = set()
    excluded = set()
    frame_counts = {}
    seqmap = dict(seqs)
    for name, seq in seqs:
        b = best.get(name)
        if b is None:
            if args.drop_unresolved:
                excluded.add(name)
                rows.append((name, len(seq), "NA", "NA", "NA", 0, 0, "NA",
                             "NA", "NA", "NA", "NA", "excluded"))
                continue
            # The HMM did not frame it, but the sequence is still a sequence.
            # Falling back to the reference-free ORF call keeps the taxon set
            # intact -- the gate promises not to drop taxa behind the user's back.
            st, fn, off, _aa = pick_orf(seq, table, args.stop_char)
            fallback.add(name)
            evidence = "orf-fallback"
            score = fsc = fev = iev = "NA"
        else:
            evidence = "hmm" if args.mode == "hmm" else "orf"
            if args.mode == "hmm":
                score, st, fn, _off, fsc, fev, iev = b
                off = {1: 0, 2: 1, 3: 2}[fn]
            else:
                score, st, fn, off = b[0], b[1], b[2], b[3]
                fsc = fev = iev = "NA"
            resolved += 1
        frame_counts[(st, fn)] = frame_counts.get((st, fn), 0) + 1
        src = seq if st == "p" else revcomp(seq)
        orf = src[off:]
        orf = orf[: len(orf) - (len(orf) % 3)]
        aa = translate(orf, table, args.stop_char)
        rows.append((name, len(seq), st, fn, off, len(orf), len(aa),
                     stop_count(orf, table), score, fsc, fev, iev, evidence))

    # write table
    if args.out_tsv:
        d = os.path.dirname(os.path.abspath(args.out_tsv))
        if d:
            os.makedirs(d, exist_ok=True)
        with open(args.out_tsv, "w", encoding="utf-8", newline="\n") as fh:
            fh.write("seq_id\tnt_len\tstrand\tframe\toffset\torf_nt_len\taa_len\t"
                     "stop_codons\tscore\thmm_full_score\thmm_full_evalue\t"
                     "hmm_dom_ievalue\tevidence\n")
            for r in rows:
                fh.write("\t".join(str(x) for x in r) + "\n")

    norm = []
    aarec = []
    for r in rows:
        name, _ntlen, st, fn, off = r[0], r[1], r[2], r[3], r[4]
        if r[12] == "excluded":
            continue
        src = seqmap[name] if st == "p" else revcomp(seqmap[name])
        orf = src[off:]
        orf = orf[: len(orf) - (len(orf) % 3)]
        norm.append((name, orf))
        aarec.append((name, translate(orf, table, args.stop_char)))
    if args.out_nt:
        write_fasta(norm, args.out_nt)
    if args.out_aa:
        write_fasta(aarec, args.out_aa)

    print()
    print("mode                       : %s" % args.mode)
    print("sequences                  : %d" % len(seqs))
    print("frame/strand resolved      : %d (%.1f%%)"
          % (resolved, 100.0 * resolved / len(seqs)))
    print("not resolved by the search : %d" % (len(seqs) - resolved))
    if fallback:
        print("   kept via ORF fallback   : %d  (evidence=orf-fallback in the table; "
              "no taxon was dropped)" % len(fallback))
    # ---- internal stop codons: the first signal of a numt or a frameshift ----
    # This is computed per sequence anyway (the evidence table has a stop_codons
    # column); it is reported here because a translated ORF with an INTERNAL stop is
    # not a plausible protein-coding sequence, and silently carrying it into a tree
    # turns an artefact into a "novel clade".  A stop at the very END is normal for a
    # complete CDS and is not counted.
    internal = []
    for r in rows:
        if r[12] == "excluded":
            continue
        src = seqmap[r[0]] if r[2] == "p" else revcomp(seqmap[r[0]])
        orf = src[r[4]:]
        orf = orf[: len(orf) - (len(orf) % 3)]
        aa = translate(orf, table, args.stop_char)
        if args.stop_char in aa[:-1]:
            internal.append(r[0])
    if internal:
        print("internal stop codons       : %d sequence(s) (%.1f%%) -- the first "
              "signal of a NUMT or a frameshift" % (len(internal),
                                                     100.0 * len(internal) / len(seqs)))
        print("   !! e.g. %s" % ", ".join(sorted(internal)[:8]))
        print("   A translated ORF with an INTERNAL stop is not a plausible")
        print("   protein-coding sequence.  Do NOT read it as a novel clade: either")
        print("   it is a numt/pseudogene (organelle markers: frequent), or the wrong")
        print("   frame was picked.  Options: inspect those sequences by eye, drop")
        print("   them as a STUDY decision under D4 (report the counts), or keep them")
        print("   and declare the artefact in the methods.")
        if args.out_tsv:
            lp = os.path.splitext(args.out_tsv)[0] + ".internal_stop.txt"
            with open(lp, "w", encoding="utf-8", newline="\n") as fh:
                fh.write("\n".join(sorted(internal)) + "\n")
            print("   list written to         : %s (LIST ONLY -- no taxon was "
                  "dropped)" % lp)
    if excluded:
        print("   EXCLUDED on request     : %d  (--drop-unresolved)" % len(excluded))
        print("   !! they are absent from --out-nt / --out-aa.  If your id list still "
              "contains them, G4's tip-set check WILL fail -- deliberately, so the")
        print("      exclusion cannot pass unnoticed.  Narrow the taxon list to match, "
              "and say so in the methods.")
    print("winning frame distribution :")
    for k in sorted(frame_counts, key=lambda x: -frame_counts[x]):
        print("    %s%d : %d" % (k[0], k[1], frame_counts[k]))
    if args.mode == "hmm":
        strands = {}
        for (st, fn), c in frame_counts.items():
            strands[st] = strands.get(st, 0) + c
        print("strand split               : %s" % ", ".join(
            "%s=%d" % (k, v) for k, v in sorted(strands.items())))
        print("  (all one strand is the expected result for amplicon data; a mixed "
              "split means the input mixes orientations)")
    if args.out_tsv:
        print("evidence table             : %s" % args.out_tsv)
    if args.out_nt:
        print("frame-normalised nt        : %s (%d records)" % (args.out_nt,
                                                                  len(norm)))
    if args.out_aa:
        print("translated aa              : %s (%d records)" % (args.out_aa,
                                                                  len(aarec)))
    print()
    print("Now align the NORMALISED sequences, then check:")
    print("  python scripts/aln_qc.py --aln aln.fasta --fasta-raw %s"
          % (args.out_nt or "norm.fasta"))
    return 0


def main():
    ap = argparse.ArgumentParser(
        description="Reading frame + strand normalisation (protein-coding markers).")
    sub = ap.add_subparsers(dest="cmd", required=True)

    t = sub.add_parser("translate", help="write all six translated frames")
    t.add_argument("--fasta", required=True)
    t.add_argument("--outdir", required=True)
    t.add_argument("--frames", default=None,
                   help="comma list to restrict, e.g. p1,p2,p3,m1,m2,m3")
    t.add_argument("--stop-char", default="X",
                   help="character for stop/unknown codons (default X, which "
                        "HMMER accepts; use '*' only if your tools expect it)")
    t.add_argument("--genetic-code", type=int, default=1, choices=sorted(CODES),
                   help="NCBI translation table: 1 standard, 2 vertebrate "
                        "mitochondrial, 5 invertebrate mitochondrial, 6 ciliate, "
                        "11 bacterial/plastid (default 1)")
    t.set_defaults(func=cmd_translate)

    p = sub.add_parser("pick", help="choose the winning frame/strand per sequence")
    p.add_argument("--fasta", required=True)
    p.add_argument("--frames-dir", default=None,
                   help="directory written by 'translate' (used to resolve "
                        "relative domtblout paths)")
    p.add_argument("--mode", choices=["hmm", "orf"], default="orf")
    p.add_argument("--domtblout-dir", default=None,
                   help="--mode hmm: a DIRECTORY of hmmsearch domtblout files named "
                        "so the frame is visible (frame_p1.domtblout ... "
                        "frame_m3.domtblout).  Preferred over --domtblout-glob: the "
                        "shell does not always expand a glob (Windows cmd/PowerShell "
                        "do not), so a directory is the portable form.")
    p.add_argument("--domtblout-glob", default=None,
                   help="--mode hmm: glob for hmmsearch domtblout files whose "
                        "FILENAMES contain the frame, e.g. 'frames/frame_*.domtblout'")
    p.add_argument("--drop-unresolved", action="store_true",
                   help="EXCLUDE sequences the HMM could not frame, instead of "
                        "keeping them via the ORF fallback.  This changes the taxon "
                        "set, so G4's tip-set check will complain until you narrow "
                        "the id list too.  Off by default: the gate does not drop "
                        "taxa silently.")
    p.add_argument("--min-score", type=float, default=0.0,
                   help="--mode hmm: ignore domain scores below this.  If you used "
                        "hmmsearch --cut_tc the files already contain only trusted "
                        "hits, so leave this at 0.")
    p.add_argument("--stop-char", default="X")
    p.add_argument("--genetic-code", type=int, default=1, choices=sorted(CODES))
    p.add_argument("--out-tsv", default=None)
    p.add_argument("--out-nt", default=None)
    p.add_argument("--out-aa", default=None)
    p.set_defaults(func=cmd_pick)

    args = ap.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
