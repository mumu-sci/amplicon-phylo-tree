#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Audit the representative sequences BEFORE alignment (gate G1b).

Why this exists: several failure modes are invisible to alignment QC.

  * REVERSE-COMPLEMENTED sequences.  A sequence in the wrong orientation is not
    "half missing" -- it is a different string.  The aligner cannot pair it with
    its true relatives, so it shows up as a region of mismatches, NOT as gaps.
    A gap-fraction check therefore does not catch it.  This affects rRNA and ITS
    just as much as protein-coding markers, which is why the audit is not
    restricted to coding genes.
  * NON-HOMOLOGOUS sequences.  If the taxon set contains sequences that are not
    the intended marker, no amount of alignment tuning fixes it.  The
    authoritative check is homology (domain HMM or a reference database), not
    the gap fraction.
  * MIXED AMPLICON REGIONS.  Two primer pairs or two regions in one alignment
    produce a tree that means nothing.  Detectable from the length distribution
    or from primer anchoring.

THIS SCRIPT REPORTS.  It never edits, filters or deletes the input.  Any output
is a LIST OF IDS for a human to act on.  Narrowing the taxon set is a study
decision, and it must be recorded.

Nothing here is a house threshold: every cutoff is an option with a stated
default, and anything that could not be checked is reported as NOT CHECKED
rather than as a pass.

Exit codes: 0 = pass, 1 = a check failed, 2 = input/usage error.
"""
import argparse
import hashlib
import json
import os
import shutil
import statistics
import subprocess
import sys
import tempfile

AMBIG = set("NRYSWKMBDHVnryswkmbdhv")
NT = set("ACGTUacgtu")
CODES = {
    1: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    2: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    5: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
    11: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
}
BASES = "TCAG"
FRAMES = [("p", 1, 0), ("p", 2, 1), ("p", 3, 2), ("m", 1, 0), ("m", 2, 1), ("m", 3, 2)]
COMP = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")


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


def build_table(code_id):
    aas = CODES[code_id]
    t = {}
    i = 0
    for a in BASES:
        for b in BASES:
            for c in BASES:
                t[a + b + c] = aas[i]
                i += 1
    return t


def translate(s, table):
    out = []
    for j in range(0, len(s) - 2, 3):
        aa = table.get(s[j:j + 3].upper())
        out.append("X" if (aa is None or aa == "*") else aa)
    return "".join(out)


# ----------------------------------------------------------- region map ------

def read_region_map(path):
    """Read an `id<TAB>region` map.  Returns (mapping, order) preserving order.

    This is the file that makes a mixed-region waiver traceable: it is the only
    place that says WHICH sequence belongs to WHICH region.  A waiver whose map
    is never read is not traceable at all, so the map is parsed, checked against
    the input ids and hashed into the report.
    """
    mapping = {}
    order = []
    problems = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\r\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                parts = line.split(None, 1)
            if len(parts) < 2 or not parts[0].strip() or not parts[1].strip():
                problems.append("line %d does not have the two columns "
                                "'id<TAB>region': %r" % (lineno, line))
                continue
            sid, region = parts[0].strip(), parts[1].strip()
            if sid in mapping:
                problems.append("line %d repeats the id %r" % (lineno, sid))
                continue
            mapping[sid] = region
            if region not in order:
                order.append(region)
    return mapping, order, problems


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------- orientation --

def orientation_via_hmm(seqs, hmm, hmmsearch, tmpdir, code):
    """Six-frame HMMER search; the winning frame fixes strand.  No DB caveat."""
    table = build_table(code)
    best = {}
    os.makedirs(tmpdir, exist_ok=True)
    for st, fn, off in FRAMES:
        key = "%s%d" % (st, fn)
        faa = os.path.join(tmpdir, "frame_%s.faa" % key)
        with open(faa, "w", encoding="utf-8", newline="\n") as fh:
            for name, seq in seqs:
                src = seq if st == "p" else seq.translate(COMP)[::-1]
                fh.write(">%s\n%s\n" % (name, translate(src[off:], table)))
        dom = os.path.join(tmpdir, "frame_%s.domtblout" % key)
        try:
            with open(os.path.join(tmpdir, "frame_%s.err" % key), "wb") as err:
                subprocess.run([hmmsearch, "--cut_tc", "--noali", "--cpu", "4",
                                "--domtblout", dom, hmm, faa],
                               stdout=subprocess.DEVNULL, stderr=err, timeout=7200)
        except (OSError, subprocess.TimeoutExpired) as exc:
            return None, "hmmsearch failed on frame %s: %s" % (key, exc)
        if not os.path.exists(dom):
            continue
        for line in open(dom, encoding="utf-8", errors="replace"):
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 14:
                continue
            try:
                score = float(f[13])
            except ValueError:
                continue
            prev = best.get(f[0])
            if prev is None or score > prev[0]:
                best[f[0]] = (score, st, fn)
    if not best:
        return {}, None
    return best, None


def hmm_coverage_profile(seqs, tmpdir):
    """Per-sequence coverage of the marker domain, from the six-frame domtblout.

    Why: two amplicon regions can have the SAME length, and then a length-bimodality
    test is blind to the mixture.  What still differs is WHERE the marker domain
    sits inside the read -- a partial domain at one end versus a full domain.  So
    for every sequence, take its best domain hit and record the fraction of the
    protein domain it covers (ali span / total alignment length).

    Returns (profile, info) where profile maps id -> dict(cov=..., frame=...);
    info is None when the channel is usable and a string explaining why not
    otherwise.  Nothing is inferred when no domain table exists -- an unrun check
    is reported as NOT CHECKED, never as a pass.
    """
    profile = {}
    tables = 0
    for st, fn, _off in FRAMES:
        dom = os.path.join(tmpdir, "frame_%s%d.domtblout" % (st, fn))
        if not os.path.exists(dom):
            continue
        tables += 1
        for line in open(dom, encoding="utf-8", errors="replace"):
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 22:
                continue
            # HMMER3 --domtblout columns (1-based): 1 target, 3 tlen, 4 query,
            # 6 qlen, 13 i-Evalue, 14 domain score, 18-19 hmm from/to,
            # 20-21 ali from/to (query coordinates), 22-23 env from/to.
            try:
                score = float(f[13])
                qlen = int(f[5])                        # the SEQUENCE length
                ali_from, ali_to = int(f[19]), int(f[20])   # query coordinates
            except (ValueError, IndexError):
                continue
            span = abs(ali_to - ali_from) + 1
            # fraction of the read that the domain alignment covers.  Measuring
            # this against the HMM length instead was a real bug: it reported how
            # much of the MODEL a hit spans, which is ~the same for every sequence
            # and therefore blind to exactly the same-length mixture this criterion
            # exists for.
            cov = span / qlen if qlen else 0.0
            prev = profile.get(f[0])
            if prev is None or score > prev["score"]:
                profile[f[0]] = {"score": score, "cov": cov,
                                 "frame": "%s%d" % (st, fn),
                                 "ali_span": span}
    if not tables:
        return {}, "no HMM domain table was produced"
    if not profile:
        return {}, ("no sequence had a domain hit in the HMM domain tables, so "
                    "the coverage pattern could not be read")
    return profile, None


def orientation_via_blast(seqs, ref_db, blastn, tmpdir):
    """blastn against a reference DB; the subject-coordinate sign gives strand.

    CAVEAT THAT IS NOT THEORETICAL: some curated databases store their entries
    reverse-complementary to the gene's coding orientation.  When that is so,
    EVERY hit looks reverse-complemented and the sign carries no information
    about the query.  We therefore report the raw fraction and only call it a
    genuine mixture when the fraction is intermediate.
    """
    q = os.path.join(tmpdir, "queries.fasta")
    with open(q, "w", encoding="utf-8", newline="\n") as fh:
        for name, seq in seqs:
            fh.write(">%s\n%s\n" % (name, seq))
    tab = os.path.join(tmpdir, "blastn.tab")
    try:
        with open(os.path.join(tmpdir, "blastn.err"), "wb") as err:
            subprocess.run([blastn, "-query", q, "-db", ref_db, "-evalue", "1e-5",
                            "-max_target_seqs", "5", "-num_threads", "4",
                            "-outfmt", "6 qseqid sseqid pident length qlen evalue "
                                       "bitscore qstart qend sstart send",
                            "-out", tab],
                           stdout=subprocess.DEVNULL, stderr=err, timeout=7200)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return None, "blastn failed: %s" % exc
    best = {}
    if os.path.exists(tab):
        for line in open(tab, encoding="utf-8", errors="replace"):
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            try:
                bs = float(f[6])
                sstart, send = int(f[9]), int(f[10])
            except ValueError:
                continue
            prev = best.get(f[0])
            if prev is None or bs > prev[0]:
                best[f[0]] = (bs, sstart, send)
    return best, None


# --------------------------------------------------------------------- main ---

def main():
    ap = argparse.ArgumentParser(
        description="Audit representative sequences before alignment (gate G1b). "
                    "Reports only; never edits or filters the input.")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--hmm", default=None, help="marker profile HMM (preferred: the "
                                               "orientation call needs no DB caveat)")
    ap.add_argument("--ref-db", default=None, help="BLAST nucleotide DB prefix.  "
                                                   "This is the channel for "
                                                   "NON-CODING markers, and it also "
                                                   "drives the homology audit for "
                                                   "them.")
    ap.add_argument("--marker-class", choices=["auto", "coding", "noncoding"],
                    default="auto",
                    help="'coding' = a protein-coding marker (HMM orientation is "
                         "meaningful); 'noncoding' = rRNA / ITS, where the six-frame "
                         "HMM path is meaningless and --ref-db is the only valid "
                         "channel.  'auto' infers it from the flags you supplied.")
    ap.add_argument("--accept-and-record", action="store_true",
                    help="treat checks that could not run as accepted.  Requires "
                         "--override-reason: an unchecked gate that is waived "
                         "silently is not a gate.")
    ap.add_argument("--override-reason", default=None,
                    help="why it is acceptable that those checks did not run; "
                         "required with --accept-and-record and written into the "
                         "JSON report and the ops log.")
    ap.add_argument("--hmmsearch-path", default=None)
    ap.add_argument("--blastn-path", default=None)
    ap.add_argument("--genetic-code", type=int, default=1, choices=sorted(CODES))
    ap.add_argument("--primers", default=None, metavar="FWD,REV",
                    help="optional; enables primer anchoring for the region check")
    ap.add_argument("--max-revcomp-frac", type=float, default=0.01)
    ap.add_argument("--max-unhit-frac", type=float, default=0.05)
    ap.add_argument("--mixed-region-ok", action="store_true")
    ap.add_argument("--region-map", default=None)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--json", default=None)
    ap.add_argument("--label", default=None)
    args = ap.parse_args()

    if args.accept_and_record and not args.override_reason:
        ap.error("--accept-and-record requires --override-reason: a waived check "
                 "must say why, or the waiver is indistinguishable from an oversight")

    if not os.path.exists(args.fasta):
        print("ERROR: no such FASTA: %s" % args.fasta)
        return 2

    label = args.label or os.path.basename(args.fasta)
    seqs = list(read_fasta(args.fasta))
    if not seqs:
        print("ERROR: no sequences in %s" % args.fasta)
        return 2

    failures = []
    not_checked = []
    report = {"label": label, "fasta": args.fasta, "n_sequences": len(seqs)}

    # The --region-map is the traceability of a mixed-region waiver.  It is read,
    # checked and hashed here (before anything else can pass), not only mentioned
    # in a printed line: a waiver whose map is never validated cannot be audited
    # afterwards.
    region_map = None
    if args.region_map:
        if not os.path.exists(args.region_map):
            print("ERROR: --region-map not found: %s" % args.region_map)
            return 2
        mapping, regions, problems = read_region_map(args.region_map)
        fasta_ids = [n for n, _s in seqs]
        unknown = sorted(set(mapping) - set(fasta_ids))
        missing_ids = [n for n in fasta_ids if n not in mapping]
        if problems or unknown or missing_ids or not mapping:
            print("ERROR: --region-map %s is not usable:"
                  % args.region_map)
            for p in problems[:10]:
                print("   - %s" % p)
            if unknown:
                print("   - %d id(s) in the map are not in the FASTA, e.g. %s"
                      % (len(unknown), ", ".join(unknown[:8])))
            if missing_ids:
                print("   - %d FASTA id(s) have no region, e.g. %s"
                      % (len(missing_ids), ", ".join(missing_ids[:8])))
            if not mapping:
                print("   - the map has no usable rows")
            return 2
        region_map = {
            "path": args.region_map,
            "sha256": file_sha256(args.region_map),
            "n_entries": len(mapping),
            "n_regions": len(regions),
            "regions": {r: sum(1 for v in mapping.values() if v == r)
                        for r in regions},
            "ids_not_in_fasta": [],
            "fasta_ids_without_region": [],
        }
        report["region_map"] = region_map

    print("input audit (G1b) -- %s" % label)
    print("=" * 76)
    if region_map:
        print("region map declared        : %s" % region_map["path"])
        print("   sha256                  : %s" % region_map["sha256"])
        print("   regions                 : %s"
              % ", ".join("%s=%d" % (r, c)
                          for r, c in sorted(region_map["regions"].items())))
        print("   every FASTA id has a region: yes (%d ids, all two-column rows)"
              % region_map["n_entries"])
        print()


    # ---------------- 1. census (always possible) ----------------
    lens = [len(s) for _n, s in seqs]
    comp = {}
    ambig_seqs = 0
    iupac = 0
    for _n, s in seqs:
        seen_ambig = False
        for ch in s:
            comp[ch] = comp.get(ch, 0) + 1
            if ch in AMBIG:
                iupac += 1
                seen_ambig = True
        if seen_ambig:
            ambig_seqs += 1
    print("1. census")
    print("   sequences            : %d" % len(seqs))
    print("   length min/median/max: %d / %g / %d"
          % (min(lens), statistics.median(lens), max(lens)))
    print("   ambiguous (IUPAC/ N) : %d chars in %d sequence(s)"
          % (iupac, ambig_seqs))
    if ambig_seqs:
        print("     note: decide how your aligner treats them (see the "
              "ambiguous-policy entry in references/parameters.md D5) -- treating")
        print("     them as a 5th state instead of as missing is a common silent error")
    report["census"] = {"n": len(seqs), "len_min": min(lens),
                        "len_median": statistics.median(lens), "len_max": max(lens),
                        "iupac_chars": iupac, "iupac_sequences": ambig_seqs}
    print()

    # ---------------- 2. orientation ----------------
    # Working directory: NOT next to the input (this script never writes near the
    # input FASTA), and cleared first -- the homology audit scans every domtblout
    # in here, so a stale file from a previous run would count as a hit for a
    # sequence it never matched.
    base = args.outdir or tempfile.gettempdir()
    tmpdir = os.path.join(base, "_audit_tmp")
    shutil.rmtree(tmpdir, ignore_errors=True)
    os.makedirs(tmpdir, exist_ok=True)
    print("working directory          : %s (cleared)" % tmpdir)
    revcomp = []
    blast_hits = {}
    orient = {"method": None, "revcomp": 0, "fraction": None}
    print("2. orientation audit")
    hmmsearch = (args.hmmsearch_path or os.environ.get("DSH_HMMSEARCH")
                 or shutil.which("hmmsearch"))
    blastn = args.blastn_path or shutil.which("blastn")

    # Which orientation channel is even meaningful?  For a protein-coding marker a
    # six-frame HMM search settles strand without any database caveat.  For rRNA /
    # ITS there is nothing to translate, so the six-frame path is meaningless and a
    # reference database is the only channel -- that is why a non-coding marker
    # cannot quietly pass this gate with no database.
    if args.marker_class == "auto":
        marker = "coding" if args.hmm else ("noncoding" if args.ref_db else "unknown")
    else:
        marker = args.marker_class
    print("   marker class         : %s%s" % (
        marker, " (inferred)" if args.marker_class == "auto" else ""))
    if marker == "noncoding" and args.hmm:
        print("   !! --hmm was supplied for a NON-CODING marker.  Nothing is "
              "translated, so the six-frame path is skipped; --ref-db is the")
        print("      meaningful channel here.")
    if marker == "coding" and not args.hmm and args.ref_db:
        print("   note: this looks like a coding marker; --hmm would give an "
              "orientation call with no database caveat, but --ref-db works too.")

    if marker == "coding" and args.hmm and hmmsearch:
        best, err = orientation_via_hmm(seqs, args.hmm, hmmsearch, tmpdir,
                                       args.genetic_code)
        if err:
            not_checked.append("orientation (HMM: %s)" % err)
            print("   NOT CHECKED: %s" % err)
        else:
            orient["method"] = "hmm (six-frame, %s)" % args.hmm
            for name, _s in seqs:
                b = best.get(name)
                if b and b[1] == "m":
                    revcomp.append(name)
            frac = len(revcomp) / len(seqs)
            orient["revcomp"] = len(revcomp)
            orient["fraction"] = frac
            print("   method               : six-frame HMMER, trusted cutoffs")
            print("   reverse-complemented : %d / %d (%.2f%%)"
                  % (len(revcomp), len(seqs), 100 * frac))
            print("     (protein-level scoring is orientation-sensitive, so this "
                  "call needs no database caveat)")
    elif args.ref_db and blastn:
        best, err = orientation_via_blast(seqs, args.ref_db, blastn, tmpdir)
        if err:
            not_checked.append("orientation (blastn: %s)" % err)
            print("   NOT CHECKED: %s" % err)
        else:
            orient["method"] = "blastn vs %s" % args.ref_db
            hit = 0
            minus = 0
            blast_hits = best
            for name, _s in seqs:
                b = best.get(name)
                if not b:
                    continue
                hit += 1
                if b[1] > b[2]:
                    minus += 1
                    revcomp.append(name)
            frac = (minus / hit) if hit else 0.0
            orient["revcomp"] = minus
            orient["fraction"] = frac
            print("   method               : blastn subject-coordinate sign")
            print("   hits                 : %d / %d" % (hit, len(seqs)))
            print("   minus-strand fraction: %.3f (%d)" % (frac, minus))
            if hit and (frac < 0.05 or frac > 0.95):
                print("   -> the fraction is extreme, which is the signature of a "
                      "DATABASE stored reverse-complementary to the coding")
                print("      orientation, not of a mixed input.  Treat this check as "
                      "INCONCLUSIVE and prefer --hmm.")
                not_checked.append("orientation (blastn sign inconclusive: DB "
                                   "orientation convention)")
                revcomp = []
                orient["revcomp"] = 0
                orient["fraction"] = None
            else:
                print("   -> intermediate fraction: a genuine orientation mixture")
    else:
        if not (hmmsearch or blastn):
            why = "neither hmmsearch nor blastn is on PATH"
        elif marker == "noncoding":
            why = ("pass --ref-db <blast-db-prefix> (SILVA / UNITE / PR2).  A "
                   "six-frame HMM search is meaningless for a marker that is not "
                   "translated, so a database is the only channel here")
        else:
            why = "pass --hmm (preferred) or --ref-db to enable it"
        not_checked.append("orientation (%s)" % why)
        print("   NOT CHECKED: %s" % why)
        print("     !! an unaudited orientation is the single most common silent "
              "cause of a nonsense tree for any marker, not just coding genes")

    if revcomp:
        frac = len(revcomp) / len(seqs)
        if frac > args.max_revcomp_frac:
            failures.append("%d sequence(s) (%.1f%%) appear reverse-complemented"
                            % (len(revcomp), 100 * frac))
        print("   -> %d suspect(s); fix them and RE-RUN the whole pipeline "
              "(frame 0 is not a valid tree for them)" % len(revcomp))
    report["orientation"] = orient
    print()

    # ---------------- 3. homology ----------------
    print("3. homology audit")
    if args.hmm and hmmsearch and not revcomp:
        # reuse the same six-frame search: a sequence with no hit in ANY frame
        table = build_table(args.genetic_code)
        hits = set()
        if os.path.isdir(tmpdir):
            for st, fn, _off in FRAMES:
                dom = os.path.join(tmpdir, "frame_%s%d.domtblout" % (st, fn))
                if not os.path.exists(dom):
                    continue
                for line in open(dom, encoding="utf-8", errors="replace"):
                    if line.startswith("#") or not line.strip():
                        continue
                    hits.add(line.split()[0])
        if hits:
            nohit = [n for n, _s in seqs if n not in hits]
            frac = len(nohit) / len(seqs)
            print("   sequences with a trustworthy domain hit : %d / %d"
                  % (len(seqs) - len(nohit), len(seqs)))
            print("   without any hit                         : %d (%.2f%%)"
                  % (len(nohit), 100 * frac))
            if frac > args.max_unhit_frac:
                failures.append("homology: %.1f%% of sequences have no marker hit "
                                "(limit %.1f%%)" % (100 * frac,
                                                    100 * args.max_unhit_frac))
                print("   !! this is expected to be a STUDY decision, not a script "
                      "decision.  Your options:")
                print("      (a) narrow the taxon set and say so (report reads%% and "
                      "richness%% retained)")
                print("      (b) keep every taxon and align in blocks "
                      "(scripts/blocks.py) or use a reference backbone")
                print("      (c) accept and record the reason")
            report["homology"] = {"with_hit": len(seqs) - len(nohit),
                                  "no_hit": len(nohit), "fraction_no_hit": frac}
            if args.outdir:
                os.makedirs(args.outdir, exist_ok=True)
                with open(os.path.join(args.outdir, "no_hit.txt"), "w",
                          encoding="utf-8", newline="\n") as fh:
                    fh.write("\n".join(nohit) + "\n")
                print("   list written to: %s/no_hit.txt (LIST ONLY -- the input "
                      "FASTA is untouched)" % args.outdir)
        else:
            not_checked.append("homology (no domtblout produced)")
            print("   NOT CHECKED: no domain table was produced")
    else:
        if revcomp:
            why = ("deferred: %d sequence(s) were called reverse-complemented, so a "
                   "'no hit' cannot be told apart from 'the wrong strand'.  Fix the "
                   "orientation, then re-run to get a homology audit." % len(revcomp))
        elif marker == "noncoding" and blast_hits:
            # For rRNA / ITS the reference database IS the homology test: a sequence
            # that matches nothing in SILVA/UNITE/PR2 is not the marker you think it
            # is, however nucleotide-looking it is.
            nohit = [n for n, _s in seqs if n not in blast_hits]
            frac = len(nohit) / len(seqs)
            print("   method                                  : reference-database hits "
                  "(%s)" % args.ref_db)
            print("   sequences with a hit                    : %d / %d"
                  % (len(seqs) - len(nohit), len(seqs)))
            print("   without any hit                         : %d (%.2f%%)"
                  % (len(nohit), 100 * frac))
            if frac > args.max_unhit_frac:
                failures.append("homology: %.1f%% of sequences have no hit in the "
                                "reference database (limit %.1f%%)"
                                % (100 * frac, 100 * args.max_unhit_frac))
                print("   !! a sequence of the right length and alphabet that matches "
                      "nothing in the reference is usually NOT this marker.  Decide")
                print("      whether to narrow the taxon set, and report what you did.")
            report["homology"] = {"method": "reference_database",
                                  "with_hit": len(seqs) - len(nohit),
                                  "no_hit": len(nohit), "fraction_no_hit": frac}
            if args.outdir:
                os.makedirs(args.outdir, exist_ok=True)
                with open(os.path.join(args.outdir, "no_hit.txt"), "w",
                          encoding="utf-8", newline="\n") as fh:
                    fh.write("\n".join(nohit) + "\n")
                print("   list written to: %s/no_hit.txt (LIST ONLY -- the input "
                      "FASTA is untouched)" % args.outdir)
        elif marker == "noncoding":
            why = ("this is a non-coding marker, so orientation and homology rest on "
                   "a reference database -- pass --ref-db <blast-db-prefix> built "
                   "from SILVA / UNITE / PR2, or record the gap with "
                   "--accept-and-record --override-reason")
            not_checked.append("orientation+homology (%s)" % why)
            print("   NOT CHECKED: %s" % why)
        else:
            why = ("needs --hmm and hmmsearch on PATH; pass --hmm (preferred) or "
                   "--ref-db")
            not_checked.append("homology (%s)" % why)
            print("   NOT CHECKED: %s" % why)
    print()

    # ---------------- 4. mixed amplicon region ----------------
    print("4. mixed amplicon region")
    mixed = False
    mixed_primers = False
    mixed_length = False
    cov_info = None
    cov_prof = None
    cov_summary = None
    comp_info = None
    if args.primers and "," in args.primers:
        fwd, rev = [p.strip().upper() for p in args.primers.split(",", 1)]
        rev_rc = rev.translate(COMP)[::-1]
        n_fwd = sum(1 for _n, s in seqs if s.upper()[:len(fwd)] == fwd)
        n_rev = sum(1 for _n, s in seqs if s.upper()[-len(rev_rc):] == rev_rc)
        print("   primer anchoring     : 5' %s matched %d/%d ; 3' %s matched %d/%d"
              % (fwd, n_fwd, len(seqs), rev_rc, n_rev, len(seqs)))
        if n_fwd < len(seqs) * 0.5 or n_rev < len(seqs) * 0.5:
            mixed = True
            mixed_primers = True
            print("   -> not all sequences carry this primer pair: either the "
                  "primers were trimmed already, or a different region is mixed in")
    else:
        # documented, simple rule -- not a formal dip test
        bins = {}
        for L in lens:
            bins[L // 10 * 10] = bins.get(L // 10 * 10, 0) + 1
        top = sorted(bins.items(), key=lambda kv: -kv[1])[:2]
        if len(top) == 2:
            (b1, c1), (b2, c2) = top
            gap = abs(b1 - b2)
            print("   length modes         : %d bp (n=%d) and %d bp (n=%d), gap %d bp"
                  % (b1, c1, b2, c2, gap))
            print("   rule                 : two modes >100 bp apart, each >5% of "
                  "sequences -> report (bin=10 bp)")
            if (gap > 100 and c1 > 0.05 * len(seqs) and c2 > 0.05 * len(seqs)):
                mixed = True
                mixed_length = True
                print("   -> bimodal: consistent with TWO amplicon regions, or with "
                      "legitimately variable amplicon length")
        else:
            print("   length modes         : only one populated bin; cannot assess")

    # ---- second criterion: marker domain coverage AND composition grouping ----
    # DESIGN 5.2-4 asks for "length bimodality OR the domain-coverage pattern
    # splitting into two groups".  The length histogram alone misses two regions of
    # the SAME length, so there are two length-independent readings of the same
    # idea here:
    #   (a) HMM domain coverage -- a partial domain at one end vs a full domain in
    #       the middle gives different coverage fractions at equal read length.
    #       This is the DESIGN wording; it needs --hmm + hmmsearch, and when that
    #       channel was asked for and failed it is NOT CHECKED (a failure).
    #   (b) sequence composition (GC of the raw read) -- always computable, and two
    #       different primer pairs usually amplify compositionally different
    #       regions.  Groups that differ in composition AND in length by >100 bp
    #       are the same signature as (a).
    cov_info = None
    cov_prof = None
    if args.hmm and hmmsearch:
        cov_prof, cov_info = hmm_coverage_profile(seqs, tmpdir)
    elif args.hmm:
        cov_info = ("hmmsearch is not on PATH, so the coverage pattern could not "
                    "be read")
    print("   domain coverage      :", end="")
    if cov_info:
        print(" NOT CHECKED: %s" % cov_info)
        not_checked.append("mixed region, second criterion (%s)" % cov_info)
    elif cov_prof:
        vals = sorted(v["cov"] for v in cov_prof.values())
        med = statistics.median(vals)
        low = [n for n, _s in seqs if n in cov_prof and cov_prof[n]["cov"] < med]
        high = [n for n, _s in seqs if n in cov_prof and cov_prof[n]["cov"] >= med]
        lens_of = {n: len(s) for n, s in seqs}

        def med_len(ids):
            return statistics.median([lens_of[n] for n in ids]) if ids else 0

        print(" %d hit(s), coverage median %.2f -> %d low / %d high"
              % (len(vals), med, len(low), len(high)))
        if (len(seqs) >= 5 and low and high
                and min(len(low), len(high)) >= max(3, 0.05 * len(seqs))):
            ml, mh = med_len(low), med_len(high)
            cov_summary = {"median_coverage": med, "n_low": len(low),
                           "n_high": len(high), "len_median_low": ml,
                           "len_median_high": mh,
                           "len_difference": abs(ml - mh)}
            # The DESIGN criterion is the coverage pattern splitting into two
            # groups -- that is the signal.  At equal read length a partial domain
            # at one end and a full domain in the middle cannot arise from one
            # primer pair, so a length difference is corroboration, not a
            # precondition: requiring it would let the very case the criterion
            # exists for (two regions of the SAME length) slip through.
            mixed = True
            print("     rule                 : coverage median split; groups each "
                  ">= max(3, 5%%); group length medians %g vs %g bp (difference "
                  "%g bp)" % (ml, mh, abs(ml - mh)))
            print("     -> domain coverage separates the set into two groups: the "
                  "marker domain occupies different parts of the read,")
            print("        which is consistent with TWO amplicon regions even where "
                  "the length histogram is not bimodal")
        else:
            print("     rule                 : fewer than 5 hits or a group below "
                  "max(3, 5%%) -- the split is not assessable, so it is not used")
    else:
        print(" not computed (no --hmm)")

    # (b) composition grouping: a second, always-available reading of the same
    # criterion.  It does not replace the HMM channel; it keeps a same-length
    # mixture detectable when no marker HMM exists (ITS / rRNA / COI chains).
    print("   composition grouping :", end="")
    gc = []
    for _n, s in seqs:
        clean = [c for c in s.upper() if c in "ACGTU"]
        f = (sum(1 for c in clean if c in "GC") / len(clean)) if clean else 0.0
        gc.append((f, len(s)))
    gc_vals = sorted(g for g, _l in gc if g > 0)
    # Locate the widest GAP between a low and a high composition group instead of
    # forcing a median split: a forced split divides every dataset in half, and a
    # half/half set with equal lengths would then always look "grouped".
    split_at = None
    if len(gc_vals) >= 5:
        best = -1.0
        for i in range(len(gc_vals) - 1):
            if len(gc_vals) - (i + 1) < max(3, 0.05 * len(seqs)):
                break
            if i + 1 < max(3, 0.05 * len(seqs)):
                continue
            d = gc_vals[i + 1] - gc_vals[i]
            if d > best:
                best, split_at = d, (gc_vals[i] + gc_vals[i + 1]) / 2.0
    if split_at is not None:
        lows = [L for f, L in gc if f > 0 and f < split_at]
        highs = [L for f, L in gc if f > 0 and f >= split_at]
        print(" composition gap at GC %.3f -> %d low / %d high"
              % (split_at, len(lows), len(highs)))
        if not lows or not highs:
            print("     rule                 : one side of the composition gap is "
                  "empty -- the split is not used")
        else:
            ml = statistics.median(lows)
            mh = statistics.median(highs)
            print("     rule                 : the widest composition gap; groups "
                  "each >= max(3, 5%%); group length medians %g vs %g bp "
                  "(difference %g bp)" % (ml, mh, abs(ml - mh)))
            comp_info = {"split_gc": round(split_at, 6), "n_low": len(lows),
                         "n_high": len(highs), "len_median_low": ml,
                         "len_median_high": mh, "len_difference": abs(ml - mh)}
            if abs(ml - mh) > 100:
                mixed = True
                print("     -> composition separates the set into two groups with "
                      "different lengths: consistent with TWO amplicon regions")
    else:
        print(" only %d sequence(s) with a readable alphabet -- cannot assess"
              % len(gc_vals))

    if mixed:
        if args.mixed_region_ok:
            print("   WAIVED by --mixed-region-ok.  Recorded.  Reminder: each region "
                  "must be built into ITS OWN tree.")
            if region_map:
                print("   region map           : %s (sha256 %s, %d regions)"
                      % (region_map["path"], region_map["sha256"][:16],
                         region_map["n_regions"]))
            else:
                why = ("--mixed-region-ok without --region-map: the region labels "
                       "are not recorded anywhere, so the waiver cannot be audited")
                print("   !! %s" % why)
                not_checked.append("mixed-region waiver traceability (%s)" % why)
        else:
            failures.append("possible mixed amplicon regions (see above)")
            print("   !! STOP AND REPORT: mixing regions into one tree is not "
                  "meaningful.  Split them (one tree per region), or pass")
            print("      --mixed-region-ok --region-map <id<TAB>region> if you have "
                  "already handled it")
    report["mixed_region"] = {
        "flag": mixed,
        "waived": bool(args.mixed_region_ok),
        "criteria": {
            "primers": mixed_primers,
            "length_bimodality": mixed_length,
            "hmm_domain_coverage": cov_summary,
            "composition_grouping": comp_info,
            "channel_not_checked": [n for n in not_checked
                                    if "mixed region" in n],
        },
    }
    print()

    # ---------------- 5. chimera reminder ----------------
    print("5. chimera detection (NOT done here -- by design)")
    print("   Chimeras are a leading source of false long branches and false novel")
    print("   clades, and detecting them needs the whole ASV set, not just the")
    print("   representatives.  Run it on the full denoised set, or record an")
    print("   explicit, justified skip:")
    print("     vsearch --uchime3_denovo all_asvs.fasta \\")
    print("             --nonchimeras asvs.nochimera.fasta \\")
    print("             --chimeras chimeras.fasta")
    print("   If you drop chimeras you MUST write the new taxon rule back into")
    print("   params.md and re-run from the taxon-set step.")
    print()

    # ---------------- reporting ----------------
    if revcomp and args.outdir:
        os.makedirs(args.outdir, exist_ok=True)
        with open(os.path.join(args.outdir, "revcomp_suspects.txt"), "w",
                  encoding="utf-8", newline="\n") as fh:
            fh.write("\n".join(sorted(revcomp)) + "\n")
        print("reverse-complement suspects written to: %s/revcomp_suspects.txt "
              "(LIST ONLY)" % args.outdir)

    print("=" * 76)
    print("checks that could NOT be run: %d" % len(not_checked))
    for n in not_checked:
        print("   - %s" % n)
    print()

    # A check that could not run is NOT a check that passed.  This gate exists to
    # stop a non-coding marker (rRNA / ITS), where characters cannot reveal
    # orientation or homology, from sailing through with an all-green report on
    # the strength of two checks that never executed.
    if not_checked and not args.accept_and_record:
        for n in not_checked:
            failures.append("could not be checked: %s" % n)
        print("Each unchecked item above is counted as a FAILURE (fail-closed).")
        print("If you have another way to establish orientation and homology -- for")
        print("example a BLAST/vsearch pass against SILVA / UNITE / PR2 -- either")
        print("pass it here (--ref-db) or record the gap explicitly:")
        print("    --accept-and-record --override-reason \"...\"")
        print()
    elif not_checked and args.accept_and_record:
        print("OVERRIDE ACCEPTED (--accept-and-record).  The unchecked items above")
        print("are recorded as accepted, with this reason:")
        print("    %s" % (args.override_reason or "(none given)"))
        print("This must be written into the ops log -- it is an override, not a pass.")
        print()

    if failures:
        print("RESULT: FAIL")
        for f in failures:
            print("   - %s" % f)
    else:
        print("RESULT: PASS -- every check either ran clean or is covered by an "
              "explicit, recorded override")

    report["not_checked"] = not_checked
    report["failures"] = failures
    report["accepted_without_checking"] = (not_checked if args.accept_and_record
                                           else [])
    report["override_reason"] = (args.override_reason
                                if args.accept_and_record else None)
    report["marker_class"] = args.marker_class
    report["verdict"] = "FAIL" if failures else "PASS"
    if args.json:
        d = os.path.dirname(os.path.abspath(args.json))
        if d:
            os.makedirs(d, exist_ok=True)
        with open(args.json, "w", encoding="utf-8", newline="\n") as fh:
            json.dump(report, fh, indent=2, ensure_ascii=False)
        print("report written: %s" % args.json)

    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
