#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Probe the phylogenetics toolchain and record exact versions.

Why this exists: a tree is only reproducible if the tool versions that produced it
are recorded next to it.  Papers get rejected over missing/ambiguous software
versions, and silent version differences (e.g. FastTree precision variants) change
branch lengths.  Run this first and keep the JSON next to the tree.

Cross-platform: tools are located on PATH (or overridden with --set NAME=PATH).
Nothing is installed, nothing is downloaded.

Exit codes: 0 = all required tools present, 1 = something required is missing.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

# name -> argv that makes the tool print its version/usage.
# Tools whose version banner is only reachable via a usage flag are noted.
PROBES = {
    "mafft": ["--version"],
    "fasttree": ["-help"],
    "FastTreeMP": ["-help"],
    "iqtree2": ["--version"],
    "iqtree": ["-v"],
    "RAxML": ["-v"],
    "raxml-ng": ["--version"],
    "hmmsearch": ["-h"],
    "hmmalign": ["-h"],
    "blastn": ["-version"],
    "blastx": ["-version"],
    "tblastn": ["-version"],
    "makeblastdb": ["-version"],
    "Rscript": ["--version"],
    "python3": ["--version"],
    "vsearch": ["--version"],
    "cd-hit": ["-h"],
    "usearch": ["-version"],
    "trimal": ["--version"],
    "seqkit": ["version"],
    "ITSx": ["-h"],
    "muscle": ["-version"],
    "clipkit": ["--version"],
}

# Which of the above are plausibly required, grouped by role.  A run needs one
# aligner, one tree program, and (for protein-coding markers) one domain search.
ROLES = {
    "aligner": ["mafft", "muscle"],
    "tree": ["fasttree", "FastTreeMP", "iqtree2", "iqtree", "RAxML", "raxml-ng"],
    "domain_search": ["hmmsearch", "hmmalign"],
    "blast": ["blastn", "blastx", "makeblastdb"],
    "r": ["Rscript"],
}


def first_meaningful_line(text, limit=160):
    for line in (text or "").splitlines():
        line = line.strip()
        if line:
            return line[:limit]
    return ""


def probe_one(name, argv, resolved, timeout=25):
    if not resolved:
        return {"found": False, "path": None, "version_line": None, "exit": None}
    cmd = [resolved] + argv
    try:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           timeout=timeout)
        out = p.stdout.decode("utf-8", errors="replace")
        rc = p.returncode
    except subprocess.TimeoutExpired:
        return {"found": True, "path": resolved, "version_line": "TIMEOUT",
                "exit": None}
    except OSError as exc:
        return {"found": True, "path": resolved,
                "version_line": "could not execute: %s" % exc, "exit": None}
    # Many tools print usage to stderr and exit non-zero on a version flag; that
    # is not a failure here -- we only want the banner text.
    return {"found": True, "path": resolved,
            "version_line": first_meaningful_line(out), "exit": rc}


def main():
    ap = argparse.ArgumentParser(
        description="Probe phylogenetics tools and record versions.")
    ap.add_argument("--require", default="",
                    help="comma list of tool names that must be present, e.g. "
                         "mafft,fasttree,Rscript.  Roles may also be used: "
                         "aligner,tree,domain_search,blast,r")
    ap.add_argument("--set", action="append", default=[], metavar="NAME=PATH",
                    help="override the executable used for NAME (repeatable)")
    ap.add_argument("--json", default=None, help="write the record to this JSON path")
    ap.add_argument("--timeout", type=int, default=25)
    args = ap.parse_args()

    overrides = {}
    for item in args.set:
        if "=" not in item:
            ap.error("--set expects NAME=PATH, got %r" % item)
        k, v = item.split("=", 1)
        overrides[k.strip()] = v.strip()

    report = {"tools": {}, "roles": {}, "missing_required": []}

    for name, argv in PROBES.items():
        resolved = overrides.get(name) or shutil.which(name)
        report["tools"][name] = probe_one(name, argv, resolved,
                                          timeout=args.timeout)

    # role availability
    for role, members in ROLES.items():
        present = [m for m in members if report["tools"].get(m, {}).get("found")]
        report["roles"][role] = present

    print("%-14s %-6s %s" % ("TOOL", "OK", "VERSION / NOTE"))
    print("-" * 78)
    for name in PROBES:
        info = report["tools"][name]
        mark = "yes" if info["found"] else "-"
        note = info["version_line"] or "not found"
        if info["found"] and info["path"] and info["path"] != name:
            note = "%s  [%s]" % (note, info["path"])
        print("%-14s %-6s %s" % (name, mark, note))

    print()
    print("role coverage")
    print("-" * 78)
    for role, members in report["roles"].items():
        print("  %-14s %s" % (role, ", ".join(members) if members else
                              "NONE AVAILABLE"))

    # The chain each marker class actually needs.  Spelled out because the
    # scripts live in separate files and it is easy to miss a step -- notably
    # that the six-frame output of frame_strand.py translate must be FILTERED by
    # `pick` before aligning, never aligned as-is.
    print()
    print("recommended call chain by marker class")
    print("-" * 78)
    for name, chain in [
        ("rrna / its (no reading frame)",
         "probe_tools -> define_taxa -> audit_inputs -> align (mafft) -> aln_qc -> "
         "infer (iqtree2 | fasttree) -> verify_tree -> root_tree.R"),
        ("protein-coding (needs frame + strand)",
         "probe_tools -> define_taxa -> audit_inputs -> frame_strand.py translate "
         "-> hmmsearch --cut_tc per frame -> frame_strand.py pick -> align -> "
         "aln_qc -> infer -> verify_tree -> root_tree.R"),
        ("fragmented set (diagnosed by aln_qc)",
         "reference backbone (hmmalign / papara / MAFFT --add --keeplength) "
         "[preferred]   else   blocks.py plan -> align each block -> "
         "blocks.py concat -> aln_qc -> infer"),
        ("short reads, or marker with a published reference tree",
         "fragment insertion (SEPP / EPA-ng) into the reference tree; declare that "
         "branch lengths come from insertion and are not comparable with de novo "
         "branch lengths"),
    ]:
        print("  %s" % name)
        print("      %s" % chain)
    print()
    print("  after changing frame_strand.py, re-run its contract self-test:")
    print("      python scripts/selftest.py")

    # resolve --require (names or roles)
    required = []
    for token in [t.strip() for t in args.require.split(",") if t.strip()]:
        if token in ROLES:
            if not report["roles"][token]:
                required.append("role:%s" % token)
        else:
            if not report["tools"].get(token, {}).get("found"):
                required.append(token)
    report["missing_required"] = required

    print()
    if required:
        print("MISSING REQUIRED: %s" % ", ".join(required))
        print("Nothing was installed.  Install the missing tool(s) and re-run, or "
              "pass --set NAME=PATH if it lives outside PATH.")
    else:
        print("All required tools present.")

    if args.json:
        d = os.path.dirname(os.path.abspath(args.json))
        if d:
            os.makedirs(d, exist_ok=True)
        with open(args.json, "w", encoding="utf-8", newline="\n") as fh:
            json.dump(report, fh, indent=2, ensure_ascii=False)
        print("record written: %s" % args.json)

    return 1 if required else 0


if __name__ == "__main__":
    sys.exit(main())
