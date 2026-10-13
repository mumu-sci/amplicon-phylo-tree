---
name: Bug report
about: A script errors out, or a documented command does not work
title: "[bug] "
labels: bug
---

<!--
Please fill in every field. Without the version and the tool log we usually cannot
reproduce the problem, and we would rather ask than guess.
-->

## Skill version

<!-- from SKILL.md frontmatter: metadata.version -->

## What you ran

```bash
# the exact command, verbatim, with every flag
```

## What happened

<!-- paste the FULL output, including the RESULT line and the exit code -->

```
```

## What you expected

## Input shape (no data needed!)

- number of sequences:
- marker / amplicon region (16S V4, ITS1, nifH, COI 313 bp, ...):
- alignment columns (if alignments are involved):
- abundance table, FASTA only, or an existing tree:

## Environment

- OS and version:
- Python version (`python --version`):
- MAFFT / FastTree / IQ-TREE / HMMER / BLAST+ / vsearch / R versions:

```
# python scripts/probe_tools.py --json tools.json   (attach tools.json)
```

## Checksums or logs

<!--
If the failure is in verification (verify_tree.py) or rooting (root_tree.R), attach
the report/log file. Do NOT attach your data if it is unpublished -- a --
synthetic reproduction is enough.
-->
