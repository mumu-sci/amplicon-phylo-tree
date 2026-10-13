---
name: Wrong tree
about: The pipeline ran and passed, but the tree is wrong or indefensible
title: "[tree] "
labels: wrong-tree
---

<!--
This template exists for the case the gates are supposed to catch and did not: the
run reported PASS but the result does not mean what it appears to mean. That is the
most valuable report we can get, so please be concrete.

If you have not run the acceptance gate yet, do that first -- it may already tell
you which gate should have stopped this:

    python scripts/verify_tree.py --rooted tree.rooted.nwk --unrooted tree.unrooted.nwk \
        --ids kept.ids --tool-log qc/tree.log --seqtype nt --require-rooted \
        --require-binary --check-same-topology --rooting-report qc/rooting.txt
-->

## Skill version

<!-- from SKILL.md frontmatter: metadata.version -->

## What is wrong with the tree

<!-- e.g. a known non-target taxon is nested inside the marker; a clade that should
     be monophyletic is not; branch lengths are absurd; support values look copied
     from a previous run -->

## What the gates said

| Gate | Verdict from your run |
|---|---|
| G1 (frame / strand) | |
| G1b (input audit) | |
| G2 (alignment QC) | |
| G3 (model vs data type) | |
| G4 (tree and tip set) | |
| G5 (traceability) | |

## The numbers

- taxa in the taxon list / tips in the tree:
- alignment: sequences x columns, gap fraction:
- alignment algorithm actually used:
- tree tool and version, model, support values:
- rooting method (outgroup / midpoint), and whether the sensitivity test ran:

## Commands and logs

```bash
# the alignment command, the tree command, and the rooting command, verbatim
```

<!-- attach qc/aln.json, qc/verify.json, qc/rooting.txt, qc/tools.json if you can -->

## What you expected instead, and why

<!-- a reference tree or a known relationship is enough; no unpublished data needed -->
