# CHANGELOG

All notable changes to this skill. Versions are for the skill folder
(`skills/amplicon-phylo-tree`); the repository layout, READMEs and this file are not
versioned separately.

## v1.0.1 — adversarial-review release

Two rounds of adversarial scientific review (design layer, then implementation
layer) found 10 blocking + 13 major + 9 minor defects in v1.0.0. This release fixes
all 10 blocking issues and the major ones that had machine-checkable consequences,
adds a 36-test black-box regression suite, and adds `check_consistency.py`, which
turns "the documentation promises X, the script does Y" into a failing test.

### The three "false passes" this release exists to eliminate

These are the reason v1.0.0 was not safe to follow as written. Each one made a run
**look verified when nothing had been verified**.

1. **Model assertions silently passed on Windows (B3).** PowerShell writes redirected
   logs as UTF-16LE. The verifier read them as UTF-8, so every wording pattern
   failed to match and the model check degraded to "no veto found" -> PASS. A log
   that actually contained the veto string
   (`WARNING! 100.0% NUCLEOTIDE CHARACTERS`) produced a green light.
   *Fix:* detect BOM/NUL before decoding; a supplied-but-unparsable log is now
   `UNVERIFIED` and fails, and `--require-tool-log` is on by default with
   `--allow-no-tool-log` as the explicit opt-out.

2. **A gate that could not run counted as a pass (B6).** For rRNA/ITS the input
   audit has no meaningful orientation or homology channel (there is nothing to
   translate), yet it printed `RESULT: PASS, exit 0`. The gate was decorative for
   exactly the markers that need it most.
   *Fix:* an unrun check is now a **failure**; `--accept-and-record` requires
   `--override-reason` and records the waiver; non-coding markers get a real channel
   (`--ref-db`) and a library-hit homology audit.

3. **`--require-binary` rejected every legitimate unrooted tree (B2).** An unrooted
   binary tree is read with a basal trifurcation (root degree 3), which was counted
   as a multifurcation. The acceptance command printed in the workflow therefore
   failed on every normal FastTree/IQ-TREE output, telling the user their correct
   tree was unresolved.
   *Fix:* the root is treated as rooted (degree 2) or unrooted (degree 3); only
   interior nodes must be bifurcating.

### Also fixed (blocking)

- **B1** rooted/unrooted files with the SAME topology were reported as stale: the
  split key used one side's hash. Now canonicalised; both sides key alike.
- **B4** the documented `--domtblout-dir` did not exist (`pick` only accepted a
  glob). Added; a missing directory is a clean exit 2.
- **B5** `pick` printed "unresolved sequences still enter the tree" and then dropped
  them (`if b is None: continue`). They are now kept via a reference-free ORF
  fallback, marked `orf-fallback` in the evidence table, and `--drop-unresolved`
  excludes them **loudly**.
- **B7** the ITS route printed "never build an ITS tree" and then printed the full
  MAFFT + IQ-TREE command; an agent executes what it is shown. The branch now prints
  no de novo command at all.
- **B8** the workflow listed a codon `%3` check that did not exist. Implemented,
  with a back-translation-evidence flag that fails when the amino-acid alignment
  contains `X` (only MACSE can map those).
- **B9** >100,000 tips were still routed to exact ML with full UFBoot. Now FastTree,
  with a loud scale warning; the residue clause can no longer override a huge N.
- **B10** the back-translation contract was self-contradictory (`frame_strand`
  records stops as `X`; PAL2NAL refuses CDS with stops; `trimAl -backtrans` cannot
  map `X`; only MACSE tolerates them). The contracts are now written down in D2b and
  in the workflow step.

### Major issues addressed

- **M1** the correction of a false measured fact: uppercase `BLOSUM` occurs **0**
  times in an IQ-TREE log, case-insensitive `Blosum` **36** times (it was recorded
  as "10"). The conclusion (never blacklist BLOSUM) stands and is stronger. Every
  word-list match must now state its case sensitivity.
- **M2** an amino-acid alignment made only of letters that are also IUPAC nucleotide
  codes was granted a data-side PASS. It is now `AMBIGUOUS`, and a waiver cannot
  resolve ambiguity without evidence.
- **M3** the promised mixed-alphabet hard stop did not exist and
  `--ambiguous-policy` was documented but unimplemented. Both implemented, plus the
  diagnostic matrix the design promised in `--json`.
- **M4** `--accept-and-record` left its trace only on stdout, although the whole
  justification of the override rests on traceability. It now lands in
  `gate_overrides` in the JSON, and requires `--override-reason`.
- **M5** fixed thresholds treated ITS and length-heterogeneous amplicons as if they
  were conserved rRNA. Defaults now follow the marker family, an explicit threshold
  still wins, and a ratio FAIL prints the copy-pasteable recorded-override command.
- **M6/M7** "you must give rooting-sensitivity evidence" had no enforcement point:
  the sensitivity test defaulted off and the gate still passed. It is now on by
  default, one file per candidate root is written, and
  `verify_tree.py --rooting-report` fails a report whose sensitivity section says
  `NOT RUN`. The design document no longer claims downstream metrics (UniFrac,
  betaNTI) are compared — the skill reports the root split only.
- **M8** the audit's temporary directory was a fixed name, never cleaned, and could
  land next to the input FASTA. It is removed first and defaults to the system temp.
- **M9** the mixed-region test only looked at length bimodality, so two regions of
  the same length slipped through. Domain-coverage grouping and composition grouping
  were added.
- **M10** `--region-map` was decorative: never read, never validated, never
  recorded. It is now parsed, checked against the input ids, and recorded with a
  SHA256; a waiver without a map counts as an unchecked item.
- **M11** a consistency self-check did not exist. `check_consistency.py` now asserts
  the D1–D11 five-field protocol, the scaling grid against the recommender's own
  grid, every `--flag` in the documentation against the scripts' argparse options,
  and the documented script paths.
- **M12** the grey zone ("build both routes and compare") had no usable definition.
  `p_aa` has a written definition (domain alignment, all pairs, gapped sites
  excluded, median **and** sum-based), the 85/90 numbers are labelled as convention
  anchors rather than a citable criterion, and the RF recipe states what must be
  equal before two trees may be compared (same tip set, same region; RF is
  topological; branch lengths only after back-translation).
- **M13** the `organelle` family lacked its interpretive ceiling row. Added, with
  the corrected COI direction (a 313 bp within-genus comparison usually goes
  nucleotide + MACSE, not amino-acid), and internal stop codons are now reported as
  the first numt/frameshift signal instead of being computed and ignored.

### Minor issues addressed

`blocks.py concat --gap-char` implemented (and an invalid filler refused);
`define_taxa.py` accepts the documented `nonzero>=N` presence form and reports an
unparseable rule instead of raising; the gap-fraction denominator is fixed at
`n_seq x width` and the two possible denominators are both reported for a ragged
alignment; a supplied-but-uninformative tool log is reported as "no positive
evidence" rather than "no log"; the no-table route states in the log and JSON that
no abundance matrix is involved; the two-dimensional table in `scaling.md` is now
actually two-dimensional and matches the recommender; and the alpha constant in the
memory heuristic is labelled unverified and contradicted by a new measurement.

### Verified in this release

```
regression_tests.py   36/36 PASS
selftest.py           PASS (coordinate contract, stops, genetic codes, real back-translation round trip)
check_consistency.py  PASS (5 check groups)
end-to-end smoke      baseline 20 -> 15 kept / 556 columns / 46.04% gaps, all gates PASS
```

### Language

`SKILL.md` and all six reference documents are in English. `SKILL.md`'s body contains
no Chinese at all; its `description` is English-first and keeps a single line of
Chinese trigger words (建树 / 系统发育树 / 进化树 / 16S 树 / ITS 树 / nifH 树 / ASV 树 /
比对建树 / 树挂了 / 空位太多 / 为什么树不对) so that Chinese-language requests still
match the skill. **`references/outputs.md` keeps its Chinese methods template in
Chinese on purpose** -- being able to paste a finished Chinese Methods paragraph is a
feature, so that template and its English counterpart are kept side by side and are
never merged.

## v1.0.0 — first release

Workflow document (Chinese) plus ten helper scripts, a decision protocol D1–D11, six
reference documents, and an `example-mini` dataset with a byte-reproducible
generator.
