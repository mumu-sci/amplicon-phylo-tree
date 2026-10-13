---
name: amplicon-phylo-tree
description: Build defensible, publishable phylogenetic trees from second-generation amplicon data (16S/18S, ITS, nifH/amoA/COI and other marker genes): parameters are fixed per study (no hard-wired thresholds), then input audit, alignment, tree inference, rooting and acceptance, producing the tree, the alignment, the QC record and bilingual (Chinese and English) methods paragraphs. Use when the user asks to build a phylogenetic tree from amplicon data, run MAFFT / FastTree / IQ-TREE, make a phylo tree from ASVs or OTUs, check or root an existing tree, or diagnose a broken one. Chinese trigger words: 建树, 系统发育树, 进化树, 16S 树, ITS 树, nifH 树, ASV 树, 比对建树, 树挂了, 空位太多, 为什么树不对. Hard gates: the model MUST match the data type (protein vs nucleotide) as decided independently from the alignment's own character set; the tip set MUST equal the input ids exactly; a structural error or a model mismatch is a HARD STOP, and you MUST NOT deliver a tree with a known defect.
compatibility: Requires MAFFT and FastTree or IQ-TREE; HMMER and/or BLAST+ for the input audit of protein-coding markers; R 4.x with ape and phangorn; Python 3.9+. Reference databases must be supplied by the user -- this skill neither downloads nor installs anything.
license: MIT
metadata:
  author: Mumu (github.com/mumu-sci)
  version: 1.0.1
  category: bioinformatics
  tags: [phylogenetics, amplicon, 16S, ITS, nifH, ecology, MAFFT, FastTree, IQ-TREE]
---

# Phylogenetic trees from second-generation amplicons

Turn "an abundance table plus representative sequences" into "a tree fit for a
paper" -- and **write down every decision** you made to get there.

This skill applies to Illumina amplicon data (16S / 18S / 25-28S, ITS,
nifH / amoA / nirK / pmoA / nosZ / COI / rbcL and other functional genes, and any
custom amplicon). It **presets no thresholds**: every parameter is fixed for the
study at hand, with its rationale and its deviation conditions written in
`references/parameters.md`.

---

## Iron rules (read first; breaking one is a failure)

1. **Parameters are decided per study and written down; implicit defaults are
   forbidden.** Each of D1-D11 (`references/parameters.md`) MUST state a **value plus
   its rationale**. Using a default is allowed, but the parameter table must say that
   it is the default. "Just run it the usual way" is not allowed.

2. **The model MUST match the data type, and it is INJECTED rather than asserted
   afterwards.**
   The data type is decided independently from the alignment's own character set,
   and the tree command is **generated** from that decision (nucleotide
   `-nt -gtr` / protein `-lg -gamma` / IQ-TREE `-st DNA|AA`).
   The most common real accident: a protein model (`-wag`/`-lg`/BLOSUM) run on pure
   ACGT data. The tool computes happily and prints a plausible likelihood.
   See G3 and `pitfalls.md`.

3. **An alignment is examined before it is used; a broken one MUST NOT be carried
   forward.** Structural errors (not rectangular, mixed alphabets, codon mode whose
   width is not a multiple of three) are a **HARD STOP**;
   ratio metrics (width ratio, gap fraction) default to FAIL and can be overridden
   with `--accept-and-record`, but the override MUST be written into the ops log.

4. **The taxon rule MUST be explicit, reproducible, and the author's decision.**
   Write it as a formula over the abundance table (e.g. `sum>=10 and prev>=2`) and
   record it as an id list.
   **When the rule changes, alignment -> tree -> rooting -> verification MUST all be
   re-run**; subsetting an existing tree and reusing it is FORBIDDEN (it silently
   changes branch lengths and the root position).
   This skill **never deletes taxa automatically** and never changes the taxon set
   silently.

5. **Nothing ships before it passes the acceptance gate.** G4 MUST pass: the tip set
   equals the input ids in both directions, the Newick is complete, and the tree is
   rooted as its purpose requires.
   A run manifest (full commands + tool versions + artefact SHA256) is produced
   alongside.

---

## When it applies / when it does not

**Applies**: you have an ASV/OTU abundance table plus representative sequences (or
only the representative-sequence FASTA) and want a de novo tree, or you want to
check or root an existing tree.

**Does not apply** (the skill says so explicitly and offers an alternative route):
specialised handling of long-read third-generation data (PacBio / Nanopore);
**executing** upstream denoising (DADA2 and friends) or OTU clustering; downstream
ecological analysis itself (computing PD / UniFrac / NRI-NTI / betaNTI).

**Input assumptions (MUST be confirmed with the user, never assumed)**: sequences
have already been **denoised and chimera-checked**, and their **orientation has been
normalised**. If that is uncertain, this skill's **G1b input audit** checks
orientation and homology for you; chimera detection needs `vsearch`, and the skill
prints the command or requires you to record an explicit skip.

---

## Three questions first (they decide every route that follows)

1. **What is the tree FOR?** -- this decides whether the taxon set may be narrowed,
   and it matters **far more than how large N is**.
   - **Distance skeleton** (UniFrac / Faith's PD / NRI-NTI / betaNTI): **every taxon
     in the analysis matrix MUST be kept** (the distance matrix is defined over all
     of them; dropping taxa redefines the metric's domain); but the **highest
     precision alignment is not needed** -- spend the effort on **the root + region
     homology + the reference-database version + tip consistency**.
   - **The conclusion itself** (locating a novel clade, functional-gene diversity
     patterns, within-clade relationships): **narrowing first is right**, but you
     MUST report the retained reads% and richness%, and run one sensitivity analysis
     over the narrowing rule.
   - **Never allowed**: build a full tree, then subset it and compute metrics from
     the subset.
2. **Which marker family?** (`rrna` / `its` / `functional-nuclear` / `organelle` /
   `other`) -- this decides the data type, the model family, whether a reading frame
   is needed, and the reference database.
3. **What scale?** N x alignment length -- this decides the alignment algorithm and
   the tree software (see the two-dimensional table in `references/scaling.md`).

**Triage question**: is your downstream a metric such as UniFrac/PD, or a de novo
tree you will publish? If it is the former and the marker is ITS, or the sequences
are short (V4 ~250 bp, COI 313 bp), **fragment insertion (SEPP / EPA-ng) or "no tree
at all" is often the more correct choice** -- see "Four routes".

---

## Four routes

| Route | When to use it | Output |
|---|---|---|
| **Fast route** (default) | ordinary cases | de novo alignment + tree |
| **Reference-backbone route** | a marker reference alignment exists; or the **first-choice remedy** after a fragmented alignment | align to the reference alignment (`hmmalign` / `papara` / `MAFFT --add --keeplength`), preserving comparability and codon continuity |
| **Insertion route** (routing + documentation) | short reads, or a marker with a mature reference tree (16S/18S/COI), or ITS whose downstream needs UniFrac/PD | SEPP / EPA-ng insertion into a reference tree; the methods **MUST state that branch lengths come from insertion and are not directly comparable with de novo branch lengths in other studies** |
| **No-tree route** (routing + documentation) | you only need names (UNITE dynamic-clustering match), or you only need the ASV/OTU table | **state explicitly that not building a tree is legitimate and often more correct** |

When the alignment is badly fragmented, priority is **reference backbone > block
concatenation**. Block concatenation (`scripts/blocks.py`) creates false signal at
the block boundaries and, for protein-coding markers, **breaks codon continuity**, so
it is a last resort; if you use it you MUST state that **branch lengths across blocks
carry no biological meaning**.

---

## Workflow

### ① Probe the environment
```bash
python scripts/probe_tools.py --require mafft,tree --set fasttree=/opt/FastTree --json qc/tools.json
```
A missing required tool exits 1 and says which one. **Versions MUST be recorded** --
different FastTree versions behave differently by default.

### ② Fill in the parameter table (D1-D11)
Work through `references/parameters.md` and write `params.md`. Three items are not
optional:
- **D1 region**: `rrna` MUST give a "region <-> reference database" pair; **different
  regions or primer pairs MUST NEVER be combined into one tree**.
- **D2b data-type criterion** (protein-coding markers): within-domain pairwise
  amino-acid identity `p_aa`
  < 85% -> **amino acids are mandatory, and the alignment MUST then be back-translated
  into a codon alignment** (PAL2NAL / `trimAl -backtrans` / MACSE);
  85-90% -> **build BOTH**, report topological agreement (RF / conflicting branches),
  then choose one;
  > 90% -> a nucleotide tree is acceptable, but you MUST state that the third codon
  position is unsaturated.
  The thresholds are restatements of a domain convention, not a citable criterion;
  they can be overridden per study, but **"build both" in the grey zone is MANDATORY**.
- **D9 rooting decision rule**: an outgroup or reference sister lineage available ->
  **the outgroup is mandatory**; the purpose is a community metric or reference
  placement -> outgroup or reference tree, otherwise **rooting-sensitivity evidence
  is mandatory**; only "single marker + homogeneous taxa + no outgroup + display
  only" permits the midpoint method.

### ③ Define the taxon rule (skippable)
```bash
python scripts/define_taxa.py --table asv_table.tsv --rule 'sum>=10 and prev>=2' \
       --fasta rep_set.fasta --out-ids kept.ids --out-fasta kept.fasta --log qc/taxa.log
```
Rule grammar: `sum` (total abundance) / `prev` (number of samples present in;
`prev>=1` is "presence") / `nonzero` (same as `prev`) / `sum_frac` (percent of the
whole library's abundance) / `prev_frac` (percent of samples present in), combined
with `and` / `or` (`or` binds loosest).
With only a FASTA use `--fasta-only rep_set.fasta --out-ids kept.ids` (and then you
**MUST NOT claim alignment with any abundance table**).
The script detects "this is actually relative abundance, not reads" and refuses to
run blind. **Once the rule changes, re-run from ④ to the end**.

### ④ Input audit (G1b)
```bash
python scripts/audit_inputs.py --fasta kept.fasta --hmm marker.hmm --outdir qc --json qc/audit.json
```
Reports: reverse-complement hits, the fraction with no significant hit, length
bimodality (= two amplicon regions mixed), and a chimera-step reminder.
Reverse complements above the threshold -> a fix list and **re-run the whole
pipeline**. A mixed region detected -> **stop and report** (`--mixed-region-ok
--region-map` is the explicit waiver, and the map is read, validated and hashed).
**A waiver without a map fails**: `--accept-and-record --override-reason` is
required to record a map-less waiver, because an untraceable waiver is an unchecked
gate, not a pass. The same applies to the second mixed-region criterion (domain
coverage / composition grouping): if that channel cannot run, it is `NOT CHECKED`
and therefore a failure.

### ⑤ Reading frame and strand normalisation (**protein-coding markers only**)
```bash
python scripts/frame_strand.py translate --fasta kept.fasta --outdir frames --genetic-code 1
for f in frames/frame_*.faa; do b=$(basename $f .faa); \
  hmmsearch --cut_tc --noali --domtblout frames/$b.domtblout marker.hmm $f > /dev/null; done
python scripts/frame_strand.py pick --fasta kept.fasta --domtblout-dir frames --mode hmm \
       --out-nt kept.frame.fasta --out-aa kept.faa --out-tsv qc/frames.tsv
```
**Genetic code**: COI and other mitochondrial coding genes MUST use
`--genetic-code 5` (invertebrates) or `2` (vertebrates); the default table treats TGA
as a stop, so **the wrong frame gets picked**. 12S is an rRNA gene and is **not
translated**.
With no reference HMM use `--mode orf` (the frame is still normalised; the strand call
is weaker).
The six-frame protein FASTA that `translate` writes **contains non-coding frame
junk and MUST NOT be aligned as-is** -- `pick` first.
After changing `frame_strand.py`, run `python scripts/selftest.py` (it asserts the
aa<->nt coordinate contract, the genetic code tables, and a real back-translation
round trip).

**Back-translation tools are NOT interchangeable** (on the amino-acid route, see
`references/parameters.md` D2b):
`translate` **does not truncate** and records a stop as `X`, therefore --
**PAL2NAL** requires in-frame CDS **without stop codons** (it errors if it finds
one); **`trimAl -backtrans`** **cannot map `X`**;
only **MACSE** is codon-aware and **tolerates stops and frameshifts** (marking them
with `!`).
-> **when the aa alignment contains `X`, only MACSE can be used**; otherwise remove
the sequences with internal stops under D4 (with the counts recorded) before
back-translating.
`verify_tree.py --backtranslation-evidence ... --aa-aln ...` asserts this (see ⑨/⑪).

### ⑥ Alignment
Choose the algorithm from the two-dimensional table in `references/scaling.md`
(**not by N alone**):
```bash
mafft --localpair --maxiterate 1000 --thread 8 kept.frame.fasta > aln.fasta   # small, high accuracy
mafft --auto --quiet --thread 8 kept.frame.fasta > aln.fasta                  # medium
mafft --parttree --retree 2 --quiet --thread 8 kept.frame.fasta > aln.fasta   # very large (declare the quality cost)
```
> ⚠ `--thread -1` is the form MAFFT actually documents. ⚠ Alignment output is wrapped
> at 60 characters, so every length statistic must reassemble the records first.

### ⑦ Alignment QC (G2)
```bash
python scripts/aln_qc.py --aln aln.fasta --fasta-raw kept.frame.fasta --label <study> \
       --codon --marker-class functional-nuclear --ambiguous-policy gap --json qc/aln.json
```
The gap convention is fixed as **gap characters / (N_seq x N_col)**, and the report
prints that definition.
**Structural errors are a HARD STOP and cannot be overridden**: not rectangular;
codon mode whose width is not a multiple of 3; **illegal characters**; and **mixed
alphabets** (non-gap symbols contain both nucleotide-only `ACGTU` and protein-only
`EFILPQZJO*` -- no single substitution model can describe both).
**Ratio metrics** (width ratio, gap fraction) default to FAIL and may be overridden
with `--accept-and-record --override-reason "..."`; **the override is written into
`qc/aln.json`** under `gate_overrides` (gate / metric / value / limit / reason), not
only printed. The FAIL report prints the copy-pasteable override command.
The default gap threshold follows the **marker family** (`--marker-class`):
0.80 for rrna/organelle, 0.85 for functional-nuclear, 0.90 for its/other; an explicit
`--max-gap*` always wins and the report says which is in force.
The **ambiguity policy** `--ambiguous-policy {gap,missing,drop-col}` is D5:
`drop-col` **only writes a list of suggested columns** and never edits the alignment
file (columns above `--ambiguous-max-frac` ambiguous, default 0.10).
`qc/aln.json` also carries the **diagnostic matrix**: fraction of gapped columns,
pairwise gap quantiles, columns before/after trimming terminal gaps, ATGC fraction,
IUPAC/N counts and affected columns, and the codon third-position difference rate in
codon mode.
On failure, work the diagnosis (non-homologous contamination / two regions /
inconsistent reading frames / ITS length heterogeneity), then take the **reference
backbone** or `blocks.py`.

### ⑧ Tree inference
Choose the tool on axis C; **use IQ-TREE if you can wait**, otherwise FastTree and
state that it is approximate ML.
```bash
iqtree2 -s aln.fasta -m MFP -B 1000 -alrt 1000 -T 8 --prefix tree -redo -pre <run-id>
fasttree -nt -gtr -gamma aln.fasta > tree.nwk 2> qc/tree.log      # nucleotide
fasttree -lg  -gamma aln.faa   > tree.nwk 2> qc/tree.log          # protein
```
IQ-TREE MUST be given `-pre`, otherwise the default output names **silently
overwrite** the previous run.

### ⑨ Model-consistency assertion (G3)
```bash
python scripts/verify_tree.py --aln aln.fasta --seqtype nt --tool iqtree --model GTR+F+I+G4 \
       --tool-log qc/tree.log --iqtree-report tree.iqtree --ids kept.ids \
       --unrooted tree.nwk --json qc/verify.json
```
Four layers of evidence, **fail-closed**: data-side character-set call (authoritative)
-> whether the invoked parameters agree with that call -> tool-specific measured
evidence -> no positive evidence anywhere means `UNVERIFIED` and a stop.
Do **not** use bare `BLOSUM` as a blacklist: IQ-TREE legitimately prints it while
choosing a protein model (measured: uppercase `BLOSUM` **0** times, case-insensitive
`Blosum` **36** times, all in the candidate model table), so it misfires.
**Any word-list match MUST state whether it is case-sensitive.**
FastTree's wording (`Nucleotide distances` and friends) **MUST NOT** be required of
IQ-TREE (measured 0 occurrences).

### ⑩ Rooting
```bash
Rscript scripts/root_tree.R --in tree.nwk --rooted tree.rooted.nwk --unrooted tree.unrooted.nwk \
        --report qc/rooting.txt --label <study> [--outgroup outgroup.txt] [--sensitivity-candidates 3]
```
**Only `tree.rooted.nwk` goes downstream**; `tree.unrooted.nwk` is for inspection and
figures (support values exist only on the unrooted tree).
The midpoint method assumes roughly clock-like rates and is indefensible when the
tree spans widely divergent clades; the report always contains a **rooting
sensitivity** section.
**The sensitivity test is ON by default** (`--sensitivity-candidates 2`;
`--no-sensitivity` turns it off and the ops log MUST say why): with no outgroup it is
the only evidence about whether the root is decided by the data or by the method.
Every compared candidate root is also written as `<rooted>.candN.nwk` (the first line
is a comment naming the candidate and how the root was placed) so you can re-check
the claim or compute downstream metrics yourself.
This skill **reports the root-split sensitivity only, and does NOT compute UniFrac /
betaNTI for you** -- do that on the same candidate trees.

### ⑪ Acceptance (G4)
```bash
python scripts/verify_tree.py --rooted tree.rooted.nwk --unrooted tree.unrooted.nwk \
       --ids kept.ids --tool-log qc/tree.log --seqtype nt --require-rooted --require-binary \
       --check-same-topology --rooting-report qc/rooting.txt --json qc/verify.json
```
Checks: the tip-set difference from the input ids is 0 in both directions; the Newick
ends with `;` and parses; tips are unique; root degree; branch-length range;
**zero-length branch fraction** (reported when >0); negative branch lengths (FAIL);
number of support labels; multifurcations; duplicate sequences; codon width %3;
the rooted and unrooted files **share the same unrooted topology** (this catches "the
file was overwritten by the previous run"); and
**with `--rooting-report` and `--require-rooted`, that the report's sensitivity
section actually ran** -- a report saying `NOT RUN` is a FAIL (you MUST NOT deliver an
untested root when there is no outgroup).

### ⑫ Delivery
Tree + unrooted tree + alignment + QC files + a **bilingual Chinese/English methods
paragraph** (the templates in `references/outputs.md`, with slots for the
limitations statement and the **support-value declaration**) + run manifest (commands,
versions, SHA256) + ops log (including every deviation and override and its reason).

---

## Gate summary

| Gate | What it covers | On failure |
|---|---|---|
| **G1** | reading-frame/strand resolvable fraction (coding markers only) | **reported only**; unresolvable sequences still enter the tree; the author decides |
| **G1b** | input audit: orientation, homology, mixed regions, chimera step | orientation mixture -> fix and re-run; insufficient homology or mixed regions -> **stop and report** |
| **G2** | alignment QC | structural errors **HARD STOP**; ratio metrics default to FAIL, overridable with `--accept-and-record` and recorded |
| **G3** | model vs data type (four layers + fail-closed) | **HARD STOP** |
| **G4** | tree structure and tip set | **HARD STOP** |
| **G5** | traceability (run manifest + SHA256) | report and stop |

---

## Common pitfalls (full table in `references/pitfalls.md`)

| Symptom | Root cause | Action |
|---|---|---|
| The likelihood looks fine but the tree is garbage | a protein model on nucleotide data (`-wag`/`-lg`/BLOSUM) | G3's four layers; `WARNING! 100.0% NUCLEOTIDE CHARACTERS` in the tool log is a veto |
| Thousands of sequences, 90%+ gaps | non-homologous sequences mixed in / two amplicon regions combined / inconsistent reading frames | G1b audit + G2 QC + `blocks.py`; prefer the **reference backbone** |
| Back-translation dies at PAL2NAL / trimAl | the aa alignment contains `X`, or the CDS contains a stop (D2b) | switch to **MACSE**, or remove the internal-stop sequences under D4 with the counts recorded |
| The column count is wrong | MAFFT wraps at 60 characters and the line length is taken as the column count | every script reassembles records first |
| The tree has fewer tips than the table | a downstream step silently dropped taxa with no tree tip | G4's two-way difference check; `--table` compares directly against the table |
| Support values come from a previous tree | IQ-TREE's default output names were overwritten | force `-pre`; G4's `--check-same-topology` also catches it |
| Support values are not comparable after switching tools | SH-like / UFBoot / SH-aLRT mean different things | re-run the **whole tree**; cross-tool thresholds are FORBIDDEN; use the declaration template in `outputs.md` |
| Nobody ever tested the root position | with no outgroup the sensitivity test was skipped and the report still passed | `root_tree.R` runs `--sensitivity-candidates 2` by default; ⑪ asserts it with `--rooting-report` (M6) |
| COI is wrong from the very first step | the standard genetic code was used (TGA read as a stop) | `--genetic-code 5`; `selftest.py` locks the distinguishing codons of every table |

---

## Smoke test (run this first on a new machine)

`assets/example-mini/` holds 20 synthetic sequences and a small abundance table
(deliberately containing four traps: a `#OTU ID` header, `X.0` floats, a trailing
taxonomy column, 60-character wrapping).

```bash
P=assets/example-mini
python scripts/probe_tools.py --require mafft,tree
python scripts/selftest.py
python scripts/define_taxa.py --table $P/asv_table.tsv --rule 'sum>=10' \
       --fasta $P/rep_set.fasta --out-ids /tmp/kept.ids --out-fasta /tmp/kept.fasta
mafft --auto --quiet /tmp/kept.fasta > /tmp/aln.fasta
python scripts/aln_qc.py --aln /tmp/aln.fasta --fasta-raw /tmp/kept.fasta
fasttree -nt -gtr -gamma /tmp/aln.fasta > /tmp/tree.nwk 2> /tmp/tree.log
Rscript scripts/root_tree.R --in /tmp/tree.nwk --rooted /tmp/tree.rooted.nwk \
        --unrooted /tmp/tree.unrooted.nwk --report /tmp/rooting.txt --label smoke
python scripts/verify_tree.py --rooted /tmp/tree.rooted.nwk --unrooted /tmp/tree.unrooted.nwk \
        --aln /tmp/aln.fasta --seqtype nt --tool fasttree --tool-log /tmp/tree.log \
        --ids /tmp/kept.ids --require-rooted --rooting-report /tmp/rooting.txt
```

---

## Reference files

| File | When to read it |
|---|---|
| `references/parameters.md` | **Required reading.** The D1-D11 decision protocol (five fields each) plus a fill-in parameter-table template |
| `references/marker-types.md` | when fixing the marker family. Per-family model family, genetic code, reference-database channels, and **interpretive ceiling** |
| `references/scaling.md` | when choosing the alignment algorithm / tree software / compute. A two-dimensional table + measured anchors + a resource heuristic + **support-value declaration templates** |
| `references/inputs.md` | when the input format is non-standard. Conversion commands for QIIME2/BIOM/mothur/phyloseq/FASTA-only |
| `references/pitfalls.md` | when something goes wrong. The full symptom / root cause / action table |
| `references/outputs.md` | at delivery time. Run manifest, **bilingual Chinese/English methods templates**, ops-log template |
