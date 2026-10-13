# Parameter Decision Protocol D1–D11

This file is the carrier of "thresholds are not hard-wired". Every item MUST state the **five fields** clearly: What to decide / How to check / Default / When to deviate / Where to record it.
Copy this file into your research project as `params.md` and fill in every item -- **fill it in before you start**.

---

## D1 Marker gene and amplicon region

- **What to decide**: marker name + amplicon region + expected length range.
- **How to check**: ask the user; or infer it from the representative-sequence length distribution (make a length histogram first).
- **Default**: No default -- it MUST be stated explicitly.
- **When to deviate**: if the length distribution is **bimodal** (peak separation > 100 bp and each peak > 5%), first use G1b to decide whether two amplicon regions were mixed; **mixed amplicon regions MUST NOT be built into the same tree**.
- **Where to record it**: `params.md`; state the amplicon region and primers in the Methods.

| Amplicon region | Expected length (rough) | Suggested reference-database pairing |
|---|---|---|
| 16S V1–V3 | ~450–600 bp | SILVA / Greengenes |
| 16S V3–V4 | ~420–470 bp | SILVA |
| 16S V3–V5 | ~600 bp | SILVA |
| 16S V4 | ~250 bp | SILVA (**for short fragments, consider fragment insertion first**) |
| 16S V4–V5 | ~300 bp | SILVA |
| 16S full length | ~1,500 bp | SILVA |
| 18S V4 / V9 | ~380 / ~130 bp | PR2 (**incompatible with the SILVA taxonomy -- pick one and use it to the end**) |
| Fungal ITS1 / ITS2 | 150–350 / 200–450 bp | UNITE (**match for naming first, not for tree inference**) |
| nifH (Zehr/PolF region) | ~230–460 bp | FunGene / ARB nifH; HMM = Fer4_NifH (PF00142) |
| amoA (bacterial / archaeal) | ~490 / ~630 bp | FunGene; the two types **MUST NOT be built together** |
| COI (Leray 313 bp) | ~313 bp | BOLD (**no HMM**; use a reference backbone or `--mode orf`); **genetic code table 5** |

## D2 Marker gene class

- **What to decide**: one of `rrna` / `its` / `functional-nuclear` / `organelle` / `other`.
- **How to check**: judge it from the marker name. **Note that 12S is a mitochondrial rRNA gene, behaves as `rrna`, and is not translated**; COI/rbcL/matK are `organelle`.
- **Default**: No default -- it MUST be stated explicitly.
- **When to deviate**: several markers in one study → a separate parameter set and a separate tree for each marker.
- **Where to record it**: `params.md`.

## D2b Data type: nucleotide / amino acid / codon-aware (**MANDATORY for protein-coding marker genes**)

- **What to decide**: whether the final tree is built on nucleotides or amino acids (if amino acids, whether to back-translate to a codon alignment).
- **How to check**: compute `p_aa` using the **hard-wired rule** below (otherwise the number cannot be reproduced, nor compared with others):

  **Definition of `p_aa` (the only rule)**
  1. **Which segment**: inside the marker HMM's **domain alignment** (`hmmalign <hmm> <fasta> > aln.sto`, take the columns marked as
     match by `#=GC RF`; **excluding** insert columns). Not computed on the raw full-length sequences, and not on the MAFFT full-length alignment.
  2. **How to pair**: **all sequences paired two by two** (n(n−1)/2 pairs), computed pair by pair.
  3. **Gapped sites**: **count only sites where neither sequence is a gap**; if the pair has **nothing comparable**
     at those sites, skip that pair, and write in the record how many pairs were skipped (too many skipped pairs means the domain alignment itself is problematic -- look at G1b first).
  4. **Pooled or mean?** Report the **median**, and report the mean as well. Method:
     `p_aa = 1 − (sum of identical sites over all pairs / sum of comparable sites over all pairs)`.
     The thresholds 85/90 are **anchor examples** (see the "source" section of D2b), whereas **a median–mean difference > 5 percentage points**
     is itself a hint that "the dataset spans functional clusters / deep lineages" → follow the most conservative side (i.e. go with amino acids).
  5. **Where to record it**: write the **pooled result, the median, the mean and the number of skipped pairs** -- all four numbers -- into the D2b row of `params.md`.

- **Default** (field convention):

  | `p_aa` | Route |
  |---|---|
  | < 85% | **amino acids MANDATORY**; after alignment the alignment **MUST be back-translated into a codon alignment** (PAL2NAL / `trimAl -backtrans` / MACSE) |
  | 85–90% | **build both**, report topological concordance (RF / conflicting branches), then choose one as the main tree (**MANDATORY**) |
  | > 90% | a nucleotide tree is acceptable, but you MUST declare "nucleotide identity >90%, third position not saturated" |

  **⚠ The two numbers 85/90 are a "convention-based starting point / anchor example", not a criterion with a single authoritative source.**
  They come from common practice in functional-gene phylogenetics (distant taxa go amino-acid, close relatives may go nucleotide); **no** paper defines
  "85% is the boundary". Therefore:
  - they **MUST NOT** be the sole basis; across functional clusters (see below) always go amino-acid, regardless of the numbers;
  - overriding them for a study is allowed, but the reason for the override MUST be written in `params.md`;
  - if a reviewer challenges it, the alternative criterion you can compute yourself is **third-codon-position saturation** (reported by `aln_qc.py --codon` as
    `third_position_difference_rate`: the third-position pairwise difference rate approaching the saturation value 0.75 → the third position is uninformative, go amino-acid).
  This item is this skill's exception to its general stance on "hard-wired thresholds": **the numbers give a starting point; the criterion rests on a quantity you can compute yourself.**

  Also judge by whether the taxa span functional clusters: across nifH Cluster I–IV, bacterial vs archaeal amoA, nirK vs nirS → always go amino-acid.

  **The rule for the grey zone "build both, report topological concordance" (you MUST follow it, otherwise the RF comparison itself is a methodological error)**
  - **Same tip set**: the two trees MUST come from the **same `kept.ids`** (not one more, not one less); different tips →
    resolve the rule problem first; comparing RF is **not allowed**.
  - **Same region**: the same amplicon region/primer pair, the same reference-database version.
  - **RF compares topology only**: `phangorn::RF.dist` by default **ignores branch lengths** (it compares the split sets), which is exactly what is wanted for
    "are the topologies concordant"; **branch lengths are not comparable** -- a nucleotide tree and an amino-acid back-translated codon tree have different branch-length dimensions, and **to compare branch lengths you MUST first back-translate to
    the same codon alignment**, otherwise you are comparing numbers on two different scales.
  - **Reproducible recipe** (R, both routes, run after building both):
    ```r
    library(ape); library(phangorn)
    a <- read.tree("tree.nt.nwk")          # the nucleotide tree
    b <- read.tree("tree.aa.codon.nwk")    # the back-translated codon tree, amino-acid route
    stopifnot(setequal(a$tip.label, b$tip.label))   # SAME tip set, or stop here
    rf  <- phangorn::RF.dist(a, b, normalize = TRUE)  # 0 = topologies identical
    # conflicting branches: compare support only on splits present in BOTH trees;
    # never apply a cross-tool threshold
    a_u <- ape::unroot(a); b_u <- ape::unroot(b)
    cat("normalised RF =", rf, "\n")
    ```
  - **Decision**: `normalize = TRUE` RF ≤ 0.10 → the two topologies are essentially concordant, and you may choose one as the main tree on the basis of information content / defensibility;
    otherwise **the main tree MUST be the most conservative one (amino acids + back-translation)**, and the Methods must report the RF value and the number of conflicting branches.
  - **Where to record it**: the "build both" row of D2b in `params.md` + the Methods (including the RF value, which one was used as the main tree, and why).

  **Back-translation tool contract (B10: check the input requirements first, then choose the tool, otherwise the back-translation step breaks outright)**

  This skill's `frame_strand.py translate` **does not truncate** and records stop codons as `X`; this is deliberate
  (truncation would silently change the aa↔nt coordinates, see the contract locked by `selftest.py`). Therefore the input requirements of the three back-translation tools
  **are not equivalent**, and you MUST choose according to the table below:

  | Tool | Input contract | Incompatibility with this skill's outputs |
  |---|---|---|
  | **PAL2NAL** | in-frame CDS, and it **must not contain stop codons**; it errors out on encountering a stop | sequences with internal stops / stops recorded as `X` will fail outright |
  | **`trimAl -backtrans`** | aa alignment + corresponding nt; maps codons position by position | **cannot map `X`** (`X` has no corresponding codon) |
  | **MACSE** | codon-aware alignment | **tolerates** stops and frameshifts, marks them with `!`, and the product is a codon alignment |

  **Rules (linked to the assertion of `verify_tree.py --aa-aln`; already enforced on the code side)**:
  1. if the aa alignment contains `X` → **only MACSE can be used** (both PAL2NAL and trimAl will fail or drop sites);
  2. when not using MACSE, first record under D4 the **removal of sequences containing internal stops**, stating the count and the abundance fraction, and only then back-translate;
  3. whichever route you take, it MUST be written into `params.md` and the operation log -- the back-translation tool determines **which sites exist** in the codon alignment,
     and that is a conclusion-level choice.
  The corresponding code-side assertion: `verify_tree.py --backtranslation-evidence <codon.fasta> --aa-aln <aa.fasta>`
  **FAILs** when the aa contains `X`, and reports that "only MACSE can handle this".

- **When to deviate**: if the goal is codon-level selection analysis (dN/dS) → nucleotides are MANDATORY and a codon model must be used; if the sequences are < 300 bp and the taxa are homogeneous → amino acids carry too little information and nucleotides are better.
- **Where to record it**: `params.md` + the Methods (state which criterion set was used and the result).

## D3 Taxon rule

- **What to decide**: an **explicit formula** applied to the abundance table.
- **How to check**: look at the abundance distribution and at the distribution of the number of samples of occurrence first, then choose a rule; trial-compute with `define_taxa.py` and look at the retained reads / richness fraction.
- **Default**: **No default** -- it MUST be stated explicitly (this is the step most often left unrecorded in the whole workflow, and the one that most affects conclusions).
- **When to deviate**: with FASTA only, use `--fasta-only`, in which case you **MUST NOT claim alignment with any abundance table**.
- **Where to record it**: `params.md`; the log of `define_taxa.py` (verbatim rule + retention fraction + the three-part discarded set); the Methods.
- **Hard constraint**: after the rule changes, **alignment → tree inference → rooting → acceptance MUST all be rerun**; pruning an existing tree and reusing it is FORBIDDEN.

## D4 Whether/how to narrow the taxa

- **What to decide**: whether to narrow by phylogenetic or taxonomic depth; whether to remove known paralogs (e.g. a homolog-cluster map).
- **How to check**: see "ask these three things first" in SKILL.md -- **look at the tree's role first**: a distance backbone MUST retain the full set; the conclusion itself should be based on a set narrowed and validated first.
- **Default**: **no automatic narrowing** (this skill does not delete taxa).
- **When to deviate**: if you remove paralogs using an external cluster map, you MUST state that map's source and reliability, and **record the count and the abundance fraction removed**.
- **Where to record it**: `params.md`; the Methods MUST report the retained reads fraction and the richness fraction.

## D5 Alignment strategy + escalation ladder + ambiguous-base policy

- **What to decide**: the alignment algorithm; which level to escalate to when G2 is not passed; how ambiguous bases/`N` are handled.
- **How to check**: use the **N × L two-dimensional table** in `scaling.md`; ambiguous bases per the table below.
- **Default**: axis C two-dimensional table; `--ambiguous-policy gap`.
- **When to deviate**: long unalignable regions → E-INS-i; reference alignment available → **reference backbone first**; high `N` fraction → `drop-col` and report the number of affected columns.
- **Where to record it**: `params.md`; commands written into the run manifest.

| `--ambiguous-policy` | What it does | When to use |
|---|---|---|
| `gap` (default) | ambiguous characters treated as gaps | general; **MUST be declared in the Methods** |
| `missing` | treated as missing data | more rigorous when the tool supports it |
| `drop-col` | drop columns whose content of ambiguous characters exceeds the threshold | high `N` fraction, or after reference sequences have been merged in |

## D6 Reference backbone / block alignment

- **What to decide**: whether alignment to a reference alignment is possible; if not, whether to split into blocks.
- **How to check**: whether the marker has an authoritative reference alignment (nifH/FunGene, 16S/SILVA, COI/BOLD); whether G2 is fragmented.
- **Default**: **reference backbone first** (`hmmalign` / `papara` / `MAFFT --add --keeplength`).
- **When to deviate**: use `blocks.py` blocking only when no reference backbone is available and the diagnostics confirm the existence of several groups that cannot be aligned to one another.
- **Where to record it**: `params.md`; the Methods MUST declare that **branch lengths across blocks have no biological meaning**.

## D7 Substitution model

- **What to decide**: the model string.
- **How to check**: determined by the data type from D2b; with IQ-TREE use `-m MFP` to select automatically.
- **Default**: nucleotide `-nt -gtr -gamma` (FastTree) / `-m MFP` (IQ-TREE); protein `-lg -gamma` (FastTree) / `-m MFP` (IQ-TREE).
- **When to deviate**: codon models (`-st CODON`) for dN/dS; known taxa → you may specify HKY/GTR+Γ etc.
- **Where to record it**: `params.md`; run manifest; **G3 verifies it independently using the alignment character set**.
- **Hard constraint**: model inconsistent with the data type → HARD STOP.

## D8 Branch support

- **What to decide**: which branch support to use, and the reporting threshold.
- **How to check**: determined by the tree-inference tool (IQ-TREE → UFBoot2 + SH-aLRT; FastTree → SH-like).
- **Default**: IQ-TREE UFBoot2 1000 + SH-aLRT 1000 (strong support suggested at UFBoot2 ≥ 95 and SH-aLRT ≥ 80); FastTree SH-like 1000.
- **When to deviate** (each has a clear applicable situation, so "none" is **not** an answer for this item):
  - **N very large** (> 10,000 tips, or a single UFBoot run no longer finishes) → reduce the number of UFBoot replicates (e.g. `-B 1000` → `-B 100`) **and state it**, or switch to IQ-TREE's **`-fast`** (approximate UFBoot, slow tree search);
  - **topology only, no support** (internal exploration, figure-making) → FastTree `-nosupport`, or IQ-TREE `-B 0`; in that case you **MUST NOT** claim in the Methods that branch support was assessed;
  - **local support only** (a cheap alternative for large datasets) → FastTree `-local 1000`, and declare that this is a **local NNI test** and **not** bootstrapping;
  - **switching to raxml-ng's standard non-parametric bootstrap** → you MUST **rerun the whole tree** and declare that the scales differ (bootstrap proportion ≠ UFBoot2 ≠ SH-aLRT);
  - **Supplementary Information and main text use the same tree set** → a single rule; if two sets are genuinely needed, list them separately and declare that they are not comparable.
- **Where to record it**: `params.md`; the Methods **MUST use the declaration template in `outputs.md`**.
- **Hard constraint**: changing tool means **rerunning the whole tree**; **thresholds across tools are FORBIDDEN** ("UFBoot ≥95 or SH-like ≥0.95" side by side is the kind of misstatement reviewers catch most easily); within one paper (including Supplementary Information) the rule is unique.

## D9 Rooting

- **What to decide**: midpoint rooting or outgroup.
- **How to check**: whether a usable outgroup/reference sister lineage exists; whether the use case is sensitive to the root.
- **Default**: a **decision rule**, not a list of alternatives:

  | Condition | Requirement |
  |---|---|
  | a usable outgroup/reference sister lineage exists | **outgroup MANDATORY**, midpoint FORBIDDEN |
  | use = reference-sequence placement or community phylogenetic metrics | outgroup or reference tree; when only midpoint is possible, **rooting-sensitivity evidence MUST be provided** |
  | needs polarization/ancestral-state inference, directional evolution, or statements about the order of divergence across taxa | **outgroup MANDATORY** |
  | taxa span more than one phylum, or the alignment is > ~1,000 columns with clear rate heterogeneity | the main-tree root MUST be given by an outgroup or a reference tree |
  | single marker + homogeneous taxa + no outgroup + display only or within-group distances | midpoint allowed, and state it |
- **When to deviate**: see the table above.
- **Where to record it**: `params.md`; the `root_tree.R` report contains a **rooting sensitivity** section; the Methods state either "using outgroup X" or "no outgroup, so midpoint rooting was used, and N candidate root positions were compared".

## D10 Compute resources

- **What to decide**: local machine or server.
- **How to check**: use the estimation formulas in `scaling.md` to estimate peak memory and wall clock; **you MUST run a subsampling benchmark once to confirm**.
- **Default**: estimate first, then decide.
- **When to deviate**: estimated peak > available RAM × 0.8, or estimated wall clock > the acceptable ceiling → go to a server or downgrade.
  The first move is not to change machines, but to **go back and check the taxon rule of D3** (narrowing is usually the cheaper and more defensible choice).
- **Where to record it**: `params.md`; the benchmark commands and measured numbers go into the operation log.

## D11 Reference database and version

- **What to decide**: which reference database/HMM, which version, which threshold.
- **How to check**: see the channels in `marker-types.md`; record the version number and the file MD5.
- **Default**: `--ref-db-name` / `--ref-db-version` / `--ref-db-md5` are mandatory; annotation MUST carry `--min-id` (default 90%) and `--max-evalue`.
- **When to deviate**: low-threshold results are all labelled `unassigned` and **MUST NOT be relied on**.
- **Where to record it**: `params.md`; the Methods MUST print the database version.
- **Hard constraint**: **using reference-database annotation as evidence of phylogenetic placement is FORBIDDEN** -- it is only a supporting statement; topological conclusions can only come from the tree itself.

---

## Parameters of this study (template; copy and fill in row by row)

| ID | Value | Rationale | Reason for deviating from the default (write "using the default" if not deviating) |
|---|---|---|---|
| D1 |  |  |  |
| D2 |  |  |  |
| D2b |  |  |  |
| D3 |  |  |  |
| D4 |  |  |  |
| D5 |  |  |  |
| D6 |  |  |  |
| D7 |  |  |  |
| D8 |  |  |  |
| D9 |  |  |  |
| D10 |  |  |  |
| D11 |  |  |  |
