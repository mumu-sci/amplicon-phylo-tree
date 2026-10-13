# Marker gene classes (axis B)

The difference between these classes is not "a different database" -- it is **the data type, whether a reading frame is required, whether alignment can hold at all, and how deep the conclusions can be pushed**.

---

## 1. `rrna` — 16S / 18S / 25–28S

| Item | Notes |
|---|---|
| Data type | Nucleic acid (**no translation**) |
| Reading frame | Not required |
| Model | `-nt -gtr -gamma` / IQ-TREE `-m MFP` |
| Alignment | Sequence lengths within an amplicon region are relatively consistent; for large N see `scaling.md` |
| Reference database | SILVA / Greengenes (16S); **PR2** (18S/eukaryotes) |
| Interpretive ceiling | The amplicon region is limited, **resolution within genus/species is low**; conclusions should mainly sit at higher taxonomic ranks. 16S is multi-copy with minor differences between copies, which introduces noise into the genealogy |

**Two real obstacles (they MUST be handled; it is not enough to write "use database X")**:
1. **Amplicon region and reference database MUST be paired**. 18S V4 and V9 differ greatly in length; **PR2 and SILVA taxonomies are incompatible** -- pick one and use it all the way through; do not mix.
2. **Deep-branching taxa are mixed in**. Marine 18S commonly contains parasitic/symbiotic groups (Syndiniales, MAST groups) mixed with host sequences, and they are extremely deeply divergent on 18S.
   Once they are mixed in, the default of "retain the full taxon set" will produce a tree with **fragmented alignment**, and G2 will not necessarily flag it.
3. **Taxonomy changes drastically between reference database versions**, so the version MUST be locked (D11); otherwise the analysis cannot be reproduced across years.

**Arbitration rule ("stratified tree building" vs "retain all taxa" -- choose one, no vagueness allowed)**:
These two statements contradict each other in the older documentation -- stratification = dropping taxa = changing the domain of the metric. The rule is as follows:

| Purpose of the tree (axis A) | Rule |
|---|---|
| **Distance skeleton / community phylogenetic metric** (PD / UniFrac / NRI-NTI / βNTI) | It is **FORBIDDEN to stratify and FORBIDDEN to drop any taxon**. Solve the deep-branching taxon problem with **reference backbone alignment** (`papara` / `MAFFT --add --keeplength` / `hmmalign`) or with **block concatenation**; both retain all tips. If it is still fragmented, **stop and report** -- do not silently drop taxa |
| **The conclusion itself / display** (placement of a novel clade, within-taxon relationships, figure production) | Stratifying **by reference taxonomic depth or building a separate tree is allowed**, but you MUST: ① write the stratification rule and its basis into `params.md`; ② report the retained **reads% and richness%** for each stratum; ③ declare explicitly in the methods that "this tree does not contain taxon group X, so it cannot be used for community metrics" |
| Unsure | Take the strictest option (retain all taxa) and state in the operation log why |

> In one sentence: **under the metrics framing "stratification" is the wrong thing to do, while under the conclusion framing it is legitimate narrowing** --
> the distinction is not in the method, it is in what you use this tree to compute. See D4 and the "ask these three things first" section of SKILL.md.

---

## 2. `its` — fungal ITS1 / ITS2 / full-length ITS, plant ITS

> ### ⚠ ITS **should not be the primary marker for a community phylogenetic tree**. Not building a tree is the default; building one is the exception.

**Three reasons**:
1. Length heterogeneity is extreme -- **cross-genus and above cannot be reliably aligned**; a single ITS1 alignment cannot simultaneously accommodate phylum A and phylum B.
   (Note the qualification: **within genus**, with the taxon set already narrowed, ITS trees are common in the literature; what this skill opposes is treating ITS
   as the primary marker for **cross-taxon / community metrics**, not opposing all ITS trees.)
2. **Intragenomic variation + incomplete concerted evolution**: the ITS of a single ASV may not represent one genomic locus, and the "species" on the tree is not a biological species.
3. The field consensus is to name taxa using **UNITE species hypotheses (SH) + match percentage**, and **not** using ITS topology.

**Priority (by degree of recommendation)**:

| Order | Approach | When it applies |
|---|---|---|
| ① **Do not build a tree** | UNITE dynamic clustering + `--min-id 97/99` match-based naming; route diversity/community metrics through **18S or LSU D1/D2** trees instead | Default; especially when the downstream use is community metrics |
| ② **Fragment insertion** | SEPP / EPA-ng insertion into a UNITE/SILVA reference tree | When "one comparable tree" is needed; preserves comparability and avoids having to solve cross-phylum alignment |
| ③ **Narrow the taxon set + secondary-structure-constrained alignment** | The taxon set has already been narrowed to a single genus / closely related species complex; use ITS2 secondary-structure-constrained alignment (or at least per-family local alignment) | **Supplementary material only** |

**Supplementary markers beyond ITS** (reviewers will ask "then why didn't you use another locus"):
**RPB2** (RNA polymerase II second largest subunit, one of the best resolvers between fungal genera/species), **TEF1α** (translation elongation factor 1α),
**β-tubulin** (`benA`). All three are protein-coding with a readable frame and **MUST be handled under D2b** -- they are better suited than ITS for between-genus phylogenetics.

**Tools for secondary-structure-constrained alignment (used by ③)**: **ITS2 database** (ready-made structure models + alignments), **4SALE**,
**RNAsalsa** (structure alignment), **ITSx** (cut ITS1/ITS2 apart first; do not align them mixed together),
**Gblocks + rRNA structure model** (remove unreliable regions, provided there is structural evidence).

**Hard constraint (qualified version)**: when "community phylogenetic metric (PD/UniFrac/NRI-NTI/βNTI)" is selected in axis A,
**for cross-genus and above, building an ITS1/ITS2 tree alone is FORBIDDEN** -- that combination cannot be defended methodologically. Route it through 18S/LSU instead, or use fragment insertion.
**Within genus**, if an ITS tree is genuinely needed: the taxon set MUST first be narrowed, secondary-structure/insertion constraints MUST be added, and it may be used **only as supplementary material**; the main text MUST NOT use it to draw cross-taxon conclusions.

③ If it is to appear in the main text/supplementary material, the methods and the figure caption MUST both state:
"ITS length heterogeneity and intragenomic variation limit the resolution and topological confidence of this tree; it is used only to display within-group relationships and
is not used for cross-taxon naming or community phylogenetic metrics."

---

## 3. `functional-nuclear` — nifH / amoA / nirK / nirS / pmoA / nosZ …

| Item | Notes |
|---|---|
| Data type | **Determined by the D2b criterion** (nucleic acid or amino acid; see `parameters.md`) |
| Reading frame | **Required**: unify the reading frame and strand orientation first (`frame_strand.py`), otherwise the alignment will necessarily fragment |
| Model | Nucleic acid `-nt -gtr -gamma`; protein `-lg -gamma`; IQ-TREE `-m MFP` |
| Reference database / HMM | Pfam/TIGRFAM HMM (e.g. **Fer4_NifH = PF00142**), **FunGene**, ARB nifH 2017 |
| Genetic code | Usually the standard table (1); **confirm per taxon** |
| Interpretive ceiling | **Deep nodes across functional clusters MUST NOT be used as evolutionary conclusions** (nifH Cluster I–IV, bacterial-type vs archaeal-type amoA, nirK vs nirS are all ancient paralogous lineages); a ~300 bp amplicon at the nucleotide level cannot resolve deep cross-cluster relationships |

**Why the reading frame MUST be normalized**: with mixed reading frames/strand orientations, the aligner cannot establish correspondences,
so it can only place them side by side → matrix explosion (measured: 10,687 un-normalized sequences gave 3,933 columns / 92.06% gaps;
after normalization, the homologous subset of the same batch gave 455 columns / 28.73% gaps). This is a **technical necessity**, not a judgment call.

**Why cross-cluster comparison MUST go through amino acids**: the third codon position is saturated between distant taxa, and the nucleotide tree will return long branches caused by long-branch attraction.
Building both trees and reporting topological concordance (RF / conflicting branches) is MANDATORY in the grey zone.

---

## 4. `organelle` — COI / rbcL / matK (protein-coding); 12S belongs to rRNA

| Item | Notes |
|---|---|
| Data type | Determined by D2b; for COI the **direction is empirical, not conventional**: in the Leray 313 bp **within-genus/between-species** alignment the third position happens to be the most informative layer, so the usual route is **nucleic acid + MACSE codon alignment**; only cross-order/cross-class comparisons favor amino acids. `-st CODON` is the choice for dN/dS or saturation diagnostics |
| Reading frame | **Required** (except 12S: it is a mitochondrial rRNA gene, **no translation**, handle it as `rrna`) |
| **Genetic code table** | **COI MUST be specified**: invertebrates **5**, vertebrates **2** (see `frame_strand.py --genetic-code`) |
| Reference database | **BOLD** / GenBank |
| Reference backbone | COI (Leray 313 bp) **has no HMM**: `--mode hmm` runs idle, so use `--mode orf` or a **reference backbone** |
| **Interpretive ceiling** | 313 bp ≈ **104 codons**. At the species level/nearby genera the **third position is the most informative**, which is why short fragments work well at species level; but a **short-fragment de novo tree and a reference tree cannot be compared by branch length** (the branch lengths of the insertion route come from the insertion), and **deep nodes are unreliable**. `frame_strand.py pick` already computes internal stops per sequence (`stop_codons`) -- this is the **first signal of numts / frameshift**; when it appears, read the WARN list for that column and do not treat it directly as a "novel clade" |

> ⚠ **The default genetic code table gets COI wrong from the very first step**: the standard table treats `TGA` as a stop, whereas the mitochondrial code reads `TGA = Trp`.
> The result is that a wrong reading frame is selected (`TGA` positions are recorded as `X`/truncated), and G1 only reports this without blocking.
> `scripts/selftest.py` already locks the key codons of each genetic code table (including the "standard table vs table 5" comparison); it MUST be run after any code change.

**Short-fragment warning**: for 313 bp COI / 250 bp 16S V4, de novo alignment loses comparability with the reference tree --
**fragment insertion** (SEPP / EPA-ng) is usually the better-accepted choice.

---

## 5. `other` — any custom amplicon

- Use the **character-set determination** of `verify_tree.py` to first confirm whether it is nucleic acid or protein.
- Reference resources are supplied by the user; when there is no reference, prefer `--mode orf` to normalize the reading frame (if it is a coding gene).
- If diagnostics show multiple groups that cannot be aligned against each other, split into blocks with `blocks.py` (last resort).

---

## Reference resource channels (this skill does not download)

| Marker gene | HMM | Reference database |
|---|---|---|
| nifH / vnfH / anfH | Pfam **PF00142** (Fer4_NifH) | FunGene, ARB nifH 2017 |
| amoA (bacterial/archaeal) | Pfam PF02461 (bacterial) / PF12998 approximate | FunGene |
| nirK / nirS | PF02239 / PF13423 (approximate) | FunGene |
| pmoA / nosZ | PF02461 family / PF03232 | FunGene |
| 16S / 18S | — (no domain search) | SILVA, Greengenes, PR2 |
| ITS | — | UNITE (+ ITSx splitting) |
| COI / 12S | — (no HMM route for COI) | BOLD, GenBank |

> The purpose of a domain search is to **determine the reading frame and strand orientation** and to **report homology**, not to supply topological conclusions.
> **It is FORBIDDEN to treat reference database annotation as phylogenetic affiliation evidence** (see `parameters.md` D11).
