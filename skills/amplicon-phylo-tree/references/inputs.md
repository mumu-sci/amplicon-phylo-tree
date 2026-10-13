# Input forms and format conversion

This skill needs two things (**the second one may be absent**):

1. **Abundance table**: rows = ASV/OTU, columns = samples, the first column is the ID; values should be **counts (reads)**, and the last column is often a taxonomy string.
2. **Representative sequences FASTA**: IDs consistent with the IDs in the abundance table.

If you have only the representative sequences FASTA, use `define_taxa.py --fasta-only`; in that case you **MUST NOT claim that the result is aligned with any abundance matrix** —
the tip-set assertion in G4 loses its point of reference, and this must also be stated in the Methods.

---

## Input preconditions (confirm with the user first; do not assume)

| Precondition | What happens if it is not met | What this skill does |
|---|---|---|
| Denoising already done (DADA2 / Deblur) | The ASV table itself is wrong | Not checked; **this belongs upstream** |
| **Chimera removal** already done | Chimeras produce false long branches and false new clades; they are the number-one contamination source | Does not run it itself; `audit_inputs.py` will flag this and give the `vsearch` command, or require an **explicitly recorded skip** |
| **Orientation already normalized** | Reverse-complement sequences show up as large blocks of mismatches (not gaps), and alignment QC does not necessarily raise an alarm | **G1b checks this** (an important link this skill adds) |
| Region/primer pair **uniform** | Mixing them yields a meaningless tree | **G1b checks this** (bimodal length / primer anchoring) |
| IDs consistent between the table and the FASTA | Taxa are silently dropped | `define_taxa.py` assertion + G4 two-way set difference |

---

## Common formats → inputs of this skill

### QIIME2 (`.qza`)

```bash
qiime tools export --input-path table.qza   --output-path table_export
biom convert -i table_export/feature-table.biom -o asv_table.tsv --to-tsv
qiime tools export --input-path rep-seqs.qza --output-path seqs_export
# result: asv_table.tsv (note that it carries a '#OTU ID' header and a 'taxonomy' last column),
#         seqs_export/dna-sequences.fasta
```
> The first row of a table exported from QIIME2 starts with `#OTU ID` — `define_taxa.py` already handles this kind of header,
> but **do NOT** use a rule like "skip every line starting with `#`" (that would eat the header along with them).

### BIOM (`.biom`)

```bash
biom convert -i table.biom -o asv_table.tsv --to-tsv            # add --header-key taxonomy to keep the taxonomy column
```

### mothur

- `.shared` → use `mothur "#shared.seqs(...)"`, or convert it into a tab table;
- `.fasta` can be used directly as representative sequences; `.count_table` → tab table;
- `.taxonomy` can be merged into the table as the last column.

### phyloseq / RDS

```r
library(phyloseq)
ps <- readRDS("ps.rds")
write.table(as.data.frame(t(otu_table(ps))), "asv_table.tsv",
            sep = "\t", quote = FALSE, col.names = NA)
Biostrings::writeXStringSet(refseq(ps), "rep_set.fasta")
```
> Note the orientation of `otu_table`: only after `t()` are **the rows ASVs and the columns samples**, which is what this skill's convention requires.

### Only representative sequences FASTA / only an OTU table

- FASTA only: `--fasta-only` (see above).
- Table only, no sequences: **no tree can be built**. You must go back upstream to obtain the representative sequences.

---

## Points to watch in the table format (every one of these pitfalls has actually occurred)

| Symptom | Explanation and how to handle it |
|---|---|
| Header starts with `#` | `#OTU ID` is a legitimate header. Header-detection rule: the first cell matches `#?OTU`/`ASV`/`Feature`/`id`… **or** more than half of the cells after it are non-numeric |
| Values written as `123.0` | These are floats written out by pandas, **but they are actually integer reads**. This skill's approach is to **parse with `float()` only, never round**: `--rule 'sum>=10'` compares float values (`10.0` passes, `9.9` does not). So the rule threshold must be written in the real unit; `sum>=10` will not admit an extra record because of any "rounding" |
| `NA` / empty string / `-` | Always treated as 0 |
| Looks like relative abundance | `define_taxa.py` detects "row sum ≈ 100 or ≈ 1" and **refuses to run blind**; pass `--assume-counts` when they really are counts |
| Last column is taxonomy | Detected and excluded automatically; `--drop-last-col` can also force it |
| ASV IDs sparse/non-contiguous | Normal (e.g. `ASV1`…`ASV354855`). **NEVER** assume they are contiguous or use the numbering as an index |
| Non-sample columns mixed in | Use `--sample-cols START:END` (0-based half-open interval) to specify the sample-column range explicitly |
| Delimiter uncertain | `--format auto` compares the number of tabs versus commas in the first line |
| Sequence FASTA line-wrapped at 60 characters | All scripts reassemble the records first, then compute the lengths |

---

## Alignment conventions on the output side

The `--out-ids` written by `define_taxa.py` is the **authoritative id list for the whole pipeline**:
alignment, tree building, `verify_tree.py --ids`, and downstream metrics should all be based on it.

**Ordering discipline**: define the taxon rule → fetch the sequences → align → build the tree.
**Once the taxon rule changes, everything after it MUST be rerun**; do not take an old tree, subset it, and use it as the tree for a new taxon rule (this silently changes the branch lengths and the position of the root — a conclusion-level error).
