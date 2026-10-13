# Pitfalls table (Symptom / Root cause / Action)

Entries marked ★ are those that **lead to wrong conclusions**; prefer to look at these first.

---

## Tree building and models

| Symptom | Root cause | Action |
|---|---|---|
| ★ Likelihood values look normal, but the tree is garbage | **A protein model applied to nucleotide data** (`-wag` / `-lg` / `-jtt` / BLOSUM used on pure ACGT). The tool reports no error; it will run to completion | G3 four-layer assertions. FastTree's `WARNING! 100.0% NUCLEOTIDE CHARACTERS` and `Amino acid distances: BLOSUM45` are a veto |
| ★ Assertions misjudged on IQ-TREE | FastTree's log wording (`Nucleotide distances`, `Generalized Time-Reversible`, `Whelan-And-Goldman`) **occurs 0 times in an IQ-TREE log** (measured) | Determine the data type independently from the alignment character set (tool-agnostic); take evidence from a per-tool table; do not demand that IQ-TREE speak FastTree's words |
| ★ Using `BLOSUM` as a blacklist causes false positives | When IQ-TREE picks a protein model it **legitimately** prints candidates such as `Blosum62+I+G4` (measured: uppercase `BLOSUM` **0** times, case-insensitive `Blosum` **36** times) | Do not use `BLOSUM` as a blacklist; use the tool's native signals instead (`ModelFinder ... DNA/protein models`, `specified sequence type is different from the detected one`). **Any word-list match MUST state whether it is case-sensitive** |
| ★ COI is wrong from the very first step | The standard genetic code table was used; in mitochondria `TGA = Trp`, while the standard table treats it as a stop → **the wrong reading frame is selected** | `frame_strand.py --genetic-code 5` (invertebrates) or `2` (vertebrates); `selftest.py` already locks down the key codons of each table |
| Branch support values are not comparable across studies | SH-like / UFBoot2 / SH-aLRT mean three different things | When switching tools **rerun the entire tree**; cross-tool thresholds are FORBIDDEN; use the declaration template in `scaling.md` |
| Branch support comes from the previous tree | IQ-TREE's default output names (`*.iqtree`, `*.treefile`) are **silently overwritten** | MUST force `-pre <run-id>` at invocation; `verify_tree.py --check-same-topology` also catches rooted/unrooted inconsistency |
| The same analysis is not reproducible on another machine | FastTree version/precision builds differ, and default behaviour differs | `probe_tools.py` captures the version and archives it; the Methods cite the version number |

## Alignment

| Symptom | Root cause | Action |
|---|---|---|
| ★ Thousands of sequences, gap rate 90%+ | Non-homologous sequences mixed in / **two amplicon regions mixed into one build** / **reading frame and strand not unified**. The aligner cannot establish correspondence and can only place them side by side | G1b audit + G2 health check. Prefer a **reference backbone** (`hmmalign` / `papara` / `MAFFT --add --keeplength`); `blocks.py` is the last resort |
| ★ The alignment looks normal but is in fact misaligned | Sequences are short and locally conserved, every metric is green (normal column count, low gap rate), yet homology does not actually hold | **Detection of non-homology is the job of the G1b homology audit**; do not expect the gap rate to do it for you |
| Column count computed wrong | MAFFT output is **wrapped at 60 characters**; treating the length of each line as the column count gives absurd results | All scripts reassemble records first. If someone added `--clustalout`, `aln_qc.py --format clustal` supports that too |
| The alignment ran for two days | L-INS-i was used on 2,000 sequences (its practical ceiling is **a few hundred**) | Choose the algorithm from the two-dimensional table in `scaling.md`; for large N use FFT-NS-i or `--parttree` |
| `--auto` suddenly becomes very slow at large N | `--auto` builds an **all-pairs guide tree** (O(N²)) | Switch to `--parttree --retree 2` and declare the quality cost |
| `--thread -1` errors out or has no effect | Mistakenly assuming you must write the number of cores | `--thread -1` **is** the correct form accepted by MAFFT (help text: `if unsure, --thread -1`) |
| Strange long branches appear after block concatenation | Concatenation boundaries create discontinuities in rate and base composition; for coding genes they also break codon continuity | Prefer a reference backbone; after concatenation you MUST declare that **branch lengths across blocks carry no biological meaning** |

## Input and taxa

| Symptom | Root cause | Action |
|---|---|---|
| ★ Fewer tips than in the table; taxa silently dropped downstream | Table IDs disagree with FASTA IDs / tree tips; or someone subset-pruned an old tree and reused it | G4 two-way difference set is 0; `verify_tree.py --table` checks against the table directly; **if the taxon set changes, everything MUST be rerun** |
| ★ Reverse-complemented sequences mixed in | Strand orientation not normalised. It shows up as large stretches of **mismatches rather than gaps**, and G2 may not flag it — it exists for rRNA/ITS as well | G1b strand audit; after detection, produce a correction list and rerun the whole workflow |
| ★ Chimeras cause spurious long branches and false novel clades | Chimera removal not performed | `vsearch --uchime3_denovo` is mandatory, or the skip MUST be **explicitly recorded**; after removal, write the taxon set back |
| The header is swallowed and the row count is one short | The header was handled by "skip every line starting with `#`" | Header-detection rule: the first cell matches `#?OTU`/`ASV`/`Feature`/`id`, **or** more than half of the following cells are non-numeric |
| `reads>=10` filters out almost everything | The table holds **relative abundance**, not counts | The script detects "row sums ≈ 100 or ≈ 1" and refuses to run blindly; when the values are confirmed to be counts, pass `--assume-counts` |
| ★ Reference-database annotation treated as a topological conclusion | Annotation treated as evidence of phylogenetic affiliation | D11 forbids this explicitly; annotation is supporting text only, and MUST carry `--min-id` / `--max-evalue`; below threshold everything is labelled `unassigned` |
| Fully duplicate sequences appear after reference sequences are merged in | Multiple references for the same species, or multiple ASVs hitting the same reference | G4 checks for duplicate tip names and duplicate sequences; declare explicitly whether multiple tips are kept or collapsed |
| Not reproducible across years | Taxonomy shifts drastically between reference-database versions (SILVA/PR2/UNITE) | D11 locks the version + records the MD5, and writes it into the Methods |

## Tree files and downstream

| Symptom | Root cause | Action |
|---|---|---|
| ★ Rooting method silently swapped | "There must be a root" was implemented as "use midpoint rooting" | D9 is a **decision rule**: when an outgroup is available you MUST use it; a community-metric or reference-placement purpose needs an outgroup or sensitivity evidence; the report MUST contain a **rooting sensitivity** section |
| The tree is not binary and downstream R code behaves strangely | FastTree/IQ-TREE may output multifurcations | The community-metric route **MUST be binary**; use `multi2di` and declare the random resolution |
| Many zero-length branches | Fully duplicate or unresolved sequences (a normal phenomenon) | Report the proportion; zero-length branches affect the midpoint root position and `ape::is.binary` |
| Newick parsing fails | The file is truncated (no trailing `;`) | G4 asserts that it ends with `;` and parses |
| Negative branch lengths | Some optimisation paths produce them | G4 reports FAIL |

## Tools and environment

| Symptom | Root cause | Action |
|---|---|---|
| External tools on Windows report `Cannot read` / `exit 127` | Paths containing Chinese characters or spaces get mangled across the WSL↔Windows boundary | **Copy the binary to an ASCII path** before invoking it; or work inside WSL throughout. Use a pure-ASCII working directory |
| `sudo` asks for a password | That WSL instance has no passwordless sudo configured | `wsl -d <distro> -u root -- ...` runs directly as root |
| `apt-get install` reports 404 | The system is noble but `sources.list` points at jammy (two source sets mixed) | **Do NOT modify the user's configuration**; use a temporary sources file plus an override of apt's list paths:<br>`apt-get -o Dir::Etc::sourcelist=/tmp/only.list -o Dir::Etc::sourceparts=/dev/null -o APT::Get::List-Cleanup=0 install ...` |
| PowerShell kills the upstream process | The output is piped into `Select-Object -First N`, which terminates the upstream once it has enough | Remove `-First`, or write the output to a file and read it afterwards |
| No pip in Python | Ubuntu's python3 ships without pip, and `ensurepip` has been stripped too | Install it with the system package manager, or switch to a Python environment that bundles pip |
| `ape` has no `midpoint()` | ape 5.8-1 indeed does not have it; `phangorn` 12.1 does | `root_tree.R` uses `phangorn::midpoint`; do **NOT** use `ape::cophenetic()` — on 69,668 tips it was measured to require about 144.6 GB |
| An R package reports "built under R version 4.5.3" | It is only a warning | Ignore it; but record the R and package versions actually used |
