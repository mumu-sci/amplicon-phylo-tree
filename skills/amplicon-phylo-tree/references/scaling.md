# Scale → Algorithm → Compute (two-dimensional table)

## Why two dimensions

Cost grows approximately with **N (number of taxa) × L (alignment columns)**, not with N alone.
The burden of 1,500 columns × 20,000 sequences can exceed that of 300 columns × 100,000 sequences.
Binning by N alone gives wrong advice on "long region + moderate N" (e.g. 1,900 sequences × 4,000 bp blows memory).

**Every threshold in the table below is an empirical starting point, not a conclusion.** It is based on the scaling notes shipped with each program and on the measured anchors below.

**How to read the two-dimensional table (this file and `scripts/recommend_strategy.py` use the same set of boundaries)**:

- **N bins**: `≲200` / `200–3,000` / `3,000–100,000` / `>100,000`
- **Two L bins**: `L ≤ 2,000` ("L small") and `L > 2,000` ("L large").
  The 2,000 boundary comes from the measured lesson below: **1,900 sequences × 4,000 bp already blows memory** —
  at the same 1,900 sequences, change the order of magnitude of L and the set of usable algorithms changes too. It is not a hard criterion, it is the starting point for asking "should I change bins?".

## Alignment (MAFFT)

| N | L ≤ 2,000 | L > 2,000 |
|---|---|---|
| ≲ 200 | `mafft --auto` (it **chooses** L-INS-i itself) | `mafft --localpair --maxiterate 1000` (explicit L-INS-i; estimate machine time first, long sequences are very slow) |
| 200–3,000 | **L-INS-i**: `--localpair --maxiterate 1000` | **E-INS-i**: `--genafpair --maxiterate 1000` (long sequences usually contain unalignable segments; **you MUST estimate machine time first**) |
| 3,000–100,000 | FFT-NS-i: `--maxiterate 1000` (`--auto` also takes this path) | As at left, **and you MUST run a benchmark first**; only promote to `--parttree` if measurement proves it too slow |
| > 100,000 | `--parttree --retree 2` (**you MUST declare the alignment-quality cost in the Methods**) | `--parttree --retree 2`; and first go back and check the taxon scope in D3 |

Key points:

- **The practical ceiling for L-INS-i is a few hundred sequences**, not a few thousand. MAFFT's `--auto` selects it only up to ≲200 sequences;
  that is not conservatism, it is MAFFT's own scaling judgement. Cramming 2,000 sequences into L-INS-i will run for a very long time.
- **Do not reach for `--parttree` too early.** MAFFT's own `--auto` switches to parttree only once N reaches about 200,000;
  below that (within roughly 100,000) it uses the same progressive strategy with a memsave guide tree.
  parttree sacrifices alignment quality, and its dependence on a guide tree makes the alignment's **systematic error distribute unevenly across taxa** —
  far more dangerous than being slower. So: **choose from the table above first, run a benchmark, and upgrade only if measurement shows it is too slow**, and record the reason for the upgrade.
- On very large N, `--auto` builds an **all-pairs guide tree** (O(N²)); that is why it becomes slow in that regime.
- **`--thread -1`** is the correct form MAFFT recognises (help text verbatim: `if unsure, --thread -1`).
- Applicability of the three indel modes: **L-INS-i** (a single alignable domain, the most common case for amplicons) / **E-INS-i** (contains long unalignable regions, and also the first thing to try when the alignment fragments) / **G-INS-i** (similar lengths, globally alignable).

## Tree inference

| N | L ≤ 2,000 | L > 2,000 |
|---|---|---|
| ≲ 20,000 | **Estimate machine time first; if you can wait, use IQ-TREE** (`-m MFP -B 1000 -alrt 1000`, plus `-fast` if needed) | **Again use IQ-TREE (estimate machine time first)**; but first estimate memory generously: N × L is where it really bites (1,900 × 4,000 already blows memory) |
| 20,000–100,000 | IQ-TREE is still feasible (UFBoot2 was designed for this), but estimate first; otherwise FastTree, **and declare that it is approximate ML** | **FastTree**, unless N×L is very small (this file's criterion: N×L ≤ 3×10⁷ still allows IQ-TREE) and you have measured memory yourself |
| > 100,000 | FastTree (`-nt -gtr -gamma` / `-lg -gamma`), SH-like support values | FastTree (`-nt -gtr -gamma` / `-lg -gamma`), SH-like support values |

**Invariant**: FastTree is **approximate** maximum likelihood (GTR+CAT approximation, NNI/SPR local search); IQ-TREE is exact ML.
In the Methods you MUST NOT present the two as equivalent, and MUST NOT compare their support values side by side (see below).

---

## Measured anchors

The numbers below come from real runs during the development of this skill (the same machine and a 16-core / 64 GB cloud host).
Their purpose is to **calibrate orders of magnitude**, not to make promises.
**⚠ These are "developer-environment measurements, for order-of-magnitude reference only"** — the machine, the versions, and single/double-precision builds all differ,
so **do NOT** write the numbers here into a Methods section as fact; run a benchmark on your own machine (see below).

| N × columns | Step | Measured time | Measured peak memory |
|---|---|---|---|
| 3,874 × 455 nt | MAFFT `--auto` | 10 s | very low |
| 3,874 × 455 nt | FastTree `-nt -gtr -gamma` | **66 s** | < 0.3 GB |
| 10,687 × 4,258 nt | FastTree `-nt -gtr -gamma` | **1,543.7 s (25.7 min)** | **≈ 0.96 GB** |
| 20,000 × 455 nt | MAFFT `--auto --thread 4` | **1,149 s (19.2 min, 4 threads)** | **204.7 MB (≈ 200 MB)** |
| 30,258 × 1,560 nt | FastTree `-nt -gtr -gamma` | 3,180.9 s (53 min) | 2–4 GB |
| 69,668 × 1,560 nt | FastTree | **did not finish** on this machine | > 13.9 GB (see the warning below) |
| 190,032 × 1,560 nt | MAFFT `--parttree --retree 2` | **158 s** | **not recorded** (only the completion time was logged at the time, no peak-memory sampling — do not use this row to estimate memory) |
| 190,032 × 1,560 nt | FastTree (16 cores / 64 GB) | 2,609.8 s (43.5 min) | on the order of several GB |
| 10,687 × 3,933 nt (not normalised, contains non-homologous sequences) | MAFFT `--auto` | 147 s | low |

> ⚠ **About the 69,668 row**: that failure on this machine was **a confounded measurement** — 4 stale processes were holding about 5.8 GB at the same time,
> and the swap file was being repeatedly swapped in and out. So it only shows that "that run did not finish"; it **MUST NOT** be taken as a clean memory requirement.
> This is also why the protocol below requires "benchmark first": the true numbers on someone else's machine may differ greatly.

**Three things the anchors show**:
1. **Memory grows faster than linearly with N** (at equal column counts, going from 30,000 to 70,000 in N crossed a line this machine could not hold).
2. **Alignment column count has a very large effect on cost**: 3,874 × 455 takes only 66 seconds, while 10,687 × 4,258 takes 25.7 minutes —
   the latter has more than 9 times the columns. So **putting unalignable sequences into their own blocks** (`blocks.py`) is not only a quality matter, it is also a compute matter.
3. **MAFFT's memory is far smaller than the α constant below**: 20,000 × 455 (= 9.1M cells, N in the moderately-large range) measured a peak
   of only **204.7 MB ≈ 22 B/cell**, almost identical to FastTree's 21 B/cell on the same machine.
   So the α≈50–200 B/cell rule **has no measurement support** (see below); do not use it to conclude "not enough memory";
   its only remaining use is for **extremely long alignments** (L in the thousands of columns) — and that bin MUST be measured yourself.

---

## Resource estimation (a planning heuristic, always confirm with a benchmark)

```
Alignment peak memory ≈ α × N × L            α ≈ 50–200 B —— an unverified wide estimate: on this machine 20,000 × 455
                                              MAFFT --auto measured only ≈ 22 B/cell (204.7 MB)
ML search memory  ≈ β × N × L × number of states    plus the overhead of bootstrap/UFBoot
                                                      FastTree measured 10,687 × 4,258 at about 0.96 GB
                                                      → about 21 B/cell, but the coefficient rises when N is large

Decision: estimated peak > available RAM × 0.8  OR  estimated wall clock > acceptable ceiling
      → move to a server, or step down a bin, or (preferred) go back and check the taxon scope in D3
```

**These constants are wide, and the α one has already been contradicted by measurement** (see point 3 above): they are only for judging order of magnitude;
**do NOT** write them into a Methods section as fact, and do NOT use them in place of a benchmark.

### Benchmark protocol (mandatory)

1. Take a subsample of about 2,000 sequences (**the rule MUST be reproducible**, e.g. `seqkit sample -n 2000 -s 11`; if you simply take the first 2,000 records of the file, write down the sentence "the first 2,000 records were taken").
2. On the subsample, run the **real** alignment command and the **real** tree-inference command, and record **peak memory** (not just time):
   - Linux: `/usr/bin/time -v <command>` (look at `Maximum resident set size`)
   - macOS: `/usr/bin/time -l <command>`
   - Windows: run inside WSL so you can use `/usr/bin/time`, or sample the process RSS periodically
3. Time grows at least linearly with the number of taxa and approximately linearly with columns; **memory grows faster than linearly with the number of taxa**. When extrapolating, **round up**.
4. If the extrapolation does not fit: **the first priority is the D3 taxon scope** (a study decision), the second is a **reference skeleton / blocking**, and only the third is a bigger machine.

---

## Tool versions

`probe_tools.py` captures versions and writes them to `qc/tools.json`. They **MUST be archived with the products** and cited in the Methods —
different FastTree versions (and single-precision/double-precision builds) have different default behaviour, which makes the same analysis non-reproducible.
Calling IQ-TREE MUST include **`-pre <run-id>`**; otherwise the default output names (`*.iqtree`, `*.treefile`, `*.contree`)
**silently overwrite** the results of the previous run, so that support-value labels come from a different tree.

---

## Branch support: meaning, templates and red lines

**The three have different meanings and their values MUST NOT be compared across tools**:
- FastTree **SH-like**: a Shimodaira–Hasegawa test over local NNI (`-local 1000`).
- IQ-TREE **UFBoot2**: based on **site-likelihood resampling** (not resampling of sites).
- IQ-TREE **SH-aLRT**: a likelihood-ratio test.

### Methods template (English, ready to adapt)

> **IQ-TREE route**: The maximum-likelihood tree was inferred with IQ-TREE 2.x, with the substitution model selected by
> ModelFinder (selected model: `<model>`). Branch support was assessed with both UFBoot2 (1,000 replicates) and
> SH-aLRT (1,000 replicates), with ≥ 95 and ≥ 80 respectively taken as strong support. Note that UFBoot2 is based
> on site-likelihood resampling, so its values cannot be compared directly with standard non-parametric bootstrap;
> SH-aLRT is a likelihood-ratio test and likewise differs in meaning from a bootstrap proportion.
>
> **FastTree route**: The tree was inferred with FastTree `<version>` (GTR+CAT approximation, NNI/SPR local search)
> and is **approximate** maximum likelihood. Branch support values are SH-like local test values from 1,000 resamples.
> These values are **not equivalent to non-parametric bootstrap and MUST NOT be compared with IQ-TREE's UFBoot2 / SH-aLRT values**.
> All branch lengths and support values in this study come from a single pipeline; no results from other programs were combined.

### English counterpart (for submission)

> The maximum-likelihood tree was inferred with IQ-TREE 2.x, with the substitution
> model selected by ModelFinder (`<model>`). Branch support was assessed with both
> UFBoot2 (1,000 replicates) and SH-aLRT (1,000 replicates), requiring UFBoot2 >= 95
> and SH-aLRT >= 80 for strong support. UFBoot2 resamples site likelihoods and its
> values are **not directly comparable to standard non-parametric bootstrap**;
> SH-aLRT is a likelihood-ratio test and likewise differs in interpretation.
>
> Alternatively (FastTree): the tree was inferred with FastTree `<version>`
> (GTR+CAT approximation, NNI/SPR local search), i.e. **approximate** ML. Support
> values are SH-like local tests from 1,000 resamples; they are **not equivalent to
> non-parametric bootstrap and cannot be compared with IQ-TREE UFBoot2/SH-aLRT
> values**. All branch lengths and support values in this study come from a single
> pipeline; no results from different programs were combined.

### Three red lines

1. **Switching tools means the whole tree MUST be re-run** — putting one tool's support-value labels onto another tool's topology is FORBIDDEN.
2. **One study (including all supplementary material) uses one convention**, or the results are explicitly listed separately and declared non-comparable.
3. **Cross-tool thresholds are FORBIDDEN**: do not write "UFBoot ≥ 95 or SH-like ≥ 0.95 is taken as strong support" —
   placing two different quantities side by side is the statement error reviewers catch most easily.
