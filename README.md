# amplicon-phylo-tree

An agent skill for building **defensible, publishable phylogenetic trees from
second-generation amplicon data** (16S / 18S / ITS / nifH / amoA / COI and other
marker genes).

It is a **workflow document plus parameterised helper scripts**, not a black box:
every threshold is chosen for the study at hand and written down, every gate has an
exit code, and every override leaves a trace.

> The skill itself lives in [`skills/amplicon-phylo-tree/`](skills/amplicon-phylo-tree/).
> Point your agent at `SKILL.md` inside it.

[![verify](https://github.com/mumu-sci/amplicon-phylo-tree/actions/workflows/verify.yml/badge.svg)](https://github.com/mumu-sci/amplicon-phylo-tree/actions/workflows/verify.yml)
[![skills.sh](https://skills.sh/b/mumu-sci/amplicon-phylo-tree)](https://skills.sh/mumu-sci/amplicon-phylo-tree)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

---

## Install

**Any Agent Skills client** (the [skills CLI](https://github.com/vercel-labs/skills); this
repository is indexed by [skills.sh](https://skills.sh)):

```bash
npx skills add mumu-sci/amplicon-phylo-tree --skill amplicon-phylo-tree
```

**Claude Code plugin marketplace:**

```
/plugin marketplace add mumu-sci/amplicon-phylo-tree
```

**By hand** — copy the folder into your agent's skills directory:

```bash
cp -r skills/amplicon-phylo-tree <your-skills-dir>/amplicon-phylo-tree
```

The folder **must** stay named `amplicon-phylo-tree`: that name has to match the
`name` field in `SKILL.md`.

**Then verify the copy you installed** (three commands, no external tools needed):

```bash
cd <your-skills-dir>/amplicon-phylo-tree
python scripts/selftest.py            # coordinate, stop and genetic-code contract
python scripts/check_consistency.py   # documentation vs scripts agree
python scripts/regression_tests.py    # 36 black-box regression tests
```

### DeepSeek Harness (DSH)

DSH has no skill marketplace: it discovers skills by scanning fixed directories for
`<name>/SKILL.md`. So installation means placing the folder in one of those roots.
The most convenient one-liners:

```bash
# user-level (rank 500, the default place for personal skills)
git clone --depth 1 https://github.com/mumu-sci/amplicon-phylo-tree /tmp/aph \
  && cp -r /tmp/aph/skills/amplicon-phylo-tree ~/.agents/skills/ \
  && rm -rf /tmp/aph/skills/amplicon-phylo-tree/__pycache__
```

```powershell
# the same thing through the bundled installer, which also verifies what it copied
git clone --depth 1 https://github.com/mumu-sci/amplicon-phylo-tree $env:TEMP\aph
pwsh -File $env:TEMP\aph\install.ps1           # or:  bash install.sh
```

The installer (`install.sh` / `install.ps1`, both at the repository root) picks the
right root, refuses to overwrite a copy you have edited, strips `__pycache__` (a
`.pyc` embeds the absolute path of the machine that created it), and then runs the
three verification commands above on what it installed.

**Where DSH looks, in priority order** (the loader scans these roots and the first
match wins):

| Rank | Source | Path |
|---|---|---|
| 100 | project `.dsh` | `<project root>/.dsh/skills` |
| 200 | project `.agents` | `<project root>/.agents/skills` |
| 300 | custom | your `customSkillDirs` configuration |
| 400 | user `.dsh` | `~/.dsh/skills` |
| 500 | user `.agents` | `~/.agents/skills` (respects `DSH_AGENTS_HOME`) |
| 600 | bundled | ships with the client (`bundledSkillDir` / `DSH_BUNDLED_SKILL_DIR`) |

A project root is the nearest ancestor containing `.git`. Use the project rows to
pin a specific version to one project (`bash install.sh --project`), and the user
rows to make the skill available everywhere.

> Rank 600 is the only way a skill reaches *every* DSH user without them doing
> anything, and that list is maintained inside the client distribution — it is not
> something a third party can add to. If you think this skill belongs there, say so
> in an issue here and we can take it to the client maintainers together.

## What it is, and what it is not

It is a **workflow** with **gates**, not a one-shot script: it refuses to hand you a
tree whose model does not match the data, whose tip set drifted from the taxon list,
or whose alignment broke — and when it *cannot* check something it says so and fails
instead of passing quietly.

It is **not** a pipeline that decides science for you: it never deletes a taxon on its
own, never judges which topology is "true", and never computes your ecological
metrics. It tells you what you must state in the Methods, and gives you the artefacts
to back it up.

## Repository layout

```
amplicon-phylo-tree/
├── README.md / README.zh.md      this file, in English and Chinese
├── LICENSE                       MIT
├── CHANGELOG.md                  what changed, including the defects that were fixed
├── PRE_UPLOAD_CHECKLIST.md       the checklist this repository is released against
├── .github/workflows/verify.yml  CI: self-test, consistency, regression, smoke test
├── .claude-plugin/plugin.json    Claude Code plugin metadata
└── skills/amplicon-phylo-tree/   the skill: SKILL.md, scripts/, references/, assets/
```

## Contributing

Run the three verify commands above before opening a pull request; CI runs them again
on 3.9 and 3.12. If you change a documented default, change it in **both** the
reference document and the script — `check_consistency.py` fails the build when the
two disagree, and that property is deliberate.

---

## What it does

```
probe the environment -> fix the taxon rule -> audit the inputs -> normalise frame/strand
   -> align -> QC the alignment -> infer the tree -> root it -> verify -> deliver
```

At each step it produces the artefacts a reviewer will ask for:

| Gate | What it checks | On failure |
|---|---|---|
| **G1** | reading frame / strand resolvable proportion (coding markers) | reports only; unresolved sequences still enter the tree |
| **G1b** | input audit: orientation, homology, mixed amplicon regions, chimera step | orientation mixture -> fix and re-run; mixed regions or no homology -> **stop and report** |
| **G2** | alignment QC (fragmentation, alphabets, ambiguity policy) | structural errors **hard stop**; ratio metrics default FAIL, overridable with a recorded reason |
| **G3** | substitution model vs data type (four layers, fail-closed) | **hard stop** |
| **G4** | tree structure and tip set (bipartition identity, root degree, branch lengths) | **hard stop** |
| **G5** | traceability (run manifest, tool versions, SHA256) | reported and stop |

## What it deliberately does NOT do

- **No automatic taxon deletion.** The skill never silently drops a taxon. When an
  alignment fragments, the remedies are a reference backbone or block alignment --
  both keep every tip.
- **No truth adjudication.** It does not decide which topology is "correct"; it
  decides whether the analysis you ran is defensible as stated.
- **No downstream ecology.** PD / UniFrac / NRI-NTI / betaNTI are **not** computed.
  It produces the tree and the statements you must make about it.
- **No upstream denoising or OTU clustering.** DADA2 and friends are upstream of it.
- **No bundled binaries or databases.** Reference databases are yours to supply.
- **No fragment insertion implementation.** SEPP / EPA-ng routing and documentation
  only (an explicit design decision).

## Requirements

- **Python 3.9+** (standard library only).
- **MAFFT**, and **FastTree** or **IQ-TREE**.
- **HMMER** and/or **BLAST+** for the input audit of coding markers.
- **R 4.x** with **ape** and **phangorn** for rooting (`scripts/root_tree.R`).
- **vsearch** (recommended) for chimera screening.

`python skills/amplicon-phylo-tree/scripts/probe_tools.py --require mafft,tree`
reports what is missing and records versions.

## Known limitations

- **FastTree is approximate ML.** Its support values are SH-like local tests and are
  not comparable with UFBoot2 / SH-aLRT. The skill forces you to say so.
- **Midpoint rooting is a convention, not evidence.** With no outgroup the skill runs
  a rooting-sensitivity test by default and refuses to let a `NOT RUN` report pass
  the acceptance gate.
- **ITS is not a community-phylogenetics marker.** Cross-genus ITS trees are refused
  for metric use; within-genus trees need narrowing plus structure-aware alignment
  and belong in supplementary material.
- **The memory heuristics are planning aids.** One of them (alpha, B/cell) has been
  contradicted by measurement; the documents say so rather than hiding it.
- **Measured anchors are from a developer machine** and are labelled as such: they
  calibrate an order of magnitude, nothing more. Run the benchmark.

## Licence

MIT. See [LICENSE](LICENSE).

## Version

1.0.1 — see [CHANGELOG.md](CHANGELOG.md).
