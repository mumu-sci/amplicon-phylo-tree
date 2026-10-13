# Pre-upload checklist

Run this before publishing the repository or handing the skill to someone else.
Every box is something that has actually gone wrong at least once in this project's
history; none is hypothetical.

Legend: **[auto]** = a command in `skills/amplicon-phylo-tree` verifies it;
**[manual]** = you must read and decide.

---

## A. Correctness (the code does what it says)

- [ ] **[auto]** `python scripts/regression_tests.py` -> `RESULT: PASS`, and the
      count matches what `CHANGELOG.md` claims.
- [ ] **[auto]** `python scripts/selftest.py` -> `RESULT: PASS`.
- [ ] **[auto]** `python scripts/check_consistency.py` -> `RESULT: PASS`
      (five check groups: D1–D11 fields, scaling grid, documented `--flag`s,
      script paths, routing rules).
- [ ] **[auto]** The end-to-end smoke run reproduces the documented baseline
      (20 -> 15 taxa kept, 556 columns, 46.04% gaps). Any change to those numbers
      is a regression unless it is explained here.
- [ ] **[auto]** Every "will be checked" promise in `SKILL.md` names a check that a
      script really performs. If `check_consistency.py` cannot see it, add a
      regression test that can.
- [ ] **[manual]** For each gate in the summary table: can it FAIL? A gate that can
      only pass is worse than no gate, because it reads as assurance. Grep the gate
      for its exit path and confirm it is reachable.

## B. Content compliance (it says only what it can support)

- [ ] **[manual]** No claim of having done something the skill does not do
      (downstream metrics, insertion, denoising). Search for the words
      `UniFrac`, `betaNTI`, `PD`, `SEPP`, `EPA-ng` and check the sentence around
      each: routing and documentation yes, performed no.
- [ ] **[manual]** "No default" statements are honest. D1/D2/D3 must not have an
      implicit default hiding in the scripts; D8 must not claim there are no
      deviation cases.
- [ ] **[manual]** Thresholds that came from convention are labelled as convention,
      not as citable findings (85/90% in D2b, 95/80 support cut-offs in D8).
- [ ] **[manual]** Every measured number carries its provenance and its limits: the
      machine, the version, and whether the measurement was clean. The confounded
      69,668-tip measurement must stay labelled as confounded.
- [ ] **[manual]** Limitations that weaken the skill are still present (ITS cannot
      be used for community metrics; midpoint rooting is a convention; FastTree is
      approximate ML). This is the failure mode reviewers actually catch, and the
      temptation is always to delete these lines.

## C. Release form (the repository is usable)

- [ ] **[manual]** `README.md` (EN) and `README.zh.md` describe the same behaviour;
      neither promises a feature the other denies.
- [ ] **[manual]** The README's dependencies, install line and verify commands run
      as written on a clean checkout.
- [ ] **[manual]** `LICENSE` (MIT) is present, and `SKILL.md` frontmatter carries
      `license: MIT`.
- [ ] **[manual]** `CHANGELOG.md` lists the version in `SKILL.md` metadata and is
      written honestly (including the false passes).
- [ ] **[manual]** `.github/ISSUE_TEMPLATE/` has a bug template and a "wrong tree"
      template, both requiring the version and the tool log.
- [ ] **[auto]** The skill folder contains no `README.md` (the skill validator
      forbids it; human docs live at the repository root) and its folder name equals
      the `name` in `SKILL.md`.

## D. Declarations and boundaries

- [ ] **[auto]** Identity scan: no developer path, user name, machine name, project
      name, data file or personal e-mail inside `skills/amplicon-phylo-tree`.
      Run it **after** deleting `__pycache__`: a `.pyc` embeds the absolute path of
      the checkout, so bytecode left behind by a local run is a real leak even
      though no source file mentions it. Purge every `__pycache__` in **both** the
      working and the installed copy.
- [ ] **[manual]** No third-party binary or reference database is bundled.
- [ ] **[manual]** The development warning at the top of `SKILL.md` has been removed
      (it must not ship), and nothing that replaced it weakens the invariants below.
- [ ] **[manual]** The hard invariants are still stated, and still have teeth:
      **no automatic taxon deletion**, **no silent "unchecked = passed"**, **no
      brittle equality against tool wording**, **model must match the data type**.

---

## The four hard invariants

These have each been broken once by a "reasonable-looking" change. Check them by
hand, not by memory:

1. **A taxon is never dropped silently.** Unresolved sequences are written out and
   marked; exclusion requires an explicit flag and prints a warning.
2. **An unchecked gate never passes.** `NOT CHECKED` affects the exit code;
   waivers require a written reason and land in the machine-readable record.
3. **Tool evidence is matched per tool, case-sensitively, and never as a blacklist
   alone.** "The word did not appear" is not positive evidence.
4. **The declared data type, the substitution model and the alignment characters
   agree**, decided by the characters themselves, and a mismatch stops the run.
