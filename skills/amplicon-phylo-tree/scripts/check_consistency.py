#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Check that the DOCUMENTATION and the SCRIPTS still agree.

Why this exists: the most expensive defects found by adversarial review were not
crashes -- they were promises the workflow made that the scripts did not keep:

  * SKILL.md documented a ``--domtblout-dir`` flag that ``pick`` did not accept,
    so following the documented command exited 2 (R2's B4).
  * ``parameters.md`` D5 documented ``--ambiguous-policy``, which existed nowhere
    in the code (R2's M3).
  * ``scaling.md`` advertised a two-dimensional table while ``alignment_plan()``
    ignored L entirely (R2's 7.1).

Each of those is mechanically findable.  A reviewer reading 60 KB by hand is not.
So this script asserts four things and exits non-zero on any failure:

  1. D1-D11 in ``references/parameters.md`` each carry all five protocol fields
     (要决定什么 / 怎么查 / 默认 / 何时偏离 / 记到哪里) with non-trivial content.
  2. The N x L grid in ``references/scaling.md`` matches what
     ``recommend_strategy.py`` actually returns, cell by cell.
  3. Every ``--flag`` written in SKILL.md / references/ is a flag some script
     really accepts (the B4 / M3 class of drift).
  4. ``scripts/*.py`` and ``scripts/*.R`` files named in the docs exist, and the
     skill's own smoke-test path is coherent.

Run:  python check_consistency.py            (from the skill folder or anywhere)
Exit: 0 = consistent, 1 = drift found, 2 = usage/IO error.
"""
import argparse
import ast
import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
REFS = os.path.join(SKILL, "references")
DOCS = [os.path.join(SKILL, "SKILL.md")]
if os.path.isdir(REFS):
    DOCS += [os.path.join(REFS, f) for f in sorted(os.listdir(REFS))
             if f.endswith(".md")]

# Scripts whose argument parser is part of the public contract.
SCRIPT_NAMES = ["aln_qc.py", "audit_inputs.py", "blocks.py", "define_taxa.py",
                "frame_strand.py", "probe_tools.py", "recommend_strategy.py",
                "root_tree.R", "selftest.py", "verify_tree.py"]

# ---------------------------------------------------------------- helpers ----

def read(path):
    with io.open(path, encoding="utf-8", errors="replace") as fh:
        return fh.read()


def script_options(path):
    """Every long option a script's parser accepts.

    Python: parse the AST and collect add_argument("--x") literals, plus options
    declared through an add_argument_group or a choices list.  R: collect the
    `--name` string literals in root_tree.R's argument loop.
    """
    text = read(path)
    opts = set()
    if path.endswith(".R"):
        opts |= set(re.findall(r'"(--[A-Za-z0-9-]+)"', text))
        opts |= set(re.findall(r"'(--[A-Za-z0-9-]+)'", text))
        return opts
    try:
        tree = ast.parse(text)
    except SyntaxError as exc:
        raise ValueError("%s does not parse: %s" % (path, exc))
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        fn = node.func
        name = getattr(fn, "attr", None) or getattr(fn, "id", None)
        if name != "add_argument":
            continue
        for arg in node.args:
            if isinstance(arg, ast.Constant) and isinstance(arg.value, str) \
                    and arg.value.startswith("-"):
                opts.add(arg.value)
        for kw in node.keywords:
            if kw.arg == "choices" and isinstance(kw.value, (ast.List, ast.Tuple)):
                for el in kw.value.elts:
                    if isinstance(el, ast.Constant) and isinstance(el.value, str) \
                            and el.value.startswith("-"):
                        opts.add(el.value)
    return opts


def all_options():
    table = {}
    for name in SCRIPT_NAMES:
        path = os.path.join(HERE, name)
        if not os.path.exists(path):
            table[name] = None            # missing script: reported elsewhere
            continue
        table[name] = script_options(path)
    return table


# The five-field protocol.  Each field is matched in either language, so the check
# works on the Chinese original and on the English translation alike -- and a
# translation that drops a field still fails.
FIELDS = [
    ("要决定什么", "What to decide"),
    ("怎么查", "How to check"),
    ("默认", "Default"),
    ("何时偏离", "When to deviate"),
    ("记到哪里", "Where to record it"),
]
# A field may legitimately be a bullet OR a table; what is not legitimate is a
# placeholder ("无", "n/a", "-") where a decision rule is required.
TRIVIAL = {"无", "n/a", "na", "none", "-", "待定", "tbd", "—"}


def check_parameters(text):
    """Return (ok, [messages]) for the D1-D11 five-field protocol."""
    msgs = []
    ids = re.findall(r"^##\s+(D\d+b?)\s", text, re.M)
    if not ids:
        return False, ["parameters.md declares no D-numbered sections"]
    expected = ["D1", "D2", "D2b", "D3", "D4", "D5", "D6", "D7", "D8", "D9",
                "D10", "D11"]
    missing = [d for d in expected if d not in ids]
    if missing:
        msgs.append("parameters.md is missing sections: %s" % ", ".join(missing))

    # split the document at each section heading
    parts = re.split(r"^##\s+(D\d+b?)\s.*$", text, flags=re.M)
    bodies = {}
    for i in range(1, len(parts) - 1, 2):
        bodies[parts[i]] = parts[i + 1]

    for did in expected:
        body = bodies.get(did)
        if body is None:
            continue
        for zh, en in FIELDS:
            m = re.search(r"\*\*(?:%s|%s)\*\*\s*[:：]?(.*)"
                          % (re.escape(zh), re.escape(en)), body)
            if not m:
                msgs.append("%s has no **%s / %s** field" % (did, zh, en))
                continue
            rest = m.group(1).strip().strip("。.")
            # a table right under the field is content too
            tail = body[m.end():m.end() + 400]
            has_table = re.search(r"^\s*\|", tail, re.M) is not None
            if (rest and rest not in TRIVIAL) or has_table:
                continue
            if rest.startswith("（") and "见" in rest and len(rest) > 4:
                continue                      # e.g. "（见上表）" with a table
            msgs.append("%s **%s** is empty or a placeholder (%r)"
                        % (did, field, rest[:40]))
    return (not msgs), msgs


# ------------------------------------------------------- 2. scaling grid -----

# The 比对/建树 tables.  The headings are matched in either language so the check
# survives the Chinese->English translation without being weakened.
ALIGN_HEADINGS = ("比对（MAFFT）", "Alignment (MAFFT)")
TREE_HEADINGS = ("建树", "Tree inference")

# The N x L grid the checker probes.  The expected values are READ OUT OF
# scaling.md's two tables, so the assertion compares DOCUMENT against CODE -- never
# a hand-copied constant against itself.
GRID = [(50, 500), (50, 4000), (1500, 500), (1500, 4000), (20000, 500),
        (20000, 4000), (190000, 500), (190000, 4000)]


def scaling_expected(text):
    """Parse the two tables of scaling.md into {(n, L): (alignment, tree)}.

    Cells repeat one value across an N row; the parser takes the row's two cells
    as the L<=2000 and L>2000 answers.  A cell written "同左 / 同上 / 同左；..." means
    "same algorithm as the row above, with extra caveats", which is resolved before
    comparison so the documentation can stay readable without weakening the check.
    """
    def rows_after(headings):
        for heading in headings:
            m = re.search(r"^##\s+%s\s*$" % re.escape(heading), text, re.M)
            if not m:
                continue
            chunk = text[m.end():]
            out, seen_table = [], False
            # The table is one contiguous block: stop at the first non-table line
            # after it, otherwise the NEXT section's table is read into this one.
            for line in chunk.splitlines():
                s = line.strip()
                if s.startswith("|"):
                    seen_table = True
                    cells = [c.strip() for c in s.strip("|").split("|")]
                    if len(cells) < 3:
                        continue
                    if cells[0] in ("N", "") or set(cells[0]) <= set("-: "):
                        continue
                    out.append(cells)
                elif seen_table:
                    break
            # Resolve an inheritance marker ("同左/同上" or "As at left / Same as
            # above") against the row above, appending whatever extra caveat
            # follows so nothing is lost.
            for i in range(1, len(out)):
                for col in (1, 2):
                    cell = out[i][col]
                    if not (re.match(r"^(同[左上右]|见上)", cell)
                            or re.match(r"^(as at left|same as above|as above|"
                                        r"as in the row above)", cell, re.I)):
                        continue
                    extra = re.sub(r"^(同[左上右]|见上)[；;，,、]?\s*", "", cell)
                    extra = re.sub(r"^(as at left|same as above|as above|as in the "
                                   r"row above)[,:;]?\s*", "", extra, flags=re.I)
                    inherited = out[i - 1][col]
                    if inherited.startswith("同"):
                        continue        # never inherit a marker; later pass
                    out[i][col] = inherited + ("; " + extra if extra else "")
            return out
        return []

    align_rows = rows_after(ALIGN_HEADINGS)
    tree_rows = rows_after(TREE_HEADINGS)

    def pick_answer(rows, n, L):
        """Rows say '≲ 200' / '200–3,000' / '3,000–100,000' / '> 100,000'."""
        for cells in rows:
            label = cells[0]
            nums = [int(x.replace(",", "")) for x in re.findall(r"[\d,]+", label)]
            small = label.startswith("≲") or label.startswith("≤") or label.startswith("<")
            greater = label.startswith(">")
            matched = False
            if small and nums and n <= nums[0]:
                matched = True
            elif greater and nums and n > nums[0]:
                matched = True
            elif not small and not greater and len(nums) == 2 and nums[0] < n <= nums[1]:
                matched = True
            if matched:
                col = 1 if L <= 2000 else 2
                if col >= len(cells):
                    return None
                return cells[col]
        return None

    expected = {}
    for n, L in GRID:
        a = pick_answer(align_rows, n, L)
        t = pick_answer(tree_rows, n, L)
        if a is not None and t is not None:
            expected[(n, L)] = (a, t)
    return expected, align_rows, tree_rows


# Order matters: a cell often mentions more than one algorithm, sometimes to say
# "escalate to this ONLY if slow" ("FFT-NS-i ...; 只有过慢才升 --parttree").  The
# recommended action is therefore matched LAST and only against the part of the
# cell before the first conditional clause marker.
ALIGN_KEYS = ["FFT-NS-i", "E-INS-i", "L-INS-i", "parttree", "--auto"]
CONDITIONAL = re.compile(r"只有|才升|才升到|再升|除非|if it is slow|only if")


def classify_alignment(cells):
    """Name the algorithm a doc cell or a code return value RECOMMENDS.

    The escalation caveat ("parttree only if the benchmark proves it necessary") is
    not the recommendation, so the conditional tail is cut before matching.
    """
    head = CONDITIONAL.split(cells)[0] if CONDITIONAL.search(cells) else cells
    for k in ALIGN_KEYS:
        if k.lower() in head.lower():
            return "auto" if k == "--auto" else k
    # the recommendation may be stated only in the conditional tail
    for k in ALIGN_KEYS:
        if k.lower() in cells.lower():
            return "auto" if k == "--auto" else k
    return None


def classify_tree(cells):
    low = cells.lower()
    if "iq-tree" in low:
        return "IQ-TREE"
    if "fasttree" in low:
        return "FastTree"
    return None


def check_scaling(doc_text):
    msgs = []
    expected, arows, trows = scaling_expected(doc_text)
    if not arows or not trows:
        return False, ["scaling.md: could not parse the 比对/建树 tables"]

    # import the implementation's own grid function
    sys.path.insert(0, HERE)
    import recommend_strategy as rs                        # noqa: E402

    mism = 0
    for n, L in GRID:
        doc_cells = expected.get((n, L))
        if doc_cells is None:
            msgs.append("scaling.md covers no cell for N=%d L=%d" % (n, L))
            continue
        doc_align, doc_tree = doc_cells
        aid = classify_alignment(doc_align)
        tid = classify_tree(doc_tree)
        m, _cmd, _why, _cost = rs.alignment_plan(n, L)
        code_align = classify_alignment(m)
        tool, _t, _s, _note = rs.tree_plan(n, L)
        code_tree = classify_tree(tool)
        if aid != code_align:
            mism += 1
            msgs.append("N=%d L=%d: scaling.md says %r, code returns %r"
                        % (n, L, aid or doc_align[:30], code_align or m[:30]))
        if tid != code_tree:
            mism += 1
            msgs.append("N=%d L=%d: scaling.md says %r, code returns %r"
                        % (n, L, tid or doc_tree[:30], code_tree or tool[:30]))
    return (mism == 0), msgs


# --------------------------------------------------------- 3. flags in docs --

FLAG_RE = re.compile(r"(?<![\w-])(--[A-Za-z][A-Za-z0-9-]*)")

# Flags that legitimately appear in this skill's PROSE but belong to tools it
# drives, or are placeholders it deliberately leaves open.  Kept explicit so that
# adding a new external tool is a conscious edit rather than a silent widening.
EXTERNAL_FLAGS = {
    "--clustalout", "--add", "--keeplength", "--preserve", "--uchime3_denovo",
    "--uchime_ref", "--nonchimeras", "--chimeras", "--usearch_global", "--uc",
    "--strand", "--maxaccepts", "--maxrejects", "--db", "--query", "--evalue",
    "--max_target_seqs", "--num_threads", "--outfmt", "--out", "--domtblout",
    "--cut_tc", "--noali", "--cpu", "--localpair", "--maxiterate", "--auto",
    "--parttree", "--retree", "--genafpair", "--thread", "--quiet", "--reorder",
    "--local", "--nosupport", "--redo", "--prefix", "--seed", "--id",
    "--model", "--outdir",
}

# Flag-looking tokens in this skill's own prose that are NOT command-line options
# but names of a parameter-table field or of another tool's option.
VIEWER_FLAGS = {
    "--min-id", "--max-evalue", "--ref-db-name", "--ref-db-version",
    "--ref-db-md5", "--cut", "--uchime3", "--input-path", "--output-path",
    "--to-tsv", "--header-key",
}


def doc_prose(path):
    """The document with fenced code-block BODIES removed (fences kept).

    Options written inside a fenced command block belong to whatever tool that
    block invokes (MAFFT, FastTree, IQ-TREE, HMMER, vsearch, ITSx, QIIME2, biom,
    seqkit ...).  Only prose and inline code are this skill's own contract, so only
    they are checked.
    """
    out, in_fence = [], False
    for line in read(path).splitlines():
        if line.strip().startswith("```"):
            in_fence = not in_fence
            out.append(line)
            continue
        if in_fence and line.startswith((" ", "\t")) and line.strip():
            out.append("")                    # body line: belongs to another tool
        elif in_fence:
            out.append(line)                  # e.g. continuation arrows, `mafft --auto`
        else:
            out.append(line)
    return "\n".join(out)


def check_flags(option_table):
    msgs = []
    known = set()
    for opts in option_table.values():
        if opts:
            known |= opts
    placeholder = re.compile(r"--[xX]+$")
    for path in DOCS:
        if not os.path.exists(path):
            msgs.append("documented file is missing: %s" % path)
            continue
        text = doc_prose(path)
        for m in FLAG_RE.finditer(text):
            flag = m.group(1)
            if flag in known or flag in EXTERNAL_FLAGS or flag in VIEWER_FLAGS:
                continue
            if placeholder.match(flag) or flag.endswith("-path"):
                continue
            line = text[:m.start()].count("\n") + 1
            msgs.append("%s:%d: %s is not a flag any script accepts"
                        % (os.path.relpath(path, SKILL), line, flag))
    return (not msgs), msgs


# -------------------------------------------------- 4. script references -----

SCRIPT_RE = re.compile(r"(?:scripts/)([A-Za-z0-9_.-]+\.(?:py|R))")


def check_scripts_exist():
    msgs = []
    for path in DOCS:
        if not os.path.exists(path):
            continue
        text = read(path)
        for m in SCRIPT_RE.finditer(text):
            target = os.path.join(HERE, m.group(1))
            if not os.path.exists(target):
                line = text[:m.start()].count("\n") + 1
                msgs.append("%s:%d: scripts/%s does not exist"
                            % (os.path.relpath(path, SKILL), line, m.group(1)))
    return (not msgs), msgs


# --------------------------------------------------- 5. route distinctness ---

ROUTE_CASES = [
    ("rrna", "skeleton"), ("rrna", "conclusion"), ("its", "skeleton"),
    ("its", "placement"), ("its", "conclusion"), ("functional-nuclear", "skeleton"),
    ("functional-nuclear", "conclusion"), ("organelle", "skeleton"),
    ("organelle", "conclusion"), ("other", "figure"),
]


def check_routes():
    """The routing invariants that the documents promise must hold.

    What is asserted (each one is written down somewhere, so it is checkable):

      * every marker-class x tree-role combination produces a non-empty route;
      * ITS in a community/placement role prints NO de novo alignment or tree
        command (B7) -- an agent executes what it is shown;
      * every other combination DOES print both;
      * the recommendation is genuinely two-dimensional: N changes it and L
        changes it (axes A and B are documented, by design, as being about what you
        KEEP and what you must PROVE, not about silently swapping the algorithm, so
        they are not asserted to change the tool).
    """
    msgs = []
    sys.path.insert(0, HERE)
    import recommend_strategy as rs                        # noqa: E402
    import subprocess                                      # noqa: E402
    src = read(os.path.join(HERE, "recommend_strategy.py"))
    if not re.search(r'marker_class == "its" and args\.tree_role', src):
        msgs.append("recommend_strategy.py no longer special-cases ITS in a "
                    "community/placement role (B7 regression)")
    if "alignment_plan(n, L)" not in src:
        msgs.append("recommend_strategy.py no longer passes L to alignment_plan: "
                    "the two-dimensional table is decorative again")

    # N changes the answer
    labels = {rs.alignment_plan(n, 500)[0] for n in (50, 1500, 20000, 190000)}
    if len(labels) < 4:
        msgs.append("alignment_plan ignores N: %d distinct answers for 4 sizes"
                    % len(labels))
    # L changes the answer (this is the 7.1 regression)
    if rs.alignment_plan(1500, 500)[0] == rs.alignment_plan(1500, 4000)[0]:
        msgs.append("alignment_plan ignores L at N=1,500 -- the two-dimensional "
                    "table is decorative again")
    if rs.tree_plan(30000, 100)[0] == rs.tree_plan(30000, 4000)[0]:
        msgs.append("tree_plan ignores L at N=30,000 (the N x L residue clause is "
                    "not being applied)")

    # every marker x role combination produces a route, and ITS in a community or
    # placement role produces NO de novo command (B7)
    denovo = ("iqtree2", "fasttree", "mafft --localpair", "mafft --genafpair",
              "mafft --maxiterate", "mafft --parttree", "mafft --auto")
    for marker, role in ROUTE_CASES:
        p = subprocess.run([sys.executable,
                            os.path.join(HERE, "recommend_strategy.py"),
                            "--n-seq", "8000", "--median-len", "300",
                            "--marker-class", marker, "--tree-role", role],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        out = p.stdout.decode("utf-8", errors="replace")
        if p.returncode != 0:
            msgs.append("recommend_strategy.py exits %d for %s/%s"
                        % (p.returncode, marker, role))
            continue
        low = out.lower()
        has = [d for d in denovo if d in low]
        if marker == "its" and role in ("skeleton", "placement"):
            if has:
                msgs.append("ITS/%s forbids an ITS tree and then prints %s"
                            % (role, has[0]))
        elif not has:
            msgs.append("%s/%s prints no alignment or tree command at all"
                        % (marker, role))
    return (not msgs), msgs


# --------------------------------------------------- 5. frontmatter is YAML --

def check_frontmatter():
    """Frontmatter must survive a STRICT YAML parser, not just a naive one.

    Why this is its own check: `description` and `compatibility` are free text, and
    a colon inside them ("Chinese trigger words: ...") makes the block invalid YAML
    while a line-based parser happily accepts it.  The universal skills CLI uses a
    real parser and skipped this whole skill with

        YAML parse error: Nested mappings are not allowed in compact mappings

    so the skill was invisible to every client that goes through that tool.  This
    check is dependency-free on purpose: tooling must not need PyYAML installed to
    catch it.
    """
    msgs = []
    path = os.path.join(SKILL, "SKILL.md")
    if not os.path.exists(path):
        return False, ["SKILL.md is missing"]
    text = read(path)
    if not text.startswith("---"):
        return False, ["SKILL.md does not start with YAML frontmatter"]
    parts = text.split("---", 2)
    if len(parts) < 3:
        return False, ["SKILL.md frontmatter is not closed"]
    block = parts[1]

    plain = re.compile(r"^[A-Za-z_][A-Za-z0-9_./-]*$")
    seen = {}
    for lineno, raw in enumerate(block.splitlines(), 1):
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        if raw[0] in " \t":
            continue                       # indented: part of a mapping/list
        if ":" not in raw:
            msgs.append("frontmatter line %d has no key: %r" % (lineno, raw[:60]))
            continue
        key, _, val = raw.partition(":")
        key, val = key.strip(), val.strip()
        if not plain.match(key):
            msgs.append("frontmatter key %r is not a plain scalar" % key)
        seen[key] = val
        if not val:
            continue                       # a mapping header, e.g. 'metadata:'
        if val[0] in "\"'":
            if len(val) < 2 or not val.endswith(val[0]):
                msgs.append("frontmatter %r has an unterminated quoted value" % key)
            continue
        # an unquoted scalar: reject the constructs that silently change meaning
        if ": " in val or val.endswith(":"):
            msgs.append("frontmatter %r contains an unquoted colon (': ') -- quote "
                        "the value, or a strict YAML parser turns the tail into a "
                        "nested mapping" % key)
        if " #" in val:
            msgs.append("frontmatter %r contains an unquoted ' #', which starts a "
                        "comment" % key)
        if val.startswith(("[", "{")) and not val.endswith(("]", "}")):
            msgs.append("frontmatter %r looks like an unclosed flow collection" % key)
    for required in ("name", "description"):
        if not seen.get(required):
            msgs.append("frontmatter is missing a non-empty %r" % required)
    desc = seen.get("description", "").strip("\"'")
    if desc and len(desc) > 1024:
        msgs.append("description is %d chars; the Agent Skills spec allows 1024"
                    % len(desc))
    return (not msgs), msgs


# ------------------------------------------------------------------- main ----

def run_checks(verbose=True):
    """Run every check; return (n_failed, report_lines)."""
    lines = []
    options = all_options()
    params = read(os.path.join(REFS, "parameters.md")) \
        if os.path.exists(os.path.join(REFS, "parameters.md")) else ""
    scaling = read(os.path.join(REFS, "scaling.md")) \
        if os.path.exists(os.path.join(REFS, "scaling.md")) else ""

    checks = [
        ("parameters.md D1-D11 five fields", check_parameters(params)),
        ("scaling.md grid vs recommend_strategy.py", check_scaling(scaling)),
        ("documented --flags vs argparse options", check_flags(options)),
        ("documented script paths exist", check_scripts_exist()),
        ("axis A x axis B routes are distinct", check_routes()),
        ("SKILL.md frontmatter parses as YAML", check_frontmatter()),
    ]
    failed = 0
    for title, (ok, msgs) in checks:
        lines.append("%-44s %s" % (title, "ok" if ok else "FAIL"))
        if not ok:
            failed += 1
            for m in msgs[:25]:
                lines.append("    - %s" % m)
            if len(msgs) > 25:
                lines.append("    ... and %d more" % (len(msgs) - 25))
    return failed, lines


def main():
    ap = argparse.ArgumentParser(
        description="Assert that the documentation and the scripts agree.")
    ap.add_argument("--verbose", action="store_true", default=True)
    args = ap.parse_args()
    print("consistency check -- docs vs scripts")
    print("=" * 74)
    failed, lines = run_checks()
    for l in lines:
        print(l)
    print("=" * 74)
    if failed:
        print("RESULT: FAIL -- %d check group(s) found drift" % failed)
        return 1
    print("RESULT: PASS -- documentation and scripts agree")
    return 0


if __name__ == "__main__":
    sys.exit(main())
