# Installing this skill

Every client that supports Agent Skills reads a `SKILL.md` inside its own skills
directory. This file lists those directories per client, plus the details that are
specific to DeepSeek Harness (DSH).

The short version: `.agents/skills/` is read by most clients, so putting the folder
there once covers them all.

```bash
npx skills add mumu-sci/amplicon-phylo-tree@amplicon-phylo-tree -g
```

---

## Where each client looks

| Client | Project-level | User-level |
|---|---|---|
| **OpenCode** | `.opencode/skills/` or `.agents/skills/` | `~/.config/opencode/skills/` or `~/.agents/skills/` |
| **Claude Code** | `.claude/skills/` | `~/.claude/skills/` |
| **Codex** | `.agents/skills/` | `~/.codex/skills/` |
| **Cursor** | `.agents/skills/` | `~/.cursor/skills/` |
| **Cline** | `.agents/skills/` | `~/.agents/skills/` |
| **Windsurf** | `.windsurf/skills/` | `~/.codeium/windsurf/skills/` |
| **DSH** | `.dsh/skills/` or `.agents/skills/` | `~/.dsh/skills/` or `~/.agents/skills/` |

`.agents/skills/` is the shared location: OpenCode, Codex, Cursor, Cline and DSH all
read it.

---

## The installer

Both scripts ship at the repository root and do the same work:

```bash
bash install.sh                 # auto-detect a skills root
bash install.sh --project       # ./.agents/skills (project-level)
bash install.sh --dest <root>   # an explicit skills root
bash install.sh --dry-run
```

```powershell
pwsh -File install.ps1
pwsh -File install.ps1 -Project
pwsh -File install.ps1 -Dest <root>
```

They pick the target root, refuse to overwrite a copy you have edited, strip
`__pycache__` (a `.pyc` embeds the absolute path of the machine that created it),
and then run the verification commands on what they installed.

---

## DeepSeek Harness (DSH)

DSH has no skill marketplace: its loader scans a fixed list of directories, in
priority order, and the first match for a name wins.

| Rank | Source | Path |
|---|---|---|
| 100 | project `.dsh` | `<project root>/.dsh/skills` |
| 200 | project `.agents` | `<project root>/.agents/skills` |
| 300 | custom | your `customSkillDirs` configuration |
| 400 | user `.dsh` | `~/.dsh/skills` |
| 500 | user `.agents` | `~/.agents/skills` (respects `DSH_AGENTS_HOME`) |
| 600 | bundled | ships with the client (`bundledSkillDir` / `DSH_BUNDLED_SKILL_DIR`) |

A project root is the nearest ancestor containing `.git`. Use ranks 100-200 to pin a
version to one project, and 400-500 to make the skill available everywhere.

Additional DSH details worth knowing:

- A skill is `<name>/SKILL.md` **at the top level of a scanned root**; nested
  `**/SKILL.md` files are deliberately not discovered. The folder name must equal
  the `name` in the frontmatter.
- DSH accepts the optional frontmatter keys `whenToUse`, `metadata`,
  `disable-model-invocation` and `user-invocable`; unknown keys are ignored.
- Rank 600 is the only way a skill reaches every DSH user without them doing
  anything, and that list lives inside the client distribution, so a third party
  cannot add to it. If you think this skill belongs there, open an issue here.
