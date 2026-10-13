#!/usr/bin/env bash
# Install this skill into a DeepSeek Harness (or any Agent Skills client) skills root.
#
# Why a script: DSH discovers skills by scanning fixed directories for
# <name>/SKILL.md -- there is no marketplace to publish to.  So "distributing" a
# skill means "getting the folder into the right root", and this is that step,
# spelled out and checked instead of described in prose.
#
# Usage:
#   bash install.sh                     # auto-detect the target root
#   bash install.sh --dest <dir>        # explicit root (must be the PARENT of the
#                                       # skill folder, e.g. ~/.agents/skills)
#   bash install.sh --project           # install into ./.agents/skills (rank 200,
#                                       # highest priority after .dsh/skills)
#   bash install.sh --url <git-url>     # clone from a repository instead of using
#                                       # the folder this script lives in
#   bash install.sh --dry-run
#
# Exit codes: 0 = installed (or already correct), 1 = failure, 2 = usage error.
set -euo pipefail

SKILL_NAME="amplicon-phylo-tree"
MODE="auto"
DEST=""
URL=""
DRY=0

while [ $# -gt 0 ]; do
  case "$1" in
    --dest)    DEST="${2:-}"; shift 2 ;;
    --project) MODE="project"; shift ;;
    --url)     URL="${2:-}"; shift 2 ;;
    --dry-run) DRY=1; shift ;;
    -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1 (see --help)" >&2; exit 2 ;;
  esac
done

say() { printf '%s\n' "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

# ---------------------------------------------------------------- pick a root --
if [ -n "$DEST" ]; then
  ROOT="$DEST"
elif [ "$MODE" = "project" ]; then
  ROOT="$PWD/.agents/skills"
else
  # Mirrors the loader's own defaults: DSH_AGENTS_HOME (or ~/.agents), with the
  # user .dsh root as a fallback if that is where skills already live.
  AGENTS_HOME="${DSH_AGENTS_HOME:-$HOME/.agents}"
  if [ -d "$AGENTS_HOME/skills" ]; then
    ROOT="$AGENTS_HOME/skills"
  elif [ -d "$HOME/.dsh/skills" ]; then
    ROOT="$HOME/.dsh/skills"
  else
    ROOT="$AGENTS_HOME/skills"
  fi
fi

say "target root : $ROOT"
TARGET="$ROOT/$SKILL_NAME"

# ------------------------------------------------------------- source folder --
# This script ships at the repository root, so the skill folder can be either
# beside it (someone copied the skill out on its own) or under skills/<name>
# (the repository layout).  Look in both, in that order.
HERE="$(cd "$(dirname "$0")" && pwd)"
SRC=""
if [ -f "$HERE/SKILL.md" ]; then
  SRC="$HERE"
elif [ -f "$HERE/skills/$SKILL_NAME/SKILL.md" ]; then
  SRC="$HERE/skills/$SKILL_NAME"
elif [ -z "$URL" ]; then
  die "no SKILL.md beside this script ($HERE) or under skills/$SKILL_NAME, and no --url given"
fi

if [ "$DRY" -eq 1 ]; then
  say "DRY RUN: would install"
  if [ -n "$URL" ]; then
    say "  from $URL"
  else
    say "  from $SRC"
  fi
  say "  into $TARGET"
  exit 0
fi

mkdir -p "$ROOT"

# --------------------------------------------------------------- install it ---
if [ -n "$URL" ]; then
  command -v git >/dev/null || die "git is not on PATH (needed for --url)"
  TMP="$(mktemp -d)"
  trap 'rm -rf "$TMP"' EXIT
  say "cloning     : $URL"
  git clone --depth 1 --quiet "$URL" "$TMP/repo"
  # the skill may be at the repository root or under skills/
  if [ -f "$TMP/repo/SKILL.md" ]; then
    SRC="$TMP/repo"
  elif [ -f "$TMP/repo/skills/$SKILL_NAME/SKILL.md" ]; then
    SRC="$TMP/repo/skills/$SKILL_NAME"
  else
    die "no SKILL.md found in $URL (looked at the root and at skills/$SKILL_NAME)"
  fi
fi

# Refuse to overwrite a DIFFERENT skill that happens to share the name.
if [ -e "$TARGET/SKILL.md" ] && ! cmp -s "$SRC/SKILL.md" "$TARGET/SKILL.md"; then
  say "existing install differs from this source:"
  say "  $TARGET/SKILL.md"
  say "Refusing to overwrite a skill you may have edited.  Move it aside first:"
  say "  mv '$TARGET' '$TARGET.backup-$(date +%Y%m%d%H%M%S)'"
  exit 1
fi

rm -rf "$TARGET"
mkdir -p "$TARGET"
# copy the skill folder's contents (cp -R of "." keeps this portable to BSD cp too)
(cd "$SRC" && tar cf - .) | (cd "$TARGET" && tar xf -)

say "installed   : $TARGET"

# ------------------------------------------------------------------ verify ----
PY="$(command -v python3 || command -v python || true)"
if [ -n "$PY" ]; then
  say ""
  say "verifying the installed copy:"
  ( cd "$TARGET" && "$PY" scripts/selftest.py | tail -1 )
  ( cd "$TARGET" && "$PY" scripts/check_consistency.py | tail -1 )
  ( cd "$TARGET" && "$PY" scripts/regression_tests.py | tail -1 )
else
  say "note: no python3 on PATH, skipped the verification step"
fi

# Strip byte-code AFTER verification: running the scripts is what creates it, and a
# .pyc embeds the absolute path of the machine that made it.
find "$TARGET" -name '__pycache__' -type d -prune -exec rm -rf {} + 2>/dev/null || true
find "$TARGET" -name '*.pyc' -delete 2>/dev/null || true

say ""
say "Done. The skill is discovered on the next catalog scan (reload or restart the"
say "client). The folder must keep its name '$SKILL_NAME'."
