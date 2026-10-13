<#
.SYNOPSIS
  Install this skill into a DeepSeek Harness (or any Agent Skills client) skills root.

.DESCRIPTION
  A client discovers skills by scanning fixed directories for <name>/SKILL.md, so
  installing one means putting the folder in the right root.  This script does
  that, checks what it is about to overwrite, strips byte-code that would leak the
  installing machine's paths, and verifies the installed copy.

.PARAMETER Dest
  The skills ROOT (the parent of the skill folder), e.g. ~/.agents/skills.

.PARAMETER Project
  Install into ./.agents/skills instead (higher priority than the user roots).

.PARAMETER Url
  Clone from this git URL instead of using the folder this script lives in.

.PARAMETER DryRun
  Show what would happen and exit.

.EXAMPLE
  pwsh -File install.ps1
.EXAMPLE
  pwsh -File install.ps1 -Url https://github.com/mumu-sci/amplicon-phylo-tree
#>
param(
  [string]$Dest = "",
  [switch]$Project,
  [string]$Url = "",
  [switch]$DryRun
)

$ErrorActionPreference = "Stop"
$SkillName = "amplicon-phylo-tree"

# ------------------------------------------------------------- source folder --
# This script ships at the repository root, so the skill folder can be either
# beside it (someone copied the skill out on its own) or under skills/<name>
# (the repository layout).  Look in both, in that order.
$Here = Split-Path -Parent $MyInvocation.MyCommand.Path
$Src = ""
foreach ($cand in @($Here, (Join-Path $Here "skills/$SkillName"))) {
  if (Test-Path (Join-Path $cand "SKILL.md")) { $Src = $cand; break }
}
if (-not $Src -and -not $Url) {
  Write-Error "no SKILL.md beside this script ($Here) or under skills/$SkillName, and no -Url given"
  exit 1
}

# ---------------------------------------------------------------- pick a root --
if ($Dest) {
  $Root = $Dest
} elseif ($Project) {
  $Root = Join-Path (Get-Location) ".agents/skills"
} else {
  # Mirrors the loader's defaults: DSH_AGENTS_HOME (or ~/.agents), falling back to
  # the user .dsh root when that is where skills already live.
  $AgentsHome = if ($env:DSH_AGENTS_HOME) { $env:DSH_AGENTS_HOME }
                else { Join-Path $env:USERPROFILE ".agents" }
  if (Test-Path (Join-Path $AgentsHome "skills")) {
    $Root = Join-Path $AgentsHome "skills"
  } elseif (Test-Path (Join-Path $env:USERPROFILE ".dsh/skills")) {
    $Root = Join-Path $env:USERPROFILE ".dsh/skills"
  } else {
    $Root = Join-Path $AgentsHome "skills"
  }
}

$Target = Join-Path $Root $SkillName
Write-Host "target root : $Root"

if ($DryRun) {
  Write-Host "DRY RUN: would install"
  if ($Url) { Write-Host "  from $Url" } else { Write-Host "  from $Src" }
  Write-Host "  into $Target"
  exit 0
}

# --------------------------------------------------------------- install it ---
New-Item -ItemType Directory -Force -Path $Root | Out-Null

$Tmp = $null
if ($Url) {
  if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
    Write-Error "git is not on PATH (needed for -Url)"; exit 1
  }
  $Tmp = Join-Path ([IO.Path]::GetTempPath()) ("aph-" + [guid]::NewGuid().ToString("N"))
  Write-Host "cloning     : $Url"
  git clone --depth 1 --quiet $Url $Tmp
  if ($LASTEXITCODE -ne 0) { Write-Error "git clone failed"; exit 1 }
  if (Test-Path (Join-Path $Tmp "SKILL.md")) {
    $Src = $Tmp
  } elseif (Test-Path (Join-Path $Tmp "skills/$SkillName/SKILL.md")) {
    $Src = Join-Path $Tmp "skills/$SkillName"
  } else {
    Write-Error "no SKILL.md found in $Url"
    exit 1
  }
}

try {
  # Refuse to overwrite a DIFFERENT skill that happens to share the name.
  $targetSkill = Join-Path $Target "SKILL.md"
  if (Test-Path $targetSkill) {
    $a = (Get-FileHash (Join-Path $Src "SKILL.md") -Algorithm SHA256).Hash
    $b = (Get-FileHash $targetSkill -Algorithm SHA256).Hash
    if ($a -ne $b) {
      Write-Host "existing install differs from this source:"
      Write-Host "  $targetSkill"
      Write-Host "Refusing to overwrite a skill you may have edited.  Move it aside first:"
      Write-Host "  Move-Item '$Target' '$Target.backup-'+(Get-Date -Format yyyyMMddHHmmss)"
      exit 1
    }
  }

  if (Test-Path $Target) { Remove-Item $Target -Recurse -Force }
  New-Item -ItemType Directory -Force -Path $Target | Out-Null
  Copy-Item (Join-Path $Src "*") -Destination $Target -Recurse -Force

  Write-Host "installed   : $Target"

  # ---------------------------------------------------------------- verify ----
  $Py = (Get-Command python -ErrorAction SilentlyContinue).Source
  if (-not $Py) { $Py = (Get-Command python3 -ErrorAction SilentlyContinue).Source }
  if ($Py) {
    Write-Host ""
    Write-Host "verifying the installed copy:"
    foreach ($script in @("selftest.py", "check_consistency.py", "regression_tests.py")) {
      Push-Location $Target
      try { & $Py "scripts/$script" | Select-Object -Last 1 }
      finally { Pop-Location }
    }
  } else {
    Write-Host "note: no python on PATH, skipped the verification step"
  }

  # Strip byte-code AFTER verification: running the scripts is what creates it, and
  # a .pyc embeds the absolute path of the machine that made it.
  Get-ChildItem $Target -Recurse -Directory -Filter "__pycache__" -ErrorAction SilentlyContinue |
    ForEach-Object { Remove-Item $_.FullName -Recurse -Force -ErrorAction SilentlyContinue }
  Get-ChildItem $Target -Recurse -File -Filter "*.pyc" -ErrorAction SilentlyContinue |
    ForEach-Object { Remove-Item $_.FullName -Force -ErrorAction SilentlyContinue }

  Write-Host ""
  Write-Host "Done. The skill is discovered on the next catalog scan (reload or restart"
  Write-Host "the client). The folder must keep its name '$SkillName'."
} finally {
  if ($Tmp -and (Test-Path $Tmp)) { Remove-Item $Tmp -Recurse -Force }
}
