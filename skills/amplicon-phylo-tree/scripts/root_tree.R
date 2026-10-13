#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
## ===========================================================================
## root_tree.R -- root a phylogenetic tree, and show how much of the root
##                placement is a decision rather than a result.
##
## WHY THE SCRIPT EXISTS
##
## A root is not a property of a tree: it is a hypothesis you add.  The same
## branch lengths and the same unrooted topology give different ancestral-state
## reconstructions, different ingroup/outgroup statements and different
## root-to-tip statistics depending on where the root was placed -- and no tree
## file records which choice was made.  So the choice has to be (a) made
## deliberately, (b) written down together with its caveats, and (c) tested for
## how much the answer depends on it.
##
##   * An outgroup is a hypothesis about relationships and therefore evidence.
##     ape::root() is used, and every outgroup label is checked against the tip
##     set first, because a tree and an outgroup list can disagree about
##     quoting, whitespace or a truncated label without anything failing loudly
##     on its own.
##   * Midpoint rooting is a CONVENTION, not evidence.  It assumes an
##     approximately clock-like tree; it is unreliable when the tree spans
##     widely divergent clades, because the longest tip-to-tip path can run
##     between two fast-evolving lineages rather than between the two most
##     distant clades.  It is NOT acceptable evidence-free when the tree feeds
##     root-sensitive metrics (UniFrac, NRI/NTI, betaNTI).  The report says all
##     of this in words that can be adapted into a methods section.
##   * Whether the root is data-decided is measurable: root at several plausible
##     positions and compare the split across the root's two children.  If the
##     split changes, the placement is a decision, and that must be discussed in
##     the methods.  --sensitivity-candidates turns that test on.
##
## Multichotomies are resolved with ape::multi2di, because a tree with
## unresolved nodes has no single answer for "the two children of the root" and
## most downstream tools require a bifurcating tree.  The resolution is
## ARBITRARY: ape resolves multichotomies with ZERO-LENGTH branches (documented
## behaviour of ape::multi2di -- "transforms all multichotomies into a series of
## dichotomies with one (or several) branch(es) of length zero"), so the
## resolution carries no support of its own.  The report states this, and it
## must be stated in the methods.
##
## NOTE on multi2di and unrooted input: multi2di resolves EVERY node with more
## than two children, and ape reads an unrooted tree with a basal trifurcation
## (the root has three children).  On an unrooted input, multi2di therefore also
## converts that basal trifurcation into a bifurcation joined by a zero-length
## branch.  This adds no bipartition and no branch length, but the written
## unrooted file then reads as "root degree 2".  The report prints the root
## degree that was actually written, so nothing is hidden.
##
## NOTE on packages: ape 5.8-1 has NO midpoint() function; phangorn 2.12.1 does.
## Midpoint rooting therefore goes through phangorn::midpoint().  Nothing else
## from phangorn is used.  Nothing is downloaded, nothing is installed, and no
## file other than the requested outputs is touched.
##
## USAGE
##
##   Rscript root_tree.R --in <tree.nwk> --rooted <out_rooted.nwk> \
##                       --unrooted <out_unrooted.nwk> --report <report.txt> \
##                       --label "<study label>" \
##                       [--outgroup <file.txt>] \
##                       [--sensitivity-candidates <n> | --no-sensitivity] [--seed <int>]
##
##   --in                      input Newick tree; exactly one tree
##   --rooted / --unrooted     where to write the rooted / the resolved tree
##   --report                  where to write the plain-text report
##   --label                   study label, used in the report
##   --outgroup <file.txt>     one tip label per line; blank lines and lines
##                             starting with '#' are ignored.  Rooting then uses
##                             ape::root(outgroup = <those tips>, resolve.root =
##                             TRUE)
##   --sensitivity-candidates  n = number of candidate root positions to try
##                             (DEFAULT 2; up to four exist: the primary method,
##                             the midpoint, the longest-branch tip, the longest
##                             interior edge).  This is the test that says whether
##                             the root is decided by the data or by the method,
##                             so it is ON by default; a tree that feeds
##                             root-sensitive metrics must report it.
##   --no-sensitivity          turn the sensitivity test OFF (equivalent to
##                             --sensitivity-candidates 0).  Use it only when the
##                             tree cannot feed root-sensitive metrics, and say so.
##   --seed                    seed for the only stochastic step, multi2di's
##                             random resolution; default 20261012, the same
##                             fixed seed the Python scripts in this folder use,
##                             so even a run without --seed is reproducible.
##   -h / --help               usage
##
## Flags are written `--name value` or `--name=value`.  There are deliberately
## NO positional arguments: positional order is an implicit convention, and a
## silent swap (the rooted tree written over the unrooted one) is invisible in
## the result.
##
## Exit codes: 0 = PASS
##             1 = FAIL (no defensible rooting was produced: a non-monophyletic
##                 outgroup, a failed ape::root()/phangorn::midpoint() call, ...)
##             2 = input or usage error (missing file, unknown flag, a required
##                 flag missing, an outgroup label that is not on the tree)
## ===========================================================================

## The same fixed seed the Python scripts in this folder use, so that a run
## without --seed is still reproducible.
DEFAULT_SEED <- 20261012L

## --------------------------------------------------------------- reporting --

## Every report line goes to stdout as it is produced AND into this buffer, so
## the --report file is written from exactly the lines the user saw.  An
## environment is used because it has reference semantics: say() modifies the
## buffer without copying it.
REPORT <- new.env(parent = emptyenv())
REPORT$lines <- character(0)

say <- function(...) {
  txt <- paste0(...)
  REPORT$lines <- c(REPORT$lines, txt)
  cat(txt, "\n", sep = "")
}

## Progress lines: stdout only, not part of the report.
note <- function(...) cat(paste0("  .. ", ...), "\n", sep = "")

## Every fatal path prints an explicit RESULT line and exits with the code that
## belongs to the class of problem (see the header).
die <- function(..., status = 1L) {
  cat(paste0("ERROR: ", ...), "\n", sep = "")
  cat("RESULT: FAIL\n")
  quit(save = "no", status = status)
}

usage <- function() {
  say("root_tree.R -- root a tree and report how much the root placement is a decision")
  say("")
  say("  Rscript root_tree.R --in <tree.nwk> --rooted <out_rooted.nwk> \\")
  say("                      --unrooted <out_unrooted.nwk> --report <report.txt> \\")
  say("                      --label \"<study label>\" [--outgroup <file.txt>] \\")
  say("                      [--sensitivity-candidates <n> | --no-sensitivity] [--seed <int>]")
  say("")
  say("  --in                      input Newick tree (exactly one tree)")
  say("  --rooted                  output: rooted tree (Newick)")
  say("  --unrooted                output: resolved tree, no root hypothesis applied")
  say("  --report                  output: plain-text report")
  say("  --label                   study label used in the report")
  say("  --outgroup <file.txt>     one tip label per line; '#' starts a comment")
  say("  --sensitivity-candidates  n = candidate root positions to compare")
  say("                            (DEFAULT 2; 0 = off)")
  say("  --no-sensitivity          turn the sensitivity test off (same as 0)")
  say("  --seed <int>              seed for multi2di's random resolution")
  say("                            (default 20261012)")
  say("  -h, --help                this text")
  say("")
  say("Exit codes: 0 = PASS, 1 = FAIL, 2 = input/usage error.")
}

## Collect warnings and errors from an expression instead of letting them
## escape: a warning from ape::root() or phangorn::midpoint() is evidence about
## the run and belongs in the report rather than on stderr only.
collect <- function(expr) {
  ws <- character(0)
  err <- NULL
  val <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      err <<- conditionMessage(e)
      NULL
    }),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(value = val, warnings = ws, error = err)
}

## Write text as UTF-8 with "\n" line endings on every platform.  A text-mode R
## connection would translate "\n" to "\r\n" on Windows, and a report that diffs
## differently per platform is not reproducible.
write_lines_utf8 <- function(path, lines) {
  d <- dirname(path)
  if (nzchar(d) && !dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeLines(enc2utf8(lines), con, sep = "\n", useBytes = TRUE)
}

## Read a small text file (the outgroup list) as UTF-8, dropping the BOM, blank
## lines and '#' comments.
read_lines_utf8 <- function(path) {
  x <- readLines(path, encoding = "UTF-8", warn = FALSE)
  if (length(x)) x[1L] <- sub("^\uFEFF", "", x[1L])
  x <- x[-grep("^\\s*#", x)]
  x <- trimws(x)
  x[nzchar(x)]
}

fmt <- function(x, digits = 6L) {
  if (is.null(x) || length(x) != 1L || is.na(x)) return("NA")
  format(x, digits = digits, scientific = FALSE, trim = TRUE)
}

pct <- function(a, b) if (b > 0) sprintf("%.2f%%", 100 * a / b) else "n/a"

pkg_version <- function(p) {
  v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  if (is.na(v)) "not installed" else v
}

## ------------------------------------------------------------ argument parse --

parse_args <- function(argv) {
  if (length(argv) == 1L && argv %in% c("-h", "--help")) {
    usage()
    quit(save = "no", status = 0L)
  }
  ## Fold "--name=value" into "--name value" so the loop below handles one shape.
  argv2 <- character(0)
  for (a in argv) {
    if (grepl("^--[A-Za-z0-9-]+=", a)) {
      argv2 <- c(argv2, sub("=.*$", "", a), sub("^[^=]*=", "", a))
    } else {
      argv2 <- c(argv2, a)
    }
  }

  o <- list(infile = NA_character_, rooted = NA_character_,
            unrooted = NA_character_, report = NA_character_,
            label = NA_character_, outgroup = NA_character_,
            sensitivity = "2", seed = as.character(DEFAULT_SEED))

  i <- 1L
  while (i <= length(argv2)) {
    a <- argv2[i]
    getv <- function() {
      if (i + 1L > length(argv2) || grepl("^--", argv2[i + 1L]))
        die(sprintf("%s needs a value", a), status = 2L)
      i <<- i + 1L
      argv2[i]
    }
    if (a == "--in") o$infile <- getv()
    else if (a == "--rooted") o$rooted <- getv()
    else if (a == "--unrooted") o$unrooted <- getv()
    else if (a == "--report") o$report <- getv()
    else if (a == "--label") o$label <- getv()
    else if (a == "--outgroup") o$outgroup <- getv()
    else if (a == "--sensitivity-candidates") o$sensitivity <- getv()
    else if (a == "--no-sensitivity") o$sensitivity <- "0"
    else if (a == "--seed") o$seed <- getv()
    else if (a == "-h" || a == "--help") { usage(); quit(save = "no", status = 0L) }
    else die(sprintf("unknown argument: %s (run with --help)", a), status = 2L)
    i <- i + 1L
  }

  for (nm in c("infile", "rooted", "unrooted", "report", "label")) {
    if (is.na(o[[nm]]) || !nzchar(o[[nm]]))
      die(sprintf("--%s is required (run with --help)",
                  if (nm == "infile") "in" else nm), status = 2L)
  }

  ## Integers are parsed strictly: as.integer("x") is NA with a warning, which
  ## would silently become a bad seed or a negative candidate count.
  as_int <- function(txt, flag, min_value) {
    v <- suppressWarnings(as.integer(txt))
    if (is.na(v) || !grepl("^[+-]?[0-9]+$", trimws(txt)))
      die(sprintf("%s expects an integer, got '%s'", flag, txt), status = 2L)
    if (v < min_value)
      die(sprintf("%s must be >= %d, got %d", flag, min_value, v), status = 2L)
    v
  }
  o$sensitivity <- as_int(o$sensitivity, "--sensitivity-candidates", 0L)
  o$seed <- as_int(o$seed, "--seed", 0L)

  if (!file.exists(o$infile))
    die(sprintf("input tree not found: %s", o$infile), status = 2L)
  if (!is.na(o$outgroup) && !file.exists(o$outgroup))
    die(sprintf("outgroup file not found: %s", o$outgroup), status = 2L)
  o
}

## ----------------------------------------------------------- tree utilities --

## The root is the only node that is never a child.  Deriving it from the edge
## matrix avoids assuming that the root is numbered n_tip + 1, which is not
## guaranteed for every function that returns a phylo object.
root_node_of <- function(phy) {
  v <- setdiff(unique(phy$edge[, 1L]), phy$edge[, 2L])
  if (length(v) != 1L) return(NA_integer_)
  v[1L]
}

child_counts <- function(phy) tabulate(phy$edge[, 1L], nbins = max(phy$edge))

root_degree_of <- function(phy) {
  rn <- root_node_of(phy)
  if (is.na(rn)) return(NA_integer_)
  as.integer(sum(phy$edge[, 1L] == rn))
}

## ape's definition, spelled out: a tree is rooted if exactly two edges meet at
## the root, or if it carries a root.edge.
is_rooted_here <- function(phy) {
  if (!is.null(phy$root.edge)) return(TRUE)
  rd <- root_degree_of(phy)
  !is.na(rd) && rd == 2L
}

## Multichotomies, counted the way THIS tree is read.
## A rooted binary tree has at most two children per node.  An UNROOTED binary
## tree is read by ape with a basal trifurcation: the root has three children,
## and that third branch is the unrooted representation, not an unresolved node.
## Counting "more than two children" blindly would call every unrooted binary
## tree unresolved, so the root is treated according to how the tree was read.
multichotomies <- function(phy, rooted) {
  rn <- root_node_of(phy)
  cnt <- child_counts(phy)
  bad <- integer(0)
  for (nd in unique(phy$edge[, 1L])) {
    limit <- 2L
    if (!rooted && !is.na(rn) && nd == rn) limit <- 3L
    if (cnt[nd] > limit) bad <- c(bad, as.integer(nd))
  }
  list(nodes = bad, degrees = if (length(bad)) as.integer(cnt[bad]) else integer(0))
}

## Children as a list, built once: the naive phy$edge[phy$edge[, 1] == nd, 2]
## inside a traversal is O(edges) per node and becomes quadratic on exactly the
## large trees this skill targets.
children_list <- function(phy) {
  ch <- vector("list", max(phy$edge))
  for (i in seq_len(nrow(phy$edge))) {
    p <- phy$edge[i, 1L]
    ch[[p]] <- c(ch[[p]], phy$edge[i, 2L])
  }
  ch
}

tips_below <- function(phy, node, ch, n_tip) {
  seen <- logical(n_tip)
  stack <- node
  while (length(stack)) {
    nd <- stack[length(stack)]
    stack <- stack[-length(stack)]
    if (nd <= n_tip) seen[nd] <- TRUE else stack <- c(stack, ch[[nd]])
  }
  phy$tip.label[seen]
}

## The split across the root's children, as a canonical key.  A bipartition has
## no left or right side, so the two sides are ordered by their sorted label
## lists; two root positions give the same key exactly when they give the same
## split, which counts of tips on each side alone could not decide.
split_key <- function(sides) {
  keys <- vapply(sides, function(s) paste(sort(s), collapse = "|"), character(1))
  min(keys)
}

n_node_labels <- function(phy) {
  if (is.null(phy$node.label)) return(0L)
  as.integer(sum(!is.na(phy$node.label) & nzchar(phy$node.label)))
}

## ------------------------------------------------------------ rooting moves --

root_outgroup <- function(tr, tips) ape::root(tr, outgroup = tips, resolve.root = TRUE)

root_midpoint <- function(tr) phangorn::midpoint(tr, node.labels = "support")

## The tip with the longest terminal branch, and the interior edge with the
## longest branch.  Both are only defined when the tree has branch lengths.
longest_tip <- function(tr) {
  if (is.null(tr$edge.length)) return(NA_character_)
  n_tip <- length(tr$tip.label)
  i <- which(tr$edge[, 2L] <= n_tip)
  if (!length(i)) return(NA_character_)
  tr$tip.label[tr$edge[i[which.max(tr$edge.length[i])], 2L]]
}

longest_interior_child <- function(tr) {
  if (is.null(tr$edge.length)) return(NA_integer_)
  n_tip <- length(tr$tip.label)
  i <- which(tr$edge[, 2L] > n_tip)
  if (!length(i)) return(NA_integer_)
  as.integer(tr$edge[i[which.max(tr$edge.length[i])], 2L])
}

## ===========================================================================
## main
## ===========================================================================

main <- function(argv) {
  args <- parse_args(argv)
  label <- args$label

  say("root_tree.R -- rooting report")
  note("input      : ", args$infile)
  note("label      : ", label)

  ## ------------------------------------------------------------- packages --
  if (!requireNamespace("ape", quietly = TRUE))
    die("the R package 'ape' is not installed (install.packages(\"ape\"))",
        status = 2L)
  ape_ver <- pkg_version("ape")
  phangorn_ok <- requireNamespace("phangorn", quietly = TRUE)
  phangorn_ver <- if (phangorn_ok) pkg_version("phangorn") else "not installed"

  ## ---------------------------------------------------- 1. read the tree ---
  out <- collect(ape::read.tree(args$infile))
  if (!is.null(out$error))
    die(sprintf("ape::read.tree could not read %s: %s", args$infile, out$error),
        status = 2L)
  tr <- out$value
  if (is.null(tr))
    die(sprintf("ape::read.tree returned nothing for %s", args$infile), status = 2L)
  if (inherits(tr, "multiPhylo"))
    die(sprintf(paste0("the input holds more than one tree (%d).  This script ",
                       "roots exactly one tree; split the file and re-run."),
                length(tr)), status = 2L)
  if (!inherits(tr, "phylo"))
    die(sprintf("the input is not a phylo object (class: %s)",
                paste(class(tr), collapse = "/")), status = 2L)

  ## Surrounding quotes: a Newick writer quotes labels that contain characters
  ## with a meaning in Newick, and a quoted label that is not stripped here will
  ## not match the label in an abundance table or an outgroup file -- which is a
  ## silent join failure further downstream, not an error here.
  raw_labels <- tr$tip.label
  clean <- trimws(gsub("^['\"]|['\"]$", "", raw_labels))
  n_quoted <- as.integer(sum(clean != raw_labels))
  tr$tip.label <- clean
  if (anyDuplicated(tr$tip.label))
    die(paste0("the tree has duplicated tip labels after de-quoting; two tips ",
               "cannot be told apart, so the tip set is not usable"),
        status = 2L)

  n_tip <- length(tr$tip.label)
  n_node <- nrow(tr$edge) - n_tip + 1L
  n_edge <- nrow(tr$edge)
  rd_in <- root_degree_of(tr)
  rooted_in <- is_rooted_here(tr)
  labels_before <- n_node_labels(tr)

  has_bl <- !is.null(tr$edge.length)
  n_na_bl <- if (has_bl) as.integer(sum(is.na(tr$edge.length))) else 0L
  bl_min <- bl_max <- NA_real_
  n_zero <- n_neg <- NA_integer_
  if (has_bl) {
    ok <- !is.na(tr$edge.length)
    if (any(ok)) {
      bl_min <- min(tr$edge.length[ok])
      bl_max <- max(tr$edge.length[ok])
    }
    n_zero <- as.integer(sum(ok & tr$edge.length == 0))
    n_neg <- as.integer(sum(ok & tr$edge.length < 0))
  }

  say("")
  say(strrep("=", 78))
  say(sprintf("rooting report -- %s", label))
  say(sprintf("root_tree.R; ape %s, phangorn %s, seed %d",
              ape_ver, phangorn_ver, args$seed))
  say(strrep("=", 78))
  say("")
  say("1  TREE AS READ")
  say(sprintf("  input tree                : %s", args$infile))
  say(sprintf("  tips                      : %d", n_tip))
  say(sprintf("  internal nodes            : %d", n_node))
  say(sprintf("  edges (branches)          : %d", n_edge))
  say(sprintf("  rooted as read            : %s (root degree %s)",
              if (rooted_in) "yes" else "no",
              if (is.na(rd_in)) "?" else as.character(rd_in)))
  say(sprintf("  rooting method requested  : %s",
              if (!is.na(args$outgroup)) "outgroup (ape::root)"
              else "midpoint (phangorn::midpoint)"))
  say(sprintf("  node labels (non-empty)   : %d%s", labels_before,
              if (labels_before > 0)
                "  (support values live on the UNROOTED tree)" else "  (none present)"))
  if (has_bl) {
    say(sprintf("  branch length min / max   : %s / %s", fmt(bl_min), fmt(bl_max)))
    say(sprintf("  zero-length branches      : %d (%s of %d)",
                n_zero, pct(n_zero, n_edge), n_edge))
    if (n_neg > 0)
      say(sprintf("  !! negative branch lengths: %d -- fix the tree before use", n_neg))
    if (n_na_bl > 0)
      say(sprintf("  !! missing (NA) lengths   : %d", n_na_bl))
  } else {
    say("  branch lengths            : ABSENT (the file has no edge lengths)")
    say("    !! without branch lengths, midpoint rooting and the longest-branch")
    say("       candidates are impossible: supply branch lengths, or root with")
    say("       --outgroup")
  }
  say(sprintf("  tip labels de-quoted      : %d", n_quoted))
  if (n_quoted > 0L)
    note("stripped surrounding quotes from ", n_quoted, " tip label(s)")

  ## ------------------------------------------------ 2. resolve multichotomies
  mc_before <- multichotomies(tr, rooted_in)
  note("multichotomies as read: ", length(mc_before$nodes))

  set.seed(args$seed)
  out <- collect(ape::multi2di(tr, random = TRUE))
  if (!is.null(out$error))
    die(sprintf("ape::multi2di failed: %s", out$error), status = 1L)
  res <- out$value

  rd_u <- root_degree_of(res)
  n_node_after <- nrow(res$edge) - length(res$tip.label) + 1L
  added_nodes <- n_node_after - n_node
  added_edges <- nrow(res$edge) - n_edge
  mc_after <- multichotomies(res, is_rooted_here(res))
  zero_after <- if (is.null(res$edge.length)) NA_integer_
                else as.integer(sum(res$edge.length == 0, na.rm = TRUE))

  ## Write the resolved tree.  This is the file that carries the support values:
  ## rooting moves the root, and where a support value ends up after re-rooting
  ## depends on the tool (see ape::root's 'edgelabel' argument and Czech et al.
  ## 2017), so the pre-rooting tree is the safe place to read them.
  out <- collect(ape::write.tree(res, file = args$unrooted))
  if (!is.null(out$error))
    die(sprintf("could not write %s: %s", args$unrooted, out$error), status = 1L)

  say("")
  say("2  MULTICHOTOMY RESOLUTION (declared; arbitrary in detail)")
  say(sprintf("  multichotomies as read    : %d%s", length(mc_before$nodes),
              if (length(mc_before$nodes))
                sprintf(" (child counts: %s)",
                        paste(mc_before$degrees, collapse = ", ")) else ""))
  say(sprintf("  ape::multi2di applied     : yes (random = TRUE, set.seed(%d))",
              args$seed))
  say(sprintf("  nodes / edges added       : %d / %d", added_nodes, added_edges))
  say(sprintf("  multichotomies remaining  : %d", length(mc_after$nodes)))
  if (!is.na(zero_after))
    say(sprintf("  zero-length branches now  : %d", zero_after))
  say("  ape resolves multichotomies with ZERO-LENGTH branches, so the resolution")
  say("  is arbitrary in detail and carries no support of its own.  It must be")
  say("  declared in the methods.")
  say(sprintf("  unrooted file root degree : %s", if (is.na(rd_u)) "?" else as.character(rd_u)))
  say("  For an unrooted tree, root degree 3 is the normal representation.  If it")
  say("  is 2 here, multi2di converted the basal trifurcation of an unrooted")
  say("  input into a bifurcation joined by a zero-length branch: that adds no")
  say("  bipartition and no branch length, but the file must be treated as")
  say("  unrooted by content, not by root degree.")
  say(sprintf("  written                   : %s", args$unrooted))

  ## ------------------------------------------------------------ 3. root ----
  ## The primary rooting is the declared method: the outgroup if one was given,
  ## the midpoint otherwise.
  og_used <- character(0)
  mono_note <- "not applicable"
  if (!is.na(args$outgroup)) {
    primary_name <- "outgroup"
    og <- read_lines_utf8(args$outgroup)
    og <- trimws(gsub("^['\"]|['\"]$", "", og))
    og <- unique(og[nzchar(og)])
    if (!length(og))
      die(sprintf("the outgroup file %s has no usable labels", args$outgroup),
          status = 2L)
    missing <- setdiff(og, res$tip.label)
    if (length(missing))
      die(sprintf(paste0("these outgroup labels are not on the tree: %s.  ",
                         "Check spelling and case, and that the label was not ",
                         "truncated at the first space.  (%d of %d labels ",
                         "missing.)"),
                  paste(utils::head(missing, 10L), collapse = ", "),
                  length(missing), length(og)), status = 2L)
    ## ape's monophyly test is the right one for an unrooted tree: with its
    ## default, reroot = !is.rooted(phy), the group is asked whether it could be
    ## a clade under any rooting, rather than under the root this file happens to
    ## carry.  ape::root() is still the authority -- it is the call that has to
    ## succeed -- so a negative result here is recorded and reported rather than
    ## used to reject a case ape itself would accept.
    mono <- collect(ape::is.monophyletic(res, og))
    og_mono <- isTRUE(mono$value)
    mono_note <- if (!is.null(mono$error))
                   sprintf("could not be tested (%s)", mono$error)
                 else if (og_mono) "yes" else "NO -- the named tips do not form a clade"
    note("outgroup: ", length(og), " tip(s), ape::is.monophyletic = ", mono_note)
    out <- collect(root_outgroup(res, og))
    if (!is.null(out$error))
      die(sprintf(paste0("ape::root(outgroup = ...) failed: %s.  Outgroup ",
                         "monophyly (ape::is.monophyletic): %s.  The named taxa ",
                         "must form a clade on this tree, otherwise they cannot ",
                         "define a root.  Narrow the outgroup to a clade, fix ",
                         "the tree, or root by an explicitly stated ",
                         "alternative."),
                  out$error, mono_note), status = 1L)
    rooted <- out$value
    root_warnings <- out$warnings
    if (!og_mono)
      root_warnings <- c(root_warnings, paste0(
        "ape::is.monophyletic reports that the outgroup is NOT a clade on this ",
        "tree, yet ape::root() returned a tree: the root was then placed at the ",
        "MRCA of the named tips, which is an arbitrary position.  Check the ",
        "rooted tree before using it."))
    og_used <- og
  } else {
    primary_name <- "midpoint"
    if (!phangorn_ok)
      die(paste0("midpoint rooting needs the R package 'phangorn' ",
                 "(phangorn::midpoint); ape 5.8-1 has no midpoint().  ",
                 "install.packages(\"phangorn\"), or pass --outgroup."),
          status = 2L)
    if (!has_bl)
      die(paste0("midpoint rooting needs branch lengths and this tree has none; ",
                 "pass --outgroup, or supply a tree with branch lengths."),
          status = 2L)
    out <- collect(root_midpoint(res))
    if (!is.null(out$error))
      die(sprintf("phangorn::midpoint failed: %s", out$error), status = 1L)
    rooted <- out$value
    root_warnings <- out$warnings
  }
  if (is.null(rooted))
    die("rooting produced no tree", status = 1L)

  out <- collect(ape::write.tree(rooted, file = args$rooted))
  if (!is.null(out$error))
    die(sprintf("could not write %s: %s", args$rooted, out$error), status = 1L)

  rd_r <- root_degree_of(rooted)
  n_node_rooted <- nrow(rooted$edge) - length(rooted$tip.label) + 1L
  labels_after <- n_node_labels(rooted)

  say("")
  say("3  ROOTING")
  say(sprintf("  method used               : %s", primary_name))
  if (primary_name == "outgroup") {
    say(sprintf("  outgroup file             : %s", args$outgroup))
    say(sprintf("  outgroup tips (%d)        : %s", length(og_used),
                if (length(og_used) <= 6) paste(og_used, collapse = ", ")
                else paste(c(utils::head(og_used, 6L), "..."), collapse = ", ")))
    say(sprintf("  outgroup monophyletic     : %s", mono_note))
    say("  call                      : ape::root(res, outgroup = <those tips>,")
    say("                                           resolve.root = TRUE)")
    say("  An outgroup is a hypothesis about relationships, so this root is")
    say("  evidence-based and has to be justified in the methods: why those taxa,")
    say("  and on what reference.")
  } else {
    say("  call                      : phangorn::midpoint(res, node.labels = \"support\")")
    say("  MIDPOINT ROOTING IS A CONVENTION, NOT A HYPOTHESIS.  It places the root")
    say("  at the midpoint of the longest tip-to-tip path, which assumes an")
    say("  approximately clock-like tree.  It is unreliable when the tree spans")
    say("  widely divergent clades, because that longest path can run between two")
    say("  fast-evolving lineages rather than between the two most distant clades.")
    say("  It is NOT acceptable evidence-free when the tree feeds root-sensitive")
    say("  metrics (UniFrac, NRI/NTI, betaNTI): state the assumption, and test how")
    say("  much the result depends on the root (section 4).")
    say("  phangorn::midpoint assumes by default (node.labels = \"support\") that")
    say("  node labels are split support values; that is what was used here.")
  }
  say(sprintf("  rooted tree root degree   : %s",
              if (is.na(rd_r)) "?" else as.character(rd_r)))
  say(sprintf("  rooted tree nodes / edges : %d / %d", n_node_rooted, nrow(rooted$edge)))
  if (!is.na(rd_r) && rd_r != 2L) {
    say("  !! the rooted tree does not have a root of degree 2; check that the")
    say("     rooting call did what you expected before using this tree")
  }
  say(sprintf("  node labels retained      : %d (of %d before rooting)",
              labels_after, labels_before))
  say(sprintf("  tip set unchanged         : %s",
              if (setequal(rooted$tip.label, res$tip.label)) "yes" else "NO (!)"))
  say(sprintf("  written                   : %s", args$rooted))
  if (length(root_warnings))
    for (w in root_warnings) say(sprintf("  warning from the root call: %s", w))
  say("  ape::write.tree rounds branch lengths to 10 significant digits (its")
  say("  default): nothing meaningful is lost, but the written file is not")
  say("  byte-identical to the input.")

  ## ------------------------------------------------- 4. sensitivity test ---
  rooted_sides <- function(phy) {
    n_tip <- length(phy$tip.label)
    rn <- root_node_of(phy)
    if (is.na(rn)) return(NULL)
    kids <- phy$edge[phy$edge[, 1L] == rn, 2L]
    ch <- children_list(phy)
    lapply(kids, function(k) tips_below(phy, k, ch, n_tip))
  }

  primary_sides <- rooted_sides(rooted)
  if (is.null(primary_sides) || length(primary_sides) < 2L)
    die("the rooted tree has no bipartition at the root, so no split can be reported",
        status = 1L)
  primary_key <- split_key(primary_sides)
  primary_counts <- vapply(primary_sides, length, integer(1))
  n_diff <- NA_integer_
  n_cands_used <- 0L

  say("")
  say("4  ROOTING SENSITIVITY")
  cand_paths <- character(0)
  if (args$sensitivity < 1L) {
    say("  NOT RUN (--no-sensitivity / --sensitivity-candidates 0).")
    say("  This is the test that says whether the root you chose is decided by the")
    say("  data or by the method, and it is ON by default (2 candidates).  Leaving")
    say("  it off is a recorded decision: say why in the ops log.  To run it:")
    say(sprintf("    Rscript root_tree.R --in %s --rooted %s --unrooted %s \\",
                args$infile, args$rooted, args$unrooted))
    say(sprintf("                        --report %s --label \"%s\"%s \\",
                args$report, label,
                if (!is.na(args$outgroup)) sprintf(" --outgroup %s", args$outgroup)
                else ""))
    say("                        --sensitivity-candidates 4")
  } else {
    cands <- list()
    cands[[length(cands) + 1L]] <- list(
      name = primary_name, how = sprintf("primary method (%s)", primary_name),
      fn = function() if (primary_name == "outgroup") root_outgroup(res, og_used)
                      else root_midpoint(res))
    if (primary_name != "midpoint" && phangorn_ok && has_bl)
      cands[[length(cands) + 1L]] <- list(
        name = "midpoint", how = "phangorn::midpoint", fn = function() root_midpoint(res))
    lt <- if (has_bl) longest_tip(res) else NA_character_
    if (!is.na(lt))
      cands[[length(cands) + 1L]] <- list(
        name = "longest tip branch",
        how = sprintf("ape::root at tip %s", lt),
        fn = function() root_outgroup(res, lt))
    li <- if (has_bl) longest_interior_child(res) else NA_integer_
    if (!is.na(li))
      cands[[length(cands) + 1L]] <- list(
        name = "longest interior edge",
        how = sprintf("ape::root(node = %d) -- see the note below", li),
        fn = function() ape::root(res, node = li, resolve.root = TRUE))

    if (length(cands) > args$sensitivity)
      cands <- cands[seq_len(args$sensitivity)]
    n_cands_used <- length(cands)
    say(sprintf("  candidates compared       : %d (up to 4 exist), seed %d",
                length(cands), args$seed))
    say(sprintf("  primary root split        : %s tips [%s]",
                paste(primary_counts, collapse = " | "),
                paste(sort(primary_sides[[1L]]), collapse = ",")))
    say("")
    say(sprintf("  %-21s %-13s %-12s %s", "candidate root", "split L|R", "same split?",
                "how the root was placed"))
    say(sprintf("  %s", strrep("-", 98)))
    n_diff <- 0L
    for (ci in seq_along(cands)) {
      cd <- cands[[ci]]
      out <- collect(cd$fn())
      # One file per candidate root, so that "the root moved" can be checked
      # afterwards instead of taken on trust.  The candidate name and how the root
      # was placed are written as a leading comment, because a bare Newick file
      # does not record either.
      cand_path <- sprintf("%s.cand%d.nwk", args$rooted, ci)
      header <- sprintf("# %s -- candidate %d of %d: %s",
                        basename(cand_path), ci, length(cands), cd$name)
      if (!is.null(out$error)) {
        say(sprintf("  %-21s %-13s %-12s %s", cd$name, "-", "-",
                    sprintf("FAILED: %s", out$error)))
        write_lines_utf8(cand_path, c(sprintf("%s -- FAILED: %s", header,
                                              out$error),
                                      sprintf("# %s", cd$how)))
        next
      }
      sides <- rooted_sides(out$value)
      if (is.null(sides) || length(sides) < 2L) {
        say(sprintf("  %-21s %-13s %-12s %s", cd$name, "-", "-",
                    sprintf("root degree %s: not a bipartition",
                            if (is.null(sides)) "?" else length(sides))))
        write_lines_utf8(cand_path, c(sprintf("%s -- NOT a bipartition at the root",
                                              header),
                                      sprintf("# %s", cd$how)))
        next
      }
      same <- identical(split_key(sides), primary_key)
      if (!same) n_diff <- n_diff + 1L
      say(sprintf("  %-21s %-13s %-12s %s", cd$name,
                  paste(vapply(sides, length, integer(1)), collapse = " | "),
                  if (same) "yes" else "NO", cd$how))
      wt <- collect(ape::write.tree(out$value))
      if (is.null(wt$error) && is.character(wt$value) && length(wt$value) == 1L) {
        write_lines_utf8(cand_path, c(
          sprintf("%s -- split %s the primary root", header,
                  if (same) "MATCHES" else "DIFFERS FROM"),
          sprintf("# how the root was placed: %s", cd$how),
          wt$value))
        cand_paths <- c(cand_paths, cand_path)
        note("candidate root written: ", cand_path)
      } else {
        say(sprintf("  !! candidate %d tree could not be written: %s", ci,
                    if (is.null(wt$error)) "unexpected value" else wt$error))
      }
    }
    say("")
    if (n_diff > 0L) {
      say(sprintf("  %d of %d candidate root positions give a DIFFERENT split across",
                  n_diff, length(cands)))
      say("  the root.  A differing split means the root placement is NOT decided")
      say("  by the data: two defensible placements put different taxa on")
      say("  different sides of the root.  Anything that depends on the root")
      say("  (ancestral states, statements about the ingroup, root-to-tip")
      say("  statistics) inherits that ambiguity, and the methods must say so.")
    } else {
      say(sprintf("  All %d candidate root positions give the SAME split across the",
                  length(cands)))
      say("  root.  That is evidence that this bipartition is robust to the root")
      say("  placement; it does not turn midpoint rooting into a hypothesis, and")
      say("  the assumption behind the method still has to be stated.")
    }
    say("  Each candidate is one plausible convention, not an exhaustive search:")
    say("  the point is to show whether the answer moves, not to find a best root.")
    say("  Note on ape::root(node = ...): when a node instead of an outgroup is")
    say("  given, ape has to pick which side is the ingroup, and it resolves that")
    say("  ambiguity arbitrarily (the clade on the right of the node in the")
    say("  default plot).  The longest-interior-edge candidate therefore also")
    say("  illustrates how much of a 'rooting' can be an implementation detail.")
  }

  ## ---------------------------------------------------- 5. methods text ----
  say("")
  say("5  STATE THIS IN THE METHODS (adapt the wording, keep the statements)")
  if (primary_name == "outgroup") {
    say(sprintf(paste0("  The tree was rooted on the outgroup %s, which was ",
                       "required to be monophyletic on the tree (checked with ",
                       "ape::is.monophyletic); rooting used ape::root(..., ",
                       "resolve.root = TRUE)."),
                if (length(og_used) <= 4) paste(og_used, collapse = ", ")
                else sprintf("%d taxa", length(og_used))))
  } else {
    say(paste0("  The tree was rooted with the midpoint criterion ",
               "(phangorn::midpoint).  Midpoint rooting is a convention, not a ",
               "hypothesis: it assumes that evolutionary rates are roughly ",
               "clock-like across the whole tree.  That assumption was not ",
               "tested here and may not hold across the clades in this dataset."))
  }
  say(paste0("  Multichotomies were resolved with ape::multi2di before rooting. ",
             "ape resolves multichotomies with zero-length branches, so the ",
             "resolution is arbitrary in detail and carries no support of its own."))
  if (args$sensitivity >= 1L && !is.na(n_diff) && n_cands_used >= 2L) {
    if (n_diff > 0L)
      say(paste0("  Rooting sensitivity was tested by re-rooting at ",
                 n_cands_used, " candidate positions (outgroup / midpoint / ",
                 "longest terminal branch / longest interior branch); ",
                 n_diff, " of them changed the bipartition across the root.  ",
                 "The root placement is therefore not decided by the data, and ",
                 "results that depend on the root are reported with that ",
                 "limitation."))
    else
      say(paste0("  Rooting sensitivity was tested by re-rooting at ",
                 n_cands_used, " candidate positions (outgroup / midpoint / ",
                 "longest terminal branch / longest interior branch); all of ",
                 "them produced the same bipartition across the root."))
  } else {
    say(paste0("  The root placement was not tested for sensitivity to ",
               "alternative positions.  If this tree feeds root-sensitive ",
               "metrics (UniFrac, NRI/NTI, betaNTI), run that test with ",
               "--sensitivity-candidates and report the outcome."))
  }
  say(paste0("  Support values are read from the UNROOTED tree (",
             args$unrooted, "): re-rooting moves where a support value sits, ",
             "and how a toolkit re-attaches it differs between tools."))
  say(paste0("  Software: R ", paste(R.version$major, R.version$minor, sep = "."),
             " with ape ", ape_ver, " and phangorn ", phangorn_ver,
             "; a multi2di seed of ", args$seed, " was set for reproducibility."))

  ## ------------------------------------------------------------ 6. files ---
  say("")
  say("6  OUTPUT FILES")
  say(sprintf("  rooted tree               : %s", args$rooted))
  say(sprintf("  resolved unrooted tree    : %s", args$unrooted))
  say(sprintf("  report                    : %s", args$report))
  if (length(cand_paths)) {
    say(sprintf("  candidate root trees (%d)  :", length(cand_paths)))
    for (p in cand_paths) say(sprintf("    %s", p))
    say("  Each candidate file carries a leading comment with its name and how the")
    say("  root was placed.  Use them to re-check the sensitivity claim, or to")
    say("  compute the downstream metrics (UniFrac / betaNTI) yourself on the SAME")
    say("  candidate trees: this skill reports the root split only, it does not")
    say("  compute those metrics for you.")
  }

  say("")
  say("RESULT: PASS")

  ## The report file is written last, from exactly the lines that were printed.
  write_lines_utf8(args$report, REPORT$lines)
  invisible(0L)
}

## An unexpected error must still exit non-zero with an explicit RESULT line, and
## must not leave a half-written report behind claiming a run happened.
tryCatch(
  main(commandArgs(trailingOnly = TRUE)),
  error = function(e) {
    cat(sprintf("ERROR: unexpected failure: %s\n", conditionMessage(e)))
    cat("RESULT: FAIL\n")
    quit(save = "no", status = 1L)
  })
