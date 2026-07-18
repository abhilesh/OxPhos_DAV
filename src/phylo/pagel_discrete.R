#!/usr/bin/env Rscript
# src/phylo/pagel_discrete.R
#
# Pagel's discrete test for correlated binary character evolution.
# Called by 01_find_compensating_partners.py.
#
# Three modes:
#
# Chunk mode (2 args) — for HPC job arrays:
#   Rscript pagel_discrete.R <chunk.tsv> <out_file>
#   chunk.tsv: manifest TSV (identity_key  tree_nwk  trait1_txt  trait2_txt)
#              paths are RELATIVE to the working directory (cd to pagel_jobs/ first)
#   out_file:  results written here (not stdout)
#   Output file format: identity_key\tpagel_p\t<p_value_or_NA>
#   Checkpointing: if out_file already exists, already-processed pairs are
#   skipped and new results are appended. Safe to resubmit cancelled jobs.
#
# Batch mode (1 arg = manifest TSV) — for local single-session run:
#   Rscript pagel_discrete.R <manifest.tsv>
#   Output to stdout: pair_id\tpagel_p\t<p_value_or_NA>
#
# Single mode (3 args):
#   Rscript pagel_discrete.R <tree.nwk> <trait1.txt> <trait2.txt>
#   Output: pagel_p\t<p_value>
#
# Input trait files:
#   species<TAB>0/1  (one per line, no header)

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
})

# Per-pair wall-time limit. Normal pairs finish in ~60-90s on HPC; this
# catches genuinely non-converging fitPagel calls without killing normal ones.
PAGEL_TIMEOUT_SEC <- 300L

# ── In-memory tree cache keyed by file path ───────────────────────────────────
.tree_cache <- new.env(hash = TRUE, parent = emptyenv())

get_tree <- function(tree_file) {
  if (exists(tree_file, envir = .tree_cache)) {
    return(get(tree_file, envir = .tree_cache))
  }
  tree <- read.tree(tree_file)
  # Force ultrametric: floating-point drift in MamPhy branch lengths causes
  # is.ultrametric() to return FALSE; phytools::fitPagel requires ultrametric.
  if (!is.ultrametric(tree)) {
    tree <- suppressMessages(force.ultrametric(tree, method = "extend"))
  }
  assign(tree_file, tree, envir = .tree_cache)
  tree
}

read_trait <- function(path) {
  d <- read.table(path, header = FALSE, sep = "\t",
                  col.names = c("species", "state"),
                  stringsAsFactors = FALSE)
  v <- as.integer(d$state)
  names(v) <- d$species
  v
}

run_pagel <- function(tree_file, trait1_file, trait2_file) {
  tryCatch({
    tree <- get_tree(tree_file)
    if (is.null(tree) || length(tree$tip.label) < 4) return(NA_real_)

    t1_raw <- read_trait(trait1_file)
    t2_raw <- read_trait(trait2_file)

    common <- intersect(tree$tip.label, intersect(names(t1_raw), names(t2_raw)))
    if (length(common) < 4) return(NA_real_)

    tree_use <- keep.tip(tree, common)
    t1 <- t1_raw[common]
    t2 <- t2_raw[common]

    if (length(unique(t1)) < 2 || length(unique(t2)) < 2) return(NA_real_)

    setTimeLimit(elapsed = PAGEL_TIMEOUT_SEC, transient = TRUE)
    fit <- fitPagel(tree_use, t1, t2)
    setTimeLimit(elapsed = Inf, transient = TRUE)

    p_val <- fit$P
    if (is.na(p_val) || !is.finite(p_val)) return(NA_real_)
    return(p_val)

  }, error = function(e) {
    setTimeLimit(elapsed = Inf, transient = TRUE)
    NA_real_
  })
}

# ── Checkpoint helpers ────────────────────────────────────────────────────────

# Read already-completed pair_ids from an existing output file.
read_done <- function(out_file) {
  if (!file.exists(out_file)) return(character(0))
  d <- tryCatch(
    read.table(out_file, header = FALSE, sep = "\t",
               col.names = c("pair_id", "field", "value"),
               stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(d) || nrow(d) == 0) return(character(0))
  unique(d$pair_id)
}

run_batch <- function(manifest, output_con, report_every = 50L,
                      done_ids = character(0)) {
  n_total  <- nrow(manifest)
  n_done   <- 0L
  n_skip   <- 0L
  t_start  <- proc.time()[["elapsed"]]

  for (i in seq_len(n_total)) {
    row <- manifest[i, ]

    if (row$pair_id %in% done_ids) {
      n_skip <- n_skip + 1L
      next
    }

    p_val <- run_pagel(row$tree_file, row$trait1_file, row$trait2_file)
    p_str <- if (is.na(p_val)) "NA" else sprintf("%.6g", p_val)
    cat(sprintf("%s\tpagel_p\t%s\n", row$pair_id, p_str), file = output_con)
    flush(output_con)

    n_done <- n_done + 1L
    if (n_done %% report_every == 0L) {
      elapsed <- proc.time()[["elapsed"]] - t_start
      rate    <- if (elapsed > 0) n_done / elapsed else NA_real_
      remain  <- if (!is.na(rate) && rate > 0) {
        (n_total - n_skip - n_done) / rate
      } else NA_real_
      message(sprintf(
        "  Pagel: %d / %d  (%.2f s/pair, ~%.0f s remaining)",
        n_done + n_skip, n_total,
        if (!is.na(rate) && rate > 0) 1 / rate else NA_real_,
        if (!is.na(remain)) remain else Inf))
      flush(stderr())
    }
  }

  if (n_skip > 0L)
    message(sprintf("  Skipped (already done): %d", n_skip))
}

# ── Dispatch ──────────────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 2) {
  # ── Chunk mode (HPC job array) with checkpointing ────────────────────────────
  chunk <- read.table(args[1], header = FALSE, sep = "\t",
                      col.names = c("pair_id", "tree_file",
                                    "trait1_file", "trait2_file"),
                      stringsAsFactors = FALSE)

  done_ids <- read_done(args[2])
  n_remaining <- nrow(chunk) - length(done_ids)
  message(sprintf("Chunk: %s  (%d pairs, %d already done, %d remaining) → %s",
                  args[1], nrow(chunk), length(done_ids), n_remaining, args[2]))

  out_con <- file(args[2], open = "a")   # append mode — safe for checkpointing
  run_batch(chunk, out_con, report_every = 50L, done_ids = done_ids)
  close(out_con)

} else if (length(args) == 1) {
  # ── Batch mode (local, stdout) ────────────────────────────────────────────────
  manifest <- read.table(args[1], header = FALSE, sep = "\t",
                         col.names = c("pair_id", "tree_file",
                                       "trait1_file", "trait2_file"),
                         stringsAsFactors = FALSE)
  run_batch(manifest, stdout(), report_every = 500L)

} else if (length(args) == 3) {
  # ── Single mode ──────────────────────────────────────────────────────────────
  p_val <- run_pagel(args[1], args[2], args[3])
  p_str <- if (is.na(p_val)) "NA" else sprintf("%.6g", p_val)
  cat(sprintf("pagel_p\t%s\n", p_str))

} else {
  cat("pagel_p\tNA\n")
}