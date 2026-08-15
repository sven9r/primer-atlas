#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(readr))

primers <- read_csv("data/primers.csv", show_col_types = FALSE)
pair_ids <- unique(primers$pair_id)
overall <- read_csv(
  "data/derived/claimed_primer_overall_summary.csv",
  show_col_types = FALSE
)

pair_templates <- nanoparquet::read_parquet(
  "data/pinned/COI/2026-08-14/pair_templates.parquet"
)

stopifnot(
  length(pair_ids) == 25L,
  setequal(overall$pair_id, pair_ids),
  nrow(overall) == length(pair_ids),
  setequal(unique(pair_templates$pair_id), pair_ids),
  all(overall$n_centroids == 67352L),
  all(overall$pair_scorable_fraction >= 0),
  all(overall$pair_scorable_fraction <= 1)
)

required_score_columns <- c(
  "accession",
  "pair_id",
  "pair_label",
  "order",
  "family",
  "genus",
  "pair_scorable",
  "pair_penalty",
  "forward_binding_sequence_primer_oriented",
  "reverse_binding_sequence_primer_oriented"
)
required_pair_template_columns <- setdiff(
  required_score_columns,
  c("accession", "family", "genus")
)
stopifnot(all(required_pair_template_columns %in% names(pair_templates)))

app_source <- paste(readLines("app.R", warn = FALSE), collapse = "\n")
function_source <- paste(
  readLines("R/functions.R", warn = FALSE),
  collapse = "\n"
)
builder_source <- paste(
  readLines("scripts/build_claimed_primer_scores.R", warn = FALSE),
  collapse = "\n"
)

stopifnot(
  grepl(
    "input.detail_pair && input.detail_pair != 'BEEPRIME'",
    app_source,
    fixed = TRUE
  ),
  grepl(
    "score_primer_pair_expanded_reference(",
    app_source,
    fixed = TRUE
  ),
  grepl(
    "expanded_scores_path",
    app_source,
    fixed = TRUE
  ),
  grepl(
    "score_primer_pair_expanded_reference <- function",
    function_source,
    fixed = TRUE
  ),
  grepl(
    "default_pairs <- unique(primers$pair_id)",
    builder_source,
    fixed = TRUE
  )
)

message(
  "Universal curated and session-custom exact-sequence drill-down checks passed."
)
