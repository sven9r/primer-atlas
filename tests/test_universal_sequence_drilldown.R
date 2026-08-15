#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(readr))

primers <- read_csv("data/primers.csv", show_col_types = FALSE)
pair_ids <- unique(primers$pair_id)
overall <- read_csv(
  "data/derived/claimed_primer_overall_summary.csv",
  show_col_types = FALSE
)

score_paths <- file.path(
  "data",
  "derived",
  paste0("claimed_", tolower(pair_ids), "_centroid_scores.csv.gz")
)

stopifnot(
  length(pair_ids) == 25L,
  setequal(overall$pair_id, pair_ids),
  nrow(overall) == length(pair_ids),
  all(file.exists(score_paths)),
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
for (path in score_paths) {
  header <- names(read.csv(
    gzfile(path),
    nrows = 1L,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ))
  stopifnot(all(required_score_columns %in% header))
}

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
