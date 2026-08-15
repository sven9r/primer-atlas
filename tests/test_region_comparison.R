#!/usr/bin/env Rscript

source("R/functions.R")

sequences <- data.frame(
  country_or_territory = c("A", "A", "B", "B"),
  forward_binding_sequence_primer_oriented = c("ACGT", "ACGT", "ATGT", "ATGT"),
  reverse_binding_sequence_primer_oriented = c("TGCA", "TGCA", "TGTA", "TGTA"),
  stringsAsFactors = FALSE
)

a <- reference_binding_profile(sequences, "A", "forward", "ACGT")
b <- reference_binding_profile(sequences, "B", "forward", "ACGT")

stopifnot(
  nrow(a) == 4L,
  nrow(b) == 4L,
  a$consensus[2] == "C",
  b$consensus[2] == "T",
  a$incompatible_fraction[2] == 0,
  b$incompatible_fraction[2] == 1,
  all(a$n_total == 2L),
  all(b$n_available == 2L)
)

app_source <- paste(readLines("app.R", warn = FALSE), collapse = "\n")
stopifnot(
  grepl("reference sequences—not regional amplification probability", app_source, fixed = TRUE),
  grepl("summarize_region_orders", app_source, fixed = TRUE),
  grepl("reference_binding_profile", app_source, fixed = TRUE),
  grepl('DTOutput("region_sequences")', app_source, fixed = TRUE)
)

message("Region comparison sequence-profile and interpretation checks passed.")
