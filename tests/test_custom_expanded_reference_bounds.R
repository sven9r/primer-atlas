#!/usr/bin/env Rscript

.libPaths(c(".Rlib", .libPaths()))
suppressPackageStartupMessages(library(readr))
source("R/functions.R")

forward_sequence <- "CWAATCAYARRGATATTGG"
reverse_sequence <- "CHACATAATAHGTRTCATG"
reference <- clean_sequence(
  parse_fasta("data/reference/coi_reference_NC_001322.1.fasta")[[1]]
)
located <- locate_primer_pair(
  reference_sequence = reference,
  forward_sequence = forward_sequence,
  reverse_sequence = reverse_sequence
)

stopifnot(
  located$best$start_forward == 23L,
  located$best$end_forward == 41L,
  located$best$start_reverse == 1096L,
  located$best$end_reverse == 1114L
)

primer_rows <- data.frame(
  pair_id = rep("CUSTOM_BOUNDS", 2L),
  pair_label = rep("Custom bounds regression", 2L),
  direction = c("forward", "reverse"),
  sequence = c(forward_sequence, reverse_sequence),
  stringsAsFactors = FALSE
)
geometry <- data.frame(
  forward_start = located$best$start_forward,
  forward_end = located$best$end_forward,
  reverse_start = located$best$start_reverse,
  reverse_end = located$best$end_reverse
)
taxonomy <- read_csv(
  "data/derived/gurten2026_centroid_taxonomy.csv",
  col_types = cols(accession = col_character()),
  show_col_types = FALSE
)

result <- score_primer_pair_expanded_reference(
  alignment_path = "data/external/gurten2026/ClusteredReferences.fasta",
  taxonomy = taxonomy,
  primer_rows = primer_rows,
  geometry = geometry,
  alignment_offset = 40L
)

stopifnot(
  nrow(result$scores) == nrow(taxonomy),
  !isTRUE(result$overall$reference_suitable),
  all(!result$scores$pair_scorable),
  grepl(
    "Reverse primer requests alignment columns 1136–1154",
    result$overall$reference_limitation_reason,
    fixed = TRUE
  ),
  grepl(
    "does not span both primer sites",
    result$overall$reference_limitation_reason,
    fixed = TRUE
  )
)

message("Expanded-reference out-of-bounds custom-pair regression passed.")
