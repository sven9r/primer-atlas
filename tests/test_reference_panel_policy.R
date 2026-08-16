#!/usr/bin/env Rscript

policy <- read.csv("data/catalog/reference_panel_policy.csv", stringsAsFactors = FALSE)
orders <- read.csv("data/orders_22.csv", stringsAsFactors = FALSE)
unite <- read.csv("data/catalog/unite_primers.csv", stringsAsFactors = FALSE)
oligos <- read.csv("data/catalog/oligos.csv", stringsAsFactors = FALSE)

general <- policy[policy$target_group == "ALL_ARTHROPOD_GROUPS", , drop = FALSE]
bee <- policy[policy$exception_pair_id == "BEEPRIME", , drop = FALSE]
stopifnot(
  nrow(general) == 1L,
  general$baseline_sequences == 1000L,
  abs(general$monthly_growth_fraction - 0.01) < 1e-12,
  general$update_mode == "append_unseen_accessions",
  general$scoring_input == "full_reference_coordinate_alignment",
  nrow(bee) == 1L,
  bee$scoring_input == "author_supplied_99.5_percent_centroids",
  sum(orders$order != "Bivalvia") == 22L
)

update_source <- paste(readLines("scripts/update_order_reference_panels.R", warn = FALSE), collapse = "\n")
alignment_source <- paste(readLines("scripts/build_22_group_alignments.R", warn = FALSE), collapse = "\n")
release_source <- paste(readLines("scripts/build_release.R", warn = FALSE), collapse = "\n")
stopifnot(
  grepl("ceiling\\(previous_n \\* \\(1 \\+ growth\\)\\)", update_source),
  grepl("setdiff\\(as.character\\(search\\$ids\\), existing\\$entrez_uid\\)", update_source),
  grepl('EVALUATION_MODE", "full', alignment_source, fixed = TRUE),
  grepl("full_order_", release_source, fixed = TRUE),
  grepl("exact_beeprime", release_source, fixed = TRUE),
  !grepl('pattern = "\\^claimed_\\.\\*', release_source)
)

its_oligos <- oligos[oligos$marker_id == "ITS_FUNGAL", , drop = FALSE]
stopifnot(
  nrow(unite) >= 120L,
  all(c("primary_reference_key", "primary_reference_url", "primary_reference_status") %in% names(unite)),
  all(nzchar(unite$primary_reference_key)),
  all(nzchar(unite$primary_reference_url)),
  all(its_oligos$source_key != "unite_primers"),
  all(grepl("^its_primary_", its_oligos$source_key))
)

message("Reference-panel growth, evidence separation, and primary ITS citation checks passed.")
