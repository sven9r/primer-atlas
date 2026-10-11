#!/usr/bin/env Rscript
.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
source("R/functions.R")
source("R/data_layer.R")
source("R/coi_retention.R")
source("R/coi_scoring.R")
source("R/coi_reference.R")

fixture_root <- tempfile("coi-runtime-")
dir.create(file.path(fixture_root, "data/reference"), recursive = TRUE)
on.exit(unlink(fixture_root, recursive = TRUE), add = TRUE)
reference_path <- file.path(fixture_root, "data/reference/coi_reference_NC_001322.1.fasta")
stopifnot(file.copy("data/reference/coi_reference_NC_001322.1.fasta", reference_path))
reference <- parse_fasta(reference_path)
write_fasta <- function(value, path) writeLines(unlist(lapply(seq_along(value), function(i) c(paste0(">", names(value)[i]), value[i])), use.names = FALSE), path)
artifact <- function(path, rows = NULL) list(local_path = path, sha256 = atlas_sha256(path), rows = rows)
csv_artifact <- function(value, name) {
  path <- file.path(fixture_root, paste0(name, ".csv"))
  write.csv(value, path, row.names = FALSE, na = "")
  artifact(path, nrow(value))
}
groups <- c("Araneae", "Diptera")
write.csv(data.frame(order = groups), file.path(fixture_root, "data/orders_22.csv"), row.names = FALSE)
artifacts <- list()
raw_values <- list()
for (group in groups) {
  raw_dir <- file.path(fixture_root, "data/primerminer/raw", group)
  dir.create(raw_dir, recursive = TRUE)
  raw <- setNames(c(unname(reference), paste0("T", substring(reference, 2))),
                  paste0(group, c("_A.1", "_B.1"), " fixture"))
  raw_values[[group]] <- raw
  raw_path <- file.path(raw_dir, paste0(group, "_GB.fasta"))
  write_fasta(raw, raw_path)
  alignment_path <- file.path(fixture_root, paste0(group, "_aligned.fasta"))
  write_fasta(c(reference, raw), alignment_path)
  input_path <- paste0(alignment_path, ".inputs.json")
  inputs <- coi_alignment_inputs(raw_path, reference_path, "fixture")
  inputs$output_sha256 <- atlas_sha256(alignment_path)
  jsonlite::write_json(inputs, input_path, auto_unbox = TRUE, pretty = TRUE)
  artifacts[[paste0("raw_", group)]] <- artifact(raw_path, 2L)
  artifacts[[paste0("alignment_", group)]] <- artifact(alignment_path, 3L)
  artifacts[[paste0("alignment_inputs_", group)]] <- artifact(input_path)
}
records <- coi_retained_records(fixture_root, groups)
artifacts$records <- csv_artifact(records, "records")
taxonomy <- records[c("target_group", "accession")]
taxonomy$order <- taxonomy$target_group
taxonomy$organism <- "Example organism"
taxonomy$family <- c("Example family", "", "Example family", "")
taxonomy$genus <- ""
taxonomy$taxonomy_status <- "unresolved"
artifacts$panel_taxonomy <- csv_artifact(taxonomy, "taxonomy")
qa <- data.frame(target_group = groups, previous_retained = 0L, added_new = 2L,
  retained_after = 2L, target_met = TRUE, verified_marker_records = 2L)
artifacts$growth_qa <- csv_artifact(qa, "growth")
pairs <- c("MCO", "FOLMER")
overview <- expand.grid(order = groups, pair_id = pairs, stringsAsFactors = FALSE)
overview$n_templates <- 2L
scores <- do.call(rbind, lapply(groups, function(group) {
  do.call(rbind, lapply(pairs, function(pair) data.frame(order = group, pair_id = pair,
    pair_label = pair, template = names(raw_values[[group]]), accession = coi_accession(names(raw_values[[group]])),
    forward_primerminer_penalty = 0, reverse_primerminer_penalty = 0,
    forward_terminal_3_mismatches = 0, reverse_terminal_3_mismatches = 0,
    forward_scorable = TRUE, reverse_scorable = TRUE, pair_scorable = TRUE,
    pair_penalty = 0, perfect_match_pair = TRUE, stringsAsFactors = FALSE)))
}))
state <- list(manifest = list(marker_id = "COI", release_version = "runtime-fixture",
  artifacts = list(overview = csv_artifact(overview, "overview"), pair_templates = csv_artifact(scores, "scores")),
  coi_retained_state = list(retained_records = 4L, target_groups = groups,
    group_counts = list(Araneae = 2L, Diptera = 2L), artifacts = artifacts,
    policy = list(minimum_per_group = 1L, minimum_total = 4L, monthly_growth = .01))),
  source = "fixture", stale = FALSE)
loaded <- coi_reference_load(state, fixture_root, minimum_per_group = 1L, minimum_total = 4L)
stopifnot(loaded$mode == "full_release", loaded$n_references == 4L, loaded$n_verified_marker == 4L,
  sum(vapply(loaded$alignments, length, integer(1))) == 4L,
  nrow(loaded$pair_templates) == 8L, nrow(loaded$taxonomy) == 4L,
  is.na(loaded$taxonomy$family[2]), all(is.na(loaded$taxonomy$genus)),
  grepl("4 retained COI records", coi_reference_status(loaded)),
  grepl("4 additions", coi_reference_status(loaded)),
  all(is.na(coi_reference_geography(loaded)$country_or_territory)))

expect_invalid <- function(modified, pattern) {
  result <- tryCatch(.coi_reference_full_release(modified, fixture_root, 1L, 4L), error = identity)
  stopifnot(inherits(result, "error"), grepl(pattern, conditionMessage(result)))
}
bad_counts <- state
bad_counts$manifest$coi_retained_state$group_counts$Araneae <- 3L
expect_invalid(bad_counts, "group counts")
missing_group <- state
missing_group$manifest$coi_retained_state$target_groups <- "Araneae"
expect_invalid(missing_group, "catalog target groups")
bad_checksum <- state
bad_checksum$manifest$coi_retained_state$artifacts$alignment_Araneae$sha256 <- strrep("0", 64)
expect_invalid(bad_checksum, "checksum")
bad_scores <- state
bad_scores$manifest$artifacts$pair_templates <- csv_artifact(scores[-1, ], "short_scores")
expect_invalid(bad_scores, "different denominators")
bad_qa <- state
wrong_qa <- qa
wrong_qa$added_new[1] <- 1L
bad_qa$manifest$coi_retained_state$artifacts$growth_qa <- csv_artifact(wrong_qa, "wrong_growth")
expect_invalid(bad_qa, "growth proof")
bad_width <- state
truncated <- c(reference, raw_values[["Araneae"]])
truncated[2] <- substring(truncated[2], 2)
bad_path <- file.path(fixture_root, "wrong_width.fasta")
write_fasta(truncated, bad_path)
bad_width$manifest$coi_retained_state$artifacts$alignment_Araneae <- artifact(bad_path, 3L)
expect_invalid(bad_width, "output differs|width changed")

# Pair detail never reads a global artifact with millions of rows, and a
# corrupt partition cannot fall back to a legacy table from another cohort.
position_rows <- data.frame(order = "Araneae", pair_id = "MCO", template = "Araneae_A.1 fixture",
  unavailable_reason = "", stringsAsFactors = FALSE)
loaded$state$manifest$artifacts$full_order_mco_positions <- csv_artifact(position_rows, "mco_positions")
loaded$state$manifest$artifacts$primer_positions <- list(rows = 11000000L,
  local_path = file.path(fixture_root, "global-must-not-load.csv"), sha256 = strrep("0", 64))
positions <- coi_reference_pair_positions(loaded, "MCO")
stopifnot(nrow(positions) == 1L, positions$pair_id == "MCO",
  nrow(coi_reference_pair_positions(loaded, "FOLMER")) == 0L)
corrupt_partition <- loaded
corrupt_partition$cache <- new.env(parent = emptyenv())
corrupt_partition$state$manifest$artifacts$full_order_mco_positions$sha256 <- strrep("0", 64)
result <- tryCatch(coi_reference_pair_positions(corrupt_partition, "MCO"), error = identity)
stopifnot(inherits(result, "error"), grepl("checksum", conditionMessage(result)))
foreign <- position_rows
foreign$template <- "OUTSIDE_PANEL.1 fixture"
outside <- loaded
outside$cache <- new.env(parent = emptyenv())
outside$state$manifest$artifacts$full_order_mco_positions <- csv_artifact(foreign, "foreign_positions")
result <- tryCatch(coi_reference_pair_positions(outside, "MCO"), error = identity)
stopifnot(inherits(result, "error"), grepl("outside", conditionMessage(result)))
normalized <- coi_reference_pair_taxonomy(loaded, scores)
stopifnot(all(normalized$cluster_size == 1L), all(normalized$forward_penalty == 0),
  all(!normalized$terminal_3_mismatch_pair), all(normalized$reference_version == "runtime-fixture"))

# A complete local build may activate before publication, with an explicit
# label. Loose full-alignment files alone cannot replace the chosen release.
dir.create(file.path(fixture_root, "data/provenance"), recursive = TRUE)
dir.create(file.path(fixture_root, "data/derived/coi_reference_panel"), recursive = TRUE)
dir.create(file.path(fixture_root, "data/primerminer/full_aligned_reference"), recursive = TRUE)
write.csv(data.frame(order = groups), file.path(fixture_root, "data/orders_22.csv"), row.names = FALSE)
write.csv(records, file.path(fixture_root, "data/provenance/coi_reference_accession_ledger.csv"), row.names = FALSE)
stopifnot(file.copy("data/primers.csv", file.path(fixture_root, "data/primers.csv")),
  file.copy("data/derived/primer_bindings.csv", file.path(fixture_root, "data/derived/primer_bindings.csv")))
for (group in groups) {
  destination <- file.path(fixture_root, "data/primerminer/full_aligned_reference", paste0(group, "_COX1_full_reference_aligned.fasta"))
  stopifnot(file.copy(artifacts[[paste0("alignment_", group)]]$local_path, destination),
    file.copy(artifacts[[paste0("alignment_inputs_", group)]]$local_path, paste0(destination, ".inputs.json")))
}
complete <- coi_scoring_inputs(fixture_root, "full", groups, 1L, 4L)
scoring_manifest <- complete$manifest
scoring_manifest$primer_catalog_sha256 <- atlas_sha256(file.path(fixture_root, "data/primers.csv"))
scoring_manifest$binding_coordinates_sha256 <- atlas_sha256(file.path(fixture_root, "data/derived/primer_bindings.csv"))
scoring_manifest$scored_at_utc <- "2026-10-11 00:00:00 UTC"
write.csv(scoring_manifest, file.path(fixture_root, "data/provenance/order_scoring_manifest.csv"), row.names = FALSE)
local_overview <- overview
local_overview$evidence_class <- "full_order_reference_alignment"
write.csv(local_overview, file.path(fixture_root, "data/derived/order_pair_summary.csv"), row.names = FALSE)
write.csv(scores, file.path(fixture_root, "data/derived/pair_template_scores.csv"), row.names = FALSE)
write.csv(taxonomy, file.path(fixture_root, "data/derived/coi_reference_panel/taxonomy.csv"), row.names = FALSE, na = "")
write.csv(qa, file.path(fixture_root, "data/provenance/coi_reference_growth_qa.csv"), row.names = FALSE)
for (pair in pairs) {
  prefix <- file.path(fixture_root, "data/derived", paste0("full_order_", tolower(pair)))
  connection <- gzfile(paste0(prefix, "_sequence_scores.csv.gz"), "wt")
  write.csv(scores[scores$pair_id == pair, ], connection, row.names = FALSE)
  close(connection)
  connection <- gzfile(paste0(prefix, "_primer_positions.csv.gz"), "wt")
  positions <- position_rows
  positions$pair_id <- pair
  write.csv(positions, connection, row.names = FALSE)
  close(connection)
}
local <- .coi_reference_local_full(state, fixture_root, 1L, 4L)
stopifnot(local$mode == "full_local", local$n_references == 4L,
  grepl("publication pending", coi_reference_status(local)), nrow(coi_reference_pair_positions(local, "MCO")) == 1L)
scoring_manifest$alignment_sha256[1] <- strrep("0", 64)
write.csv(scoring_manifest, file.path(fixture_root, "data/provenance/order_scoring_manifest.csv"), row.names = FALSE)
result <- tryCatch(.coi_reference_local_full(state, fixture_root, 1L, 4L), error = identity)
stopifnot(inherits(result, "error"), grepl("stale", conditionMessage(result)))

# A corrupt frozen release changes the entire selected cohort to a named
# fallback. Its good tables cannot be combined with another panel's rows.
dir.create(file.path(fixture_root, "data/releases/COI"), recursive = TRUE)
dir.create(file.path(fixture_root, "data/primerminer/aligned_reference"), recursive = TRUE)
for (group in groups) {
  stopifnot(file.copy(artifacts[[paste0("alignment_", group)]]$local_path,
    file.path(fixture_root, "data/primerminer/aligned_reference", paste0(group, "_COX1_reference_aligned.fasta"))))
}
fallback <- state$manifest
fallback$coi_retained_state <- NULL
fallback$release_version <- "legacy-fixture"
jsonlite::write_json(fallback, file.path(fixture_root, "data/releases/COI/latest.json"), auto_unbox = TRUE, pretty = TRUE)
fallback_loaded <- coi_reference_load(bad_checksum, fixture_root, 1L, 4L)
stopifnot(fallback_loaded$mode == "legacy_preview", fallback_loaded$version == "legacy-fixture",
  fallback_loaded$state$manifest$release_version != "runtime-fixture", grepl("checksum", fallback_loaded$issue))

# The checked-in preview remains a distinct dataset even when a real new full
# local build is installed by another worker during this test run.
preview <- .coi_reference_preview(getwd())
pinned <- nanoparquet::read_parquet("data/pinned/COI/2026-08-14/overview.parquet")
stopifnot(sum(preview$counts) == 1431L, sum(.coi_reference_counts(pinned)) == 1431L)
cat("COI runtime checks passed: coherent full cohort, frozen identities/counts, growth proof, no stale mixtures, pair-lazy positions and legacy labeling.\n")
