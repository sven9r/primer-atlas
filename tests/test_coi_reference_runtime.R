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
position_rows <- data.frame(order = scores$order[scores$pair_id == "MCO"], pair_id = "MCO", template = scores$template[scores$pair_id == "MCO"],
  unavailable_reason = "", stringsAsFactors = FALSE)
loaded$state$manifest$artifacts$full_order_mco_positions <- csv_artifact(position_rows, "mco_positions")
loaded$state$manifest$artifacts$primer_positions <- list(rows = 11000000L,
  local_path = file.path(fixture_root, "global-must-not-load.csv"), sha256 = strrep("0", 64))
positions <- coi_reference_pair_positions(loaded, "MCO")
stopifnot(nrow(positions) == 4L, all(positions$pair_id == "MCO"),
  nrow(coi_reference_pair_positions(loaded, "FOLMER")) == 0L)
corrupt_partition <- loaded
corrupt_partition$cache <- new.env(parent = emptyenv())
corrupt_partition$state$manifest$artifacts$full_order_mco_positions$sha256 <- strrep("0", 64)
result <- tryCatch(coi_reference_pair_positions(corrupt_partition, "MCO"), error = identity)
stopifnot(inherits(result, "error"), grepl("checksum", conditionMessage(result)))
foreign <- position_rows
foreign$template[1L] <- "OUTSIDE_PANEL.1 fixture"
outside <- loaded
outside$cache <- new.env(parent = emptyenv())
outside$state$manifest$artifacts$full_order_mco_positions <- csv_artifact(foreign, "foreign_positions")
result <- tryCatch(coi_reference_pair_positions(outside, "MCO"), error = identity)
stopifnot(inherits(result, "error"), grepl("outside", conditionMessage(result)))
incomplete <- loaded
incomplete$cache <- new.env(parent = emptyenv())
partial_positions <- rbind(position_rows[-1L, ], position_rows[2L, ])
partial_positions$alignment_position <- c(1L, 1L, 1L, 2L)
incomplete$state$manifest$artifacts$full_order_mco_positions <- csv_artifact(partial_positions, "incomplete_positions")
result <- tryCatch(coi_reference_pair_positions(incomplete, "MCO"), error = identity)
stopifnot(inherits(result, "error"), grepl("every template", conditionMessage(result)))
normalized <- coi_reference_pair_taxonomy(loaded, scores)
stopifnot(all(normalized$cluster_size == 1L), all(normalized$forward_penalty == 0),
  all(!normalized$terminal_3_mismatch_pair), all(normalized$reference_version == "runtime-fixture"))

# Full partitions must match the requested pair and the complete selected
# cohort; BeePrime's named publication scope remains independent.
pair_rows <- scores[scores$pair_id == "MCO", ]
checked <- coi_reference_check_pair_scores(loaded, pair_rows, "MCO")
stopifnot(nrow(checked) == 4L)
wrong_pair <- pair_rows
wrong_pair$pair_id[1L] <- "FOLMER"
foreign_pair <- pair_rows
foreign_pair$order[1L] <- "Foreign"
wrong_accession <- pair_rows
wrong_accession$accession[1L] <- "Foreign.1"
for (bad in list(wrong_pair, foreign_pair, wrong_accession, pair_rows[-1L, ])) {
  result <- tryCatch(coi_reference_check_pair_scores(loaded, bad, "MCO"), error = identity)
  stopifnot(inherits(result, "error"))
}
publication_rows <- pair_rows[1L, ]
publication_rows$pair_id <- "BEEPRIME"
publication_rows$order <- "Independent publication group"
stopifnot(identical(coi_reference_check_pair_scores(loaded, publication_rows, "BEEPRIME"), publication_rows))

normalized_pair <- coi_reference_pair_taxonomy(loaded, checked)
families <- coi_reference_taxon_summary(normalized_pair, "Araneae", "family")
unknown_family <- families[families$taxon == "family unresolved", ]
unknown_rows <- coi_reference_taxon_filter(normalized_pair[normalized_pair$order == "Araneae", ], "family", unknown_family$taxon)
stopifnot(sum(families$n_centroids) == 2L, unknown_family$n_centroids == 1L,
          nrow(unknown_rows) == 1L, is.na(unknown_rows$family),
          unknown_rows$organism_label == "Example organism", unknown_rows$scientific_name_ncbi == "Example organism")
exact <- coi_reference_selected_alignments(loaded, unknown_rows$accession)
stopifnot(length(exact) == 1L, nchar(exact) == 1536L,
          coi_accession(names(exact)) == unknown_rows$accession)
result <- tryCatch(coi_reference_selected_alignments(loaded, "Foreign.1"), error = identity)
stopifnot(inherits(result, "error"))

# Exercise the actual main-app reactive branches with a small frozen cohort,
# including selection of an unresolved family and its exact organisms.
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(dplyr))
app_ast <- parse("app.R")
assignment_name <- function(value) {
  if (is.call(value) && identical(value[[1L]], as.name("<-")) && is.symbol(value[[2L]])) as.character(value[[2L]]) else ""
}
app_assignment <- function(name) app_ast[[which(vapply(app_ast, assignment_name, character(1)) == name)]]
runtime_environment <- new.env(parent = globalenv())
runtime_environment$coi_reference_panel <- loaded
runtime_environment$coi_release_state <- loaded$state
runtime_environment$coi_release_state$manifest$artifacts$full_order_mco <- csv_artifact(pair_rows, "mco_sequence_partition")
runtime_environment$coi_reference_panel$state <- runtime_environment$coi_release_state
runtime_environment$coi_reference_current <- function() runtime_environment$coi_reference_panel
eval(app_assignment("pair_sequence_evidence_spec"), runtime_environment)
eval(app_assignment("read_pair_sequence_evidence"), runtime_environment)
server_statements <- as.list(app_assignment("server")[[3L]][[3L]])[-1L]
callback_names <- c("custom_detail_active", "full_general_detail", "full_general_raw_scores", "claimed_pair_overall",
                    "targeted_zbj_active", "claimed_selected_scope", "claimed_taxon_data", "claimed_selected_orders",
                    "claimed_taxon_selection", "claimed_taxon_selection_observer", "claimed_selected_taxon",
                    "claimed_raw_scores", "claimed_selected_sequence_data")
callbacks <- server_statements[vapply(server_statements, assignment_name, character(1)) %in% callback_names]
fixture_server <- function(input, output, session) {
  beeprime_reference_mode <- function() "publication"
  for (statement in callbacks) eval(statement, environment())
}
environment(fixture_server) <- runtime_environment
shiny::testServer(fixture_server, {
  session$setInputs(detail_pair = "MCO", detail_order = "Araneae", claimed_taxon_rank = "family", claimed_taxon_sort = "taxon")
  stopifnot(full_general_detail(), !targeted_zbj_active(), claimed_pair_overall()$n_centroids == 4L,
            claimed_selected_scope()$n_centroids == 2L, sum(claimed_taxon_data()$n_centroids) == 2L)
  unknown_index <- which(claimed_taxon_data()$taxon == "family unresolved")
  session$setInputs(claimed_taxon_table_rows_selected = unknown_index)
  selected <- claimed_selected_sequence_data()
  stopifnot(nrow(selected) == 1L, selected$accession == unknown_rows$accession,
            selected$organism_label == "Example organism", is.na(selected$family))
})

# Human-readable display projections preserve the exact source row indices;
# filtered downloads keep all scientific fields and the frozen gapped bases.
lineage_names <- c("lineage_geography", "lineage_scope", "lineage_summary_table", "lineage_raw_scores",
  "lineage_family_selection", "lineage_family_selection_observer", "lineage_selected_family", "lineage_selected_family_order", "lineage_sequence_data", "lineage_summary_display", "lineage_sequence_display",
  "lineage_selected_alignments", "write_lineage_sequence_csv", "write_lineage_alignment_fasta")
lineage_callbacks <- server_statements[vapply(server_statements, assignment_name, character(1)) %in% lineage_names]
lineage_server <- function(input, output, session) {
  beeprime_reference_mode <- function() "publication"
  claimed_primer_families <- claimed_primer_subfamilies <- claimed_primer_genera <- tibble()
  for (statement in lineage_callbacks) eval(statement, environment())
}
environment(lineage_server) <- runtime_environment
shiny::testServer(lineage_server, {
  session$setInputs(lineage_marker = "COI", lineage_pair = "MCO", lineage_target = "Araneae",
                    lineage_country = "All", lineage_locality = "All")
  source_summary <- lineage_summary_table(claimed_primer_families, "family")
  display_summary <- lineage_summary_display(source_summary, "family")
  stopifnot(identical(display_summary$Family, source_summary$family),
            identical(display_summary$References, source_summary$n_sequences),
            identical(display_summary$`Scorable %`, round(source_summary$pair_scorable_fraction * 100, 1)),
            !any(c("reference_version", "pair_id", "pair_label") %in% names(display_summary)))
  index <- which(source_summary$family == "family unresolved")
  session$setInputs(lineage_families_rows_selected = index)
  selected <- lineage_sequence_data()
  display_sequences <- lineage_sequence_display(selected)
  stopifnot(nrow(selected) == 1L, identical(selected$accession, unknown_rows$accession),
            identical(display_sequences$Accession, selected$accession), display_sequences$Family == "Unresolved",
            sum(names(display_sequences) == "Organism") == 1L, !"reference_version" %in% names(display_sequences))
  csv_path <- tempfile(fileext = ".csv")
  fasta_path <- tempfile(fileext = ".fasta")
  write_lineage_sequence_csv(csv_path)
  write_lineage_alignment_fasta(fasta_path)
  exported <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
  exported_fasta <- parse_fasta(fasta_path)
  stopifnot(identical(names(exported), names(selected)), identical(exported$accession, selected$accession),
            isTRUE(all.equal(exported$pair_penalty, selected$pair_penalty)),
            identical(exported_fasta, lineage_selected_alignments()), nchar(exported_fasta) == 1536L)
  unlink(c(csv_path, fasta_path))
})

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
  grepl("publication pending", coi_reference_status(local)), nrow(coi_reference_pair_positions(local, "MCO")) == 4L)

# A local cohort must meet the same growth and marker-proof contract as an
# immutable hosted release; plausible score counts cannot rescue stale QA.
growth_path <- file.path(fixture_root, "data/provenance/coi_reference_growth_qa.csv")
python_qa <- qa
python_qa$target_met <- "True"
write.csv(python_qa, growth_path, row.names = FALSE)
stopifnot(.coi_reference_local_full(state, fixture_root, 1L, 4L)$mode == "full_local")
invalid_growth <- list(
  missing_group = qa[-1L, ],
  foreign_group = transform(qa, target_group = c("Foreign", groups[2L])),
  stale_count = transform(qa, retained_after = retained_after + 1L),
  failed_target = transform(qa, target_met = FALSE),
  python_false = transform(qa, target_met = "False"),
  invalid_flag = transform(qa, target_met = "unknown"),
  missed_growth = transform(qa, previous_retained = retained_after, added_new = 0L),
  unverified = transform(qa, verified_marker_records = 0L),
  inflated_verified = transform(qa, verified_marker_records = retained_after + 1L),
  missing_value = transform(qa, verified_marker_records = NA_integer_)
)
for (bad in invalid_growth) {
  write.csv(bad, growth_path, row.names = FALSE)
  result <- tryCatch(.coi_reference_local_full(state, fixture_root, 1L, 4L), error = identity)
  stopifnot(inherits(result, "error"), grepl("growth proof", conditionMessage(result)))
}
write.csv(qa, growth_path, row.names = FALSE)
stopifnot(.coi_reference_local_full(state, fixture_root, 1L, 4L)$mode == "full_local")
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
