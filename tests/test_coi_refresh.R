#!/usr/bin/env Rscript
.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
source("R/functions.R")
source("R/data_layer.R")
source("R/coi_retention.R")
source("R/coi_scoring.R")
source("R/coi_reference.R")
source("R/coi_refresh.R")
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(dplyr))

root <- tempfile("coi-refresh-fixture-")
dir.create(file.path(root, "data/reference"), recursive = TRUE)
reference_file <- file.path(root, "data/reference/coi_reference_NC_001322.1.fasta")
stopifnot(file.copy("data/reference/coi_reference_NC_001322.1.fasta", reference_file))
reference <- parse_fasta(reference_file)
groups <- c("Araneae", "Diptera")
write.csv(data.frame(order = groups), file.path(root, "data/orders_22.csv"), row.names = FALSE)
write_fasta <- function(value, path) writeLines(unlist(lapply(seq_along(value), function(i) c(paste0(">", names(value)[i]), value[i]))), path)
descriptor <- function(path, rows = NULL) list(local_path = path, sha256 = atlas_sha256(path), rows = rows)

make_release <- function(version, per_group, previous = 0L) {
  directory <- file.path(root, version)
  dir.create(directory)
  csv <- function(value, name) {
    path <- file.path(directory, paste0(name, ".csv"))
    write.csv(value, path, row.names = FALSE, na = "")
    descriptor(path, nrow(value))
  }
  retained <- list()
  records <- scores <- positions <- list()
  for (group in groups) {
    bases <- vapply(seq_len(per_group), function(i) {
      if (i == 1L) unname(reference) else paste0(paste(rep("T", i), collapse = ""), substring(reference, i + 1L))
    }, character(1))
    names(bases) <- paste0(group, "_", LETTERS[seq_len(per_group)], ".1 fixture")
    raw_path <- file.path(directory, paste0(group, "_raw.fasta"))
    aligned_path <- file.path(directory, paste0(group, "_aligned.fasta"))
    write_fasta(bases, raw_path)
    write_fasta(c(reference, bases), aligned_path)
    input_path <- paste0(aligned_path, ".inputs.json")
    jsonlite::write_json(list(evaluation_mode = "full", reference_sha256 = atlas_sha256(reference_file),
      raw_sha256 = atlas_sha256(raw_path), output_sha256 = atlas_sha256(aligned_path)), input_path, auto_unbox = TRUE)
    retained[[paste0("raw_", group)]] <- descriptor(raw_path, per_group)
    retained[[paste0("alignment_", group)]] <- descriptor(aligned_path, per_group + 1L)
    retained[[paste0("alignment_inputs_", group)]] <- descriptor(input_path)
    accessions <- coi_accession(names(bases))
    records[[group]] <- data.frame(target_group = group, accession = accessions,
      sequence_sha256 = vapply(unname(bases), digest::digest, character(1), algo = "sha256", serialize = FALSE),
      header_sha256 = vapply(names(bases), digest::digest, character(1), algo = "sha256", serialize = FALSE))
    scores[[group]] <- data.frame(order = group, pair_id = "MCO", pair_label = "MCO",
      template = names(bases), accession = accessions, forward_primerminer_penalty = 0,
      reverse_primerminer_penalty = 0, forward_terminal_3_mismatches = 0,
      reverse_terminal_3_mismatches = 0, forward_scorable = TRUE, reverse_scorable = TRUE,
      pair_scorable = TRUE, pair_penalty = 0, perfect_match_pair = TRUE)
    positions[[group]] <- data.frame(order = group, pair_id = "MCO", template = names(bases),
      alignment_position = 1L, template_base = substring(bases, 1, 1), position_penalty = 0)
  }
  records <- bind_rows(records)
  scores <- bind_rows(scores)
  positions <- bind_rows(positions)
  taxonomy <- records[c("target_group", "accession")]
  taxonomy$order <- taxonomy$target_group
  taxonomy$organism <- paste("Organism", taxonomy$accession)
  taxonomy$family <- rep(c("Example family", "", rep("AAA new family", max(0L, per_group - 2L))), length(groups))
  retained$records <- csv(records, "records")
  retained$panel_taxonomy <- csv(taxonomy, "taxonomy")
  retained$growth_qa <- csv(data.frame(target_group = groups, previous_retained = previous,
    added_new = per_group - previous, retained_after = per_group, target_met = "True",
    verified_marker_records = per_group), "growth")
  overview <- data.frame(order = groups, pair_id = "MCO", n_templates = per_group)
  publication <- scores[1L, ]
  publication$pair_id <- "BEEPRIME"
  publication$order <- "Independent publication group"
  list(manifest = list(marker_id = "COI", release_version = version,
    artifacts = list(overview = csv(overview, "overview"), pair_templates = csv(scores, "pair_templates"),
      full_order_mco = csv(scores, "sequence_partition"), full_order_mco_primer_positions = csv(positions, "positions"),
      exact_beeprime = csv(publication, "publication")),
    coi_retained_state = list(schema_version = "1.0.0", retained_records = nrow(records),
      target_groups = groups, group_counts = as.list(setNames(rep(per_group, 2L), groups)),
      artifacts = retained, policy = list(minimum_per_group = 1L, minimum_total = 4L, monthly_growth = .01))),
    stale = FALSE, source = "mock latest pointer")
}
first_state <- make_release("refresh-v1", 2L)
second_state <- make_release("refresh-v2", 3L, 2L)
third_state <- make_release("refresh-v3", 4L, 3L)
loader <- function(state) .coi_reference_full_release(state, root, 1L, 4L)
first <- loader(first_state)
latest <- first_state
reader <- function() latest
coordinator <- coi_reference_refresh_coordinator(first, root, reader, loader, now = function() "fixture UTC")
other_session <- coi_reference_refresh_coordinator(first, root, reader, loader)
stopifnot(!identical(coordinator$current()$cache, first$cache),
  !identical(coordinator$current()$cache, other_session$current()$cache),
  !coordinator$check()$changed, coordinator$status()$phase == "current",
  grepl("bundled reference", coordinator$status()$message), !grepl("latest verified hosted", coordinator$status()$message))
hosted_pointer <- first_state
hosted_pointer$source <- "https://example.org/releases/COI/latest.json"
hosted_coordinator <- coi_reference_refresh_coordinator(first, root, function() hosted_pointer, loader)
stopifnot(!hosted_coordinator$check()$changed, grepl("latest verified hosted", hosted_coordinator$status()$message))
latest <- second_state
updated <- coordinator$check()
stopifnot(updated$changed, updated$reference$version == "refresh-v2", updated$reference$n_references == 6L,
  coordinator$current()$version == "refresh-v2", other_session$current()$version == "refresh-v1", first$n_references == 4L)
rows <- coi_reference_pair_sequences(coordinator$current(), "MCO")
trace <- coi_reference_pair_positions(coordinator$current(), "MCO")
stopifnot(nrow(rows) == 6L, all(rows$reference_version == "refresh-v2"), nrow(trace) == 6L,
  setequal(trace$template, rows$template), length(coi_reference_selected_alignments(coordinator$current(), rows$accession)) == 6L)
stopifnot(identical(lapply(coi_reference_pair_sequences(first, "BEEPRIME"), as.character),
  lapply(coi_reference_pair_sequences(coordinator$current(), "BEEPRIME"), as.character)))

assert_retained <- function(state, pattern) {
  latest <<- state
  before <- coordinator$current()
  result <- coordinator$check()
  stopifnot(!result$changed, result$status$phase == "retained", grepl(pattern, result$status$message),
    identical(coordinator$current(), before), nrow(coi_reference_pair_sequences(coordinator$current(), "MCO")) == 6L)
}
corrupt <- third_state
corrupt$manifest$coi_retained_state$artifacts$alignment_Araneae$sha256 <- strrep("0", 64L)
assert_retained(corrupt, "checksum")
corrupt_trace <- third_state
corrupt_trace$manifest$artifacts$full_order_mco_primer_positions$sha256 <- strrep("0", 64L)
assert_retained(corrupt_trace, "checksum")
corrupt_sequence <- third_state
corrupt_sequence$manifest$artifacts$full_order_mco$sha256 <- strrep("0", 64L)
assert_retained(corrupt_sequence, "checksum")
dual_trace <- third_state
dual_trace$manifest$artifacts$full_order_mco_positions <- dual_trace$manifest$artifacts$full_order_mco_primer_positions
dual_trace$manifest$artifacts$full_order_mco_positions$sha256 <- strrep("0", 64L)
assert_retained(dual_trace, "checksum")
dual_read <- tryCatch(coi_reference_pair_positions(loader(dual_trace), "MCO"), error = identity)
stopifnot(inherits(dual_read, "error"), grepl("checksum", conditionMessage(dual_read)))
missing <- third_state
missing$manifest$artifacts$full_order_mco_primer_positions <- NULL
assert_retained(missing, "partition")
mutation <- second_state
mutation$manifest$artifacts$full_order_mco$sha256 <- strrep("0", 64L)
assert_retained(mutation, "immutable")
changed_record <- third_state
path <- changed_record$manifest$coi_retained_state$artifacts$records$local_path
records <- read.csv(path, stringsAsFactors = FALSE)
records$sequence_sha256[1L] <- strrep("a", 64L)
path <- file.path(root, "changed-records.csv")
write.csv(records, path, row.names = FALSE)
changed_record$manifest$coi_retained_state$artifacts$records <- descriptor(path, nrow(records))
assert_retained(changed_record, "previously retained sequence_sha256")
unavailable <- second_state
unavailable$stale <- TRUE
assert_retained(unavailable, "pointer")

# Run the real app's selection/read/download callbacks around the real refresh
# module. Added family AAA sorts ahead of the selected unresolved family, so
# retaining a row index would silently download the wrong accession after v2.
ast <- parse("app.R")
assignment_name <- function(x) if (is.call(x) && identical(x[[1L]], as.name("<-")) && is.symbol(x[[2L]])) as.character(x[[2L]]) else ""
assignment <- function(name) ast[[which(vapply(ast, assignment_name, character(1)) == name)]]
statements <- as.list(assignment("server")[[3L]][[3L]])[-1L]
names_needed <- c("current_order_scores", "current_pair_template_scores", "read_order_position_scores", "read_pair_sequence_evidence",
  "custom_location", "custom_pair_collection", "custom_pair_counter", "selected_custom_id", "its_catalog_pending",
  "archived_coi_custom_pairs", "custom_cohort_version", "coi_custom_refresh_observer",
  "order_lens_scores", "order_lens_template_scores",
  "custom_detail_active", "full_general_detail", "full_general_raw_scores", "claimed_pair_overall", "targeted_zbj_active",
  "claimed_selected_scope", "claimed_taxon_data", "claimed_selected_orders", "claimed_taxon_selection",
  "claimed_taxon_selection_observer", "claimed_selected_taxon", "claimed_raw_scores", "claimed_selected_sequence_data",
  "lineage_geography", "lineage_scope", "lineage_summary_table", "lineage_raw_scores", "lineage_family_selection",
  "lineage_family_selection_observer", "lineage_selected_family", "lineage_selected_family_order", "lineage_sequence_data", "lineage_selected_alignments",
  "write_lineage_sequence_csv", "write_lineage_alignment_fasta")
callbacks <- statements[vapply(statements, assignment_name, character(1)) %in% names_needed]
runtime <- new.env(parent = globalenv())
runtime$coi_reference_panel <- first
for (name in c("pair_sequence_evidence_spec", "read_pair_sequence_evidence", "read_order_position_scores",
               "coi_read_pair_sequence_evidence", "coi_read_order_position_scores")) eval(assignment(name), runtime)
app_fixture <- function(input, output, session) {
  refresh <- coiReferenceRefreshServer("coi", first, root, interval_ms = 600000L, reader = reader, loader = loader)
  coi_reference_current <- refresh$reference
  beeprime_reference_mode <- function() "publication"
  claimed_primer_families <- claimed_primer_subfamilies <- claimed_primer_genera <- tibble()
  for (statement in callbacks) eval(statement, environment())
}
environment(app_fixture) <- runtime
latest <- first_state
shiny::testServer(app_fixture, {
  session$setInputs(detail_pair = "MCO", detail_order = "Araneae", claimed_taxon_rank = "family", claimed_taxon_sort = "taxon",
    lineage_marker = "COI", lineage_pair = "MCO", lineage_target = "Araneae", lineage_country = "All", lineage_locality = "All")
  first_family_index <- which(lineage_summary_table(tibble(), "family")$family == "family unresolved")
  first_taxon_index <- which(claimed_taxon_data()$taxon == "family unresolved")
  session$setInputs(lineage_families_rows_selected = first_family_index, claimed_taxon_table_rows_selected = first_taxon_index)
  selected_accession <- lineage_sequence_data()$accession
  stopifnot(length(selected_accession) == 1L, nrow(full_general_raw_scores()) == 4L,
    sum(current_order_scores()$n_templates) == 4L, all(claimed_selected_sequence_data()$accession == selected_accession))
  session$setInputs(lineage_target = "All", lineage_families_rows_selected = integer())
  all_families <- lineage_summary_table(tibble(), "family")
  stopifnot(sum(all_families$family == "family unresolved") == 2L)
  first_family_index <- which(all_families$order == "Araneae" & all_families$family == "family unresolved")
  session$setInputs(lineage_families_rows_selected = first_family_index)
  stopifnot(lineage_selected_family_order() == "Araneae", nrow(lineage_sequence_data()) == 1L,
    identical(lineage_sequence_data()$accession, selected_accession))
  custom_pair_collection(list(CUSTOM_1 = list(marker_id = "COI", coi_reference_version = "refresh-v1",
    order_scores = data.frame(pair_id = "CUSTOM_1", order = groups, n_templates = 2L),
    pair_template_scores = tibble(), primer_position_scores = tibble())))
  custom_location(list(marker_id = "COI"))
  selected_custom_id("CUSTOM_1")
  stopifnot("CUSTOM_1" %in% order_lens_scores()$pair_id)
  latest <<- second_state
  session$setInputs(`coi-check_latest` = 1L)
  stopifnot(coi_reference_current()$version == "refresh-v2", nrow(full_general_raw_scores()) == 6L,
    sum(current_order_scores()$n_templates) == 6L, nrow(current_pair_template_scores()) == 6L,
    nrow(read_order_position_scores("MCO")) == 6L, claimed_pair_overall()$n_centroids == 6L,
    claimed_selected_scope()$n_centroids == 3L,
    first_taxon_index != which(claimed_taxon_data()$taxon == "family unresolved"),
    first_family_index != which(lineage_summary_table(tibble(), "family")$order == "Araneae" &
      lineage_summary_table(tibble(), "family")$family == "family unresolved"),
    lineage_selected_family() == "family unresolved", identical(lineage_sequence_data()$accession, selected_accession),
    lineage_selected_family_order() == "Araneae",
    identical(claimed_selected_sequence_data()$accession, selected_accession),
    all(lineage_sequence_data()$reference_version == "refresh-v2"))
  stopifnot(!length(custom_pair_collection()), length(archived_coi_custom_pairs()) == 1L,
    is.null(custom_location()), is.null(selected_custom_id()), !"CUSTOM_1" %in% order_lens_scores()$pair_id,
    nrow(order_lens_template_scores()) == 6L)
  csv_path <- tempfile(fileext = ".csv")
  fasta_path <- tempfile(fileext = ".fasta")
  write_lineage_sequence_csv(csv_path)
  write_lineage_alignment_fasta(fasta_path)
  stopifnot(identical(read.csv(csv_path)$accession, selected_accession),
    identical(coi_accession(names(parse_fasta(fasta_path))), selected_accession))
  unlink(c(csv_path, fasta_path))
  latest <<- corrupt
  session$setInputs(`coi-check_latest` = 2L)
  stopifnot(coi_reference_current()$version == "refresh-v2", refresh$status()$phase == "retained",
    identical(lineage_sequence_data()$accession, selected_accession), nrow(read_order_position_scores("MCO")) == 6L)
  # A periodic check also detects a new pointer in the same live session.
  latest <<- third_state
  session$elapse(600001L)
  session$flushReact()
  stopifnot(coi_reference_current()$version == "refresh-v3", nrow(full_general_raw_scores()) == 8L,
    identical(lineage_sequence_data()$accession, selected_accession))
})
latest <- second_state
shiny::testServer(coiReferenceRefreshServer, args = list(initial = first, root = root, reader = reader, loader = loader), {
  session$flushReact()
  stopifnot(reference()$version == "refresh-v2", reference()$n_references == 6L)
})
unlink(root, recursive = TRUE)
message("COI refresh checks passed: session isolation, atomic verified updates, invalid-release retention, stable taxon identities, refreshed fit/trace/download readers, periodic and new-session checks.")
