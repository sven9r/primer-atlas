#!/usr/bin/env Rscript

root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(root, ".Rlib"), .libPaths()))
suppressPackageStartupMessages({library(jsonlite); library(digest); library(nanoparquet)})
source("R/functions.R")
source("R/data_layer.R")
source("R/coi_retention.R")
source("R/coi_scoring.R")

marker_id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(marker_id) || !nzchar(marker_id)) stop("Usage: validate_release.R MARKER_ID")
release_id <- Sys.getenv("RELEASE_ID", "")
if (!nzchar(release_id)) {
  candidates <- list.files(file.path("data", "releases", marker_id), pattern = "^[0-9].*", full.names = FALSE)
  release_id <- tail(sort(candidates), 1)
}
manifest_path <- file.path("data", "releases", marker_id, release_id, "manifest.json")
stopifnot(file.exists(manifest_path))
manifest <- read_json(manifest_path, simplifyVector = FALSE)
stopifnot(manifest$schema_version == "1.0.0", manifest$marker_id == marker_id)

pairs <- read.csv("data/catalog/primer_pairs.csv", stringsAsFactors = FALSE)
links <- read.csv("data/catalog/pair_oligos.csv", stringsAsFactors = FALSE)
oligos <- read.csv("data/catalog/oligos.csv", stringsAsFactors = FALSE)
claims <- read.csv("data/catalog/claims.csv", stringsAsFactors = FALSE)
sources <- read.csv("data/catalog/sources.csv", stringsAsFactors = FALSE)
facets <- read.csv("data/catalog/primer_pair_facets.csv", stringsAsFactors = FALSE)
active_pairs <- pairs[pairs$marker_id == marker_id, , drop = FALSE]
stopifnot(nrow(active_pairs) > 0L)
stopifnot(all(active_pairs$pair_id %in% links$pair_id))
active_oligos <- merge(links[links$pair_id %in% active_pairs$pair_id, ], oligos, by = "oligo_id")
stopifnot(all(active_oligos$direction %in% c("forward", "reverse")))
stopifnot(all(grepl("^[ACGTRYSWKMBDHVNI]+$", active_oligos$sequence)))
stopifnot(all(nzchar(active_oligos$source_key)))
source_keys <- unique(unlist(strsplit(active_oligos$source_key, "[|]")))
stopifnot(all(source_keys %in% sources$key), all(nzchar(sources$url[match(source_keys, sources$key)])))
if (nrow(claims[claims$marker_id == marker_id, , drop = FALSE])) {
  stopifnot(all(nzchar(claims$source_keys[claims$marker_id == marker_id])))
}
stopifnot(all(c("barcoding", "bulk_community", "edna", "diet") %in% facets$facet_value))

for (artifact in manifest$artifacts) {
  path <- artifact$local_path
  stopifnot(file.exists(path))
  stopifnot(identical(
    digest(file = path, algo = "sha256", serialize = FALSE), artifact$sha256
  ))
  if (grepl("[.]fasta$", path, ignore.case = TRUE)) {
    stopifnot(sum(startsWith(readLines(path, warn = FALSE), ">")) == artifact$rows)
  } else if (grepl("[.]parquet$", path, ignore.case = TRUE)) {
    stopifnot(nrow(nanoparquet::read_parquet(path)) == artifact$rows)
  } else if (grepl("[.]csv$", path, ignore.case = TRUE)) {
    stopifnot(nrow(read.csv(path, stringsAsFactors = FALSE)) == artifact$rows)
  } else if (grepl("[.]json$", path, ignore.case = TRUE)) {
    jsonlite::read_json(path, simplifyVector = FALSE)
  }
}

if (marker_id == "ITS_FUNGAL") {
  source("R/fungal_versions.R")
  stopifnot(length(manifest$fungal_snapshots) >= 2L)
  for (descriptor in manifest$fungal_snapshots) {
    snapshot <- read_json(descriptor$manifest$local_path, simplifyVector = FALSE)
    stopifnot(identical(snapshot$id, descriptor$id), snapshot$sequence_count >= 1L)
    required <- c("catalog", "raw_catalog", "sequences", "alignment", "taxonomy", "geometry", "landmarks", "pairs", "pair_sequences", "sites")
    stopifnot(all(required %in% names(snapshot$artifacts)))
    values <- lapply(required, function(name) fungal_snapshot_artifact(snapshot, name, dirname(descriptor$manifest$local_path)))
    names(values) <- required
    stopifnot(!anyDuplicated(names(values$sequences)), !anyDuplicated(values$taxonomy$accession))
    stopifnot(setequal(sub(" .*", "", names(values$sequences)), values$taxonomy$accession))
    stopifnot(length(unique(nchar(values$alignment))) == 1L)
    stopifnot(setequal(sub("^_R_", "", sub(" .*", "", names(values$alignment))), names(values$sequences)))
    stopifnot(snapshot$sequence_count == length(values$sequences))
    if (snapshot$version_number >= 2L) stopifnot(length(values$sequences) >= 2L, nrow(values$pairs) >= 10L)
  }
}

beeprime_expanded <- manifest$reference_panels$beeprime_expanded
if (!is.null(beeprime_expanded)) {
  stopifnot(
    identical(beeprime_expanded$pair_id, "BEEPRIME"),
    identical(beeprime_expanded$reference_mode, "expanded"),
    identical(beeprime_expanded$taxonomic_scope, "bees"),
    identical(beeprime_expanded$publication_baseline, "gurten2026_suppl5"),
    beeprime_expanded$added_sequences > 0L
  )
  expanded_scores <- nanoparquet::read_parquet(
    beeprime_expanded$sequence_scores$local_path
  )
  publication_scores <- read.csv(
    "data/derived/beeprime_hymenoptera_centroid_scores.csv",
    stringsAsFactors = FALSE
  )
  publication_accessions <- publication_scores$accession[
    publication_scores$is_bee %in% TRUE
  ]
  stopifnot(
    nrow(expanded_scores) == length(publication_accessions) + beeprime_expanded$added_sequences,
    all(expanded_scores$pair_id == "BEEPRIME"),
    all(expanded_scores$order == "Hymenoptera"),
    all(expanded_scores$is_bee %in% TRUE),
    !anyDuplicated(expanded_scores$accession),
    setequal(
      expanded_scores$accession[expanded_scores$reference_origin == "publication"],
      publication_accessions
    ),
    sum(expanded_scores$reference_origin == "added") == beeprime_expanded$added_sequences,
    nrow(nanoparquet::read_parquet(beeprime_expanded$geography$local_path)) == nrow(expanded_scores)
  )
  aligned_headers <- readLines(beeprime_expanded$alignment$local_path, warn = FALSE)
  aligned_sequences <- parse_fasta(beeprime_expanded$alignment$local_path)
  aligned_accessions <- sub("^_R_", "", sub(" .*", "", names(aligned_sequences)))
  publication_alignment_width <- unique(nchar(parse_fasta(
    "data/external/gurten2026/ClusteredReferences.fasta"
  )))
  stopifnot(
    !anyDuplicated(aligned_accessions),
    setequal(aligned_accessions, expanded_scores$accession),
    length(publication_alignment_width) == 1L,
    all(nchar(aligned_sequences) == publication_alignment_width)
  )
}

previous_release_id <- if (marker_id == "COI") Sys.getenv("COI_PREVIOUS_RELEASE_ID", "") else ""
previous_state <- tryCatch(atlas_read_manifest(marker_id, root,
  release_id = if (nzchar(previous_release_id)) previous_release_id else NULL), error = function(e) NULL)
if (!is.null(previous_state)) {
  previous <- previous_state$manifest
  for (name in intersect(names(manifest$artifacts), names(previous$artifacts))) {
    old_rows <- previous$artifacts[[name]]$rows
    new_rows <- manifest$artifacts[[name]]$rows
    stopifnot(new_rows >= 0.90 * old_rows)
  }
  for (metric in c("taxonomy_resolved_fraction", "geography_resolved_fraction")) {
    old <- previous$completeness[[metric]]
    new <- manifest$completeness[[metric]]
    if (!is.null(old) && !is.null(new) && is.finite(old) && is.finite(new)) stopifnot(new >= old - 0.05)
  }
  old_geometry <- tryCatch(atlas_read_artifact(previous$artifacts$pair_geometry, root), error = function(e) NULL)
  if (!is.null(old_geometry)) {
    new_geometry <- nanoparquet::read_parquet(manifest$artifacts$pair_geometry$local_path)
    coordinates <- c("pair_id", "forward_start", "forward_end", "reverse_start", "reverse_end")
    shared <- intersect(old_geometry$pair_id, new_geometry$pair_id)
    stopifnot(identical(
      old_geometry[match(shared, old_geometry$pair_id), coordinates, drop = FALSE],
      new_geometry[match(shared, new_geometry$pair_id), coordinates, drop = FALSE]
    ))
  }
}
if (!is.null(manifest$artifacts$taxonomy)) {
  taxonomy <- nanoparquet::read_parquet(manifest$artifacts$taxonomy$local_path)
  stopifnot(!anyDuplicated(taxonomy$accession))
}
if (marker_id == "COI") stopifnot(nrow(active_pairs) == 25L)
if (marker_id == "COI") {
  state <- manifest$coi_retained_state
  stopifnot(!is.null(state), state$schema_version == "1.0.0")
  checked <- coi_validate_retained_state(root)
  stopifnot(nrow(checked$records) == state$retained_records, state$retained_records >= 10000L,
    length(state$target_groups) == 22L, all(unlist(state$group_counts) >= 100L))
  scoring <- coi_scoring_inputs(root, mode = "full", minimum_per_group = 100L, minimum_total = 10000L)
  stopifnot(sum(scoring$manifest$n_templates) == nrow(checked$records))
  stopifnot(all(paste0("alignment_", unlist(state$target_groups)) %in% names(state$artifacts)),
    all(c("records", "ledger", "growth_qa", "panel_panel", "panel_panel_manifest", "panel_taxonomy", "panel_run_manifest", "scoring_manifest") %in% names(state$artifacts)))
  if (!is.null(previous_state$manifest$coi_retained_state$artifacts$records)) {
    previous_records <- as.data.frame(atlas_read_artifact(previous_state$manifest$coi_retained_state$artifacts$records, root))
    coi_assert_retention(previous_records, checked$records)
  }
  published_membership <- file.path(root, "data/provenance/coi_published_baseline.csv")
  if (file.exists(published_membership)) {
    prior <- read.csv(published_membership, stringsAsFactors = FALSE)
    stopifnot(all(paste(prior$target_group, prior$accession) %in% paste(checked$records$target_group, checked$records$accession)))
  }
  growth_qa <- read.csv(state$artifacts$growth_qa$local_path, stringsAsFactors = FALSE)
  stopifnot(
    nrow(growth_qa) == 22L,
    all(growth_qa$target_retained >= 100L),
    all(growth_qa$target_met),
    all(growth_qa$retained_after == as.integer(state$group_counts[growth_qa$target_group])),
    all(growth_qa$retained_after >= ceiling(growth_qa$previous_retained * 1.01)),
    all(growth_qa$verified_marker_records >= 100L),
    sum(growth_qa$verified_marker_records) >= 10000L
  )
  ledger <- checked$ledger
  panel_taxonomy <- read.csv(state$artifacts$panel_taxonomy$local_path, stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(panel_taxonomy$accession), setequal(panel_taxonomy$accession, ledger$accession))
  partition_keys <- paste0("full_order_", tolower(active_pairs$pair_id), "_primer_positions")
  stopifnot(all(partition_keys %in% names(manifest$artifacts)))
  scored_overview <- nanoparquet::read_parquet(manifest$artifacts$overview$local_path)
  expected_counts <- setNames(scoring$manifest$n_templates, scoring$manifest$order)
  stopifnot(nrow(scored_overview) == length(expected_counts) * nrow(active_pairs),
    !anyDuplicated(paste(scored_overview$order, scored_overview$pair_id)),
    setequal(unique(scored_overview$order), names(expected_counts)),
    setequal(unique(scored_overview$pair_id), active_pairs$pair_id),
    all(scored_overview$n_templates == expected_counts[scored_overview$order]),
    all(scored_overview$evidence_class == "full_order_reference_alignment"))
  scored_templates <- nanoparquet::read_parquet(manifest$artifacts$pair_templates$local_path)
  stopifnot(nrow(scored_templates) == nrow(checked$records) * nrow(active_pairs),
    !anyDuplicated(paste(scored_templates$order, scored_templates$pair_id, scored_templates$accession)),
    all(scored_templates$evidence_class == "full_order_reference_alignment"))
  retained_keys <- paste(checked$records$target_group, checked$records$accession)
  primer_rows <- read.csv("data/primers.csv", stringsAsFactors = FALSE)
  for (pair_id in active_pairs$pair_id) {
    pair_rows <- scored_templates[scored_templates$pair_id == pair_id, , drop = FALSE]
    stopifnot(setequal(paste(pair_rows$order, pair_rows$accession), retained_keys))
    pair_primers <- primer_rows[primer_rows$pair_id == pair_id, , drop = FALSE]
    stopifnot(nrow(pair_primers) == 2L, setequal(pair_primers$direction, c("forward", "reverse")))
    expected_positions <- nrow(checked$records) * sum(nchar(vapply(pair_primers$sequence, clean_sequence, character(1))))
    artifact <- manifest$artifacts[[paste0("full_order_", tolower(pair_id), "_primer_positions")]]
    stopifnot(artifact$rows == expected_positions)
    sequence_artifact <- manifest$artifacts[[paste0("full_order_", tolower(pair_id))]]
    stopifnot(!is.null(sequence_artifact), sequence_artifact$rows == nrow(checked$records))
  }
}
message("Release gates passed for ", marker_id, " ", manifest$release_version)
