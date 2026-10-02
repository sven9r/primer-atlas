#!/usr/bin/env Rscript

root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(root, ".Rlib"), .libPaths()))
suppressPackageStartupMessages({library(jsonlite); library(digest); library(nanoparquet)})
source("R/functions.R")
source("R/data_layer.R")

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
  } else {
    stopifnot(nrow(nanoparquet::read_parquet(path)) == artifact$rows)
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

previous_state <- tryCatch(atlas_read_manifest(marker_id, root), error = function(e) NULL)
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
if (marker_id == "COI" && file.exists("data/provenance/coi_reference_growth_qa.csv")) {
  growth_qa <- read.csv("data/provenance/coi_reference_growth_qa.csv", stringsAsFactors = FALSE)
  stopifnot(
    nrow(growth_qa) >= 22L,
    all(growth_qa$target_retained >= 1000L),
    all(growth_qa$target_met),
    all(growth_qa$target_retained >= ceiling(growth_qa$previous_retained * 1.01))
  )
}
message("Release gates passed for ", marker_id, " ", manifest$release_version)
