#!/usr/bin/env Rscript

root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(root, ".Rlib"), .libPaths()))
suppressPackageStartupMessages({library(jsonlite); library(nanoparquet); library(digest)})
source("R/data_layer.R")

release_id <- Sys.getenv("RELEASE_ID", format(Sys.time(), tz = "UTC", "%Y%m%dT%H%M%SZ"))
if (!grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", release_id)) stop("Invalid release identity")
data_base <- sub("/+$", "", Sys.getenv("ATLAS_DATA_BASE_URL", "https://data.example.invalid/primer-atlas"))
markers <- read.csv("data/catalog/markers.csv", stringsAsFactors = FALSE)

artifact_specs <- list(
  COI = c(
    overview = "data/derived/order_pair_summary.csv",
    pair_geometry = "data/derived/marker_pair_geometry.csv",
    pair_templates = "data/derived/pair_template_scores.csv",
    primer_positions = "data/derived/primer_position_scores.csv",
    reference_growth_qa = "data/provenance/coi_reference_growth_qa.csv",
    taxonomy = "data/derived/gurten2026_centroid_taxonomy.csv",
    geography = "data/derived/reference_geography.csv",
    families = "data/derived/claimed_primer_family_summary.csv",
    subfamilies = "data/derived/claimed_primer_subfamily_summary.csv",
    genera = "data/derived/claimed_primer_genus_summary.csv"
  ),
  ITS_FUNGAL = c(
    pair_geometry = "data/derived/marker_pair_geometry.csv",
    primer_source = "data/catalog/unite_primers.csv"
  )
)

write_marker <- function(marker_id) {
  release_dir <- file.path("data", "releases", marker_id, release_id)
  if (dir.exists(release_dir) && length(list.files(release_dir, all.files = FALSE))) {
    stop("Refusing to overwrite an existing immutable release: ", release_id)
  }
  dir.create(release_dir, recursive = TRUE, showWarnings = FALSE)
  artifacts <- list()
  specs <- artifact_specs[[marker_id]]
  if (marker_id == "COI") {
    beeprime <- "data/derived/claimed_beeprime_centroid_scores.csv.gz"
    if (file.exists(beeprime)) specs <- c(specs, exact_beeprime = beeprime)
    full <- list.files("data/derived", pattern = "^full_order_.*_sequence_scores[.]csv[.]gz$", full.names = TRUE)
    names(full) <- sub("_sequence_scores[.]csv[.]gz$", "", basename(full))
    specs <- c(specs, full)
    expanded_specs <- c(
      beeprime_expanded_scores = "data/derived/beeprime_expanded_scores.csv.gz",
      beeprime_expanded_geography = "data/derived/beeprime_expanded_geography.csv",
      beeprime_expanded_taxonomy_audit = "data/provenance/beeprime_expanded_taxonomy_audit.csv",
      beeprime_expanded_growth_qa = "data/provenance/beeprime_expanded_growth_qa.csv",
      beeprime_expanded_alignment = "data/derived/beeprime_expanded_alignment.fasta",
      beeprime_expanded_added_sequences = "data/derived/beeprime_expanded_added_sequences.fasta"
    )
    expanded_specs <- expanded_specs[file.exists(expanded_specs)]
    specs <- c(specs, expanded_specs)
  }
  for (name in names(specs)) {
    source <- specs[[name]]
    if (!file.exists(source)) next
    if (grepl("[.]fasta$", source, ignore.case = TRUE)) {
      destination <- file.path(release_dir, paste0(name, ".fasta"))
      if (!file.copy(source, destination, overwrite = TRUE)) {
        stop("Could not stage reference FASTA artifact: ", source)
      }
      artifact_rows <- sum(startsWith(readLines(destination, warn = FALSE), ">"))
    } else {
      x <- read.csv(source, stringsAsFactors = FALSE, check.names = FALSE)
      if (name == "pair_geometry") x <- x[x$marker_id == marker_id, , drop = FALSE]
      destination <- file.path(release_dir, paste0(name, ".parquet"))
      nanoparquet::write_parquet(x, destination, compression = "gzip")
      artifact_rows <- nrow(x)
    }
    relative <- file.path("releases", marker_id, release_id, basename(destination))
    artifacts[[name]] <- list(
      url = paste(data_base, relative, sep = "/"),
      local_path = destination,
      bytes = unname(file.info(destination)$size),
      rows = artifact_rows,
      sha256 = digest::digest(file = destination, algo = "sha256", serialize = FALSE)
    )
  }
  fungal_snapshots <- list()
  if (marker_id == "ITS_FUNGAL") {
    source("R/fungal_versions.R")
    index <- fungal_read_index(root, remote = FALSE)
    if (!nrow(index)) stop("Build frozen fungal snapshots before creating an ITS release")
    snapshot_root <- fungal_version_root(root)
    add_snapshot_artifact <- function(name, source_path, relative_path, rows = NULL) {
      target <- file.path(release_dir, relative_path)
      dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
      if (!file.copy(source_path, target, overwrite = FALSE)) stop("Could not stage fungal snapshot artifact")
      artifact <- list(
        url = paste(data_base, "releases", marker_id, release_id, relative_path, sep = "/"),
        local_path = target, bytes = unname(file.info(target)$size),
        sha256 = atlas_sha256(target)
      )
      if (!is.null(rows)) artifact$rows <- rows
      artifacts[[name]] <<- artifact
      artifact
    }
    add_snapshot_artifact("fungal_snapshot_index", file.path(snapshot_root, "index.json"), "fungal_versions/index.json")
    for (i in seq_len(nrow(index))) {
      id <- index$id[i]
      source_manifest <- file.path(snapshot_root, index$manifest[i])
      snapshot <- read_json(source_manifest, simplifyVector = FALSE)
      snapshot_artifacts <- list()
      for (name in names(snapshot$artifacts)) {
        # This verifies the frozen input/output before packaging it.
        fungal_snapshot_artifact(snapshot, name, dirname(source_manifest))
        descriptor <- snapshot$artifacts[[name]]
        relative <- file.path("fungal_versions", id, descriptor$path)
        snapshot_artifacts[[name]] <- add_snapshot_artifact(
          paste("fungal", id, name, sep = "_"),
          file.path(dirname(source_manifest), descriptor$path), relative, descriptor$rows
        )
      }
      manifest_artifact <- add_snapshot_artifact(
        paste0("fungal_", id, "_manifest"), source_manifest,
        file.path("fungal_versions", id, "manifest.json")
      )
      fungal_snapshots[[id]] <- list(
        id = id, label = snapshot$label, description = snapshot$description,
        created_at = snapshot$created_at, input_sha256 = snapshot$input_sha256,
        manifest = manifest_artifact, artifacts = snapshot_artifacts
      )
    }
  }
  marker <- markers[markers$marker_id == marker_id, , drop = FALSE]
  previous_pointer <- file.path("data", "releases", marker_id, "latest.json")
  previous_state <- tryCatch(atlas_read_manifest(marker_id, root), error = function(e) NULL)
  previous <- if (!is.null(previous_state)) previous_state$manifest$release_version else NULL
  reference_panels <- list()
  if (marker_id == "COI" && all(c(
    "beeprime_expanded_scores", "beeprime_expanded_geography",
    "beeprime_expanded_taxonomy_audit", "beeprime_expanded_growth_qa",
    "beeprime_expanded_alignment", "beeprime_expanded_added_sequences"
  ) %in% names(artifacts))) {
    expanded_scores <- read.csv("data/derived/beeprime_expanded_scores.csv.gz",
                                stringsAsFactors = FALSE, check.names = FALSE)
    reference_panels$beeprime_expanded <- list(
      pair_id = "BEEPRIME",
      reference_mode = "expanded",
      taxonomic_scope = "bees",
      publication_baseline = "gurten2026_suppl5",
      version = paste0(release_id, "-BeePrime"),
      added_sequences = sum(expanded_scores$reference_origin == "added"),
      sequence_scores = artifacts$beeprime_expanded_scores,
      geography = artifacts$beeprime_expanded_geography,
      alignment = artifacts$beeprime_expanded_alignment,
      added_sequences_fasta = artifacts$beeprime_expanded_added_sequences,
      taxonomy_audit = artifacts$beeprime_expanded_taxonomy_audit,
      growth_qa = artifacts$beeprime_expanded_growth_qa
    )
  }
  manifest <- list(
    schema_version = "1.0.0", release_version = release_id,
    marker_id = marker_id, released_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    upstream_sources = list(list(name = marker$source_name, url = marker$source_url)),
    coordinate_reference = marker$coordinate_reference,
    reference_version = marker$reference_version,
    completeness = list(
      taxonomy_resolved_fraction = if (marker_id == "COI") mean(nzchar(read.csv("data/derived/gurten2026_centroid_taxonomy.csv")$order), na.rm = TRUE) else NA,
      geography_resolved_fraction = if (marker_id == "COI") mean(nzchar(read.csv("data/derived/reference_geography.csv")$country_or_territory), na.rm = TRUE) else NA
    ),
    previous_successful_release = previous,
    artifacts = artifacts,
    reference_panels = reference_panels,
    fungal_snapshots = fungal_snapshots,
    code_commit = tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[1], error = function(e) "unknown"),
    builder_sha256 = atlas_sha256("scripts/build_release.R"),
    snapshot_builder_sha256 = if (marker_id == "ITS_FUNGAL") atlas_sha256("scripts/build_fungal_snapshots.R") else NULL
  )
  versioned_manifest <- file.path(release_dir, "manifest.json")
  write_json(manifest, versioned_manifest, auto_unbox = TRUE, pretty = TRUE, na = "null")
  message(marker_id, ": ", length(artifacts), " artifacts in ", release_dir)
}

requested <- commandArgs(trailingOnly = TRUE)
if (!length(requested)) requested <- intersect(c("COI", "ITS_FUNGAL"), markers$marker_id)
for (marker_id in requested) {
  if (!marker_id %in% names(artifact_specs)) stop("No release specification for ", marker_id)
  write_marker(marker_id)
}
