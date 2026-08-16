#!/usr/bin/env Rscript

root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(root, ".Rlib"), .libPaths()))
suppressPackageStartupMessages({library(jsonlite); library(nanoparquet); library(digest)})
source("R/data_layer.R")

release_id <- Sys.getenv("RELEASE_ID", format(Sys.Date(), "%Y-%m-%d"))
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
  dir.create(release_dir, recursive = TRUE, showWarnings = FALSE)
  artifacts <- list()
  specs <- artifact_specs[[marker_id]]
  if (marker_id == "COI") {
    beeprime <- "data/derived/claimed_beeprime_centroid_scores.csv.gz"
    if (file.exists(beeprime)) specs <- c(specs, exact_beeprime = beeprime)
    full <- list.files("data/derived", pattern = "^full_order_.*_sequence_scores[.]csv[.]gz$", full.names = TRUE)
    names(full) <- sub("_sequence_scores[.]csv[.]gz$", "", basename(full))
    specs <- c(specs, full)
  }
  for (name in names(specs)) {
    source <- specs[[name]]
    if (!file.exists(source)) next
    x <- read.csv(source, stringsAsFactors = FALSE, check.names = FALSE)
    if (name == "pair_geometry") x <- x[x$marker_id == marker_id, , drop = FALSE]
    destination <- file.path(release_dir, paste0(name, ".parquet"))
    nanoparquet::write_parquet(x, destination, compression = "gzip")
    relative <- file.path("releases", marker_id, release_id, basename(destination))
    artifacts[[name]] <- list(
      url = paste(data_base, relative, sep = "/"),
      local_path = destination,
      bytes = unname(file.info(destination)$size),
      rows = nrow(x),
      sha256 = digest::digest(file = destination, algo = "sha256", serialize = FALSE)
    )
  }
  marker <- markers[markers$marker_id == marker_id, , drop = FALSE]
  previous_pointer <- file.path("data", "releases", marker_id, "latest.json")
  previous_state <- tryCatch(atlas_read_manifest(marker_id, root), error = function(e) NULL)
  previous <- if (!is.null(previous_state)) previous_state$manifest$release_version else NULL
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
    artifacts = artifacts
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
