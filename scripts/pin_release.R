#!/usr/bin/env Rscript

.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
suppressPackageStartupMessages(library(jsonlite))
marker_id <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(marker_id) || !nzchar(marker_id)) stop("Usage: pin_release.R MARKER_ID")
release_root <- file.path("data", "releases", marker_id)
releases <- list.files(release_root, pattern = "^[0-9].*", full.names = FALSE)
release_id <- tail(sort(releases), 1)
manifest_path <- file.path(release_root, release_id, "manifest.json")
manifest <- read_json(manifest_path, simplifyVector = FALSE)
compact <- c("overview", "pair_geometry", "pair_templates", "primer_positions", "taxonomy", "geography", "families", "subfamilies", "genera", "primer_source")
destination <- file.path("data", "pinned", marker_id, release_id)
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
for (name in names(manifest$artifacts)) {
  artifact <- manifest$artifacts[[name]]
  if (name %in% compact) {
    target <- file.path(destination, basename(artifact$local_path))
    if (!file.copy(artifact$local_path, target, overwrite = TRUE)) stop("Could not pin ", name)
    artifact$local_path <- target
  } else {
    artifact$local_path <- NULL
  }
  manifest$artifacts[[name]] <- artifact
}
dir.create(release_root, recursive = TRUE, showWarnings = FALSE)
write_json(manifest, file.path(release_root, "latest.json"), auto_unbox = TRUE, pretty = TRUE, na = "null")
message("Pinned compact fallback for ", marker_id, " ", release_id)
