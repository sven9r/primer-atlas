#!/usr/bin/env Rscript

.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
source("R/data_layer.R")

old_base <- Sys.getenv("ATLAS_MANIFEST_BASE_URL", unset = NA_character_)
on.exit({
  if (is.na(old_base)) Sys.unsetenv("ATLAS_MANIFEST_BASE_URL") else Sys.setenv(ATLAS_MANIFEST_BASE_URL = old_base)
}, add = TRUE)
Sys.setenv(ATLAS_MANIFEST_BASE_URL = "http://127.0.0.1:9/unavailable")
fallback <- atlas_read_manifest("COI")
stopifnot(fallback$stale, fallback$manifest$marker_id == "COI")

artifact <- fallback$manifest$artifacts$overview
overview <- atlas_read_artifact(artifact)
stopifnot(nrow(overview) == artifact$rows)

build_source <- paste(readLines("scripts/build_release.R", warn = FALSE), collapse = "\n")
publish_source <- paste(readLines("scripts/publish_r2.sh", warn = FALSE), collapse = "\n")
gate_source <- paste(readLines("scripts/validate_release.R", warn = FALSE), collapse = "\n")
stopifnot(
  !grepl("write_json(manifest, previous_pointer", build_source, fixed = TRUE),
  grepl("latest.json", publish_source, fixed = TRUE),
  grepl("new_rows >= 0.90 * old_rows", gate_source, fixed = TRUE),
  grepl("new >= old - 0.05", gate_source, fixed = TRUE),
  grepl("anyDuplicated(taxonomy$accession)", gate_source, fixed = TRUE)
)

message("R2 outage fallback, lazy artifact, and atomic release checks passed.")
