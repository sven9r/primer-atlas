#!/usr/bin/env Rscript

# Exercise the application reader against this candidate, before publishing it.
# Merely starting the app can succeed while it falls back to an older release.
root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(root, ".Rlib"), .libPaths()))
source("R/functions.R")
source("R/data_layer.R")
source("R/coi_retention.R")
source("R/coi_scoring.R")
source("R/coi_reference.R")

release_id <- Sys.getenv("RELEASE_ID", "")
if (!nzchar(release_id)) stop("RELEASE_ID is required for candidate runtime verification")
path <- atlas_manifest_path("COI", root, release_id)
state <- list(manifest = jsonlite::read_json(path, simplifyVector = FALSE),
              source = path, stale = FALSE)
reference <- coi_reference_load(state, root)
stopifnot(identical(reference$mode, "full_release"),
          identical(reference$version, release_id),
          identical(reference$alignment_mode, "full"),
          reference$n_references == state$manifest$coi_retained_state$retained_records,
          reference$n_verified_marker >= 10000L,
          length(reference$counts) == 22L,
          all(reference$counts >= 100L))

pair_id <- unique(reference$overview$pair_id)[1]
positions <- coi_reference_pair_positions(reference, pair_id)
primers <- read.csv(file.path(root, "data/primers.csv"), stringsAsFactors = FALSE)
primers <- primers[primers$pair_id == pair_id, , drop = FALSE]
expected <- reference$n_references * sum(nchar(vapply(primers$sequence, clean_sequence, character(1))))
stopifnot(nrow(positions) == expected,
          setequal(unique(positions$order), names(reference$counts)),
          all(positions$pair_id == pair_id))
message("Candidate application reader passed: ", release_id, "; ", reference$n_references,
        " retained / ", reference$n_verified_marker, " verified COI records; ",
        nrow(positions), " checked position rows for ", pair_id)
