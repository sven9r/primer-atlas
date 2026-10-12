#!/usr/bin/env Rscript
# Exercise the actual restore CLI against bundled annotation history and an
# older/corrupt mocked hosted boundary. No primary source or pointer is changed.
root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(root, ".Rlib"), .libPaths()))
source("R/data_layer.R")
source("R/fungal_versions.R")
shipped_index <- fungal_read_index(remote = FALSE)
manifests <- lapply(shipped_index$manifest, function(path) {
  jsonlite::read_json(file.path(fungal_version_root(), path), simplifyVector = FALSE)
})
annotations <- which(vapply(manifests, function(value) !is.null(value$annotation_parent_id), logical(1)))
stopifnot(length(annotations) > 0L)
latest <- annotations[which.max(vapply(manifests[annotations], `[[`, integer(1), "version_number"))]
v3_id <- shipped_index$id[latest]
v3 <- fungal_read_snapshot(v3_id, remote = FALSE)
v2_id <- v3$manifest$annotation_parent_id
v2_row <- match(v2_id, shipped_index$id)
stopifnot(!is.na(v2_row))
original_files <- list.files(fungal_version_root(), recursive = TRUE, full.names = TRUE)
original_hashes <- vapply(original_files, atlas_sha256, character(1))
worker <- tempfile("fungal-restore-selection-worker-")
dir.create(worker)
for (path in c("R/functions.R", "R/data_layer.R", "R/fungal_versions.R", "scripts/build_fungal_snapshots.R")) {
  destination <- file.path(worker, path)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(file.path(root, path), destination))
}
dir.create(file.path(worker, "data/pinned/ITS_FUNGAL"), recursive = TRUE)
stopifnot(file.copy(fungal_version_root(), file.path(worker, "data/pinned/ITS_FUNGAL"), recursive = TRUE))
fixture <- file.path(worker, "hosted_fixture")
dir.create(fixture)
served <- file.path(fixture, "v2")
dir.create(served)
v2_directory <- dirname(file.path(fungal_version_root(), shipped_index$manifest[v2_row]))
stopifnot(file.copy(list.files(v2_directory, full.names = TRUE), served, recursive = TRUE))
hosted <- list(id = v2_id, label = manifests[[v2_row]]$label,
  description = manifests[[v2_row]]$description, created_at = manifests[[v2_row]]$created_at,
  manifest = list(url = "https://restore-fixture.invalid/v2/manifest.json",
    sha256 = atlas_sha256(file.path(served, "manifest.json"))), artifacts = list())
paths <- list()
paths[[hosted$manifest$url]] <- file.path(served, "manifest.json")
for (name in names(manifests[[v2_row]]$artifacts)) {
  descriptor <- manifests[[v2_row]]$artifacts[[name]]
  url <- paste0("https://restore-fixture.invalid/v2/", descriptor$path)
  hosted$artifacts[[name]] <- list(url = url, sha256 = descriptor$sha256)
  paths[[url]] <- file.path(served, descriptor$path)
}
jsonlite::write_json(paths, file.path(fixture, "paths.json"), auto_unbox = TRUE)
# Mock only transport. The CLI, immutable manifest/artifact checks, history
# installation, typing declarations and selected prior-panel copies are real.
cat('\n\natlas_read_manifest <- function(...) {\n  state <- jsonlite::read_json("hosted_fixture/state.json", simplifyVector = FALSE)\n  list(manifest = state$manifest, stale = state$stale, source = "mocked hosted boundary")\n}\natlas_cached_artifact <- function(url, sha256 = NULL, ...) {\n  paths <- jsonlite::read_json("hosted_fixture/paths.json", simplifyVector = FALSE)\n  path <- paths[[url]]\n  if (is.null(path) || !file.exists(path) || !identical(atlas_sha256(path), sha256)) stop("Mocked hosted artifact checksum failed")\n  path\n}\n',
  file = file.path(worker, "R/data_layer.R"), append = TRUE)
set_state <- function(descriptors = list(), stale = FALSE) {
  jsonlite::write_json(list(manifest = list(marker_id = "ITS_FUNGAL", fungal_snapshots = descriptors), stale = stale),
    file.path(fixture, "state.json"), auto_unbox = TRUE)
}
run_restore <- function(expected_success = TRUE) {
  previous <- getwd()
  setwd(worker)
  output <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    c("scripts/build_fungal_snapshots.R", "restore"),
    env = c(paste0("R_LIBS=", shQuote(paste(.libPaths(), collapse = .Platform$path.sep))),
            "ATLAS_MANIFEST_BASE_URL=https://restore-fixture.invalid/releases"), stdout = TRUE, stderr = TRUE))
  setwd(previous)
  success <- is.null(attr(output, "status")) || attr(output, "status") == 0L
  if (!identical(success, expected_success)) stop(paste(output, collapse = "\n"))
  output
}
panel <- file.path(worker, "data/derived/fungal_previous_reference_panel")
assert_v3_panel <- function() {
  mapping <- c(sequences = "panel.fasta", taxonomy = "panel_manifest.csv", landmarks = "landmarks.csv", run_manifest = "run_manifest.json")
  directory <- dirname(file.path(fungal_version_root(worker), shipped_index$manifest[latest]))
  for (name in names(mapping)) {
    descriptor <- v3$manifest$artifacts[[name]]
    stopifnot(identical(atlas_sha256(file.path(panel, mapping[[name]])), descriptor$sha256),
              identical(readBin(file.path(panel, mapping[[name]]), "raw", n = file.info(file.path(panel, mapping[[name]]))$size),
                readBin(file.path(directory, descriptor$path), "raw", n = file.info(file.path(directory, descriptor$path))$size)))
  }
  descriptor <- v3$manifest$artifacts$landmarks
  descriptor$local_path <- file.path(panel, "landmarks.csv")
  restored_landmarks <- fungal_snapshot_read_artifact(descriptor)
  domain <- restored_landmarks$landmark_type %in% "divergent_domain"
  stopifnot(sum(domain) == sum(v3$landmarks$landmark_type %in% "divergent_domain"),
    identical(restored_landmarks$parent_landmark[domain], v3$landmarks$parent_landmark[domain]),
    is.character(restored_landmarks$parent_landmark), is.numeric(restored_landmarks$start),
    is.numeric(restored_landmarks$source_feature_start),
    all(restored_landmarks$evidence_status[domain] == "mapped_reference_domain"),
    all(restored_landmarks$accession %in% v3$taxonomy$accession))
}

# An older hosted V2 must not take precedence over bundled reviewed V3.
set_state(setNames(list(hosted), v2_id))
output <- run_restore()
stopifnot(any(grepl(v3_id, output, fixed = TRUE)))
assert_v3_panel()

# A catalog-only latest still retains the complete bundled V3 panel, including
# all provenance columns, numeric boundaries and accession-scoped domains.
unlink(panel, recursive = TRUE)
set_state()
output <- run_restore()
stopifnot(any(grepl(v3_id, output, fixed = TRUE)))
assert_v3_panel()

# Restore a missing older snapshot from checked hosted bytes, then choose V3.
index <- fungal_read_index(worker, remote = FALSE)
index <- index[index$id != v2_id, , drop = FALSE]
jsonlite::write_json(list(schema_version = "1.0.0", marker_id = "ITS_FUNGAL",
  versions = lapply(seq_len(nrow(index)), function(i) as.list(index[i, , drop = FALSE]))),
  file.path(fungal_version_root(worker), "index.json"), auto_unbox = TRUE)
unlink(file.path(fungal_version_root(worker), v2_id), recursive = TRUE)
set_state(setNames(list(hosted), v2_id))
invisible(run_restore())
stopifnot(v2_id %in% fungal_read_index(worker, remote = FALSE)$id)
assert_v3_panel()

# Corrupt newer bundled evidence cannot silently fall back to older hosted V2.
landmark_path <- file.path(fungal_version_root(worker), v3_id, v3$manifest$artifacts$landmarks$path)
bytes <- readBin(landmark_path, "raw", n = file.info(landmark_path)$size)
cat("corrupt bytes", file = landmark_path, append = TRUE)
output <- run_restore(FALSE)
stopifnot(any(grepl("checksum", output, ignore.case = TRUE)))
writeBin(bytes, landmark_path)
assert_v3_panel()

# Coherent checksums do not excuse a landmark assigned to a foreign accession.
# This mutation exists only inside the disposable fixture, then is restored.
manifest_path <- file.path(fungal_version_root(worker), v3_id, "manifest.json")
manifest_bytes <- readBin(manifest_path, "raw", n = file.info(manifest_path)$size)
malformed <- read.csv(landmark_path, check.names = FALSE, stringsAsFactors = FALSE)
malformed$accession[which(malformed$landmark_type %in% "divergent_domain")[1]] <- "FOREIGN.1"
write.csv(malformed, landmark_path, row.names = FALSE, na = "")
manifest <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
manifest$artifacts$landmarks$sha256 <- atlas_sha256(landmark_path)
jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE)
output <- run_restore(FALSE)
stopifnot(any(grepl("accession ownership", output)))
writeBin(bytes, landmark_path)
writeBin(manifest_bytes, manifest_path)
assert_v3_panel()

# Stale hosted state remains fail-closed even if bundled evidence is available.
set_state(stale = TRUE)
output <- run_restore(FALSE)
stopifnot(any(grepl("unreachable", output)))
assert_v3_panel()
stopifnot(identical(vapply(original_files, atlas_sha256, character(1)), original_hashes))
unlink(worker, recursive = TRUE)
message("Restore CLI retains bundled annotation history over older/catalog-only hosting, restores checked history, preserves provenance typing/ownership, and rejects corrupt or stale state.")
