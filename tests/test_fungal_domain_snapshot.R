#!/usr/bin/env Rscript
# Execute the real annotation-snapshot CLI in a clean worker. No primary
# snapshot, panel, index, release pointer or scientific source is modified.
root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(root, ".Rlib"), .libPaths()))
source("R/data_layer.R")
source("R/fungal_versions.R")
worker <- tempfile("fungal-domain-snapshot-worker-")
dir.create(worker)
for (path in c("R/data_layer.R", "R/fungal_versions.R", "scripts/build_fungal_domain_snapshot.R", "scripts/build_fungal_reference_panel.py")) {
  target <- file.path(worker, path)
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(file.path(root, path), target))
}
version_root <- fungal_version_root(worker)
dir.create(version_root, recursive = TRUE)
write_fixture <- function(number) {
  id <- paste0("v", number, "-", sprintf("%012d", number))
  directory <- file.path(version_root, id)
  dir.create(directory)
  artifact <- list()
  frames <- list(catalog = data.frame(primer_name = "TEST", sequence = "ACGT"),
    taxonomy = data.frame(accession = "TEST.1", class = "Synthetic"),
    geometry = data.frame(pair_id = "TEST_PAIR"),
    landmarks = data.frame(accession = "TEST.1", landmark_id = "28S", label = "28S",
      start = 6L, end = 20L, source_note = "Synthetic complete gene fixture"),
    pairs = data.frame(pair_id = "TEST_PAIR", pair_label = "Synthetic fixture"),
    pair_sequences = data.frame(pair_id = "TEST_PAIR", accession = "TEST.1", amplicon_bp = 12L),
    sites = data.frame(accession = "TEST.1", primer_name = "TEST", start = 2L, end = 5L))
  for (name in names(frames)) {
    path <- file.path(directory, paste0(name, ".csv"))
    write.csv(frames[[name]], path, row.names = FALSE)
    artifact[[name]] <- list(path = basename(path), sha256 = atlas_sha256(path), rows = nrow(frames[[name]]))
  }
  for (name in c("sequences", "alignment")) {
    path <- file.path(directory, paste0(name, ".fasta"))
    sequence <- if (name == "sequences") "ACGTACGTACGTACGTACGT" else "ACGT--ACGTACGT--ACGTACGT"
    writeLines(c(">TEST.1", sequence), path)
    artifact[[name]] <- list(path = basename(path), sha256 = atlas_sha256(path), rows = 1L)
  }
  path <- file.path(directory, "run_manifest.json")
  jsonlite::write_json(list(source = "synthetic fixture"), path, auto_unbox = TRUE)
  artifact$run_manifest <- list(path = basename(path), sha256 = atlas_sha256(path))
  manifest <- list(schema_version = "1.0.0", marker_id = "ITS_FUNGAL", id = id,
    version_number = number, reference_accession = "TEST.1", sequence_count = 1L,
    label = paste0("V", number, " fixture"), description = "Synthetic fixture",
    created_at = "2026-10-12", input_sha256 = paste0("fixture", number), artifacts = artifact)
  jsonlite::write_json(manifest, file.path(directory, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
  list(id = id, label = manifest$label, description = manifest$description,
    created_at = manifest$created_at, manifest = file.path(id, "manifest.json"))
}
entries <- lapply(1:2, write_fixture)
jsonlite::write_json(list(schema_version = "1.0.0", marker_id = "ITS_FUNGAL", versions = entries),
  file.path(version_root, "index.json"), auto_unbox = TRUE, pretty = TRUE)
base <- fungal_read_snapshot(entries[[2]]$id, worker, remote = FALSE)
old_files <- list.files(version_root, recursive = TRUE, full.names = TRUE)
old_files <- old_files[basename(old_files) != "index.json"]
old_hashes <- vapply(old_files, atlas_sha256, character(1))
annotation_dir <- file.path(worker, "data/derived/fungal_domain_annotations")
dir.create(annotation_dir, recursive = TRUE)
annotations <- read.csv(file.path(version_root, entries[[2]]$id, "landmarks.csv"), stringsAsFactors = FALSE)
annotations$evidence_status <- "synthetic_gene"
annotations$parent_region <- "28S"
annotations$source_url <- "synthetic test fixture"
annotations$source_start <- NA_integer_
annotations$projection_identity <- NA_real_
domain <- annotations[1, , drop = FALSE]
domain$landmark_id <- "D1"
domain$label <- "D1 synthetic fixture"
domain$start <- 7L
domain$end <- 10L
domain$source_note <- "Synthetic structural-reference projection fixture"
domain$evidence_status <- "reference_projection"
domain$source_start <- 101L
domain$projection_identity <- 0.94
annotations <- rbind(annotations, domain)
write_annotations <- function(value, generated_at = "first fixture build") {
  write.csv(value, file.path(annotation_dir, "landmarks.csv"), row.names = FALSE, na = "")
  jsonlite::write_json(list(source = "synthetic fixture", base_snapshot = entries[[2]]$id, generated_at = generated_at,
    artifacts_sha256 = list(landmarks.csv = atlas_sha256(file.path(annotation_dir, "landmarks.csv"))),
    source_artifacts = list(projection_code = list(repo_path = "scripts/projection_fixture.R",
      sha256 = atlas_sha256(file.path(worker, "scripts/projection_fixture.R"))))),
    file.path(annotation_dir, "run_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
}
writeLines("# Synthetic projection-source fixture", file.path(worker, "scripts/projection_fixture.R"))
writeLines("<synthetic_reference>fixture</synthetic_reference>", file.path(annotation_dir, "reference.xml"))
writeBin(as.raw(c(0L, 1L, 255L)), file.path(annotation_dir, "binary_source.bin"))
write_annotations(annotations)
run_builder <- function(expect_success = TRUE) {
  previous <- getwd()
  setwd(worker)
  output <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), "scripts/build_fungal_domain_snapshot.R",
    env = paste0("R_LIBS=", shQuote(paste(.libPaths(), collapse = .Platform$path.sep))), stdout = TRUE, stderr = TRUE))
  setwd(previous)
  successful <- is.null(attr(output, "status")) || attr(output, "status") == 0L
  if (!identical(successful, expect_success)) stop(paste(output, collapse = "\n"))
  invisible(output)
}
run_builder()
index <- fungal_read_index(worker, remote = FALSE)
stopifnot(nrow(index) == 3L)
updated <- fungal_read_snapshot(tail(index$id, 1), worker, remote = FALSE)
stopifnot(updated$manifest$version_number == 3L,
  identical(updated$manifest$annotation_parent_id, entries[[2]]$id),
  identical(updated$manifest$annotation_retention_pipeline_sha256,
    atlas_sha256(file.path(worker, "scripts/build_fungal_reference_panel.py"))),
  identical(updated$landmarks$evidence_status[updated$landmarks$landmark_id == "D1"], "reference_projection"),
  identical(updated$landmarks$parent_region, rep("28S", 2)),
  is.numeric(updated$landmarks$source_start), is.numeric(updated$landmarks$projection_identity),
  updated$landmarks$source_start[updated$landmarks$landmark_id == "D1"] == 101,
  identical(updated$sequences, base$sequences), identical(updated$alignment, base$alignment),
  isTRUE(all.equal(updated$pairs, base$pairs, check.attributes = FALSE)),
  isTRUE(all.equal(updated$pair_sequences, base$pair_sequences, check.attributes = FALSE)))
for (name in setdiff(names(base$manifest$artifacts), "landmarks")) {
  stopifnot(identical(updated$manifest$artifacts[[name]], base$manifest$artifacts[[name]]))
}
binary <- names(Filter(function(x) x$encoding == "base64", updated$manifest$annotation_source_archives))
stopifnot(length(binary) == 1L,
  identical(jsonlite::base64_dec(updated[[binary]]), as.raw(c(0L, 1L, 255L))))
stopifnot(any(vapply(updated$manifest$annotation_source_archives, function(x) x$original_name == "reference.xml", logical(1))))

# The hosted read path must take typing rules from the frozen inner manifest,
# even when the enclosing release descriptor has no typing metadata.
saved_cache <- atlas_cached_artifact
saved_remote <- Sys.getenv("ATLAS_MANIFEST_BASE_URL", unset = NA_character_)
Sys.setenv(ATLAS_MANIFEST_BASE_URL = "https://annotation-reader-fixture.invalid/releases")
remote_directory <- dirname(file.path(version_root, index$manifest[nrow(index)]))
remote_manifest <- file.path(remote_directory, "manifest.json")
url_paths <- list()
descriptor <- list(id = updated$manifest$id,
  manifest = list(url = "https://annotation-reader-fixture.invalid/manifest.json", sha256 = atlas_sha256(remote_manifest)),
  artifacts = list())
url_paths[[descriptor$manifest$url]] <- remote_manifest
for (name in names(updated$manifest$artifacts)) {
  artifact <- updated$manifest$artifacts[[name]]
  url <- paste0("https://annotation-reader-fixture.invalid/artifacts/", name)
  descriptor$artifacts[[name]] <- list(url = url, sha256 = artifact$sha256)
  url_paths[[url]] <- file.path(remote_directory, artifact$path)
}
atlas_cached_artifact <- function(url, sha256, ...) {
  path <- url_paths[[url]]
  stopifnot(!is.null(path), identical(atlas_sha256(path), sha256))
  path
}
assign(fungal_remote_key(root, updated$manifest$id), descriptor, envir = .fungal_remote_snapshots)
hosted <- fungal_read_snapshot(updated$manifest$id)
stopifnot(identical(hosted$manifest$transport_source, "verified hosted snapshot"),
  identical(hosted$landmarks$parent_region, rep("28S", 2)), is.numeric(hosted$landmarks$start),
  is.numeric(hosted$landmarks$source_start),
  identical(hosted$landmarks$evidence_status[hosted$landmarks$landmark_id == "D1"], "reference_projection"))
atlas_cached_artifact <- saved_cache
if (is.na(saved_remote)) Sys.unsetenv("ATLAS_MANIFEST_BASE_URL") else Sys.setenv(ATLAS_MANIFEST_BASE_URL = saved_remote)
bad_manifest <- updated$manifest
bad_manifest$artifacts$landmarks$character_columns <- "missing_column"
failure <- tryCatch({fungal_snapshot_artifact(bad_manifest, "landmarks", remote_directory); NULL}, error = identity)
stopifnot(inherits(failure, "error"), grepl("headers", conditionMessage(failure)))

run_builder()
stopifnot(identical(fungal_read_index(worker, remote = FALSE)$id, index$id))
write_annotations(annotations, "later replay timestamp")
run_builder()
stopifnot(identical(fungal_read_index(worker, remote = FALSE)$id, index$id))

# New reviewed biological coordinates allocate V4 while retaining both earlier
# scientific snapshots and the first annotation revision byte-for-byte.
annotations$end[annotations$landmark_id == "D1"] <- 11L
write_annotations(annotations)
run_builder()
new_index <- fungal_read_index(worker, remote = FALSE)
stopifnot(nrow(new_index) == 4L,
  fungal_read_snapshot(tail(new_index$id, 1), worker, remote = FALSE)$manifest$version_number == 4L)
hash <- atlas_sha256(file.path(version_root, "index.json"))
bad <- annotations
bad$end[bad$landmark_id == "D1"] <- 21L
write_annotations(bad)
run_builder(FALSE)
stopifnot(identical(atlas_sha256(file.path(version_root, "index.json")), hash))
bad <- annotations
bad$source_note[bad$landmark_id == "28S"] <- "Changed historical annotation"
write_annotations(bad)
run_builder(FALSE)
stopifnot(identical(atlas_sha256(file.path(version_root, "index.json")), hash))
write_annotations(annotations)
provenance_file <- file.path(annotation_dir, "run_manifest.json")
wrong_base <- jsonlite::read_json(provenance_file, simplifyVector = FALSE)
wrong_base$base_snapshot <- entries[[1]]$id
jsonlite::write_json(wrong_base, provenance_file, auto_unbox = TRUE, pretty = TRUE)
run_builder(FALSE)
stopifnot(identical(atlas_sha256(file.path(version_root, "index.json")), hash))
write_annotations(annotations)
old_alignment <- file.path(version_root, entries[[2]]$id, "alignment.fasta")
saved_bytes <- readBin(old_alignment, "raw", n = file.info(old_alignment)$size)
cat("corrupt source alignment", file = old_alignment, append = TRUE)
run_builder(FALSE)
writeBin(saved_bytes, old_alignment)
stopifnot(identical(atlas_sha256(file.path(version_root, "index.json")), hash))
cat("changed source bytes", file = file.path(worker, "scripts/projection_fixture.R"), append = TRUE)
run_builder(FALSE)
stopifnot(identical(atlas_sha256(file.path(version_root, "index.json")), hash),
  identical(vapply(old_files, atlas_sha256, character(1)), old_hashes))
unlink(worker, recursive = TRUE)
message("Annotation CLI creates immutable V3/V4, preserves scientific artifact bytes, archives source/code evidence, reuses timestamp-only replays, and rejects invalid coordinates or changed history/source proofs.")
