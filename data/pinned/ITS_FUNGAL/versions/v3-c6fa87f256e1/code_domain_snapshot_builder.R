#!/usr/bin/env Rscript
# Add reviewed annotations to a new immutable snapshot; never reevaluate or
# modify the source snapshot's sequences, alignment, catalog, or primer results.
root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(root, ".Rlib"), .libPaths()))
source("R/data_layer.R")
source("R/fungal_versions.R")

options <- list(base_id = Sys.getenv("FUNGAL_DOMAIN_BASE_ID", ""),
  annotation_dir = Sys.getenv("FUNGAL_DOMAIN_ANNOTATION_DIR", "data/derived/fungal_domain_annotations"),
  label = Sys.getenv("FUNGAL_DOMAIN_SNAPSHOT_LABEL", "LSU domain annotations"),
  description = Sys.getenv("FUNGAL_DOMAIN_SNAPSHOT_DESCRIPTION", ""))
args <- commandArgs(trailingOnly = TRUE)
while (length(args)) {
  if (length(args) < 2L || !args[1] %in% c("--base-id", "--annotation-dir", "--label", "--description")) {
    stop("Use --base-id ID, --annotation-dir DIR, --label TEXT or --description TEXT")
  }
  key <- gsub("-", "_", sub("^--", "", args[1]))
  options[[key]] <- args[2]
  args <- args[-c(1L, 2L)]
}

version_root <- fungal_version_root(root)
index_path <- file.path(version_root, "index.json")
index <- fungal_read_index(root, remote = FALSE)
if (!nrow(index)) stop("No frozen fungal source snapshot is installed")
index_hash <- atlas_sha256(index_path)
manifests <- lapply(index$manifest, function(path) {
  jsonlite::read_json(file.path(version_root, path), simplifyVector = FALSE)
})
if (!nzchar(options$base_id)) {
  eligible <- which(vapply(manifests, function(value) {
    is.null(value$annotation_parent_id) && value$version_number > 1L
  }, logical(1)))
  if (!length(eligible)) stop("No expanded source snapshot is available; select --base-id explicitly")
  options$base_id <- index$id[tail(eligible, 1L)]
}
base_row <- match(options$base_id, index$id)
if (is.na(base_row)) stop("Unknown source snapshot: ", options$base_id)
base_path <- file.path(version_root, index$manifest[base_row])
base_directory <- dirname(base_path)
base <- fungal_read_snapshot(options$base_id, root, remote = FALSE)
if (is.null(base$pairs) || is.null(base$pair_sequences) || is.null(base$sites)) {
  stop("Annotation-only snapshots require already frozen primer results")
}
base_hash <- atlas_sha256(base_path)
annotation_dir <- normalizePath(options$annotation_dir, mustWork = TRUE)
landmark_path <- file.path(annotation_dir, "landmarks.csv")
provenance_path <- file.path(annotation_dir, "run_manifest.json")
if (!file.exists(landmark_path) || !file.exists(provenance_path)) {
  stop("Reviewed annotation inputs require landmarks.csv and run_manifest.json")
}
annotations <- read.csv(landmark_path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("accession", "landmark_id", "label", "start", "end", "source_note")
if (!all(required %in% names(annotations))) stop("Annotation table lacks accession, boundaries, labels or source notes")
if (!nrow(annotations)) stop("Reviewed annotation table is empty")
if (anyDuplicated(names(annotations))) stop("Annotation CSV repeats a column name")
if (anyNA(annotations$label) || anyNA(annotations$source_note) ||
    any(!nzchar(annotations$label)) || any(!nzchar(annotations$source_note))) {
  stop("Every annotation requires a label and source or unavailable-evidence note")
}
provenance <- jsonlite::read_json(provenance_path, simplifyVector = FALSE)
if (!is.null(provenance$base_snapshot) && !identical(provenance$base_snapshot, options$base_id)) {
  stop("Reviewed annotations belong to a different base snapshot")
}
for (relative in names(provenance$artifacts_sha256)) {
  if (grepl("(^/|(^|/)\\.\\.(/|$))", relative)) stop("Invalid reviewed annotation artifact path")
  path <- file.path(annotation_dir, relative)
  if (!file.exists(path) || !identical(atlas_sha256(path), provenance$artifacts_sha256[[relative]])) {
    stop("Annotation inputs differ from their reviewed provenance: ", relative)
  }
}

accession <- function(x) sub("^_R_", "", sub("[[:space:]].*$", "", x))
sequences <- base$sequences
names(sequences) <- accession(names(sequences))
key <- function(x) paste(x$accession, x$landmark_id, sep = "|")
if (anyNA(annotations$accession) || anyNA(annotations$landmark_id) ||
    any(!nzchar(annotations$accession)) || any(!nzchar(annotations$landmark_id)) ||
    anyDuplicated(key(annotations)) || any(!annotations$accession %in% names(sequences))) {
  stop("Annotation identities must be unique and belong to the frozen source sequences")
}
old_landmarks <- read.csv(file.path(base_directory, base$manifest$artifacts$landmarks$path),
  stringsAsFactors = FALSE, check.names = FALSE)
if (!"accession" %in% names(old_landmarks)) old_landmarks$accession <- base$manifest$reference_accession
matched <- match(key(old_landmarks), key(annotations))
if (anyNA(matched) || !all(names(old_landmarks) %in% names(annotations))) {
  stop("Annotation update must retain every earlier gene/region annotation")
}
for (column in names(old_landmarks)) {
  before <- as.character(old_landmarks[[column]])
  after <- as.character(annotations[[column]][matched])
  before[is.na(before)] <- ""
  after[is.na(after)] <- ""
  if (!identical(before, after)) stop("Earlier frozen landmark changed: ", column)
}
for (i in seq_len(nrow(annotations))) {
  row <- annotations[i, , drop = FALSE]
  bounds <- suppressWarnings(as.numeric(c(row$start, row$end)))
  if (all(is.na(bounds))) {
    raw_bounds <- as.character(c(row$start, row$end))
    if (any(!is.na(raw_bounds) & !trimws(raw_bounds) %in% c("", "NA"))) {
      stop("Non-numeric annotation boundaries cannot be treated as unavailable")
    }
    statuses <- unlist(row[intersect(c("evidence_status", "annotation_status", "status"), names(row))], use.names = FALSE)
    if (!any(tolower(statuses) %in% c("unavailable", "unresolved", "not_mapped")) ||
        is.na(row$source_note) || !nzchar(row$source_note)) {
      stop("Missing annotation boundaries require an explicit unavailable status and reason")
    }
  } else if (any(!is.finite(bounds)) || any(bounds != floor(bounds)) ||
             bounds[1] < 1L || bounds[2] < bounds[1] ||
             bounds[2] > nchar(sequences[[row$accession]])) {
    stop("Annotation boundary is outside its exact raw source template: ", row$accession, " / ", row$landmark_id)
  }
}

# Freeze local source proofs and any generator code explicitly named by the
# reviewed provenance. Unknown text formats are stored as .txt so the existing
# reader can inspect them; binary sources are recoverable base64 text archives.
sources <- list()
add_source <- function(path, original_name, expected = NULL) {
  path <- normalizePath(path, mustWork = TRUE)
  hash <- atlas_sha256(path)
  if (!is.null(expected) && !identical(tolower(hash), tolower(expected))) stop("Annotation source checksum failed: ", original_name)
  sources[[length(sources) + 1L]] <<- list(path = path, original_name = original_name, sha256 = hash)
}
files <- sort(list.files(annotation_dir, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE))
files <- files[!file.info(files)$isdir]
for (path in files) {
  relative <- substring(path, nchar(annotation_dir) + 2L)
  if (!relative %in% c("landmarks.csv", "run_manifest.json")) add_source(path, relative)
}
for (name in names(provenance$source_artifacts)) {
  descriptor <- provenance$source_artifacts[[name]]
  relative <- if (!is.null(descriptor$repo_path)) descriptor$repo_path else descriptor$path
  if (is.null(relative) || !nzchar(relative) || grepl("(^/|(^|/)\\.\\.(/|$))", relative)) stop("Annotation source paths must be safe relative paths")
  source_root <- if (!is.null(descriptor$repo_path)) root else annotation_dir
  add_source(file.path(source_root, relative), name, descriptor$sha256)
}
if (length(sources)) {
  identities <- vapply(sources, function(x) paste(x$original_name, x$sha256, sep = "|"), character(1))
  sources <- sources[!duplicated(identities)]
  sources <- sources[order(vapply(sources, `[[`, character(1), "original_name"))]
}
stable_provenance <- function(x) {
  if (!is.list(x)) return(x)
  volatile <- c("created_at", "built_at", "generated_at", "acquired_at", "checked_at", "timestamp", "elapsed_seconds")
  if (!is.null(names(x))) x <- x[!names(x) %in% volatile]
  lapply(x, stable_provenance)
}
script_path <- file.path(root, "scripts/build_fungal_domain_snapshot.R")
source_identity <- digest::digest(jsonlite::toJSON(stable_provenance(provenance), auto_unbox = TRUE, null = "null"), algo = "sha256", serialize = FALSE)
input_hash <- digest::digest(paste(c(base_hash, atlas_sha256(landmark_path), source_identity,
  vapply(sources, function(x) paste(x$original_name, x$sha256, sep = "|"), character(1))), collapse = "|"), algo = "sha256", serialize = FALSE)
existing <- which(vapply(manifests, function(x) identical(x$input_sha256, input_hash), logical(1)))
if (length(existing)) {
  fungal_read_snapshot(index$id[existing[1]], root, remote = FALSE)
  message("Identical reviewed annotation snapshot already exists: ", index$id[existing[1]])
  quit(status = 0L)
}
number <- max(vapply(manifests, function(x) as.integer(x$version_number), integer(1))) + 1L
id <- paste0("v", number, "-", substr(input_hash, 1L, 12L))
destination <- file.path(version_root, id)
if (dir.exists(destination)) stop("Refusing to overwrite immutable annotation snapshot: ", id)
staging <- tempfile(".domain-snapshot-", tmpdir = version_root)
dir.create(staging)
published <- FALSE
cleanup <- function() {
  unlink(staging, recursive = TRUE)
  if (!published && dir.exists(destination)) unlink(destination, recursive = TRUE)
}
tryCatch({
  artifacts <- base$manifest$artifacts
  for (name in names(artifacts)) {
    old <- file.path(base_directory, artifacts[[name]]$path)
    new <- file.path(staging, artifacts[[name]]$path)
    dir.create(dirname(new), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(old, new) || !identical(atlas_sha256(new), artifacts[[name]]$sha256)) stop("Could not preserve source artifact bytes: ", name)
  }
  add_artifact <- function(name, path, filename, rows = NULL) {
    target <- file.path(staging, filename)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(path, target, overwrite = TRUE)) stop("Could not freeze reviewed annotation artifact")
    descriptor <- list(path = filename, sha256 = atlas_sha256(target), bytes = unname(file.info(target)$size))
    if (!is.null(rows)) descriptor$rows <- rows
    artifacts[[name]] <<- descriptor
  }
  add_artifact("landmarks", landmark_path, artifacts$landmarks$path, nrow(annotations))
  artifacts$landmarks$character_columns <- names(annotations)[
    vapply(annotations, is.character, logical(1)) & !names(annotations) %in% c("start", "end")]
  add_artifact("domain_annotation_provenance", provenance_path, "domain_annotation_run_manifest.json")
  add_artifact("code_domain_snapshot_builder", script_path, "code_domain_snapshot_builder.R")
  add_artifact("code_domain_snapshot_reader", file.path(root, "R/fungal_versions.R"), "code_domain_snapshot_reader.R")
  pipeline_path <- file.path(root, "scripts/build_fungal_reference_panel.py")
  if (file.exists(pipeline_path)) add_artifact("code_domain_retention_pipeline", pipeline_path, "code_domain_retention_pipeline.py")
  source_descriptors <- list()
  for (i in seq_along(sources)) {
    source <- sources[[i]]
    bytes <- readBin(source$path, "raw", n = file.info(source$path)$size)
    known <- grepl("[.](csv|json|fasta|fa|fna|R|py|txt)$", source$path)
    binary <- any(bytes == as.raw(0))
    extension <- if (known && !binary) tools::file_ext(source$path) else "txt"
    filename <- paste0("annotation_sources/source_", sprintf("%03d", i), ".", extension)
    artifact_name <- paste0("domain_annotation_source_", sprintf("%03d", i))
    if (binary) {
      encoded <- tempfile(tmpdir = staging)
      writeLines(jsonlite::base64_enc(bytes), encoded, useBytes = TRUE)
      add_artifact(artifact_name, encoded, filename)
      unlink(encoded)
    } else add_artifact(artifact_name, source$path, filename)
    source_descriptors[[artifact_name]] <- list(original_name = source$original_name,
      original_sha256 = source$sha256, encoding = if (binary) "base64" else "original_bytes")
  }
  manifest <- base$manifest
  manifest$transport_source <- NULL
  manifest$id <- id
  manifest$version_number <- number
  manifest$label <- paste0("V", number, " · ", options$label)
  manifest$description <- if (nzchar(options$description)) options$description else paste0(
    "Sequences, alignment, catalog and primer results are unchanged from ", options$base_id,
    "; reviewed LSU domain annotations and their source evidence are frozen in this version.")
  manifest$created_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  manifest$input_sha256 <- input_hash
  manifest$annotation_parent_id <- options$base_id
  manifest$annotation_parent_manifest_sha256 <- base_hash
  manifest$annotation_source_identity_sha256 <- source_identity
  manifest$annotation_provenance_sha256 <- artifacts$domain_annotation_provenance$sha256
  manifest$annotation_builder_sha256 <- artifacts$code_domain_snapshot_builder$sha256
  manifest$annotation_reader_sha256 <- artifacts$code_domain_snapshot_reader$sha256
  if (!is.null(artifacts$code_domain_retention_pipeline)) manifest$annotation_retention_pipeline_sha256 <- artifacts$code_domain_retention_pipeline$sha256
  manifest$annotation_base_commit <- tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[1], error = function(e) "unknown")
  manifest$annotation_working_tree_modified <- length(system2("git", c("status", "--porcelain"), stdout = TRUE, stderr = FALSE)) > 0L
  manifest$annotation_source_archives <- source_descriptors
  manifest$annotation_policy <- "Per-accession raw coordinates; source evidence distinguishes structural-reference projection from direct gene annotation. Existing primer results are not recomputed."
  manifest$artifacts <- artifacts
  jsonlite::write_json(manifest, file.path(staging, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
  for (name in names(artifacts)) fungal_snapshot_artifact(manifest, name, staging)
  if (!identical(atlas_sha256(index_path), index_hash)) stop("Version index changed during annotation build; retry from the updated history")
  if (!file.rename(staging, destination)) stop("Could not finalize annotation snapshot")
  entry <- data.frame(id = id, label = manifest$label, description = manifest$description,
    created_at = manifest$created_at, manifest = file.path(id, "manifest.json"), stringsAsFactors = FALSE)
  updated <- rbind(index, entry)
  temporary_index <- tempfile(tmpdir = version_root)
  jsonlite::write_json(list(schema_version = "1.0.0", marker_id = "ITS_FUNGAL",
    versions = lapply(seq_len(nrow(updated)), function(i) as.list(updated[i, , drop = FALSE]))),
    temporary_index, auto_unbox = TRUE, pretty = TRUE)
  if (!file.rename(temporary_index, index_path)) stop("Could not append annotation snapshot index")
  published <- TRUE
}, error = function(e) { cleanup(); stop(conditionMessage(e), call. = FALSE) })
cleanup()
message("Created immutable annotation-only snapshot ", id, "; source ", options$base_id,
  "; unchanged ", length(base$sequences), " sequences and primer results; ", nrow(annotations), " landmark rows.")
