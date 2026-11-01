#!/usr/bin/env Rscript

# Capture V1 before updating source files, then build V2 from a reviewed panel.
# Content-derived identities and refusal to overwrite keep each release stable.
root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(root, ".Rlib"), .libPaths()))
source("R/functions.R")
source("R/data_layer.R")
source("R/fungal_versions.R")
if (file.exists("R/fungal_explorer.R")) source("R/fungal_explorer.R")
version_root <- fungal_version_root(root)
dir.create(version_root, recursive = TRUE, showWarnings = FALSE)
baseline_root <- file.path(version_root, ".baseline-inputs")

capture_fungal_baseline <- function() {
  dir.create(baseline_root, recursive = TRUE, showWarnings = FALSE)
  paths <- c(
    raw_catalog.csv = "data/catalog/unite_primers.csv",
    panel.fasta = "data/reference/its_fungal_reference_FN812768.2.fasta",
    geometry.csv = "data/derived/marker_pair_geometry.csv",
    landmarks.csv = "data/catalog/marker_landmarks.csv",
    primer_pairs.csv = "data/catalog/primer_pairs.csv",
    oligos.csv = "data/catalog/oligos.csv",
    pair_oligos.csv = "data/catalog/pair_oligos.csv"
  )
  existing_index <- fungal_read_index(root, remote = FALSE)
  baseline_ids <- existing_index$id[startsWith(existing_index$id, "v1-")]
  if (length(baseline_ids)) {
    baseline <- fungal_read_snapshot(baseline_ids[1], root, remote = FALSE)
    base_manifest_path <- file.path(version_root, existing_index$manifest[match(baseline_ids[1], existing_index$id)])
    snapshot_names <- c("raw_catalog", "sequences", "geometry", "landmarks", "primer_pairs", "oligos", "pair_oligos")
    paths <- setNames(vapply(snapshot_names, function(name) {
      file.path(dirname(base_manifest_path), baseline$manifest$artifacts[[name]]$path)
    }, character(1)), names(paths))
  }
  for (name in names(paths)) {
    target <- file.path(baseline_root, name)
    if (!file.exists(target) && !file.copy(paths[[name]], target)) stop("Could not capture V1 input: ", name)
  }
}

snapshot_fasta <- function(path) atlas_read_artifact(list(local_path = path))
write_snapshot_fasta <- function(sequences, path) {
  lines <- unlist(lapply(seq_along(sequences), function(i) {
    c(paste0(">", names(sequences)[i]), unname(sequences[[i]]))
  }), use.names = FALSE)
  writeLines(lines, path, useBytes = TRUE)
}
read_snapshot_csv <- function(path) read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
write_snapshot_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "")

build_fungal_snapshot <- function(number, panel_dir = Sys.getenv("FUNGAL_PANEL_DIR", "data/derived/fungal_reference_panel")) {
  if (number == 1L) {
    existing <- fungal_read_index(root, remote = FALSE)
    ids <- existing$id[startsWith(existing$id, "v1-")]
    if (length(ids)) {
      fungal_read_snapshot(ids[1], root, remote = FALSE)
      message("Retained original immutable baseline: ", ids[1])
      return(invisible(ids[1]))
    }
  }
  if (!exists("fungal_evaluate_panel", mode = "function")) {
    stop("R/fungal_explorer.R must exist to freeze reproducible explorer outputs; V1 source inputs are already captured")
  }
  frozen_catalog <- file.path(baseline_root, "raw_catalog.csv")
  catalog_source <- if (number == 1L) frozen_catalog else "data/catalog/unite_primers.csv"
  catalog <- fungal_catalog(catalog_source)
  raw_catalog <- read_snapshot_csv(catalog_source)
  geometry <- read_snapshot_csv(file.path(baseline_root, "geometry.csv"))
  geometry <- geometry[geometry$marker_id == "ITS_FUNGAL", , drop = FALSE]
  if (number == 1L) {
    sequence_path <- file.path(baseline_root, "panel.fasta")
    alignment_path <- sequence_path
    sequences <- snapshot_fasta(sequence_path)
    taxonomy <- data.frame(
      accession = sub(" .*", "", names(sequences)),
      organism = "Uncultured fungus", scientific_name = "Uncultured fungus",
      phylum = NA_character_, class = NA_character_, order = NA_character_,
      family = NA_character_, genus = NA_character_, taxid = NA_character_,
      taxonomy_status = "unresolved_baseline", source = "NCBI FN812768.2 FASTA header",
      stringsAsFactors = FALSE
    )
    landmarks <- read_snapshot_csv(file.path(baseline_root, "landmarks.csv"))
    landmarks <- landmarks[landmarks$marker_id == "ITS_FUNGAL", , drop = FALSE]
    description <- "Original FN812768.2 reference, captured UNITE catalog, eight mapped pairs and pilot landmarks. Taxonomic class and clade are unresolved."
    run_manifest <- list(source = "NCBI nuccore FN812768.2", reference_accession = "FN812768.2", scope = "one reference")
  } else {
    required <- file.path(panel_dir, c("panel.fasta", "panel_manifest.csv", "alignment.fasta", "landmarks.csv", "run_manifest.json"))
    if (any(!file.exists(required))) stop("Missing reviewed V2 inputs: ", paste(basename(required[!file.exists(required)]), collapse = ", "))
    sequence_path <- required[1]
    alignment_path <- required[3]
    sequences <- snapshot_fasta(sequence_path)
    taxonomy <- read_snapshot_csv(required[2])
    landmarks <- read_snapshot_csv(required[4])
    run_manifest <- jsonlite::read_json(required[5], simplifyVector = FALSE)
    for (filename in names(run_manifest$artifacts_sha256)) {
      path <- file.path(panel_dir, filename)
      if (!file.exists(path) || !identical(atlas_sha256(path), run_manifest$artifacts_sha256[[filename]])) {
        stop("Fungal panel differs from its acquisition manifest: ", filename)
      }
    }
    description <- "Bounded NCBI fungal rDNA reference panel with frozen accession taxonomy, alignment, landmark evidence and exact primer-site results."
  }
  names(sequences) <- sub(" .*", "", names(sequences))
  if (anyDuplicated(names(sequences))) stop("Fungal panel repeats an accession")
  if (anyDuplicated(taxonomy$accession) || !setequal(names(sequences), taxonomy$accession)) stop("Fungal sequence and taxonomy identities differ")
  if ("sequence_sha256" %in% names(taxonomy)) {
    sequence_hashes <- vapply(unname(sequences), digest::digest, character(1), algo = "sha256", serialize = FALSE)
    expected <- taxonomy$sequence_sha256[match(names(sequences), taxonomy$accession)]
    if (!identical(unname(sequence_hashes), unname(expected))) stop("Fungal accession sequence hash differs from taxonomy ledger")
  }
  alignment <- snapshot_fasta(alignment_path)
  aligned_accessions <- sub("^_R_", "", sub(" .*", "", names(alignment)))
  if (anyDuplicated(aligned_accessions) || !setequal(aligned_accessions, names(sequences)) || length(unique(nchar(alignment))) != 1L) {
    stop("Fungal alignment does not retain the complete sequence panel")
  }
  evaluated <- fungal_evaluate_panel(catalog, sequences, taxonomy)
  staging <- tempfile(pattern = ".snapshot-", tmpdir = version_root)
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE), add = TRUE)
  frames <- c(list(raw_catalog = raw_catalog, catalog = catalog, taxonomy = taxonomy,
                   geometry = geometry, landmarks = landmarks),
              evaluated[c("pairs", "pair_sequences", "sites")])
  # Retain the original documented-pair records as provenance, independent
  # of whether today's selected references span both binding sites.
  for (name in c("primer_pairs", "oligos", "pair_oligos")) {
    x <- read_snapshot_csv(file.path(baseline_root, paste0(name, ".csv")))
    if ("marker_id" %in% names(x)) x <- x[x$marker_id == "ITS_FUNGAL", , drop = FALSE]
    frames[[name]] <- x
  }
  artifacts <- list()
  for (name in names(frames)) {
    filename <- paste0(name, ".csv")
    path <- file.path(staging, filename)
    write_snapshot_csv(frames[[name]], path)
    artifacts[[name]] <- list(path = filename, sha256 = atlas_sha256(path), bytes = unname(file.info(path)$size), rows = nrow(frames[[name]]))
  }
  for (name in c("sequences", "alignment")) {
    path <- file.path(staging, paste0(name, ".fasta"))
    value <- if (name == "sequences") sequences else alignment
    write_snapshot_fasta(value, path)
    artifacts[[name]] <- list(path = basename(path), sha256 = atlas_sha256(path), bytes = unname(file.info(path)$size), rows = length(value))
  }
  jsonlite::write_json(run_manifest, file.path(staging, "run_manifest.json"), auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
  artifacts$run_manifest <- list(path = "run_manifest.json", sha256 = atlas_sha256(file.path(staging, "run_manifest.json")), bytes = unname(file.info(file.path(staging, "run_manifest.json"))$size))
  code_paths <- c(code_functions = "R/functions.R", code_explorer = "R/fungal_explorer.R",
                  code_acquisition = "scripts/build_fungal_reference_panel.py",
                  code_snapshot_builder = "scripts/build_fungal_snapshots.R",
                  code_version_reader = "R/fungal_versions.R", code_data_layer = "R/data_layer.R")
  for (name in names(code_paths)) {
    source_path <- code_paths[[name]]
    if (!file.exists(source_path)) next
    filename <- paste0(name, ".", tools::file_ext(source_path))
    path <- file.path(staging, filename)
    if (!file.copy(source_path, path)) stop("Could not freeze scientific code")
    artifacts[[name]] <- list(path = filename, sha256 = atlas_sha256(path), bytes = unname(file.info(path)$size), rows = length(readLines(path, warn = FALSE)))
  }
  scientific_code <- c("code_functions", "code_explorer", "code_acquisition")
  tooling_code <- c("code_snapshot_builder", "code_version_reader", "code_data_layer")
  code_sha256 <- digest::digest(paste(vapply(artifacts[intersect(scientific_code, names(artifacts))], `[[`, character(1), "sha256"), collapse = "|"), algo = "sha256", serialize = FALSE)
  tooling_sha256 <- digest::digest(paste(vapply(artifacts[intersect(tooling_code, names(artifacts))], `[[`, character(1), "sha256"), collapse = "|"), algo = "sha256", serialize = FALSE)
  # Acquisition time and retained/new counters document the run, but cannot
  # turn an unchanged biological panel into a new version on offline replay.
  stable_source <- run_manifest[intersect(c("source", "sampling_policy", "searches", "raw_responses", "taxonomy_snapshot_sha256", "mafft_version", "alignment_policy"), names(run_manifest))]
  if (!is.null(stable_source$searches)) {
    stable_source$searches <- lapply(stable_source$searches, function(search) {
      search[intersect(c("clade", "query", "count_available", "selected", "sha256"), names(search))]
    })
  }
  source_identity <- digest::digest(jsonlite::toJSON(stable_source, auto_unbox = TRUE, null = "null"), algo = "sha256", serialize = FALSE)
  identity_artifacts <- artifacts[setdiff(names(artifacts), c("run_manifest", tooling_code))]
  input_sha256 <- digest::digest(paste(c(vapply(identity_artifacts, `[[`, character(1), "sha256"), code_sha256, source_identity), collapse = "|"), algo = "sha256", serialize = FALSE)
  index <- fungal_read_index(root, remote = FALSE)
  existing_manifests <- lapply(index$manifest, function(relative) {
    jsonlite::read_json(file.path(version_root, relative), simplifyVector = FALSE)
  })
  matching <- which(vapply(existing_manifests, function(existing) {
    identical(existing$input_sha256, input_sha256)
  }, logical(1)))
  if (length(matching)) {
    id <- existing_manifests[[matching[1]]]$id
    fungal_read_snapshot(id, root, remote = FALSE)
    message("Fungal snapshot content already exists: ", id)
    return(invisible(id))
  }
  if (number > 1L) {
    previous_numbers <- vapply(existing_manifests, function(existing) {
      as.integer(existing$version_number)
    }, integer(1))
    number <- if (length(previous_numbers)) max(previous_numbers) + 1L else 2L
  }
  id <- paste0("v", number, "-", substr(input_sha256, 1L, 12L))
  label <- if (number == 1L) "V1 · original fungal reference" else paste0("V", number, " · NCBI fungal reference panel")
  destination <- file.path(version_root, id)
  if (dir.exists(destination)) {
    existing <- jsonlite::read_json(file.path(destination, "manifest.json"), simplifyVector = FALSE)
    if (!identical(existing$input_sha256, input_sha256)) stop("Immutable fungal snapshot identity collision")
    for (name in names(existing$artifacts)) fungal_snapshot_artifact(existing, name, destination)
    message("Fungal snapshot already exists: ", id)
    return(invisible(id))
  }
  commit <- tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[1], error = function(e) "unknown")
  manifest <- list(
    schema_version = "1.0.0", marker_id = "ITS_FUNGAL", id = id,
    version_number = number, label = label, description = description,
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    sequence_count = length(sequences), reference_accession = if (number == 1L) "FN812768.2" else run_manifest$reference_accession,
    input_sha256 = input_sha256, code_sha256 = code_sha256,
    tooling_sha256 = tooling_sha256,
    base_commit = commit, code_commit = commit, code_commit_role = "base checkout; frozen code artifacts are authoritative",
    working_tree_modified = length(system2("git", c("status", "--porcelain"), stdout = TRUE, stderr = FALSE)) > 0L,
    source_identity_sha256 = source_identity,
    engine = "exact_iupac_site_search", source_manifest = run_manifest,
    artifacts = artifacts
  )
  jsonlite::write_json(manifest, file.path(staging, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
  if (!file.rename(staging, destination)) stop("Could not finalize fungal snapshot")
  entry <- data.frame(id = id, label = label, description = description, created_at = manifest$created_at,
                      manifest = file.path(id, "manifest.json"), stringsAsFactors = FALSE)
  index <- rbind(index, entry)
  entries <- lapply(seq_len(nrow(index)), function(i) as.list(index[i, , drop = FALSE]))
  temporary <- tempfile(tmpdir = version_root)
  jsonlite::write_json(list(schema_version = "1.0.0", marker_id = "ITS_FUNGAL", versions = entries), temporary, auto_unbox = TRUE, pretty = TRUE)
  if (!file.rename(temporary, file.path(version_root, "index.json"))) stop("Could not update fungal snapshot index")
  message("Created immutable fungal snapshot ", id, ": ", length(sequences), " sequences, ", nrow(evaluated$pairs), " candidate pairs")
  invisible(id)
}

requested <- commandArgs(trailingOnly = TRUE)
if (!length(requested)) requested <- "all"
if ("restore" %in% requested) {
  previous <- atlas_read_manifest("ITS_FUNGAL", root)
  if (isTRUE(previous$stale) && nzchar(Sys.getenv("ATLAS_MANIFEST_BASE_URL", ""))) {
    stop("Previous fungal release is unreachable; refusing to reset retained sequences")
  }
  descriptors <- previous$manifest$fungal_snapshots
  if (!length(descriptors)) {
    message("Previous release is catalog-only; no fungal sequence panel to restore")
  } else {
    # A clean git runner must restore every earlier immutable snapshot, not
    # only the latest sequence panel, before allocating the next ordinal.
    local_index <- fungal_read_index(root, remote = FALSE)
    for (historical in descriptors) {
      if (!grepl("^v[0-9]+-[a-f0-9]{12}$", historical$id)) stop("Invalid prior fungal snapshot identity")
      remote_manifest <- atlas_cached_artifact(historical$manifest$url, historical$manifest$sha256)
      frozen <- jsonlite::read_json(remote_manifest, simplifyVector = FALSE)
      if (!identical(frozen$id, historical$id)) stop("Prior fungal manifest identity differs")
      destination_snapshot <- file.path(version_root, historical$id)
      if (dir.exists(destination_snapshot)) {
        if (!identical(atlas_sha256(file.path(destination_snapshot, "manifest.json")), historical$manifest$sha256)) stop("Prior immutable fungal identity changed")
        for (name in names(frozen$artifacts)) fungal_snapshot_artifact(frozen, name, destination_snapshot)
      } else {
        temporary_snapshot <- tempfile(pattern = ".restore-", tmpdir = version_root)
        dir.create(temporary_snapshot)
        for (name in names(frozen$artifacts)) {
          artifact <- historical$artifacts[[name]]
          inner <- frozen$artifacts[[name]]
          if (is.null(artifact) || !identical(artifact$sha256, inner$sha256) ||
              grepl("(^/|(^|/)\\.\\.(/|$))", inner$path)) stop("Prior fungal artifact identity differs")
          source_path <- atlas_cached_artifact(artifact$url, artifact$sha256)
          target <- file.path(temporary_snapshot, inner$path)
          dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
          if (!file.copy(source_path, target)) stop("Could not restore prior fungal snapshot artifact")
        }
        if (!file.copy(remote_manifest, file.path(temporary_snapshot, "manifest.json")) ||
            !file.rename(temporary_snapshot, destination_snapshot)) stop("Could not restore immutable fungal snapshot")
      }
      if (!historical$id %in% local_index$id) {
        local_index <- rbind(local_index, data.frame(id = historical$id, label = frozen$label,
          description = frozen$description, created_at = frozen$created_at,
          manifest = file.path(historical$id, "manifest.json"), stringsAsFactors = FALSE))
      }
    }
    ordinals <- vapply(local_index$manifest, function(relative) {
      as.integer(jsonlite::read_json(file.path(version_root, relative))$version_number)
    }, integer(1))
    local_index <- local_index[order(ordinals), , drop = FALSE]
    entries <- lapply(seq_len(nrow(local_index)), function(i) as.list(local_index[i, , drop = FALSE]))
    temporary_index <- tempfile(tmpdir = version_root)
    jsonlite::write_json(list(schema_version = "1.0.0", marker_id = "ITS_FUNGAL", versions = entries), temporary_index, auto_unbox = TRUE, pretty = TRUE)
    if (!file.rename(temporary_index, file.path(version_root, "index.json"))) stop("Could not restore fungal history index")
    descriptor <- descriptors[[length(descriptors)]]
    destination <- "data/derived/fungal_previous_reference_panel"
    dir.create(destination, recursive = TRUE, showWarnings = FALSE)
    mapping <- c(sequences = "panel.fasta", taxonomy = "panel_manifest.csv", landmarks = "landmarks.csv", run_manifest = "run_manifest.json")
    for (name in names(mapping)) {
      artifact <- descriptor$artifacts[[name]]
      if (is.null(artifact)) stop("Previous fungal panel lacks retained input: ", name)
      source_path <- artifact$local_path
      if (is.null(source_path) || !file.exists(source_path)) source_path <- atlas_cached_artifact(artifact$url, artifact$sha256)
      if (!identical(atlas_sha256(source_path), artifact$sha256)) stop("Previous fungal input checksum differs")
      if (!file.copy(source_path, file.path(destination, mapping[[name]]), overwrite = TRUE)) stop("Could not restore fungal input")
    }
    message("Restored complete fungal sequence, taxonomy and landmark state: ", descriptor$id)
  }
} else if ("capture" %in% requested) {
  capture_fungal_baseline()
  message("Captured original fungal inputs")
} else {
  capture_fungal_baseline()
  if (any(requested %in% c("all", "v1"))) build_fungal_snapshot(1L)
  if (any(requested %in% c("all", "v2"))) build_fungal_snapshot(2L)
  unlink(baseline_root, recursive = TRUE)
}
