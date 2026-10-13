#!/usr/bin/env Rscript

project_root <- normalizePath(getwd(), mustWork = TRUE)
local_library <- file.path(project_root, ".Rlib")
.libPaths(c(local_library, .libPaths()))

suppressPackageStartupMessages({
  library(PrimerMiner)
  library(rentrez)
})

source(file.path(project_root, "R", "functions.R"))
source(file.path(project_root, "R", "data_layer.R"))
source(file.path(project_root, "R", "coi_retention.R"))

orders_path <- file.path(project_root, "data", "orders_22.csv")
orders <- read.csv(orders_path, stringsAsFactors = FALSE, check.names = FALSE)
orders <- orders[orders$order != "Bivalvia", , drop = FALSE]
requested_raw <- Sys.getenv("ORDERS", "")
if (nzchar(requested_raw)) {
  requested <- trimws(strsplit(requested_raw, ",", fixed = TRUE)[[1]])
  missing_orders <- setdiff(requested, orders$order)
  if (length(missing_orders)) {
    stop("Unknown ORDERS value(s): ", paste(missing_orders, collapse = ", "))
  }
  orders <- orders[match(requested, orders$order), , drop = FALSE]
}

evaluation_mode <- tolower(Sys.getenv("EVALUATION_MODE", "full"))
if (!evaluation_mode %in% c("full", "centroid_preview")) {
  stop("EVALUATION_MODE must be 'full' or 'centroid_preview'.")
}
force <- identical(tolower(Sys.getenv("FORCE", "false")), "true")
set.seed(20260728)

vsearch <- Sys.which("vsearch")
mafft <- Sys.which("mafft")
if (!nzchar(vsearch)) stop("VSEARCH was not found on PATH.")
if (!nzchar(mafft)) stop("MAFFT was not found on PATH.")

reference_dir <- file.path(project_root, "data", "reference")
raw_root <- file.path(project_root, "data", "primerminer", "raw")
cluster_root <- file.path(project_root, "data", "primerminer", "clustered")
alignment_root <- file.path(project_root, "data", "primerminer", "aligned_reference")
full_alignment_root <- file.path(project_root, "data", "primerminer", "full_aligned_reference")
provenance_dir <- file.path(project_root, "data", "provenance")
dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_root, recursive = TRUE, showWarnings = FALSE)
dir.create(cluster_root, recursive = TRUE, showWarnings = FALSE)
dir.create(alignment_root, recursive = TRUE, showWarnings = FALSE)
dir.create(full_alignment_root, recursive = TRUE, showWarnings = FALSE)
dir.create(provenance_dir, recursive = TRUE, showWarnings = FALSE)

reference_accession <- "NC_001322.1"
reference_start <- 1474L
reference_stop <- 3009L
reference_path <- file.path(reference_dir, "coi_reference_NC_001322.1.fasta")

write_fasta_text <- function(text, path) {
  text <- sub("^>[^\\r\\n]*", paste0(
    ">NC_001322.1_COX1_reference coordinates=", reference_start, "-", reference_stop
  ), text)
  writeLines(text, path, useBytes = TRUE)
}

if (!file.exists(reference_path)) stop("Install the frozen versioned COX1 coordinate reference before aligning; this builder does not acquire sequences")
reference_sequences <- parse_fasta(reference_path)
if (length(reference_sequences) != 1L || nchar(reference_sequences[[1]]) != 1536L) {
  stop("Reference FASTA must contain exactly one 1,536-bp sequence.")
}

count_fasta <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  sum(startsWith(readLines(path, warn = FALSE), ">"))
}

run_mafft <- function(reference, fragments, output, log_path) {
  temporary <- tempfile(pattern = "coi_alignment_", tmpdir = dirname(output))
  on.exit(unlink(temporary), add = TRUE)
  args <- c(
    "--quiet",
    "--adjustdirectionaccurately",
    "--keeplength",
    "--addfragments", shQuote(fragments),
    shQuote(reference)
  )
  status <- system2(mafft, args, stdout = temporary, stderr = log_path)
  if (!identical(status, 0L)) stop("MAFFT failed; see ", log_path)

  aligned <- parse_fasta(temporary)
  widths <- nchar(aligned)
  if (!length(aligned) || any(widths != 1536L)) {
    stop("Reference-coordinate alignment is not exactly 1,536 columns: ", output)
  }
  if (!file.rename(temporary, output)) stop("Could not install completed COI alignment")
}

run_primerminer_clustering <- function(raw_path, destination_dir) {
  scratch <- tempfile(pattern = "primerminer_cluster_")
  dir.create(scratch, recursive = TRUE)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  scratch_input <- file.path(scratch, basename(raw_path))
  file.copy(raw_path, scratch_input, overwrite = TRUE)

  PrimerMiner::Clustering(
    file = basename(scratch_input),
    vsearchpath = vsearch,
    id = 0.97,
    threshold = "Majority",
    setwd = scratch
  )

  generated <- file.path(
    scratch,
    paste0(sub("\\.fasta$", "", basename(raw_path)), "_cons_cluster_Majority.fasta")
  )
  if (!file.exists(generated)) {
    stop("PrimerMiner clustering did not create ", basename(generated))
  }
  file.copy(generated, destination_dir, overwrite = TRUE)
  if (file.exists(file.path(scratch, "log.txt"))) {
    file.copy(file.path(scratch, "log.txt"), destination_dir, overwrite = TRUE)
  }
  if (dir.exists(file.path(scratch, "Vsearch"))) {
    unlink(file.path(destination_dir, "Vsearch"), recursive = TRUE)
    dir.create(file.path(destination_dir, "Vsearch"), showWarnings = FALSE)
    file.copy(
      list.files(file.path(scratch, "Vsearch"), full.names = TRUE),
      file.path(destination_dir, "Vsearch"),
      recursive = TRUE,
      overwrite = TRUE
    )
  }
}

manifest_rows <- vector("list", nrow(orders))
aligner_version <- system2(mafft, "--version", stdout = TRUE, stderr = TRUE)

# Source labels describe the frozen acquisition input. This stage only aligns
# those records; it never submits a new query or applies an acquisition subset.
coi_alignment_acquisition_provenance <- function(root, group) {
  relative <- "data/derived/coi_reference_panel/run_manifest.json"
  path <- file.path(root, relative)
  result <- list(query = "", candidate_limit = NA_integer_,
    acquisition_run_manifest = "", acquisition_run_manifest_sha256 = "",
    acquisition_note = "Alignment stage acquires no sequences; no acquisition manifest is available for this legacy input.")
  if (!file.exists(path)) return(result)
  manifest <- jsonlite::read_json(path, simplifyVector = FALSE)
  if (!identical(manifest$marker_id, "COI")) stop("Alignment acquisition manifest is not a COI panel")
  matches <- Filter(function(row) identical(row$target_group, group), manifest$searches)
  if (length(matches) != 1L) stop("Acquisition provenance needs exactly one current search for ", group)
  row <- matches[[1]]
  if (!is.character(row$query) || length(row$query) != 1L || !nzchar(row$query) ||
      !is.numeric(row$fetched_candidate_limit) || length(row$fetched_candidate_limit) != 1L ||
      is.na(row$fetched_candidate_limit) || row$fetched_candidate_limit < 1 ||
      row$fetched_candidate_limit != as.integer(row$fetched_candidate_limit)) {
    stop("Acquisition search provenance is incomplete for ", group)
  }
  result$query <- row$query
  result$candidate_limit <- as.integer(row$fetched_candidate_limit)
  result$acquisition_run_manifest <- relative
  result$acquisition_run_manifest_sha256 <- atlas_sha256(path)
  result$acquisition_note <- paste("Alignment stage acquires no sequences. Query and candidate_limit describe the current acquisition search before QC;",
    "retained earlier records can have different source searches recorded in their provenance.")
  result
}

for (i in seq_len(nrow(orders))) {
  order_name <- orders$order[i]
  message("\n=== ", order_name, " ===")
  raw_dir <- file.path(raw_root, order_name)
  cluster_dir <- file.path(cluster_root, order_name)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)

  raw_path <- file.path(raw_dir, paste0(order_name, "_GB.fasta"))
  cluster_path <- file.path(
    cluster_dir, paste0(order_name, "_GB_cons_cluster_Majority.fasta")
  )
  aligned_path <- if (evaluation_mode == "full") {
    file.path(full_alignment_root, paste0(order_name, "_COX1_full_reference_aligned.fasta"))
  } else {
    file.path(alignment_root, paste0(order_name, "_COX1_reference_aligned.fasta"))
  }
  mafft_log <- file.path(dirname(aligned_path), paste0(order_name, "_mafft.log"))
  input_manifest_path <- paste0(aligned_path, ".inputs.json")
  alignment_inputs <- NULL

  status <- "complete"
  note <- ""
  tryCatch({
    if (!file.exists(raw_path)) stop("Retained raw input is missing; run the feature-aware append-only COI acquisition builder first")
    if (!file.exists(raw_path) || count_fasta(raw_path) < 1L) {
      stop("No GenBank sequences were downloaded.")
    }

    cluster_input_path <- paste0(cluster_path, ".raw.sha256")
    current_raw_hash <- atlas_sha256(raw_path)
    cluster_current <- file.exists(cluster_input_path) && identical(readLines(cluster_input_path, warn = FALSE), current_raw_hash)
    if (evaluation_mode == "centroid_preview" && (!file.exists(cluster_path) || force || !cluster_current)) {
      unlink(file.path(cluster_dir, "Vsearch"), recursive = TRUE)
      unlink(file.path(cluster_dir, "log.txt"))
      run_primerminer_clustering(raw_path, cluster_dir)
      generated <- file.path(
        cluster_dir,
        paste0(sub("\\.fasta$", "", basename(raw_path)), "_cons_cluster_Majority.fasta")
      )
      if (!identical(generated, cluster_path) && file.exists(generated)) {
        file.rename(generated, cluster_path)
      }
      writeLines(current_raw_hash, cluster_input_path)
    }
    if (evaluation_mode == "centroid_preview" && (!file.exists(cluster_path) || count_fasta(cluster_path) < 1L)) {
      stop("PrimerMiner clustering did not produce consensus sequences.")
    }

    alignment_input <- if (evaluation_mode == "full") raw_path else cluster_path
    alignment_inputs <- coi_alignment_inputs(alignment_input, reference_path, aligner_version,
      evaluation_mode, file.path(project_root, "scripts/build_22_group_alignments.R"))
    cached_inputs <- if (file.exists(input_manifest_path)) tryCatch(jsonlite::read_json(input_manifest_path), error = function(e) NULL) else NULL
    cache_valid <- !force && file.exists(aligned_path) && count_fasta(aligned_path) > 0L &&
      !is.null(cached_inputs) && identical(cached_inputs$input_sha256, alignment_inputs$input_sha256) &&
      identical(cached_inputs$output_sha256, atlas_sha256(aligned_path))
    if (cache_valid && evaluation_mode == "full") {
      cache_valid <- isTRUE(tryCatch({coi_validate_full_alignment(aligned_path, raw_path, reference_path, cached_inputs); TRUE}, error = function(e) FALSE))
    }
    if (!cache_valid) {
      run_mafft(reference_path, alignment_input, aligned_path, mafft_log)
      if (evaluation_mode == "full") coi_validate_full_alignment(aligned_path, raw_path, reference_path, alignment_inputs)
      alignment_inputs$output_sha256 <- atlas_sha256(aligned_path)
      alignment_inputs$built_at_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
      alignment_inputs$raw_sequences <- count_fasta(alignment_input)
      alignment_inputs$aligned_templates <- count_fasta(aligned_path) - 1L
      jsonlite::write_json(alignment_inputs, input_manifest_path, auto_unbox = TRUE, pretty = TRUE)
    }
  }, error = function(e) {
    status <<- "failed"
    note <<- conditionMessage(e)
    message("FAILED: ", note)
  })

  acquisition <- coi_alignment_acquisition_provenance(project_root, order_name)
  manifest_rows[[i]] <- data.frame(
    order = order_name,
    environment_group = orders$environment_group[i],
    query = acquisition$query,
    requested_subset = NA_integer_,
    candidate_limit = acquisition$candidate_limit,
    acquisition_run_manifest = acquisition$acquisition_run_manifest,
    acquisition_run_manifest_sha256 = acquisition$acquisition_run_manifest_sha256,
    acquisition_note = acquisition$acquisition_note,
    downloaded_sequences = count_fasta(raw_path),
    consensus_clusters = count_fasta(cluster_path),
    aligned_sequences_including_reference = count_fasta(aligned_path),
    reference_accession = reference_accession,
    reference_feature_coordinates = paste0(reference_start, "-", reference_stop),
    reference_length_bp = 1536L,
    clustering_identity = 0.97,
    evaluation_mode = evaluation_mode,
    scoring_sequences = if (evaluation_mode == "full") count_fasta(raw_path) else count_fasta(cluster_path),
    raw_sha256 = if (file.exists(raw_path)) atlas_sha256(raw_path) else NA_character_,
    reference_sha256 = atlas_sha256(reference_path),
    alignment_input_sha256 = if (!is.null(alignment_inputs)) alignment_inputs$input_sha256 else NA_character_,
    alignment_sha256 = if (file.exists(aligned_path)) atlas_sha256(aligned_path) else NA_character_,
    alignment_input_manifest = sub(paste0("^", project_root, "/"), "", input_manifest_path),
    status = status,
    note = note,
    raw_fasta = sub(paste0("^", project_root, "/"), "", raw_path),
    cluster_consensus_fasta = sub(paste0("^", project_root, "/"), "", cluster_path),
    reference_aligned_fasta = sub(paste0("^", project_root, "/"), "", aligned_path),
    built_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
}

manifest_path <- file.path(provenance_dir, "order_alignment_manifest.csv")
new_manifest <- do.call(rbind, manifest_rows)
if (file.exists(manifest_path) && nzchar(requested_raw)) {
  old_manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  old_manifest <- old_manifest[!old_manifest$order %in% new_manifest$order, , drop = FALSE]
  all_columns <- union(names(old_manifest), names(new_manifest))
  for (column in setdiff(all_columns, names(old_manifest))) old_manifest[[column]] <- NA
  for (column in setdiff(all_columns, names(new_manifest))) new_manifest[[column]] <- NA
  old_manifest <- old_manifest[, all_columns, drop = FALSE]
  new_manifest <- new_manifest[, all_columns, drop = FALSE]
  new_manifest <- rbind(old_manifest, new_manifest)
}
new_manifest <- new_manifest[order(new_manifest$order), , drop = FALSE]
write.csv(new_manifest, manifest_path, row.names = FALSE, na = "")
message("\nWrote ", manifest_path)
print(new_manifest[, c(
  "order", "downloaded_sequences", "consensus_clusters",
  "aligned_sequences_including_reference", "status", "note"
)], row.names = FALSE)
if (any(manifest_rows |> vapply(function(x) x$status != "complete", logical(1)))) stop("One or more retained COI groups could not be aligned; partial builds cannot be released")
