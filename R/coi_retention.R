# Monthly COI state is part of an immutable release. A ledger alone is never
# sufficient evidence that its sequences survived a clean-worker rebuild.

coi_accession <- function(header) sub("[|].*$", "", sub("^_R_", "", sub("[[:space:]].*$", "", header)))

coi_read_fasta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  starts <- which(startsWith(lines, ">"))
  if (!length(starts)) stop("COI FASTA contains no records: ", basename(path))
  ends <- c(starts[-1] - 1L, length(lines))
  values <- vapply(seq_along(starts), function(i) {
    if (ends[i] <= starts[i]) return("")
    toupper(paste0(gsub("[[:space:]]", "", lines[seq.int(starts[i] + 1L, ends[i])]), collapse = ""))
  }, character(1))
  names(values) <- sub("^>", "", lines[starts])
  if (any(!nzchar(values)) || any(!grepl("^[ACGTURYSWKMBDHVNI]+$", values))) stop("Invalid retained COI sequence")
  accessions <- coi_accession(names(values))
  if (any(!nzchar(accessions)) || anyDuplicated(accessions)) stop("Duplicate or missing retained COI accession")
  values
}

coi_retained_records <- function(root = getwd(), groups = NULL) {
  if (is.null(groups)) {
    orders <- read.csv(file.path(root, "data/orders_22.csv"), stringsAsFactors = FALSE)
    groups <- orders$order[orders$order != "Bivalvia"]
  }
  blocks <- lapply(groups, function(group) {
    path <- file.path(root, "data/primerminer/raw", group, paste0(group, "_GB.fasta"))
    if (!file.exists(path)) stop("Retained raw COI FASTA is missing for ", group)
    seqs <- coi_read_fasta(path)
    data.frame(target_group = group, accession = coi_accession(names(seqs)),
      header = names(seqs), sequence_length = nchar(seqs),
      sequence_sha256 = vapply(unname(seqs), digest::digest, character(1), algo = "sha256", serialize = FALSE),
      header_sha256 = vapply(names(seqs), digest::digest, character(1), algo = "sha256", serialize = FALSE),
      stringsAsFactors = FALSE)
  })
  records <- do.call(rbind, blocks)
  if (is.null(records) || !nrow(records)) stop("No retained COI records are installed")
  if (anyDuplicated(paste(records$target_group, records$accession))) stop("Duplicate retained group/accession")
  records
}

coi_validate_retained_state <- function(root = getwd(), groups = NULL) {
  records <- coi_retained_records(root, groups)
  path <- file.path(root, "data/provenance/coi_reference_accession_ledger.csv")
  if (!file.exists(path)) stop("COI accession ledger is missing")
  ledger <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("target_group", "accession", "sequence_sha256", "header_sha256") %in% names(ledger))) stop("COI ledger lacks retained identities or byte hashes")
  keys <- paste(records$target_group, records$accession)
  retained <- ledger[ledger$target_group %in% unique(records$target_group), , drop = FALSE]
  ledger_keys <- paste(retained$target_group, retained$accession)
  if (anyDuplicated(ledger_keys) || !setequal(keys, ledger_keys)) stop("COI ledger and actual retained FASTA identities differ")
  for (column in c("sequence_sha256", "header_sha256")) {
    expected <- retained[[column]][match(keys, ledger_keys)]
    if (anyNA(expected) || any(!grepl("^[A-Fa-f0-9]{64}$", expected)) || !identical(tolower(expected), tolower(records[[column]]))) {
      stop("COI retained bytes differ from the accession ledger: ", column)
    }
  }
  list(records = records, ledger = ledger)
}

coi_assert_retention <- function(previous, current) {
  columns <- c("target_group", "accession", "sequence_sha256", "header_sha256")
  if (!all(columns %in% names(previous)) || !all(columns %in% names(current))) stop("Retained-record proof lacks identity or byte hashes")
  old_keys <- paste(previous$target_group, previous$accession)
  new_keys <- paste(current$target_group, current$accession)
  match_new <- match(old_keys, new_keys)
  if (anyNA(match_new)) stop("A previously retained COI accession was dropped")
  for (column in columns) {
    if (!identical(as.character(previous[[column]]), as.character(current[[column]][match_new]))) {
      stop("Previously retained COI record changed: ", column)
    }
  }
  invisible(TRUE)
}

coi_alignment_inputs <- function(raw_path, reference_path, aligner_version,
                                  evaluation_mode = "full", builder_path = NULL) {
  value <- list(raw_sha256 = atlas_sha256(raw_path), reference_sha256 = atlas_sha256(reference_path),
    evaluation_mode = evaluation_mode, reference_columns = 1536L,
    aligner = "MAFFT", aligner_version = paste(aligner_version, collapse = " "),
    parameters = c("--quiet", "--adjustdirectionaccurately", "--keeplength", "--addfragments"))
  if (!is.null(builder_path)) value$builder_sha256 <- atlas_sha256(builder_path)
  value$input_sha256 <- digest::digest(jsonlite::toJSON(value, auto_unbox = TRUE), algo = "sha256", serialize = FALSE)
  value
}

coi_validate_full_alignment <- function(alignment_path, raw_path, reference_path,
                                        input_manifest = NULL, minimum = 0L) {
  raw <- coi_read_fasta(raw_path)
  if (length(raw) < minimum) stop("Full COI panel is below its per-group minimum")
  aligned <- {
    lines <- readLines(alignment_path, warn = FALSE)
    starts <- which(startsWith(lines, ">")); ends <- c(starts[-1] - 1L, length(lines))
    if (!length(starts)) stop("Full COI alignment contains no rows")
    values <- vapply(seq_along(starts), function(i) {
      if (ends[i] <= starts[i]) return("")
      toupper(paste0(gsub("[[:space:]]", "", lines[seq.int(starts[i] + 1L, ends[i])]), collapse = ""))
    }, character(1))
    names(values) <- sub("^>", "", lines[starts]); values
  }
  reference <- grepl("^(_R_)?NC_001322[.]1_COX1_reference([[:space:]]|$)", names(aligned))
  if (sum(reference) != 1L || any(nchar(aligned) != 1536L)) stop("Full COI alignment must contain one reference and1536 columns")
  accessions <- coi_accession(names(aligned)[!reference])
  if (anyDuplicated(accessions) || !setequal(accessions, coi_accession(names(raw))) || length(accessions) != length(raw)) {
    stop("Full COI alignment does not contain every retained raw accession exactly once")
  }
  if (!identical(unname(aligned[reference]), unname(coi_read_fasta(reference_path)))) stop("Aligned COI coordinate reference changed")
  if (!is.null(input_manifest)) {
    if (!identical(input_manifest$raw_sha256, atlas_sha256(raw_path)) ||
        !identical(input_manifest$reference_sha256, atlas_sha256(reference_path))) stop("COI alignment input hash changed")
    if (!is.null(input_manifest$output_sha256) && !identical(input_manifest$output_sha256, atlas_sha256(alignment_path))) stop("COI alignment output checksum failed")
  }
  invisible(list(raw_sequences = length(raw), aligned_templates = length(accessions), columns = 1536L))
}

coi_stage_retained_state <- function(root, release_dir, release_id, data_base, groups = NULL) {
  state <- coi_validate_retained_state(root, groups)
  records_path <- file.path(root, "data/provenance/coi_reference_records.csv")
  write.csv(state$records, records_path, row.names = FALSE, na = "")
  sources <- c(
    records = records_path,
    ledger = file.path(root, "data/provenance/coi_reference_accession_ledger.csv"),
    growth_qa = file.path(root, "data/provenance/coi_reference_growth_qa.csv"),
    alignment_manifest = file.path(root, "data/provenance/order_alignment_manifest.csv"),
    scoring_manifest = file.path(root, "data/provenance/order_scoring_manifest.csv"),
    published_baseline = file.path(root, "data/provenance/coi_published_baseline.csv"),
    seed_restore = file.path(root, "data/provenance/coi_seed_restore.json")
  )
  panel_dir <- file.path(root, "data/derived/coi_reference_panel")
  # Freeze the acquisition metadata and the extracted/legacy panel as well as
  # the per-group input. Previous provenance is restored with sequence bytes.
  panel_names <- c("panel.fasta", "taxonomy.csv", "ledger.csv", "panel_manifest.csv", "growth_qa.csv", "rejections.csv", "run_manifest.json")
  panel_files <- file.path(panel_dir, panel_names)
  panel_files <- panel_files[file.exists(panel_files)]
  if (length(panel_files)) sources <- c(sources, setNames(panel_files, paste0("panel_", sub("[.].*$", "", basename(panel_files)))))
  archives <- unlist(lapply(c("raw", "legacy_raw"), function(folder) {
    list.files(file.path(panel_dir, folder), recursive = TRUE, full.names = TRUE)
  }), use.names = FALSE)
  archives <- archives[!dir.exists(archives)]
  if (length(archives)) sources <- c(sources, setNames(archives, paste0("source_", seq_along(archives))))
  for (group in unique(state$records$target_group)) {
    sources[[paste0("raw_", group)]] <- file.path(root, "data/primerminer/raw", group, paste0(group, "_GB.fasta"))
    alignment <- file.path(root, "data/primerminer/full_aligned_reference", paste0(group, "_COX1_full_reference_aligned.fasta"))
    if (file.exists(alignment)) sources[[paste0("alignment_", group)]] <- alignment
    provenance <- paste0(alignment, ".inputs.json")
    if (file.exists(provenance)) sources[[paste0("alignment_inputs_", group)]] <- provenance
  }
  sources <- sources[file.exists(sources)]
  descriptors <- list()
  for (name in names(sources)) {
    source <- sources[[name]]
    restore_path <- substring(normalizePath(source), nchar(normalizePath(root)) + 2L)
    relative <- file.path("coi_state", restore_path)
    target <- file.path(release_dir, relative)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (file.exists(target) || !file.copy(source, target)) stop("Could not stage immutable COI state")
    descriptors[[name]] <- list(
      url = paste(sub("/+$", "", data_base), "releases/COI", release_id, relative, sep = "/"),
      local_path = target, restore_path = restore_path,
      bytes = unname(file.info(target)$size), sha256 = atlas_sha256(target)
    )
    if (grepl("[.]fasta$", source)) descriptors[[name]]$rows <- sum(startsWith(readLines(source, warn = FALSE), ">"))
    if (grepl("[.]csv$", source)) descriptors[[name]]$rows <- nrow(read.csv(source, stringsAsFactors = FALSE))
  }
  counts <- table(state$records$target_group)
  list(schema_version = "1.0.0", retained_records = nrow(state$records),
       target_groups = unique(state$records$target_group), group_counts = as.list(counts),
       policy = list(minimum_per_group = as.integer(Sys.getenv("COI_MIN_PER_GROUP", "100")),
         minimum_total = as.integer(Sys.getenv("COI_TARGET_TOTAL", "10000")),
         monthly_growth = as.numeric(Sys.getenv("COI_MONTHLY_GROWTH", "0.01"))),
       artifacts = descriptors)
}

coi_recover_published_membership <- function(manifest, root = getwd()) {
  artifact <- manifest$artifacts$full_order_anml
  if (is.null(artifact)) artifact <- manifest$artifacts$pair_templates
  if (is.null(artifact)) stop("Legacy COI release has no recoverable template membership artifact")
  rows <- as.data.frame(atlas_read_artifact(artifact, root), stringsAsFactors = FALSE)
  group_column <- if ("target_group" %in% names(rows)) "target_group" else if ("order" %in% names(rows)) "order" else NULL
  template_column <- if ("template" %in% names(rows)) "template" else if ("accession" %in% names(rows)) "accession" else NULL
  if (is.null(group_column) || is.null(template_column)) stop("Legacy COI template artifact lacks group/accession identities")
  orders <- read.csv(file.path(root, "data/orders_22.csv"), stringsAsFactors = FALSE)
  groups <- orders$order[orders$order != "Bivalvia"]
  result <- data.frame(target_group = rows[[group_column]], accession = coi_accession(rows[[template_column]]),
    template = rows[[template_column]], published_release_id = manifest$release_version,
    publication_artifact_sha256 = artifact$sha256,
    original_raw_bytes_available = FALSE, stringsAsFactors = FALSE)
  result <- result[result$target_group %in% groups & nzchar(result$accession) &
    !grepl("^NC_001322[.]1_COX1_reference", result$accession), , drop = FALSE]
  result <- result[!duplicated(paste(result$target_group, result$accession)), , drop = FALSE]
  if (!nrow(result) || !setequal(unique(result$target_group), groups)) stop("Legacy published COI membership is incomplete")
  dir.create(file.path(root, "data/provenance"), recursive = TRUE, showWarnings = FALSE)
  write.csv(result, file.path(root, "data/provenance/coi_published_baseline.csv"), row.names = FALSE)
  result
}

coi_restore_retained_state <- function(root = getwd(), manifest_state = NULL,
                                       allow_bootstrap = FALSE, groups = NULL) {
  if (is.null(manifest_state)) manifest_state <- atlas_read_manifest("COI", root)
  if (isTRUE(manifest_state$stale)) stop("Monthly COI restore requires the verified hosted latest manifest; stale fallback is unsafe")
  manifest <- manifest_state$manifest
  if (!identical(manifest$marker_id, "COI")) stop("Retained state manifest belongs to a different marker")
  state <- manifest$coi_retained_state
  if (is.null(state) || !length(state$artifacts)) {
    if (isTRUE(allow_bootstrap)) return(invisible(list(restored = FALSE, reason = "legacy release has no retained sequence state")))
    stop("Previous COI release has no retained sequence state. An explicit bootstrap is required")
  }
  if (!all(c("records", "ledger") %in% names(state$artifacts))) stop("Incomplete retained COI release contract")
  scratch <- tempfile("coi_restore_")
  dir.create(scratch); on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  staged <- list()
  for (name in names(state$artifacts)) {
    artifact <- state$artifacts[[name]]
    relative <- artifact$restore_path
    if (is.null(relative) || !grepl("^data/(primerminer/|provenance/|derived/coi_reference_panel/)", relative) ||
        grepl("(^|/)\\.\\.(/|$)|[\\\\]", relative)) stop("Unsafe COI state restore path")
    source <- artifact$local_path
    if (!is.null(source) && !startsWith(source, "/")) source <- file.path(root, source)
    if (is.null(source) || !file.exists(source)) source <- atlas_cached_artifact(artifact$url, artifact$sha256)
    if (!identical(tolower(atlas_sha256(source)), tolower(artifact$sha256))) stop("Retained COI artifact checksum failed: ", name)
    target <- file.path(scratch, relative)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(source, target)) stop("Could not prepare COI restore")
    staged[[relative]] <- target
  }
  restored_groups <- unlist(state$target_groups, use.names = FALSE)
  check <- coi_validate_retained_state(scratch, restored_groups)
  archived <- read.csv(file.path(scratch, state$artifacts$records$restore_path), stringsAsFactors = FALSE)
  coi_assert_retention(archived, check$records)
  if (nrow(check$records) != state$retained_records) stop("Restored COI sequence count differs from manifest")
  # Refuse to overwrite local records which were appended since the chosen
  # release. All validation completes before any working state is replaced.
  existing <- lapply(restored_groups, function(group) {
    path <- file.path(root, "data/primerminer/raw", group, paste0(group, "_GB.fasta"))
    if (file.exists(path)) coi_retained_records(root, group) else NULL
  })
  existing <- do.call(rbind, existing)
  if (!is.null(existing) && nrow(existing)) coi_assert_retention(existing, check$records)
  for (relative in names(staged)) {
    destination <- file.path(root, relative)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(staged[[relative]], destination, overwrite = TRUE)) stop("Could not install retained COI state")
  }
  provenance <- list(previous_release_id = manifest$release_version,
    restored_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    retained_records = nrow(check$records), manifest_source = manifest_state$source,
    artifact_sha256 = lapply(state$artifacts, `[[`, "sha256"))
  jsonlite::write_json(provenance, file.path(root, "data/provenance/coi_reference_restore.json"), auto_unbox = TRUE, pretty = TRUE)
  invisible(list(restored = TRUE, previous_release_id = manifest$release_version, retained_records = nrow(check$records)))
}
