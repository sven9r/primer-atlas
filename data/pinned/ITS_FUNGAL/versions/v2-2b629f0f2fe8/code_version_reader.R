# Fungal snapshots keep catalogs, sequences, taxonomy, alignments and results
# together. Reading an older snapshot never fills it with today's data.
fungal_version_root <- function(root = getwd()) {
  file.path(root, "data", "pinned", "ITS_FUNGAL", "versions")
}

.fungal_remote_snapshots <- new.env(parent = emptyenv())
fungal_remote_key <- function(root, id) paste(normalizePath(root), Sys.getenv("ATLAS_MANIFEST_BASE_URL", ""), id, sep = "|")
fungal_remote_enabled <- function(root, remote) {
  isTRUE(remote) && nzchar(Sys.getenv("ATLAS_MANIFEST_BASE_URL", "")) &&
    identical(normalizePath(root), normalizePath(getwd()))
}

fungal_read_index <- function(root = getwd(), remote = TRUE) {
  fallback_reason <- NULL
  if (fungal_remote_enabled(root, remote)) {
    state <- tryCatch(atlas_read_manifest("ITS_FUNGAL", root), error = function(e) NULL)
    descriptors <- if (!is.null(state) && !isTRUE(state$stale)) state$manifest$fungal_snapshots else NULL
    if (length(descriptors)) {
      index_artifact <- state$manifest$artifacts$fungal_snapshot_index
      if (is.null(index_artifact)) stop("Hosted fungal release lacks a checksummed index")
      index_artifact$local_path <- NULL
      index_path <- tryCatch(atlas_cached_artifact(index_artifact$url, index_artifact$sha256), error = function(e) NULL)
      if (is.null(index_path)) return(fungal_read_index(root, remote = FALSE))
      hosted_index <- jsonlite::read_json(index_path, simplifyVector = FALSE)
      indexed_ids <- vapply(hosted_index$versions, `[[`, character(1), "id")
      descriptor_ids <- vapply(descriptors, `[[`, character(1), "id")
      if (!setequal(indexed_ids, descriptor_ids)) stop("Hosted fungal index and release identities differ")
      rows <- lapply(descriptors, function(descriptor) {
        if (!grepl("^v[0-9]+-[a-f0-9]{12}$", descriptor$id) || is.null(descriptor$manifest$sha256)) {
          stop("Hosted fungal snapshot has an invalid identity")
        }
        key <- fungal_remote_key(root, descriptor$id)
        if (exists(key, envir = .fungal_remote_snapshots, inherits = FALSE) &&
            !identical(get(key, envir = .fungal_remote_snapshots)$manifest$sha256, descriptor$manifest$sha256)) {
          stop("Hosted immutable fungal snapshot changed its manifest")
        }
        assign(key, descriptor, envir = .fungal_remote_snapshots)
        data.frame(id = descriptor$id, label = descriptor$label, description = descriptor$description,
                   created_at = descriptor$created_at, manifest = descriptor$manifest$url,
                   stringsAsFactors = FALSE)
      })
      result <- do.call(rbind, rows)
      if (anyDuplicated(result$id)) stop("Hosted fungal snapshot index repeats an identity")
      attr(result, "source") <- state$source
      attr(result, "stale") <- FALSE
      attr(result, "hosted") <- TRUE
      return(result)
    }
    fallback_reason <- if (is.null(state) || isTRUE(state$stale)) "Hosted release unavailable" else "Hosted release has no sequence snapshots"
  }
  path <- file.path(fungal_version_root(root), "index.json")
  columns <- c("id", "label", "description", "created_at", "manifest")
  if (!file.exists(path)) return(as.data.frame(setNames(rep(list(character()), length(columns)), columns)))
  index <- jsonlite::read_json(path, simplifyVector = FALSE)
  entries <- index$versions
  if (is.null(entries)) entries <- index
  if (!length(entries)) return(as.data.frame(setNames(rep(list(character()), length(columns)), columns)))
  rows <- lapply(entries, function(entry) {
    values <- lapply(columns, function(column) {
      value <- entry[[column]]
      if (is.null(value)) "" else as.character(value)
    })
    as.data.frame(setNames(values, columns), stringsAsFactors = FALSE)
  })
  result <- do.call(rbind, rows)
  if (anyDuplicated(result$id)) stop("Fungal snapshot index repeats an identity")
  attr(result, "source") <- path
  attr(result, "stale") <- !is.null(fallback_reason)
  attr(result, "hosted") <- FALSE
  attr(result, "fallback_reason") <- fallback_reason
  result
}

fungal_snapshot_artifact <- function(manifest, name, directory) {
  artifact <- manifest$artifacts[[name]]
  if (is.null(artifact)) return(NULL)
  path <- artifact$path
  if (is.null(path) || !nzchar(path) || grepl("(^/|(^|/)\\.\\.(/|$))", path)) {
    stop("Invalid path for fungal snapshot artifact: ", name)
  }
  if (is.null(artifact$sha256) || !grepl("^[a-fA-F0-9]{64}$", artifact$sha256)) {
    stop("Fungal snapshot artifact lacks SHA-256: ", name)
  }
  artifact$local_path <- file.path(directory, path)
  value <- atlas_read_artifact(artifact)
  actual_rows <- if (is.data.frame(value)) nrow(value) else if (is.character(value)) length(value) else NULL
  if (!is.null(actual_rows) && !is.null(artifact$rows) && actual_rows != artifact$rows) {
    stop("Fungal snapshot row count differs: ", name)
  }
  value
}

fungal_read_snapshot <- function(version_id, root = getwd(), remote = TRUE) {
  if (fungal_remote_enabled(root, remote)) {
    key <- fungal_remote_key(root, version_id)
    if (!exists(key, envir = .fungal_remote_snapshots, inherits = FALSE)) fungal_read_index(root)
    if (exists(key, envir = .fungal_remote_snapshots, inherits = FALSE)) {
      descriptor <- get(key, envir = .fungal_remote_snapshots)
      read_remote <- function() {
        manifest_path <- atlas_cached_artifact(descriptor$manifest$url, descriptor$manifest$sha256)
        manifest <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
        if (!identical(manifest$id, version_id) || !identical(manifest$marker_id, "ITS_FUNGAL")) stop("Hosted fungal manifest identity differs")
        result <- list(manifest = manifest)
        for (name in names(manifest$artifacts)) {
          artifact <- descriptor$artifacts[[name]]
          inner <- manifest$artifacts[[name]]
          if (is.null(artifact) || !identical(artifact$sha256, inner$sha256)) stop("Hosted fungal artifact differs from its frozen manifest")
          artifact$local_path <- NULL
          value <- atlas_read_artifact(artifact)
          rows <- if (is.data.frame(value)) nrow(value) else if (is.character(value)) length(value) else NULL
          if (!is.null(rows) && !is.null(inner$rows) && rows != inner$rows) stop("Hosted fungal artifact row count differs")
          result[[name]] <- value
        }
        required <- c("catalog", "sequences", "taxonomy", "geometry", "landmarks", "pairs", "pair_sequences", "sites")
        if (!all(required %in% names(result))) stop("Incomplete hosted fungal snapshot")
        result$manifest$transport_source <- "verified hosted snapshot"
        result
      }
      result <- tryCatch(read_remote(), error = identity)
      if (!inherits(result, "error")) return(result)
      # Fallback must be the exact selected immutable version, never a newer
      # local catalog or another snapshot assembled table by table.
      local_index <- fungal_read_index(root, remote = FALSE)
      local_row <- match(version_id, local_index$id)
      if (!is.na(local_row)) {
        local_manifest <- file.path(fungal_version_root(root), local_index$manifest[local_row])
        if (file.exists(local_manifest) && identical(atlas_sha256(local_manifest), descriptor$manifest$sha256)) {
          result <- fungal_read_snapshot(version_id, root, remote = FALSE)
          result$manifest$transport_source <- "exact pinned fallback"
          return(result)
        }
      }
      stop(conditionMessage(result))
    }
  }
  index <- fungal_read_index(root, remote = FALSE)
  if (length(version_id) != 1L || !version_id %in% index$id) stop("Unknown fungal snapshot identity")
  relative <- index$manifest[match(version_id, index$id)]
  if (grepl("(^/|(^|/)\\.\\.(/|$))", relative)) stop("Invalid fungal snapshot manifest path")
  path <- file.path(fungal_version_root(root), relative)
  manifest <- jsonlite::read_json(path, simplifyVector = FALSE)
  if (!identical(manifest$id, version_id) || !identical(manifest$marker_id, "ITS_FUNGAL")) {
    stop("Fungal snapshot manifest identity differs from its index")
  }
  result <- list(manifest = manifest)
  for (name in names(manifest$artifacts)) {
    result[[name]] <- fungal_snapshot_artifact(manifest, name, dirname(path))
  }
  required <- c("catalog", "sequences", "taxonomy", "geometry", "landmarks")
  if (!all(required %in% names(result))) stop("Incomplete fungal snapshot")
  # Early captured baselines may precede the explorer. Evaluation uses only
  # that snapshot's frozen inputs; its code identity is disclosed separately.
  if (!all(c("pairs", "pair_sequences", "sites") %in% names(result))) {
    if (!exists("fungal_evaluate_panel", mode = "function")) stop("Snapshot has no frozen explorer outputs")
    evaluated <- fungal_evaluate_panel(result$catalog, result$sequences, result$taxonomy)
    result[c("pairs", "pair_sequences", "sites")] <- evaluated[c("pairs", "pair_sequences", "sites")]
    result$manifest$results_derived_at_read <- TRUE
  }
  result$manifest$transport_source <- "bundled pinned snapshot"
  result
}
