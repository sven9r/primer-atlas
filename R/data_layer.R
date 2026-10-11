atlas_sha256 <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

atlas_manifest_path <- function(marker_id, root = getwd(), release_id = NULL) {
  if (is.null(release_id) || identical(release_id, "latest")) {
    return(file.path(root, "data", "releases", marker_id, "latest.json"))
  }
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", release_id)) stop("Invalid release identity")
  file.path(root, "data", "releases", marker_id, release_id, "manifest.json")
}

atlas_read_manifest <- function(marker_id, root = getwd(), release_id = NULL) {
  remote_base <- sub("/+$", "", Sys.getenv("ATLAS_MANIFEST_BASE_URL", ""))
  fallback <- atlas_manifest_path(marker_id, root, release_id)
  remote_error <- NULL
  if (nzchar(remote_base)) {
    suffix <- if (is.null(release_id) || identical(release_id, "latest")) {
      "latest.json"
    } else paste0(release_id, "/manifest.json")
    url <- paste0(remote_base, "/", marker_id, "/", suffix)
    response <- tryCatch(httr2::request(url) |> httr2::req_timeout(8) |> httr2::req_perform(), error = identity)
    if (!inherits(response, "error") && httr2::resp_status(response) == 200L) {
      return(list(
        manifest = jsonlite::fromJSON(httr2::resp_body_string(response), simplifyVector = FALSE),
        stale = FALSE, source = url
      ))
    }
    remote_error <- if (inherits(response, "error")) conditionMessage(response) else paste("HTTP", httr2::resp_status(response))
  }
  if (!file.exists(fallback)) stop("No release manifest is available for marker ", marker_id)
  list(
    manifest = jsonlite::read_json(fallback, simplifyVector = FALSE),
    stale = nzchar(remote_base), source = fallback, remote_error = remote_error
  )
}

atlas_cached_artifact <- function(url, sha256 = NULL, cache_dir = file.path(tempdir(), "primer-atlas")) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  # The URL identifies marker and release; content identity also separates
  # artifacts after a publication error without sharing a same-named cache.
  identity <- digest::digest(paste(url, sha256, sep = "|"), algo = "sha256", serialize = FALSE)
  destination <- file.path(cache_dir, identity, basename(sub("[?].*$", "", url)))
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(destination)) {
    response <- httr2::request(url) |> httr2::req_timeout(120) |> httr2::req_perform()
    temporary <- tempfile(tmpdir = dirname(destination))
    on.exit(unlink(temporary), add = TRUE)
    writeBin(httr2::resp_body_raw(response), temporary)
    if (!is.null(sha256) && nzchar(sha256) &&
        !identical(tolower(atlas_sha256(temporary)), tolower(sha256))) {
      stop("Artifact checksum failed for ", basename(destination))
    }
    if (!file.rename(temporary, destination)) stop("Could not cache artifact")
  }
  if (!is.null(sha256) && nzchar(sha256)) {
    actual <- digest::digest(file = destination, algo = "sha256", serialize = FALSE)
    if (!identical(tolower(actual), tolower(sha256))) {
      unlink(destination)
      stop("Artifact checksum failed for ", basename(destination))
    }
  }
  destination
}

atlas_read_artifact <- function(artifact, root = getwd()) {
  path <- artifact$local_path
  if (!is.null(path) && !grepl("^/", path)) path <- file.path(root, path)
  if (is.null(path) || !file.exists(path)) {
    if (is.null(artifact$url) || !nzchar(artifact$url)) stop("Artifact has no available local path or URL")
    path <- atlas_cached_artifact(artifact$url, artifact$sha256)
  }
  if (!is.null(artifact$sha256) && nzchar(artifact$sha256) &&
      !identical(tolower(atlas_sha256(path)), tolower(artifact$sha256))) {
    stop("Artifact checksum failed for ", basename(path))
  }
  if (grepl("[.]parquet$", path, ignore.case = TRUE)) return(nanoparquet::read_parquet(path))
  if (grepl("[.]json$", path, ignore.case = TRUE)) return(jsonlite::read_json(path, simplifyVector = TRUE))
  if (grepl("[.](R|py|txt)$", path)) return(readLines(path, warn = FALSE))
  if (grepl("[.](fasta|fa|fna)$", path, ignore.case = TRUE)) {
    lines <- readLines(path, warn = FALSE)
    headers <- which(startsWith(lines, ">"))
    if (!length(headers)) stop("Artifact contains no FASTA records")
    ends <- c(headers[-1] - 1L, length(lines))
    sequences <- vapply(seq_along(headers), function(i) {
      if (ends[i] <= headers[i]) return("")
      toupper(paste0(gsub("[[:space:]]", "", lines[seq.int(headers[i] + 1L, ends[i])]), collapse = ""))
    }, character(1))
    names(sequences) <- sub("^>", "", lines[headers])
    return(sequences)
  }
  readr::read_csv(path, show_col_types = FALSE)
}
