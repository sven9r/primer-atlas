atlas_sha256 <- function(path) {
  unname(tools::md5sum(path))
}

atlas_manifest_path <- function(marker_id, root = getwd()) {
  file.path(root, "data", "releases", marker_id, "latest.json")
}

atlas_read_manifest <- function(marker_id, root = getwd()) {
  remote_base <- sub("/+$", "", Sys.getenv("ATLAS_MANIFEST_BASE_URL", ""))
  fallback <- atlas_manifest_path(marker_id, root)
  remote_error <- NULL
  if (nzchar(remote_base)) {
    url <- paste0(remote_base, "/", marker_id, "/latest.json")
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
  destination <- file.path(cache_dir, basename(sub("[?].*$", "", url)))
  if (!file.exists(destination)) {
    response <- httr2::request(url) |> httr2::req_timeout(120) |> httr2::req_perform()
    writeBin(httr2::resp_body_raw(response), destination)
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
  if (grepl("[.]parquet$", path, ignore.case = TRUE)) return(nanoparquet::read_parquet(path))
  if (grepl("[.]json$", path, ignore.case = TRUE)) return(jsonlite::read_json(path, simplifyVector = TRUE))
  readr::read_csv(path, show_col_types = FALSE)
}
