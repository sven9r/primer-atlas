# Register application images explicitly, including when a sourced app object
# is launched without shinyAppDir's automatic www handler.
atlas_register_assets <- function(asset_directory = "www", prefix = "atlas-assets") {
  if (length(prefix) != 1L || is.na(prefix) ||
      !grepl("^[A-Za-z][A-Za-z0-9_-]*$", prefix)) {
    stop("Application asset prefix must be one safe path segment.")
  }
  if (length(asset_directory) != 1L || is.na(asset_directory) ||
      !dir.exists(asset_directory)) {
    stop("Application image directory is unavailable.")
  }
  directory <- normalizePath(asset_directory, winslash = "/", mustWork = TRUE)
  shiny::addResourcePath(prefix, directory)
  invisible(prefix)
}

atlas_asset_url <- function(filename, prefix = "atlas-assets") {
  if (length(filename) != 1L || is.na(filename) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", filename)) {
    stop("Application image filename must be one safe path segment.")
  }
  directory <- shiny::resourcePaths()[[prefix]]
  if (is.null(directory) || !file.exists(file.path(directory, filename))) {
    stop("Application image is not available in the registered asset directory: ", filename)
  }
  paste0(prefix, "/", filename)
}
