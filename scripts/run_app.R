#!/usr/bin/env Rscript

project_root <- normalizePath(getwd(), mustWork = TRUE)
local_r_library <- file.path(project_root, ".Rlib")
if (dir.exists(local_r_library)) {
  .libPaths(c(local_r_library, .libPaths()))
}

required <- c(
  "shiny", "bslib", "readr", "dplyr", "stringr", "DT", "PrimerMiner"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop(
    "Install missing app package(s) first: ",
    paste(missing, collapse = ", "),
    "\nRun: install.packages(c(",
    paste(sprintf('"%s"', missing), collapse = ", "),
    "))"
  )
}

shiny::runApp(
  appDir = project_root,
  host = "127.0.0.1",
  port = as.integer(Sys.getenv("PORT", "3838")),
  launch.browser = TRUE
)
