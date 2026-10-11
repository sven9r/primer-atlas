#!/usr/bin/env Rscript

# Exercise the source()+runApp(app-object) launch that lacks automatic www
# serving, using the actual application banner and a clean copied asset folder.
check_app_assets <- function() {
  root <- normalizePath(getwd(), mustWork = TRUE)
  .libPaths(c(file.path(root, ".Rlib"), .libPaths()))
  required <- c("shiny", "httpuv", "httr2", "jsonlite", "digest")
  stopifnot(all(vapply(required, requireNamespace, logical(1), quietly = TRUE)))
  scratch <- tempfile("atlas-assets-http-")
  dir.create(scratch)
  child_pid <- NA_integer_
  on.exit({
    if (is.finite(child_pid)) try(tools::pskill(child_pid, signal = 15L), silent = TRUE)
    unlink(scratch, recursive = TRUE)
  }, add = TRUE)

  manifest <- jsonlite::read_json(file.path(root, "manifest.json"), simplifyVector = FALSE)
  images <- c("evolab-berkeley.png", "fulbright.png", "walter-benjamin-dfg.png", "uc-berkeley.png")
  image_paths <- paste0("www/", images)
  stopifnot(all(image_paths %in% names(manifest$files)))
  stopifnot("R/app_assets.R" %in% names(manifest$files),
            identical(unname(tools::md5sum(file.path(root, "R/app_assets.R"))),
                      manifest$files[["R/app_assets.R"]]$checksum))
  dir.create(file.path(scratch, "www"))
  for (path in image_paths) {
    stopifnot(identical(unname(tools::md5sum(file.path(root, path))), manifest$files[[path]]$checksum),
              file.copy(file.path(root, path), file.path(scratch, path)))
  }
  stopifnot(file.copy(file.path(root, "R/app_assets.R"), file.path(scratch, "app_assets.R")))
  expressions <- parse(file.path(root, "app.R"))
  is_banner <- vapply(expressions, function(expr) is.call(expr) &&
      identical(expr[[1L]], as.name("<-")) && identical(expr[[2L]], as.name("institution_banner")), logical(1))
  stopifnot(sum(is_banner) == 1L)
  saveRDS(expressions[[which(is_banner)]], file.path(scratch, "banner.rds"))
  saveRDS(.libPaths(), file.path(scratch, "libraries.rds"))
  port <- httpuv::randomPort()
  writeLines(as.character(port), file.path(scratch, "port.txt"))
  child_script <- file.path(scratch, "serve.R")
  writeLines(c(
    "args <- commandArgs(trailingOnly = TRUE)",
    "setwd(args[1L])",
    ".libPaths(readRDS('libraries.rds'))",
    "library(shiny)",
    "source('app_assets.R')",
    "atlas_register_assets('www')",
    "eval(readRDS('banner.rds'))",
    "app <- shinyApp(fluidPage(institution_banner), function(input, output, session) {})",
    "stopifnot(is.null(app$staticPaths))",
    "writeLines(as.character(Sys.getpid()), 'pid.txt')",
    "runApp(app, host = '127.0.0.1', port = as.integer(readLines('port.txt')), launch.browser = FALSE)"
  ), child_script)
  log <- file.path(scratch, "server.log")
  system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(child_script), shQuote(scratch)),
          stdout = log, stderr = log, wait = FALSE)
  base_url <- paste0("http://127.0.0.1:", port)
  page <- NULL
  deadline <- Sys.time() + 20
  while (Sys.time() < deadline) {
    if (file.exists(file.path(scratch, "pid.txt"))) child_pid <- as.integer(readLines(file.path(scratch, "pid.txt")))
    page <- tryCatch(httr2::request(base_url) |> httr2::req_timeout(1) |> httr2::req_perform(),
                     error = function(e) NULL)
    if (!is.null(page)) break
    Sys.sleep(.1)
  }
  if (is.null(page)) stop("Asset HTTP fixture did not start: ", paste(readLines(log), collapse = "\n"))
  html <- httr2::resp_body_string(page)
  stopifnot(httr2::resp_status(page) == 200L)
  png_signature <- as.raw(c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L))
  for (image in images) {
    path <- paste0("atlas-assets/", image)
    stopifnot(grepl(paste0('src="', path, '"'), html, fixed = TRUE))
    response <- httr2::request(paste0(base_url, "/", path)) |>
      httr2::req_url_query(uncached = basename(scratch)) |>
      httr2::req_headers(`Cache-Control` = "no-cache") |>
      httr2::req_timeout(5) |> httr2::req_perform()
    bytes <- httr2::resp_body_raw(response)
    stopifnot(httr2::resp_status(response) == 200L,
              grepl("^image/png", httr2::resp_header(response, "content-type")),
              identical(bytes[seq_along(png_signature)], png_signature),
              identical(digest::digest(bytes, algo = "md5", serialize = FALSE),
                        manifest$files[[paste0("www/", image)]]$checksum))
  }
  message("All four copied institution images return uncached HTTP200 image/png with original bytes from a sourced Shiny app object.")
}

check_app_assets()
