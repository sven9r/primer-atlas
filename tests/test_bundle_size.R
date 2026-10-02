#!/usr/bin/env Rscript

tracked_roots <- c("app.R", "R", "www", "data/catalog", "data/reference", "data/pinned", "data/releases/COI/latest.json", "data/releases/ITS_FUNGAL/latest.json")
files <- unique(unlist(lapply(tracked_roots, function(path) {
  if (!file.exists(path)) return(character())
  if (dir.exists(path)) list.files(path, recursive = TRUE, full.names = TRUE) else path
})))
files <- files[file.exists(files) & !dir.exists(files)]
bytes <- sum(file.info(files)$size, na.rm = TRUE)
stopifnot(bytes < 100 * 1024^2)

manifest <- jsonlite::read_json("manifest.json", simplifyVector = FALSE)
runtime_packages <- c(
  "shiny", "bslib", "readr", "dplyr", "stringr", "DT",
  "httr2", "jsonlite", "digest", "nanoparquet"
)
stopifnot(
  all(runtime_packages %in% names(manifest$packages)),
  !"readxl" %in% names(manifest$packages),
  !"PrimerMiner" %in% names(manifest$packages)
)
required_files <- c(
  "app.R", "data/catalog/primer_pair_organisms.csv",
  "data/catalog/primer_pair_diet_uses.csv", "data/catalog/primer_pair_facets.csv"
)
stopifnot(all(required_files %in% names(manifest$files)))
for (path in required_files) {
  stopifnot(identical(unname(tools::md5sum(path)), manifest$files[[path]]$checksum))
}
stopifnot(!any(grepl("its_fungal_pilot|regional_coi_audit", names(manifest$files))))
message(sprintf("Deployable application data is %.1f MB (<100 MB gate).", bytes / 1024^2))
