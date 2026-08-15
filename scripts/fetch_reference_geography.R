#!/usr/bin/env Rscript

project_root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(project_root, ".Rlib"), .libPaths()))
suppressPackageStartupMessages(library(rentrez))
source(file.path(project_root, "R", "functions.R"))
source(file.path(project_root, "R", "data_layer.R"))

taxonomy_path <- file.path(
  project_root, "data", "derived", "gurten2026_centroid_taxonomy.csv"
)
output_path <- file.path(
  project_root, "data", "derived", "reference_geography.csv"
)
cache_path <- file.path(
  project_root, "data", "provenance", "reference_geography_fetched.csv"
)
taxonomy <- read.csv(taxonomy_path, stringsAsFactors = FALSE, check.names = FALSE)
cache <- if (file.exists(cache_path)) {
  read.csv(cache_path, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  data.frame(
    accession = character(), geo_loc_name_raw = character(),
    lat_lon_raw = character(), stringsAsFactors = FALSE
  )
}
if (!nrow(cache)) {
  prior <- tryCatch(atlas_read_manifest("COI", project_root), error = function(e) NULL)
  if (!is.null(prior) && !is.null(prior$manifest$artifacts$geography)) {
    prior_geography <- tryCatch(
      atlas_read_artifact(prior$manifest$artifacts$geography, project_root),
      error = function(e) NULL
    )
    if (!is.null(prior_geography)) {
      prior_geography <- prior_geography[
        !is.na(prior_geography$geo_loc_name_raw) | !is.na(prior_geography$lat_lon_raw),
        c("accession", "geo_loc_name_raw", "lat_lon_raw"), drop = FALSE
      ]
      cache <- as.data.frame(prior_geography, stringsAsFactors = FALSE)
      message("Seeded ", nrow(cache), " previously resolved geography records.")
    }
  }
}

limit <- suppressWarnings(as.integer(Sys.getenv("N_GEO_RECORDS", "500")))
if (!is.finite(limit) || limit < 0L) stop("N_GEO_RECORDS must be a non-negative integer")
pending <- setdiff(taxonomy$accession, cache$accession)
pending <- head(pending, limit)

extract_qualifier <- function(record, names) {
  for (name in names) {
    pattern <- paste0("/", name, "=\"([^\"]+)\"")
    hit <- regexec(pattern, record, perl = TRUE)
    value <- regmatches(record, hit)[[1]]
    if (length(value) == 2L) return(value[2])
  }
  NA_character_
}

if (length(pending)) {
  batches <- split(pending, ceiling(seq_along(pending) / 100L))
  new_rows <- list()
  for (i in seq_along(batches)) {
    message("Fetching geography batch ", i, "/", length(batches))
    text <- entrez_fetch(
      db = "nuccore", id = batches[[i]], rettype = "gb", retmode = "text"
    )
    records <- strsplit(text, "\\n//\\n", perl = TRUE)[[1]]
    rows <- lapply(records, function(record) {
      version <- regmatches(
        record, regexec("(?:^|\\n)VERSION[[:space:]]+([A-Z0-9_.]+)", record, perl = TRUE)
      )[[1]]
      if (length(version) != 2L) return(NULL)
      data.frame(
        accession = trimws(version[2]),
        geo_loc_name_raw = extract_qualifier(record, c("geo_loc_name", "country")),
        lat_lon_raw = extract_qualifier(record, "lat_lon"),
        stringsAsFactors = FALSE
      )
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows)) new_rows[[length(new_rows) + 1L]] <- do.call(rbind, rows)
    Sys.sleep(0.4)
  }
  if (length(new_rows)) {
    cache <- rbind(cache, do.call(rbind, new_rows))
    cache <- cache[!duplicated(cache$accession, fromLast = TRUE), , drop = FALSE]
    write.csv(cache, cache_path, row.names = FALSE, na = "")
  }
}

matched <- match(taxonomy$accession, cache$accession)
raw_geo <- cache$geo_loc_name_raw[matched]
raw_lat_lon <- cache$lat_lon_raw[matched]
normalized <- normalize_geo_loc_name(raw_geo)
coordinates <- parse_lat_lon(raw_lat_lon)
geography <- cbind(
  data.frame(
    marker_id = "COI", accession = taxonomy$accession,
    lat_lon_raw = raw_lat_lon, stringsAsFactors = FALSE
  ),
  normalized,
  coordinates
)
write.csv(geography, output_path, row.names = FALSE, na = "")
message(
  "Wrote ", nrow(geography), " geography rows; ",
  sum(!is.na(geography$country_or_territory)), " have a resolved country or territory."
)
