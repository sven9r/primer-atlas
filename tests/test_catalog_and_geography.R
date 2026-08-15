#!/usr/bin/env Rscript

source("R/functions.R")

markers <- read.csv("data/catalog/markers.csv", stringsAsFactors = FALSE)
pairs <- read.csv("data/catalog/primer_pairs.csv", stringsAsFactors = FALSE)
facets <- read.csv("data/catalog/primer_pair_facets.csv", stringsAsFactors = FALSE)
unite <- read.csv("data/catalog/unite_primers.csv", stringsAsFactors = FALSE)
geometry <- read.csv("data/derived/marker_pair_geometry.csv", stringsAsFactors = FALSE)

stopifnot(
  setequal(markers$marker_id, c("COI", "ITS_FUNGAL", "18S", "16S_PROK", "16S_MT", "28S")),
  sum(pairs$marker_id == "COI") == 25L,
  sum(pairs$marker_id == "ITS_FUNGAL") >= 8L,
  all(c("application", "target_taxon", "environment", "design_intent") %in% facets$facet_type),
  setequal(
    unique(facets$facet_value[facets$facet_type == "application"]),
    c("barcoding", "bulk_community", "edna", "diet")
  ),
  nrow(unite) >= 100L,
  all(grepl("^[ACGTRYSWKMBDHVNI]+$", unite$sequence)),
  all(c("COI", "ITS_FUNGAL") %in% geometry$marker_id)
)

geo <- normalize_geo_loc_name(c("USA:Hawaii", "Morocco: Atlas Mountains", NA))
coordinates <- parse_lat_lon(c("21.3069 N 157.8583 W", "31.7917 N 7.0926 W", ""))
stopifnot(
  geo$normalized_location[1] == "USA · Hawaii",
  geo$country_or_territory[2] == "Morocco",
  geo$geographic_resolution[3] == "unresolved",
  abs(coordinates$latitude[1] - 21.3069) < 1e-8,
  abs(coordinates$longitude[1] + 157.8583) < 1e-8,
  is.na(coordinates$latitude[3])
)

message("Marker catalog, multi-facet, ITS import, and geography checks passed.")

geography <- read.csv("data/derived/reference_geography.csv", stringsAsFactors = FALSE, na.strings = c("", "NA"))
stopifnot(
  nrow(geography) == nrow(read.csv("data/derived/gurten2026_centroid_taxonomy.csv")),
  sum(!is.na(geography$country_or_territory)) >= 500L,
  all(geography$marker_id == "COI")
)
