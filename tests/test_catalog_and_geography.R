#!/usr/bin/env Rscript

source("R/functions.R")

markers <- read.csv("data/catalog/markers.csv", stringsAsFactors = FALSE)
pairs <- read.csv("data/catalog/primer_pairs.csv", stringsAsFactors = FALSE)
facets <- read.csv("data/catalog/primer_pair_facets.csv", stringsAsFactors = FALSE)
unite <- read.csv("data/catalog/unite_primers.csv", stringsAsFactors = FALSE)
geometry <- read.csv("data/derived/marker_pair_geometry.csv", stringsAsFactors = FALSE)
marker_groups <- read.csv("data/catalog/marker_organism_groups.csv", stringsAsFactors = FALSE)
pr2_primers <- read.csv("data/catalog/pr2_18s_primers.csv", stringsAsFactors = FALSE)
pr2_sets <- read.csv("data/catalog/pr2_18s_primer_sets.csv", stringsAsFactors = FALSE)

stopifnot(
  setequal(markers$marker_id, c("COI", "12S_MT", "ITS_FUNGAL", "18S", "16S_PROK", "16S_MT", "28S")),
  sum(pairs$marker_id == "COI") == 25L,
  sum(pairs$marker_id == "ITS_FUNGAL") >= 8L,
  all(c("application", "target_taxon", "environment", "design_intent") %in% facets$facet_type),
  setequal(
    unique(facets$facet_value[facets$facet_type == "application"]),
    c("barcoding", "bulk_community", "edna", "diet")
  ),
  nrow(unite) == 121L,
  sum(unite$direction == "forward") == 45L,
  sum(unite$direction == "reverse") == 76L,
  !anyDuplicated(paste(unite$primer_name, unite$direction, unite$sequence, sep = "|")),
  all(grepl("^[ACGTRYSWKMBDHVNI]+$", unite$sequence)),
  all(c("COI", "ITS_FUNGAL") %in% geometry$marker_id),
  all(c("12S_MT", "COI", "16S_MT") %in% marker_groups$marker_id[marker_groups$group_id == "FISH"]),
  all(c("FUNGI", "PROTISTS", "PHYTOPLANKTON") %in% marker_groups$group_id[marker_groups$marker_id == "28S"])
)

pr2_mappable <- with(
  pr2_sets,
  !is.na(forward_start) & !is.na(forward_end) &
    !is.na(reverse_start) & !is.na(reverse_end) &
    nzchar(forward_sequence) & nzchar(reverse_sequence) &
    reverse_start > forward_end
)
stopifnot(
  markers$status[markers$marker_id == "18S"] == "pilot",
  markers$reference_length[markers$marker_id == "18S"] == 1799L,
  nrow(pr2_primers) == 321L,
  nrow(pr2_sets) == 123L,
  sum(!is.na(pr2_primers$start_yeast)) == 244L,
  sum(pr2_mappable) == 93L,
  all(c("TAReuk454FWD1", "1391F", "E572F", "Uni18SF") %in% pr2_sets$forward_primer),
  all(nzchar(pr2_sets$reference)),
  all(pr2_sets$source_version == "2.1.1")
)

app_source <- paste(readLines("app.R", warn = FALSE), collapse = "\n")
stopifnot(
  grepl('"its_catalog_forward"', app_source, fixed = TRUE),
  grepl('"its_catalog_reverse"', app_source, fixed = TRUE),
  grepl('"its_add_catalog_pair"', app_source, fixed = TRUE),
  grepl('DTOutput("its_map_catalog")', app_source, fixed = TRUE),
  grepl("source-supported pair", app_source, fixed = TRUE),
  grepl("catalog-built session combination", app_source, fixed = TRUE)
  ,grepl('"PR2_18S_008"', app_source, fixed = TRUE)
  ,grepl('"PR2_18S_017"', app_source, fixed = TRUE)
  ,grepl('"PR2_18S_027"', app_source, fixed = TRUE)
  ,grepl('"PR2_18S_040"', app_source, fixed = TRUE)
  ,grepl("Documented 18S combinations", app_source, fixed = TRUE)
  ,grepl("Browse individual fungal-marker primers", app_source, fixed = TRUE)
  ,grepl("Deselect primers to reduce documented combinations", app_source, fixed = TRUE)
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

geography <- nanoparquet::read_parquet(
  "data/pinned/COI/2026-08-14/geography.parquet"
)
taxonomy <- nanoparquet::read_parquet(
  "data/pinned/COI/2026-08-14/taxonomy.parquet"
)
stopifnot(
  nrow(geography) == nrow(taxonomy),
  sum(!is.na(geography$country_or_territory)) >= 500L,
  all(geography$marker_id == "COI")
)
