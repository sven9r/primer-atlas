#!/usr/bin/env Rscript

source("R/functions.R")
pairs <- read.csv("data/catalog/primer_pairs.csv", stringsAsFactors = FALSE)
facets <- read.csv("data/catalog/primer_pair_facets.csv", stringsAsFactors = FALSE)
landmarks <- read.csv("data/catalog/marker_landmarks.csv", stringsAsFactors = FALSE)
marker_groups <- read.csv("data/catalog/marker_organism_groups.csv", stringsAsFactors = FALSE)

or_result <- filter_pair_facets(
  pairs$pair_id, facets,
  list(application = c("edna", "diet"))
)
and_result <- filter_pair_facets(
  pairs$pair_id, facets,
  list(application = c("edna", "diet"), design_intent = "exclusion_blocking")
)
stopifnot(length(or_result) >= length(and_result), length(and_result) > 0L)
stopifnot(all(and_result %in% or_result))
stopifnot(setequal(
  landmarks$label[landmarks$marker_id == "ITS_FUNGAL"],
  c("18S", "ITS1", "5.8S", "ITS2", "28S")
))

coi <- clean_sequence(parse_fasta("data/reference/coi_reference_NC_001322.1.fasta")[[1]])
its <- clean_sequence(parse_fasta("data/reference/its_fungal_reference_FN812768.2.fasta")[[1]])
coi_pair <- read.csv("data/primers.csv", stringsAsFactors = FALSE)
coi_pair <- coi_pair[coi_pair$pair_id == "FOLMER", ]
cross_marker <- tryCatch({
  locate_primer_pair(
    its,
    coi_pair$sequence[coi_pair$direction == "forward"],
    coi_pair$sequence[coi_pair$direction == "reverse"]
  )
  FALSE
}, error = function(e) TRUE)
stopifnot(cross_marker)

app_source <- paste(readLines("app.R", warn = FALSE), collapse = "\n")
stopifnot(
  all(c("FISH", "BIRDS") %in% marker_groups$group_id[marker_groups$marker_id == "12S_MT"]),
  grepl('"Region comparison"', app_source, fixed = TRUE),
  grepl('"region_position_differences"', app_source, fixed = TRUE),
  grepl('"region_sequences"', app_source, fixed = TRUE),
  grepl('"map_organism"', app_source, fixed = TRUE)
)

message("Marker switching, landmark, facet, and cross-marker rejection checks passed.")
