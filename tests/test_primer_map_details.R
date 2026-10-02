#!/usr/bin/env Rscript

app <- source("app.R")$value
metadata <- catalog_download_rows()
pr2_download <- metadata[metadata$marker == "18S", ]
pr2_reported <- pr2_map_sets$amplicon_size[match(pr2_download$pair_id, pr2_map_sets$pair_id)]
stopifnot(
  nrow(metadata) == nrow(pair_meta),
  all(c("map_start", "map_end", "forward_sequence_5to3", "source_urls",
        "community_source_urls", "dietary_source_urls", "dietary_predator",
        "dietary_sample", "dietary_evidence_status") %in% names(metadata)),
  metadata$dietary_predator[metadata$pair_id == "ZBJ_ART"] == "bats",
  metadata$dietary_sample[metadata$pair_id == "NOSPID"] == "whole-spider DNA extracts",
  grepl("10.1111/j.1755-0998.2010.02920.x",
        metadata$dietary_source_urls[metadata$pair_id == "ZBJ_ART"], fixed = TRUE),
  !any(grepl("primer_sheet", metadata$source_keys, fixed = TRUE)),
  !"primer_sheet" %in% public_citations$key,
  all(pr2_download$catalog_status == "documented_PR2_set"),
  identical(pr2_download$reported_amplicon_bp, as.numeric(pr2_reported)),
  !anyDuplicated(metadata$pair_id),
  all(nzchar(metadata$forward_sequence_5to3)),
  all(nzchar(metadata$reverse_sequence_5to3)),
  all(nzchar(metadata$source_urls)),
  all(unlist(strsplit(pair_organisms$source_keys, "|", fixed = TRUE)) %in% citations$key),
  all(unlist(strsplit(diet_uses$source_key, "|", fixed = TRUE)) %in% citations$key)
)

shiny::testServer(app$serverFuncSource(), {
  fungal_18s <- pair_choices_for_marker("18S", "FUNGI", apply_filters = FALSE)
  stopifnot(setequal(unname(fungal_18s), c("PR2_18S_092", "PR2_18S_149")))

  session$setInputs(
    map_organism = "ALL", map_marker = "COI",
    pair_select = "FISH_F2_R1", amplicon_range = c(110, 710),
    sodium = 50, ta_offset = 4
  )

  labels <- names(pair_choices_for_marker("COI"))
  stopifnot(
    any(grepl("FishF2 + FishR1 · Fishes target; arthropod off-target check only", labels, fixed = TRUE)),
    any(grepl("MollCOI253 · Marine mollusks target; arthropod off-target check only", labels, fixed = TRUE)),
    any(grepl("Vertebrate COI · Vertebrata target; arthropod off-target check only", labels, fixed = TRUE))
  )

  map_html <- paste(as.character(output$primer_map), collapse = " ")
  stopifnot(
    grepl("Estimated primer Tm (50 mM salt): FishF2 52.6", map_html, fixed = TRUE),
    grepl("Suggested pair Ta (starting estimate): 47.6–49.6 °C", map_html, fixed = TRUE),
    grepl("Reported Ta: 52 °C", map_html, fixed = TRUE),
    grepl("Arthropod reference: off-target binding check only", map_html, fixed = TRUE),
    grepl("<circle cx=\"16\"", map_html, fixed = TRUE),
    grepl("class=\"primer-tooltip-target\" data-tooltip=\"Map row 1: FishF2", map_html, fixed = TRUE)
  )

  session$setInputs(map_organism = "ARTHROPODS")
  arthropod_labels <- names(pair_choices_for_marker("COI"))
  stopifnot(
    !any(grepl("FishF2", arthropod_labels, fixed = TRUE)),
    any(grepl("ZBJ-Art", arthropod_labels, fixed = TRUE)),
    any(grepl("Spidprey", arthropod_labels, fixed = TRUE)),
    grepl("No primer pairs match", paste(as.character(output$primer_map), collapse = " "), fixed = TRUE)
  )

  session$setInputs(map_organism = "ALL")
  stopifnot(grepl("FishF2 52.6", paste(as.character(output$primer_map), collapse = " "), fixed = TRUE))

  session$setInputs(sodium = 75)
  stopifnot(grepl(
    "Estimated primer Tm (75 mM salt):",
    paste(as.character(output$primer_map), collapse = " "),
    fixed = TRUE
  ))
})

shiny::testServer(app$serverFuncSource(), {
  sent <- list()
  session$sendInputMessage <- function(inputId, message) sent[[inputId]] <<- message
  eligible <- pair_organisms$pair_id[pair_organisms$group_id == "ARTHROPODS"]
  dietary <- c("ANML", "FWH1", "ZBJ_ART", "NOSPID")
  session$setInputs(
    map_organism = "ARTHROPODS", map_marker = "COI",
    pair_select = eligible, amplicon_range = c(110, 710)
  )
  session$setInputs(application_filter = "diet")
  stopifnot(setequal(sent$pair_select$value, dietary))
  session$setInputs(pair_select = dietary, environment_filter = "marine")
  stopifnot(length(sent$pair_select$value) == 0L)
  session$setInputs(
    pair_select = character(), environment_filter = character(),
    application_filter = character()
  )
  stopifnot(setequal(sent$pair_select$value, eligible))
  manually_selected <- setdiff(eligible, "BF2_BR2")
  session$setInputs(pair_select = manually_selected, map_sort = "binding_site")
  stopifnot(setequal(sent$pair_select$value, manually_selected))
})

cat("Primer-map usage and temperature details passed.\n")
