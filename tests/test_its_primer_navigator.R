#!/usr/bin/env Rscript

env <- new.env()
source("app.R", local = env)

stopifnot(
  nrow(env$unite_catalog) == 121L,
  nrow(env$its_documented_members) == 9L,
  nrow(env$its_documented_pair_index) == 8L,
  all(c("gene_locus", "target", "remarks", "reference") %in%
        names(env$its_documented_members))
)

shiny::testServer(env$server, {
  session$setInputs(
    map_organism = "FUNGI",
    map_marker = "ITS_FUNGAL",
    map_sort = "binding_site",
    amplicon_range = c(100, 800),
    pair_select = env$its_documented_pair_ids,
    its_combo_forward = unique(env$its_documented_pair_index$forward_primer),
    its_combo_reverse = unique(env$its_documented_pair_index$reverse_primer),
    unite_region_filter = "ALL",
    application_filter = character(),
    environment_filter = character(),
    intent_filter = character(),
    show_custom = TRUE,
    lineage_marker = "COI",
    lineage_country = "All"
  )
  session$flushReact()

  stopifnot(length(its_filtered_pair_ids()) == 8L)
  stopifnot(nrow(combination_list_data()) == 8L)
  stopifnot(nrow(filtered_pairs()) == 8L)
  stopifnot(grepl(
    "ITS pairs are controlled by the primer checkboxes above the map",
    output$pair_selector_ui$html,
    fixed = TRUE
  ))
  stopifnot(grepl(
    "Deselect primers to reduce documented combinations",
    output$combination_browser_ui$html,
    fixed = TRUE
  ))

  session$setInputs(its_combo_reverse = c("ITS4B", "LR21"))
  session$flushReact()
  stopifnot(
    setequal(
      its_filtered_pair_ids(),
      c("ITS1F_ITS4B", "ITS1F_LR21")
    ),
    nrow(combination_list_data()) == 2L,
    nrow(filtered_pairs()) == 2L
  )
})

message("ITS primer-first navigator and combination-reduction checks passed.")
