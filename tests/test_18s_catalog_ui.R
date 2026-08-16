#!/usr/bin/env Rscript

env <- new.env()
source("app.R", local = env)

stopifnot(
  nrow(env$pr2_18s_primers) == 321L,
  nrow(env$pr2_18s_sets) == 123L,
  sum(env$pair_meta$marker_id == "18S") == 93L,
  all(c(
    "PR2_18S_008", "PR2_18S_017",
    "PR2_18S_027", "PR2_18S_040"
  ) %in% env$pair_meta$pair_id)
)

shiny::testServer(env$server, {
  session$setInputs(
    map_organism = "ALL",
    map_marker = "18S",
    map_sort = "binding_site",
    amplicon_range = c(100, 1800),
    pair_select = c(
      "PR2_18S_008", "PR2_18S_017",
      "PR2_18S_027", "PR2_18S_040"
    ),
    application_filter = character(),
    environment_filter = character(),
    intent_filter = character(),
    show_custom = TRUE,
    lineage_marker = "COI",
    lineage_country = "All"
  )
  session$flushReact()

  stopifnot(nrow(filtered_pairs()) == 4L)
  stopifnot(nrow(combination_list_data()) == 123L)
  stopifnot(grepl(
    "Documented 18S combinations (123)",
    output$combination_browser_ui$html,
    fixed = TRUE
  ))
  stopifnot(grepl("primer-map", output$primer_map$html, fixed = TRUE))
  stopifnot(grepl("FU970071", output$marker_coordinate_note$html, fixed = TRUE))
})

message("Accessible PR2 18S combination list and map checks passed.")
