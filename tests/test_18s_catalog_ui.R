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
stopifnot(
  setequal(env$marker_rows_for_organism("ALL", released_only = TRUE)$marker_id, c("COI", "ITS_FUNGAL", "18S")),
  setequal(env$marker_rows_for_organism("FUNGI", released_only = TRUE)$marker_id, c("ITS_FUNGAL", "18S")),
  identical(env$marker_rows_for_organism("PROTISTS", released_only = TRUE)$marker_id, "18S")
)
ui <- as.character(env$ui)
first_start_card <- regmatches(ui, regexpr('<div class="start-card">[\\s\\S]*?</div>', ui, perl = TRUE))
stopifnot(grepl('id="go_18s"', ui, fixed = TRUE),
          grepl('id="go_18s_start"', ui, fixed = TRUE),
          grepl("Open 18S primers", ui, fixed = TRUE),
          grepl('id="go_18s_start"', first_start_card, fixed = TRUE),
          length(regmatches(ui, gregexpr('id="go_18s_start"', ui, fixed = TRUE))[[1L]]) == 1L)

# MockShinySession does not store output suspension settings. Capture the real
# server configuration call while still delegating rendering to Shiny.
configured_output_options <- list()
env$outputOptions <- function(x, name, ...) {
  configured_output_options[[name]] <<- list(...)
  shiny::outputOptions(x, name, ...)
}

shiny::testServer(env$server, {
  # Pair-selection inputs must be created even before a user opens Filters.
  # Hidden selector suspension previously left both fresh COI and 18S maps empty.
  stopifnot(identical(configured_output_options$pair_selector_ui$suspendWhenHidden, FALSE),
            is.null(input$pair_select), is.null(input$pair_select_18s))
  # Supply the inputs present in the static page, but not the pair inputs that
  # the dynamic selector output is responsible for creating.
  session$setInputs(map_organism = "ALL", map_marker = "COI", map_sort = "binding_site",
    amplicon_range = c(100, 1800), application_filter = character(),
    environment_filter = character(), intent_filter = character(), show_custom = TRUE,
    lineage_marker = "COI", lineage_country = "All")
  session$flushReact()
  initial_selectors <- output$pair_selector_ui$html
  stopifnot(grepl('id="pair_select"', initial_selectors, fixed = TRUE),
            grepl('id="pair_select_18s"', initial_selectors, fixed = TRUE),
            grepl('value="ANML" checked', initial_selectors, fixed = TRUE),
            grepl('<option value="PR2_18S_008" selected>', initial_selectors, fixed = TRUE))
  session$setInputs(
    map_organism = "ALL",
    map_marker = "18S",
    map_sort = "binding_site",
    amplicon_range = c(100, 1800),
    pair_select_18s = c(
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

  # A visible entry restores the existing pilot from another organism/marker
  # and stale COI facets. Capture actual client messages rather than assuming
  # a test-session input changes automatically when update*Input is called.
  sent <- list()
  session$sendInputMessage <- function(inputId, message) sent[[inputId]] <<- message
  session$setInputs(
    map_organism = "ARTHROPODS", map_marker = "COI", pair_select = "ANML",
    application_filter = "diet", environment_filter = "freshwater",
    intent_filter = "exclusion_blocking", amplicon_range = c(100, 350)
  )
  session$flushReact()
  session$setInputs(go_18s = 1L)
  session$flushReact()
  defaults <- c("PR2_18S_008", "PR2_18S_017", "PR2_18S_027", "PR2_18S_040")
  stopifnot(sent$map_organism$value == "ALL", sent$map_marker$value == "18S",
            grepl('value="18S"', sent$map_marker$options, fixed = TRUE),
            sent$main_nav$value == "find", sent$map_sort$value == "binding_site",
            length(sent$application_filter$value) == 0L,
            length(sent$environment_filter$value) == 0L,
            length(sent$intent_filter$value) == 0L,
            setequal(sent$pair_select_18s$value, defaults),
            length(sent$pair_select_18s$value) == 4L,
            identical(map_navigation_marker(), "18S"))
  # Slider updates serialize numbers as strings; the browser returns numeric
  # values. Model that conversion when applying captured updates in this test.
  limits <- as.numeric(sent$amplicon_range$value)
  bp <- env$pair_meta$aligned_amplicon_bp[env$pair_meta$marker_id == "18S"]
  stopifnot(length(limits) == 2L, limits[1L] <= min(bp), limits[2L] >= max(bp))

  # The group reset can arrive before the marker reset. The old COI value must
  # not undo the explicit 18S request or select all 93 mapped sets at once.
  session$setInputs(map_organism = "ALL", map_marker = "COI")
  session$flushReact()
  stopifnot(sent$map_marker$value == "18S", setequal(sent$pair_select_18s$value, defaults),
            length(sent$pair_select_18s$value) == 4L, is.null(map_navigation_marker()))
  session$setInputs(map_marker = "18S", pair_select_18s = defaults,
    amplicon_range = limits, application_filter = character(),
    environment_filter = character(), intent_filter = character())
  session$flushReact()
  stopifnot(nrow(filtered_pairs()) == 4L, all(filtered_pairs()$marker_id == "18S"),
            nrow(combination_list_data()) == 123L,
            grepl("FU970071", output$marker_coordinate_note$html, fixed = TRUE),
            grepl("primer-map", output$primer_map$html, fixed = TRUE))

  # The landing button shares the same route and remains usable when already
  # on the all-group 18S pilot; it cannot leave a pending navigation override.
  session$setInputs(go_18s_start = 1L)
  session$flushReact()
  stopifnot(sent$map_marker$value == "18S", sent$main_nav$value == "find",
            setequal(sent$pair_select_18s$value, defaults), is.null(map_navigation_marker()))
})

message("Accessible PR2 18S entry, clean scope/filter reset, combination list and map checks passed.")
