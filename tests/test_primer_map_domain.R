#!/usr/bin/env Rscript

project_root <- normalizePath(getwd(), mustWork = TRUE)
source(file.path(project_root, "R", "functions.R"))

pairs <- read.csv(
  file.path(project_root, "data", "derived", "pair_geometry.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
folmer <- read.csv(
  file.path(project_root, "data", "derived", "folmer_region.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
bindings <- read.csv(
  file.path(project_root, "data", "derived", "primer_bindings.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

inside <- pairs[pairs$pair_id %in% c("NOSPID", "NOSPI2_LAURELIN", "NOPLANT"), ]
inside_domain <- primer_map_domain(inside, folmer$start, folmer$end)
stopifnot(
  identical(unname(inside_domain), c(0, 700))
)

default_ids <- c(
  "MCO", "LERAY_XT", "ZBJ_ART", "FWH2", "BF2_BR2", "ANML",
  "BEEPRIME", "NOSPID", "NOSPI2_LAURELIN"
)
default_domain <- primer_map_domain(
  pairs[pairs$pair_id %in% default_ids, ],
  folmer$start,
  folmer$end
)
stopifnot(
  identical(unname(default_domain), c(0, 800))
)

mollusk_pair <- pairs[pairs$pair_id == "MOLLCOI253", ]
mollusk_bindings <- bindings[bindings$pair_id == "MOLLCOI253", ]
stopifnot(
  nrow(mollusk_pair) == 1L,
  nrow(mollusk_bindings) == 2L,
  mollusk_pair$forward_start < mollusk_pair$reverse_start,
  abs(mollusk_pair$aligned_amplicon_bp - 253L) <= 5L,
  all(
    mollusk_bindings$placement_mode ==
      "pair-aware published-size constraint"
  ),
  all(mollusk_bindings$pair_size_delta_bp <= 5L)
)

late_pair <- data.frame(
  forward_start = 1181,
  forward_end = 1205,
  reverse_start = 1276,
  reverse_end = 1300
)
late_domain <- primer_map_domain(late_pair, folmer$start, folmer$end)
stopifnot(
  identical(unname(late_domain), c(0, 1300))
)

app_source <- paste(
  readLines(file.path(project_root, "app.R"), warn = FALSE),
  collapse = "\n"
)
positions <- c(
  marker = regexpr('"map_marker", "Marker"', app_source, fixed = TRUE)[1],
  application = regexpr('"application_filter", "Applications"', app_source, fixed = TRUE)[1],
  map_sort = regexpr('"map_sort", "Order rows by"', app_source, fixed = TRUE)[1],
  select_all = regexpr('"select_all", "Select all"', app_source, fixed = TRUE)[1],
  amplicon_range = regexpr(
    '"amplicon_range", "Amplicon + primers (bp)"',
    app_source,
    fixed = TRUE
  )[1],
  pair_select = regexpr('"pair_select", "Primer pairs"', app_source, fixed = TRUE)[1]
)
stopifnot(
  all(positions > 0),
  positions["marker"] < positions["application"],
  positions["application"] < positions["map_sort"],
  positions["map_sort"] < positions["select_all"],
  positions["select_all"] < positions["amplicon_range"],
  positions["amplicon_range"] < positions["pair_select"],
  grepl("primer_map_domain(", app_source, fixed = TRUE),
  grepl('paste0("View ", axis_min, "–", axis_max, " bp")', app_source, fixed = TRUE),
  grepl("sorted_filtered_pairs <- reactive({", app_source, fixed = TRUE),
  grepl("sorted_pair_choices <- reactive({", app_source, fixed = TRUE),
  grepl(
    "number_primer_choice_labels(\n      sorted_pair_choices(),",
    app_source,
    fixed = TRUE
  ),
  grepl("number_primer_choice_labels(", app_source, fixed = TRUE),
  grepl('updateCheckboxGroupInput(', app_source, fixed = TRUE),
  grepl("number_x <- 16", app_source, fixed = TRUE),
  grepl("left <- 112", app_source, fixed = TRUE),
  grepl("primer-hover-tooltip", app_source, fixed = TRUE),
  grepl("pinnedPrimerTooltipTarget", app_source, fixed = TRUE),
  grepl('class = "primer-tooltip-target"', app_source, fixed = TRUE),
  grepl('`data-tooltip` = tooltip', app_source, fixed = TRUE),
  grepl('"\\nTags: "', app_source, fixed = TRUE),
  grepl('"Map row ", i', app_source, fixed = TRUE),
  grepl("Sorting and filtering recalculate this number.", app_source, fixed = TRUE),
  grepl('class = "legend-number"', app_source, fixed = TRUE)
)

sorting_fixture <- data.frame(
  pair_label = c("Zulu", "Alpha", "Middle"),
  pair_start = c(300, 100, 200),
  pair_end = c(450, 500, 260),
  folmer_overlap_bp = c(100, 300, 50),
  aligned_amplicon_bp = c(151, 401, 61),
  stringsAsFactors = FALSE
)
stopifnot(
  identical(
    sort_primer_map_rows(sorting_fixture, "binding_site")$pair_label,
    c("Alpha", "Middle", "Zulu")
  ),
  identical(
    sort_primer_map_rows(sorting_fixture, "alphabetical")$pair_label,
    c("Alpha", "Middle", "Zulu")
  ),
  identical(
    sort_primer_map_rows(sorting_fixture, "folmer_overlap")$pair_label,
    c("Alpha", "Zulu", "Middle")
  ),
  identical(
    sort_primer_map_rows(sorting_fixture, "amplicon_length")$pair_label,
    c("Middle", "Zulu", "Alpha")
  )
)

tag_fixture <- primer_use_tags(
  "vertebrate diet metabarcoding / eDNA",
  "Vertebrata"
)
stopifnot(
  "diet metabarcoding" %in% tag_fixture,
  "eDNA metabarcoding" %in% tag_fixture,
  "target: Vertebrata" %in% tag_fixture
)

choice_fixture <- data.frame(
  pair_id = c("A", "B", "C"),
  pair_label = c("Alpha", "Beta", "Gamma"),
  stringsAsFactors = FALSE
)
display_fixture <- data.frame(
  pair_id = c("C", "A"),
  stringsAsFactors = FALSE
)
numbered_choices <- number_primer_choice_labels(
  choice_fixture,
  display_fixture
)
stopifnot(
  identical(
    names(numbered_choices),
    c("2 · Alpha", "Beta", "1 · Gamma")
  ),
  identical(unname(numbered_choices), c("A", "B", "C"))
)

alphabetical_choice_fixture <- choice_fixture[c(3, 1, 2), , drop = FALSE]
alphabetical_choices <- number_primer_choice_labels(
  alphabetical_choice_fixture,
  display_fixture
)
stopifnot(
  identical(
    names(alphabetical_choices),
    c("1 · Gamma", "2 · Alpha", "Beta")
  ),
  identical(unname(alphabetical_choices), c("C", "A", "B"))
)

message("Dynamic primer-map domain and control-order checks passed.")
