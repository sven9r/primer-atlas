project_root <- normalizePath(getwd(), mustWork = TRUE)
app_path <- file.path(project_root, "app.R")
invisible(parse(file = app_path))
app_source <- paste(readLines(app_path, warn = FALSE), collapse = "\n")

required_fragments <- c(
  "heatmap-scroll-top",
  "heatmap-scroll-body",
  "wireHeatmapScrollers",
  "heatmap_cell_click",
  "scroll_order_detail",
  "function(message)",
  "heatmap-reading-guide",
  "How to read each cell",
  "How to investigate a warning",
  "Loading selected comparison",
  'class = "comparison-trace-shell"',
  'class = "order-detail-card order-inspector-card"',
  ".order-inspector-card .selectize-dropdown-content",
  "max-height:360px!important",
  "overflow-y:auto!important",
  "cancelOutput = TRUE",
  'freezeReactiveValue(input, "detail_order")',
  'freezeReactiveValue(input, "detail_pair")',
  'class = "order-detail-card"',
  'class = "deep-detail-card"',
  'class = "detail-table-scroll"',
  'DTOutput("claimed_taxon_table", fill = FALSE)',
  'DTOutput("claimed_sequence_table", fill = FALSE)',
  'scrollY = "520px"',
  'scrollY = "380px"'
)
stopifnot(vapply(
  required_fragments,
  function(fragment) grepl(fragment, app_source, fixed = TRUE),
  logical(1)
))
stopifnot(
  grepl("% available", app_source, fixed = TRUE),
  !grepl("% covered", app_source, fixed = TRUE)
)

cat("Order-lens scrolling, legend, and detail-navigation checks passed.\n")
