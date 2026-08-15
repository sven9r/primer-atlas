project_root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(project_root, ".Rlib"), .libPaths()))
source(file.path(project_root, "R", "functions.R"))

primers <- read.csv(
  file.path(project_root, "data", "primers.csv"),
  stringsAsFactors = FALSE
)
geometry <- read.csv(
  file.path(project_root, "data", "derived", "pair_geometry.csv"),
  stringsAsFactors = FALSE
)
stored_order <- read.csv(
  file.path(project_root, "data", "derived", "order_pair_summary.csv"),
  stringsAsFactors = FALSE
)
stored_templates <- read.csv(
  file.path(project_root, "data", "derived", "pair_template_scores.csv"),
  stringsAsFactors = FALSE,
  colClasses = c(template = "character")
)
stored_positions <- read.csv(
  file.path(project_root, "data", "derived", "primer_position_scores.csv"),
  stringsAsFactors = FALSE,
  colClasses = c(template = "character")
)

pair_id <- "BEEPRIME"
order_name <- "Hymenoptera"
runtime <- score_primer_pair_primerminer(
  file.path(
    project_root,
    "data", "primerminer", "aligned_reference",
    paste0(order_name, "_COX1_reference_aligned.fasta")
  ),
  primers[primers$pair_id == pair_id, ],
  geometry[geometry$pair_id == pair_id, ]
)

observed <- runtime$order_scores
expected <- stored_order[
  stored_order$order == order_name & stored_order$pair_id == pair_id,
]
summary_fields <- c(
  "n_templates", "n_pair_scorable", "pair_scorable_fraction",
  "forward_median_penalty", "reverse_median_penalty",
  "pair_median_penalty", "pair_p90_penalty",
  "terminal_3_mismatch_fraction"
)
stopifnot(all.equal(
  observed[summary_fields],
  expected[summary_fields],
  tolerance = 1e-10,
  check.attributes = FALSE
))
stopifnot(
  nrow(runtime$pair_template_scores) ==
    nrow(stored_templates[
      stored_templates$order == order_name &
        stored_templates$pair_id == pair_id,
    ]),
  nrow(runtime$primer_position_scores) ==
    nrow(stored_positions[
      stored_positions$order == order_name &
        stored_positions$pair_id == pair_id,
    ])
)

app_source <- paste(
  readLines(file.path(project_root, "app.R"), warn = FALSE),
  collapse = "\n"
)
stopifnot(
  grepl('"order_pair_select"', app_source, fixed = TRUE),
  grepl("order_lens_scores <- reactive", app_source, fixed = TRUE),
  grepl("order_lens_position_scores()", app_source, fixed = TRUE)
)

cat("Runtime custom-pair Order-lens scoring checks passed.\n")
