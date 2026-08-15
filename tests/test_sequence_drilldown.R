#!/usr/bin/env Rscript

project_root <- normalizePath(getwd(), mustWork = TRUE)

spidprey_path <- file.path(
  project_root,
  "data", "derived", "claimed_nospid_centroid_scores.csv.gz"
)
targeted_path <- file.path(
  project_root,
  "data", "derived", "targeted_zbj_sequence_scores.csv.gz"
)
app_path <- file.path(project_root, "app.R")

stopifnot(file.exists(spidprey_path), file.exists(targeted_path), file.exists(app_path))

spidprey <- read.csv(
  gzfile(spidprey_path),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
ozyptila <- spidprey[
  spidprey$order == "Araneae" &
    spidprey$family == "Thomisidae" &
    spidprey$genus == "Ozyptila",
  ,
  drop = FALSE
]

stopifnot(
  nrow(ozyptila) == 20L,
  sum(ozyptila$cluster_size) == 69L,
  length(unique(ozyptila$accession)) == 20L,
  all(ozyptila$pair_scorable),
  all(nzchar(ozyptila$organism_label)),
  all(nzchar(ozyptila$forward_binding_sequence_primer_oriented)),
  all(nzchar(ozyptila$reverse_binding_sequence_primer_oriented)),
  all(c(
    "Ozyptila americana",
    "Ozyptila atomaria",
    "Ozyptila brevipes",
    "Ozyptila claveata",
    "Ozyptila distans",
    "Ozyptila praticola",
    "Ozyptila trux"
  ) %in% ozyptila$organism_label)
)

targeted <- read.csv(
  gzfile(targeted_path),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
stopifnot(
  nrow(targeted) > 1000L,
  all(c(
    "forward_binding_sequence_primer_oriented",
    "reverse_binding_sequence_primer_oriented",
    "accession",
    "species",
    "pair_penalty"
  ) %in% names(targeted)),
  all(nzchar(targeted$forward_binding_sequence_primer_oriented)),
  all(nzchar(targeted$reverse_binding_sequence_primer_oriented))
)

app_source <- paste(readLines(app_path, warn = FALSE), collapse = "\n")
stopifnot(
  grepl('DTOutput\\("beeprime_taxon_table"', app_source),
  grepl('DTOutput\\("beeprime_sequence_table"', app_source),
  grepl('DTOutput\\("claimed_taxon_table"', app_source),
  grepl('DTOutput\\("claimed_sequence_table"', app_source),
  grepl("beeprime_selected_sequence_data", app_source, fixed = TRUE),
  grepl("download_beeprime_sequences_fasta", app_source, fixed = TRUE),
  grepl("download_claimed_sequences_fasta", app_source, fixed = TRUE),
  grepl("Why the penalty columns contain", app_source, fixed = TRUE)
)

invisible(parse(file = app_path))
message("Filterable taxon-table and exact sequence drill-down checks passed.")
