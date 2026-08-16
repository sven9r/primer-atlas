#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !dir.exists(args[[1]])) {
  stop(
    "Usage: Rscript scripts/import_pr2_18s_primers.R /path/to/pr2-primers-shiny",
    call. = FALSE
  )
}
if (!requireNamespace("qs2", quietly = TRUE)) {
  stop("The qs2 package is required to read the upstream PR2 snapshot.", call. = FALSE)
}

source_dir <- normalizePath(args[[1]])
source_commit <- system2(
  "git", c("-C", shQuote(source_dir), "rev-parse", "HEAD"),
  stdout = TRUE
)
snapshot_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
source_url <- "https://github.com/pr2database/pr2-primers-shiny"
source_version <- "2.1.1"

primers <- qs2::qs_read(file.path(source_dir, "data", "primers.qs2"))
sets <- qs2::qs_read(file.path(source_dir, "data", "primer_sets.qs2"))

primers <- primers[primers$gene == "18S rRNA", , drop = FALSE]
sets <- sets[sets$gene == "18S rRNA", , drop = FALSE]

primer_export <- data.frame(
  marker_id = "18S",
  primer_id = primers$primer_id,
  primer_name = primers$name,
  synonyms = primers$synonyms,
  direction = ifelse(primers$direction == "fwd", "forward", "reverse"),
  sequence = primers$sequence,
  start_yeast = primers$start_yeast,
  end_yeast = primers$end_yeast,
  specificity = primers$specificity,
  reference = primers$reference,
  doi = trimws(primers$doi),
  remarks = primers$remark,
  source_version = source_version,
  source_commit = source_commit,
  source_url = source_url,
  source_snapshot_utc = snapshot_utc,
  stringsAsFactors = FALSE
)

set_export <- data.frame(
  marker_id = "18S",
  primer_set_id = sets$primer_set_id,
  primer_set_name = sets$primer_set_name,
  gene_region = sets$gene_region,
  specificity = sets$specificity,
  tested = sets$tested,
  used_for = sets$used_for,
  forward_primer = sets$fwd_name,
  forward_sequence = sets$fwd_seq,
  forward_start = sets$fwd_start,
  forward_end = sets$fwd_end,
  reverse_primer = sets$rev_name,
  reverse_sequence = sets$rev_seq,
  reverse_start = sets$rev_start,
  reverse_end = sets$rev_end,
  amplicon_size = sets$amplicon_size,
  reference = sets$reference,
  doi = trimws(sets$doi),
  metabarcoding_reference = sets$metabarcoding_reference,
  metabarcoding_doi = trimws(sets$metabarcoding_doi),
  remarks = sets$remark,
  source_version = source_version,
  source_commit = source_commit,
  source_url = source_url,
  source_snapshot_utc = snapshot_utc,
  stringsAsFactors = FALSE
)

dir.create(file.path("data", "catalog"), recursive = TRUE, showWarnings = FALSE)
write.csv(
  primer_export,
  file.path("data", "catalog", "pr2_18s_primers.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  set_export,
  file.path("data", "catalog", "pr2_18s_primer_sets.csv"),
  row.names = FALSE,
  na = ""
)

message(
  "Imported ", nrow(primer_export), " 18S primers and ",
  nrow(set_export), " documented primer sets from PR2-primer ",
  source_version, " (", substr(source_commit, 1, 12), ")."
)
