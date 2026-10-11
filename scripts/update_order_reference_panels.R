#!/usr/bin/env Rscript
# Acquire through the feature-aware, append-only Entrez builder. Its panel
# manifest, original headers and sequence bytes are the retention authority.
root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(root, ".Rlib"), .libPaths()))
source("R/data_layer.R")
source("R/coi_retention.R")

panel_dir <- Sys.getenv("COI_PANEL_DIR", "data/derived/coi_reference_panel")
target_total <- as.integer(Sys.getenv("COI_TARGET_TOTAL", "10000"))
minimum <- as.integer(Sys.getenv("COI_MIN_PER_GROUP", "100"))
growth <- as.numeric(Sys.getenv("COI_MONTHLY_GROWTH", "0.01"))
release_id <- Sys.getenv("RELEASE_ID", format(Sys.Date(), "%Y-%m-%d"))
if (!is.finite(target_total) || target_total < 22L * minimum || !is.finite(minimum) || minimum < 100L) stop("COI target must support at least100 records per arthropod target group")
if (!is.finite(growth) || growth < 0.01) stop("COI monthly growth must be at least0.01")
python <- Sys.getenv("COI_PYTHON", Sys.which("python3"))
if (!nzchar(python)) stop("Python3 is required for the COI Entrez builder")
args <- c("scripts/build_coi_reference_panel.py", "--output-dir", panel_dir,
          "--target-total", target_total, "--minimum-per-group", minimum,
          "--growth", growth, "--release-id", release_id)
previous_dir <- Sys.getenv("COI_PREVIOUS_PANEL_DIR", panel_dir)
if (file.exists(file.path(previous_dir, "panel.fasta"))) {
  args <- c(args, "--previous-dir", previous_dir)
} else if (dir.exists("data/primerminer/raw")) {
  args <- c(args, "--legacy-raw-dir", "data/primerminer/raw")
}
if (identical(tolower(Sys.getenv("COI_OFFLINE", "false")), "true")) args <- c(args, "--offline")
if (identical(tolower(Sys.getenv("COI_REFRESH", "false")), "true")) args <- c(args, "--refresh")
published_baseline <- Sys.getenv("COI_PUBLISHED_BASELINE", "data/provenance/coi_published_baseline.csv")
if (!file.exists(published_baseline)) published_baseline <- ""
if (nzchar(published_baseline)) args <- c(args, "--published-baseline", published_baseline)
use_existing <- identical(tolower(Sys.getenv("COI_USE_EXISTING_PANEL", "false")), "true")
if (identical(tolower(Sys.getenv("DRY_RUN", "false")), "true")) {
  message("COI acquisition dry run; no sequences or provenance were changed")
  quit(status = 0L)
}
if (!use_existing) {
  status <- system2(python, vapply(as.character(args), shQuote, character(1)))
  if (status != 0L) stop("Feature-aware COI acquisition failed; existing retained inputs were not replaced")
}
required <- c("panel.fasta", "panel_manifest.csv", "ledger.csv", "taxonomy.csv", "growth_qa.csv", "run_manifest.json")
if (!all(file.exists(file.path(panel_dir, required)))) stop("Incomplete acquired COI panel")
sequences <- coi_read_fasta(file.path(panel_dir, "panel.fasta"))
manifest <- read.csv(file.path(panel_dir, "panel_manifest.csv"), stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c("target_group", "accession", "sequence_sha256", "header_sha256", "first_release_id", "added_at_utc")
if (!all(required_columns %in% names(manifest))) stop("COI panel manifest lacks immutable retention provenance")
accessions <- coi_accession(names(sequences))
if (anyDuplicated(manifest$accession) || !setequal(accessions, manifest$accession)) stop("COI panel manifest and sequence identities differ")
manifest <- manifest[match(accessions, manifest$accession), , drop = FALSE]
actual_hash <- vapply(unname(sequences), digest::digest, character(1), algo = "sha256", serialize = FALSE)
header_hash <- vapply(names(sequences), digest::digest, character(1), algo = "sha256", serialize = FALSE)
if (!identical(tolower(actual_hash), tolower(manifest$sequence_sha256)) ||
    !identical(tolower(header_hash), tolower(manifest$header_sha256))) stop("COI panel bytes differ from acquired provenance")
orders <- read.csv("data/orders_22.csv", stringsAsFactors = FALSE)
groups <- orders$order[orders$order != "Bivalvia"]
if (!setequal(unique(manifest$target_group), groups)) stop("Acquired COI panel does not contain the22 arthropod target groups")
# Verify all old records before replacing any group file.
for (group in groups) {
  path <- file.path("data/primerminer/raw", group, paste0(group, "_GB.fasta"))
  if (!file.exists(path)) next
  old <- coi_retained_records(root, group)
  current <- manifest[manifest$target_group == group, , drop = FALSE]
  coi_assert_retention(old, current)
}
write_fasta <- function(seqs, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lines <- unlist(lapply(seq_along(seqs), function(i) {
    start <- seq(1L, nchar(seqs[[i]]), by = 80L)
    c(paste0(">", names(seqs)[i]), substring(seqs[[i]], start, pmin(start + 79L, nchar(seqs[[i]]))))
  }), use.names = FALSE)
  temporary <- tempfile(tmpdir = dirname(path))
  writeLines(lines, temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) stop("Could not install retained COI group input")
}
for (group in groups) {
  seqs <- sequences[manifest$target_group == group]
  write_fasta(seqs, file.path("data/primerminer/raw", group, paste0(group, "_GB.fasta")))
}
dir.create("data/provenance", recursive = TRUE, showWarnings = FALSE)
if (!file.copy(file.path(panel_dir, "ledger.csv"), "data/provenance/coi_reference_accession_ledger.csv", overwrite = TRUE) ||
    !file.copy(file.path(panel_dir, "growth_qa.csv"), "data/provenance/coi_reference_growth_qa.csv", overwrite = TRUE)) stop("Could not install acquired COI ledger and QA")
state <- coi_validate_retained_state(root, groups)
message("Installed ", nrow(state$records), " retained COI records across", length(groups), " arthropod target groups; historicalBivalvia inputs remain archived separately")
