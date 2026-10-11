#!/usr/bin/env Rscript
.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
source("R/data_layer.R")
source("R/coi_retention.R")
source("R/coi_scoring.R")

root <- tempfile("coi_scoring_test_")
dir.create(root)
on.exit(unlink(root, recursive = TRUE), add = TRUE)
groups <- c("OrderA", "OrderB")
ref <- file.path(root, "data/reference/coi_reference_NC_001322.1.fasta")
dir.create(dirname(ref), recursive = TRUE)
writeLines(c(">NC_001322.1_COX1_reference", strrep("A", 1536)), ref)
for (group in groups) {
  raw <- file.path(root, "data/primerminer/raw", group, paste0(group, "_GB.fasta"))
  aligned <- file.path(root, "data/primerminer/full_aligned_reference", paste0(group, "_COX1_full_reference_aligned.fasta"))
  dir.create(dirname(raw), recursive = TRUE)
  dir.create(dirname(aligned), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(paste0(">", group, "001.1"), strrep("C", 600),
               paste0(">", group, "002.1"), strrep("G", 600)), raw)
  writeLines(c(readLines(ref), paste0(">", group, "001.1"), paste0(strrep("C", 600), strrep("-", 936)),
               paste0(">", group, "002.1"), paste0(strrep("G", 600), strrep("-", 936))), aligned)
  inputs <- coi_alignment_inputs(raw, ref, "test-only")
  inputs$output_sha256 <- atlas_sha256(aligned)
  jsonlite::write_json(inputs, paste0(aligned, ".inputs.json"), auto_unbox = TRUE)
}
ledger_path <- file.path(root, "data/provenance/coi_reference_accession_ledger.csv")
dir.create(dirname(ledger_path), recursive = TRUE)
fixture_ledger <- coi_retained_records(root, groups)
write.csv(fixture_ledger, ledger_path, row.names = FALSE)
selected <- coi_scoring_inputs(root, "full", groups, 2L, 4L)
stopifnot(length(selected$files) == 2L, sum(selected$manifest$n_templates) == 4L,
          selected$output_prefix == "full_order_", selected$evaluation_mode == "full")
expect_failure <- function(expression) stopifnot(inherits(tryCatch(force(expression), error = identity), "error"))
expect_failure(coi_scoring_inputs(root, "full", groups, 3L, 4L))
expect_failure(coi_scoring_inputs(root, "full", groups, 2L, 5L))
unlink(ledger_path)
expect_failure(coi_scoring_inputs(root, "full", groups, 2L, 4L))
write.csv(fixture_ledger, ledger_path, row.names = FALSE)

# The same record count does not prove that old identities survived.
alignment <- selected$files[1]
original <- readLines(alignment)
writeLines(sub("OrderA002.1", "OrderA999.1", original, fixed = TRUE), alignment)
expect_failure(coi_scoring_inputs(root, "full", groups, 2L, 4L))
writeLines(original, alignment)
raw <- file.path(root, "data/primerminer/raw/OrderA/OrderA_GB.fasta")
writeLines(sub(strrep("C", 600), strrep("T", 600), readLines(raw), fixed = TRUE), raw)
expect_failure(coi_scoring_inputs(root, "full", groups, 2L, 4L))

# A partially installed full panel cannot fall back to a different denominator.
unlink(selected$files[2])
expect_failure(coi_scoring_inputs(root, "auto", groups, 2L, 4L))
unlink(file.path(root, "data/primerminer/full_aligned_reference"), recursive = TRUE)
preview <- file.path(root, "data/primerminer/aligned_reference/OrderA_COX1_reference_aligned.fasta")
dir.create(dirname(preview), recursive = TRUE)
writeLines(original, preview)
legacy <- coi_scoring_inputs(root, "auto", groups)
stopifnot(legacy$evaluation_mode == "centroid_preview", legacy$output_prefix == "preview_order_",
          legacy$evidence_class == "legacy_centroid_preview")
message("COI scoring rejects stale, incomplete, or below-target full panels and labels previews explicitly.")
