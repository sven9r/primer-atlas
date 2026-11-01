#!/usr/bin/env Rscript
.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
source("R/functions.R")
source("R/data_layer.R")
source("R/fungal_explorer.R")
source("R/fungal_versions.R")
source("R/fungal_fit.R")
source("R/unite_reference.R")

# Replay the independently frozen R engine through the broad panel's sparse
# index. This checks all 884 pair counts and all 18,564 template identities.
old <- fungal_read_snapshot("v3-c6fa87f256e1", remote = FALSE)
meta <- old$taxonomy
meta$record_id <- meta$accession
meta$class_name <- meta$class
meta$order_name <- meta$order
old_fit <- fungal_fit_summary(old, rank = "class", pair_ids = old$pairs$pair_id[1L])
class_lookup <- old_fit$details[, c("accession", "group_id")]
meta$class_group_id <- class_lookup$group_id[match(meta$accession, class_lookup$accession)]
meta$raw_length_bp <- nchar(old$sequences[meta$accession])
coverage <- old_fit$summary[, c("rank", "taxon", "parent_path", "group_id", "taxon_resolved", "n_panel")]
names(coverage)[names(coverage) == "n_panel"] <- "n_selected"
coverage$n_source <- coverage$n_selected
sites <- old$sites[old$sites$status != "no_exact_binding_site", ]
sites$record_id <- sites$accession
data <- list(manifest = list(id = "unite-replay-fixture", n_source = 21L, n_selected = 21L),
    metadata = meta, coverage = coverage, catalog = old$catalog, pairs = old$pairs,
    sites = sites, sequences = old$sequences)
data$site_keys <- paste(sites$record_id, sites$catalog_id, sep = "\n")
for (pair_id in old$pairs$pair_id) {
  actual <- unite_reference_pair_fit(data, pair_id)
  expected <- old$pair_sequences[old$pair_sequences$pair_id == pair_id, ]
  expected <- expected[match(actual$details$record_id, expected$accession), ]
  stopifnot(identical(actual$details$amplicon_bp, as.numeric(expected$amplicon_bp)),
      sum(actual$summary$n_exact_products) == sum(expected$status == "scorable"),
      all(rowSums(actual$summary[c("n_exact_products", "n_single_sites", "n_no_exact_sites", "n_unresolved", "n_unavailable", "n_withheld")]) == actual$summary$n_selected))
}

# A genuine region alignment must contain precisely the extracted source
# bases. A lost base, wrong owner, or unequal row widths is rejected.
scratch <- tempfile("unite-alignment-")
dir.create(scratch)
file <- file.path(scratch, "its2.fasta")
writeLines(c(">A", "AC-GT", ">B", "ACGGT"), file)
fixture <- list(manifest = list(id = "unite-region-good"), directory = scratch,
    metadata = data.frame(record_id = c("A", "B", "C"), accession = c("AA", "BB", "CC"), class_group_id = "class-A"),
    sequences = c(A = "ACGT", B = "ACGGT", C = "TTTT"),
    regions = data.frame(record_id = c("A", "B"), region = "ITS2", start = 1L, end = c(4L, 5L), boundary_status = "complete"),
    alignment_index = data.frame(group_id = "class-A", region = "ITS2", file = "its2.fasta"))
aligned <- unite_reference_alignment(fixture, "class-A", "ITS2", limit = 1L)
stopifnot(length(aligned$alignment) == 1L, aligned$n_aligned == 2L, aligned$n_regions == 2L,
    aligned$n_total == 3L, aligned$n_missing_boundaries == 1L, aligned$n_columns == 5L,
    aligned$missing_metadata$record_id == "C")
fixture$manifest$id <- "unite-region-bad"
writeLines(c(">A", "AT-GT", ">B", "ACGGT"), file)
stopifnot(inherits(try(unite_reference_alignment(fixture, "class-A"), silent = TRUE), "try-error"))
stopifnot(inherits(try(.unite_safe_path(scratch, "../other"), silent = TRUE), "try-error"))
stopifnot(is.null(unite_reference_read(tempdir())))
message("UNITE sparse evidence matches all frozen R pairs; regional bases and coverage validated.")
