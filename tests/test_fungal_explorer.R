#!/usr/bin/env Rscript

source("R/functions.R")
source("R/fungal_explorer.R")

catalog <- fungal_catalog()
stopifnot(
  nrow(catalog) == 112L,
  all(c("source_direction", "direction_note", "catalog_id") %in% names(catalog)),
  !any(catalog$gene_locus %in% c("trnL", "rbcL", "IGS")),
  catalog$direction[catalog$primer_name == "ITS2"] == "reverse",
  catalog$source_direction[catalog$primer_name == "ITS2"] == "forward",
  !anyDuplicated(catalog$catalog_id)
)

reference <- parse_fasta("data/reference/its_fungal_reference_FN812768.2.fasta")[[1]]
panel <- fungal_evaluate_panel(catalog, setNames(reference, "FN812768.2"))
documented <- c("ITS1_ITS4", "ITS3_ITS4", "FITS7_ITS4", "GITS7_ITS4", "ITS86F_ITS4")
stopifnot(
  all(documented %in% panel$pairs$pair_id[panel$pairs$n_scorable == 1L]),
  all(panel$pairs$evidence_status[match(documented, panel$pairs$pair_id)] == "documented_pair"),
  sum(panel$pairs$n_scorable > 0L) >= 20L,
  panel$pairs$amplicon_median_bp[match("ITS1_ITS4", panel$pairs$pair_id)] == 737,
  panel$pairs$amplicon_median_bp[match("ITS3_ITS4", panel$pairs$pair_id)] == 391,
  any(panel$pairs$forward_primer == "ITS1" & panel$pairs$reverse_primer == "ITS2" & panel$pairs$n_scorable == 1L)
)

unresolved <- c("ITS1F_ITS4", "ITS1F_ITS4B", "ITS1F_LR21")
platform_rows <- fungal_platform_summary(panel$pair_sequences, panel$pairs)
stopifnot(
  all(unresolved %in% panel$pairs$pair_id),
  all(panel$pairs$n_scorable[match(unresolved, panel$pairs$pair_id)] == 0L),
  all(is.na(panel$pairs$amplicon_median_bp[match(unresolved, panel$pairs$pair_id)])),
  all(platform_rows$length_tier[match(unresolved, platform_rows$pair_id)] == "Unresolved"),
  all(platform_rows$n_length_eligible[match(unresolved, platform_rows$pair_id)] == 0L)
)

forward <- "AACCGGTTGACCTGACTGAC"
reverse <- "TCGTAGGCTACGATCGACTG"
small_catalog <- data.frame(
  primer_name = c("ExampleF", "ExampleR"), sequence = c(forward, reverse),
  direction = c("forward", "reverse"), gene_locus = c("ITS", "ITS"),
  catalog_id = c("F", "R"), target = c("Fungi", "Fungi"),
  primary_reference_status = c("publication_cited", "publication_cited"),
  stringsAsFactors = FALSE
)
make_template <- function(bp) {
  paste0(forward, strrep("A", bp - nchar(forward) - nchar(reverse)), reverse_complement(reverse))
}

# Same accession and inward nonoverlap are necessary; matching components on
# separate templates never provide pair support.
split_panel <- fungal_evaluate_panel(
  small_catalog, c(F_only = paste0(forward, strrep("A", 150)),
                   R_only = paste0(strrep("A", 150), reverse_complement(reverse)))
)
wrong_order <- fungal_evaluate_panel(
  small_catalog, c(reversed = paste0(reverse_complement(reverse), strrep("A", 150), forward))
)
stopifnot(nrow(split_panel$pairs) == 0L, nrow(wrong_order$pairs) == 0L)

# A second exact binding site stays ambiguous instead of being chosen to give
# a desired amplicon length.
ambiguous <- fungal_evaluate_panel(
  small_catalog,
  c(repeated = paste0(forward, strrep("A", 150), forward,
                     strrep("A", 150), reverse_complement(reverse)))
)
stopifnot(
  nrow(ambiguous$pairs) == 0L,
  ambiguous$sites$status[ambiguous$sites$direction == "forward"] == "ambiguous_binding_site",
  ambiguous$sites$n_matches[ambiguous$sites$direction == "forward"] == 2L
)

# Primer degeneracy is accepted, while ambiguity in the template binding site
# cannot establish exact compatibility.
degenerate_catalog <- small_catalog
degenerate_catalog$sequence[1] <- sub("A", "R", forward)
degenerate <- fungal_evaluate_panel(degenerate_catalog, c(exact = make_template(500)))
ambiguous_base <- fungal_evaluate_panel(
  small_catalog, c(uncertain = paste0("N", substring(make_template(500), 2)))
)
stopifnot(nrow(degenerate$pairs) == 1L, nrow(ambiguous_base$pairs) == 0L)

# Product length counts real bases, inclusive of primers, without alignment
# gaps. Metadata follows the sequence identifier rather than row order.
gapped <- fungal_evaluate_panel(
  small_catalog, c(B = paste0(forward, "--", strrep("A", 460), reverse_complement(reverse)),
                   A = make_template(501)),
  taxonomy = data.frame(accession = c("A", "B"), class = c("ClassA", "ClassB"),
                        phylum = c("PhylumA", "PhylumB"), organism = c("Species A", "Species B"))
)
stopifnot(
  identical(gapped$pair_sequences$amplicon_bp, c(500L, 501L)),
  identical(gapped$pair_sequences$class, c("ClassB", "ClassA")),
  gapped$pairs$n_panel == 2L, gapped$pairs$n_scorable == 2L
)
partial <- fungal_platform_summary(gapped$pair_sequences, gapped$pairs, "Illumina")
stopifnot(partial$n_length_eligible == 1L, partial$length_tier == "PacBio",
          partial$length_range_crosses_tiers,
          grepl("Illumina (1/2)", partial$platforms, fixed = TRUE))
one_class <- fungal_subset_summary(gapped$pair_sequences, gapped$pairs, "B")
no_class <- fungal_subset_summary(gapped$pair_sequences, gapped$pairs, character())
stopifnot(one_class$n_panel == 1L, one_class$n_scorable == 1L,
          one_class$amplicon_median_bp == 500,
          no_class$n_panel == 0L, no_class$n_scorable == 0L,
          is.na(no_class$amplicon_min_bp))
crossing_rows <- data.frame(pair_id = c("crossing", "crossing"), accession = c("C", "D"),
                            amplicon_bp = c(400, 600), status = "scorable")
crossing_pair <- data.frame(pair_id = "crossing", n_panel = 2L)
crossing <- fungal_platform_summary(crossing_rows, crossing_pair, "Illumina")
stopifnot(crossing$length_tier == "Illumina", crossing$length_range_crosses_tiers,
          crossing$length_tier_range == "Illumina → PacBio",
          crossing$length_eligible_fraction == 0.5, crossing$panel_eligible_fraction == 0.5)

# Inclusive user planning ceilings at every relevant boundary.
lengths <- c(452, 453, 500, 501, 2000, 2001, 3500, 3501)
threshold_pairs <- data.frame(pair_id = paste0("boundary_", lengths))
threshold_sequences <- data.frame(pair_id = threshold_pairs$pair_id,
                                  amplicon_bp = lengths, status = "scorable")
threshold_result <- fungal_platform_summary(threshold_sequences, threshold_pairs)
stopifnot(identical(threshold_result$length_tier,
                    c("Illumina", "Illumina", "Illumina", "PacBio", "PacBio", "Nanopore", "Nanopore", "Beyond configured limits")))
strict_illumina <- fungal_platform_summary(threshold_sequences, threshold_pairs,
                                           platform = "Illumina", illumina_max = 452)
stopifnot(identical(strict_illumina$n_length_eligible, c(1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L)))
pacbio <- fungal_platform_summary(threshold_sequences, threshold_pairs, "PacBio")
nanopore <- fungal_platform_summary(threshold_sequences, threshold_pairs, "Nanopore")
stopifnot(identical(pacbio$n_length_eligible, c(1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L)),
          identical(nanopore$n_length_eligible, c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L)))

# Old low-confidence coordinates cannot acquire platform colors just because
# a numeric product length is present.
untrusted_sequences <- data.frame(pair_id = "untrusted", amplicon_bp = 284,
                                   status = "reference_specific_low_confidence")
untrusted <- fungal_platform_summary(untrusted_sequences, data.frame(pair_id = "untrusted"))
stopifnot(untrusted$length_tier == "Unresolved", untrusted$n_length_eligible == 0L)

expect_error <- function(expr) inherits(tryCatch(force(expr), error = identity), "error")
stopifnot(
  expect_error(fungal_evaluate_panel(small_catalog, setNames(c(make_template(500), make_template(501)), c("same", "same")))),
  expect_error(fungal_platform_summary(threshold_sequences, threshold_pairs, illumina_max = 2100)),
  expect_error(fungal_platform_summary(threshold_sequences, threshold_pairs, "unknown"))
)

message("Fungal discovery, exact-site evidence, product lengths, and platform boundary checks passed.")
