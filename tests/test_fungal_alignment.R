#!/usr/bin/env Rscript

source("R/functions.R")
source("R/fungal_alignment.R")

expect_error <- function(expression, pattern) {
  value <- tryCatch(force(expression), error = identity)
  stopifnot(inherits(value, "error"), grepl(pattern, conditionMessage(value)))
}

raw <- c(A = "ACGTTGCAACGT", B = "TTGGACCAATGC")
gap_template <- function(sequence) paste0("--", substring(sequence, 1, 2), "-",
                                        substring(sequence, 3, 8), "--",
                                        substring(sequence, 9), "-")
alignment <- c(A = gap_template(raw[["A"]]),
               `_R_B` = gap_template(reverse_complement(raw[["B"]])))
index <- fungal_alignment_index(alignment, raw)
stopifnot(
  index$n_columns == 18L,
  identical(names(index$alignment), c("A", "B")),
  identical(index$metadata$status, c("verified", "verified")),
  identical(index$metadata$orientation, c("forward", "reverse")),
  identical(unname(index$maps$A), c(3L, 4L, 6L, 7L, 8L, 9L, 10L, 11L, 14L, 15L, 16L, 17L)),
  identical(unname(index$maps$B), rev(unname(index$maps$A))),
  sum(index$occupancy$n_bases[index$occupancy$accession == "A"]) == nchar(raw[["A"]]),
  identical(index$occupancy$start_col[index$occupancy$accession == "A"], c(3L, 6L, 14L)),
  identical(index$occupancy$end_col[index$occupancy$accession == "A"], c(4L, 11L, 17L))
)

evidence <- data.frame(pair_id = "example", accession = c("A", "B"),
                       forward_start = 2L, forward_end = 4L,
                       reverse_start = 9L, reverse_end = 11L,
                       amplicon_bp = 10L, status = "scorable",
                       forward_status = "exact_compatible_site",
                       reverse_status = "exact_compatible_site", stringsAsFactors = FALSE)
mapped <- fungal_alignment_pair_rows(pair_sequences = evidence, index = index)
stopifnot(
  identical(mapped$alignment_status, c("mapped", "mapped")),
  identical(mapped$forward_start_col, c(4L, 16L)),
  identical(mapped$forward_end_col, c(7L, 14L)),
  identical(mapped$reverse_start_col, c(14L, 7L)),
  identical(mapped$reverse_end_col, c(16L, 4L)),
  identical(mapped$forward_direction, c(1L, -1L)),
  identical(mapped$reverse_direction, c(-1L, 1L)),
  all(mapped$product_start_col == 4L), all(mapped$product_end_col == 16L),
  all(mapped$amplicon_bp == 10L), all(mapped$product_alignment_span == 13L),
  all(mapped$n_gap_columns_in_product == 3L),
  identical(fungal_alignment_pair_rows(index, evidence), mapped),
  nrow(fungal_alignment_pair_rows(pair_sequences = evidence, pair_ids = character(), index = index)) == 0L,
  nrow(fungal_alignment_pair_rows(pair_sequences = evidence, accessions = character(), index = index)) == 0L,
  nrow(fungal_alignment_pair_rows(pair_sequences = evidence, accessions = "B", index = index)) == 1L
)

# A verified exact component remains visible when the other is missing, but
# ambiguous/low-specificity coordinates and unresolved product spans do not.
partial <- evidence[1, ]
partial$status <- "reverse_site_missing"
partial$reverse_start <- partial$reverse_end <- partial$amplicon_bp <- NA_integer_
partial$reverse_status <- "no_exact_binding_site"
partial_mapped <- fungal_alignment_pair_rows(pair_sequences = partial, index = index)
stopifnot(partial_mapped$alignment_status == "partial_exact_sites",
          partial_mapped$forward_alignment_status == "mapped",
          partial_mapped$forward_start_col == 4L,
          is.na(partial_mapped$reverse_start_col), is.na(partial_mapped$product_start_col))
low <- evidence[1, ]
low$status <- "low_specificity_site"
low$forward_status <- low$reverse_status <- "low_specificity_site"
low_mapped <- fungal_alignment_pair_rows(pair_sequences = low, index = index)
stopifnot(low_mapped$alignment_status == "unresolved_binding_evidence",
          is.na(low_mapped$forward_start_col), is.na(low_mapped$product_start_col))
missing_status <- evidence[1, ]
missing_status$status <- NA_character_
missing_status$forward_status <- missing_status$reverse_status <- "no_exact_binding_site"
stopifnot(fungal_alignment_pair_rows(pair_sequences = missing_status, index = index)$alignment_status == "unresolved_binding_evidence")

# Alignment evidence and raw products must agree. No product is drawn from a
# missing template, mismatched sequence, outward sites, or inconsistent length.
no_raw <- fungal_alignment_index(alignment)
stopifnot(all(fungal_alignment_pair_rows(pair_sequences = evidence, index = no_raw)$alignment_status == "raw_template_unavailable"))
wrong_raw <- raw
wrong_raw[["A"]] <- paste0("T", substring(raw[["A"]], 2))
mismatch <- fungal_alignment_index(alignment, wrong_raw)
stopifnot(fungal_alignment_pair_rows(pair_sequences = evidence[1, ], index = mismatch)$alignment_status == "alignment_template_mismatch")
unknown <- evidence[1, ]
unknown$accession <- "missing"
stopifnot(fungal_alignment_pair_rows(pair_sequences = unknown, index = index)$alignment_status == "template_not_in_alignment")
outward <- evidence[1, ]
outward$reverse_start <- 3L
stopifnot(fungal_alignment_pair_rows(pair_sequences = outward, index = index)$alignment_status == "invalid_raw_coordinates")
invalid <- evidence[1, ]
invalid$forward_start <- 0L
invalid_mapped <- fungal_alignment_pair_rows(pair_sequences = invalid, index = index)
stopifnot(invalid_mapped$forward_alignment_status == "invalid_raw_coordinates",
          is.na(invalid_mapped$product_start_col))
length_mismatch <- evidence[1, ]
length_mismatch$amplicon_bp <- 11L
stopifnot(fungal_alignment_pair_rows(pair_sequences = length_mismatch, index = index)$alignment_status == "raw_product_length_mismatch")
inconsistent <- evidence[1, ]
inconsistent$forward_status <- "ambiguous_binding_site"
stopifnot(fungal_alignment_pair_rows(pair_sequences = inconsistent, index = index)$alignment_status == "invalid_binding_site_evidence")
expect_error(fungal_alignment_index(c(A = "ACG", B = "AC"), raw), "different column")
expect_error(fungal_alignment_index(c(A = "ACG", `_R_A` = "ACG"), raw), "repeats an accession")
expect_error(fungal_alignment_index(c("ACG")), "named character")

# A nucleotide window reports the actual frozen columns, including gaps.
# Overlay bounds stay in the global alignment coordinate system.
window <- fungal_alignment_window(index, mapped, start = 5L, max_columns = 5L)
stopifnot(identical(window$columns, 5:9),
          identical(unname(window$rows$sequence), c("-GTTG", "-ATTG")),
          identical(window$bases["A", ], setNames(c("-", "G", "T", "T", "G"), as.character(5:9))),
          all(window$overlays$product_visible),
          window$overlays$forward_visible[1], !window$overlays$reverse_visible[1],
          window$overlays$forward_window_start_col[1] == 5L,
          window$overlays$forward_window_end_col[1] == 7L,
          fungal_alignment_window(index, start = 100L, max_columns = 5L)$start_col == 14L,
          fungal_alignment_window(index, start = -10L, max_columns = 5L)$start_col == 1L,
          length(fungal_alignment_window(index, max_columns = 100L)$columns) == 18L)
partial_window <- fungal_alignment_window(index, partial_mapped, start = 4L, max_columns = 5L)
stopifnot(partial_window$overlays$forward_visible,
          !partial_window$overlays$reverse_visible, !partial_window$overlays$product_visible)
empty <- fungal_alignment_index(setNames(character(), character()), setNames(character(), character()))
stopifnot(empty$n_columns == 0L, nrow(empty$metadata) == 0L,
          ncol(fungal_alignment_window(empty)$bases) == 0L,
          nrow(fungal_alignment_window(index, accessions = character())$bases) == 0L)

# Validate the shipped frozen panel, rather than replacing its alignment with
# a new alignment or mapping a consensus/median position onto every template.
version_index <- jsonlite::read_json("data/pinned/ITS_FUNGAL/versions/index.json", simplifyVector = FALSE)
version_directory <- file.path("data/pinned/ITS_FUNGAL/versions", tail(version_index$versions, 1)[[1]]$id)
frozen_alignment <- parse_fasta(file.path(version_directory, "alignment.fasta"))
frozen_sequences <- parse_fasta(file.path(version_directory, "sequences.fasta"))
frozen_pairs <- read.csv(file.path(version_directory, "pair_sequences.csv"), stringsAsFactors = FALSE)
frozen_index <- fungal_alignment_index(frozen_alignment, frozen_sequences)
frozen_mapped <- fungal_alignment_pair_rows(pair_sequences = frozen_pairs, index = frozen_index)
supported <- frozen_mapped[frozen_mapped$status == "scorable", , drop = FALSE]
stopifnot(all(frozen_index$metadata$status == "verified"),
          nrow(frozen_mapped) == nrow(frozen_pairs),
          all(supported$alignment_status == "mapped"),
          all(supported$product_alignment_span >= supported$amplicon_bp),
          all(is.na(frozen_mapped$product_start_col[frozen_mapped$status != "scorable"])))
for (accession in frozen_index$metadata$accession) {
  columns <- frozen_index$maps[[accession]]
  projected <- paste0(strsplit(frozen_index$alignment[[accession]], "", fixed = TRUE)[[1]][columns], collapse = "")
  expected <- frozen_sequences[[accession]]
  if (frozen_index$metadata$strand[match(accession, frozen_index$metadata$accession)] == -1L) {
    # Descending columns retrieve the complemented raw sequence, without
    # reversing its raw-coordinate order a second time.
    projected <- reverse_complement(paste0(rev(strsplit(projected, "", fixed = TRUE)[[1]]), collapse = ""))
  }
  stopifnot(identical(projected, expected))
}

# The renderer uses these mappings and physical lengths without inventing an
# unresolved product, reversing primer arrows, or coloring alignment gaps as
# nucleotides inside a primer or product.
source("R/fungal_alignment_ui.R")
ceilings <- c(500L, 2000L, 3500L)
fixture_pairs <- data.frame(pair_id = "example", pair_label = "ExampleF + ExampleR", stringsAsFactors = FALSE)
fixture_taxonomy <- data.frame(accession = names(raw), organism = c("Fixture A", "Fixture B"),
                               class = "Fixture class", stringsAsFactors = FALSE)
scale <- function(column) column * 10L
span <- as.character(.fungal_alignment_ui_span(mapped[1, ], c(1L, 18L), scale, 40L, ceilings, "ExampleF + ExampleR"))
partial_span <- as.character(.fungal_alignment_ui_span(partial_mapped, c(1L, 18L), scale, 40L, ceilings))
reverse_span <- as.character(.fungal_alignment_ui_span(mapped[2, ], c(1L, 18L), scale, 40L, ceilings))
stopifnot(grepl('data-product-bp="10"', span, fixed = TRUE),
          grepl('data-start-col="4"', span, fixed = TRUE),
          grepl('data-end-col="16"', span, fixed = TRUE),
          grepl('class="fungal-alignment-product"', span, fixed = TRUE),
          grepl('data-alignment-status="partial_exact_sites"', partial_span, fixed = TRUE),
          grepl('data-primer-role="forward"', partial_span, fixed = TRUE),
          !grepl('data-primer-role="reverse"', partial_span, fixed = TRUE),
          !grepl('data-product-bp=', partial_span, fixed = TRUE),
          !grepl('class="fungal-alignment-product"', partial_span, fixed = TRUE),
          grepl('data-primer-role="forward" data-direction="-1"', reverse_span, fixed = TRUE),
          grepl('data-primer-role="reverse" data-direction="1"', reverse_span, fixed = TRUE))
fixture_bases <- as.character(fungal_alignment_bases(index, mapped, fixture_pairs, fixture_taxonomy,
                                                     accessions = "A", start = 4L, width = 5L,
                                                     active_id = "example", ceilings = ceilings))
gap_cell <- regmatches(fixture_bases, regexpr('<td[^>]*data-alignment-col="5"[^>]*>[^<]*</td>', fixture_bases))
stopifnot(grepl('data-template-count="1"', fixture_bases, fixed = TRUE),
          grepl('data-raw-position="3"', fixture_bases, fixed = TRUE),
          grepl('class="fungal-alignment-base gap"', gap_cell, fixed = TRUE),
          grepl('>-</td>', gap_cell, fixed = TRUE),
          !grepl('data-raw-position=', gap_cell, fixed = TRUE),
          !grepl('style=', gap_cell, fixed = TRUE))
partial_bases <- as.character(fungal_alignment_bases(index, partial_mapped, fixture_pairs, fixture_taxonomy,
                                                     accessions = "A", start = 4L, width = 5L,
                                                     active_id = "example", ceilings = ceilings))
stopifnot(grepl('class="fungal-alignment-base forward"', partial_bases, fixed = TRUE),
          !grepl('class="fungal-alignment-base product', partial_bases, fixed = TRUE),
          !grepl('class="fungal-alignment-base reverse', partial_bases, fixed = TRUE))
frozen_pair_metadata <- read.csv(file.path(version_directory, "pairs.csv"), stringsAsFactors = FALSE)
frozen_taxonomy <- read.csv(file.path(version_directory, "taxonomy.csv"), stringsAsFactors = FALSE)
chosen <- frozen_pair_metadata$pair_id[which.max(frozen_pair_metadata$n_scorable)]
overview <- as.character(fungal_alignment_overview(frozen_index, frozen_mapped, head(frozen_pair_metadata, 20L),
                         frozen_taxonomy, names(frozen_sequences)[1], names(frozen_sequences),
                         compare_ids = chosen, active_id = chosen, domain = c(1L, frozen_index$n_columns),
                         ceilings = ceilings, ns = function(value) paste0("fixture-", value)))
full_bases <- as.character(fungal_alignment_bases(frozen_index, frozen_mapped, frozen_pair_metadata,
                           frozen_taxonomy, names(frozen_sequences), start = 100L, width = 120L,
                           active_id = chosen, ceilings = ceilings))
stopifnot(grepl(sprintf('data-template-count="%d"', length(frozen_sequences)), overview, fixed = TRUE),
          grepl('data-coordinate-system="alignment_columns"', overview, fixed = TRUE),
          grepl('data-start-col="100"', full_bases, fixed = TRUE),
          grepl('data-end-col="219"', full_bases, fixed = TRUE),
          grepl(sprintf('data-template-count="%d"', length(frozen_sequences)), full_bases, fixed = TRUE))
cat(sprintf("Fungal alignment checks passed: %d frozen templates, %d columns, %d exact product mappings.\n",
            length(frozen_sequences), frozen_index$n_columns, nrow(supported)))
