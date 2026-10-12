#!/usr/bin/env Rscript
# Physical amplicon bp and frozen alignment columns are separate quantities.
.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
source("R/functions.R")
source("R/fungal_alignment.R")
source("R/fungal_alignment_ui.R")

find_tags <- function(node, class = NULL, name = NULL) {
  found <- list()
  visit <- function(item) {
    if (inherits(item, "shiny.tag")) {
      classes <- strsplit(if (is.null(item$attribs$class)) "" else item$attribs$class, " ", fixed = TRUE)[[1L]]
      if ((is.null(class) || class %in% classes) && (is.null(name) || item$name == name)) {
        found[[length(found) + 1L]] <<- item
      }
      lapply(item$children, visit)
    } else if (is.list(item)) lapply(item, visit)
    invisible(NULL)
  }
  visit(node)
  found
}
one_tag <- function(node, class) {
  tags <- find_tags(node, class)
  stopifnot(length(tags) == 1L)
  tags[[1L]]
}
near <- function(actual, expected) stopifnot(isTRUE(all.equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-10)))

raw <- c(A = "ACGTTGCAACGT", B = "TTGGACCAATGC")
gapped <- function(sequence) paste0("--", substring(sequence, 1, 2), "-",
  substring(sequence, 3, 8), "--", substring(sequence, 9), "-")
index <- fungal_alignment_index(c(A = gapped(raw[["A"]]), `_R_B` = gapped(reverse_complement(raw[["B"]]))), raw)
evidence <- data.frame(pair_id = "example", accession = c("A", "B"),
  forward_start = 2L, forward_end = 4L, reverse_start = 9L, reverse_end = 11L,
  amplicon_bp = 10L, status = "scorable", forward_status = "exact_compatible_site",
  reverse_status = "exact_compatible_site", stringsAsFactors = FALSE)
mapped <- fungal_alignment_pair_rows(pair_sequences = evidence, index = index)
pairs <- data.frame(pair_id = "example", pair_label = "ExampleF + ExampleR")
taxonomy <- data.frame(accession = c("A", "B"), organism = c("Fixture A", "Fixture B"), class = "Fixture class")
landmarks <- data.frame(accession = c("A", "B"), landmark_id = c("ITS1", "B_ONLY"),
  label = c("ITS1", "Other template feature"), start = 2L, end = 9L,
  source_note = "NCBI annotated complete forward feature")
regions <- fungal_alignment_landmark_rows(index, landmarks)
ceilings <- c(500L, 2000L, 3500L)
view <- function(reference = "A", domain = c(1L, 12L), rows = mapped, annotations = regions,
                 enabled = TRUE, alignment_index = index, pair_metadata = pairs, active = "example") {
  fungal_alignment_template_bp_view(alignment_index, rows, pair_metadata, taxonomy,
    reference_accession = reference, active_id = active, domain = domain, ceilings = ceilings,
    ns = function(value) paste0("fixture-", value), landmarks = annotations, show_regions = enabled)
}

# An inclusive 2..11 product occupies ten complete base cells. It is wider than
# 2..10, excludes all three alignment gaps, and keeps its full raw bounds.
forward_view <- view()
reverse_view <- view("B")
for(class_name in c("fungal-template-bp-window-regions","fungal-template-bp-accounting")) {
  secondary <- one_tag(forward_view,class_name)
  stopifnot(secondary$name=="details",is.null(secondary$attribs$open))
}
forward_product <- one_tag(forward_view, "fungal-template-bp-product")
reverse_product <- one_tag(reverse_view, "fungal-template-bp-product")
near(forward_product$attribs$x, 290 + 710 / 12)
near(forward_product$attribs$width, 10 * 710 / 12)
near(reverse_product$attribs$x, forward_product$attribs$x)
near(reverse_product$attribs$width, forward_product$attribs$width)
stopifnot(forward_view$attribs$`data-coordinate-system` == "raw_bases",
          forward_view$attribs$`data-template-bases` == 12L,
          forward_view$attribs$`data-start-base` == 1L, forward_view$attribs$`data-end-base` == 12L,
          one_tag(forward_view, "fungal-template-bp-pair")$attribs$`data-product-bp` == 10L,
          one_tag(forward_view, "fungal-template-bp-pair")$attribs$`data-start-base` == 2L,
          one_tag(forward_view, "fungal-template-bp-pair")$attribs$`data-end-base` == 11L)
for (rendered in list(forward_view, reverse_view)) {
  primers <- find_tags(rendered, "fungal-template-bp-primer")
  arrows <- find_tags(rendered, "fungal-template-bp-arrow")
  stopifnot(length(primers) == 2L, length(arrows) == 2L,
    identical(vapply(arrows, function(tag) tag$attribs$`data-primer-role`, character(1)), c("forward", "reverse")),
    identical(vapply(arrows, function(tag) tag$attribs$`data-direction`, integer(1)), c(1L, -1L)))
  near(vapply(primers, function(tag) tag$attribs$width, numeric(1)), rep(3 * 710 / 12, 2L))
  account <- one_tag(rendered, "fungal-template-bp-accounting-row")
  stopifnot(account$attribs$`data-product-bp` == 10L,
            account$attribs$`data-product-alignment-span` == 13L, account$attribs$`data-gap-columns` == 3L,
            account$attribs$`data-clipped` == "false", account$attribs$`data-visible-bp` == 10L,
            grepl("10 bp total · raw bases 2–11 (F 2–4; R 9–11) · 13 alignment columns including 3 gap columns", as.character(account), fixed = TRUE))
}
# MAFFT's reversed orientation reverses aligned arrows, but not raw-axis arrows.
stopifnot(mapped$forward_direction[2L] == -1L, mapped$reverse_direction[2L] == 1L)

# Genes use their own raw interval and preserve both source and aligned bounds.
# They must never borrow another template's annotation in either summary/track.
for (reference in c("A", "B")) {
  rendered <- view(reference)
  region <- one_tag(rendered, "fungal-template-bp-region")
  summary <- one_tag(rendered, "fungal-template-bp-window-region")
  rectangle <- find_tags(region, name = "rect")[[1L]]
  near(rectangle$attribs$x, 290 + 710 / 12)
  near(rectangle$attribs$width, 8 * 710 / 12)
  expected <- regions[regions$accession == reference, , drop = FALSE]
  stopifnot(region$attribs$`data-accession` == reference, summary$attribs$`data-accession` == reference,
            region$attribs$`data-start-base` == 2L, region$attribs$`data-end-base` == 9L,
            region$attribs$`data-start-col` == expected$left_col,
            region$attribs$`data-end-col` == expected$right_col,
            summary$attribs$`data-clipped` == "false",
            grepl("Raw bases 2–9 (8 bp)", summary$attribs$title, fixed = TRUE),
            grepl("NCBI annotated complete forward feature", summary$attribs$title, fixed = TRUE))
}
only_other <- view(annotations = regions[regions$accession == "B", , drop = FALSE])
hidden <- view(enabled = FALSE)
stopifnot(length(find_tags(only_other, "fungal-template-bp-region")) == 0L,
          length(find_tags(only_other, "fungal-template-bp-window-region")) == 0L,
          length(find_tags(hidden, "fungal-template-bp-region")) == 0L,
          length(find_tags(hidden, "fungal-template-bp-window-region")) == 0L,
          length(find_tags(hidden, "fungal-template-bp-product")) == 1L)
unresolved_regions <- regions
unresolved_regions$alignment_status <- "invalid_raw_coordinates"
stopifnot(length(find_tags(view(annotations = unresolved_regions), "fungal-template-bp-region")) == 0L)

# Clipping changes only what is drawn, not the full inclusive length/accounting.
clipped <- view(domain = c(3L, 7L))
near(one_tag(clipped, "fungal-template-bp-product")$attribs$width, 710)
near(find_tags(one_tag(clipped, "fungal-template-bp-region"), name = "rect")[[1L]]$attribs$width, 710)
stopifnot(one_tag(clipped, "fungal-template-bp-window-region")$attribs$`data-clipped` == "true",
          one_tag(clipped, "fungal-template-bp-pair")$attribs$`data-product-bp` == 10L,
          one_tag(clipped, "fungal-template-bp-product")$attribs$`data-visible-bp` == 5L,
          one_tag(clipped, "fungal-template-bp-product")$attribs$`data-clipped` == "true",
          grepl("10 bp total", as.character(clipped), fixed = TRUE),
          grepl("5 bp · partly in view", as.character(clipped), fixed = TRUE),
          grepl("5 bp in this window (partly in view)", as.character(clipped), fixed = TRUE),
          length(find_tags(clipped, "fungal-template-bp-primer")) == 1L)
near(one_tag(clipped, "fungal-template-bp-primer")$attribs$width, 2 * 710 / 5)
single_base <- view(domain = c(5L, 5L))
near(one_tag(single_base, "fungal-template-bp-product")$attribs$width, 710)
stopifnot(single_base$attribs$`data-start-base` == 5L, single_base$attribs$`data-end-base` == 5L)
outside <- view(domain = c(12L, 12L))
stopifnot(length(find_tags(outside, "fungal-template-bp-product")) == 0L,
          length(find_tags(outside, "fungal-template-bp-region")) == 0L,
          grepl("Outside this raw-base window", as.character(outside), fixed = TRUE),
          one_tag(outside, "fungal-template-bp-accounting-row")$attribs$`data-product-bp` == 10L,
          one_tag(outside, "fungal-template-bp-accounting-row")$attribs$`data-visible-bp` == 0L,
          grepl("0 bp in this window (outside window)", as.character(outside), fixed = TRUE))

# Exact individual components are useful without assigning a product or tier.
partial <- evidence[1L, , drop = FALSE]
partial$status <- "reverse_site_missing"
partial$reverse_start <- partial$reverse_end <- partial$amplicon_bp <- NA_integer_
partial$reverse_status <- "no_exact_binding_site"
partial_mapped <- fungal_alignment_pair_rows(pair_sequences = partial, index = index)
partial_view <- view(rows = partial_mapped)
stopifnot(length(find_tags(partial_view, "fungal-template-bp-product")) == 0L,
          length(find_tags(partial_view, "fungal-template-bp-primer")) == 1L,
          one_tag(partial_view, "fungal-template-bp-primer")$attribs$`data-primer-role` == "forward",
          is.null(one_tag(partial_view, "fungal-template-bp-accounting-row")$attribs$`data-product-bp`),
          grepl("Product unresolved", as.character(partial_view), fixed = TRUE))
ambiguous <- evidence[1L, , drop = FALSE]
ambiguous$status <- "low_specificity_site"
ambiguous$forward_status <- ambiguous$reverse_status <- "low_specificity_site"
ambiguous_mapped <- fungal_alignment_pair_rows(pair_sequences = ambiguous, index = index)
ambiguous_view <- view(rows = ambiguous_mapped)
stopifnot(length(find_tags(ambiguous_view, "fungal-template-bp-product")) == 0L,
          length(find_tags(ambiguous_view, "fungal-template-bp-primer")) == 0L,
          is.null(one_tag(ambiguous_view, "fungal-template-bp-accounting-row")$attribs$`data-gap-columns`))
duplicate_view <- view(rows = rbind(mapped, mapped))
stopifnot(length(find_tags(duplicate_view, "fungal-template-bp-product")) == 0L,
          grepl("Ambiguous pair evidence", as.character(duplicate_view), fixed = TRUE))

# A gap-only aligned window has no corresponding raw domain. It must stay empty
# instead of silently resetting to the whole template. Stale mapped rows cannot
# override an unverified template/index either.
gap_only <- view(domain = integer())
stopifnot(gap_only$attribs$class == "callout fungal-template-bp-empty-window",
          grepl("contains no bases from the chosen template", as.character(gap_only), fixed = TRUE),
          length(find_tags(gap_only, name = "svg")) == 0L)
stopifnot(view(domain = NULL)$attribs$`data-end-base` == 12L)
bad_raw <- raw
bad_raw[["A"]] <- paste0("T", substring(bad_raw[["A"]], 2L))
unverified <- view(alignment_index = fungal_alignment_index(c(A = gapped(raw[["A"]]),
  `_R_B` = gapped(reverse_complement(raw[["B"]]))), bad_raw))
stopifnot(length(find_tags(unverified, name = "svg")) == 0L,
          grepl("verified raw template is required", as.character(unverified), fixed = TRUE))

# Candidate limits retain the active pair and namespaced buttons. Empty filters
# do not invent a candidate or borrow one from the binding rows.
many_pairs <- data.frame(pair_id = paste0("p", 1:25), pair_label = paste("Pair", 1:25))
many_view <- view(pair_metadata = many_pairs, active = "p25")
buttons <- find_tags(many_view, "fungal-pair-choice")
stopifnot(length(buttons) == 20L, buttons[[1L]]$attribs$`data-pair` == "p25",
          buttons[[1L]]$attribs$`data-input` == "fixture-pair_click",
          length(find_tags(view(pair_metadata = pairs[FALSE, ], active = character()), "fungal-template-bp-accounting-row")) == 0L)

# In bp mode, the lower chart still uses a shared alignment axis. Its original
# column geometry is preserved, and the upper candidate chart is absent.
lower <- fungal_alignment_overview(index, mapped, pairs, taxonomy, "A", c("A", "B"),
  compare_ids = "example", active_id = "example", domain = c(1L, 18L), ceilings = ceilings,
  ns = identity, landmarks = regions, show_candidate_tracks = FALSE)
stopifnot(lower$attribs$`data-coordinate-system` == "alignment_columns",
          length(find_tags(lower, name = "svg")) == 1L,
          length(find_tags(lower, "fungal-pair-choice")) == 0L,
          length(find_tags(lower, "fungal-alignment-window-regions")) == 0L,
          grepl("Shared alignment columns", as.character(lower), fixed = TRUE),
          grepl("Alignment columns (including gaps)", as.character(lower), fixed = TRUE),
          !grepl("Candidate combinations on", as.character(lower), fixed = TRUE))
default_aligned <- fungal_alignment_overview(index, mapped, pairs, taxonomy, "A", c("A", "B"),
  compare_ids = "example", active_id = "example", domain = c(1L, 18L), ceilings = ceilings, ns = identity)
stopifnot(length(find_tags(default_aligned, name = "svg")) == 2L,
          length(find_tags(default_aligned, "fungal-pair-choice")) == 1L)

# The challenged real example has 335 physical bases and 1573 frozen columns.
# Count the actual aligned characters independently, then verify both reported
# quantities and the 335-base physical width; no snapshot files are rewritten.
versions <- jsonlite::read_json("data/pinned/ITS_FUNGAL/versions/index.json", simplifyVector = FALSE)$versions
directory <- file.path("data/pinned/ITS_FUNGAL/versions", versions[[2L]]$id)
paths <- file.path(directory, c("alignment.fasta", "sequences.fasta", "pair_sequences.csv", "landmarks.csv"))
hashes <- tools::md5sum(paths)
frozen_index <- fungal_alignment_index(parse_fasta(paths[1L]), parse_fasta(paths[2L]))
frozen_evidence <- read.csv(paths[3L], stringsAsFactors = FALSE)
example <- frozen_evidence[frozen_evidence$pair_id == "FUNGAL_ITS1_ITS2" &
  frozen_evidence$accession == "PZ962702.1", , drop = FALSE]
real_mapped <- fungal_alignment_pair_rows(pair_sequences = example, index = frozen_index)
stopifnot(nrow(real_mapped) == 1L, real_mapped$alignment_status == "mapped",
          real_mapped$forward_start == 1775L, real_mapped$reverse_end == 2109L,
          real_mapped$amplicon_bp == 335L, real_mapped$product_start_col == 7995L,
          real_mapped$product_end_col == 9567L, real_mapped$product_alignment_span == 1573L,
          real_mapped$n_gap_columns_in_product == 1238L)
characters <- strsplit(substr(frozen_index$alignment[["PZ962702.1"]], 7995L, 9567L), "", fixed = TRUE)[[1L]]
stopifnot(length(characters) == 1573L, sum(characters == "-") == 1238L,
          sum(characters != "-") == 335L, 2109L - 1775L + 1L == 335L)
real_regions <- fungal_alignment_landmark_rows(frozen_index, read.csv(paths[4L], stringsAsFactors = FALSE))
real_pairs <- read.csv(file.path(directory, "pairs.csv"), stringsAsFactors = FALSE)
real_pairs <- real_pairs[real_pairs$pair_id == example$pair_id, , drop = FALSE]
real_taxonomy <- read.csv(file.path(directory, "taxonomy.csv"), stringsAsFactors = FALSE)
real_view <- fungal_alignment_template_bp_view(frozen_index, real_mapped, real_pairs, real_taxonomy,
  "PZ962702.1", active_id = example$pair_id, domain = c(1774L, 2110L), ceilings = ceilings,
  ns = identity, landmarks = real_regions)
near(one_tag(real_view, "fungal-template-bp-product")$attribs$width, 335 * 710 / 337)
real_account <- one_tag(real_view, "fungal-template-bp-accounting-row")
stopifnot(real_account$attribs$`data-product-bp` == 335L,
          real_account$attribs$`data-product-alignment-span` == 1573L,
          real_account$attribs$`data-gap-columns` == 1238L,
          grepl("335 bp total · raw bases 1775–2109", as.character(real_account), fixed = TRUE),
          grepl("1,573 alignment columns including 1,238 gap columns", as.character(real_account), fixed = TRUE),
          identical(tools::md5sum(paths), hashes))
real_zoom <- fungal_alignment_template_bp_view(frozen_index, real_mapped, real_pairs, real_taxonomy,
  "PZ962702.1", active_id = example$pair_id, domain = c(2057L, 2231L), ceilings = ceilings,
  ns = identity, landmarks = real_regions)
zoom_product <- one_tag(real_zoom, "fungal-template-bp-product")
zoom_account <- one_tag(real_zoom, "fungal-template-bp-accounting-row")
near(zoom_product$attribs$width, 53 * 710 / 175)
stopifnot(zoom_product$attribs$`data-visible-bp` == 53L,
          zoom_product$attribs$`data-clipped` == "true",
          zoom_account$attribs$`data-visible-bp` == 53L, zoom_account$attribs$`data-product-bp` == 335L,
          zoom_account$attribs$`data-product-alignment-span` == 1573L,
          zoom_account$attribs$`data-gap-columns` == 1238L,
          grepl("335 bp total", as.character(real_zoom), fixed = TRUE),
          grepl("53 bp · partly in view", as.character(real_zoom), fixed = TRUE),
          grepl("53 bp in this window (partly in view)", as.character(zoom_account), fixed = TRUE),
          identical(tools::md5sum(paths), hashes))
message("Fungal template-bp view checks passed: inclusive geometry, raw strand arrows, exact partial sites, scoped genes, gap-only windows and independently counted 335-bp frozen example.")
