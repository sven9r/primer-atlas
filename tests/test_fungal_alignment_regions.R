#!/usr/bin/env Rscript
# Frozen gene/region annotations must follow their own raw-template bases.
# Their alignment spans can contain gaps without changing the physical length.
.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
source("R/functions.R")
source("R/fungal_alignment.R")

expect_error <- function(expression, pattern) {
  value <- tryCatch(force(expression), error = identity)
  stopifnot(inherits(value, "error"), grepl(pattern, conditionMessage(value)))
}

raw <- c(A = "ACGTTGCAACGT", B = "TTGGACCAATGC")
gapped <- function(sequence) paste0("--", substring(sequence, 1, 2), "-",
  substring(sequence, 3, 8), "--", substring(sequence, 9), "-")
alignment <- c(A = gapped(raw[["A"]]), `_R_B` = gapped(reverse_complement(raw[["B"]])))
index <- fungal_alignment_index(alignment, raw)
landmarks <- data.frame(accession = c("A", "B"), landmark_id = "ITS1",
  label = "ITS1", start = 2L, end = 9L,
  source_note = "NCBI annotated complete forward feature", stringsAsFactors = FALSE)
mapped <- fungal_alignment_landmark_rows(index, landmarks)
stopifnot(
  identical(mapped$alignment_status, c("mapped", "mapped")),
  identical(mapped$start_col, c(4L, 16L)), identical(mapped$end_col, c(14L, 7L)),
  identical(mapped$left_col, c(4L, 7L)), identical(mapped$right_col, c(14L, 16L)),
  identical(mapped$direction, c(1L, -1L)), all(mapped$raw_bp == 8L),
  identical(mapped$right_col - mapped$left_col + 1L, c(11L, 10L)),
  identical(mapped$source_note, landmarks$source_note),
  identical(mapped$evidence_status, rep("annotated_reference", 2L)))

# First/last raw bases map to non-gap cells, excluding terminal padding.
ends <- landmarks
ends$start <- 1L
ends$end <- 12L
mapped_ends <- fungal_alignment_landmark_rows(index, ends)
stopifnot(identical(mapped_ends$start_col, c(3L, 17L)),
          identical(mapped_ends$end_col, c(17L, 3L)), all(mapped_ends$raw_bp == 12L))

# A reference-scoped historical table belongs only to the manifest reference.
# The currently displayed different template cannot inherit its architecture.
unscoped <- landmarks[1L, setdiff(names(landmarks), "accession"), drop = FALSE]
unscoped$source_note <- "A feature architecture; boundaries curated for the pilot display"
unscoped$landmark_type <- "conserved_anchor"
pilot <- fungal_alignment_landmark_rows(index, unscoped, reference_accession = "A")
stopifnot(nrow(pilot) == 1L, pilot$accession == "A", pilot$alignment_status == "mapped",
          pilot$evidence_status == "curated_reference", pilot$landmark_type == "conserved_anchor",
          !any(pilot$accession == "B"), pilot$start_col == 4L)
no_scope <- fungal_alignment_landmark_rows(index, unscoped)
stopifnot(nrow(no_scope) == 1L, no_scope$alignment_status != "mapped",
          is.na(no_scope$start_col), is.na(no_scope$end_col))
explicit_missing <- landmarks[1L, , drop = FALSE]
explicit_missing$accession <- NA_character_
stopifnot(is.na(fungal_alignment_landmark_rows(index, explicit_missing, "A")$start_col))

# A valid interval still needs exact raw-template and alignment agreement.
no_raw <- fungal_alignment_index(alignment)
stopifnot(all(fungal_alignment_landmark_rows(no_raw, landmarks)$alignment_status == "raw_template_unavailable"))
wrong_raw <- raw
wrong_raw[["A"]] <- paste0("T", substring(wrong_raw[["A"]], 2L))
wrong_index <- fungal_alignment_index(alignment, wrong_raw)
wrong <- fungal_alignment_landmark_rows(wrong_index, landmarks)
stopifnot(wrong$alignment_status[1L] == "alignment_template_mismatch", is.na(wrong$start_col[1L]),
          wrong$alignment_status[2L] == "mapped")
outside <- landmarks[1L, , drop = FALSE]
outside$accession <- "A.2"
stopifnot(fungal_alignment_landmark_rows(index, outside)$alignment_status == "template_not_in_alignment")

# Invalid annotations remain visible as unresolved rows, with no coordinates.
invalid <- landmarks[rep(1L, 6L), , drop = FALSE]
invalid$start <- c(0, 2.5, NA, Inf, 9, 2)
invalid$end <- c(9, 9, 9, 9, 2, 13)
invalid_rows <- fungal_alignment_landmark_rows(index, invalid)
stopifnot(nrow(invalid_rows) == 6L, all(invalid_rows$alignment_status == "invalid_raw_coordinates"),
          all(is.na(invalid_rows$start_col)), all(is.na(invalid_rows$left_col)), all(is.na(invalid_rows$raw_bp)))
numeric_text <- landmarks
numeric_text$start <- as.character(numeric_text$start)
numeric_text$end <- as.character(numeric_text$end)
stopifnot(identical(fungal_alignment_landmark_rows(index, numeric_text)$start_col, mapped$start_col))
empty <- fungal_alignment_landmark_rows(index, landmarks[FALSE, , drop = FALSE])
stopifnot(nrow(empty) == 0L, all(c("alignment_status", "start_col", "end_col", "raw_bp") %in% names(empty)),
          nrow(fungal_alignment_landmark_rows(index, NULL)) == 0L)
expect_error(fungal_alignment_landmark_rows(list(), landmarks), "index must")
expect_error(fungal_alignment_landmark_rows(index, landmarks[setdiff(names(landmarks), "start")]), "missing raw")

# The shipped V1 and V2 keep distinct immutable annotation populations. V1's
# five curated regions/anchors do not silently augment V2's nine NCBI features.
versions_root <- "data/pinned/ITS_FUNGAL/versions"
versions <- jsonlite::read_json(file.path(versions_root, "index.json"), simplifyVector = FALSE)$versions
stopifnot(length(versions) >= 2L)
load_version <- function(version) {
  directory <- file.path(versions_root, version$id)
  paths <- file.path(directory, c("manifest.json", "sequences.fasta", "alignment.fasta", "landmarks.csv"))
  hashes <- tools::md5sum(paths)
  manifest <- jsonlite::read_json(paths[1L], simplifyVector = FALSE)
  sequences <- parse_fasta(paths[2L])
  frozen_index <- fungal_alignment_index(parse_fasta(paths[3L]), sequences)
  frozen_landmarks <- read.csv(paths[4L], stringsAsFactors = FALSE)
  rows <- fungal_alignment_landmark_rows(frozen_index, frozen_landmarks, manifest$reference_accession)
  stopifnot(identical(tools::md5sum(paths), hashes))
  list(index = frozen_index, landmarks = frozen_landmarks, rows = rows)
}
v1 <- load_version(versions[[1L]])
v2 <- load_version(versions[[2L]])
stopifnot(nrow(v1$rows) == 5L, all(v1$rows$accession == "FN812768.2"),
          all(v1$rows$alignment_status == "mapped"), all(v1$rows$evidence_status == "curated_reference"),
          identical(v1$rows$start_col, as.integer(v1$landmarks$start)),
          sum(v1$rows$landmark_type == "conserved_anchor") == 2L,
          nrow(v2$rows) == 9L, all(v2$rows$alignment_status == "mapped"),
          all(v2$rows$evidence_status == "annotated_reference"),
          identical(v2$rows$accession, v2$landmarks$accession),
          !"FN812768.2" %in% v2$rows$accession,
          setequal(unique(v2$rows$accession), c("PZ962702.1", "PZ276773.1", "FJ430740.1")))
for (i in seq_len(nrow(v2$rows))) {
  row <- v2$rows[i, , drop = FALSE]
  map <- v2$index$maps[[row$accession]]
  stopifnot(row$start_col == unname(map[row$start]), row$end_col == unname(map[row$end]),
            row$raw_bp == row$end - row$start + 1L)
}
# Independently inspected V2 annotations make a strong real-world regression:
# FJ430740's 168 raw 5.8S bases span 918 frozen columns, including 750 gaps.
# Other templates' 5.8S endpoints must not replace these with a common band.
gene <- v2$rows[v2$rows$landmark_id == "5_8S", , drop = FALSE]
gene <- gene[match(c("PZ962702.1", "PZ276773.1", "FJ430740.1"), gene$accession), , drop = FALSE]
stopifnot(identical(gene$left_col, c(9518L, 9519L, 8717L)),
          identical(gene$right_col, c(9688L, 9672L, 9634L)),
          identical(gene$raw_bp, c(158L, 153L, 168L)),
          gene$right_col[3L] - gene$left_col[3L] + 1L == 918L,
          gene$right_col[3L] - gene$left_col[3L] + 1L - gene$raw_bp[3L] == 750L)

# Rendering clips the displayed span to the visible alignment window while
# keeping full raw and mapped coordinates in its metadata and source tooltip.
source("R/fungal_alignment_ui.R")
scale <- function(column) column * 10
clipped <- as.character(.fungal_alignment_ui_regions(mapped, "A", c(5L, 8L), scale, 40L))
reverse_clipped <- as.character(.fungal_alignment_ui_regions(mapped, "B", c(8L, 10L), scale, 40L))
stopifnot(grepl('data-accession="A"', clipped, fixed = TRUE),
          grepl('data-start-col="4"', clipped, fixed = TRUE),
          grepl('data-end-col="14"', clipped, fixed = TRUE),
          grepl('<rect x="50"', clipped, fixed = TRUE),
          grepl('width="30"', clipped, fixed = TRUE),
          grepl("(8 bp)", clipped, fixed = TRUE),
          grepl('data-start-col="7"', reverse_clipped, fixed = TRUE),
          grepl('data-end-col="16"', reverse_clipped, fixed = TRUE),
          grepl('<rect x="80"', reverse_clipped, fixed = TRUE),
          grepl('width="20"', reverse_clipped, fixed = TRUE))
for (absent in list(
    .fungal_alignment_ui_regions(mapped, "A", c(15L, 18L), scale, 40L),
    .fungal_alignment_ui_regions(mapped, "C", c(1L, 18L), scale, 40L),
    .fungal_alignment_ui_regions(invalid_rows, "A", c(1L, 18L), scale, 40L))) {
  html <- as.character(absent)
  stopifnot(grepl("No frozen gene / region boundary", html, fixed = TRUE),
            !grepl('class="fungal-alignment-region"', html, fixed = TRUE))
}
pilot_html <- as.character(.fungal_alignment_ui_regions(v1$rows, "FN812768.2", c(1L, 739L), scale, 40L))
stopifnot(grepl("18S anchor (partial)", pilot_html, fixed = TRUE),
          grepl("28S anchor (partial)", pilot_html, fixed = TRUE),
          grepl('stroke-dasharray="3 2"', pilot_html, fixed = TRUE),
          grepl("Curated pilot boundaries", pilot_html, fixed = TRUE),
          !grepl("Frozen NCBI annotation", pilot_html, fixed = TRUE))

# Region header colspans include alignment gaps, but the nucleotide tooltip
# labels only actual annotated raw bases, never a gap or another template.
cells <- .fungal_alignment_ui_region_cells(mapped, "A", 2:16)
stopifnot(length(cells) == 3L,
          identical(vapply(cells, function(cell) as.integer(cell$attribs$colspan), integer(1)), c(2L, 11L, 2L)),
          identical(vapply(cells, function(cell) as.integer(cell$attribs$`data-start-col`), integer(1)), c(2L, 4L, 15L)),
          identical(vapply(cells, function(cell) as.integer(cell$attribs$`data-end-col`), integer(1)), c(3L, 14L, 16L)))
empty_pairs <- data.frame(pair_id = character(), pair_label = character())
empty_binding <- data.frame(pair_id = character(), accession = character())
taxonomy <- data.frame(accession = c("A", "B"), organism = c("Fixture A", "Fixture B"), class = "Fixture class")
render_bases <- function(annotations, enabled = TRUE) as.character(fungal_alignment_bases(
  index, empty_binding, empty_pairs, taxonomy, accessions = c("A", "B"),
  start = 4L, width = 5L, active_id = character(), ceilings = c(500, 2000, 3500),
  landmarks = annotations, show_regions = enabled))
base_html <- render_bases(mapped[mapped$accession == "A", , drop = FALSE])
region_rows <- regmatches(base_html, gregexpr('<tr class="fungal-alignment-region-row"[^>]*>', base_html))[[1]]
stopifnot(length(region_rows) == 1L, grepl('data-accession="A"', region_rows, fixed = TRUE))
gap_cell <- regmatches(base_html, regexpr('<td[^>]*data-alignment-col="5"[^>]*>[^<]*</td>', base_html))
base_cell <- regmatches(base_html, regexpr('<td[^>]*data-alignment-col="4"[^>]*>[^<]*</td>', base_html))
stopifnot(grepl("gap; no raw base", gap_cell, fixed = TRUE), !grepl("ITS1", gap_cell, fixed = TRUE),
          !grepl('data-raw-position=', gap_cell, fixed = TRUE),
          grepl("ITS1", base_cell, fixed = TRUE), grepl('data-raw-position="2"', base_cell, fixed = TRUE))
hidden <- render_bases(mapped, FALSE)
stopifnot(!grepl('class="fungal-alignment-region-row"', hidden, fixed = TRUE),
          !grepl(" · ITS1", hidden, fixed = TRUE), grepl('data-alignment-col="4"', hidden, fixed = TRUE))
overview <- function(enabled, annotations = mapped, domain = c(1L, 18L), reference = "A") as.character(fungal_alignment_overview(
  index, empty_binding, empty_pairs, taxonomy, reference_accession = reference,
  selected_accessions = c("A", "B"), compare_ids = character(), active_id = character(),
  domain = domain, ceilings = c(500, 2000, 3500), ns = identity,
  landmarks = annotations, show_regions = enabled))
stopifnot(grepl('class="fungal-alignment-region"', overview(TRUE), fixed = TRUE),
          !grepl('class="fungal-alignment-region"', overview(FALSE), fixed = TRUE))

# Compact summaries remain outside the scrolling chart. Even when another
# selected template has an overlapping feature, the reference's summary may
# contain only its own mapped annotations that overlap the displayed window.
summary_fixture <- mapped
summary_fixture$landmark_id[summary_fixture$accession == "B"] <- "B_ONLY"
summary_fixture$label[summary_fixture$accession == "B"] <- "Other template feature"
outside_landmark <- landmarks[1L, , drop = FALSE]
outside_landmark$landmark_id <- "A_OUTSIDE"
outside_landmark$label <- "Outside feature"
outside_landmark$start <- 10L
outside_landmark$end <- 12L
outside_mapped <- fungal_alignment_landmark_rows(index, outside_landmark)
unresolved_landmark <- invalid_rows[1L, , drop = FALSE]
unresolved_landmark$landmark_id <- "A_UNRESOLVED"
unresolved_landmark$label <- "Unresolved feature"
# Stale coordinates must not override an explicitly unresolved mapping status.
unresolved_landmark$left_col <- 5L
unresolved_landmark$right_col <- 7L
summary_fixture <- rbind(summary_fixture, outside_mapped, unresolved_landmark)
summary_tags <- function(html) regmatches(html, gregexpr(
  '<span class="fungal-alignment-window-region [^"]*"[^>]*>', html, perl = TRUE))[[1L]]
summary_html <- overview(TRUE, summary_fixture, c(5L, 8L))
secondary_summary <- regmatches(summary_html,regexpr('<details class="fungal-alignment-window-regions"[^>]*>',summary_html))
stopifnot(length(secondary_summary)==1L,nzchar(secondary_summary),!grepl(" open",secondary_summary,fixed=TRUE),
  grepl("Region coordinates and sources",summary_html,fixed=TRUE))
summary <- summary_tags(summary_html)
stopifnot(length(summary) == 1L,
          grepl('data-accession="A"', summary, fixed = TRUE),
          grepl('data-landmark-id="ITS1"', summary, fixed = TRUE),
          grepl('data-start-col="4"', summary, fixed = TRUE),
          grepl('data-end-col="14"', summary, fixed = TRUE),
          grepl('data-clipped="true"', summary, fixed = TRUE),
          grepl("Raw bases 2–9 (8 bp)", summary, fixed = TRUE),
          grepl("NCBI annotated complete forward feature", summary, fixed = TRUE),
          grepl("columns 4–14 · partly in view", summary_html, fixed = TRUE),
          regexpr('class="fungal-alignment-window-regions"', summary_html, fixed = TRUE)[1L] <
            regexpr('class="fungal-alignment-scroll"', summary_html, fixed = TRUE)[1L],
          grepl('class="tiny fungal-alignment-pan-hint"', summary_html, fixed = TRUE))
# Exact boundary overlap is included, and a fully visible feature is not marked
# clipped. Reversed raw orientation still reports the ordered alignment bounds.
complete_summary <- summary_tags(overview(TRUE, summary_fixture, c(4L, 14L)))
boundary_summary <- summary_tags(overview(TRUE, summary_fixture, c(1L, 4L)))
reverse_summary <- summary_tags(overview(TRUE, summary_fixture, c(8L, 10L), "B"))
stopifnot(length(complete_summary) == 1L,
          grepl('data-clipped="false"', complete_summary, fixed = TRUE),
          length(boundary_summary) == 1L,
          grepl('data-landmark-id="ITS1"', boundary_summary, fixed = TRUE),
          length(reverse_summary) == 1L,
          grepl('data-accession="B"', reverse_summary, fixed = TRUE),
          grepl('data-landmark-id="B_ONLY"', reverse_summary, fixed = TRUE),
          grepl('data-start-col="7"', reverse_summary, fixed = TRUE),
          grepl('data-end-col="16"', reverse_summary, fixed = TRUE),
          grepl('data-clipped="true"', reverse_summary, fixed = TRUE))
absent_summary <- overview(TRUE, summary_fixture[summary_fixture$accession == "B", , drop = FALSE], c(5L, 8L))
hidden_summary <- overview(FALSE, summary_fixture, c(5L, 8L))
stopifnot(length(summary_tags(absent_summary)) == 0L,
          grepl("No frozen gene or region boundary overlaps this window.", absent_summary, fixed = TRUE),
          length(summary_tags(hidden_summary)) == 0L,
          !grepl('class="fungal-alignment-window-regions"', hidden_summary, fixed = TRUE),
          !grepl("Regions in this window", hidden_summary, fixed = TRUE))
message("Fungal region checks passed: own-template gaps/strand, scoped curated anchors, frozen versions, clipped tracks/cells, visible-window summaries and annotation toggle.")
