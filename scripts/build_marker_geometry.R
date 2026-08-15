#!/usr/bin/env Rscript

project_root <- normalizePath(getwd(), mustWork = TRUE)
source(file.path(project_root, "R", "functions.R"))

catalog_dir <- file.path(project_root, "data", "catalog")
derived_dir <- file.path(project_root, "data", "derived")
pairs <- read.csv(file.path(catalog_dir, "primer_pairs.csv"), stringsAsFactors = FALSE)
oligos <- read.csv(file.path(catalog_dir, "oligos.csv"), stringsAsFactors = FALSE)
links <- read.csv(file.path(catalog_dir, "pair_oligos.csv"), stringsAsFactors = FALSE)
facets <- read.csv(file.path(catalog_dir, "primer_pair_facets.csv"), stringsAsFactors = FALSE)

coi <- read.csv(
  file.path(derived_dir, "pair_geometry.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
coi$marker_id <- "COI"
coi$placement_status <- "sequence_aligned"

its_reference_path <- file.path(
  project_root, "data", "reference", "its_fungal_reference_FN812768.2.fasta"
)
if (!file.exists(its_reference_path)) {
  stop("Run scripts/download_marker_references.R before building marker geometry.")
}
its_reference <- clean_sequence(parse_fasta(its_reference_path)[[1]])

pair_primers <- function(pair_id) {
  pair_links <- links[links$pair_id == pair_id, , drop = FALSE]
  result <- merge(pair_links, oligos, by = "oligo_id", sort = FALSE)
  result[match(c("forward", "reverse"), result$direction), , drop = FALSE]
}

facet_text <- function(pair_id, type) {
  values <- facets$facet_value[
    facets$pair_id == pair_id & facets$facet_type == type
  ]
  paste(unique(values), collapse = " / ")
}

locate_reference_pair <- function(pair_row) {
  primer_rows <- pair_primers(pair_row$pair_id)
  forward_sequence <- primer_rows$sequence[primer_rows$direction == "forward"]
  reverse_sequence <- primer_rows$sequence[primer_rows$direction == "reverse"]
  forward <- primer_candidate_table(its_reference, forward_sequence, "forward")
  reverse <- primer_candidate_table(its_reference, reverse_sequence, "reverse")
  candidates <- merge(forward, reverse, by = NULL, suffixes = c("_forward", "_reverse"))
  candidates <- candidates[
    candidates$start_forward < candidates$start_reverse,
    , drop = FALSE
  ]
  if (!nrow(candidates)) return(NULL)
  candidates$rank_score <-
    candidates$mismatch_count_forward + candidates$mismatch_count_reverse +
    2 * (
      candidates$terminal_3_mismatches_forward +
        candidates$terminal_3_mismatches_reverse
    )
  candidates <- candidates[order(
    candidates$rank_score,
    candidates$localization_penalty_forward + candidates$localization_penalty_reverse,
    candidates$start_forward,
    candidates$start_reverse
  ), , drop = FALSE]
  best <- candidates[1, , drop = FALSE]
  credible <- isTRUE(assess_primer_binding(
    forward_sequence,
    data.frame(
      compatible_identity = best$compatible_identity_forward,
      terminal_3_mismatches = best$terminal_3_mismatches_forward,
      longest_compatible_run = best$longest_compatible_run_forward,
      random_match_evalue = best$random_match_evalue_forward
    )
  )$credible) && isTRUE(assess_primer_binding(
    reverse_sequence,
    data.frame(
      compatible_identity = best$compatible_identity_reverse,
      terminal_3_mismatches = best$terminal_3_mismatches_reverse,
      longest_compatible_run = best$longest_compatible_run_reverse,
      random_match_evalue = best$random_match_evalue_reverse
    )
  )$credible)
  targeted_start <- best$end_forward + 1L
  targeted_end <- best$start_reverse - 1L
  data.frame(
    pair_id = pair_row$pair_id,
    pair_label = pair_row$pair_label,
    use_case = facet_text(pair_row$pair_id, "application"),
    target = facet_text(pair_row$pair_id, "target_taxon"),
    forward_primer = primer_rows$primer_name[primer_rows$direction == "forward"],
    reverse_primer = primer_rows$primer_name[primer_rows$direction == "reverse"],
    forward_start = best$start_forward,
    forward_end = best$end_forward,
    reverse_start = best$start_reverse,
    reverse_end = best$end_reverse,
    pair_start = best$start_forward,
    pair_end = best$end_reverse,
    aligned_amplicon_bp = best$end_reverse - best$start_forward + 1L,
    targeted_region_start = targeted_start,
    targeted_region_end = targeted_end,
    targeted_region_bp = max(0L, targeted_end - targeted_start + 1L),
    folmer_overlap_bp = 0L,
    folmer_coverage_fraction = NA_real_,
    reference_accession = "FN812768.2",
    sources = "unite_primers",
    folmer_overlap_rank = NA_integer_,
    marker_id = "ITS_FUNGAL",
    placement_status = if (credible) "sequence_aligned" else "reference_specific_low_confidence",
    forward_identity = best$compatible_identity_forward,
    reverse_identity = best$compatible_identity_reverse,
    stringsAsFactors = FALSE
  )
}

its_pairs <- pairs[pairs$marker_id == "ITS_FUNGAL", , drop = FALSE]
its_rows <- lapply(seq_len(nrow(its_pairs)), function(i) {
  locate_reference_pair(its_pairs[i, , drop = FALSE])
})
its_rows <- its_rows[!vapply(its_rows, is.null, logical(1))]
its <- if (length(its_rows)) do.call(rbind, its_rows) else data.frame()

all_columns <- union(names(coi), names(its))
fill_columns <- function(x) {
  for (column in setdiff(all_columns, names(x))) x[[column]] <- NA
  x[, all_columns, drop = FALSE]
}
geometry <- rbind(fill_columns(coi), fill_columns(its))
write.csv(
  geometry, file.path(derived_dir, "marker_pair_geometry.csv"),
  row.names = FALSE, na = ""
)
message(
  "Wrote geometry for ", nrow(geometry), " pairs across ",
  length(unique(geometry$marker_id)), " markers."
)
