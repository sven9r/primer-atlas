#!/usr/bin/env Rscript

project_root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(project_root, ".Rlib"), .libPaths()))
suppressPackageStartupMessages(library(PrimerMiner))
source(file.path(project_root, "R", "functions.R"))

primer_path <- file.path(project_root, "data", "primers.csv")
reference_path <- file.path(
  project_root, "data", "reference", "coi_reference_NC_001322.1.fasta"
)
output_dir <- file.path(project_root, "data", "derived")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

primers <- read.csv(primer_path, stringsAsFactors = FALSE, check.names = FALSE)
reference <- parse_fasta(reference_path)
if (length(reference) != 1L) stop("Expected one sequence in the COI reference FASTA.")
reference_sequence <- clean_sequence(reference[[1]])
reference_length <- nchar(reference_sequence)
primerminer_reference <- tempfile(fileext = ".fasta")
writeLines(
  c(">reference_copy_1", reference_sequence, ">reference_copy_2", reference_sequence),
  primerminer_reference
)
on.exit(unlink(primerminer_reference), add = TRUE)

candidate_tables <- lapply(seq_len(nrow(primers)), function(i) {
  primer_candidate_table(
    reference_sequence, primers$sequence[i], primers$direction[i]
  )
})
selected_candidate_rows <- rep(1L, nrow(primers))
placement_mode <- rep("individual sequence best", nrow(primers))
pair_expected_amplicon_bp <- rep(NA_integer_, nrow(primers))
pair_full_amplicon_bp <- rep(NA_integer_, nrow(primers))
pair_informative_bp <- rep(NA_integer_, nrow(primers))
pair_size_delta_bp <- rep(NA_integer_, nrow(primers))

# A lineage-focused primer can have a lower-mismatch accidental site on an
# unrelated reference. Keep sequence matching as the default, but resolve an
# impossible or size-inconsistent pair jointly using the published product
# size. The size can refer either to the primer-inclusive product or to the
# informative region, both conventions used in the source literature.
for (pair_id in unique(primers$pair_id)) {
  pair_indices <- which(primers$pair_id == pair_id)
  forward_index <- pair_indices[primers$direction[pair_indices] == "forward"]
  reverse_index <- pair_indices[primers$direction[pair_indices] == "reverse"]
  expected_sizes <- unique(stats::na.omit(
    suppressWarnings(as.integer(primers$reported_amplicon_bp[pair_indices]))
  ))
  if (
    length(forward_index) != 1L ||
      length(reverse_index) != 1L ||
      length(expected_sizes) != 1L
  ) {
    next
  }

  expected_size <- expected_sizes[[1]]
  forward_candidates <- head(candidate_tables[[forward_index]], 100L)
  reverse_candidates <- head(candidate_tables[[reverse_index]], 100L)
  combinations <- expand.grid(
    forward_row = seq_len(nrow(forward_candidates)),
    reverse_row = seq_len(nrow(reverse_candidates))
  )
  combinations$forward_start <-
    forward_candidates$start[combinations$forward_row]
  combinations$forward_end <-
    forward_candidates$end[combinations$forward_row]
  combinations$reverse_start <-
    reverse_candidates$start[combinations$reverse_row]
  combinations$reverse_end <-
    reverse_candidates$end[combinations$reverse_row]
  combinations <- combinations[
    combinations$reverse_start > combinations$forward_end,
    ,
    drop = FALSE
  ]
  if (!nrow(combinations)) next

  combinations$full_amplicon_bp <-
    combinations$reverse_end - combinations$forward_start + 1L
  combinations$informative_bp <-
    combinations$reverse_start - combinations$forward_end - 1L
  combinations$size_delta_bp <- pmin(
    abs(combinations$full_amplicon_bp - expected_size),
    abs(combinations$informative_bp - expected_size)
  )
  combinations$terminal_mismatches <-
    forward_candidates$terminal_5_mismatches[combinations$forward_row] +
    reverse_candidates$terminal_5_mismatches[combinations$reverse_row]
  combinations$total_mismatches <-
    forward_candidates$mismatch_count[combinations$forward_row] +
    reverse_candidates$mismatch_count[combinations$reverse_row]
  combinations$total_penalty <-
    forward_candidates$localization_penalty[combinations$forward_row] +
    reverse_candidates$localization_penalty[combinations$reverse_row]
  combinations <- combinations[
    order(
      combinations$size_delta_bp,
      combinations$terminal_mismatches,
      combinations$total_mismatches,
      combinations$total_penalty
    ),
    ,
    drop = FALSE
  ]

  independent_forward <- candidate_tables[[forward_index]][1, ]
  independent_reverse <- candidate_tables[[reverse_index]][1, ]
  independent_order_valid <-
    independent_reverse$start > independent_forward$end
  independent_full <- if (independent_order_valid) {
    independent_reverse$end - independent_forward$start + 1L
  } else {
    NA_integer_
  }
  independent_informative <- if (independent_order_valid) {
    independent_reverse$start - independent_forward$end - 1L
  } else {
    NA_integer_
  }
  independent_delta <- if (independent_order_valid) {
    min(
      abs(independent_full - expected_size),
      abs(independent_informative - expected_size)
    )
  } else {
    Inf
  }
  allowed_delta <- max(5L, as.integer(round(expected_size * 0.03)))

  if (!independent_order_valid || independent_delta > allowed_delta) {
    selected <- combinations[1, ]
    selected_candidate_rows[forward_index] <- selected$forward_row
    selected_candidate_rows[reverse_index] <- selected$reverse_row
    placement_mode[pair_indices] <-
      "pair-aware published-size constraint"
  }

  selected_forward <-
    candidate_tables[[forward_index]][selected_candidate_rows[forward_index], ]
  selected_reverse <-
    candidate_tables[[reverse_index]][selected_candidate_rows[reverse_index], ]
  selected_full <- selected_reverse$end - selected_forward$start + 1L
  selected_informative <-
    selected_reverse$start - selected_forward$end - 1L
  selected_delta <- min(
    abs(selected_full - expected_size),
    abs(selected_informative - expected_size)
  )
  pair_expected_amplicon_bp[pair_indices] <- expected_size
  pair_full_amplicon_bp[pair_indices] <- selected_full
  pair_informative_bp[pair_indices] <- selected_informative
  pair_size_delta_bp[pair_indices] <- selected_delta
}

binding_rows <- vector("list", nrow(primers))
candidate_rows <- vector("list", nrow(primers))

for (i in seq_len(nrow(primers))) {
  primer <- primers[i, ]
  candidates <- candidate_tables[[i]]
  best <- candidates[selected_candidate_rows[i], ]
  same_best <- candidates[
    candidates$mismatch_count == best$mismatch_count &
      candidates$terminal_5_mismatches == best$terminal_5_mismatches &
      abs(candidates$localization_penalty - best$localization_penalty) < 1e-9,
    ,
    drop = FALSE
  ]
  runner_up <- candidates[min(2L, nrow(candidates)), ]

  pm <- PrimerMiner::evaluate_primer(
    # PrimerMiner 0.22 drops matrix dimensions for a one-sequence alignment.
    # Two identical temporary copies preserve the intended calculation.
    alignment_imp = primerminer_reference,
    primer_sequ = clean_sequence(primer$sequence),
    start = best$start,
    stop = best$end,
    forward = primer$direction == "forward",
    gap_NA = TRUE,
    N_NA = TRUE,
    mm_position = "Position_v1",
    mm_type = "Type_v1",
    adjacent = 2,
    sequ_names = TRUE
  )

  binding_rows[[i]] <- data.frame(
    pair_id = primer$pair_id,
    pair_label = primer$pair_label,
    primer_name = primer$primer_name,
    direction = primer$direction,
    sequence = primer$sequence,
    alignment_start = best$start,
    alignment_end = best$end,
    primer_length = nchar(clean_sequence(primer$sequence)),
    reference_mismatches = best$mismatch_count,
    terminal_5_mismatches = best$terminal_5_mismatches,
    localization_penalty = round(best$localization_penalty, 3),
    primerminer_reference_penalty = pm$sum[1],
    equally_best_sites = nrow(same_best),
    individually_best_start = candidates$start[1],
    individually_best_end = candidates$end[1],
    placement_mode = placement_mode[i],
    pair_expected_amplicon_bp = pair_expected_amplicon_bp[i],
    pair_full_amplicon_bp = pair_full_amplicon_bp[i],
    pair_informative_bp = pair_informative_bp[i],
    pair_size_delta_bp = pair_size_delta_bp[i],
    second_candidate_start = runner_up$start,
    second_candidate_mismatches = runner_up$mismatch_count,
    second_candidate_penalty = round(runner_up$localization_penalty, 3),
    reference_window = best$reference_window,
    coordinate_method = if (placement_mode[i] ==
      "pair-aware published-size constraint") {
      paste(
        "exhaustive IUPAC-compatible search with forward/reverse order and",
        "published product-size constraint; size may be primer-inclusive or",
        "primer-excluded, matching source conventions"
      )
    } else {
      paste(
        "exhaustive IUPAC-compatible search; ranked by mismatch count,",
        "terminal-5 mismatch count, then PrimerMiner-style localization penalty;",
        "pair order and published product size audited"
      )
    },
    reference_accession = "NC_001322.1",
    reference_feature = "COX1 1474..3009",
    reference_length_bp = reference_length,
    stringsAsFactors = FALSE
  )

  displayed_candidate_rows <- unique(c(
    seq_len(min(5L, nrow(candidates))),
    selected_candidate_rows[i]
  ))
  candidate_rows[[i]] <- cbind(
    data.frame(
      pair_id = primer$pair_id,
      primer_name = primer$primer_name,
      direction = primer$direction,
      candidate_rank = displayed_candidate_rows,
      selected_for_map = displayed_candidate_rows == selected_candidate_rows[i],
      stringsAsFactors = FALSE
    ),
    candidates[displayed_candidate_rows, , drop = FALSE]
  )
}

bindings <- do.call(rbind, binding_rows)
top_candidates <- do.call(rbind, candidate_rows)

folmer_forward <- bindings[
  bindings$pair_id == "FOLMER" & bindings$direction == "forward", ,
  drop = FALSE
]
folmer_reverse <- bindings[
  bindings$pair_id == "FOLMER" & bindings$direction == "reverse", ,
  drop = FALSE
]
if (nrow(folmer_forward) != 1L || nrow(folmer_reverse) != 1L) {
  stop("The FOLMER pair must contain exactly one forward and one reverse primer.")
}
# Match the convention in the Elbrecht overview: the targeted region is the
# sequence strictly between the forward and reverse primer 3-prime ends.
# The full amplicon (including primer binding sequences) remains separate.
folmer_start <- folmer_forward$alignment_end + 1L
folmer_end <- folmer_reverse$alignment_start - 1L

pair_ids <- unique(primers$pair_id)
pair_rows <- lapply(pair_ids, function(pair_id) {
  pair_primers <- primers[primers$pair_id == pair_id, , drop = FALSE]
  pair_bindings <- bindings[bindings$pair_id == pair_id, , drop = FALSE]
  forward <- pair_bindings[pair_bindings$direction == "forward", , drop = FALSE]
  reverse <- pair_bindings[pair_bindings$direction == "reverse", , drop = FALSE]
  if (nrow(forward) != 1L || nrow(reverse) != 1L) return(NULL)

  pair_start <- forward$alignment_start
  pair_end <- reverse$alignment_end
  target_start <- forward$alignment_end + 1L
  target_end <- reverse$alignment_start - 1L
  targeted_region_bp <- max(0L, target_end - target_start + 1L)
  overlap <- max(0L, min(target_end, folmer_end) - max(target_start, folmer_start) + 1L)
  data.frame(
    pair_id = pair_id,
    pair_label = pair_primers$pair_label[1],
    use_case = pair_primers$use_case[1],
    target = paste(unique(pair_primers$target), collapse = " / "),
    forward_primer = forward$primer_name,
    reverse_primer = reverse$primer_name,
    forward_start = forward$alignment_start,
    forward_end = forward$alignment_end,
    reverse_start = reverse$alignment_start,
    reverse_end = reverse$alignment_end,
    pair_start = pair_start,
    pair_end = pair_end,
    aligned_amplicon_bp = pair_end - pair_start + 1L,
    targeted_region_start = target_start,
    targeted_region_end = target_end,
    targeted_region_bp = targeted_region_bp,
    folmer_overlap_bp = overlap,
    folmer_coverage_fraction = overlap / (folmer_end - folmer_start + 1L),
    reference_accession = "NC_001322.1",
    sources = paste(unique(pair_primers$source_key), collapse = ", "),
    stringsAsFactors = FALSE
  )
})
pair_geometry <- do.call(rbind, pair_rows)
pair_geometry <- pair_geometry[
  order(-pair_geometry$folmer_overlap_bp, pair_geometry$pair_start, pair_geometry$pair_label),
  ,
  drop = FALSE
]
pair_geometry$folmer_overlap_rank <- seq_len(nrow(pair_geometry))

folmer_region <- data.frame(
  region = "Folmer barcode region",
  start = folmer_start,
  end = folmer_end,
  length_bp = folmer_end - folmer_start + 1L,
  forward_primer = folmer_forward$primer_name,
  reverse_primer = folmer_reverse$primer_name,
  coordinate_method = paste(
    "interior strictly between alignment-derived LCO1490 and HCO2198",
    "3-prime ends; full primer-inclusive amplicon is reported separately"
  ),
  reference_accession = "NC_001322.1",
  stringsAsFactors = FALSE
)

write.csv(bindings, file.path(output_dir, "primer_bindings.csv"), row.names = FALSE)
write.csv(top_candidates, file.path(output_dir, "primer_binding_top_candidates.csv"), row.names = FALSE)
write.csv(pair_geometry, file.path(output_dir, "pair_geometry.csv"), row.names = FALSE)
write.csv(folmer_region, file.path(output_dir, "folmer_region.csv"), row.names = FALSE)

cat("\nAlignment-derived primer positions:\n")
print(bindings[, c(
  "pair_id", "primer_name", "direction", "alignment_start", "alignment_end",
  "reference_mismatches", "primerminer_reference_penalty", "equally_best_sites"
)], row.names = FALSE)
cat("\nAlignment-derived pair geometry, sorted by Folmer overlap:\n")
print(pair_geometry[, c(
  "pair_label", "pair_start", "pair_end", "targeted_region_bp", "aligned_amplicon_bp",
  "folmer_overlap_bp"
)], row.names = FALSE)
cat("\nFolmer region:\n")
print(folmer_region, row.names = FALSE)
