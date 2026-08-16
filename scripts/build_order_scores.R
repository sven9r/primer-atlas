#!/usr/bin/env Rscript

project_root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(project_root, ".Rlib"), .libPaths()))
suppressPackageStartupMessages(library(PrimerMiner))
source(file.path(project_root, "R", "functions.R"))

primers <- read.csv(
  file.path(project_root, "data", "primers.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
bindings <- read.csv(
  file.path(project_root, "data", "derived", "primer_bindings.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
full_alignment_dir <- file.path(project_root, "data", "primerminer", "full_aligned_reference")
alignment_dir <- if (dir.exists(full_alignment_dir) && length(list.files(full_alignment_dir, pattern = "_COX1_full_reference_aligned\\.fasta$"))) {
  full_alignment_dir
} else {
  file.path(project_root, "data", "primerminer", "aligned_reference")
}
alignment_files <- list.files(
  alignment_dir,
  pattern = "_COX1(_full)?_reference_aligned\\.fasta$",
  full.names = TRUE
)
if (!length(alignment_files)) {
  stop("No installed reference-coordinate alignments found.")
}

output_dir <- file.path(project_root, "data", "derived")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

order_from_path <- function(path) {
  sub("_COX1(_full)?_reference_aligned\\.fasta$", "", basename(path))
}

quantile_safe <- function(x, probability) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  unname(quantile(x, probability, na.rm = TRUE))
}

median_safe <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x)
}

primer_template_rows <- list()
template_index <- 1L
position_row_count <- 0L
position_output_path <- file.path(output_dir, "primer_position_scores.csv")
position_temp_path <- paste0(position_output_path, ".tmp")
if (file.exists(position_temp_path)) unlink(position_temp_path)
on.exit({
  if (file.exists(position_temp_path)) unlink(position_temp_path)
}, add = TRUE)

for (alignment_path in alignment_files) {
  order_name <- order_from_path(alignment_path)
  message("Evaluating ", order_name)

  for (i in seq_len(nrow(primers))) {
    primer <- primers[i, ]
    binding <- bindings[
      bindings$pair_id == primer$pair_id &
        bindings$primer_name == primer$primer_name &
        bindings$direction == primer$direction,
      ,
      drop = FALSE
    ]
    if (nrow(binding) != 1L) {
      stop("Expected one binding coordinate for ", primer$pair_id, "/", primer$primer_name)
    }

    evaluated <- suppressMessages(PrimerMiner::evaluate_primer(
      alignment_imp = alignment_path,
      primer_sequ = clean_sequence(primer$sequence),
      start = binding$alignment_start,
      stop = binding$alignment_end,
      forward = primer$direction == "forward",
      gap_NA = TRUE,
      N_NA = TRUE,
      mm_position = "Position_v1",
      mm_type = "Type_v1",
      adjacent = 2,
      sequ_names = TRUE
    ))
    evaluated <- evaluated[
      !grepl("^NC_001322\\.1_COX1_reference", evaluated$Template),
      ,
      drop = FALSE
    ]

    score_columns <- grep("^V[0-9]+$", names(evaluated), value = TRUE)
    score_matrix <- as.matrix(evaluated[, score_columns, drop = FALSE])
    mismatch_matrix <- score_matrix > 0
    mismatch_matrix[is.na(mismatch_matrix)] <- FALSE
    missing_matrix <- is.na(score_matrix)

    adjacent_count <- apply(mismatch_matrix, 1, function(x) {
      if (length(x) < 2L) return(0L)
      sum(x[-length(x)] & x[-1L])
    })
    distance_values <- as.integer(sub("^V", "", score_columns))
    terminal_count <- function(max_distance) {
      columns <- which(distance_values <= max_distance)
      if (!length(columns)) return(rep(0L, nrow(evaluated)))
      rowSums(mismatch_matrix[, columns, drop = FALSE])
    }

    primer_template_rows[[template_index]] <- data.frame(
      order = order_name,
      pair_id = primer$pair_id,
      pair_label = primer$pair_label,
      primer_name = primer$primer_name,
      direction = primer$direction,
      template = evaluated$Template,
      binding_sequence_primer_oriented = evaluated$sequ,
      alignment_start = binding$alignment_start,
      alignment_end = binding$alignment_end,
      primerminer_penalty = evaluated$sum,
      scorable = !is.na(evaluated$sum),
      missing_binding_positions = rowSums(missing_matrix),
      mismatch_count = rowSums(mismatch_matrix),
      terminal_1_mismatches = terminal_count(1L),
      terminal_3_mismatches = terminal_count(3L),
      terminal_5_mismatches = terminal_count(5L),
      adjacent_mismatch_pairs = adjacent_count,
      stringsAsFactors = FALSE
    )
    template_index <- template_index + 1L

    primer_bases <- split_bases(primer$sequence)
    template_bases <- lapply(evaluated$sequ, split_bases)
    position_block <- vector("list", length(score_columns))
    for (column_index in seq_along(score_columns)) {
      distance <- distance_values[column_index]
      primer_index <- length(primer_bases) - distance + 1L
      base_values <- vapply(
        template_bases,
        function(x) if (length(x) >= primer_index) x[primer_index] else NA_character_,
        character(1)
      )
      position_block[[column_index]] <- data.frame(
        order = order_name,
        pair_id = primer$pair_id,
        pair_label = primer$pair_label,
        primer_name = primer$primer_name,
        direction = primer$direction,
        template = evaluated$Template,
        alignment_position = if (primer$direction == "forward") {
          binding$alignment_start + primer_index - 1L
        } else {
          binding$alignment_end - primer_index + 1L
        },
        primer_position_5_to_3 = primer_index,
        distance_from_3prime = distance,
        primer_base = primer_bases[primer_index],
        template_base = base_values,
        position_penalty = score_matrix[, column_index],
        mismatch = mismatch_matrix[, column_index],
        unavailable_reason = ifelse(
          !missing_matrix[, column_index],
          "",
          ifelse(base_values == "-", "alignment gap", ifelse(base_values == "N", "ambiguous N", "unscorable"))
        ),
        stringsAsFactors = FALSE
      )
    }
    position_block <- do.call(rbind, position_block)
    write.table(
      position_block,
      file = position_temp_path,
      sep = ",",
      quote = TRUE,
      row.names = FALSE,
      col.names = !file.exists(position_temp_path),
      append = file.exists(position_temp_path),
      na = ""
    )
    position_row_count <- position_row_count + nrow(position_block)
  }
}

primer_template_scores <- do.call(rbind, primer_template_rows)

forward <- primer_template_scores[
  primer_template_scores$direction == "forward",
  ,
  drop = FALSE
]
reverse <- primer_template_scores[
  primer_template_scores$direction == "reverse",
  ,
  drop = FALSE
]
names(forward)[!names(forward) %in% c("order", "pair_id", "pair_label", "template")] <- paste0(
  "forward_", names(forward)[!names(forward) %in% c("order", "pair_id", "pair_label", "template")]
)
names(reverse)[!names(reverse) %in% c("order", "pair_id", "pair_label", "template")] <- paste0(
  "reverse_", names(reverse)[!names(reverse) %in% c("order", "pair_id", "pair_label", "template")]
)

pair_template_scores <- merge(
  forward,
  reverse,
  by = c("order", "pair_id", "pair_label", "template"),
  all = TRUE,
  sort = FALSE
)
pair_template_scores$pair_scorable <-
  pair_template_scores$forward_scorable %in% TRUE &
  pair_template_scores$reverse_scorable %in% TRUE
pair_template_scores$pair_penalty <- ifelse(
  pair_template_scores$pair_scorable,
  pair_template_scores$forward_primerminer_penalty +
    pair_template_scores$reverse_primerminer_penalty,
  NA_real_
)
pair_template_scores$perfect_match_pair <-
  pair_template_scores$pair_scorable &
  pair_template_scores$forward_mismatch_count == 0L &
  pair_template_scores$reverse_mismatch_count == 0L
pair_template_scores$accession <- sub(
  "[|].*$", "",
  sub(" .*", "", sub("^_R_", "", pair_template_scores$template))
)

groups <- split(
  pair_template_scores,
  interaction(pair_template_scores$order, pair_template_scores$pair_id, drop = TRUE)
)
order_pair_summary <- do.call(rbind, lapply(groups, function(x) {
  scorable <- x$pair_scorable %in% TRUE
  data.frame(
    order = x$order[1],
    pair_id = x$pair_id[1],
    pair_label = x$pair_label[1],
    n_templates = nrow(x),
    n_pair_scorable = sum(scorable),
    pair_scorable_fraction = mean(scorable),
    unavailable_fraction = mean(!scorable),
    perfect_match_fraction_among_scorable = if (sum(scorable)) {
      mean(x$perfect_match_pair[scorable])
    } else {
      NA_real_
    },
    forward_median_penalty = median_safe(x$forward_primerminer_penalty),
    forward_p90_penalty = quantile_safe(x$forward_primerminer_penalty, 0.9),
    reverse_median_penalty = median_safe(x$reverse_primerminer_penalty),
    reverse_p90_penalty = quantile_safe(x$reverse_primerminer_penalty, 0.9),
    pair_median_penalty = median_safe(x$pair_penalty),
    pair_p90_penalty = quantile_safe(x$pair_penalty, 0.9),
    terminal_1_mismatch_fraction = if (sum(scorable)) {
      mean(
        x$forward_terminal_1_mismatches[scorable] +
          x$reverse_terminal_1_mismatches[scorable] > 0
      )
    } else {
      NA_real_
    },
    terminal_3_mismatch_fraction = if (sum(scorable)) {
      mean(
        x$forward_terminal_3_mismatches[scorable] +
          x$reverse_terminal_3_mismatches[scorable] > 0
      )
    } else {
      NA_real_
    },
    adjacent_mismatch_fraction = if (sum(scorable)) {
      mean(
        x$forward_adjacent_mismatch_pairs[scorable] +
          x$reverse_adjacent_mismatch_pairs[scorable] > 0
      )
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}))
order_pair_summary <- order_pair_summary[
  order(order_pair_summary$order, order_pair_summary$pair_label),
  ,
  drop = FALSE
]

model_metadata <- data.frame(
  field = c(
    "model_name", "implementation", "position_matrix", "type_matrix",
    "adjacent_multiplier", "gap_handling", "N_handling", "interpretation"
  ),
  value = c(
    "PrimerMiner mismatch penalty",
    "PrimerMiner 0.22 evaluate_primer",
    "Position_v1",
    "Type_v1",
    "2",
    "unavailable, reported separately from penalty",
    "unavailable, reported separately from penalty",
    "relative mismatch penalty; not PCR probability and not a validated universal failure threshold"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  primer_template_scores,
  file.path(output_dir, "primer_template_scores.csv"),
  row.names = FALSE,
  na = ""
)
if (!file.rename(position_temp_path, position_output_path)) {
  stop("Could not atomically replace ", position_output_path)
}
write.csv(
  pair_template_scores,
  file.path(output_dir, "pair_template_scores.csv"),
  row.names = FALSE,
  na = ""
)
for (pair_id in unique(pair_template_scores$pair_id)) {
  pair_rows <- pair_template_scores[pair_template_scores$pair_id == pair_id, , drop = FALSE]
  path <- file.path(
    output_dir,
    paste0("full_order_", tolower(pair_id), "_sequence_scores.csv.gz")
  )
  connection <- gzfile(path, open = "wt")
  write.csv(pair_rows, connection, row.names = FALSE, na = "")
  close(connection)
}
write.csv(
  order_pair_summary,
  file.path(output_dir, "order_pair_summary.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  model_metadata,
  file.path(output_dir, "penalty_model_metadata.csv"),
  row.names = FALSE
)

cat(
  "Wrote", nrow(primer_template_scores), "primer-template rows,",
  position_row_count, "position-level rows,",
  nrow(pair_template_scores), "pair-template rows, and",
  nrow(order_pair_summary), "order-pair summaries.\n"
)
