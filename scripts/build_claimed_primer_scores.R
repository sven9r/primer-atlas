#!/usr/bin/env Rscript

# Score diet/gut-content and other target-claim primers on one shared,
# taxonomically annotated arthropod COI reference. The goal is to turn an
# order-level warning into traceable family/genus evidence, not to manufacture
# a universal PCR pass/fail threshold.

project_root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(project_root, ".Rlib"), .libPaths()))

suppressPackageStartupMessages(library(PrimerMiner))
source(file.path(project_root, "R", "functions.R"))

external_dir <- file.path(project_root, "data", "external", "gurten2026")
derived_dir <- file.path(project_root, "data", "derived")
cache_dir <- file.path(derived_dir, "expanded_primer_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

alignment_path <- file.path(external_dir, "ClusteredReferences.fasta")
taxonomy_path <- file.path(derived_dir, "gurten2026_centroid_taxonomy.csv")
bindings_path <- file.path(derived_dir, "primer_bindings.csv")
primer_path <- file.path(project_root, "data", "primers.csv")
claims_path <- file.path(project_root, "data", "primer_target_claims.csv")
reference_manifest_path <- file.path(
  derived_dir,
  "beeprime_reference_manifest.csv"
)

required <- c(
  alignment_path, taxonomy_path, bindings_path, primer_path, claims_path,
  reference_manifest_path
)
missing <- required[!file.exists(required)]
if (length(missing)) {
  stop("Missing expanded-reference input(s): ", paste(missing, collapse = ", "))
}

primers <- read.csv(primer_path, stringsAsFactors = FALSE, check.names = FALSE)
bindings <- read.csv(bindings_path, stringsAsFactors = FALSE, check.names = FALSE)
claims <- read.csv(claims_path, stringsAsFactors = FALSE, check.names = FALSE)
taxonomy <- read.csv(taxonomy_path, stringsAsFactors = FALSE, check.names = FALSE)
reference_manifest <- read.csv(
  reference_manifest_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
aligned_sequences <- parse_fasta(alignment_path)

if (length(aligned_sequences) != nrow(taxonomy)) {
  stop(
    "Alignment/taxonomy row mismatch: ", length(aligned_sequences),
    " sequences versus ", nrow(taxonomy), " taxonomy rows."
  )
}
alignment_accessions <- sub(" .*", "", names(aligned_sequences))
if (!setequal(alignment_accessions, taxonomy$accession)) {
  stop("Alignment accession labels do not match the taxonomy table.")
}

default_pairs <- unique(primers$pair_id)
pair_ids <- strsplit(
  Sys.getenv("CLAIM_PAIRS", paste(default_pairs, collapse = ",")),
  ",",
  fixed = TRUE
)[[1]]
pair_ids <- trimws(pair_ids[nzchar(trimws(pair_ids))])
unknown_pairs <- setdiff(pair_ids, unique(primers$pair_id))
if (length(unknown_pairs)) {
  stop("Unknown pair ID(s): ", paste(unknown_pairs, collapse = ", "))
}

force_scores <- identical(tolower(Sys.getenv("FORCE", "false")), "true")

# Supplement 5 carries a stable 40-column alignment offset relative to the
# NC_001322.1 COX1 coordinate used in primer_bindings.csv. Derive that offset
# from both BeePrime sites and fail if the anchors disagree.
manifest_value <- function(field) {
  reference_manifest$value[match(field, reference_manifest$field)]
}
column_start <- function(value) as.integer(sub("-.*$", "", value))
column_end <- function(value) as.integer(sub("^.*-", "", value))
bee_forward <- bindings[
  bindings$pair_id == "BEEPRIME" & bindings$direction == "forward",
  ,
  drop = FALSE
]
bee_reverse <- bindings[
  bindings$pair_id == "BEEPRIME" & bindings$direction == "reverse",
  ,
  drop = FALSE
]
forward_manifest <- manifest_value("beeprime_forward_alignment_columns")
reverse_manifest <- manifest_value("beeprime_reverse_alignment_columns")
offsets <- c(
  column_start(forward_manifest) - bee_forward$alignment_start,
  column_end(forward_manifest) - bee_forward$alignment_end,
  column_start(reverse_manifest) - bee_reverse$alignment_start,
  column_end(reverse_manifest) - bee_reverse$alignment_end
)
if (length(unique(offsets)) != 1L) {
  stop("Gurten alignment offset anchors disagree: ", paste(offsets, collapse = ", "))
}
alignment_offset <- unique(offsets)

evaluate_primer_panel <- function(primer_row, binding_row) {
  cache_name <- paste0(
    gsub("[^A-Za-z0-9_.-]+", "_", primer_row$primer_name),
    "_gurten2026_scores.rds"
  )
  cache_path <- file.path(cache_dir, cache_name)
  if (file.exists(cache_path) && !force_scores) {
    cached <- readRDS(cache_path)
    if (nrow(cached) == length(aligned_sequences)) {
      message("Using cached ", primer_row$primer_name)
      return(cached)
    }
  }

  start_column <- binding_row$alignment_start + alignment_offset
  end_column <- binding_row$alignment_end + alignment_offset
  if (
    start_column < 1L ||
      end_column > unique(nchar(aligned_sequences)) ||
      end_column - start_column + 1L != nchar(primer_row$sequence)
  ) {
    stop("Invalid expanded-reference columns for ", primer_row$primer_name)
  }

  message(
    "Evaluating ", primer_row$primer_name, " at columns ",
    start_column, "-", end_column
  )
  evaluated <- suppressMessages(PrimerMiner::evaluate_primer(
    alignment_imp = alignment_path,
    primer_sequ = clean_sequence(primer_row$sequence),
    start = start_column,
    stop = end_column,
    forward = primer_row$direction == "forward",
    gap_NA = TRUE,
    N_NA = TRUE,
    mm_position = "Position_v1",
    mm_type = "Type_v1",
    adjacent = 2,
    sequ_names = TRUE
  ))
  score_columns <- grep("^V[0-9]+$", names(evaluated), value = TRUE)
  score_matrix <- as.matrix(evaluated[, score_columns, drop = FALSE])
  mismatch_matrix <- score_matrix > 0
  mismatch_matrix[is.na(mismatch_matrix)] <- FALSE
  missing_matrix <- is.na(score_matrix)
  distances <- as.integer(sub("^V", "", score_columns))

  terminal_count <- function(max_distance) {
    columns <- which(distances <= max_distance)
    if (!length(columns)) return(rep(0L, nrow(evaluated)))
    rowSums(mismatch_matrix[, columns, drop = FALSE])
  }

  result <- data.frame(
    accession = evaluated$Template,
    binding_sequence_primer_oriented = evaluated$sequ,
    penalty = evaluated$sum,
    scorable = !is.na(evaluated$sum),
    missing_binding_positions = rowSums(missing_matrix),
    mismatch_count = rowSums(mismatch_matrix),
    terminal_1_mismatches = terminal_count(1L),
    terminal_3_mismatches = terminal_count(3L),
    terminal_5_mismatches = terminal_count(5L),
    stringsAsFactors = FALSE
  )
  saveRDS(result, cache_path, compress = "xz")
  result
}

median_safe <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) median(x) else NA_real_
}
mean_safe <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}
quantile_safe <- function(x, probability) {
  x <- x[is.finite(x)]
  if (length(x)) unname(quantile(x, probability, na.rm = TRUE)) else NA_real_
}

summarize_scores <- function(data, group_columns) {
  usable <- data
  for (column in group_columns) {
    usable[[column]][
      is.na(usable[[column]]) | !nzchar(trimws(usable[[column]]))
    ] <- paste0(column, " unresolved")
  }
  group_key <- interaction(
    usable[, group_columns, drop = FALSE],
    drop = TRUE,
    lex.order = TRUE
  )
  groups <- split(usable, group_key)
  result <- do.call(rbind, lapply(groups, function(x) {
    scorable <- x$pair_scorable %in% TRUE
    group_values <- x[1, group_columns, drop = FALSE]
    metrics <- data.frame(
      n_centroids = nrow(x),
      source_records_represented = sum(x$cluster_size, na.rm = TRUE),
      n_forward_scorable = sum(x$forward_scorable %in% TRUE),
      forward_scorable_fraction = mean(x$forward_scorable %in% TRUE),
      n_reverse_scorable = sum(x$reverse_scorable %in% TRUE),
      reverse_scorable_fraction = mean(x$reverse_scorable %in% TRUE),
      n_pair_scorable = sum(scorable),
      pair_scorable_fraction = mean(scorable),
      forward_median_penalty = median_safe(x$forward_penalty),
      reverse_median_penalty = median_safe(x$reverse_penalty),
      pair_median_penalty = median_safe(x$pair_penalty),
      pair_mean_penalty = mean_safe(x$pair_penalty),
      pair_p90_penalty = quantile_safe(x$pair_penalty, 0.9),
      perfect_match_fraction_among_scorable = if (sum(scorable)) {
        mean(x$perfect_match_pair[scorable])
      } else {
        NA_real_
      },
      terminal_3_mismatch_fraction_among_scorable = if (sum(scorable)) {
        mean(x$terminal_3_mismatch_pair[scorable])
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
    cbind(group_values, metrics)
  }))
  rownames(result) <- NULL
  result
}

score_cache <- list()
summary_overall <- list()
summary_order <- list()
summary_family <- list()
summary_genus <- list()
summary_subfamily <- list()
summary_phylogeny_group <- list()
position_manifest <- list()

acari_orders <- c(
  "Ixodida", "Mesostigmata", "Sarcoptiformes", "Trombidiformes",
  "Opilioacarida", "Holothyrida"
)
collembola_orders <- c(
  "Entomobryomorpha", "Poduromorpha", "Symphypleona", "Neelipleona"
)
butterfly_families <- c(
  "Hedylidae", "Hesperiidae", "Lycaenidae", "Nymphalidae",
  "Papilionidae", "Pieridae", "Riodinidae"
)
macroheterocera_families <- c(
  "Apatelodidae", "Bombycidae", "Drepanidae", "Erebidae", "Eupterotidae",
  "Euteliidae", "Geometridae", "Lasiocampidae", "Lemoniidae", "Noctuidae",
  "Nolidae", "Notodontidae", "Saturniidae", "Sphingidae", "Uraniidae"
)
lepidoptera_boundary_families <- c(
  "Castniidae", "Megalopygidae", "Mimallonidae", "Sesiidae",
  "Thyrididae", "Zygaenidae"
)

assign_phylogeny_group <- function(data) {
  data$display_group_type <- NA_character_
  data$display_group <- NA_character_
  data$group_definition <- NA_character_
  data$group_caveat <- NA_character_

  acari <- data$order %in% acari_orders
  data$display_group_type[acari] <- "Composite clade"
  data$display_group[acari] <- "Acari"
  data$group_definition[acari] <- paste(acari_orders, collapse = " | ")
  data$group_caveat[acari] <-
    "Acari is a composite arachnid clade here, not a single modern order."

  collembola <- data$order %in% collembola_orders
  data$display_group_type[collembola] <- "Composite clade"
  data$display_group[collembola] <- "Collembola"
  data$group_definition[collembola] <- paste(collembola_orders, collapse = " | ")
  data$group_caveat[collembola] <-
    "Collembola is a composite hexapod clade here, not a single modern order."

  lepidoptera <- data$order == "Lepidoptera"
  butterfly <- lepidoptera & data$family %in% butterfly_families
  macro <- lepidoptera & data$family %in% macroheterocera_families
  boundary <- lepidoptera & data$family %in% lepidoptera_boundary_families
  micro <- lepidoptera & !butterfly & !macro & !boundary

  data$display_group_type[butterfly] <- "Lepidoptera study lens"
  data$display_group[butterfly] <- "Butterflies (Papilionoidea)"
  data$group_definition[butterfly] <- paste(butterfly_families, collapse = " | ")
  data$group_caveat[butterfly] <-
    "Butterflies are shown separately instead of being forced into a moth size bin."

  data$display_group_type[macro] <- "Lepidoptera study lens"
  data$display_group[macro] <- "Macroheterocera (macro-moth core)"
  data$group_definition[macro] <- paste(
    macroheterocera_families,
    collapse = " | "
  )
  data$group_caveat[macro] <-
    "A phylogenetic macro-moth core; this is more defensible than a body-size cutoff."

  data$display_group_type[boundary] <- "Lepidoptera study lens"
  data$display_group[boundary] <- "Boundary / convention-sensitive moths"
  data$group_definition[boundary] <- paste(
    lepidoptera_boundary_families,
    collapse = " | "
  )
  data$group_caveat[boundary] <-
    "Large-bodied families whose micro/macro placement varies by convention."

  data$display_group_type[micro] <- "Lepidoptera study lens"
  data$display_group[micro] <- "Microlepidoptera (operational grade)"
  data$group_definition[micro] <-
    "Remaining represented Lepidoptera families after explicit butterfly, Macroheterocera, and boundary assignments."
  data$group_caveat[micro] <-
    "Microlepidoptera is an operational, non-monophyletic grade."
  data
}

for (pair_id in pair_ids) {
  pair_primers <- primers[primers$pair_id == pair_id, , drop = FALSE]
  pair_bindings <- bindings[bindings$pair_id == pair_id, , drop = FALSE]
  if (
    sum(pair_primers$direction == "forward") != 1L ||
      sum(pair_primers$direction == "reverse") != 1L ||
      nrow(pair_bindings) != 2L
  ) {
    stop("Pair ", pair_id, " does not contain exactly one forward and one reverse primer.")
  }

  directional_scores <- list()
  for (direction in c("forward", "reverse")) {
    primer_row <- pair_primers[pair_primers$direction == direction, , drop = FALSE]
    binding_row <- pair_bindings[pair_bindings$direction == direction, , drop = FALSE]
    cache_key <- primer_row$primer_name
    if (is.null(score_cache[[cache_key]])) {
      score_cache[[cache_key]] <- evaluate_primer_panel(primer_row, binding_row)
    }
    directional_scores[[direction]] <- score_cache[[cache_key]]
    position_manifest[[length(position_manifest) + 1L]] <- data.frame(
      pair_id = pair_id,
      pair_label = primer_row$pair_label,
      primer_name = primer_row$primer_name,
      direction = direction,
      reference_start = binding_row$alignment_start,
      reference_end = binding_row$alignment_end,
      expanded_alignment_start = binding_row$alignment_start + alignment_offset,
      expanded_alignment_end = binding_row$alignment_end + alignment_offset,
      alignment_offset = alignment_offset,
      stringsAsFactors = FALSE
    )
  }

  forward <- directional_scores$forward
  reverse <- directional_scores$reverse
  names(forward)[-1] <- paste0("forward_", names(forward)[-1])
  names(reverse)[-1] <- paste0("reverse_", names(reverse)[-1])
  pair_scores <- merge(forward, reverse, by = "accession", sort = FALSE)
  taxonomy_match <- match(pair_scores$accession, taxonomy$accession)
  pair_scores <- cbind(
    data.frame(
      pair_id = pair_id,
      pair_label = pair_primers$pair_label[1],
      stringsAsFactors = FALSE
    ),
    pair_scores,
    taxonomy[
      taxonomy_match,
      setdiff(names(taxonomy), "accession"),
      drop = FALSE
    ]
  )
  pair_scores$pair_scorable <-
    pair_scores$forward_scorable %in% TRUE &
    pair_scores$reverse_scorable %in% TRUE
  pair_scores$pair_penalty <- ifelse(
    pair_scores$pair_scorable,
    pair_scores$forward_penalty + pair_scores$reverse_penalty,
    NA_real_
  )
  pair_scores$perfect_match_pair <-
    pair_scores$pair_scorable &
    pair_scores$forward_mismatch_count == 0L &
    pair_scores$reverse_mismatch_count == 0L
  pair_scores$terminal_3_mismatch_pair <-
    pair_scores$pair_scorable &
    (
      pair_scores$forward_terminal_3_mismatches +
        pair_scores$reverse_terminal_3_mismatches
    ) > 0L

  # A mismatch penalty is only interpretable when the reference actually
  # spans both primer sites. Early-Folmer primers such as ZBJ can look
  # catastrophically poor on an alignment built mostly from shorter internal
  # COI fragments. Keep those records for a reference audit, but never present
  # the sparse subset as biological primer coverage.
  pair_reference_scorable_fraction <- mean(pair_scores$pair_scorable)
  pair_reference_suitable <-
    sum(pair_scores$pair_scorable) >= 100L &&
    pair_reference_scorable_fraction >= 0.50
  pair_scores$reference_suitable <- pair_reference_suitable
  pair_scores$reference_suitability <- if (pair_reference_suitable) {
    "usable"
  } else {
    "binding-region incomplete"
  }
  pair_scores$reference_limitation_reason <- if (pair_reference_suitable) {
    ""
  } else {
    paste0(
      "Only ", sum(pair_scores$pair_scorable), "/", nrow(pair_scores),
      " centroids (", round(pair_reference_scorable_fraction * 100, 1),
      "%) span both primer sites; penalties cannot be interpreted as taxon coverage."
    )
  }

  raw_path <- file.path(
    derived_dir,
    paste0("claimed_", tolower(pair_id), "_centroid_scores.csv.gz")
  )
  raw_connection <- gzfile(raw_path, open = "wt")
  write.csv(pair_scores, raw_connection, row.names = FALSE, na = "")
  close(raw_connection)

  overall <- summarize_scores(pair_scores, c("pair_id", "pair_label"))
  overall$reference_suitable <- pair_reference_suitable
  overall$reference_suitability <- unique(pair_scores$reference_suitability)
  overall$reference_limitation_reason <- unique(
    pair_scores$reference_limitation_reason
  )
  by_order <- summarize_scores(pair_scores, c("pair_id", "pair_label", "order"))
  by_order$reference_suitable <- pair_reference_suitable
  by_order$reference_suitability <- unique(pair_scores$reference_suitability)
  by_order$reference_limitation_reason <- unique(
    pair_scores$reference_limitation_reason
  )
  known_order <- pair_scores[
    !is.na(pair_scores$order) & nzchar(trimws(pair_scores$order)),
    ,
    drop = FALSE
  ]
  by_family <- summarize_scores(
    known_order,
    c("pair_id", "pair_label", "order", "family")
  )
  by_family$reference_suitable <- pair_reference_suitable
  by_family$reference_suitability <- unique(pair_scores$reference_suitability)
  known_family <- known_order[
    !is.na(known_order$family) & nzchar(trimws(known_order$family)),
    ,
    drop = FALSE
  ]
  by_genus <- summarize_scores(
    known_family,
    c("pair_id", "pair_label", "order", "family", "genus")
  )
  by_genus$reference_suitable <- pair_reference_suitable
  by_genus$reference_suitability <- unique(pair_scores$reference_suitability)
  known_subfamily <- known_family[
    !is.na(known_family$subfamily) &
      nzchar(trimws(known_family$subfamily)),
    ,
    drop = FALSE
  ]
  by_subfamily <- if (nrow(known_subfamily)) {
    subfamily_result <- summarize_scores(
      known_subfamily,
      c("pair_id", "pair_label", "order", "family", "subfamily")
    )
    subfamily_result$reference_suitable <- pair_reference_suitable
    subfamily_result$reference_suitability <- unique(
      pair_scores$reference_suitability
    )
    subfamily_result
  } else {
    NULL
  }

  summary_overall[[pair_id]] <- overall
  summary_order[[pair_id]] <- by_order
  summary_family[[pair_id]] <- by_family
  summary_genus[[pair_id]] <- by_genus
  summary_subfamily[[pair_id]] <- by_subfamily
  grouped_scores <- assign_phylogeny_group(pair_scores)
  grouped_scores <- grouped_scores[
    !is.na(grouped_scores$display_group),
    ,
    drop = FALSE
  ]
  if (nrow(grouped_scores)) {
    group_summary <- summarize_scores(
      grouped_scores,
      c("pair_id", "pair_label", "display_group_type", "display_group")
    )
    group_notes <- unique(grouped_scores[, c(
      "display_group_type", "display_group", "group_definition", "group_caveat"
    )])
    group_summary <- merge(
      group_summary,
      group_notes,
      by = c("display_group_type", "display_group"),
      all.x = TRUE,
      sort = FALSE
    )
    group_summary$reference_suitable <- pair_reference_suitable
    group_summary$reference_suitability <- unique(
      pair_scores$reference_suitability
    )
    summary_phylogeny_group[[pair_id]] <- group_summary
  }

  message(
    pair_id, ": median pair penalty ",
    round(overall$pair_median_penalty, 2),
    "; both sites scorable in ",
    round(overall$pair_scorable_fraction * 100, 1), "%"
  )
}

bind_rows_base <- function(items) {
  items <- items[!vapply(items, is.null, logical(1))]
  if (!length(items)) return(data.frame())
  do.call(rbind, items)
}
sort_rows <- function(data, columns) {
  if (!nrow(data)) return(data)
  rownames(data) <- NULL
  data[do.call(order, data[columns]), , drop = FALSE]
}

overall_output <- sort_rows(
  bind_rows_base(summary_overall),
  c("pair_label")
)
order_output <- sort_rows(
  bind_rows_base(summary_order),
  c("pair_label", "order")
)
family_output <- sort_rows(
  bind_rows_base(summary_family),
  c("pair_label", "order", "family")
)
genus_output <- sort_rows(
  bind_rows_base(summary_genus),
  c("pair_label", "order", "family", "genus")
)
subfamily_output <- sort_rows(
  bind_rows_base(summary_subfamily),
  c("pair_label", "order", "family", "subfamily")
)
phylogeny_group_output <- sort_rows(
  bind_rows_base(summary_phylogeny_group),
  c("pair_label", "display_group_type", "display_group")
)

write.csv(
  overall_output,
  file.path(derived_dir, "claimed_primer_overall_summary.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  order_output,
  file.path(derived_dir, "claimed_primer_order_summary.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  family_output,
  file.path(derived_dir, "claimed_primer_family_summary.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  genus_output,
  file.path(derived_dir, "claimed_primer_genus_summary.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  subfamily_output,
  file.path(derived_dir, "claimed_primer_subfamily_summary.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  phylogeny_group_output,
  file.path(derived_dir, "claimed_primer_phylogeny_group_summary.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  bind_rows_base(position_manifest),
  file.path(derived_dir, "claimed_primer_scoring_manifest.csv"),
  row.names = FALSE,
  na = ""
)

message("Expanded target-claim summaries written for ", length(pair_ids), " primer pairs.")
