#!/usr/bin/env Rscript

project_root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(project_root, ".Rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(PrimerMiner)
  library(jsonlite)
  library(readxl)
  library(xml2)
})

source(file.path(project_root, "R", "functions.R"))

external_dir <- file.path(project_root, "data", "external", "gurten2026")
derived_dir <- file.path(project_root, "data", "derived")
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)

alignment_path <- file.path(external_dir, "ClusteredReferences.fasta")
bee_lookup_path <- file.path(
  external_dir,
  "supplement_7_bees.xlsx"
)
arthropod_lookup_path <- file.path(
  external_dir,
  "supplement_8_arthropods.xlsx"
)

required_paths <- c(alignment_path, bee_lookup_path, arthropod_lookup_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) {
  stop("Missing Gurten et al. (2026) input(s): ", paste(missing_paths, collapse = ", "))
}

beeprime_forward <- "ATGAATTAATAATGATCANATYTATAAYWC"
beeprime_reverse <- "GGATAWACWGTTCAWCCWGTWCC"

message("Reading Gurten et al. (2026) 99.5% centroid alignment")
aligned_sequences <- parse_fasta(alignment_path)
headers <- names(aligned_sequences)
if (length(aligned_sequences) != 67352L) {
  warning(
    "Expected 67,352 centroid sequences from Supplement 5; found ",
    length(aligned_sequences), "."
  )
}
alignment_widths <- nchar(aligned_sequences)
if (length(unique(alignment_widths)) != 1L) {
  stop("Supplement 5 is not a rectangular alignment.")
}

accession <- sub(" .*", "", headers)
description <- sub("^[^ ]+ +", "", headers)
genus <- sub(" .*", "", description)
genus <- gsub("^\\[|\\]$", "", genus)
cluster_size <- ifelse(
  grepl(";size=[0-9]+;", headers),
  as.integer(sub(".*;size=([0-9]+);.*", "\\1", headers)),
  1L
)

organism_label <- sub(";size=[0-9]+;.*$", "", description)
organism_label <- sub(
  paste0(
    " (voucher|isolate|clone|haplotype|specimen|strain|sample|",
    "cytochrome|mitochondrial|COI|COX|cox1|coxI|tRNA).*$"
  ),
  "",
  organism_label,
  ignore.case = TRUE
)
organism_label <- trimws(organism_label)

bee_lookup <- as.data.frame(read_excel(
  bee_lookup_path,
  sheet = "BeesInReference"
), stringsAsFactors = FALSE)
arthropod_lookup <- as.data.frame(read_excel(
  arthropod_lookup_path,
  sheet = "SpeciesInReference"
), stringsAsFactors = FALSE)

unique_rank_for_genus <- function(values) {
  values <- unique(trimws(values[!is.na(values) & nzchar(trimws(values))]))
  if (length(values) == 1L) values else NA_character_
}

lookup_groups <- split(arthropod_lookup, arthropod_lookup$genus)
genus_lookup <- do.call(rbind, lapply(lookup_groups, function(x) {
  data.frame(
    genus = x$genus[1],
    phylum = unique_rank_for_genus(x$phylum),
    class = unique_rank_for_genus(x$class),
    order = unique_rank_for_genus(x$order),
    family = unique_rank_for_genus(x$family),
    genus_mapping_ambiguous = length(unique(
      paste(x$phylum, x$class, x$order, x$family, sep = "|")
    )) > 1L,
    stringsAsFactors = FALSE
  )
}))
rownames(genus_lookup) <- NULL

taxonomy_match <- match(genus, genus_lookup$genus)
centroid_metadata <- data.frame(
  accession = accession,
  header = headers,
  organism_label = organism_label,
  genus = genus,
  cluster_size = cluster_size,
  phylum = genus_lookup$phylum[taxonomy_match],
  class = genus_lookup$class[taxonomy_match],
  order = genus_lookup$order[taxonomy_match],
  family = genus_lookup$family[taxonomy_match],
  genus_mapping_ambiguous = genus_lookup$genus_mapping_ambiguous[taxonomy_match],
  stringsAsFactors = FALSE
)
centroid_metadata$genus_mapping_ambiguous[
  is.na(centroid_metadata$genus_mapping_ambiguous)
] <- FALSE

bee_match <- match(centroid_metadata$genus, bee_lookup$genus)
centroid_metadata$is_bee <- !is.na(bee_match)
centroid_metadata$phylum[centroid_metadata$is_bee] <- "Arthropoda"
centroid_metadata$class[centroid_metadata$is_bee] <- "Insecta"
centroid_metadata$order[centroid_metadata$is_bee] <- "Hymenoptera"
centroid_metadata$family[centroid_metadata$is_bee] <- bee_lookup$family[
  bee_match[centroid_metadata$is_bee]
]

read_esummary_map <- function(paths) {
  rows <- list()
  row_index <- 1L
  for (path in paths) {
    payload <- jsonlite::fromJSON(path, simplifyVector = FALSE)
    ids <- payload$result$uids
    for (uid in ids) {
      record <- payload$result[[uid]]
      rows[[row_index]] <- data.frame(
        accession = record$accessionversion,
        taxid = as.character(record$taxid),
        ncbi_title = record$title,
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }
  result <- do.call(rbind, rows)
  result[!duplicated(result$accession), , drop = FALSE]
}

esummary_paths <- sort(list.files(
  external_dir,
  pattern = "^bee_centroid_esummary_batch_[0-9]+\\.json$",
  full.names = TRUE
))
if (length(esummary_paths)) {
  esummary_map <- read_esummary_map(esummary_paths)
  summary_match <- match(centroid_metadata$accession, esummary_map$accession)
  centroid_metadata$taxid <- esummary_map$taxid[summary_match]
  centroid_metadata$ncbi_title <- esummary_map$ncbi_title[summary_match]
} else {
  centroid_metadata$taxid <- NA_character_
  centroid_metadata$ncbi_title <- NA_character_
}

taxonomy_value <- function(node, wanted_rank) {
  own_rank <- xml_text(xml_find_first(node, "./Rank"))
  if (identical(own_rank, wanted_rank)) {
    return(xml_text(xml_find_first(node, "./ScientificName")))
  }
  match_node <- xml_find_first(
    node,
    paste0("./LineageEx/Taxon[Rank='", wanted_rank, "']/ScientificName")
  )
  if (inherits(match_node, "xml_missing")) NA_character_ else xml_text(match_node)
}

read_taxonomy_map <- function(paths) {
  rows <- list()
  row_index <- 1L
  for (path in paths) {
    doc <- read_xml(path)
    taxa <- xml_find_all(doc, "/TaxaSet/Taxon")
    for (node in taxa) {
      rows[[row_index]] <- data.frame(
        taxid = xml_text(xml_find_first(node, "./TaxId")),
        scientific_name = xml_text(xml_find_first(node, "./ScientificName")),
        subfamily = taxonomy_value(node, "subfamily"),
        tribe = taxonomy_value(node, "tribe"),
        genus_ncbi = taxonomy_value(node, "genus"),
        family_ncbi = taxonomy_value(node, "family"),
        order_ncbi = taxonomy_value(node, "order"),
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }
  result <- do.call(rbind, rows)
  result[!duplicated(result$taxid), , drop = FALSE]
}

taxonomy_paths <- sort(list.files(
  external_dir,
  pattern = "^bee_taxonomy_batch_[0-9]+\\.xml$",
  full.names = TRUE
))
if (length(taxonomy_paths)) {
  taxonomy_map <- read_taxonomy_map(taxonomy_paths)
  taxonomy_match <- match(centroid_metadata$taxid, taxonomy_map$taxid)
  centroid_metadata$subfamily <- taxonomy_map$subfamily[taxonomy_match]
  centroid_metadata$tribe <- taxonomy_map$tribe[taxonomy_match]
  centroid_metadata$scientific_name_ncbi <- taxonomy_map$scientific_name[
    taxonomy_match
  ]
} else {
  centroid_metadata$subfamily <- NA_character_
  centroid_metadata$tribe <- NA_character_
  centroid_metadata$scientific_name_ncbi <- NA_character_
}

centroid_metadata$subfamily[
  centroid_metadata$is_bee &
    (is.na(centroid_metadata$subfamily) | !nzchar(centroid_metadata$subfamily))
] <- "Subfamily unresolved at NCBI"
centroid_metadata$tribe[
  centroid_metadata$is_bee &
    (is.na(centroid_metadata$tribe) | !nzchar(centroid_metadata$tribe))
] <- "Tribe unresolved at NCBI"

lineage_group <- rep("Other arthropod", nrow(centroid_metadata))
lineage_group[centroid_metadata$order == "Hymenoptera"] <- "Other Hymenoptera"
lineage_group[centroid_metadata$family == "Formicidae"] <- "Ant"
lineage_group[centroid_metadata$family %in% c(
  "Argidae", "Cimbicidae", "Diprionidae", "Pergidae", "Tenthredinidae",
  "Pamphiliidae", "Cephidae", "Siricidae", "Xiphydriidae"
)] <- "Sawfly or woodwasp"
lineage_group[centroid_metadata$family %in% c(
  "Chrysididae", "Crabronidae", "Mutillidae", "Pompilidae", "Scoliidae",
  "Sphecidae", "Tiphiidae", "Vespidae"
)] <- "Aculeate wasp"
lineage_group[centroid_metadata$is_bee] <- "Bee"
centroid_metadata$lineage_group <- lineage_group

reference_accession <- "OM794648.1"
reference_index <- match(reference_accession, centroid_metadata$accession)
if (is.na(reference_index)) {
  stop("Reference bee centroid ", reference_accession, " is absent.")
}
reference_aligned <- toupper(aligned_sequences[[reference_index]])
reference_ungapped <- gsub("-", "", reference_aligned, fixed = TRUE)
alignment_columns <- which(split_bases(reference_aligned) != "-")

forward_candidates <- primer_candidate_table(
  reference_ungapped,
  beeprime_forward,
  "forward"
)
reverse_candidates <- primer_candidate_table(
  reference_ungapped,
  beeprime_reverse,
  "reverse"
)
forward_best <- forward_candidates[1, , drop = FALSE]
reverse_best <- reverse_candidates[1, , drop = FALSE]
forward_start <- alignment_columns[forward_best$start]
forward_end <- alignment_columns[forward_best$end]
reverse_start <- alignment_columns[reverse_best$start]
reverse_end <- alignment_columns[reverse_best$end]

if (!identical(c(forward_start, forward_end, reverse_start, reverse_end),
  c(172L, 201L, 401L, 423L))) {
  stop(
    "Unexpected BeePrime columns in Gurten et al. alignment: ",
    paste(c(forward_start, forward_end, reverse_start, reverse_end), collapse = ", ")
  )
}

evaluate_primer_panel <- function(sequence, start, end, forward, primer_name) {
  message("Evaluating ", primer_name, " on ", length(aligned_sequences), " centroids")
  evaluated <- suppressMessages(PrimerMiner::evaluate_primer(
    alignment_imp = alignment_path,
    primer_sequ = clean_sequence(sequence),
    start = start,
    stop = end,
    forward = forward,
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

  data.frame(
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
}

forward_scores <- evaluate_primer_panel(
  beeprime_forward,
  forward_start,
  forward_end,
  TRUE,
  "BeePrimeF"
)
reverse_scores <- evaluate_primer_panel(
  beeprime_reverse,
  reverse_start,
  reverse_end,
  FALSE,
  "BeePrimeR"
)

names(forward_scores)[-1] <- paste0("forward_", names(forward_scores)[-1])
names(reverse_scores)[-1] <- paste0("reverse_", names(reverse_scores)[-1])
pair_scores <- merge(
  forward_scores,
  reverse_scores,
  by = "accession",
  all = TRUE,
  sort = FALSE
)
metadata_match <- match(pair_scores$accession, centroid_metadata$accession)
pair_scores <- cbind(
  pair_scores,
  centroid_metadata[
    metadata_match,
    setdiff(names(centroid_metadata), "accession"),
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
      n_pair_scorable = sum(scorable),
      pair_scorable_fraction = mean(scorable),
      forward_median_penalty = median_safe(x$forward_penalty),
      reverse_median_penalty = median_safe(x$reverse_penalty),
      pair_median_penalty = median_safe(x$pair_penalty),
      pair_mean_penalty = mean_safe(x$pair_penalty),
      pair_p90_penalty = quantile_safe(x$pair_penalty, 0.9),
      weighted_mean_pair_penalty = if (sum(scorable)) {
        weighted.mean(
          x$pair_penalty[scorable],
          w = pmax(x$cluster_size[scorable], 1L),
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
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

order_summary <- summarize_scores(pair_scores, "order")
family_summary <- summarize_scores(pair_scores, c("order", "family"))
hymenoptera_scores <- pair_scores[
  pair_scores$order %in% "Hymenoptera",
  ,
  drop = FALSE
]
hymenoptera_family_summary <- summarize_scores(
  hymenoptera_scores,
  c("lineage_group", "family")
)
hymenoptera_lineage_summary <- summarize_scores(
  hymenoptera_scores,
  "lineage_group"
)
bee_scores <- pair_scores[pair_scores$is_bee %in% TRUE, , drop = FALSE]
bee_family_summary <- summarize_scores(bee_scores, "family")
bee_subfamily_summary <- summarize_scores(bee_scores, c("family", "subfamily"))
bee_genus_summary <- summarize_scores(
  bee_scores,
  c("family", "subfamily", "genus")
)
scored_bee_genera_count <- length(unique(bee_scores$genus))

# Supplement 7 lists 42 bee genera, while only 30 genus names can be recovered
# directly from the centroid FASTA headers. Complete the genus output against
# the authors' catalog so those 12 mapping gaps remain visible instead of being
# silently dropped.
bee_genus_catalog <- unique(bee_lookup[, c("family", "genus", "species number")])
names(bee_genus_catalog)[names(bee_genus_catalog) == "species number"] <-
  "species_listed_in_supplement"
bee_genus_summary <- merge(
  bee_genus_catalog,
  bee_genus_summary,
  by = c("family", "genus"),
  all.x = TRUE,
  sort = FALSE
)
zero_when_unmapped <- c(
  "n_centroids", "source_records_represented", "n_pair_scorable"
)
for (column in zero_when_unmapped) {
  bee_genus_summary[[column]][is.na(bee_genus_summary[[column]])] <- 0
}
bee_genus_summary$subfamily[
  is.na(bee_genus_summary$subfamily) | !nzchar(bee_genus_summary$subfamily)
] <- "subfamily unresolved"

sort_summary <- function(x, columns) {
  x[do.call(order, x[columns]), , drop = FALSE]
}

write.csv(
  centroid_metadata,
  file.path(derived_dir, "gurten2026_centroid_taxonomy.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  pair_scores,
  file.path(derived_dir, "beeprime_centroid_scores_gurten2026.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  hymenoptera_scores,
  file.path(derived_dir, "beeprime_hymenoptera_centroid_scores.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  sort_summary(order_summary, "order"),
  file.path(derived_dir, "beeprime_order_summary_gurten2026.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  sort_summary(family_summary, c("order", "family")),
  file.path(derived_dir, "beeprime_family_summary_gurten2026.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  sort_summary(hymenoptera_lineage_summary, "lineage_group"),
  file.path(derived_dir, "beeprime_hymenoptera_lineage_summary.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  sort_summary(hymenoptera_family_summary, c("lineage_group", "family")),
  file.path(derived_dir, "beeprime_hymenoptera_family_summary.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  sort_summary(bee_family_summary, "family"),
  file.path(derived_dir, "beeprime_bee_family_summary.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  sort_summary(bee_subfamily_summary, c("family", "subfamily")),
  file.path(derived_dir, "beeprime_bee_subfamily_summary.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  sort_summary(bee_genus_summary, c("family", "subfamily", "genus")),
  file.path(derived_dir, "beeprime_bee_genus_summary.csv"),
  row.names = FALSE,
  na = ""
)

manifest <- data.frame(
  field = c(
    "source_article",
    "source_doi",
    "source_supplement",
    "source_supplement_doi",
    "source_access_date",
    "source_alignment_centroids",
    "source_alignment_width",
    "underlying_ncbi_sequences_reported_by_authors",
    "authors_ncbi_access_date",
    "clustering_identity",
    "bee_centroids_mapped",
    "bee_genera_in_supplement",
    "bee_genera_with_scored_centroid_labels",
    "bee_families_in_supplement",
    "bee_species_listed_in_supplement",
    "beeprime_forward_alignment_columns",
    "beeprime_reverse_alignment_columns",
    "penalty_model",
    "interpretation"
  ),
  value = c(
    "Gurten et al. 2026",
    "10.3897/mbmg.10.183708",
    "Supplementary material 5: ClusteredReferences.fasta",
    "10.3897/mbmg.10.183708.suppl5",
    as.character(Sys.Date()),
    length(aligned_sequences),
    unique(alignment_widths),
    316254,
    "May 2023",
    0.995,
    nrow(bee_scores),
    length(unique(bee_lookup$genus)),
    scored_bee_genera_count,
    length(unique(bee_lookup$family)),
    sum(bee_lookup[["species number"]]),
    paste0(forward_start, "-", forward_end),
    paste0(reverse_start, "-", reverse_end),
    "PrimerMiner 0.22 Position_v1 + Type_v1 + adjacency 2",
    paste(
      "Relative in-silico mismatch penalties; lower is a closer predicted fit.",
      "Not a PCR probability or universal amplification threshold."
    )
  ),
  stringsAsFactors = FALSE
)
write.csv(
  manifest,
  file.path(derived_dir, "beeprime_reference_manifest.csv"),
  row.names = FALSE,
  na = ""
)

message("\nBeePrime expanded reference summary")
message("  all centroids: ", nrow(pair_scores))
message("  Hymenoptera centroids: ", nrow(hymenoptera_scores))
message("  bee centroids: ", nrow(bee_scores))
message("  bee families: ", length(unique(bee_scores$family)))
message(
  "  bee pair median penalty: ",
  round(median_safe(bee_scores$pair_penalty), 2)
)
message(
  "  bee pair scorable fraction: ",
  round(mean(bee_scores$pair_scorable), 4)
)
print(bee_family_summary, row.names = FALSE)
