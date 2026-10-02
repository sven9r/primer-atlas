#!/usr/bin/env Rscript

root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(root, ".Rlib"), .libPaths()))
suppressPackageStartupMessages({
  library(PrimerMiner)
  library(rentrez)
  library(xml2)
})
source("R/functions.R")
source("R/data_layer.R")

derived_dir <- "data/derived"
provenance_dir <- "data/provenance"
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(provenance_dir, recursive = TRUE, showWarnings = FALSE)

publication_fasta <- "data/external/gurten2026/ClusteredReferences.fasta"
publication_scores_path <- file.path(derived_dir, "beeprime_hymenoptera_centroid_scores.csv")
current_hymenoptera_fasta <- "data/primerminer/raw/Hymenoptera/Hymenoptera_GB.fasta"
alignment_output <- file.path(derived_dir, "beeprime_expanded_alignment.fasta")
added_fasta_output <- file.path(derived_dir, "beeprime_expanded_added_sequences.fasta")
scores_output <- file.path(derived_dir, "beeprime_expanded_scores.csv.gz")
geography_output <- file.path(derived_dir, "beeprime_expanded_geography.csv")
taxonomy_audit_output <- file.path(provenance_dir, "beeprime_expanded_taxonomy_audit.csv")
growth_qa_output <- file.path(provenance_dir, "beeprime_expanded_growth_qa.csv")

required_inputs <- c(publication_fasta, publication_scores_path, current_hymenoptera_fasta)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop("Missing BeePrime expanded-reference input(s): ", paste(missing_inputs, collapse = ", "))
}

accession_from_header <- function(header) {
  token <- sub(" .*", "", sub("^>", "", header))
  token <- sub("^_R_", "", token)
  sub("[|].*$", "", token)
}

write_fasta <- function(sequences, path) {
  lines <- unlist(lapply(seq_along(sequences), function(i) {
    sequence <- toupper(gsub("[[:space:]]", "", unname(sequences[[i]])))
    starts <- seq.int(1L, nchar(sequence), by = 80L)
    c(paste0(">", names(sequences)[i]),
      substring(sequence, starts, pmin(starts + 79L, nchar(sequence))))
  }), use.names = FALSE)
  writeLines(lines, path, useBytes = TRUE)
}

checked_manifest_artifact_path <- function(artifact) {
  if (is.null(artifact) || is.null(artifact$sha256)) {
    stop("Previous expanded BeePrime release is missing a checksummed artifact.")
  }
  path <- artifact$local_path
  if (!is.null(path) && !grepl("^/", path)) path <- file.path(root, path)
  if (is.null(path) || !file.exists(path)) {
    if (is.null(artifact$url) || !nzchar(artifact$url)) {
      stop("Previous expanded BeePrime artifact has no downloadable URL.")
    }
    path <- atlas_cached_artifact(artifact$url, artifact$sha256)
  }
  actual <- digest::digest(file = path, algo = "sha256", serialize = FALSE)
  if (!identical(tolower(actual), tolower(artifact$sha256))) {
    stop("Previous expanded BeePrime artifact checksum failed: ", basename(path))
  }
  path
}

previous_state <- tryCatch(atlas_read_manifest("COI", root), error = identity)
if (inherits(previous_state, "error") && nzchar(Sys.getenv("ATLAS_MANIFEST_BASE_URL", ""))) {
  stop("Could not read the previous COI release; refusing to reset the BeePrime sequence panel: ",
       conditionMessage(previous_state))
}
if (!inherits(previous_state, "error") &&
    nzchar(Sys.getenv("ATLAS_MANIFEST_BASE_URL", "")) &&
    isTRUE(previous_state$stale)) {
  stop("The remote COI manifest could not be verified; refusing to reset the BeePrime sequence panel.")
}
previous_descriptor <- if (!inherits(previous_state, "error")) {
  previous_state$manifest$reference_panels$beeprime_expanded
} else NULL

publication_scores <- read.csv(publication_scores_path, stringsAsFactors = FALSE,
                               check.names = FALSE)
publication_bees <- publication_scores[publication_scores$is_bee %in% TRUE, , drop = FALSE]
publication_accessions <- as.character(publication_bees$accession)
if (nrow(publication_bees) != 590L || anyDuplicated(publication_accessions)) {
  stop("The Gurten publication BeePrime baseline must contain 590 unique bee accessions.")
}

previous_scores <- data.frame()
previous_added_sequences <- character()
if (!is.null(previous_descriptor)) {
  if (!identical(previous_descriptor$pair_id, "BEEPRIME") ||
      !identical(previous_descriptor$taxonomic_scope, "bees") ||
      !identical(previous_descriptor$publication_baseline, "gurten2026_suppl5")) {
    stop("The previous release's BeePrime panel has an unexpected identity or taxonomic scope.")
  }
  previous_scores <- as.data.frame(atlas_read_artifact(previous_descriptor$sequence_scores, root),
                                   stringsAsFactors = FALSE)
  sequence_path <- checked_manifest_artifact_path(previous_descriptor$added_sequences_fasta)
  previous_added_sequences <- parse_fasta(sequence_path)
  if (!nrow(previous_scores) || !length(previous_added_sequences)) {
    stop("The previous expanded BeePrime release is incomplete; refusing to replace it.")
  }
  if (!setequal(previous_scores$accession[previous_scores$reference_origin == "publication"],
                publication_accessions)) {
    stop("The previous expanded BeePrime release does not preserve the exact publication baseline.")
  }
  previous_score_additions <- as.character(
    previous_scores$accession[previous_scores$reference_origin == "added"]
  )
  previous_fasta_additions <- vapply(names(previous_added_sequences), accession_from_header, character(1))
  if (!setequal(previous_score_additions, previous_fasta_additions)) {
    stop("The previous BeePrime added-sequence FASTA does not match its score artifact.")
  }
}

publication_alignment <- parse_fasta(publication_fasta)
publication_headers <- names(publication_alignment)
publication_alignment_accessions <- vapply(publication_headers, accession_from_header, character(1))
baseline_match <- match(publication_accessions, publication_alignment_accessions)
if (anyNA(baseline_match)) {
  stop("One or more publication BeePrime accessions are absent from the publication alignment.")
}
baseline_alignment <- publication_alignment[baseline_match]
baseline_widths <- unique(nchar(baseline_alignment))
if (length(baseline_widths) != 1L || baseline_widths < 423L) {
  stop("The publication BeePrime baseline must retain a rectangular alignment spanning both primer sites.")
}

current_sequences <- parse_fasta(current_hymenoptera_fasta)
current_accessions <- vapply(names(current_sequences), accession_from_header, character(1))
if (any(!nzchar(current_accessions))) stop("Hymenoptera input includes a record without an accession.")
current_sequences <- current_sequences[!duplicated(current_accessions)]
current_accessions <- current_accessions[!duplicated(current_accessions)]

known_accessions <- publication_accessions
if (nrow(previous_scores)) {
  known_accessions <- unique(c(known_accessions, previous_scores$accession))
}
candidate_index <- which(!current_accessions %in% known_accessions)
candidate_sequences <- current_sequences[candidate_index]
candidate_accessions <- current_accessions[candidate_index]

taxonomy_for_candidates <- function(accessions) {
  empty <- data.frame(
    accession = character(), taxid = character(), scientific_name_ncbi = character(),
    order = character(), family = character(), subfamily = character(), genus = character(),
    is_bee = logical(), taxonomy_status = character(),
    stringsAsFactors = FALSE
  )
  if (!length(accessions)) return(empty)

  batch <- function(x, size = 75L) split(x, ceiling(seq_along(x) / size))
  summary_rows <- list()
  summary_index <- 1L
  field_or_na <- function(record, field) {
    value <- record[[field]]
    if (is.null(value) || !length(value)) return(NA_character_)
    as.character(value[[1]])
  }
  for (ids in batch(accessions)) {
    result <- tryCatch(rentrez::entrez_summary(db = "nuccore", id = ids, retmode = "json"),
                       error = function(e) stop("NCBI nuccore summary failed: ", conditionMessage(e)))
    records <- if (inherits(result, "esummary_list")) unclass(result) else list(result)
    for (record in records) {
      summary_rows[[summary_index]] <- data.frame(
        accession = field_or_na(record, "accessionversion"),
        taxid = field_or_na(record, "taxid"),
        scientific_name_ncbi = field_or_na(record, "organism"),
        stringsAsFactors = FALSE
      )
      summary_index <- summary_index + 1L
    }
    Sys.sleep(if (nzchar(Sys.getenv("NCBI_API_KEY", ""))) 0.12 else 0.36)
  }
  summaries <- if (length(summary_rows)) do.call(rbind, summary_rows) else empty[0, ]
  summaries <- summaries[!duplicated(summaries$accession), , drop = FALSE]
  summaries <- summaries[match(accessions, summaries$accession), , drop = FALSE]
  summaries$taxonomy_status <- ifelse(is.na(summaries$taxid) | !nzchar(summaries$taxid),
                                      "unresolved_taxid", "unresolved_lineage")
  summaries$is_bee <- FALSE
  summaries$order <- summaries$family <- summaries$subfamily <- summaries$genus <- NA_character_

  taxids <- unique(summaries$taxid[!is.na(summaries$taxid) & nzchar(summaries$taxid)])
  taxonomy_rows <- list()
  taxonomy_index <- 1L
  for (ids in batch(taxids)) {
    xml_text_value <- tryCatch(
      rentrez::entrez_fetch(db = "taxonomy", id = ids, rettype = "xml", retmode = "xml"),
      error = function(e) stop("NCBI taxonomy fetch failed: ", conditionMessage(e))
    )
    document <- tryCatch(xml2::read_xml(xml_text_value), error = function(e) {
      stop("NCBI taxonomy returned invalid XML: ", conditionMessage(e))
    })
    nodes <- xml2::xml_find_all(document, "/TaxaSet/Taxon")
    for (node in nodes) {
      lineage <- xml2::xml_find_all(node, "./LineageEx/Taxon")
      ranks <- vapply(lineage, function(x) xml2::xml_text(xml2::xml_find_first(x, "./Rank")), character(1))
      names <- vapply(lineage, function(x) xml2::xml_text(xml2::xml_find_first(x, "./ScientificName")), character(1))
      ids_in_lineage <- vapply(lineage, function(x) xml2::xml_text(xml2::xml_find_first(x, "./TaxId")), character(1))
      own_rank <- xml2::xml_text(xml2::xml_find_first(node, "./Rank"))
      own_name <- xml2::xml_text(xml2::xml_find_first(node, "./ScientificName"))
      own_id <- xml2::xml_text(xml2::xml_find_first(node, "./TaxId"))
      ranks <- c(ranks, own_rank)
      names <- c(names, own_name)
      ids_in_lineage <- c(ids_in_lineage, own_id)
      rank_value <- function(rank) {
        found <- which(ranks == rank)
        if (length(found)) names[tail(found, 1)] else NA_character_
      }
      taxid <- own_id
      taxonomy_rows[[taxonomy_index]] <- data.frame(
        taxid = taxid,
        scientific_name_ncbi = own_name,
        order = rank_value("order"), family = rank_value("family"),
        subfamily = rank_value("subfamily"), genus = rank_value("genus"),
        is_bee = "3042114" %in% ids_in_lineage,
        stringsAsFactors = FALSE
      )
      taxonomy_index <- taxonomy_index + 1L
    }
    Sys.sleep(if (nzchar(Sys.getenv("NCBI_API_KEY", ""))) 0.12 else 0.36)
  }
  if (length(taxonomy_rows)) {
    taxonomy <- do.call(rbind, taxonomy_rows)
    taxonomy <- taxonomy[!duplicated(taxonomy$taxid), , drop = FALSE]
    taxonomy_match <- match(summaries$taxid, taxonomy$taxid)
    for (column in c("order", "family", "subfamily", "genus", "is_bee")) {
      summaries[[column]] <- taxonomy[[column]][taxonomy_match]
    }
    summaries$is_bee[is.na(summaries$is_bee)] <- FALSE
    summaries$taxonomy_status <- ifelse(
      is.na(taxonomy_match), "unresolved_lineage",
      ifelse(summaries$is_bee, "included_bee", "not_bee")
    )
  }
  summaries
}

taxonomy_audit <- taxonomy_for_candidates(candidate_accessions)
write.csv(taxonomy_audit, taxonomy_audit_output, row.names = FALSE, na = "")

new_bee_accessions <- taxonomy_audit$accession[taxonomy_audit$taxonomy_status == "included_bee"]
new_bee_sequences <- candidate_sequences[match(new_bee_accessions, candidate_accessions)]
previous_added_accessions <- if (length(previous_added_sequences)) {
  vapply(names(previous_added_sequences), accession_from_header, character(1))
} else character()
new_keep <- !new_bee_accessions %in% previous_added_accessions
new_bee_accessions <- new_bee_accessions[new_keep]
new_bee_sequences <- new_bee_sequences[new_keep]

added_sequences <- c(previous_added_sequences, new_bee_sequences)
added_accessions <- if (length(added_sequences)) {
  vapply(names(added_sequences), accession_from_header, character(1))
} else character()
keep_unique <- !duplicated(added_accessions)
added_sequences <- added_sequences[keep_unique]
added_accessions <- added_accessions[keep_unique]
if (any(added_accessions %in% publication_accessions)) {
  stop("An added sequence duplicates an accession in the publication baseline.")
}

qa <- data.frame(
  baseline_publication_bees = length(publication_accessions),
  previously_retained_added_bees = length(previous_added_accessions),
  current_hymenoptera_records = length(current_sequences),
  new_records_taxonomy_checked = length(candidate_accessions),
  newly_identified_bee_records = length(new_bee_accessions),
  retained_added_bee_records = length(added_accessions),
  current_taxonomy_unresolved = sum(taxonomy_audit$taxonomy_status %in% c("unresolved_taxid", "unresolved_lineage")),
  release_id = Sys.getenv("RELEASE_ID", format(Sys.Date(), "%Y-%m-%d")),
  stringsAsFactors = FALSE
)
write.csv(qa, growth_qa_output, row.names = FALSE, na = "")
print(qa, row.names = FALSE)

if (!length(added_sequences)) {
  unlink(c(alignment_output, added_fasta_output, scores_output, geography_output))
  message("No new bee sequences are available yet; the publication alignment remains the only selectable BeePrime reference.")
  quit(save = "no", status = 0L)
}

names(added_sequences) <- vapply(seq_along(added_sequences), function(i) {
  paste(added_accessions[i], sub("^[^ ]+\\s*", "", names(added_sequences)[i]))
}, character(1))
write_fasta(added_sequences, added_fasta_output)

seed_accession <- "OM794648.1"
seed_index <- match(seed_accession, publication_alignment_accessions)
if (is.na(seed_index)) stop("The fixed BeePrime coordinate sequence is absent from the publication alignment.")
seed_alignment <- publication_alignment[seed_index]
seed_path <- tempfile(fileext = ".fasta")
on.exit(unlink(seed_path), add = TRUE)
added_path <- tempfile(fileext = ".fasta")
on.exit(unlink(added_path), add = TRUE)
new_alignment_path <- tempfile(fileext = ".fasta")
on.exit(unlink(new_alignment_path), add = TRUE)
write_fasta(seed_alignment, seed_path)
writeLines(readLines(added_fasta_output, warn = FALSE), added_path, useBytes = TRUE)
mafft <- Sys.which("mafft")
if (!nzchar(mafft)) stop("MAFFT was not found on PATH.")
mafft_log <- tempfile(fileext = ".log")
status <- system2(mafft, c("--quiet", "--adjustdirectionaccurately", "--keeplength", "--addfragments",
                           shQuote(added_path), shQuote(seed_path)),
                   stdout = new_alignment_path, stderr = mafft_log)
if (!identical(status, 0L)) stop("MAFFT failed while aligning the expanded BeePrime panel: ", mafft_log)
new_alignment <- parse_fasta(new_alignment_path)
new_alignment_accessions <- vapply(names(new_alignment), accession_from_header, character(1))
seed_output_index <- match(seed_accession, new_alignment_accessions)
if (is.na(seed_output_index) ||
    !identical(toupper(unname(new_alignment[[seed_output_index]])), toupper(unname(seed_alignment)))) {
  stop("MAFFT changed the fixed BeePrime coordinate sequence; refusing to alter publication coordinates.")
}
new_alignment_index <- match(added_accessions, new_alignment_accessions)
if (anyNA(new_alignment_index)) stop("MAFFT lost one or more newly added BeePrime accessions.")
names(new_alignment) <- vapply(names(new_alignment), accession_from_header, character(1))
expanded_alignment <- c(baseline_alignment, new_alignment[new_alignment_index])
expanded_alignment_accessions <- c(publication_accessions, added_accessions)
if (anyDuplicated(expanded_alignment_accessions) ||
    length(expanded_alignment) != length(expanded_alignment_accessions) ||
    any(nchar(expanded_alignment) != baseline_widths)) {
  stop("Expanded BeePrime alignment does not retain the publication accessions and alignment width.")
}
expanded_headers <- names(expanded_alignment)
expanded_accessions <- vapply(expanded_headers, accession_from_header, character(1))
expanded_baseline <- expanded_alignment[match(publication_accessions, expanded_accessions)]
if (anyDuplicated(expanded_accessions) ||
    !setequal(expanded_accessions, c(publication_accessions, added_accessions)) ||
    !identical(toupper(unname(expanded_baseline)), toupper(unname(baseline_alignment)))) {
  stop("MAFFT changed a publication baseline sequence; refusing to publish altered study evidence.")
}
write_fasta(expanded_alignment, alignment_output)

score_primer <- function(sequence, start, stop, forward) {
  evaluated <- suppressMessages(PrimerMiner::evaluate_primer(
    alignment_imp = alignment_output,
    primer_sequ = clean_sequence(sequence),
    start = start, stop = stop, forward = forward,
    gap_NA = TRUE, N_NA = TRUE, mm_position = "Position_v1",
    mm_type = "Type_v1", adjacent = 2, sequ_names = TRUE
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
    accession = vapply(evaluated$Template, accession_from_header, character(1)),
    binding_sequence_primer_oriented = evaluated$sequ,
    penalty = evaluated$sum, scorable = !is.na(evaluated$sum),
    mismatch_count = rowSums(mismatch_matrix),
    terminal_3_mismatches = terminal_count(3L),
    stringsAsFactors = FALSE
  )
}

forward <- score_primer("ATGAATTAATAATGATCANATYTATAAYWC", 172L, 201L, TRUE)
reverse <- score_primer("GGATAWACWGTTCAWCCWGTWCC", 401L, 423L, FALSE)
names(forward)[-1] <- paste0("forward_", names(forward)[-1])
names(reverse)[-1] <- paste0("reverse_", names(reverse)[-1])
scores <- merge(forward, reverse, by = "accession", sort = FALSE)
if (nrow(scores) != length(c(publication_accessions, added_accessions)) || anyDuplicated(scores$accession)) {
  stop("PrimerMiner output did not map one score row to each BeePrime alignment sequence.")
}

baseline_metadata <- data.frame(
  accession = publication_bees$accession,
  order = "Hymenoptera", family = publication_bees$family,
  subfamily = publication_bees$subfamily, genus = publication_bees$genus,
  is_bee = TRUE, cluster_size = publication_bees$cluster_size,
  reference_origin = "publication", organism_label = publication_bees$organism_label,
  scientific_name_ncbi = publication_bees$scientific_name_ncbi,
  taxid = publication_bees$taxid,
  taxonomy_source = "Gurten et al. Supplement 8; NCBI taxonomy reconciliation",
  stringsAsFactors = FALSE
)

added_metadata <- data.frame(
  accession = added_accessions, order = "Hymenoptera",
  family = NA_character_, subfamily = NA_character_, genus = NA_character_,
  is_bee = TRUE, cluster_size = 1L, reference_origin = "added",
  organism_label = NA_character_, scientific_name_ncbi = NA_character_,
  taxid = NA_character_,
  taxonomy_source = "NCBI Taxonomy lineage containing Anthophila taxid 3042114",
  stringsAsFactors = FALSE
)
new_metadata_match <- match(added_accessions, taxonomy_audit$accession)
new_metadata_rows <- !is.na(new_metadata_match)
if (any(new_metadata_rows)) {
  new_taxonomy <- taxonomy_audit[new_metadata_match[new_metadata_rows], , drop = FALSE]
  added_metadata$family[new_metadata_rows] <- new_taxonomy$family
  added_metadata$subfamily[new_metadata_rows] <- new_taxonomy$subfamily
  added_metadata$genus[new_metadata_rows] <- new_taxonomy$genus
  added_metadata$scientific_name_ncbi[new_metadata_rows] <- new_taxonomy$scientific_name_ncbi
  added_metadata$taxid[new_metadata_rows] <- new_taxonomy$taxid
  added_metadata$organism_label[new_metadata_rows] <- new_taxonomy$scientific_name_ncbi
  if (any(!new_taxonomy$is_bee %in% TRUE)) {
    stop("A non-bee taxon entered the BeePrime expanded reference.")
  }
}
reused_rows <- !new_metadata_rows
if (any(reused_rows)) {
  previous_rows <- previous_scores[
    previous_scores$reference_origin == "added", , drop = FALSE
  ]
  previous_match <- match(added_accessions[reused_rows], previous_rows$accession)
  if (anyNA(previous_match)) {
    stop("Prior retained BeePrime sequences are missing their taxonomy metadata.")
  }
  previous_meta <- previous_rows[previous_match, , drop = FALSE]
  for (column in c("family", "subfamily", "genus", "organism_label",
                   "scientific_name_ncbi", "taxid", "taxonomy_source")) {
    if (column %in% names(previous_meta)) added_metadata[[column]][reused_rows] <- previous_meta[[column]]
  }
}
if (any(!added_metadata$is_bee %in% TRUE)) stop("Expanded BeePrime records must all be bees.")

metadata <- rbind(
  baseline_metadata[, names(added_metadata), drop = FALSE],
  added_metadata[, names(baseline_metadata), drop = FALSE]
)
metadata_match <- match(scores$accession, metadata$accession)
if (anyNA(metadata_match)) stop("BeePrime sequence scores did not map to taxonomy metadata.")
scores <- cbind(scores, metadata[metadata_match, setdiff(names(metadata), "accession"), drop = FALSE])
names(scores)[names(scores) == "forward_penalty"] <- "forward_penalty"
names(scores)[names(scores) == "reverse_penalty"] <- "reverse_penalty"
names(scores)[names(scores) == "forward_scorable"] <- "forward_scorable"
names(scores)[names(scores) == "reverse_scorable"] <- "reverse_scorable"
scores$pair_id <- "BEEPRIME"
scores$pair_label <- "BeePrime"
scores$pair_scorable <- scores$forward_scorable %in% TRUE & scores$reverse_scorable %in% TRUE
scores$pair_penalty <- ifelse(scores$pair_scorable,
                              scores$forward_penalty + scores$reverse_penalty, NA_real_)
scores$perfect_match_pair <- scores$pair_scorable &
  scores$forward_mismatch_count == 0L & scores$reverse_mismatch_count == 0L
scores$terminal_3_mismatch_pair <- scores$pair_scorable &
  (scores$forward_terminal_3_mismatches + scores$reverse_terminal_3_mismatches) > 0L
scores$lineage_group <- "Bee"
scores$reference_mode <- "expanded"
scores$reference_version <- paste0(Sys.getenv("RELEASE_ID", format(Sys.Date(), "%Y-%m-%d")), "-BeePrime")
scores$source_records_represented <- scores$cluster_size
scores <- scores[order(match(scores$accession, c(publication_accessions, added_accessions))), , drop = FALSE]
score_connection <- gzfile(scores_output, open = "wt")
write.csv(scores, score_connection, row.names = FALSE, na = "")
close(score_connection)

geography_source <- "data/derived/reference_geography.csv"
geography_columns <- c("accession", "country_or_territory", "locality", "normalized_location",
                       "geographic_resolution", "latitude", "longitude")
geography <- data.frame(
  accession = scores$accession, country_or_territory = NA_character_, locality = NA_character_,
  normalized_location = NA_character_, geographic_resolution = "unresolved",
  latitude = NA_real_, longitude = NA_real_, stringsAsFactors = FALSE
)
if (file.exists(geography_source)) {
  supplied <- read.csv(geography_source, stringsAsFactors = FALSE, check.names = FALSE)
  if (all(geography_columns %in% names(supplied))) {
    supplied <- supplied[!duplicated(supplied$accession), geography_columns, drop = FALSE]
    m <- match(geography$accession, supplied$accession)
    found <- !is.na(m)
    geography[found, geography_columns[-1]] <- supplied[m[found], geography_columns[-1]]
    geography$geographic_resolution[is.na(geography$geographic_resolution) |
                                      !nzchar(geography$geographic_resolution)] <- "unresolved"
  }
}
write.csv(geography, geography_output, row.names = FALSE, na = "")
message("Built BeePrime reference with ", nrow(publication_bees),
        " publication centroids and ", length(added_accessions), " added bee sequences.")
