#!/usr/bin/env Rscript

# Build a ZBJ-specific site-complete validation panel. This deliberately does not
# reuse the Gurten et al. 67k centroid alignment because most of that alignment
# lacks the early forward-primer site. The four priority clades use annotated
# complete COX1 coding features extracted from complete mitochondrial records,
# aligned onto the same 1,536-bp Drosophila reference coordinate system.
# Raw sequences, accession taxonomy, alignments, and denominators are retained.

project_root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(project_root, ".Rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(jsonlite)
  library(PrimerMiner)
  library(xml2)
})
source(file.path(project_root, "R", "functions.R"))

external_dir <- file.path(project_root, "data", "external", "targeted_zbj")
derived_dir <- file.path(project_root, "data", "derived")
dir.create(external_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)

reference_path <- file.path(
  project_root, "data", "reference", "coi_reference_NC_001322.1.fasta"
)
primer_path <- file.path(project_root, "data", "primers.csv")
binding_path <- file.path(derived_dir, "primer_bindings.csv")
raw_paths <- c(
  Hymenoptera = file.path(
    external_dir, "raw", "Hymenoptera_complete_COX1.fasta"
  ),
  Lepidoptera = file.path(
    external_dir, "raw", "Lepidoptera_complete_COX1.fasta"
  ),
  Acari = file.path(
    external_dir, "raw", "Acari_complete_COX1.fasta"
  ),
  Collembola = file.path(
    external_dir, "raw", "Collembola_complete_COX1.fasta"
  ),
  Coleoptera = file.path(
    external_dir, "raw", "Coleoptera_complete_COX1.fasta"
  )
)
required <- c(reference_path, primer_path, binding_path, raw_paths)
missing <- required[!file.exists(required)]
if (length(missing)) {
  stop("Missing targeted-reference input(s): ", paste(missing, collapse = ", "))
}

mafft <- Sys.which("mafft")
if (!nzchar(mafft)) stop("MAFFT was not found on PATH.")
force <- identical(tolower(Sys.getenv("FORCE", "false")), "true")
ncbi_email <- Sys.getenv("NCBI_EMAIL", "")
ncbi_key <- Sys.getenv("NCBI_API_KEY", "")
query_suffix <- paste0(
  "&tool=coi_primer_atlas",
  if (nzchar(ncbi_email)) paste0("&email=", URLencode(ncbi_email)) else "",
  if (nzchar(ncbi_key)) paste0("&api_key=", URLencode(ncbi_key)) else ""
)

normalize_accession <- function(x) {
  x <- sub("^_R_", "", x)
  sub(" .*", "", x)
}

records <- do.call(rbind, lapply(names(raw_paths), function(group_name) {
  sequences <- parse_fasta(raw_paths[[group_name]])
  data.frame(
    accession = vapply(names(sequences), normalize_accession, character(1)),
    source_group = group_name,
    fasta_header = names(sequences),
    sequence_length = nchar(sequences),
    stringsAsFactors = FALSE
  )
}))
records <- records[!duplicated(records$accession), , drop = FALSE]

batch_values <- function(values, size = 75L) {
  split(values, ceiling(seq_along(values) / size))
}

download_checked <- function(url, destination, validator) {
  if (file.exists(destination) && !force) {
    validator(destination)
    return(invisible(destination))
  }
  temporary <- tempfile(pattern = "targeted_full_coi_")
  on.exit(unlink(temporary), add = TRUE)
  download.file(url, temporary, mode = "wb", quiet = TRUE)
  validator(temporary)
  if (!file.copy(temporary, destination, overwrite = TRUE)) {
    stop("Could not write ", destination)
  }
  Sys.sleep(if (nzchar(ncbi_key)) 0.12 else 0.4)
  invisible(destination)
}

accession_batches <- batch_values(records$accession)
summary_paths <- character(length(accession_batches))
for (i in seq_along(accession_batches)) {
  destination <- file.path(
    external_dir,
    sprintf("nuccore_esummary_batch_%02d.json", i)
  )
  url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
    "?db=nuccore&retmode=json&id=",
    URLencode(paste(accession_batches[[i]], collapse = ","), reserved = TRUE),
    query_suffix
  )
  validator <- function(path) {
    payload <- fromJSON(path, simplifyVector = FALSE)
    if (is.null(payload$result$uids) || !length(payload$result$uids)) {
      stop("NCBI ESummary returned no records for ", basename(destination))
    }
  }
  download_checked(url, destination, validator)
  summary_paths[i] <- destination
}

accession_taxids <- list()
summary_index <- 1L
for (path in summary_paths) {
  payload <- fromJSON(path, simplifyVector = FALSE)
  for (uid in payload$result$uids) {
    item <- payload$result[[uid]]
    accession_taxids[[summary_index]] <- data.frame(
      accession = item$accessionversion,
      taxid = as.character(item$taxid),
      ncbi_title = item$title,
      stringsAsFactors = FALSE
    )
    summary_index <- summary_index + 1L
  }
}
accession_taxids <- do.call(rbind, accession_taxids)
accession_taxids <- accession_taxids[!duplicated(accession_taxids$accession), ]
taxids <- unique(accession_taxids$taxid[nzchar(accession_taxids$taxid)])

taxonomy_batches <- batch_values(taxids)
taxonomy_paths <- character(length(taxonomy_batches))
for (i in seq_along(taxonomy_batches)) {
  destination <- file.path(
    external_dir,
    sprintf("taxonomy_batch_%02d.xml", i)
  )
  url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
    "?db=taxonomy&retmode=xml&id=",
    URLencode(paste(taxonomy_batches[[i]], collapse = ","), reserved = TRUE),
    query_suffix
  )
  validator <- function(path) {
    document <- read_xml(path)
    if (!length(xml_find_all(document, "/TaxaSet/Taxon"))) {
      stop("NCBI Taxonomy returned no records for ", basename(destination))
    }
  }
  download_checked(url, destination, validator)
  taxonomy_paths[i] <- destination
}

rank_value <- function(taxon_node, rank_name) {
  current_rank <- xml_text(xml_find_first(taxon_node, "./Rank"))
  if (identical(current_rank, rank_name)) {
    return(xml_text(xml_find_first(taxon_node, "./ScientificName")))
  }
  lineage_nodes <- xml_find_all(taxon_node, "./LineageEx/Taxon")
  lineage_ranks <- xml_text(xml_find_all(taxon_node, "./LineageEx/Taxon/Rank"))
  index <- which(lineage_ranks == rank_name)
  if (!length(index)) return(NA_character_)
  xml_text(xml_find_first(lineage_nodes[index[length(index)]], "./ScientificName"))
}

taxonomy_rows <- list()
taxonomy_index <- 1L
for (path in taxonomy_paths) {
  document <- read_xml(path)
  for (node in xml_find_all(document, "/TaxaSet/Taxon")) {
    taxonomy_rows[[taxonomy_index]] <- data.frame(
      taxid = xml_text(xml_find_first(node, "./TaxId")),
      scientific_name = xml_text(xml_find_first(node, "./ScientificName")),
      phylum = rank_value(node, "phylum"),
      class = rank_value(node, "class"),
      order = rank_value(node, "order"),
      family = rank_value(node, "family"),
      subfamily = rank_value(node, "subfamily"),
      genus = rank_value(node, "genus"),
      species = rank_value(node, "species"),
      stringsAsFactors = FALSE
    )
    taxonomy_index <- taxonomy_index + 1L
  }
}
taxonomy <- do.call(rbind, taxonomy_rows)
taxonomy <- taxonomy[!duplicated(taxonomy$taxid), , drop = FALSE]

records <- merge(records, accession_taxids, by = "accession", all.x = TRUE)
records <- merge(records, taxonomy, by = "taxid", all.x = TRUE)
records <- records[order(records$source_group, records$accession), ]
write.csv(
  records,
  file.path(derived_dir, "targeted_zbj_taxonomy.csv"),
  row.names = FALSE,
  na = ""
)

alignment_paths <- character(length(raw_paths))
names(alignment_paths) <- names(raw_paths)
for (group_name in names(raw_paths)) {
  output <- file.path(
    external_dir,
    paste0(group_name, "_raw_reference_aligned.fasta")
  )
  log_path <- file.path(
    external_dir,
    paste0(group_name, "_raw_reference_aligned.mafft.log")
  )
  if (!file.exists(output) || force) {
    status <- system2(
      mafft,
      c(
        "--quiet", "--adjustdirectionaccurately", "--keeplength",
        "--addfragments", shQuote(raw_paths[[group_name]]),
        shQuote(reference_path)
      ),
      stdout = output,
      stderr = log_path
    )
    if (!identical(status, 0L)) stop("MAFFT failed; see ", log_path)
  }
  aligned <- parse_fasta(output)
  if (!length(aligned) || any(nchar(aligned) != 1536L)) {
    stop("Invalid 1,536-column alignment: ", output)
  }
  alignment_paths[group_name] <- output
}

primers <- read.csv(primer_path, stringsAsFactors = FALSE, check.names = FALSE)
primers <- primers[
  primers$pair_id %in% c("ZBJ_ART", "ZBJ_ART_DEG"),
  ,
  drop = FALSE
]
bindings <- read.csv(binding_path, stringsAsFactors = FALSE, check.names = FALSE)
bindings <- bindings[
  bindings$pair_id %in% unique(primers$pair_id),
  ,
  drop = FALSE
]

evaluate_primer <- function(alignment_path, primer, binding) {
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
  distances <- as.integer(sub("^V", "", score_columns))
  terminal_count <- function(max_distance) {
    columns <- which(distances <= max_distance)
    if (!length(columns)) return(rep(0L, nrow(evaluated)))
    rowSums(mismatch_matrix[, columns, drop = FALSE])
  }
  data.frame(
    accession = vapply(evaluated$Template, normalize_accession, character(1)),
    binding_sequence_primer_oriented = evaluated$sequ,
    penalty = evaluated$sum,
    scorable = !is.na(evaluated$sum),
    missing_binding_positions = rowSums(missing_matrix),
    mismatch_count = rowSums(mismatch_matrix),
    terminal_3_mismatches = terminal_count(3L),
    stringsAsFactors = FALSE
  )
}

median_safe <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) median(x) else NA_real_
}
quantile_safe <- function(x, probability) {
  x <- x[is.finite(x)]
  if (length(x)) unname(quantile(x, probability, na.rm = TRUE)) else NA_real_
}
summarize_scores <- function(data, groups) {
  for (column in groups) {
    data[[column]][is.na(data[[column]]) | !nzchar(data[[column]])] <-
      paste0(column, " unresolved")
  }
  split_data <- split(
    data,
    interaction(data[, groups, drop = FALSE], drop = TRUE, lex.order = TRUE)
  )
  do.call(rbind, lapply(split_data, function(x) {
    scorable <- x$pair_scorable %in% TRUE
    cbind(
      x[1, groups, drop = FALSE],
      data.frame(
        n_sequences = nrow(x),
        n_species = length(unique(x$species[!is.na(x$species)])),
        n_pair_scorable = sum(scorable),
        pair_scorable_fraction = mean(scorable),
        forward_median_penalty = median_safe(x$forward_penalty),
        reverse_median_penalty = median_safe(x$reverse_penalty),
        pair_median_penalty = median_safe(x$pair_penalty),
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
    )
  }))
}

score_rows <- list()
score_index <- 1L
for (group_name in names(alignment_paths)) {
  for (pair_id in unique(primers$pair_id)) {
    directional <- list()
    for (direction in c("forward", "reverse")) {
      primer <- primers[
        primers$pair_id == pair_id & primers$direction == direction,
        ,
        drop = FALSE
      ]
      binding <- bindings[
        bindings$pair_id == pair_id & bindings$direction == direction,
        ,
        drop = FALSE
      ]
      directional[[direction]] <- evaluate_primer(
        alignment_paths[[group_name]],
        primer,
        binding
      )
    }
    forward <- directional$forward
    reverse <- directional$reverse
    names(forward)[-1] <- paste0("forward_", names(forward)[-1])
    names(reverse)[-1] <- paste0("reverse_", names(reverse)[-1])
    pair_scores <- merge(forward, reverse, by = "accession", sort = FALSE)
    pair_scores$pair_id <- pair_id
    pair_scores$pair_label <- primers$pair_label[
      match(pair_id, primers$pair_id)
    ]
    pair_scores$source_group <- group_name
    pair_scores$pair_scorable <-
      pair_scores$forward_scorable & pair_scores$reverse_scorable
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
    score_rows[[score_index]] <- pair_scores
    score_index <- score_index + 1L
  }
}
scores <- do.call(rbind, score_rows)
scores <- merge(scores, records, by = c("accession", "source_group"), all.x = TRUE)

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
scores$study_group <- scores$source_group
lepidoptera <- scores$source_group == "Lepidoptera"
family_resolved <- !is.na(scores$family) & nzchar(scores$family)
scores$study_group[
  lepidoptera & scores$family %in% butterfly_families
] <- "Butterflies (Papilionoidea)"
scores$study_group[
  lepidoptera & scores$family %in% macroheterocera_families
] <- "Macroheterocera (macro-moth core)"
scores$study_group[
  lepidoptera & scores$family %in% lepidoptera_boundary_families
] <- "Boundary / convention-sensitive moths"
scores$study_group[
  lepidoptera &
    family_resolved &
    !scores$family %in% c(
      butterfly_families, macroheterocera_families,
      lepidoptera_boundary_families
    )
] <- "Microlepidoptera (operational grade)"
scores$study_group[
  lepidoptera & !family_resolved
] <- "Lepidoptera family unresolved"

raw_connection <- gzfile(
  file.path(derived_dir, "targeted_zbj_sequence_scores.csv.gz"),
  open = "wt"
)
write.csv(scores, raw_connection, row.names = FALSE, na = "")
close(raw_connection)

overall <- summarize_scores(scores, c("pair_id", "pair_label", "source_group"))
family <- summarize_scores(
  scores[!is.na(scores$family) & nzchar(scores$family), ],
  c("pair_id", "pair_label", "source_group", "family")
)
genus <- summarize_scores(
  scores[
    !is.na(scores$family) & nzchar(scores$family) &
      !is.na(scores$genus) & nzchar(scores$genus),
  ],
  c("pair_id", "pair_label", "source_group", "family", "genus")
)
subfamily <- summarize_scores(
  scores[
    !is.na(scores$family) & nzchar(scores$family) &
      !is.na(scores$subfamily) & nzchar(scores$subfamily),
  ],
  c("pair_id", "pair_label", "source_group", "family", "subfamily")
)
study_group <- summarize_scores(
  scores,
  c("pair_id", "pair_label", "source_group", "study_group")
)

write.csv(
  overall,
  file.path(derived_dir, "targeted_zbj_overall_summary.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  family,
  file.path(derived_dir, "targeted_zbj_family_summary.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  genus,
  file.path(derived_dir, "targeted_zbj_genus_summary.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  subfamily,
  file.path(derived_dir, "targeted_zbj_subfamily_summary.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  study_group,
  file.path(derived_dir, "targeted_zbj_study_group_summary.csv"),
  row.names = FALSE,
  na = ""
)

message(
  "Wrote targeted ZBJ site-complete scores for ", nrow(scores),
  " pair-sequence comparisons across ", length(raw_paths), " study groups."
)
