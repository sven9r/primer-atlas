# Taxonomic summaries of frozen fungal exact-site evidence. Fractions describe
# retained reference templates, not PCR probabilities or whole-taxon coverage.

.fungal_fit_unresolved <- function(x) {
  value <- tolower(trimws(as.character(x)))
  is.na(value) | !nzchar(value) |
    value %in% c("unresolved", "unknown", "unclassified", "uncultured", "na", "n/a",
                 "none", "not assigned", "not available", "incertae sedis") |
    grepl("^(unresolved|unknown|unclassified|uncultured)[[:space:]]", value)
}

fungal_fit_available_ranks <- function(snapshot) {
  taxonomy <- snapshot$taxonomy
  available <- c("phylum", "class")
  if (!is.null(taxonomy)) {
    for (rank in c("order", "family")) {
      if (rank %in% names(taxonomy) && any(!.fungal_fit_unresolved(taxonomy[[rank]]))) {
        available <- c(available, rank)
      }
    }
  }
  available
}

.fungal_fit_empty_details <- function() {
  data.frame(version_id = character(), pair_id = character(), pair_label = character(),
    rank = character(), group_id = character(), taxon_key = character(), taxon = character(),
    parent_path = character(), taxon_resolved = logical(), accession = character(),
    organism = character(), phylum = character(), class = character(), order = character(),
    family = character(), taxid = character(), taxonomy_source_sha256 = character(),
    source_url = character(), fit_state = character(), reason = character(),
    status = character(), forward_status = character(), reverse_status = character(),
    forward_start = numeric(), forward_end = numeric(), reverse_start = numeric(), reverse_end = numeric(),
    forward_exact = logical(), reverse_exact = logical(), n_exact_sites = integer(),
    taxonomy_unresolved = logical(), amplicon_bp = numeric(), length_basis = character(),
    n_ambiguous_product_bases = integer(), has_ambiguous_product_bases = logical(),
    illumina_eligible = logical(), pacbio_eligible = logical(), nanopore_eligible = logical(),
    length_eligible = logical(), platform = character(), pair_evidence_status = character(),
    pair_source_citation = character(), pair_source_url = character(), stringsAsFactors = FALSE)
}

.fungal_fit_empty_summary <- function() {
  data.frame(version_id = character(), pair_id = character(), pair_label = character(), rank = character(),
    group_id = character(), taxon_key = character(), taxon = character(), parent_path = character(),
    taxon_resolved = logical(), n_panel = integer(), n_evidence_available = integer(),
    n_exact_products = integer(), n_partial_sites = integer(), n_unresolved_sites = integer(),
    n_no_exact_sites = integer(), n_product_withheld = integer(), n_unavailable = integer(),
    n_taxonomy_unresolved = integer(), n_products_with_ambiguous_bases = integer(),
    n_illumina_eligible = integer(), n_pacbio_eligible = integer(), n_nanopore_eligible = integer(),
    n_length_eligible = integer(), exact_product_fraction = numeric(), eligible_panel_fraction = numeric(),
    eligible_product_fraction = numeric(), amplicon_min_bp = numeric(), amplicon_median_bp = numeric(),
    amplicon_max_bp = numeric(), platform = character(), pair_evidence_status = character(),
    stringsAsFactors = FALSE)
}

.fungal_fit_select <- function(requested, available, label) {
  if (is.null(requested)) return(available)
  requested <- unique(as.character(requested))
  if (anyNA(requested) || any(!nzchar(requested)) || any(!requested %in% available)) {
    stop("Selected ", label, " must belong to this frozen snapshot.")
  }
  requested
}

fungal_fit_summary <- function(snapshot, rank = "class", pair_ids = NULL, accessions = NULL,
                               platform = "ALL", illumina_max = 500, pacbio_max = 2000,
                               nanopore_max = 3500) {
  ranks <- fungal_fit_available_ranks(snapshot)
  if (length(rank) != 1L || is.na(rank) || !rank %in% ranks) {
    stop("Taxonomic rank is not supported by this frozen snapshot.")
  }
  platform <- toupper(platform)
  ceilings <- c(ILLUMINA = illumina_max, PACBIO = pacbio_max, NANOPORE = nanopore_max)
  if (length(platform) != 1L || is.na(platform) || !platform %in% c("ALL", names(ceilings))) {
    stop("Platform must be ALL, Illumina, PacBio, or Nanopore.")
  }
  if (length(ceilings) != 3L || any(!is.finite(ceilings)) || any(ceilings <= 0) || any(diff(ceilings) < 0)) {
    stop("Length ceilings must be positive and ordered Illumina <= PacBio <= Nanopore.")
  }
  sequences <- snapshot$sequences
  if (!is.character(sequences) || is.null(names(sequences)) || anyNA(names(sequences)) ||
      any(!nzchar(names(sequences))) || anyDuplicated(names(sequences)) || anyNA(sequences)) {
    stop("Frozen sequences require unique accession names and nonmissing character sequences.")
  }
  sequences <- toupper(gsub("[-[:space:]]", "", sequences))
  pairs <- snapshot$pairs
  evidence <- snapshot$pair_sequences
  if (is.null(pairs) || !"pair_id" %in% names(pairs) || anyNA(pairs$pair_id) ||
      any(!nzchar(pairs$pair_id)) || anyDuplicated(pairs$pair_id)) {
    stop("Frozen pair identities must be unique.")
  }
  if (is.null(evidence)) evidence <- data.frame(pair_id = character(), accession = character())
  if (!all(c("pair_id", "accession") %in% names(evidence)) || anyNA(evidence$pair_id) ||
      anyNA(evidence$accession) || anyDuplicated(paste(evidence$pair_id, evidence$accession, sep = "\r"))) {
    stop("Frozen pair-template evidence must have unique pair and accession identities.")
  }
  if (any(!evidence$pair_id %in% pairs$pair_id) || any(!evidence$accession %in% names(sequences))) {
    stop("Frozen pair-template evidence contains an identity outside the snapshot.")
  }
  taxonomy <- snapshot$taxonomy
  if (is.null(taxonomy)) taxonomy <- data.frame(accession = character())
  if (!"accession" %in% names(taxonomy) || anyNA(taxonomy$accession) || anyDuplicated(taxonomy$accession)) {
    stop("Frozen taxonomy accession identities must be unique.")
  }
  pair_ids <- .fungal_fit_select(pair_ids, as.character(pairs$pair_id), "pair identities")
  accessions <- .fungal_fit_select(accessions, names(sequences), "accession identities")
  version_id <- snapshot$manifest$id
  if (length(version_id) != 1L || is.na(version_id) || !nzchar(version_id)) {
    stop("A frozen snapshot version identity is required.")
  }
  metadata <- list(version_id = version_id, rank = rank, platform = platform,
    n_snapshot_templates = length(sequences), n_selected_templates = length(accessions),
    n_selected_pairs = length(pair_ids), available_ranks = ranks, ceilings = ceilings,
    length_basis = "encoded_ungapped_template_bases_inclusive_primers",
    fraction_denominators = c(exact_product_fraction = "all selected templates in this taxon",
      eligible_panel_fraction = "all selected templates in this taxon",
      eligible_product_fraction = "templates with a stored exact-compatible inward product"),
    evidence_note = paste("Exact-site evidence describes the retained templates. Missing sites can reflect incomplete marker coverage.",
      "It does not estimate mismatch tolerance, PCR success, off-panel specificity, or coverage of the whole taxon."),
    length_note = paste("Length eligibility is planning eligibility from encoded reference bases, not a sequencing or PCR validation.",
      "Internal ambiguous bases and possible source gaps are flagged; encoded length does not establish an unknown assembly gap's physical length."))
  if (!length(pair_ids) || !length(accessions)) {
    return(list(summary = .fungal_fit_empty_summary(), details = .fungal_fit_empty_details(),
                metadata = metadata, ranks = ranks))
  }
  details <- expand.grid(accession = accessions, pair_id = pair_ids,
                         KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  n <- nrow(details)
  details$version_id <- rep(version_id, n)
  details$rank <- rep(rank, n)
  pair_index <- match(details$pair_id, pairs$pair_id)
  text <- function(data, column, index) {
    if (column %in% names(data)) as.character(data[[column]][index]) else rep(NA_character_, length(index))
  }
  details$pair_label <- text(pairs, "pair_label", pair_index)
  missing_label <- is.na(details$pair_label) | !nzchar(details$pair_label)
  details$pair_label[missing_label] <- details$pair_id[missing_label]
  for (column in c("evidence_status", "source_citation", "source_url")) {
    output <- c(evidence_status = "pair_evidence_status", source_citation = "pair_source_citation",
                source_url = "pair_source_url")[[column]]
    details[[output]] <- text(pairs, column, pair_index)
  }
  tax_index <- match(details$accession, taxonomy$accession)
  for (column in c("organism", "phylum", "class", "order", "family", "taxid", "taxonomy_source_sha256", "source_url")) {
    details[[column]] <- text(taxonomy, column, tax_index)
  }
  hierarchy <- c("phylum", "class", "order", "family")
  hierarchy <- hierarchy[seq_len(match(rank, hierarchy))]
  display <- lapply(hierarchy, function(column) {
    value <- details[[column]]
    missing <- .fungal_fit_unresolved(value)
    value[missing] <- paste("Unresolved", column)
    value
  })
  names(display) <- hierarchy
  details$taxon <- display[[rank]]
  details$taxon_resolved <- !.fungal_fit_unresolved(details[[rank]])
  details$taxonomy_unresolved <- !details$taxon_resolved
  encoded <- lapply(display, function(value) vapply(value, utils::URLencode, character(1), reserved = TRUE))
  details$group_id <- paste0(rank, "|", do.call(paste, c(encoded, sep = "|")))
  details$taxon_key <- details$group_id
  details$parent_path <- if (length(display) == 1L) rep("", n) else {
    do.call(paste, c(display[-length(display)], sep = " › "))
  }
  evidence_index <- match(paste(details$pair_id, details$accession, sep = "\r"),
                          paste(evidence$pair_id, evidence$accession, sep = "\r"))
  for (column in c("status", "forward_status", "reverse_status")) {
    details[[column]] <- text(evidence, column, evidence_index)
  }
  for (column in c("forward_start", "forward_end", "reverse_start", "reverse_end", "amplicon_bp")) {
    details[[column]] <- suppressWarnings(as.numeric(text(evidence, column, evidence_index)))
  }
  template_lengths <- nchar(sequences[details$accession])
  valid_site <- function(start, end) {
    is.finite(start) & is.finite(end) & start == floor(start) & end == floor(end) &
      start >= 1 & end >= start & end <= template_lengths
  }
  forward_site <- valid_site(details$forward_start, details$forward_end)
  reverse_site <- valid_site(details$reverse_start, details$reverse_end)
  details$forward_exact <- !is.na(details$forward_status) & details$forward_status == "exact_compatible_site" & forward_site
  details$reverse_exact <- !is.na(details$reverse_status) & details$reverse_status == "exact_compatible_site" & reverse_site
  details$n_exact_sites <- as.integer(details$forward_exact) + as.integer(details$reverse_exact)
  details$fit_state <- rep("evidence_unavailable", n)
  details$reason <- rep("Stored pair-template evidence is unavailable in this version.", n)
  recognized <- c("scorable", "both_sites_missing", "reverse_site_missing", "forward_site_missing",
    "ambiguous_binding_site", "low_specificity_site", "outward_or_overlapping_sites", "below_minimum_length")
  site_states <- c("exact_compatible_site", "no_exact_binding_site", "ambiguous_binding_site",
                   "low_specificity_site", "template_shorter_than_oligo")
  available <- !is.na(evidence_index) & !is.na(details$status) & details$status %in% recognized &
    !is.na(details$forward_status) & details$forward_status %in% site_states &
    !is.na(details$reverse_status) & details$reverse_status %in% site_states
  invalid_exact <- (!is.na(details$forward_status) & details$forward_status == "exact_compatible_site" & !forward_site) |
    (!is.na(details$reverse_status) & details$reverse_status == "exact_compatible_site" & !reverse_site)
  available <- available & !invalid_exact
  ambiguous <- available & (details$forward_status == "ambiguous_binding_site" | details$reverse_status == "ambiguous_binding_site")
  low_specificity <- available & !ambiguous &
    (details$forward_status == "low_specificity_site" | details$reverse_status == "low_specificity_site")
  details$fit_state[ambiguous | low_specificity] <- "unresolved_sites"
  details$reason[ambiguous] <- "Multiple exact-compatible sites prevent a unique stored product."
  details$reason[low_specificity] <- "A match fails the stored specificity guard; an exact-supported product is withheld."
  clean <- available & !ambiguous & !low_specificity
  no_sites <- clean & details$n_exact_sites == 0L & details$status == "both_sites_missing"
  details$fit_state[no_sites] <- "no_exact_sites"
  details$reason[no_sites] <- "No exact-compatible primer site was found in this retained sequence; marker coverage may be incomplete."
  partial <- clean & details$n_exact_sites == 1L &
    ((details$forward_exact & details$status == "reverse_site_missing") |
       (details$reverse_exact & details$status == "forward_site_missing"))
  details$fit_state[partial] <- "partial_sites"
  details$reason[partial & details$forward_exact] <- "Forward site is exact-compatible; no exact reverse site was found in the retained sequence."
  details$reason[partial & details$reverse_exact] <- "Reverse site is exact-compatible; no exact forward site was found in the retained sequence."
  withheld <- clean & details$n_exact_sites == 2L &
    details$status %in% c("outward_or_overlapping_sites", "below_minimum_length")
  details$fit_state[withheld] <- "product_withheld"
  details$reason[withheld & details$status == "outward_or_overlapping_sites"] <- "Both sites are exact-compatible, but outward or overlapping geometry withholds an inward product."
  details$reason[withheld & details$status == "below_minimum_length"] <- "Both sites are exact-compatible, but the stored product is below the evaluation's minimum length."
  expected_bp <- details$reverse_end - details$forward_start + 1
  exact <- clean & details$n_exact_sites == 2L & details$status == "scorable" &
    details$forward_end < details$reverse_start & is.finite(details$amplicon_bp) &
    details$amplicon_bp > 0 & details$amplicon_bp == expected_bp
  exact[is.na(exact)] <- FALSE
  details$fit_state[exact] <- "exact_product"
  details$reason[exact] <- "Two unique exact-compatible sites support a stored inward product on this retained template."
  inconsistent <- available & details$fit_state == "evidence_unavailable"
  details$reason[inconsistent | invalid_exact] <- "Stored site status, raw coordinates, or product length are incomplete or inconsistent; product support is unavailable."
  unknown <- !is.na(evidence_index) & !available & !invalid_exact
  details$reason[unknown] <- paste0("Stored evidence status is unavailable or unsupported: ",
    ifelse(is.na(details$status[unknown]), "missing status", details$status[unknown]), ".")
  details$amplicon_bp[!exact] <- NA_real_
  details$length_basis <- ifelse(exact, "encoded_ungapped_template_bases_inclusive_primers", "unavailable")
  details$n_ambiguous_product_bases <- rep(NA_integer_, n)
  for (i in which(exact)) {
    product <- substring(sequences[[details$accession[i]]], details$forward_start[i], details$reverse_end[i])
    details$n_ambiguous_product_bases[i] <- nchar(gsub("[ACGT]", "", product))
  }
  details$has_ambiguous_product_bases <- !is.na(details$n_ambiguous_product_bases) & details$n_ambiguous_product_bases > 0L
  details$illumina_eligible <- exact & !is.na(details$amplicon_bp) & details$amplicon_bp <= illumina_max
  details$pacbio_eligible <- exact & !is.na(details$amplicon_bp) & details$amplicon_bp <= pacbio_max
  details$nanopore_eligible <- exact & !is.na(details$amplicon_bp) & details$amplicon_bp <= nanopore_max
  details$length_eligible <- if (platform == "ALL") details$nanopore_eligible else {
    details[[paste0(tolower(platform), "_eligible")]]
  }
  details$platform <- rep(platform, n)
  details <- details[, names(.fungal_fit_empty_details()), drop = FALSE]
  keys <- paste(details$pair_id, details$group_id, sep = "\r")
  grouped <- split(seq_len(n), factor(keys, levels = unique(keys)))
  summary <- do.call(rbind, lapply(grouped, function(indices) {
    rows <- details[indices, , drop = FALSE]
    result <- rows[1, c("version_id", "pair_id", "pair_label", "rank", "group_id", "taxon_key",
                       "taxon", "parent_path", "taxon_resolved", "platform", "pair_evidence_status"), drop = FALSE]
    result$n_panel <- nrow(rows)
    for (state in c(exact_products = "exact_product", partial_sites = "partial_sites", unresolved_sites = "unresolved_sites",
                    no_exact_sites = "no_exact_sites", product_withheld = "product_withheld", unavailable = "evidence_unavailable")) {
      # The state-to-column mapping is explicit so these categories partition n_panel.
      column <- c(exact_product = "n_exact_products", partial_sites = "n_partial_sites", unresolved_sites = "n_unresolved_sites",
                  no_exact_sites = "n_no_exact_sites", product_withheld = "n_product_withheld", evidence_unavailable = "n_unavailable")[[state]]
      result[[column]] <- sum(rows$fit_state == state)
    }
    result$n_evidence_available <- result$n_panel - result$n_unavailable
    result$n_taxonomy_unresolved <- sum(rows$taxonomy_unresolved)
    result$n_products_with_ambiguous_bases <- sum(rows$has_ambiguous_product_bases)
    for (system in c("illumina", "pacbio", "nanopore")) {
      result[[paste0("n_", system, "_eligible")]] <- sum(rows[[paste0(system, "_eligible")]])
    }
    result$n_length_eligible <- sum(rows$length_eligible)
    result$exact_product_fraction <- result$n_exact_products / result$n_panel
    result$eligible_panel_fraction <- result$n_length_eligible / result$n_panel
    result$eligible_product_fraction <- if (result$n_exact_products) result$n_length_eligible / result$n_exact_products else NA_real_
    lengths <- rows$amplicon_bp[rows$fit_state == "exact_product"]
    result$amplicon_min_bp <- if (length(lengths)) min(lengths) else NA_real_
    result$amplicon_median_bp <- if (length(lengths)) median(lengths) else NA_real_
    result$amplicon_max_bp <- if (length(lengths)) max(lengths) else NA_real_
    result
  }))
  summary <- summary[, names(.fungal_fit_empty_summary()), drop = FALSE]
  rownames(summary) <- rownames(details) <- NULL
  list(summary = summary, details = details, metadata = metadata, ranks = ranks)
}
