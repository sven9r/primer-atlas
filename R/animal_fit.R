# Taxonomic summaries of the actual session's marker-matched sequence searches.
# A denominator is a set of selected input templates, never a whole taxon or a
# PCR success probability. Source geometry and lineage provenance stay separate.
.animal_fit_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
.animal_fit_text <- function(x, default = "") {
  if (is.null(x)) return(default)
  x <- as.character(x); x[is.na(x)] <- default; x
}
.animal_fit_or <- function(x, fallback) if (is.null(x) || !length(x)) fallback else x
.animal_fit_unresolved <- function(x) {
  value <- tolower(trimws(.animal_fit_text(x)))
  !nzchar(value) | value %in% c("na", "n/a", "none", "unknown", "unresolved", "unclassified",
    "uncultured", "incertae sedis", "not assigned", "not available") |
    grepl("^(unresolved|unknown|unclassified|uncultured)([[:space:]]|$)", value)
}
.animal_fit_tuple <- function(...) {
  columns <- list(...)
  encoded <- lapply(columns, function(x) vapply(as.character(x), utils::URLencode,
    character(1), reserved = TRUE))
  do.call(paste, c(encoded, sep = "|"))
}
.animal_fit_versions <- function(x) {
  x <- sort(unique(.animal_fit_text(x))); paste(x[nzchar(x)], collapse = "; ")
}

animal_fit_parse_taxonomy_csv <- function(text, template_ids) {
  if (length(text) != 1L || is.na(text) || !nzchar(trimws(text)))
    stop("Paste CSV containing template_id and at least one taxonomic rank.")
  if (nchar(text, type = "bytes") > 200000L) stop("Keep session taxonomy CSV below 200 kB.")
  rows <- tryCatch(utils::read.csv(text = text, stringsAsFactors = FALSE,
    check.names = FALSE, colClasses = "character"), error = function(e) stop("Could not read taxonomy CSV: ", conditionMessage(e)))
  if (!"template_id" %in% names(rows) || !any(.animal_fit_ranks %in% names(rows)))
    stop("CSV needs template_id and at least one of kingdom, phylum, class, order, family, genus or species.")
  if (anyDuplicated(names(rows))) stop("CSV column names must be unique.")
  unexpected <- setdiff(names(rows), c("template_id", .animal_fit_ranks))
  if (length(unexpected)) stop("Unexpected taxonomy CSV columns: ", paste(unexpected, collapse = ", "), ".")
  if (!nrow(rows) || nrow(rows) > length(template_ids)) stop("Supply one or more rows for existing templates only.")
  if (anyNA(rows$template_id) || any(!nzchar(rows$template_id)) || anyDuplicated(rows$template_id))
    stop("Taxonomy template_id values must be nonempty and unique.")
  if (any(!rows$template_id %in% template_ids)) stop("Every taxonomy template_id must exactly match a currently loaded template.")
  for (rank in .animal_fit_ranks) {
    if (!rank %in% names(rows)) rows[[rank]] <- rep("Unresolved", nrow(rows))
    value <- trimws(.animal_fit_text(rows[[rank]]))
    value[.animal_fit_unresolved(value)] <- "Unresolved"
    if (any(nchar(value) > 200L) || any(grepl("[\r\n]", value))) stop("Use short, single-line taxon labels.")
    rows[[rank]] <- value
  }
  rows$taxonomy_status <- "user_provided_metadata"
  rows$taxonomy_source_url <- ""
  rows$taxonomy_version <- "session_metadata"
  rows$taxonomy_source_sha256 <- digest::digest(text, algo = "sha256", serialize = FALSE)
  rows
}

.animal_fit_taxonomy <- function(templates, taxonomy = NULL, user_taxonomy = NULL) {
  rows <- as.data.frame(templates, stringsAsFactors = FALSE)
  for (field in c(.animal_fit_ranks, "taxonomy_status", "taxonomy_source_url",
      "taxonomy_version", "taxonomy_source_sha256")) {
    if (!field %in% names(rows)) rows[[field]] <- rep("", nrow(rows))
    rows[[field]] <- .animal_fit_text(rows[[field]])
  }
  if (!is.null(taxonomy) && nrow(taxonomy)) {
    if (!is.data.frame(taxonomy) || !"template_id" %in% names(taxonomy) || anyNA(taxonomy$template_id) ||
        anyDuplicated(taxonomy$template_id)) stop("Supplemental taxonomy needs unique exact template identities.")
    index <- match(rows$template_id, taxonomy$template_id)
    take <- !is.na(index)
    # Session headers can intentionally match an accession string. They never
    # inherit a frozen record's lineage from that display identity alone.
    if ("version" %in% names(rows)) take <- take & !grepl("session|user", .animal_fit_text(rows$version), ignore.case = TRUE)
    for (field in c("marker_id", "accession_version")) if (field %in% names(rows) && field %in% names(taxonomy))
      take <- take & !is.na(rows[[field]]) & !is.na(taxonomy[[field]][index]) & rows[[field]] == taxonomy[[field]][index]
    take[is.na(take)] <- FALSE
    for (field in intersect(c(.animal_fit_ranks, "taxid", "taxonomy_status", "taxonomy_source_url",
        "taxonomy_version", "taxonomy_source_sha256"), names(taxonomy))) {
      if (!field %in% names(rows)) rows[[field]] <- rep("", nrow(rows))
      missing <- .animal_fit_unresolved(rows[[field]])
      own <- take & missing
      rows[[field]][own] <- .animal_fit_text(taxonomy[[field]][index[own]])
    }
  }
  if (!is.null(user_taxonomy) && nrow(user_taxonomy)) {
    if (!is.data.frame(user_taxonomy) || !all(c("template_id", .animal_fit_ranks) %in% names(user_taxonomy)) ||
        anyDuplicated(user_taxonomy$template_id) || any(!user_taxonomy$template_id %in% rows$template_id))
      stop("User taxonomy must be validated against the loaded templates.")
    index <- match(user_taxonomy$template_id, rows$template_id)
    for (field in c(.animal_fit_ranks, "taxonomy_status", "taxonomy_source_url", "taxonomy_version", "taxonomy_source_sha256"))
      rows[[field]][index] <- .animal_fit_text(user_taxonomy[[field]])
    if ("taxid" %in% names(rows)) rows$taxid[index] <- ""
    for (field in grep("^raw_taxonomy_", names(rows), value = TRUE)) rows[[field]][index] <- ""
  }
  for (rank in .animal_fit_ranks) {
    missing <- .animal_fit_unresolved(rows[[rank]])
    rows[[rank]][missing] <- "Unresolved"
  }
  rows$taxonomy_status[!nzchar(rows$taxonomy_status)] <- "unresolved_or_unsourced_metadata"
  rows
}

animal_fit_data <- function(result, rank = "class", taxonomy = NULL, user_taxonomy = NULL) {
  if (length(rank) != 1L || !rank %in% .animal_fit_ranks) stop("Choose an explicit taxonomic rank.")
  if (!is.list(result) || !all(c("pairs", "templates", "products") %in% names(result)))
    stop("Taxonomic fit requires the actual selected pair, template and product results.")
  pairs <- as.data.frame(result$pairs, stringsAsFactors = FALSE)
  templates <- as.data.frame(result$templates, stringsAsFactors = FALSE)
  products <- as.data.frame(result$products, stringsAsFactors = FALSE)
  components <- if (is.data.frame(result$components)) result$components else data.frame()
  if (nrow(components) && (!"primer_id" %in% names(components) || anyDuplicated(components$primer_id)))
    stop("Primer provenance needs unique component identities.")
  if (!all(c("pair_id", "pair_label") %in% names(pairs)) || anyNA(pairs$pair_id) ||
      any(!nzchar(pairs$pair_id)) || anyDuplicated(pairs$pair_id))
    stop("Fit requires unique selected pair identities.")
  if (!all(c("template_id", "marker_id") %in% names(templates)) || anyNA(templates$template_id) || any(!nzchar(templates$template_id)) ||
      anyDuplicated(templates$template_id)) stop("Fit requires unique explicit template identities and marker IDs.")
  templates <- .animal_fit_taxonomy(templates, taxonomy, user_taxonomy)
  parents <- .animal_fit_ranks[seq_len(match(rank, .animal_fit_ranks))]
  labels <- lapply(parents, function(field) {
    x <- templates[[field]]; x[.animal_fit_unresolved(x)] <- paste("Unresolved", field); x
  }); names(labels) <- parents
  templates$taxon <- labels[[rank]]
  templates$taxon_resolved <- !.animal_fit_unresolved(templates[[rank]])
  templates$taxon_key <- if (nrow(templates)) paste0(rank, "|", do.call(.animal_fit_tuple, labels)) else character()
  templates$parent_path <- if (length(parents) > 1L)
    do.call(paste, c(labels[-length(labels)], sep = " › ")) else rep("Root lineage", nrow(templates))
  templates$rank <- rep(rank, nrow(templates))
  meta <- list(rank = rank, n_templates = nrow(templates), n_pairs = nrow(pairs),
    marker_id = result$marker_id, settings = result$settings, method = result$method,
    sequence_versions = .animal_fit_versions(templates$version),
    taxonomy_versions = .animal_fit_versions(templates$taxonomy_version),
    denominator = "All selected input templates in the full hierarchy-scoped taxon, including unresolved or unavailable evidence",
    scorable_denominator = "Templates with two unique eligible full-length sites; geometry-withheld sites remain explicit",
    missing_definition = "No eligible sites, withheld product geometry, unavailable/unsupported results or marker mismatch; inspect raw state for the distinction",
    note = "Sequence compatibility on explicit input records is not PCR success, whole-taxon coverage or off-panel specificity.")
  if (!nrow(pairs) || !nrow(templates)) return(list(summary = data.frame(), details = data.frame(), templates = templates, components = components, metadata = meta))
  if (nrow(products) && !all(c("pair_id", "template_id", "status") %in% names(products)))
    stop("Product evidence needs pair_id, template_id and raw status.")
  if (nrow(products) && any(!products$pair_id %in% pairs$pair_id | !products$template_id %in% templates$template_id))
    stop("Product evidence contains an identity outside the selected cohort.")
  grid <- expand.grid(template_index = seq_len(nrow(templates)), pair_index = seq_len(nrow(pairs)),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  details <- templates[grid$template_index, , drop = FALSE]
  if ("status" %in% names(details)) {
    details$template_record_status <- details$status
    details$status <- NULL
  }
  # Sequences are retained by the mixer/FASTA export. This evidence CSV keeps
  # their hashes, exact source identities and positions without repeating them.
  if ("sequence" %in% names(details)) {
    details$template_sequence_sha256 <- vapply(.animal_fit_text(details$sequence),
      digest::digest, character(1), algo = "sha256", serialize = FALSE)
    details$sequence <- NULL
  }
  details$pair_id <- pairs$pair_id[grid$pair_index]
  details$pair_label <- pairs$pair_label[grid$pair_index]
  for (field in setdiff(names(pairs), c("pair_id", "pair_label", names(details))))
    details[[field]] <- pairs[[field]][grid$pair_index]
  for (role in c("forward", "reverse")) {
    id <- paste0(role, "_id")
    if (!id %in% names(details) || !nrow(components)) next
    ci <- match(details[[id]], components$primer_id)
    for (field in intersect(c("citation", "source_url", "source_key", "source_sequence", "catalog_origin",
        "parent_primer_id", "parent_primer_name", "parent_sequence", "normalization_note", "evidence_note"), names(components)))
      details[[paste0(role, "_", field)]] <- .animal_fit_text(components[[field]][ci])
  }
  key <- .animal_fit_tuple(details$pair_id, details$template_id)
  evidence_key <- if (nrow(products)) .animal_fit_tuple(products$pair_id, products$template_id) else character()
  index <- match(key, evidence_key)
  duplicated_key <- evidence_key[duplicated(evidence_key) | duplicated(evidence_key, fromLast = TRUE)]
  for (field in setdiff(names(products), c("pair_id", "template_id", "marker_id", names(details))))
    details[[field]] <- products[[field]][index]
  if (!"status" %in% names(details)) details$status <- rep(NA_character_, nrow(details))
  details$evidence_available <- !is.na(index) & !key %in% duplicated_key
  details$status[is.na(index)] <- "evidence_unavailable"
  details$status[key %in% duplicated_key] <- "inconsistent_duplicate_evidence"
  details$status[is.na(details$status) | !nzchar(details$status)] <- "unsupported_evidence_state"
  details$raw_status <- details$status
  known_states <- c("exact_product", "mismatch_product", "one_eligible_site", "template_too_short",
    "ambiguous_sites", "no_eligible_sites", "non_inward_or_overlapping_sites", "marker_mismatch",
    "template_unavailable", "evidence_unavailable", "inconsistent_duplicate_evidence", "unsupported_evidence_state")
  details$status[!details$status %in% known_states] <- "unsupported_evidence_state"
  expected_marker <- .animal_fit_text(result$marker_id)
  matching_marker <- if (length(expected_marker) == 1L) details$marker_id == expected_marker else rep(FALSE, nrow(details))
  matching_marker[is.na(matching_marker)] <- FALSE
  details$status[!matching_marker] <- "marker_mismatch"
  details$fit_category <- "missing"
  details$fit_category[details$status == "exact_product"] <- "exact"
  details$fit_category[details$status == "mismatch_product"] <- "substitution"
  details$fit_category[details$status %in% c("one_eligible_site", "template_too_short")] <- "partial"
  details$fit_category[details$status == "ambiguous_sites"] <- "ambiguous"
  for (field in c("forward_status", "reverse_status")) {
    if (!field %in% names(details)) details[[field]] <- rep("", nrow(details))
    details[[field]] <- .animal_fit_text(details[[field]])
  }
  details$two_unique_sites <- details$evidence_available & matching_marker &
    details$forward_status %in% c("exact_site", "mismatch_site") &
    details$reverse_status %in% c("exact_site", "mismatch_site")
  details$evidence_available <- details$evidence_available & matching_marker &
    !details$status %in% c("template_unavailable", "unsupported_evidence_state")
  details$length_eligible <- if ("passes_length_filter" %in% names(details))
    !is.na(details$passes_length_filter) & details$passes_length_filter &
      details$fit_category %in% c("exact", "substitution") else rep(FALSE, nrow(details))
  details$allowed_substitutions_per_primer <- rep(.animal_fit_or(result$settings$max_mismatches, NA_integer_), nrow(details))
  details$minimum_product_bp <- rep(.animal_fit_or(result$settings$min_bp, NA_real_), nrow(details))
  details$maximum_product_bp <- rep(.animal_fit_or(result$settings$max_bp, NA_real_), nrow(details))
  rownames(details) <- NULL
  summary <- do.call(rbind, lapply(split(seq_len(nrow(details)), .animal_fit_tuple(details$pair_id, details$taxon_key)), function(ii) {
    rows <- details[ii, , drop = FALSE]
    row <- data.frame(pair_id = rows$pair_id[1L], pair_label = rows$pair_label[1L], rank = rank,
      taxon_key = rows$taxon_key[1L], taxon = rows$taxon[1L], parent_path = rows$parent_path[1L],
      taxon_resolved = rows$taxon_resolved[1L], n_panel = nrow(rows),
      n_exact = sum(rows$fit_category == "exact"), n_substitution = sum(rows$fit_category == "substitution"),
      n_partial = sum(rows$fit_category == "partial"), n_ambiguous = sum(rows$fit_category == "ambiguous"),
      n_missing = sum(rows$fit_category == "missing"), n_scorable = sum(rows$two_unique_sites),
      n_evidence_available = sum(rows$evidence_available), n_length_eligible = sum(rows$length_eligible),
      n_geometry_withheld = sum(rows$status == "non_inward_or_overlapping_sites"),
      n_no_eligible_sites = sum(rows$status == "no_eligible_sites"),
      n_unavailable = sum(!rows$evidence_available),
      n_user_taxonomy = sum(rows$taxonomy_status == "user_provided_metadata"),
      marker_id = if (length(expected_marker) == 1L) expected_marker else "Unspecified",
      sequence_versions = .animal_fit_versions(rows$version),
      taxonomy_versions = .animal_fit_versions(rows$taxonomy_version),
      stringsAsFactors = FALSE)
    row$exact_panel_fraction <- row$n_exact / row$n_panel
    row$exact_scorable_fraction <- if (row$n_scorable) row$n_exact / row$n_scorable else NA_real_
    row
  }))
  summary <- summary[order(!summary$taxon_resolved, summary$taxon, summary$parent_path, summary$pair_label), , drop = FALSE]
  rownames(summary) <- NULL
  list(summary = summary, details = details, templates = templates, components = components, metadata = meta)
}
