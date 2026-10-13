# Assemble requested combinations from a frozen, exhaustive component-site
# search. This is a working view; it never changes a UNITE manifest or snapshot.
.unite_candidates_sequence <- function(value) {
  toupper(gsub("[[:space:]]", "", as.character(value)))
}

.unite_candidates_component_key <- function(rows) {
  paste(as.character(rows$forward_primer), .unite_candidates_sequence(rows$forward_sequence),
    as.character(rows$reverse_primer), .unite_candidates_sequence(rows$reverse_sequence), sep = "\n")
}

.unite_candidates_check_search <- function(data) {
  if (!is.list(data) || !is.data.frame(data$metadata) || !is.data.frame(data$catalog) ||
      !is.data.frame(data$sites) || !is.data.frame(data$pairs)) stop("UNITE working data is incomplete")
  if (!all(c("record_id") %in% names(data$metadata)) ||
      !all(c("catalog_id", "primer_name", "sequence", "direction") %in% names(data$catalog)) ||
      !all(c("record_id", "catalog_id", "status", "start", "end") %in% names(data$sites))) {
    stop("UNITE component-site identity fields are incomplete")
  }
  proof <- data$search_manifest
  artifacts <- data$manifest$artifacts
  hashes <- c(proof$sequence_sha256, proof$catalog_artifact_sha256, proof$pair_artifact_sha256,
    proof$site_artifact_sha256, proof$record_ids_sha256)
  if (length(hashes) != 5L || anyNA(hashes) || any(!grepl("^[a-f0-9]{64}$", hashes))) {
    stop("UNITE exhaustive search checksum evidence is incomplete")
  }
  prior_ids <- data$on_demand_candidate_pair_ids
  if (is.null(prior_ids)) prior_ids <- character()
  frozen_count <- nrow(data$pairs) - sum(data$pairs$pair_id %in% prior_ids)
  checks <- c(
    isTRUE(proof$complete),
    identical(proof$scan_rule_version, "exact_iupac_unique_sites_v1"),
    identical(proof$sequence_sha256, artifacts$sequences$sha256),
    identical(proof$catalog_artifact_sha256, artifacts$catalog$sha256),
    identical(proof$pair_artifact_sha256, artifacts$pairs$sha256),
    identical(proof$site_artifact_sha256, artifacts$sites$sha256),
    identical(proof$record_ids_sha256, digest::digest(paste(sort(data$metadata$record_id), collapse = "\n"),
      algo = "sha256", serialize = FALSE)),
    proof$n_records == nrow(data$metadata),
    proof$n_catalog == nrow(data$catalog),
    proof$n_pairs == frozen_count,
    proof$n_site_searches == nrow(data$metadata) * nrow(data$catalog),
    proof$n_sparse_rows == nrow(data$sites),
    proof$n_no_exact_hits + proof$n_sparse_rows == proof$n_site_searches
  )
  if (length(checks) != 13L || anyNA(checks) || !all(checks)) {
    stop("On-demand combinations require the complete frozen all-component UNITE site-search proof")
  }
  site_keys <- paste(data$sites$record_id, data$sites$catalog_id, sep = "\n")
  if (anyNA(data$metadata$record_id) || anyDuplicated(data$metadata$record_id) ||
      anyDuplicated(data$catalog$catalog_id) || anyDuplicated(data$pairs$pair_id) ||
      anyDuplicated(site_keys) || any(!data$sites$record_id %in% data$metadata$record_id) ||
      any(!data$sites$catalog_id %in% data$catalog$catalog_id) ||
      !identical(data$site_keys, site_keys)) stop("UNITE precomputed site ownership is inconsistent")
  invisible(TRUE)
}

.unite_candidates_empty_like <- function(rows) {
  out <- rows[FALSE, , drop = FALSE]
  out$candidate_origin <- character()
  out$reason <- character()
  out
}

# Returns data plus accepted/withheld request rows. A withheld request is not an
# observed zero-product result: no taxonomic support is calculated for it.
unite_reference_with_candidates <- function(data, candidate_pairs) {
  if (is.null(candidate_pairs)) candidate_pairs <- data.frame()
  if (!is.data.frame(candidate_pairs)) stop("Candidate combinations must be a data frame")
  available <- .unite_candidates_empty_like(candidate_pairs)
  unavailable <- .unite_candidates_empty_like(candidate_pairs)
  if (!nrow(candidate_pairs)) return(list(data = data, available_pairs = available, unavailable_pairs = unavailable))
  required <- c("pair_id", "pair_label", "forward_primer", "reverse_primer", "forward_sequence", "reverse_sequence")
  if (!all(required %in% names(candidate_pairs))) stop("Candidate component definitions are missing: ",
    paste(setdiff(required, names(candidate_pairs)), collapse = ", "))
  candidate_pairs <- as.data.frame(candidate_pairs, stringsAsFactors = FALSE)
  if (anyNA(candidate_pairs$pair_id) || any(!nzchar(as.character(candidate_pairs$pair_id)))) {
    stop("Candidate pair identifiers are missing")
  }
  keys <- .unite_candidates_component_key(candidate_pairs)
  duplicate_ids <- unique(candidate_pairs$pair_id[duplicated(candidate_pairs$pair_id)])
  if (any(vapply(duplicate_ids, function(id) length(unique(keys[candidate_pairs$pair_id == id])) > 1L, logical(1)))) {
    stop("A requested pair identifier refers to conflicting component identities")
  }
  candidate_pairs <- candidate_pairs[!duplicated(candidate_pairs$pair_id), , drop = FALSE]
  if (is.null(data)) {
    unavailable <- candidate_pairs
    unavailable$candidate_origin <- rep("unavailable", nrow(unavailable))
    unavailable$reason <- rep("The frozen UNITE reference collection is unavailable", nrow(unavailable))
    return(list(data = NULL, available_pairs = available, unavailable_pairs = unavailable))
  }
  .unite_candidates_check_search(data)
  original_columns <- names(data$pairs)
  if (!all(required %in% original_columns)) stop("Frozen UNITE pair definitions are incomplete")
  working <- data
  if (is.null(working$on_demand_candidate_pair_ids)) working$on_demand_candidate_pair_ids <- character()
  accepted <- list()
  withheld <- list()
  additions <- list()
  for (i in seq_len(nrow(candidate_pairs))) {
    row <- candidate_pairs[i, , drop = FALSE]
    id <- as.character(row$pair_id)
    reason <- ""
    components <- list()
    for (direction in c("forward", "reverse")) {
      name <- as.character(row[[paste0(direction, "_primer")]])
      sequence <- .unite_candidates_sequence(row[[paste0(direction, "_sequence")]])
      matches <- which(!is.na(data$catalog$primer_name) & !is.na(data$catalog$sequence) &
        !is.na(data$catalog$direction) & data$catalog$primer_name == name &
        .unite_candidates_sequence(data$catalog$sequence) == sequence & data$catalog$direction == direction)
      if (length(matches) != 1L) {
        reason <- paste0("The ", direction, " oligo has no unique exact name, sequence and direction match in this UNITE version")
        break
      }
      component <- data$catalog[matches, , drop = FALSE]
      catalog_field <- paste0(direction, "_catalog_id")
      if (catalog_field %in% names(row) && !is.na(row[[catalog_field]]) && nzchar(as.character(row[[catalog_field]])) &&
          !identical(as.character(row[[catalog_field]]), as.character(component$catalog_id))) {
        reason <- paste0("The ", direction, " catalog identifier conflicts with its frozen component identity")
        break
      }
      components[[direction]] <- component
    }
    existing <- data$pairs[data$pairs$pair_id == id, , drop = FALSE]
    if (!nzchar(reason) && nrow(existing) && !identical(
        .unite_candidates_component_key(row), .unite_candidates_component_key(existing))) {
      reason <- "The requested pair identifier conflicts with the frozen UNITE component definition"
    }
    if (nzchar(reason)) {
      row$candidate_origin <- "unavailable"
      row$reason <- reason
      withheld[[length(withheld) + 1L]] <- row
      next
    }
    for (direction in c("forward", "reverse")) {
      row[[paste0(direction, "_sequence")]] <- as.character(components[[direction]]$sequence)
      row[[paste0(direction, "_catalog_id")]] <- as.character(components[[direction]]$catalog_id)
    }
    statuses <- vapply(components, function(component) {
      if (!"primary_reference_status" %in% names(component) || is.na(component$primary_reference_status[1L])) "" else component$primary_reference_status[1L]
    }, character(1))
    row$evidence_status <- if (nrow(existing) && "evidence_status" %in% names(existing)) as.character(existing$evidence_status) else
      if ("unpublished" %in% statuses) "candidate_unpublished_component" else
      if (all(statuses == "publication_cited")) "candidate_published_components" else "candidate_incomplete_provenance"
    row$candidate_origin <- if (nrow(existing) && !id %in% data$on_demand_candidate_pair_ids) "frozen_pair" else "on_demand_component_combination"
    row$reason <- ""
    accepted[[length(accepted) + 1L]] <- row
    if (nrow(existing)) next
    addition <- row
    # NCBI-map statistics supplied by the caller do not describe this cohort.
    addition$n_panel <- nrow(data$metadata)
    for (field in c("n_scorable", "amplicon_min_bp", "amplicon_median_bp", "amplicon_max_bp")) addition[[field]] <- NA_real_
    for (field in setdiff(original_columns, names(addition))) addition[[field]] <- NA
    addition <- addition[, original_columns, drop = FALSE]
    # Match the frozen table's column types, avoiding coercion of old values.
    for (field in original_columns) {
      if (is.integer(data$pairs[[field]])) addition[[field]] <- as.integer(addition[[field]]) else
        if (is.numeric(data$pairs[[field]])) addition[[field]] <- as.numeric(addition[[field]]) else
          if (is.logical(data$pairs[[field]])) addition[[field]] <- as.logical(addition[[field]]) else
            addition[[field]] <- as.character(addition[[field]])
    }
    additions[[length(additions) + 1L]] <- addition
  }
  if (length(additions)) {
    appended <- do.call(rbind, additions)
    working$pairs <- rbind(data$pairs, appended)
    working$on_demand_candidate_pair_ids <- unique(c(data$on_demand_candidate_pair_ids, appended$pair_id))
  }
  if (length(accepted)) available <- do.call(rbind, accepted)
  if (length(withheld)) unavailable <- do.call(rbind, withheld)
  working$on_demand_candidate_pair_origins <- data.frame(pair_id = working$on_demand_candidate_pair_ids,
    candidate_origin = rep("on_demand_component_combination", length(working$on_demand_candidate_pair_ids)), stringsAsFactors = FALSE)
  list(data = working, available_pairs = available, unavailable_pairs = unavailable)
}
