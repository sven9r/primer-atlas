# Assemble requested combinations from a frozen, exhaustive component-site
# search. This is a working view; it never changes a UNITE manifest or snapshot.
.unite_candidates_sequence <- function(value) {
  toupper(gsub("[[:space:]]", "", as.character(value)))
}

.unite_candidates_component_key <- function(rows) {
  paste(as.character(rows$forward_primer), .unite_candidates_sequence(rows$forward_sequence),
    as.character(rows$reverse_primer), .unite_candidates_sequence(rows$reverse_sequence), sep = "\n")
}

.unite_candidates_validate_parent_sites <- function(data,variant,parent) {
  query <- strsplit(if (variant$direction=="reverse") reverse_complement(variant$sequence) else variant$sequence,"",fixed=TRUE)[[1L]]
  for (j in which(parent$status=="exact_compatible_site")) {
    sequence <- data$sequences[[parent$record_id[j]]]
    bounds <- c(parent$start[j],parent$end[j])
    valid <- !is.null(sequence) && all(is.finite(bounds)) && all(bounds==floor(bounds)) && bounds[1]>=1 && bounds[2]>=bounds[1] &&
      bounds[2]<=nchar(sequence) && diff(bounds)+1==length(query)
    if ("n_matches" %in% names(parent)) valid <- valid && isTRUE(parent$n_matches[j]==1)
    if (valid) {
      observed <- strsplit(substr(toupper(sequence),bounds[1],bounds[2]),"",fixed=TRUE)[[1L]]
      valid <- all(vapply(seq_along(query),function(k) observed[k] %in% iupac[[query[k]]],logical(1)))
    }
    if (!valid) stop("Frozen UNITE parent binding coordinates do not match their raw template and derived oligo")
  }
  invisible(TRUE)
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
  frozen_catalog <- if (is.data.frame(data$frozen_component_catalog)) data$frozen_component_catalog else data$catalog
  frozen_sites <- if (is.data.frame(data$frozen_component_sites)) data$frozen_component_sites else data$sites
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
    proof$n_catalog == nrow(frozen_catalog),
    proof$n_pairs == frozen_count,
    proof$n_site_searches == nrow(data$metadata) * nrow(frozen_catalog),
    proof$n_sparse_rows == nrow(frozen_sites),
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
  frozen_keys <- paste(frozen_sites$record_id, frozen_sites$catalog_id, sep = "\n")
  if (anyDuplicated(frozen_catalog$catalog_id) || anyDuplicated(frozen_keys) ||
      any(!frozen_sites$catalog_id %in% frozen_catalog$catalog_id) ||
      any(!frozen_sites$record_id %in% data$metadata$record_id)) stop("Frozen UNITE parent site ownership is inconsistent")
  # Working derived rows do not enlarge the immutable 112-component scan.
  # Check that every original column and site row remains exactly intact.
  original_catalog <- data$catalog[match(frozen_catalog$catalog_id,data$catalog$catalog_id),names(frozen_catalog),drop=FALSE]
  original_sites <- data$sites[match(frozen_keys,site_keys),names(frozen_sites),drop=FALSE]
  for (field in names(frozen_catalog)) if (!identical(original_catalog[[field]],frozen_catalog[[field]])) stop("A working catalog changes frozen UNITE parent values")
  for (field in names(frozen_sites)) if (!identical(original_sites[[field]],frozen_sites[[field]])) stop("Working sites change frozen UNITE parent values")
  derived_ids <- data$on_demand_derived_component_ids
  if (is.null(derived_ids)) derived_ids <- character()
  if (!setequal(setdiff(data$catalog$catalog_id,frozen_catalog$catalog_id),derived_ids)) stop("UNITE working catalog has an unverified added component")
  if (length(derived_ids)) {
    if (!exists("fungal_orientation_validate",mode="function")) stop("Derived orientation verification is unavailable")
    frozen <- data;frozen$catalog <- frozen_catalog
    variants <- data$catalog[match(derived_ids,data$catalog$catalog_id),,drop=FALSE]
    fungal_orientation_validate(frozen,variants,require_source_version=FALSE)
    for (i in seq_len(nrow(variants))) {
      parent <- frozen_sites[frozen_sites$catalog_id==variants$parent_catalog_id[i],,drop=FALSE]
      mirrored <- data$sites[data$sites$catalog_id==variants$catalog_id[i],,drop=FALSE]
      if (nrow(parent)!=nrow(mirrored) || !identical(parent$record_id,mirrored$record_id)) stop("Derived UNITE sites do not preserve parent search membership")
      for (field in setdiff(names(frozen_sites),c("catalog_id","primer_name","direction","parent_catalog_id","site_evidence_origin","evidence_version_id"))) {
        if (!identical(parent[[field]],mirrored[[field]])) stop("Derived UNITE sites differ from frozen parent results")
      }
      for (field in intersect(c("primer_name","direction"),names(mirrored))) {
        if (anyNA(mirrored[[field]]) || any(mirrored[[field]]!=variants[[field]][i])) stop("Derived UNITE site identity is inconsistent")
      }
      .unite_candidates_validate_parent_sites(data,variants[i,,drop=FALSE],parent)
    }
  }
  invisible(TRUE)
}

.unite_candidates_bind <- function(original, additions) {
  if (!nrow(additions)) return(original)
  fields <- union(names(original),names(additions))
  for (field in setdiff(fields,names(original))) original[[field]] <- rep(NA,nrow(original))
  for (field in setdiff(fields,names(additions))) additions[[field]] <- rep(NA,nrow(additions))
  rbind(original[,fields,drop=FALSE],additions[,fields,drop=FALSE])
}

# Only requested, explicitly verified opposite-role RC variants are appended.
# An omitted sparse parent result remains a verified no-hit; no scan is run.
.unite_candidates_orientations <- function(data,candidate_pairs,catalog) {
  if (is.null(catalog)) return(list(data=data,reasons=character()))
  if (!is.data.frame(catalog) || !all(c("catalog_id","primer_name","sequence","direction") %in% names(catalog)) ||
      anyNA(catalog$catalog_id) || anyDuplicated(catalog$catalog_id)) stop("Explicit UNITE orientation catalog identity is incomplete or repeated")
  if (!nrow(catalog)) return(list(data=data,reasons=character()))
  requested <- rep(FALSE,nrow(catalog))
  for (direction in c("forward","reverse")) {
    keys <- paste(candidate_pairs[[paste0(direction,"_primer")]],.unite_candidates_sequence(candidate_pairs[[paste0(direction,"_sequence")]]),direction,sep="\n")
    requested <- requested | paste(catalog$primer_name,.unite_candidates_sequence(catalog$sequence),catalog$direction,sep="\n") %in% keys
    id_field <- paste0(direction,"_catalog_id")
    if (id_field %in% names(candidate_pairs)) requested <- requested | catalog$catalog_id %in% candidate_pairs[[id_field]]
  }
  catalog <- catalog[requested,,drop=FALSE]
  if (!nrow(catalog)) return(list(data=data,reasons=character()))
  frozen_catalog <- if (is.data.frame(data$frozen_component_catalog)) data$frozen_component_catalog else data$catalog
  frozen_sites <- if (is.data.frame(data$frozen_component_sites)) data$frozen_component_sites else data$sites
  working <- data;reasons <- character()
  for (i in seq_len(nrow(catalog))) {
    row <- catalog[i,,drop=FALSE];id <- as.character(row$catalog_id)
    old <- working$catalog[working$catalog$catalog_id==id,,drop=FALSE]
    original <- frozen_catalog[frozen_catalog$catalog_id==id,,drop=FALSE]
    if (nrow(original)) {
      if (!identical(as.character(row$primer_name),as.character(original$primer_name)) ||
          !identical(as.character(row$sequence),as.character(original$sequence)) ||
          !identical(as.character(row$direction),as.character(original$direction))) reasons[id] <- "A supplied orientation changes a frozen UNITE component identity"
      next
    }
    verification <- tryCatch({
      if (!exists("fungal_orientation_validate",mode="function")) stop("Derived orientation verification is unavailable")
      frozen <- data;frozen$catalog <- frozen_catalog
      fungal_orientation_validate(frozen,row,require_source_version=FALSE)
      TRUE
    },error=function(error) conditionMessage(error))
    if (!isTRUE(verification)) {reasons[id] <- paste("The derived oligo cannot reuse this UNITE version's parent evidence:",verification);next}
    if (nrow(old)) {
      if (!identical(as.character(old$sequence),as.character(row$sequence)) ||
          !identical(as.character(old$direction),as.character(row$direction))) reasons[id] <- "A derived catalog identifier conflicts with its working component identity"
      next
    }
    parent <- frozen_sites[frozen_sites$catalog_id==row$parent_catalog_id,,drop=FALSE]
    # Check every reused positive against its actual raw template and the new
    # oriented query, including I/N nucleotide-set normalization.
    .unite_candidates_validate_parent_sites(working,row,parent)
    mirrored <- parent
    mirrored$catalog_id <- rep(id,nrow(parent))
    if ("primer_name" %in% names(mirrored)) mirrored$primer_name <- rep(as.character(row$primer_name),nrow(parent))
    if ("direction" %in% names(mirrored)) mirrored$direction <- rep(as.character(row$direction),nrow(parent))
    mirrored$parent_catalog_id <- rep(as.character(row$parent_catalog_id),nrow(parent))
    mirrored$site_evidence_origin <- rep("derived_reverse_complement_parent_frozen_sites",nrow(parent))
    mirrored$evidence_version_id <- rep(data$manifest$id,nrow(parent))
    working$frozen_component_catalog <- frozen_catalog
    working$frozen_component_sites <- frozen_sites
    working$catalog <- .unite_candidates_bind(working$catalog,row)
    working$sites <- .unite_candidates_bind(working$sites,mirrored)
    working$site_keys <- paste(working$sites$record_id,working$sites$catalog_id,sep="\n")
    working$on_demand_derived_component_ids <- unique(c(working$on_demand_derived_component_ids,id))
  }
  list(data=working,reasons=reasons)
}

.unite_candidates_empty_like <- function(rows) {
  out <- rows[FALSE, , drop = FALSE]
  out$candidate_origin <- character()
  out$reason <- character()
  out
}

# Returns data plus accepted/withheld request rows. A withheld request is not an
# observed zero-product result: no taxonomic support is calculated for it.
unite_reference_with_candidates <- function(data, candidate_pairs, catalog=NULL) {
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
  prepared <- .unite_candidates_orientations(data,candidate_pairs,catalog)
  data <- prepared$data
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
      catalog_field <- paste0(direction, "_catalog_id")
      requested_id <- if (catalog_field %in% names(row)) as.character(row[[catalog_field]]) else NA_character_
      if (!is.na(requested_id) && requested_id %in% names(prepared$reasons)) {
        reason <- unname(prepared$reasons[requested_id]);break
      }
      matches <- which(!is.na(data$catalog$primer_name) & !is.na(data$catalog$sequence) &
        !is.na(data$catalog$direction) & data$catalog$primer_name == name &
        .unite_candidates_sequence(data$catalog$sequence) == sequence & data$catalog$direction == direction)
      if (length(matches) != 1L) {
        reason <- paste0("The ", direction, " oligo has no unique exact name, sequence and direction match in this UNITE version")
        break
      }
      component <- data$catalog[matches, , drop = FALSE]
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
    derived <- any(statuses=="derived_reverse_complement")
    row$evidence_status <- if (nrow(existing) && "evidence_status" %in% names(existing)) as.character(existing$evidence_status) else
      if ("unpublished" %in% statuses) "candidate_unpublished_component" else
      if (all(statuses == "publication_cited")) "candidate_published_components" else "candidate_incomplete_provenance"
    if (derived) row$evidence_status <- "candidate_incomplete_provenance"
    row$has_derived_reverse_complement <- derived
    for (direction in c("forward","reverse")) {
      component <- components[[direction]]
      row[[paste0(direction,"_orientation_origin")]] <- if ("orientation_origin" %in% names(component) && !is.na(component$orientation_origin[1L])) as.character(component$orientation_origin) else "frozen_effective_role"
      row[[paste0(direction,"_parent_catalog_id")]] <- if ("parent_catalog_id" %in% names(component) && !is.na(component$parent_catalog_id[1L])) as.character(component$parent_catalog_id) else ""
    }
    row$candidate_origin <- if (nrow(existing) && !id %in% data$on_demand_candidate_pair_ids) "frozen_pair" else
      if (derived) "on_demand_reverse_complement_combination" else "on_demand_component_combination"
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
  if (is.data.frame(data$on_demand_candidate_pair_origins)) {
    prior <- match(working$on_demand_candidate_pair_origins$pair_id,data$on_demand_candidate_pair_origins$pair_id)
    found <- !is.na(prior)
    working$on_demand_candidate_pair_origins$candidate_origin[found] <- data$on_demand_candidate_pair_origins$candidate_origin[prior[found]]
  }
  if (nrow(available)) {
    current <- match(working$on_demand_candidate_pair_origins$pair_id,available$pair_id)
    found <- !is.na(current)
    working$on_demand_candidate_pair_origins$candidate_origin[found] <- available$candidate_origin[current[found]]
  }
  list(data = working, available_pairs = available, unavailable_pairs = unavailable)
}
