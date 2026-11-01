# Broad ITS evidence has its own immutable identity. It never fills a long-rDNA
# snapshot with trimmed ITS fragments or borrows the old alignment coordinates.
.unite_reference_cache <- new.env(parent = emptyenv())

.unite_read_fasta <- function(path) {
  connection <- if (grepl("[.]gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(connection))
  lines <- readLines(connection, warn = FALSE)
  starts <- which(startsWith(lines, ">"))
  if (!length(starts)) return(setNames(character(), character()))
  ends <- c(starts[-1L]-1L, length(lines))
  ids <- sub("[[:space:]].*$", "", substring(lines[starts], 2L))
  if (anyDuplicated(ids) || any(!nzchar(ids))) stop("UNITE FASTA identities are invalid")
  values <- vapply(seq_along(starts), function(i) {
    if (ends[i] <= starts[i]) "" else paste0(lines[seq.int(starts[i]+1L, ends[i])], collapse = "")
  }, character(1))
  setNames(toupper(values), ids)
}

.unite_safe_path <- function(root, relative) {
  if (length(relative) != 1L || is.na(relative) || !nzchar(relative) ||
      grepl("(^/|(^|/)[.][.](/|$)|\\\\)", relative)) stop("Unsafe UNITE artifact path")
  file.path(root, relative)
}

unite_reference_index <- function(root = getwd()) {
  path <- file.path(root, "data", "pinned", "UNITE_ITS", "index.json")
  if (!file.exists(path)) return(data.frame(id = character(), label = character(), manifest = character(), sha256 = character(), n_selected = integer()))
  entries <- jsonlite::read_json(path, simplifyVector = FALSE)$versions
  if (!length(entries)) return(data.frame())
  rows <- do.call(rbind, lapply(entries, function(entry) data.frame(id = entry$id, label = entry$label,
      manifest = entry$manifest, sha256 = entry$sha256, n_selected = entry$n_selected, stringsAsFactors = FALSE)))
  if (anyDuplicated(rows$id)) stop("UNITE historical identity repeats")
  rows
}

unite_reference_read <- function(root = getwd(), version_id = NULL) {
  directory <- file.path(root, "data", "pinned", "UNITE_ITS")
  pointer <- file.path(directory, "latest.json")
  if (!file.exists(pointer)) return(NULL)
  latest <- jsonlite::read_json(pointer, simplifyVector = FALSE)
  versions <- unite_reference_index(root)
  if (!is.null(version_id) && nzchar(version_id)) {
    row <- versions[versions$id == version_id, , drop = FALSE]
    if (nrow(row) != 1L) stop("Unknown UNITE reference version")
    latest <- as.list(row[1L, , drop = FALSE])
  }
  path <- .unite_safe_path(directory, latest$manifest)
  if (!file.exists(path)) stop("UNITE reference manifest is missing")
  digest <- atlas_sha256(path)
  if (!identical(digest, latest$sha256)) stop("UNITE pointer checksum differs")
  key <- paste(normalizePath(root), latest$id, sep = "|")
  if (exists(key, .unite_reference_cache, inherits = FALSE)) {
    value <- get(key, .unite_reference_cache)
    if (!identical(attr(value, "manifest_sha256"), digest)) stop("Immutable UNITE identity changed")
    value$versions <- versions
    return(value)
  }
  manifest <- jsonlite::read_json(path, simplifyVector = FALSE)
  if (!identical(manifest$id, latest$id) || !identical(manifest$marker_id, "UNITE_ITS")) stop("UNITE manifest identity differs")
  base <- dirname(path)
  paths <- lapply(manifest$artifacts, function(artifact) {
    file <- .unite_safe_path(base, artifact$path)
    if (!file.exists(file) || !identical(atlas_sha256(file), artifact$sha256)) stop("UNITE artifact checksum failed: ", artifact$path)
    file
  })
  value <- list(manifest = manifest, directory = base, paths = paths, versions = versions)
  for (name in c("metadata", "coverage", "regions", "alignment_index", "sites", "catalog", "pairs")) {
    if (is.null(paths[[name]])) stop("Missing UNITE artifact: ", name)
    value[[name]] <- as.data.frame(readr::read_csv(paths[[name]], col_types = readr::cols(.default = readr::col_character()),
        show_col_types = FALSE), stringsAsFactors = FALSE)
  }
  numeric_fields <- c("raw_length_bp", "standard_length_bp", "ambiguity_count", "ambiguity_fraction", "n_source", "n_eligible_unique",
      "n_selected", "shortfall_to_100", "start", "end", "n_matches", "n_panel", "n_scorable", "amplicon_min_bp", "amplicon_median_bp",
      "amplicon_max_bp", "n_columns", "n_input", "n_source_regions")
  for (name in c("metadata", "coverage", "regions", "alignment_index", "sites", "pairs")) {
    for (column in intersect(names(value[[name]]), numeric_fields)) value[[name]][[column]] <- suppressWarnings(as.numeric(value[[name]][[column]]))
  }
  value$coverage$taxon_resolved <- tolower(value$coverage$taxon_resolved) %in% c("true", "1", "t")
  meta <- value$metadata
  if (anyNA(meta$record_id) || anyDuplicated(meta$record_id) || any(!nzchar(meta$record_id))) stop("UNITE metadata identities repeat or are missing")
  if (!"accession" %in% names(meta)) meta$accession <- meta$accession_raw
  meta$organism <- if ("organism" %in% names(meta)) meta$organism else gsub("_", " ", meta$species)
  if (!"orientation" %in% names(meta)) meta$orientation <- "source"
  value$metadata <- meta
  value$sequences <- .unite_read_fasta(paths$sequences)
  if (!setequal(names(value$sequences), meta$record_id) || nrow(meta) != manifest$n_selected) stop("UNITE retained membership differs from manifest")
  if (any(!value$sites$record_id %in% meta$record_id) || any(!value$sites$catalog_id %in% value$catalog$catalog_id) ||
      anyDuplicated(paste(value$sites$record_id, value$sites$catalog_id, sep = "\n"))) stop("UNITE site ownership differs")
  if (anyDuplicated(value$catalog$catalog_id) || anyDuplicated(value$pairs$pair_id)) stop("UNITE catalog identity repeats")
  if (any(!value$regions$record_id %in% meta$record_id)) stop("UNITE region ownership differs")
  if (is.null(paths$search_manifest)) stop("UNITE exact-site search completion proof is missing")
  proof <- jsonlite::read_json(paths$search_manifest, simplifyVector = FALSE)
  proof_checks <- c(isTRUE(proof$complete), identical(proof$scan_rule_version, "exact_iupac_unique_sites_v1"),
      identical(proof$sequence_sha256, manifest$artifacts$sequences$sha256),
      identical(proof$catalog_artifact_sha256, manifest$artifacts$catalog$sha256),
      identical(proof$pair_artifact_sha256, manifest$artifacts$pairs$sha256),
      identical(proof$site_artifact_sha256, manifest$artifacts$sites$sha256),
      identical(proof$record_ids_sha256, digest::digest(paste(sort(meta$record_id), collapse = "\n"), algo = "sha256", serialize = FALSE)),
      proof$n_records == nrow(meta), proof$n_catalog == nrow(value$catalog), proof$n_pairs == nrow(value$pairs),
      proof$n_site_searches == nrow(meta)*nrow(value$catalog), proof$n_sparse_rows == nrow(value$sites),
      proof$n_no_exact_hits + proof$n_sparse_rows == proof$n_site_searches)
  if (length(proof_checks) != 13L || anyNA(proof_checks) || !all(proof_checks)) stop("UNITE exhaustive site-search proof differs from the loaded cohort")
  value$search_manifest <- proof
  value$site_keys <- paste(value$sites$record_id, value$sites$catalog_id, sep = "\n")
  attr(value, "manifest_sha256") <- digest
  assign(key, value, .unite_reference_cache)
  value
}

unite_reference_groups <- function(data, rank = "class", search = "") {
  if (is.null(data)) return(data.frame())
  if (!rank %in% c("phylum", "class", "order", "family")) stop("Unknown UNITE taxonomic rank")
  rows <- data$coverage[data$coverage$rank == rank, , drop = FALSE]
  search <- trimws(search)
  if (nzchar(search)) rows <- rows[grepl(tolower(search), tolower(paste(rows$taxon, rows$parent_path)), fixed = TRUE), , drop = FALSE]
  rows[order(!rows$taxon_resolved, -rows$n_selected, rows$parent_path, rows$taxon), , drop = FALSE]
}

unite_reference_alignment <- function(data, group_id, region = "ITS2", limit = 12L) {
  if (!region %in% c("ITS1", "5.8S", "ITS2")) stop("Unknown UNITE alignment region")
  scope <- "Independently aligned region within one frozen UNITE clade; computational ITSx boundaries."
  empty <- list(alignment = setNames(character(), character()), metadata = data.frame(), missing_metadata = data.frame(),
      n_total = 0L, n_aligned = 0L, n_columns = 0L, reason = "No regional alignment is available for this group.", scope = scope)
  if (is.null(data) || length(group_id) != 1L || is.na(group_id)) return(empty)
  owner <- if (is.data.frame(data$coverage)) data$coverage[data$coverage$group_id == group_id & data$coverage$rank %in% c("class", "order"), , drop = FALSE] else data.frame()
  rank <- if (nrow(owner) == 1L) owner$rank else "class"
  meta <- data$metadata[data$metadata[[paste0(rank, "_group_id")]] == group_id, , drop = FALSE]
  empty$n_total <- nrow(meta)
  index <- data$alignment_index[data$alignment_index$group_id == group_id & data$alignment_index$region == region, , drop = FALSE]
  if (!nrow(index)) {
    empty$missing_metadata <- head(meta, limit)
    return(empty)
  }
  if (nrow(index) != 1L) stop("Regional alignment identity repeats")
  file <- .unite_safe_path(data$directory, index$file[1L])
  key <- paste(normalizePath(data$directory), attr(data, "manifest_sha256") %||% "", data$manifest$id, group_id, region, sep = "|")
  if (exists(key, .unite_reference_cache, inherits = FALSE)) alignment <- get(key, .unite_reference_cache) else {
    alignment <- .unite_read_fasta(file)
    if (length(unique(nchar(alignment))) > 1L || any(!names(alignment) %in% meta$record_id)) stop("UNITE regional alignment has invalid width or ownership")
    regions <- data$regions[data$regions$region == region & data$regions$boundary_status == "complete", , drop = FALSE]
    for (id in names(alignment)) {
      row <- regions[regions$record_id == id, , drop = FALSE]
      if (nrow(row) != 1L) stop("Aligned region lacks unique complete boundaries")
      expected <- substr(data$sequences[[id]], row$start, row$end)
      if (!identical(gsub("-", "", alignment[[id]], fixed = TRUE), expected)) stop("Regional alignment does not retain its source bases")
    }
    assign(key, alignment, .unite_reference_cache)
  }
  shown <- head(alignment, max(1L, as.integer(limit)))
  metadata <- meta[match(names(shown), meta$record_id), , drop = FALSE]
  complete_ids <- data$regions$record_id[data$regions$region == region & data$regions$boundary_status == "complete"]
  missing <- meta[!meta$record_id %in% complete_ids, , drop = FALSE]
  list(alignment = shown, metadata = metadata, missing_metadata = head(missing, limit), n_missing_boundaries = nrow(missing),
      n_total = nrow(meta), n_aligned = length(alignment), n_regions = sum(meta$record_id %in% complete_ids),
      n_columns = if (length(alignment)) nchar(alignment[[1L]]) else 0L, reason = "", scope = scope)
}

.unite_pair_site <- function(data, primer_name, sequence, direction) {
  catalog <- data$catalog
  row <- which(catalog$primer_name == primer_name & catalog$sequence == sequence & catalog$direction == direction)
  if (length(row) != 1L) stop("UNITE pair component identity is unresolved")
  ids <- data$metadata$record_id
  position <- match(paste(ids, catalog$catalog_id[row], sep = "\n"), data$site_keys)
  rows <- data.frame(start = data$sites$start[position], end = data$sites$end[position], status = data$sites$status[position], stringsAsFactors = FALSE)
  rows$status[is.na(position)] <- "no_exact_binding_site"
  rows
}

unite_reference_pair_fit <- function(data, pair_id, rank = "class") {
  if (is.null(data)) return(list(summary = data.frame(), details = data.frame()))
  pair <- data$pairs[data$pairs$pair_id == pair_id, , drop = FALSE]
  if (nrow(pair) != 1L) stop("Choose a pair in the frozen UNITE search catalog")
  f <- .unite_pair_site(data, pair$forward_primer, pair$forward_sequence, "forward")
  r <- .unite_pair_site(data, pair$reverse_primer, pair$reverse_sequence, "reverse")
  fa <- f$status == "exact_compatible_site"; ra <- r$status == "exact_compatible_site"
  state <- rep("no_exact_sites", nrow(f))
  state[xor(fa, ra)] <- "partial_sites"
  state[f$status %in% c("ambiguous_binding_site", "low_specificity_site") | r$status %in% c("ambiguous_binding_site", "low_specificity_site")] <- "unresolved_sites"
  state[f$status %in% c("invalid_oligo", "template_shorter_than_oligo") | r$status %in% c("invalid_oligo", "template_shorter_than_oligo")] <- "evidence_unavailable"
  both <- fa & ra
  state[both] <- "product_withheld"
  length <- r$end - f$start + 1L
  inward <- both & !is.na(f$end) & !is.na(r$start) & f$end < r$start & length >= 100
  state[inward] <- "exact_product"
  length[!inward] <- NA_real_
  reasons <- c(exact_product = "Two unique exact-compatible sites support an inward encoded product; no PCR validation inferred.",
      partial_sites = "One unique exact site found; the other site may be outside the retained fragment or contain mismatches.",
      no_exact_sites = "Neither exact site found in this record; fragment coverage and primer mismatches are not distinguished by this search.",
      unresolved_sites = "Multiple compatible hits or low site specificity prevent a unique product assignment.",
      evidence_unavailable = "The oligo or retained template does not support this site search.", product_withheld = "Exact sites are overlapping, outward, or encode a product shorter than 100 bp.")
  meta <- data$metadata
  column <- paste0(rank, "_group_id")
  if (!column %in% names(meta)) stop("UNITE rank grouping is unavailable")
  details <- data.frame(meta, version_id = data$manifest$id, pair_id = pair_id,
      group_id = meta[[column]], fit_state = state, amplicon_bp = length,
      forward_status = f$status, reverse_status = r$status, forward_start = f$start, forward_end = f$end,
      reverse_start = r$start, reverse_end = r$end, reason = unname(reasons[state]),
      length_basis = "encoded_oriented_reference_span_inclusive_primers", stringsAsFactors = FALSE)
  details$n_ambiguous_product_bases <- vapply(seq_len(nrow(details)), function(i) {
    if (state[i] != "exact_product") return(NA_integer_)
    as.integer(nchar(gsub("[ACGT]", "", substr(data$sequences[[meta$record_id[i]]], f$start[i], r$end[i]))))
  }, integer(1))
  groups <- unite_reference_groups(data, rank)
  fields <- c(n_exact_products="exact_product", n_single_sites="partial_sites", n_no_exact_sites="no_exact_sites",
      n_unresolved="unresolved_sites", n_unavailable="evidence_unavailable", n_withheld="product_withheld")
  for (field in names(fields)) {
    counted <- table(details$group_id[details$fit_state == fields[[field]]])
    groups[[field]] <- as.integer(counted[match(groups$group_id, names(counted))])
    groups[[field]][is.na(groups[[field]])] <- 0L
  }
  if (any(rowSums(groups[names(fields)]) != groups$n_selected)) stop("UNITE support states do not partition retained references")
  list(summary = groups, details = details, metadata = list(pair = pair, version_id = data$manifest$id,
      n_source = data$manifest$n_source, n_selected = nrow(meta), rank = rank))
}
