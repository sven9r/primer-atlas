# Promotion compares against the pointer that is hosted now, rather than only
# the release chosen at the beginning of a worker. No stale/local latest fallback
# is acceptable here. The workflow serializes publishers separately per marker.

coi_promotion_public_url <- function(url) {
  if (length(url) != 1L || is.na(url) || !grepl("^https://[^/?#@[:space:]]+(/[^?#[:space:]]*)?$", url)) {
    stop("COI promotion requires an anonymous public HTTPS URL")
  }
  url
}

coi_promotion_fetch_latest <- function(url) {
  coi_promotion_public_url(url)
  response <- tryCatch(httr2::request(url) |>
    httr2::req_headers(`Cache-Control` = "no-cache", Pragma = "no-cache") |>
    httr2::req_timeout(30) |> httr2::req_perform(), error = function(e) NULL)
  if (is.null(response) || httr2::resp_status(response) != 200L) {
    stop("Cannot verify the currently hosted COI latest pointer; promotion refused")
  }
  bytes <- httr2::resp_body_raw(response)
  manifest <- tryCatch(jsonlite::fromJSON(rawToChar(bytes), simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(manifest)) stop("Currently hosted COI latest pointer is invalid JSON")
  list(manifest = manifest, stale = FALSE, source = url,
       manifest_sha256 = digest::digest(bytes, algo = "sha256", serialize = FALSE))
}

coi_promotion_artifact_path <- function(artifact, root, hosted) {
  coi_promotion_public_url(artifact$url)
  path <- if (isTRUE(hosted)) NULL else artifact$local_path
  if (!is.null(path) && !startsWith(path, "/")) path <- file.path(root, path)
  if (is.null(path) || !file.exists(path)) path <- atlas_cached_artifact(artifact$url, artifact$sha256)
  path
}

coi_verify_promotion <- function(candidate, root = getwd(), latest_url = NULL,
                                 fetch_latest = coi_promotion_fetch_latest,
                                 artifact_path = coi_promotion_artifact_path) {
  if (!identical(candidate$marker_id, "COI") || is.null(candidate$coi_retained_state)) stop("Candidate lacks a complete COI retained state")
  if (is.null(latest_url)) {
    base <- sub("/+$", "", Sys.getenv("ATLAS_MANIFEST_BASE_URL", ""))
    if (!nzchar(base)) {
      first_url <- candidate$coi_retained_state$artifacts$ledger$url
      base <- paste0(strsplit(first_url, "/releases/", fixed = TRUE)[[1]][1], "/releases")
    }
    latest_url <- paste0(base, "/COI/latest.json")
  }
  coi_promotion_public_url(latest_url)
  current_state <- fetch_latest(latest_url)
  if (!identical(current_state$stale, FALSE) || !identical(current_state$source, latest_url)) {
    stop("COI promotion refuses stale or local latest-pointer fallback")
  }
  if (length(current_state$manifest_sha256) != 1L ||
      !grepl("^[a-fA-F0-9]{64}$", current_state$manifest_sha256)) stop("Hosted latest pointer lacks verified byte identity")
  current <- current_state$manifest
  if (!identical(current$marker_id, "COI") || is.null(current$release_version)) stop("Hosted latest is not a COI release")

  checked <- new.env(parent = emptyenv())
  path_for <- function(artifact, hosted) {
    if (is.null(artifact) || length(artifact$sha256) != 1L ||
        !grepl("^[a-fA-F0-9]{64}$", artifact$sha256)) stop("Promotion artifact lacks a SHA-256 checksum")
    coi_promotion_public_url(artifact$url)
    key <- paste(hosted, artifact$url, artifact$sha256, artifact$local_path, sep = "|")
    if (length(key) != 1L) key <- paste(hosted, artifact$url, artifact$sha256, sep = "|")
    if (exists(key, checked, inherits = FALSE)) return(get(key, checked, inherits = FALSE))
    path <- artifact_path(artifact, root, hosted)
    if (!file.exists(path) || !identical(tolower(atlas_sha256(path)), tolower(artifact$sha256))) {
      stop("COI promotion artifact checksum failed")
    }
    assign(key, path, checked); path
  }
  read_table <- function(artifact, hosted) {
    path <- path_for(artifact, hosted)
    if (grepl("[.]parquet$", path)) as.data.frame(nanoparquet::read_parquet(path), stringsAsFactors = FALSE)
    else read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  key_for <- function(rows) paste(rows$target_group, rows$accession, sep = "|")
  validate_state <- function(manifest, hosted) {
    state <- manifest$coi_retained_state
    artifacts <- state$artifacts
    groups <- unlist(state$target_groups, use.names = FALSE)
    if (!length(groups) || anyDuplicated(groups) ||
        !all(c("records", "ledger", paste0("raw_", groups)) %in% names(artifacts))) {
      stop("COI promotion requires the ledger and actual raw sequences for every group")
    }
    records <- do.call(rbind, lapply(groups, function(group) {
      artifact <- artifacts[[paste0("raw_", group)]]
      sequences <- coi_read_fasta(path_for(artifact, hosted))
      if (!is.null(artifact$rows) && artifact$rows != length(sequences)) stop("Hosted raw COI count differs from descriptor")
      data.frame(target_group = group, accession = coi_accession(names(sequences)),
        sequence_sha256 = vapply(unname(sequences), digest::digest, character(1), algo = "sha256", serialize = FALSE),
        header_sha256 = vapply(names(sequences), digest::digest, character(1), algo = "sha256", serialize = FALSE),
        stringsAsFactors = FALSE)
    }))
    if (anyDuplicated(key_for(records)) || nrow(records) != state$retained_records) stop("COI retained-state count/identity mismatch")
    for (group in groups) if (is.null(state$group_counts[[group]]) || state$group_counts[[group]] != sum(records$target_group == group)) stop("COI group count differs from actual raw sequences")
    ledger <- read_table(artifacts$ledger, hosted)
    archived <- read_table(artifacts$records, hosted)
    for (rows in list(ledger, archived)) {
      if (!all(c("target_group", "accession", "sequence_sha256", "header_sha256") %in% names(rows)) ||
          anyDuplicated(key_for(rows)) || !setequal(key_for(rows), key_for(records))) stop("COI promotion ledger/raw membership differs")
      coi_assert_retention(rows, records)
    }
    list(records = records, ledger = ledger, state = state)
  }
  next_state <- validate_state(candidate, FALSE)
  source_count <- 0L
  if (!is.null(current$coi_retained_state)) {
    prior <- validate_state(current, TRUE)
    coi_assert_retention(prior$records, next_state$records)
    prior_keys <- key_for(prior$ledger)
    matches <- match(prior_keys, key_for(next_state$ledger))
    # Preserve scientific source identity of every retained row. New taxonomic
    # reinterpretations need a separate versioned product, not an overwritten row.
    columns <- intersect(c("source_sha256", "taxonomy_source_sha256", "first_release_id",
      "added_at_utc", "baseline_source", "published_release_id", "published_artifact_sha256",
      "published_artifact_url"), names(prior$ledger))
    for (column in columns) {
      if (!column %in% names(next_state$ledger) ||
          !identical(as.character(prior$ledger[[column]]), as.character(next_state$ledger[[column]][matches]))) {
        stop("Previously hosted COI provenance changed: ", column)
      }
    }
    prior_artifacts <- prior$state$artifacts
    next_artifacts <- next_state$state$artifacts
    immutable_sources <- Filter(function(x) !is.null(x$restore_path) &&
      grepl("^data/derived/coi_reference_panel/(raw/content/|legacy_raw/)", x$restore_path), prior_artifacts)
    needed <- unique(c(vapply(immutable_sources, `[[`, character(1), "sha256"),
      unlist(prior$ledger[intersect(c("source_sha256", "taxonomy_source_sha256"), names(prior$ledger))], use.names = FALSE)))
    needed <- tolower(needed[!is.na(needed) & nzchar(needed)])
    if (any(!grepl("^[a-f0-9]{64}$", needed))) stop("Hosted COI scientific source hash is malformed")
    prior_available <- tolower(vapply(prior_artifacts, `[[`, character(1), "sha256"))
    if (!all(needed %in% prior_available)) stop("Hosted COI ledger refers to an unavailable scientific source archive")
    for (hash in needed) path_for(prior_artifacts[[match(hash, prior_available)]], TRUE)
    available <- tolower(vapply(next_artifacts, `[[`, character(1), "sha256"))
    if (!all(needed %in% available)) stop("A previously hosted COI source archive was dropped")
    for (hash in needed) path_for(next_artifacts[[match(hash, available)]], FALSE)
    source_count <- length(needed)
    retained <- nrow(prior$records)
    evidence <- "exact retained sequences, headers and scientific source archives"
  } else {
    artifact <- current$artifacts$full_order_anml
    if (is.null(artifact)) artifact <- current$artifacts$pair_templates
    if (is.null(artifact)) stop("Legacy hosted COI release has no recoverable accession membership")
    rows <- read_table(artifact, TRUE)
    group_column <- if ("target_group" %in% names(rows)) "target_group" else if ("order" %in% names(rows)) "order" else NULL
    accession_column <- if ("template" %in% names(rows)) "template" else if ("accession" %in% names(rows)) "accession" else NULL
    if (is.null(group_column) || is.null(accession_column)) stop("Legacy hosted COI membership lacks group/accession identities")
    membership <- data.frame(target_group = rows[[group_column]], accession = coi_accession(rows[[accession_column]]), stringsAsFactors = FALSE)
    membership <- membership[membership$target_group %in% unique(next_state$records$target_group) &
      nzchar(membership$accession) & !grepl("^NC_001322[.]1_COX1_reference", membership$accession), , drop = FALSE]
    membership <- membership[!duplicated(key_for(membership)), , drop = FALSE]
    if (!nrow(membership) || !setequal(unique(membership$target_group), unique(next_state$records$target_group))) stop("Legacy hosted membership is incomplete")
    if (!all(key_for(membership) %in% key_for(next_state$records))) stop("A currently published legacy COI accession was dropped")
    baseline <- read_table(next_state$state$artifacts$published_baseline, FALSE)
    columns <- c("target_group", "accession", "published_release_id", "published_artifact_sha256", "published_artifact_url")
    if (!all(columns %in% names(baseline)) || anyDuplicated(key_for(baseline))) stop("Candidate lacks frozen legacy publication source identity")
    matches <- match(key_for(membership), key_for(baseline))
    if (anyNA(matches) || anyNA(baseline[matches, columns]) ||
        any(baseline$published_release_id[matches] != current$release_version) ||
        any(tolower(baseline$published_artifact_sha256[matches]) != tolower(artifact$sha256)) ||
        any(baseline$published_artifact_url[matches] != artifact$url)) stop("Candidate legacy membership provenance differs from the currently hosted release")
    retained <- nrow(membership)
    evidence <- "legacy accession-version membership and publication source identity; original raw bytes unavailable"
  }
  list(candidate_release_id = candidate$release_version, current_release_id = current$release_version,
    current_manifest_sha256 = current_state$manifest_sha256, retained_from_current = retained,
    candidate_retained_records = nrow(next_state$records), preserved_source_archives = source_count,
    evidence = evidence, verified_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
}
