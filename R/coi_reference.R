# Runtime COI views read a single verified release or a complete, provenance-
# checked local scoring build. The BeePrime publication panel stays separate.

coi_reference_artifact_path <- function(artifact, root = getwd()) {
  if (is.null(artifact$sha256) || !grepl("^[a-fA-F0-9]{64}$", artifact$sha256)) {
    stop("COI reference artifact lacks a SHA-256 identity")
  }
  path <- artifact$local_path
  if (!is.null(path) && !grepl("^/", path)) path <- file.path(root, path)
  if (is.null(path) || !file.exists(path)) {
    if (is.null(artifact$url) || !nzchar(artifact$url)) stop("COI reference artifact is unavailable")
    path <- atlas_cached_artifact(artifact$url, artifact$sha256)
  }
  if (!identical(tolower(atlas_sha256(path)), tolower(artifact$sha256))) stop("COI reference artifact checksum failed")
  path
}

coi_reference_artifact <- function(artifact, root = getwd()) {
  path <- coi_reference_artifact_path(artifact, root)
  local <- artifact
  local$local_path <- path
  value <- atlas_read_artifact(local, root)
  n <- if (is.data.frame(value)) nrow(value) else if (is.character(value)) length(value) else NULL
  if (!is.null(n) && !is.null(artifact$rows) && n != artifact$rows) stop("COI reference artifact row count differs")
  value
}

.coi_reference_counts <- function(overview) {
  if (!all(c("order", "pair_id", "n_templates") %in% names(overview)) || !nrow(overview)) stop("COI overview is incomplete")
  if (anyDuplicated(paste(overview$order, overview$pair_id))) stop("COI overview repeats a group/pair")
  counts <- vapply(split(overview$n_templates, overview$order), function(n) {
    if (any(!is.finite(n)) || any(n < 1L) || length(unique(n)) != 1L) stop("COI overview has inconsistent group denominators")
    as.integer(n[1])
  }, integer(1))
  counts
}

.coi_reference_check_scores <- function(scores, overview, groups, records = NULL) {
  if (!all(c("order", "pair_id", "template") %in% names(scores))) stop("COI sequence scores lack template identities")
  relevant <- scores[scores$order %in% groups, , drop = FALSE]
  if (!nrow(relevant)) stop("COI sequence score partitions are empty")
  keys <- paste(relevant$order, relevant$pair_id)
  observed <- table(keys)
  expected <- overview[overview$order %in% groups, , drop = FALSE]
  expected_keys <- paste(expected$order, expected$pair_id)
  if (!setequal(names(observed), expected_keys) || any(as.integer(observed[expected_keys]) != expected$n_templates)) {
    stop("COI overview and sequence scores have different denominators")
  }
  accession <- if ("accession" %in% names(relevant)) relevant$accession else coi_accession(relevant$template)
  if (anyDuplicated(paste(keys, accession))) stop("COI sequence scores repeat an accession within a pair")
  if (!is.null(records) && any(!paste(relevant$order, accession) %in% paste(records$target_group, records$accession))) {
    stop("COI scores include accessions outside the frozen panel")
  }
  invisible(TRUE)
}

.coi_reference_check_growth <- function(growth, groups, counts, minimum_per_group, minimum_total) {
  columns <- c("target_group", "previous_retained", "added_new", "retained_after", "target_met", "verified_marker_records")
  if (!all(columns %in% names(growth)) || anyNA(growth[columns])) {
    stop("COI growth proof and the actual retained panel differ")
  }
  if (!is.logical(growth$target_met)) {
    flags <- tolower(trimws(as.character(growth$target_met)))
    if (any(!flags %in% c("true", "false"))) stop("COI growth proof has an invalid target flag")
    growth$target_met <- flags == "true"
  }
  if (anyDuplicated(growth$target_group) ||
      !setequal(growth$target_group, groups) || any(!growth$target_met) ||
      any(growth$retained_after != as.integer(counts[growth$target_group])) ||
      any(growth$retained_after != growth$previous_retained + growth$added_new) ||
      any(growth$previous_retained < 0L) || any(growth$added_new < 0L) ||
      any(growth$retained_after < ceiling(growth$previous_retained * 1.01)) ||
      any(growth$verified_marker_records > growth$retained_after) ||
      any(growth$verified_marker_records < minimum_per_group) || sum(growth$verified_marker_records) < minimum_total) {
    stop("COI growth proof and the actual retained panel differ")
  }
  invisible(TRUE)
}

.coi_reference_taxonomy <- function(value, groups, records = NULL) {
  if (is.null(value)) return(data.frame())
  if (!all(c("target_group", "accession") %in% names(value))) stop("COI taxonomy lacks group/accession identities")
  value <- as.data.frame(value, stringsAsFactors = FALSE)
  value <- value[value$target_group %in% groups, , drop = FALSE]
  keys <- paste(value$target_group, value$accession)
  if (anyDuplicated(keys)) stop("COI taxonomy repeats a group/accession")
  if (!is.null(records) && !setequal(keys, paste(records$target_group, records$accession))) {
    stop("COI taxonomy and retained sequences differ")
  }
  if (!"order" %in% names(value)) value$order <- ""
  value$taxonomic_order <- value$order
  # Acari and Collembola are display groups, not inferred modern orders.
  value$order <- value$target_group
  for (column in c("organism", "family", "subfamily", "genus", "species", "taxid", "taxonomy_status")) {
    if (!column %in% names(value)) value[[column]] <- rep(NA_character_, nrow(value))
    value[[column]] <- as.character(value[[column]])
    value[[column]][!is.na(value[[column]]) & !nzchar(value[[column]])] <- NA_character_
  }
  value
}

.coi_reference_preview <- function(root) {
  paths <- list.files(file.path(root, "data/primerminer/aligned_reference"),
    pattern = "_COX1_reference_aligned[.]fasta$", full.names = TRUE)
  if (!length(paths)) stop("No COI alignment preview is available")
  groups <- sub("_COX1_reference_aligned[.]fasta$", "", basename(paths))
  alignments <- lapply(paths, function(path) {
    values <- parse_fasta(path)
    values[!grepl("^(_R_)?NC_001322[.]1_COX1_reference([[:space:]]|$)", names(values))]
  })
  names(alignments) <- groups
  list(paths = setNames(paths, groups), alignments = alignments,
       counts = setNames(vapply(alignments, length, integer(1)), groups))
}

.coi_reference_full_release <- function(state, root, minimum_per_group, minimum_total) {
  manifest <- state$manifest
  retained <- manifest$coi_retained_state
  if (is.null(retained) || !length(retained$artifacts)) return(NULL)
  if (!identical(manifest$marker_id, "COI")) stop("COI reference manifest has the wrong marker")
  groups <- as.character(unlist(retained$target_groups, use.names = FALSE))
  if (!length(groups) || anyDuplicated(groups) || any(!nzchar(groups))) stop("COI retained target groups are invalid")
  orders_path <- file.path(root, "data/orders_22.csv")
  if (file.exists(orders_path)) {
    expected_groups <- read.csv(orders_path, stringsAsFactors = FALSE)$order
    expected_groups <- expected_groups[expected_groups != "Bivalvia"]
    if (!setequal(groups, expected_groups)) stop("COI retained release does not cover all catalog target groups")
  }
  artifacts <- retained$artifacts
  if (!all(c("records", "growth_qa") %in% names(artifacts))) stop("COI retained release is incomplete")
  records <- coi_reference_artifact(artifacts$records, root)
  if (!all(c("target_group", "accession", "sequence_sha256", "header_sha256") %in% names(records))) stop("COI retained record proof is incomplete")
  if (anyDuplicated(records$accession) || any(!records$target_group %in% groups)) stop("COI retained identities are inconsistent")
  for (column in c("sequence_sha256", "header_sha256")) {
    if (anyNA(records[[column]]) || any(!grepl("^[a-fA-F0-9]{64}$", records[[column]]))) stop("COI retained byte identities are invalid")
  }
  counts <- table(factor(records$target_group, levels = groups))
  if (sum(counts) != retained$retained_records || any(counts < minimum_per_group) || sum(counts) < minimum_total) {
    stop("COI retained reference does not meet its sequence targets")
  }
  if (!is.null(retained$group_counts)) {
    declared <- unlist(retained$group_counts, use.names = TRUE)
    if (!setequal(names(declared), groups) || any(as.numeric(declared[groups]) != as.numeric(counts))) stop("COI declared group counts differ from records")
  }
  reference_path <- file.path(root, "data/reference/coi_reference_NC_001322.1.fasta")
  reference <- parse_fasta(reference_path)
  reference_sha <- atlas_sha256(reference_path)
  alignments <- paths <- setNames(vector("list", length(groups)), groups)
  for (group in groups) {
    required <- paste0(c("raw_", "alignment_", "alignment_inputs_"), group)
    if (!all(required %in% names(artifacts))) stop("COI retained alignment or input proof is missing for ", group)
    artifact <- artifacts[[paste0("alignment_", group)]]
    path <- coi_reference_artifact_path(artifact, root)
    aligned <- coi_reference_artifact(artifact, root)
    inputs <- coi_reference_artifact(artifacts[[paste0("alignment_inputs_", group)]], root)
    if (!identical(inputs$evaluation_mode, "full") || !identical(inputs$reference_sha256, reference_sha) ||
        !identical(inputs$raw_sha256, artifacts[[paste0("raw_", group)]]$sha256)) stop("COI alignment inputs differ from the frozen release")
    if (!is.null(inputs$output_sha256) && !identical(inputs$output_sha256, artifact$sha256)) stop("COI alignment output differs from its input proof")
    is_reference <- grepl("^(_R_)?NC_001322[.]1_COX1_reference([[:space:]]|$)", names(aligned))
    if (sum(is_reference) != 1L || any(nchar(aligned) != 1536L) ||
        !identical(unname(aligned[is_reference]), unname(reference))) stop("COI coordinate reference or alignment width changed")
    templates <- aligned[!is_reference]
    accession <- coi_accession(names(templates))
    expected <- records$accession[records$target_group == group]
    if (anyDuplicated(accession) || !setequal(accession, expected) || length(accession) != counts[[group]]) stop("COI alignment lost retained accession identities")
    alignments[[group]] <- templates
    paths[[group]] <- path
  }
  overview <- coi_reference_artifact(manifest$artifacts$overview, root)
  scored_counts <- .coi_reference_counts(overview)
  if (any(!groups %in% names(scored_counts)) || any(scored_counts[groups] != as.integer(counts))) stop("COI overview and frozen alignments have different denominators")
  overview$reference_evidence <- ifelse(overview$order %in% groups,
    "full retained alignment", "archived release score; alignment not supplied")
  sequence_scores <- coi_reference_artifact(manifest$artifacts$pair_templates, root)
  .coi_reference_check_scores(sequence_scores, overview, groups, records)
  taxonomy <- if ("panel_taxonomy" %in% names(artifacts)) {
    .coi_reference_taxonomy(coi_reference_artifact(artifacts$panel_taxonomy, root), groups, records)
  } else data.frame()
  growth <- coi_reference_artifact(artifacts$growth_qa, root)
  .coi_reference_check_growth(growth, groups, counts, minimum_per_group, minimum_total)
  list(mode = "full_release", alignment_mode = "full", state = state,
    version = manifest$release_version, overview = overview, paths = unlist(paths), alignments = alignments,
    counts = setNames(as.integer(counts), groups), n_references = as.integer(sum(counts)),
    taxonomy = taxonomy, growth_qa = growth, records = records, pair_templates = sequence_scores,
    n_verified_marker = as.integer(sum(growth$verified_marker_records)),
    policy = retained$policy, cache = new.env(parent = emptyenv()), issue = NULL)
}

.coi_reference_local_full <- function(state, root, minimum_per_group, minimum_total) {
  path <- file.path(root, "data/provenance/order_scoring_manifest.csv")
  if (!file.exists(path)) return(NULL)
  installed <- read.csv(path, stringsAsFactors = FALSE)
  groups <- read.csv(file.path(root, "data/orders_22.csv"), stringsAsFactors = FALSE)$order
  groups <- groups[groups != "Bivalvia"]
  inputs <- coi_scoring_inputs(root, "full", groups, minimum_per_group, minimum_total)
  columns <- c("order", "n_templates", "alignment_sha256", "raw_sha256", "reference_sha256", "alignment_input_sha256")
  if (!all(columns %in% names(installed)) || anyDuplicated(installed$order) || !setequal(installed$order, groups)) stop("Local COI scoring provenance is incomplete")
  installed <- installed[match(groups, installed$order), , drop = FALSE]
  expected <- inputs$manifest[match(groups, inputs$manifest$order), , drop = FALSE]
  for (column in columns) if (!identical(as.character(installed[[column]]), as.character(expected[[column]]))) stop("Local COI scores are stale for ", column)
  for (name in c("primer_catalog_sha256", "binding_coordinates_sha256")) {
    source <- if (name == "primer_catalog_sha256") "data/primers.csv" else "data/derived/primer_bindings.csv"
    if (!name %in% names(installed) || any(installed[[name]] != atlas_sha256(file.path(root, source)))) stop("Local COI score catalog or binding coordinates changed")
  }
  overview_path <- file.path(root, "data/derived/order_pair_summary.csv")
  overview <- read.csv(overview_path, stringsAsFactors = FALSE)
  counts <- .coi_reference_counts(overview)
  if (!setequal(names(counts), groups) || any(counts[groups] != expected$n_templates)) stop("Local COI overview and alignments differ")
  if (!"evidence_class" %in% names(overview) || any(overview$evidence_class != "full_order_reference_alignment")) stop("Local COI overview is not full-sequence evidence")
  overview$reference_evidence <- "full retained alignment · local verified build"
  records <- coi_retained_records(root, groups)
  if (anyDuplicated(records$accession)) stop("Local COI reference repeats an accession across target groups")
  sequence_scores <- read.csv(file.path(root, "data/derived/pair_template_scores.csv"), stringsAsFactors = FALSE)
  .coi_reference_check_scores(sequence_scores, overview, groups, records)
  artifact <- function(file, rows = NULL) {
    list(local_path = file, sha256 = atlas_sha256(file), rows = rows)
  }
  local_state <- state
  local_state$manifest$release_version <- paste0("local-", max(installed$scored_at_utc))
  local_state$manifest$coi_retained_state <- NULL
  local_state$source <- path
  local_state$stale <- FALSE
  for (name in c("overview", "pair_templates")) {
    file <- if (name == "overview") overview_path else file.path(root, "data/derived/pair_template_scores.csv")
    if (!file.exists(file)) stop("Local COI score table is missing")
    local_state$manifest$artifacts[[name]] <- artifact(file)
  }
  # Local score partitions replace production artifacts as a complete cohort.
  old_general <- grepl("^(full_order_|primer_positions$)", names(local_state$manifest$artifacts))
  local_state$manifest$artifacts[old_general] <- NULL
  pair_ids <- unique(overview$pair_id)
  for (pair in pair_ids) {
    prefix <- paste0("full_order_", tolower(pair))
    file <- file.path(root, "data/derived", paste0(prefix, "_sequence_scores.csv.gz"))
    if (!file.exists(file)) stop("Local COI per-pair scores are incomplete")
    local_state$manifest$artifacts[[prefix]] <- artifact(file)
    positions <- file.path(root, "data/derived", paste0(prefix, "_primer_positions.csv.gz"))
    if (!file.exists(positions)) stop("Local COI per-pair position traces are incomplete")
    local_state$manifest$artifacts[[paste0(prefix, "_positions")]] <- artifact(positions)
  }
  alignments <- lapply(inputs$files, function(file) {
    values <- parse_fasta(file)
    values[!grepl("^(_R_)?NC_001322[.]1_COX1_reference([[:space:]]|$)", names(values))]
  })
  names(alignments) <- groups
  taxonomy_path <- file.path(root, "data/derived/coi_reference_panel/taxonomy.csv")
  taxonomy <- if (file.exists(taxonomy_path)) .coi_reference_taxonomy(read.csv(taxonomy_path, stringsAsFactors = FALSE), groups, records) else data.frame()
  qa_path <- file.path(root, "data/provenance/coi_reference_growth_qa.csv")
  growth <- if (file.exists(qa_path)) read.csv(qa_path, stringsAsFactors = FALSE) else data.frame()
  .coi_reference_check_growth(growth, groups, counts, minimum_per_group, minimum_total)
  list(mode = "full_local", alignment_mode = "full", state = local_state,
    version = local_state$manifest$release_version, overview = overview, paths = setNames(inputs$files, groups),
    alignments = alignments, counts = counts, n_references = as.integer(sum(counts)), taxonomy = taxonomy,
    pair_templates = sequence_scores, records = records,
    growth_qa = growth, n_verified_marker = if ("verified_marker_records" %in% names(growth)) as.integer(sum(growth$verified_marker_records)) else NA_integer_,
    policy = list(minimum_per_group = minimum_per_group, minimum_total = minimum_total, monthly_growth = .01),
    cache = new.env(parent = emptyenv()), issue = NULL)
}

coi_reference_load <- function(state, root = getwd(), minimum_per_group = 100L, minimum_total = 10000L) {
  issue <- NULL
  full <- tryCatch(.coi_reference_full_release(state, root, minimum_per_group, minimum_total), error = identity)
  if (!inherits(full, "error") && !is.null(full)) return(full)
  if (inherits(full, "error")) issue <- conditionMessage(full)
  local <- tryCatch(.coi_reference_local_full(state, root, minimum_per_group, minimum_total), error = identity)
  if (!inherits(local, "error") && !is.null(local)) {
    local$issue <- issue
    return(local)
  }
  if (inherits(local, "error")) issue <- paste(c(issue, conditionMessage(local)), collapse = "; ")
  # Invalid new releases fall back to one entire bundled release. A missing
  # artifact never silently borrows a table from a different release.
  if (!is.null(state$manifest$coi_retained_state) && !is.null(issue)) {
    fallback <- atlas_manifest_path("COI", root)
    if (!file.exists(fallback)) stop(issue)
    state <- list(manifest = jsonlite::read_json(fallback, simplifyVector = FALSE),
                  source = fallback, stale = TRUE, remote_error = issue)
    fallback_full <- tryCatch(.coi_reference_full_release(state, root, minimum_per_group, minimum_total), error = identity)
    if (!inherits(fallback_full, "error") && !is.null(fallback_full)) {
      fallback_full$issue <- issue
      return(fallback_full)
    }
    if (!is.null(state$manifest$coi_retained_state)) stop(issue)
  }
  preview <- .coi_reference_preview(root)
  overview <- tryCatch(coi_reference_artifact(state$manifest$artifacts$overview, root), error = identity)
  if (inherits(overview, "error")) stop("The chosen COI overview is unavailable: ", conditionMessage(overview))
  counts <- .coi_reference_counts(overview)
  same_preview <- setequal(names(counts), names(preview$counts)) && all(counts[names(preview$counts)] == preview$counts)
  overview$reference_evidence <- if (same_preview) "legacy centroid preview" else "published full-score overview; frozen alignment not supplied"
  list(mode = if (same_preview) "legacy_preview" else "published_scores_legacy_preview",
    alignment_mode = "legacy_centroid_preview", state = state, version = state$manifest$release_version,
    overview = overview, paths = preview$paths, alignments = preview$alignments,
    counts = counts, n_references = as.integer(sum(counts)), taxonomy = data.frame(), growth_qa = data.frame(),
    preview_count = sum(preview$counts), policy = NULL, cache = new.env(parent = emptyenv()), issue = issue)
}

coi_reference_status <- function(reference) {
  count <- format(reference$n_references, big.mark = ",")
  if (reference$mode %in% c("full_release", "full_local")) {
    scope <- if (all(c("Acari", "Collembola") %in% names(reference$counts)) && length(reference$counts) == 22L) {
      "22 target groups (20 orders plus the Acari and Collembola composite groups)"
    } else paste0(length(reference$counts), " target groups")
    text <- paste0(count, " retained COI records across ", scope,
      " · ", if (reference$mode == "full_local") "verified local build; publication pending" else paste0("release ", reference$version),
      ". Full alignments and scoring denominators agree. Minimum ", min(reference$counts), " references per group.")
    if (length(reference$n_verified_marker) && is.finite(reference$n_verified_marker)) {
      text <- paste0(text, " ", format(reference$n_verified_marker, big.mark = ","), " records have verified COI feature evidence; inherited records retain their provenance labels.")
    }
  } else if (reference$mode == "published_scores_legacy_preview") {
    text <- paste0(count, " published COI scoring references · release ", reference$version,
      ". This older release has no frozen full alignments; custom placement uses an explicitly separate ",
      format(reference$preview_count, big.mark = ","), "-template legacy centroid preview.")
  } else {
    text <- paste0(count, " legacy COI centroid preview templates · release ", reference$version,
      ". A complete full-sequence panel requires at least 10,000 marker-verified COI records overall and 100 per target group; inherited records with unresolved marker QC are reported separately.")
  }
  if (!is.null(reference$issue) && nzchar(reference$issue)) text <- paste0(text, " A newer reference could not be verified: ", reference$issue, ".")
  growth <- reference$growth_qa
  if (nrow(growth) && all(c("previous_retained", "added_new", "retained_after") %in% names(growth))) {
    text <- paste0(text, " This build reports ", format(sum(growth$previous_retained), big.mark = ","),
      " previous records and ", format(sum(growth$added_new), big.mark = ","), " additions; counts come from the retained panel.")
  }
  text
}

.coi_reference_empty_positions <- function() {
  data.frame(order = character(), pair_id = character(), pair_label = character(), primer_name = character(),
    direction = character(), template = character(), alignment_position = integer(), primer_position_5_to_3 = integer(),
    distance_from_3prime = integer(), primer_base = character(), template_base = character(),
    position_penalty = numeric(), mismatch = logical(), unavailable_reason = character(), stringsAsFactors = FALSE)
}

coi_reference_pair_sequences <- function(reference, pair_id, root = getwd()) {
  if (length(pair_id) != 1L || is.na(pair_id) || !nzchar(pair_id)) return(data.frame())
  name <- if (identical(pair_id, "BEEPRIME")) "exact_beeprime" else paste0("full_order_", tolower(pair_id))
  artifact <- reference$state$manifest$artifacts[[name]]
  if (is.null(artifact)) return(data.frame())
  key <- paste0("scores_", pair_id, "_", artifact$sha256)
  if (!identical(pair_id, "BEEPRIME") && exists(key, envir = reference$cache, inherits = FALSE)) {
    return(get(key, envir = reference$cache))
  }
  scores <- coi_reference_artifact(artifact, root)
  scores <- coi_reference_check_pair_scores(reference, scores, pair_id)
  if (!identical(pair_id, "BEEPRIME")) {
    scores <- coi_reference_pair_taxonomy(reference, scores)
    previous <- ls(reference$cache, pattern = "^scores_")
    if (length(previous) >= 2L) rm(list = previous[1L], envir = reference$cache)
    assign(key, scores, envir = reference$cache)
  }
  scores
}

coi_reference_pair_positions <- function(reference, pair_id, maximum_legacy_rows = 2000000L) {
  if (length(pair_id) != 1L || is.na(pair_id) || !nzchar(pair_id) || grepl("^CUSTOM_", pair_id)) return(.coi_reference_empty_positions())
  key <- paste0("positions_", pair_id)
  if (exists(key, envir = reference$cache, inherits = FALSE)) return(get(key, envir = reference$cache))
  artifacts <- reference$state$manifest$artifacts
  prefix <- paste0("full_order_", tolower(pair_id))
  keys <- c(paste0(prefix, "_positions"), paste0(prefix, "_primer_positions"))
  selected <- keys[keys %in% names(artifacts)]
  if (length(selected)) {
    value <- coi_reference_artifact(artifacts[[selected[1]]])
    if (!all(c("pair_id", "order", "template") %in% names(value)) || any(value$pair_id != pair_id)) stop("COI position partition has the wrong pair identity")
    if (reference$alignment_mode == "full" && nrow(reference$records)) {
      membership <- paste(reference$records$target_group, reference$records$accession)
      observed <- paste(value$order, coi_accession(value$template))
      if (any(!observed %in% membership)) stop("COI position partition is outside the selected frozen panel")
      if (!setequal(observed, membership)) stop("COI position partition does not cover every template in the selected frozen panel")
      if ("alignment_position" %in% names(value) &&
          any(!is.finite(value$alignment_position) | value$alignment_position < 1L | value$alignment_position > 1536L)) {
        stop("COI position partition contains invalid alignment columns")
      }
    }
  } else {
    legacy <- artifacts$primer_positions
    if (reference$alignment_mode == "full" || is.null(legacy) || is.null(legacy$rows) || legacy$rows > maximum_legacy_rows) {
      empty <- .coi_reference_empty_positions()
      attr(empty, "position_issue") <- "A verified position partition is unavailable for this pair in the selected release."
      return(empty)
    }
    value <- coi_reference_artifact(legacy)
    value <- value[value$pair_id == pair_id, , drop = FALSE]
  }
  if (!"unavailable_reason" %in% names(value)) value$unavailable_reason <- rep("", nrow(value))
  value$unavailable_reason[is.na(value$unavailable_reason)] <- ""
  # Bound runtime memory when a user opens many pairs in one worker.
  previous <- ls(reference$cache, pattern = "^positions_")
  if (length(previous) >= 2L) rm(list = previous[1], envir = reference$cache)
  assign(key, value, envir = reference$cache)
  value
}

coi_reference_pair_taxonomy <- function(reference, scores) {
  if (!nrow(scores)) return(scores)
  taxonomy <- reference$taxonomy
  if (nrow(taxonomy) && all(c("order", "accession") %in% names(scores))) {
    matched <- match(paste(scores$order, scores$accession), paste(taxonomy$target_group, taxonomy$accession))
    for (column in setdiff(names(taxonomy), c("target_group", "accession", "order"))) {
      if (!column %in% names(scores)) scores[[column]] <- taxonomy[[column]][matched]
    }
  }
  for (column in c("family", "subfamily", "genus", "species", "organism")) {
    if (!column %in% names(scores)) scores[[column]] <- rep(NA_character_, nrow(scores))
  }
  if (!"organism_label" %in% names(scores)) scores$organism_label <- scores$organism
  if (!"scientific_name_ncbi" %in% names(scores)) scores$scientific_name_ncbi <- scores$organism
  if (!"forward_penalty" %in% names(scores)) scores$forward_penalty <- scores$forward_primerminer_penalty
  if (!"reverse_penalty" %in% names(scores)) scores$reverse_penalty <- scores$reverse_primerminer_penalty
  if (!"cluster_size" %in% names(scores)) scores$cluster_size <- rep(1L, nrow(scores))
  if (!"terminal_3_mismatch_pair" %in% names(scores)) {
    scores$terminal_3_mismatch_pair <- (scores$forward_terminal_3_mismatches + scores$reverse_terminal_3_mismatches) > 0
  }
  scores$reference_mode <- "full"
  scores$reference_version <- reference$version
  scores
}

coi_reference_check_pair_scores <- function(reference, scores, pair_id) {
  if (!all(c("pair_id", "order", "template") %in% names(scores)) ||
      anyNA(scores$pair_id) || any(scores$pair_id != pair_id)) {
    stop("COI sequence partition has the wrong pair identity")
  }
  if (reference$alignment_mode == "full" && !identical(pair_id, "BEEPRIME")) {
    accession <- coi_accession(scores$template)
    if ("accession" %in% names(scores) && anyNA(scores$accession)) stop("COI sequence partition lacks accession identities")
    if ("accession" %in% names(scores) && any(scores$accession != accession)) stop("COI sequence partition has inconsistent accession identities")
    scores$accession <- accession
    membership <- paste(reference$records$target_group, reference$records$accession)
    if (any(!paste(scores$order, scores$accession) %in% membership)) stop("COI sequence partition is outside the selected frozen panel")
    overview <- reference$overview[reference$overview$pair_id == pair_id, , drop = FALSE]
    .coi_reference_check_scores(scores, overview, names(reference$counts), reference$records)
  }
  scores
}

coi_reference_taxon_filter <- function(scores, rank, value) {
  if (!rank %in% names(scores) || length(value) != 1L || is.na(value)) return(scores[FALSE, , drop = FALSE])
  values <- as.character(scores[[rank]])
  unresolved <- is.na(values) | !nzchar(trimws(values))
  keep <- if (identical(as.character(value), paste0(rank, " unresolved"))) unresolved else !unresolved & values == value
  scores[keep, , drop = FALSE]
}

coi_reference_taxon_summary <- function(scores, group, rank) {
  if (!rank %in% c("family", "subfamily", "genus")) stop("Unsupported COI taxonomy rank")
  rows <- scores[scores$order == group, , drop = FALSE]
  if (!nrow(rows)) return(data.frame())
  groups <- c("pair_id", "pair_label", "order", "family")
  if (rank != "family") groups <- c(groups, rank)
  summary <- summarize_expanded_primer_scores(rows, groups)
  summary$taxon <- summary[[rank]]
  summary$context <- if (rank == "family") summary$order else paste(summary$order, summary$family, sep = " · ")
  summary
}

coi_reference_selected_alignments <- function(reference, accessions) {
  if (reference$alignment_mode != "full") stop("Frozen full COI alignments are unavailable")
  sequences <- do.call(c, unname(reference$alignments))
  keys <- coi_accession(names(sequences))
  selected <- match(unique(as.character(accessions)), keys)
  if (anyNA(selected)) stop("Selected COI accessions are outside the frozen alignment")
  sequences[selected]
}

coi_reference_geography <- function(reference) {
  accessions <- unique(as.character(reference$taxonomy$accession))
  data.frame(accession = accessions, country_or_territory = rep(NA_character_, length(accessions)),
    locality = rep(NA_character_, length(accessions)), normalized_location = rep(NA_character_, length(accessions)),
    geographic_resolution = rep("unresolved", length(accessions)), latitude = rep(NA_real_, length(accessions)),
    longitude = rep(NA_real_, length(accessions)), stringsAsFactors = FALSE)
}
