#!/usr/bin/env Rscript
# Read-only proof across two immutable public releases. For example:
# ATLAS_MANIFEST_BASE_URL=https://public-host/releases Rscript \
#   scripts/verify_coi_hosted_retention.R first-release-id second-release-id
.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
source("R/data_layer.R")
source("R/coi_retention.R")
source("R/coi_promotion.R")
coi_hosted_artifact_path <- function(artifact, root, hosted) {
  # Even a candidate whose original GHA local_path exists must be downloaded
  # anonymously for this independent hosted-retention audit.
  coi_promotion_artifact_path(artifact, root, hosted = TRUE)
}
coi_hosted_read_artifact <- function(artifact) {
  artifact$local_path <- NULL
  atlas_read_artifact(artifact)
}
coi_hosted_retention_details <- function(first, second, read_hosted = coi_hosted_read_artifact) {
  old <- first$coi_retained_state
  new <- second$coi_retained_state
  if (is.null(old)) stop("The first release must contain full retained state for a two-build byte-retention proof")
  groups <- unlist(old$target_groups, use.names = FALSE)
  if (!setequal(groups, unlist(new$target_groups, use.names = FALSE))) stop("The two releases use different target groups")
  old_counts <- vapply(groups, function(group) old$group_counts[[group]], numeric(1))
  new_counts <- vapply(groups, function(group) new$group_counts[[group]], numeric(1))
  growth <- new$policy$monthly_growth
  if (is.null(growth) || !is.finite(growth) || growth < 0.01 || any(new_counts < ceiling(old_counts * (1 + growth)))) {
    stop("The second hosted release does not add the required real records in every group")
  }
  old_qa <- as.data.frame(read_hosted(old$artifacts$growth_qa), stringsAsFactors = FALSE)
  new_qa <- as.data.frame(read_hosted(new$artifacts$growth_qa), stringsAsFactors = FALSE)
  if (!all(c("target_group", "previous_retained", "retained_after", "verified_marker_records") %in% names(new_qa)) ||
      !setequal(new_qa$target_group, groups) || anyDuplicated(new_qa$target_group)) stop("New hosted QA lacks verified counts for each group")
  new_qa <- new_qa[match(groups, new_qa$target_group), , drop = FALSE]
  verified <- new_qa$verified_marker_records
  if (anyNA(verified) || any(!is.finite(verified)) || any(verified != as.integer(verified)) ||
      any(verified < 100L) || sum(verified) < 10000L || any(verified > new_counts) ||
      anyNA(new_qa$previous_retained) || anyNA(new_qa$retained_after) ||
      any(new_qa$previous_retained != old_counts) || any(new_qa$retained_after != new_counts)) {
    stop("Hosted growth QA does not reconcile to actual retained counts and verified marker floors")
  }
  if (!all(c("target_group", "retained_after", "verified_marker_records") %in% names(old_qa)) ||
      !setequal(old_qa$target_group, groups) || anyDuplicated(old_qa$target_group)) stop("First hosted QA lacks verified counts for each group")
  old_qa <- old_qa[match(groups, old_qa$target_group), , drop = FALSE]
  if (anyNA(old_qa$verified_marker_records) || any(!is.finite(old_qa$verified_marker_records)) ||
      any(old_qa$verified_marker_records < 100L) || sum(old_qa$verified_marker_records) < 10000L ||
      any(old_qa$verified_marker_records > old_counts) || anyNA(old_qa$retained_after) ||
      any(old_qa$retained_after != old_counts)) stop("First hosted verified QA does not reconcile to its retained counts")
  old_run <- read_hosted(old$artifacts$panel_run_manifest)
  new_run <- read_hosted(new$artifacts$panel_run_manifest)
  if (!identical(old_run$release_id, first$release_version) || !identical(new_run$release_id, second$release_version) ||
      !identical(new_run$previous_run_manifest_sha256, old$artifacts$panel_run_manifest$sha256)) {
    stop("Second acquisition does not cite the exact first frozen acquisition manifest")
  }
  if (!identical(old_run$marker_id, "COI") || !identical(new_run$marker_id, "COI") ||
      old_run$retained_count != sum(old_counts) || new_run$retained_count != sum(new_counts) ||
      old_run$verified_marker_count != sum(old_qa$verified_marker_records) || new_run$verified_marker_count != sum(verified)) {
    stop("Frozen acquisition manifests do not reconcile to actual retained/verified QA counts")
  }
  for (group in groups) {
    old_inputs <- read_hosted(old$artifacts[[paste0("alignment_inputs_", group)]])
    new_inputs <- read_hosted(new$artifacts[[paste0("alignment_inputs_", group)]])
    raw_name <- paste0("raw_", group)
    alignment_name <- paste0("alignment_", group)
    if (!identical(old_inputs$raw_sha256, old$artifacts[[raw_name]]$sha256) ||
        !identical(new_inputs$raw_sha256, new$artifacts[[raw_name]]$sha256) ||
        !identical(new_inputs$output_sha256, new$artifacts[[alignment_name]]$sha256) ||
        identical(new_inputs$input_sha256, old_inputs$input_sha256) ||
        new_inputs$aligned_templates != new_counts[[group]]) stop("Hosted alignment inputs/counts did not follow the grown raw panel: ", group)
  }
  list(real_added_records = sum(new_counts - old_counts), per_group_added = as.list(new_counts - old_counts),
    all_group_growth_targets_met = TRUE, all_alignment_input_hashes_changed = TRUE,
    acquisition_lineage_verified = TRUE, previous_panel_run_manifest_sha256 = old$artifacts$panel_run_manifest$sha256,
    first_verified_marker_records = sum(old_qa$verified_marker_records),
    second_verified_marker_records = sum(verified), second_minimum_verified_per_group = min(verified))
}

coi_hosted_retention_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2L || any(!grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", args)) || any(args == "latest")) {
    stop("Provide the first and second immutable COI release IDs")
  }
  base <- sub("/+$", "", Sys.getenv("ATLAS_MANIFEST_BASE_URL", ""))
  if (!nzchar(base)) stop("ATLAS_MANIFEST_BASE_URL is required")
  urls <- paste0(base, "/COI/", args, "/manifest.json")
  states <- lapply(urls, coi_promotion_fetch_latest)
  for (i in seq_along(args)) {
    if (!identical(states[[i]]$manifest$release_version, args[i])) stop("Public manifest release identity differs from its requested immutable path")
  }
  # The injected fetch uses exact frozen first-release bytes, while both the
  # previous and candidate artifacts always use anonymous public downloads.
  proof <- coi_verify_promotion(states[[2]]$manifest, latest_url = urls[1],
    fetch_latest = function(url) { stopifnot(identical(url, urls[1])); states[[1]] },
    artifact_path = coi_hosted_artifact_path)
  proof <- c(proof, coi_hosted_retention_details(states[[1]]$manifest, states[[2]]$manifest))
  proof$first_manifest_url <- urls[1]
  proof$second_manifest_url <- urls[2]
  proof$second_manifest_sha256 <- states[[2]]$manifest_sha256
  dir.create("logs", showWarnings = FALSE)
  destination <- file.path("logs", paste0("coi_hosted_retention_", args[1], "_to_", args[2], ".json"))
  jsonlite::write_json(proof, destination, auto_unbox = TRUE, pretty = TRUE)
  message("Two hosted COI releases retain exact prior sequence/header hashes and scientific source archives; ",
    proof$real_added_records, " real accession versions added, every group grew and acquisition lineage/alignment inputs agree. Proof: ", destination)
}
if (sys.nframe() == 0L) coi_hosted_retention_main()
