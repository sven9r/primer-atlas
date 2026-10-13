#!/usr/bin/env Rscript
source("scripts/verify_coi_hosted_retention.R")
groups <- read.csv("data/orders_22.csv", stringsAsFactors = FALSE)$order
groups <- groups[groups != "Bivalvia"]
hash <- function(x) digest::digest(x, algo = "sha256", serialize = FALSE)
make <- function(id, n, previous_n, previous_hash = NULL) {
  artifacts <- list(growth_qa = list(value = data.frame(target_group = groups,
    previous_retained = previous_n, retained_after = n, verified_marker_records = n)),
    panel_run_manifest = list(sha256 = hash(paste0(id, " manifest")), value = list(marker_id = "COI",
      release_id = id, previous_run_manifest_sha256 = previous_hash, retained_count = n * length(groups),
      verified_marker_count = n * length(groups))))
  for (group in groups) {
    raw_hash <- hash(paste(id, group, "raw"))
    output_hash <- hash(paste(id, group, "alignment"))
    artifacts[[paste0("raw_", group)]] <- list(sha256 = raw_hash)
    artifacts[[paste0("alignment_", group)]] <- list(sha256 = output_hash)
    artifacts[[paste0("alignment_inputs_", group)]] <- list(value = list(raw_sha256 = raw_hash,
      output_sha256 = output_hash, input_sha256 = hash(paste(id, group, "inputs")), aligned_templates = n))
  }
  list(marker_id = "COI", release_version = id, coi_retained_state = list(target_groups = groups,
    group_counts = setNames(as.list(rep(n, length(groups))), groups), policy = list(monthly_growth = 0.01), artifacts = artifacts))
}
first <- make("fixture-first", 500L, 490L)
second <- make("fixture-second", 505L, 500L, first$coi_retained_state$artifacts$panel_run_manifest$sha256)
read_fixture <- function(artifact) artifact$value
verify <- function(next_manifest = second, prior_manifest = first) coi_hosted_retention_details(prior_manifest, next_manifest, read_fixture)
expect_error <- function(value, pattern) {
  result <- try(force(value), silent = TRUE)
  stopifnot(inherits(result, "try-error"), grepl(pattern, as.character(result)))
}
proof <- verify()
stopifnot(proof$real_added_records == 110L, proof$acquisition_lineage_verified,
  proof$second_verified_marker_records == 11110L, proof$second_minimum_verified_per_group == 505L)
bad <- second; bad$coi_retained_state$artifacts$panel_run_manifest$value$previous_run_manifest_sha256 <- hash("wrong predecessor")
expect_error(verify(bad), "exact first frozen acquisition manifest")
bad <- second; bad$coi_retained_state$artifacts$panel_run_manifest$value$previous_run_manifest_sha256 <- NULL
expect_error(verify(bad), "exact first frozen acquisition manifest")
bad <- second; bad$coi_retained_state$artifacts$growth_qa$value$previous_retained[1] <- 499L
expect_error(verify(bad), "reconcile to actual retained counts")
bad <- second; bad$coi_retained_state$artifacts$growth_qa$value$verified_marker_records[1] <- 506L
expect_error(verify(bad), "verified marker floors")
bad <- second; bad$coi_retained_state$artifacts$growth_qa$value$verified_marker_records[1] <- NA_integer_
expect_error(verify(bad), "verified marker floors")
bad <- second; bad$coi_retained_state$artifacts$panel_run_manifest$value$verified_marker_count <- 11109L
expect_error(verify(bad), "acquisition manifests do not reconcile")
bad <- second; bad$coi_retained_state$artifacts[[paste0("alignment_inputs_", groups[1])]]$value$input_sha256 <-
  first$coi_retained_state$artifacts[[paste0("alignment_inputs_", groups[1])]]$value$input_sha256
expect_error(verify(bad), "alignment inputs/counts")
bad <- second; bad$coi_retained_state$group_counts[[groups[1]]] <- 500L
expect_error(verify(bad), "required real records")
bad <- first; bad$coi_retained_state$artifacts$growth_qa$value$verified_marker_records[1] <- NA_integer_
expect_error(verify(prior_manifest = bad), "First hosted verified QA")
# The hosted proof must ignore existing candidate local files as well as old
# release local files. The mock represents the anonymous download boundary.
fixture <- tempfile("coi_hosted_proof_"); dir.create(fixture)
on.exit(unlink(fixture, recursive = TRUE), add = TRUE)
local <- file.path(fixture, "existing_local.csv"); remote <- file.path(fixture, "downloaded_public.csv")
write.csv(data.frame(source = "local"), local, row.names = FALSE)
write.csv(data.frame(source = "public"), remote, row.names = FALSE)
original_cache <- atlas_cached_artifact
calls <- 0L
atlas_cached_artifact <- function(url, sha256 = NULL, ...) {
  stopifnot(identical(url, "https://data.example.invalid/public.csv"), identical(sha256, atlas_sha256(remote)))
  calls <<- calls + 1L
  remote
}
artifact <- list(url = "https://data.example.invalid/public.csv", local_path = local, sha256 = atlas_sha256(remote))
stopifnot(identical(coi_hosted_artifact_path(artifact, getwd(), FALSE), remote),
  identical(coi_hosted_artifact_path(artifact, getwd(), TRUE), remote),
  coi_hosted_read_artifact(artifact)$source == "public", calls == 3L)
atlas_cached_artifact <- original_cache
message("Hosted COI proof requires anonymous artifacts for both releases, exact acquisition lineage, reconciled verified QA and changed alignment inputs.")
