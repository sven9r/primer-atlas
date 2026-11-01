#!/usr/bin/env Rscript
.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
source("R/data_layer.R")
source("R/coi_retention.R")
source("R/coi_promotion.R")
fixture <- tempfile("coi_promotion_"); dir.create(fixture)
on.exit(unlink(fixture, recursive = TRUE), add = TRUE)
groups <- c("Araneae", "Diptera")
public <- "https://data.example.invalid"
latest_url <- paste0(public, "/releases/COI/latest.json")
source_bytes <- "immutable marker response"
taxonomy_bytes <- "immutable taxonomy response"
hash_text <- function(x) digest::digest(x, algo = "sha256", serialize = FALSE)
path_map <- new.env(parent = emptyenv())
write_fasta <- function(values, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(unlist(lapply(seq_along(values), function(i) c(paste0(">", names(values)[i]), values[i])), use.names = FALSE), path)
}
build <- function(id, values, baseline = NULL) {
  root <- file.path(fixture, id); dir.create(root)
  dir.create(file.path(root, "data/provenance"), recursive = TRUE)
  for (group in groups) write_fasta(values[[group]], file.path(root, "data/primerminer/raw", group, paste0(group, "_GB.fasta")))
  rows <- coi_retained_records(root, groups)
  ledger <- rows
  ledger$source_sha256 <- hash_text(source_bytes)
  ledger$taxonomy_source_sha256 <- hash_text(taxonomy_bytes)
  ledger$first_release_id <- "original-biological-version"
  ledger$added_at_utc <- "2026-10-01T00:00:00Z"
  source_dir <- file.path(root, "data/derived/coi_reference_panel/raw/content")
  dir.create(source_dir, recursive = TRUE)
  writeBin(charToRaw(source_bytes), file.path(source_dir, paste0(hash_text(source_bytes), ".xml")))
  writeBin(charToRaw(taxonomy_bytes), file.path(source_dir, paste0(hash_text(taxonomy_bytes), ".json")))
  write.csv(ledger, file.path(root, "data/provenance/coi_reference_accession_ledger.csv"), row.names = FALSE)
  if (!is.null(baseline)) write.csv(baseline, file.path(root, "data/provenance/coi_published_baseline.csv"), row.names = FALSE)
  release_dir <- file.path(root, "release"); dir.create(release_dir)
  state <- coi_stage_retained_state(root, release_dir, id, public, groups)
  for (artifact in state$artifacts) assign(artifact$url, artifact$local_path, path_map)
  list(marker_id = "COI", release_version = id, coi_retained_state = state,
       previous_release_id = "an-older-selected-base")
}
original <- list(Araneae = c("A.1 frozen provenance" = "ACGTACGT", "B.1 frozen provenance" = "ACGTACGA"),
                 Diptera = c("D.1 frozen provenance" = "ACGTACGC"))
expanded <- original
expanded$Araneae <- c(expanded$Araneae, "C.1 new provenance" = "ACGTACGG")
prior <- build("hosted-newer", original)
candidate <- build("valid-superset", expanded)
# Hosted reads ignore any supplied local_path. The fixture supplies a verified
# public-object lookup instead of a network server; all bytes still get hashed.
artifact_path <- function(artifact, root, hosted) {
  if (hosted) get(artifact$url, path_map, inherits = FALSE) else artifact$local_path
}
fetch <- function(manifest, stale = FALSE, source = latest_url) function(url) {
  stopifnot(identical(url, latest_url))
  list(manifest = manifest, stale = stale, source = source, manifest_sha256 = hash_text("current hosted manifest"))
}
verify <- function(next_manifest = candidate, prior_manifest = prior, fetch_latest = fetch(prior_manifest)) {
  coi_verify_promotion(next_manifest, latest_url = latest_url, fetch_latest = fetch_latest, artifact_path = artifact_path)
}
expect_error <- function(value, pattern) {
  result <- try(force(value), silent = TRUE)
  stopifnot(inherits(result, "try-error"), grepl(pattern, as.character(result)))
}
proof <- verify()
stopifnot(proof$retained_from_current == 3L, proof$candidate_retained_records == 4L,
          proof$preserved_source_archives == 2L, proof$current_release_id == "hosted-newer")
# A worker based on an older release must still retain a newer hosted addition.
missing <- expanded; missing$Araneae <- missing$Araneae[-2]
expect_error(verify(build("missing-hosted-accession", missing)), "accession was dropped")
# Internally coherent ledger/raw changes cannot hide replacement of old bytes.
changed <- expanded; changed$Araneae[1] <- "TCGTACGT"
expect_error(verify(build("changed-hosted-sequence", changed)), "sequence_sha256")
changed <- expanded; names(changed$Araneae)[1] <- "A.1 edited provenance"
expect_error(verify(build("changed-hosted-header", changed)), "header_sha256")
expect_error(verify(fetch_latest = function(url) stop("hosted latest unavailable")), "unavailable")
expect_error(verify(fetch_latest = fetch(prior, stale = TRUE)), "stale or local")
expect_error(verify(fetch_latest = fetch(prior, source = "/local/latest.json")), "stale or local")
expect_error(coi_promotion_public_url("https://user:secret@example.invalid/latest.json"), "public HTTPS")
expect_error(coi_promotion_public_url("https://example.invalid/latest.json?token=value"), "public HTTPS")
missing_ledger <- candidate; missing_ledger$coi_retained_state$artifacts$ledger <- NULL
expect_error(verify(missing_ledger), "ledger and actual raw")
missing_archive <- candidate
source_name <- names(missing_archive$coi_retained_state$artifacts)[vapply(missing_archive$coi_retained_state$artifacts,
  function(x) identical(x$sha256, hash_text(source_bytes)), logical(1))]
missing_archive$coi_retained_state$artifacts[source_name] <- NULL
expect_error(verify(missing_archive), "source archive was dropped")
corrupt <- build("corrupt-candidate", expanded)
cat("changed after packaging", file = corrupt$coi_retained_state$artifacts$ledger$local_path, append = TRUE)
expect_error(verify(corrupt), "artifact checksum failed")
# Initial activation from the published legacy release proves membership and
# publication identity, explicitly without an unavailable original raw hash.
legacy_path <- file.path(fixture, "legacy_membership.csv")
legacy_rows <- data.frame(order = c("Araneae", "Araneae", "Diptera", "Bivalvia"),
  template = c("A.1 original published header", "B.1 published header", "D.1 published header", "Z.1 excluded archive"))
write.csv(legacy_rows, legacy_path, row.names = FALSE)
legacy_artifact <- list(url = paste0(public, "/releases/COI/legacy/full_order_anml.csv"),
  sha256 = atlas_sha256(legacy_path), local_path = "/deliberately-unavailable-local-path")
assign(legacy_artifact$url, legacy_path, path_map)
legacy <- list(marker_id = "COI", release_version = "legacy", artifacts = list(full_order_anml = legacy_artifact))
baseline <- data.frame(target_group = c("Araneae", "Araneae", "Diptera"), accession = c("A.1", "B.1", "D.1"),
  published_release_id = "legacy", published_artifact_sha256 = legacy_artifact$sha256,
  published_artifact_url = legacy_artifact$url, original_raw_bytes_available = FALSE)
legacy_candidate <- build("legacy-activation", expanded, baseline)
proof <- verify(legacy_candidate, legacy)
stopifnot(proof$retained_from_current == 3L, grepl("original raw bytes unavailable", proof$evidence))
expect_error(verify(build("legacy-missing-member", missing, baseline), legacy), "legacy COI accession was dropped")
bad_baseline <- baseline; bad_baseline$published_artifact_sha256[1] <- hash_text("different publication")
expect_error(verify(build("wrong-legacy-source", expanded, bad_baseline), legacy), "provenance differs")
message("COI promotion compares current hosted membership, raw/header hashes and immutable source archives; stale/local latest and destructive supersets fail closed.")
