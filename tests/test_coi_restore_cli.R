#!/usr/bin/env Rscript
.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
source("R/data_layer.R")
source("R/coi_retention.R")
repo <- normalizePath(getwd())
fixture <- tempfile("coi_restore_cli_"); dir.create(fixture)
on.exit(unlink(fixture, recursive = TRUE), add = TRUE)
groups <- c("Araneae", "Diptera")
worker <- function(name, manifest = NULL) {
  root <- file.path(fixture, name); dir.create(root)
  for (file in c("R/data_layer.R", "R/coi_retention.R", "scripts/restore_coi_reference_state.R")) {
    target <- file.path(root, file); dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(file.path(repo, file), target))
  }
  if (dir.exists(file.path(repo, ".Rlib"))) stopifnot(file.symlink(file.path(repo, ".Rlib"), file.path(root, ".Rlib")))
  dir.create(file.path(root, "data"), showWarnings = FALSE)
  write.csv(data.frame(order = groups), file.path(root, "data/orders_22.csv"), row.names = FALSE)
  if (!is.null(manifest)) {
    path <- file.path(root, "data/releases/COI/fixture-first/manifest.json")
    dir.create(dirname(path), recursive = TRUE)
    jsonlite::write_json(manifest, path, pretty = TRUE, auto_unbox = TRUE)
  }
  root
}
write_fasta <- function(values, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(unlist(lapply(seq_along(values), function(i) c(paste0(">", names(values)[i]), values[i])), use.names = FALSE), path)
}
original <- worker("original")
for (i in seq_along(groups)) write_fasta(setNames(c("ACGTACGT", "ACGTACGA"),
  paste0("CLI", i, c("A.1 frozen header", "B.1 frozen header"))),
  file.path(original, "data/primerminer/raw", groups[i], paste0(groups[i], "_GB.fasta")))
records <- coi_retained_records(original, groups)
dir.create(file.path(original, "data/provenance"), recursive = TRUE)
write.csv(records, file.path(original, "data/provenance/coi_reference_accession_ledger.csv"), row.names = FALSE)
release_dir <- file.path(fixture, "immutable"); dir.create(release_dir)
state <- coi_stage_retained_state(original, release_dir, "fixture-first", "https://data.example.invalid", groups)
manifest <- list(marker_id = "COI", release_version = "fixture-first", coi_retained_state = state)
run_cli <- function(root, remote = "", bootstrap = "false") {
  previous <- getwd(); on.exit(setwd(previous), add = TRUE); setwd(root)
  output <- file.path(root, "restore_cli.log")
  env <- c(paste0("ATLAS_MANIFEST_BASE_URL=", shQuote(remote)),
    "COI_PREVIOUS_RELEASE_ID=fixture-first", paste0("COI_ALLOW_BOOTSTRAP=", bootstrap))
  status <- system2(file.path(R.home("bin"), "Rscript"), "scripts/restore_coi_reference_state.R",
    env = env, stdout = output, stderr = output)
  list(status = status, output = paste(readLines(output, warn = FALSE), collapse = "\n"))
}
restored <- worker("full-state", manifest)
result <- run_cli(restored)
stopifnot(result$status == 0L, grepl("Restored 4 retained COI sequences from fixture-first", result$output))
coi_assert_retention(records, coi_validate_retained_state(restored, groups)$records)
for (artifact in state$artifacts) stopifnot(identical(atlas_sha256(file.path(restored, artifact$restore_path)), artifact$sha256))
lineage <- jsonlite::read_json(file.path(restored, "data/provenance/coi_reference_restore.json"))
stopifnot(lineage$previous_release_id == "fixture-first", lineage$retained_records == 4L)
expect_failure <- function(name, manifest, remote = "", pattern) {
  root <- worker(name, manifest)
  result <- run_cli(root, remote)
  stopifnot(result$status != 0L, grepl(pattern, result$output),
    !dir.exists(file.path(root, "data/primerminer/raw")))
}
expect_failure("stale-fallback", manifest, "https://127.0.0.1:1", "stale fallback is unsafe")
expect_failure("missing-manifest", NULL, pattern = "No release manifest is available")
expect_failure("missing-state", list(marker_id = "COI", release_version = "fixture-first"), pattern = "explicit bootstrap is required")
ledger_only <- manifest; ledger_only$coi_retained_state$artifacts <- state$artifacts[c("records", "ledger")]
expect_failure("ledger-only", ledger_only, pattern = "raw COI FASTA is missing")
# The authorized legacy branch must parse and execute too, with its narrower
# membership-only provenance claim and no fabricated raw sequence state.
legacy_path <- file.path(fixture, "legacy_membership.csv")
write.csv(data.frame(order = records$target_group, template = records$header), legacy_path, row.names = FALSE)
legacy <- list(marker_id = "COI", release_version = "fixture-first", artifacts = list(full_order_anml = list(
  local_path = legacy_path, url = "https://data.example.invalid/legacy.csv", sha256 = atlas_sha256(legacy_path))))
bootstrap <- worker("explicit-bootstrap", legacy)
result <- run_cli(bootstrap, bootstrap = "true")
baseline <- read.csv(file.path(bootstrap, "data/provenance/coi_published_baseline.csv"), stringsAsFactors = FALSE)
stopifnot(result$status == 0L, grepl("Explicit COI bootstrap: recovered 4", result$output),
  !dir.exists(file.path(bootstrap, "data/primerminer/raw")), nrow(baseline) == 4L,
  all(!baseline$original_raw_bytes_available), all(baseline$published_artifact_sha256 == atlas_sha256(legacy_path)))
message("Actual COI restore CLI executes full-state/bootstrap branches and rejects stale, missing and ledger-only inputs before installation.")
