#!/usr/bin/env Rscript
.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
source("R/data_layer.R")
source("R/fungal_versions.R")

index <- fungal_read_index(remote = FALSE)
shipped_index_path <- file.path(fungal_version_root(), "index.json")
shipped_index_hash <- atlas_sha256(shipped_index_path)
shipped_manifest_hashes <- vapply(index$manifest, function(path) atlas_sha256(file.path(fungal_version_root(), path)), character(1))
stopifnot(nrow(index) >= 2L, !anyDuplicated(index$id))
stopifnot(!dir.exists(file.path(fungal_version_root(), ".baseline-inputs")))
v1 <- fungal_read_snapshot(index$id[startsWith(index$id, "v1-")][1], remote = FALSE)
growing_ids <- index$id[!startsWith(index$id, "v1-")]
v2 <- fungal_read_snapshot(tail(growing_ids, 1), remote = FALSE)
stopifnot(length(v1$sequences) == 1L, names(v1$sequences) == "FN812768.2")
stopifnot(length(v2$sequences) > length(v1$sequences), "FN812768.2" %in% names(v2$sequences))
stopifnot(identical(unname(v1$sequences[[1]]), unname(v2$sequences[["FN812768.2"]])))
stopifnot(nrow(v1$raw_catalog) == 121L, nrow(v1$geometry) == 8L, nrow(v2$pairs) >= 20L)
stopifnot(setequal(names(v2$sequences), v2$taxonomy$accession))
stopifnot(length(unique(nchar(v2$alignment))) == 1L)
stopifnot(!identical(v1$manifest$input_sha256, v2$manifest$input_sha256))

# Simulate a clean worker: neither current catalog nor derived inputs exists.
clean_root <- tempfile("fungal-clean-worker-")
dir.create(file.path(clean_root, "data", "pinned", "ITS_FUNGAL"), recursive = TRUE)
stopifnot(file.copy(fungal_version_root(), file.path(clean_root, "data", "pinned", "ITS_FUNGAL"), recursive = TRUE))
clean_v1 <- fungal_read_snapshot(v1$manifest$id, clean_root)
clean_v2 <- fungal_read_snapshot(v2$manifest$id, clean_root)
stopifnot(identical(clean_v1$sequences, v1$sequences), identical(as.data.frame(clean_v2$pairs), as.data.frame(v2$pairs)))
stopifnot(identical(as.data.frame(clean_v1$catalog), as.data.frame(v1$catalog)), identical(as.data.frame(clean_v2$taxonomy), as.data.frame(v2$taxonomy)))

# Corruption must be rejected even when the artifact is already local.
manifest_path <- file.path(fungal_version_root(clean_root), index$manifest[match(v2$manifest$id, index$id)])
sequence_path <- file.path(dirname(manifest_path), v2$manifest$artifacts$sequences$path)
cat("corrupt\n", file = sequence_path, append = TRUE)
failure <- tryCatch({fungal_read_snapshot(v2$manifest$id, clean_root); NULL}, error = identity)
stopifnot(inherits(failure, "error"), grepl("checksum", conditionMessage(failure)))
unlink(clean_root, recursive = TRUE)

# Historical identity must not silently resolve to latest.
stopifnot(grepl("/example-release/manifest.json$", atlas_manifest_path("ITS_FUNGAL", release_id = "example-release")))
failure <- tryCatch({atlas_manifest_path("ITS_FUNGAL", release_id = "../latest"); NULL}, error = identity)
stopifnot(inherits(failure, "error"))

# Exercise the real builder only inside a temporary checkout. Current source
# code can differ from a shipped snapshot, so establish its content identity
# before checking that acquisition timestamps alone reuse that version.
# Repository-only CI still tests frozen loads when build inputs are absent.
panel <- "data/derived/fungal_reference_panel"
if (dir.exists(panel) && file.exists("R/fungal_explorer.R")) {
  # Two independent worker builds: copy the frozen preceding release into a
  # clean checkout, then append a clearly synthetic template and its metadata.
  # Only this temporary fixture receives new versions; shipped history stays unchanged.
  worker <- tempfile("fungal-version-worker-")
  dir.create(worker)
  files <- c("R/functions.R", "R/data_layer.R", "R/fungal_versions.R", "R/fungal_explorer.R",
             "scripts/build_fungal_snapshots.R", "scripts/build_fungal_reference_panel.py",
             "data/catalog/unite_primers.csv")
  for (path in files) {
    target <- file.path(worker, path)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(path, target))
  }
  dir.create(file.path(worker, "data", "pinned", "ITS_FUNGAL"), recursive = TRUE)
  stopifnot(file.copy(fungal_version_root(), file.path(worker, "data", "pinned", "ITS_FUNGAL"), recursive = TRUE))
  worker_panel <- file.path(worker, panel)
  dir.create(worker_panel, recursive = TRUE)
  stopifnot(all(file.copy(list.files(panel, full.names = TRUE), worker_panel)))
  builder <- function() {
    old_directory <- getwd()
    on.exit(setwd(old_directory), add = TRUE)
    setwd(worker)
    output <- system2(file.path(R.home("bin"), "Rscript"), c("scripts/build_fungal_snapshots.R", "v2"),
                      env = c(paste0("R_LIBS=", shQuote(paste(.libPaths(), collapse = .Platform$path.sep))),
                              paste0("FUNGAL_PANEL_DIR=", shQuote(worker_panel))),
                      stdout = TRUE, stderr = TRUE)
    status <- attr(output, "status")
    if (!is.null(status) && status != 0L) stop(paste(output, collapse = "\n"))
  }
  builder()
  before_ids <- fungal_read_index(worker, remote = FALSE)$id
  before_manifests <- vapply(fungal_read_index(worker, remote = FALSE)$manifest,
                            function(path) atlas_sha256(file.path(fungal_version_root(worker), path)), character(1))
  run_path <- file.path(worker_panel, "run_manifest.json")
  run <- jsonlite::read_json(run_path, simplifyVector = FALSE)
  run$acquired_at <- "2099-01-01T00:00:00Z"
  run$previous_count <- run$panel_count
  run$retained_previous_count <- run$panel_count
  run$added_count <- 0L
  run$built_at <- "2099-01-01T00:00:00Z"
  run$searches <- lapply(run$searches, function(search) {
    search$acquired_at <- "2099-01-01T00:00:00Z"
    search$cache_file <- "/different/worker/cache.json"
    search
  })
  jsonlite::write_json(run, run_path, auto_unbox = TRUE, pretty = TRUE)
  builder()
  worker_index <- fungal_read_index(worker, remote = FALSE)
  stopifnot(identical(before_ids, worker_index$id))
  stopifnot(identical(before_manifests, vapply(worker_index$manifest,
    function(path) atlas_sha256(file.path(fungal_version_root(worker), path)), character(1))))
  before <- fungal_read_snapshot(tail(worker_index$id, 1), worker, remote = FALSE)
  previous_manifest_hash <- atlas_sha256(file.path(fungal_version_root(worker), tail(worker_index$manifest, 1)))
  raw <- atlas_read_artifact(list(local_path = file.path(worker_panel, "panel.fasta")))
  aligned <- atlas_read_artifact(list(local_path = file.path(worker_panel, "alignment.fasta")))
  source_accession <- names(raw)[2]
  synthetic_accession <- "SYNTHETIC_RETENTION_FIXTURE.1"
  raw[synthetic_accession] <- raw[[source_accession]]
  aligned[synthetic_accession] <- aligned[[source_accession]]
  write_fasta <- function(value, path) {
    writeLines(unlist(lapply(seq_along(value), function(i) c(paste0(">", names(value)[i]), unname(value[[i]]))), use.names = FALSE), path)
  }
  write_fasta(raw, file.path(worker_panel, "panel.fasta"))
  write_fasta(aligned, file.path(worker_panel, "alignment.fasta"))
  metadata <- read.csv(file.path(worker_panel, "panel_manifest.csv"), stringsAsFactors = FALSE, check.names = FALSE)
  extra <- metadata[match(source_accession, metadata$accession), , drop = FALSE]
  extra$accession <- synthetic_accession
  extra$organism <- "Synthetic retention test fixture"
  extra$source <- "synthetic test fixture"
  extra$source_url <- ""
  metadata <- rbind(metadata, extra)
  write.csv(metadata, file.path(worker_panel, "panel_manifest.csv"), row.names = FALSE, na = "")
  # Extended-provenance parser fixture: text '28S' must remain text while
  # numeric reference indices remain numeric in new general monthly snapshots.
  landmark_fixture <- read.csv(file.path(worker_panel, "landmarks.csv"), stringsAsFactors = FALSE)
  landmark_fixture$parser_fixture_label <- "28S"
  landmark_fixture$parser_fixture_reference_index <- 101L
  write.csv(landmark_fixture, file.path(worker_panel, "landmarks.csv"), row.names = FALSE, na = "")
  provenance_path <- file.path(worker_panel, "run_manifest.json")
  provenance <- jsonlite::read_json(provenance_path, simplifyVector = FALSE)
  provenance$panel_count <- length(raw)
  provenance$previous_count <- length(raw) - 1L
  provenance$retained_previous_count <- length(raw) - 1L
  provenance$added_count <- 1L
  for (filename in names(provenance$artifacts_sha256)) {
    provenance$artifacts_sha256[[filename]] <- atlas_sha256(file.path(worker_panel, filename))
  }
  jsonlite::write_json(provenance, provenance_path, auto_unbox = TRUE, pretty = TRUE)
  builder()
  after_index <- fungal_read_index(worker, remote = FALSE)
  after <- fungal_read_snapshot(tail(after_index$id, 1), worker, remote = FALSE)
  stopifnot(nrow(after_index) == nrow(worker_index) + 1L)
  stopifnot(after$manifest$version_number == before$manifest$version_number + 1L)
  stopifnot(startsWith(after$manifest$id, paste0("v", after$manifest$version_number, "-")))
  stopifnot(identical(after$sequences[names(before$sequences)], before$sequences))
  stopifnot(synthetic_accession %in% after$taxonomy$accession)
  stopifnot(identical(after$landmarks$parser_fixture_label, rep("28S", nrow(after$landmarks))),
    is.numeric(after$landmarks$parser_fixture_reference_index), is.numeric(after$landmarks$start))
  stopifnot(!identical(after$manifest$artifacts$alignment$sha256, before$manifest$artifacts$alignment$sha256))
  stopifnot(identical(previous_manifest_hash, atlas_sha256(file.path(fungal_version_root(worker), tail(worker_index$manifest, 1)))))
  builder()
  stopifnot(identical(fungal_read_index(worker, remote = FALSE)$id, after_index$id))

  # Serve the new version over HTTP: it must appear without a code redeploy,
  # and a clean runner must restore its complete history before appending.
  remote_checks <- function() {
  server_script <- file.path(worker, "serve_fixture.py")
  server_state <- file.path(worker, "server_state.json")
  writeLines(c(
    "import http.server,json,os,pathlib",
    paste0("os.chdir(", jsonlite::toJSON(worker, auto_unbox = TRUE), ")"),
    "server=http.server.ThreadingHTTPServer(('127.0.0.1',0),http.server.SimpleHTTPRequestHandler)",
    paste0("pathlib.Path(", jsonlite::toJSON(server_state, auto_unbox = TRUE), ").write_text(json.dumps({'port':server.server_port,'pid':os.getpid()}))"),
    "server.serve_forever()"
  ), server_script)
  system2("python3", shQuote(server_script), stdout = file.path(worker, "server.log"), stderr = file.path(worker, "server.log"), wait = FALSE)
  for (attempt in seq_len(50L)) {
    if (file.exists(server_state)) break
    Sys.sleep(0.1)
  }
  stopifnot(file.exists(server_state))
  server <- jsonlite::read_json(server_state, simplifyVector = TRUE)
  base <- paste0("http://127.0.0.1:", server$port)
  descriptor <- function(path, rows = NULL) {
    relative <- substring(path, nchar(worker) + 2L)
    value <- list(url = paste0(base, "/", relative), sha256 = atlas_sha256(path), bytes = unname(file.info(path)$size))
    if (!is.null(rows)) value$rows <- rows
    value
  }
  hosted <- list()
  for (i in seq_len(nrow(after_index))) {
    path <- file.path(fungal_version_root(worker), after_index$manifest[i])
    frozen <- jsonlite::read_json(path, simplifyVector = FALSE)
    hosted_artifacts <- lapply(frozen$artifacts, function(artifact) descriptor(file.path(dirname(path), artifact$path), artifact$rows))
    hosted[[frozen$id]] <- list(id = frozen$id, label = frozen$label, description = frozen$description,
                               created_at = frozen$created_at, manifest = descriptor(path), artifacts = hosted_artifacts)
  }
  hosted_latest <- file.path(worker, "releases", "ITS_FUNGAL", "latest.json")
  dir.create(dirname(hosted_latest), recursive = TRUE)
  jsonlite::write_json(list(marker_id = "ITS_FUNGAL", release_version = "fixture-hosted-v3", fungal_snapshots = hosted,
    artifacts = list(fungal_snapshot_index = descriptor(file.path(fungal_version_root(worker), "index.json")))),
    hosted_latest, auto_unbox = TRUE, pretty = TRUE)
  old_remote <- Sys.getenv("ATLAS_MANIFEST_BASE_URL", unset = NA_character_)
  on.exit({
    if (!is.null(server$pid)) tools::pskill(as.integer(server$pid))
    if (is.na(old_remote)) Sys.unsetenv("ATLAS_MANIFEST_BASE_URL") else Sys.setenv(ATLAS_MANIFEST_BASE_URL = old_remote)
  }, add = TRUE)
  Sys.setenv(ATLAS_MANIFEST_BASE_URL = paste0(base, "/releases"))
  remote_index <- fungal_read_index()
  stopifnot(isTRUE(attr(remote_index, "hosted")), !isTRUE(attr(remote_index, "stale")), after$manifest$id %in% remote_index$id)
  remote_v3 <- fungal_read_snapshot(after$manifest$id)
  stopifnot(identical(remote_v3$sequences, after$sequences))
  stopifnot(identical(remote_v3$manifest$transport_source, "verified hosted snapshot"))
  cached_v1 <- atlas_cached_artifact(hosted[[v1$manifest$id]]$artifacts$catalog$url, hosted[[v1$manifest$id]]$artifacts$catalog$sha256)
  cached_v3 <- atlas_cached_artifact(hosted[[after$manifest$id]]$artifacts$catalog$url, hosted[[after$manifest$id]]$artifacts$catalog$sha256)
  stopifnot(!identical(cached_v1, cached_v3))

  fresh <- tempfile("fungal-cloud-history-worker-")
  dir.create(fresh)
  for (path in files) {
    target <- file.path(fresh, path)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(path, target))
  }
  dir.create(file.path(fresh, "data", "pinned", "ITS_FUNGAL"), recursive = TRUE)
  stopifnot(file.copy(fungal_version_root(), file.path(fresh, "data", "pinned", "ITS_FUNGAL"), recursive = TRUE))
  old_directory <- getwd()
  setwd(fresh)
  restored_output <- system2(file.path(R.home("bin"), "Rscript"), c("scripts/build_fungal_snapshots.R", "restore"),
    env = paste0("R_LIBS=", shQuote(paste(.libPaths(), collapse = .Platform$path.sep))), stdout = TRUE, stderr = TRUE)
  setwd(old_directory)
  status <- attr(restored_output, "status")
  if (!is.null(status) && status != 0L) stop(paste(restored_output, collapse = "\n"))
  restored_index <- fungal_read_index(fresh, remote = FALSE)
  stopifnot(after$manifest$id %in% restored_index$id, nrow(restored_index) == nrow(after_index))
  restored <- fungal_read_snapshot(after$manifest$id, fresh, remote = FALSE)
  stopifnot(identical(restored$sequences, after$sequences))
  stopifnot(file.exists(file.path(fresh, "data/derived/fungal_previous_reference_panel/run_manifest.json")))

  # Corrupt both served and cached bytes. The selected cloud-only version must
  # fail rather than borrowing a bundled version or mixing its current catalog.
  selected_artifact <- hosted[[after$manifest$id]]$artifacts$sequences
  cached <- atlas_cached_artifact(selected_artifact$url, selected_artifact$sha256)
  served <- file.path(fungal_version_root(worker), after_index$manifest[nrow(after_index)])
  served <- file.path(dirname(served), after$manifest$artifacts$sequences$path)
  cat("corrupt\n", file = cached, append = TRUE)
  cat("corrupt\n", file = served, append = TRUE)
  failure <- tryCatch({fungal_read_snapshot(after$manifest$id); NULL}, error = identity)
  stopifnot(inherits(failure, "error"), grepl("checksum", conditionMessage(failure)))
  tools::pskill(as.integer(server$pid))
  server$pid <- NULL
  fallback_index <- fungal_read_index()
  stopifnot(!isTRUE(attr(fallback_index, "hosted")), isTRUE(attr(fallback_index, "stale")))
  stopifnot(identical(fallback_index$id, fungal_read_index(remote = FALSE)$id))
  if (is.na(old_remote)) Sys.unsetenv("ATLAS_MANIFEST_BASE_URL") else Sys.setenv(ATLAS_MANIFEST_BASE_URL = old_remote)
  unlink(fresh, recursive = TRUE)
  }
  remote_checks()
  unlink(worker, recursive = TRUE)
}
stopifnot(identical(shipped_index_hash, atlas_sha256(shipped_index_path)))
stopifnot(identical(shipped_manifest_hashes, vapply(index$manifest,
  function(path) atlas_sha256(file.path(fungal_version_root(), path)), character(1))))
message("Frozen history, clean-worker retention, ordinal version growth, replay reuse, and corruption rejection passed.")
