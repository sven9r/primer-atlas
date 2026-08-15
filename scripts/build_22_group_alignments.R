#!/usr/bin/env Rscript

project_root <- normalizePath(getwd(), mustWork = TRUE)
local_library <- file.path(project_root, ".Rlib")
.libPaths(c(local_library, .libPaths()))

suppressPackageStartupMessages({
  library(PrimerMiner)
  library(rentrez)
})

source(file.path(project_root, "R", "functions.R"))

orders_path <- file.path(project_root, "data", "orders_22.csv")
orders <- read.csv(orders_path, stringsAsFactors = FALSE, check.names = FALSE)
requested_raw <- Sys.getenv("ORDERS", "")
if (nzchar(requested_raw)) {
  requested <- trimws(strsplit(requested_raw, ",", fixed = TRUE)[[1]])
  missing_orders <- setdiff(requested, orders$order)
  if (length(missing_orders)) {
    stop("Unknown ORDERS value(s): ", paste(missing_orders, collapse = ", "))
  }
  orders <- orders[match(requested, orders$order), , drop = FALSE]
}

subset_size <- as.integer(Sys.getenv("GB_SUBSET", "50"))
if (is.na(subset_size) || subset_size < 1) stop("GB_SUBSET must be a positive integer.")
force <- identical(tolower(Sys.getenv("FORCE", "false")), "true")
set.seed(20260728)

vsearch <- Sys.which("vsearch")
mafft <- Sys.which("mafft")
if (!nzchar(vsearch)) stop("VSEARCH was not found on PATH.")
if (!nzchar(mafft)) stop("MAFFT was not found on PATH.")

reference_dir <- file.path(project_root, "data", "reference")
raw_root <- file.path(project_root, "data", "primerminer", "raw")
cluster_root <- file.path(project_root, "data", "primerminer", "clustered")
alignment_root <- file.path(project_root, "data", "primerminer", "aligned_reference")
provenance_dir <- file.path(project_root, "data", "provenance")
dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(raw_root, recursive = TRUE, showWarnings = FALSE)
dir.create(cluster_root, recursive = TRUE, showWarnings = FALSE)
dir.create(alignment_root, recursive = TRUE, showWarnings = FALSE)
dir.create(provenance_dir, recursive = TRUE, showWarnings = FALSE)

reference_accession <- "NC_001322.1"
reference_start <- 1474L
reference_stop <- 3009L
reference_path <- file.path(reference_dir, "coi_reference_NC_001322.1.fasta")

write_fasta_text <- function(text, path) {
  text <- sub("^>[^\\r\\n]*", paste0(
    ">NC_001322.1_COX1_reference coordinates=", reference_start, "-", reference_stop
  ), text)
  writeLines(text, path, useBytes = TRUE)
}

if (!file.exists(reference_path) || force) {
  message("Fetching versioned COX1 reference ", reference_accession, ":", reference_start, "-", reference_stop)
  reference_fasta <- rentrez::entrez_fetch(
    db = "nuccore",
    id = reference_accession,
    rettype = "fasta",
    seq_start = reference_start,
    seq_stop = reference_stop,
    strand = 1
  )
  write_fasta_text(reference_fasta, reference_path)
}
reference_sequences <- parse_fasta(reference_path)
if (length(reference_sequences) != 1L || nchar(reference_sequences[[1]]) != 1536L) {
  stop("Reference FASTA must contain exactly one 1,536-bp sequence.")
}

count_fasta <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  sum(startsWith(readLines(path, warn = FALSE), ">"))
}

run_mafft <- function(reference, fragments, output, log_path) {
  args <- c(
    "--quiet",
    "--adjustdirectionaccurately",
    "--keeplength",
    "--addfragments", shQuote(fragments),
    shQuote(reference)
  )
  status <- system2(mafft, args, stdout = output, stderr = log_path)
  if (!identical(status, 0L)) stop("MAFFT failed; see ", log_path)

  aligned <- parse_fasta(output)
  widths <- nchar(aligned)
  if (!length(aligned) || any(widths != 1536L)) {
    stop("Reference-coordinate alignment is not exactly 1,536 columns: ", output)
  }
}

run_primerminer_clustering <- function(raw_path, destination_dir) {
  scratch <- tempfile(pattern = "primerminer_cluster_")
  dir.create(scratch, recursive = TRUE)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  scratch_input <- file.path(scratch, basename(raw_path))
  file.copy(raw_path, scratch_input, overwrite = TRUE)

  PrimerMiner::Clustering(
    file = basename(scratch_input),
    vsearchpath = vsearch,
    id = 0.97,
    threshold = "Majority",
    setwd = scratch
  )

  generated <- file.path(
    scratch,
    paste0(sub("\\.fasta$", "", basename(raw_path)), "_cons_cluster_Majority.fasta")
  )
  if (!file.exists(generated)) {
    stop("PrimerMiner clustering did not create ", basename(generated))
  }
  file.copy(generated, destination_dir, overwrite = TRUE)
  if (file.exists(file.path(scratch, "log.txt"))) {
    file.copy(file.path(scratch, "log.txt"), destination_dir, overwrite = TRUE)
  }
  if (dir.exists(file.path(scratch, "Vsearch"))) {
    unlink(file.path(destination_dir, "Vsearch"), recursive = TRUE)
    dir.create(file.path(destination_dir, "Vsearch"), showWarnings = FALSE)
    file.copy(
      list.files(file.path(scratch, "Vsearch"), full.names = TRUE),
      file.path(destination_dir, "Vsearch"),
      recursive = TRUE,
      overwrite = TRUE
    )
  }
}

manifest_rows <- vector("list", nrow(orders))
query_mode <- tolower(Sys.getenv("QUERY_MODE", "broad"))
query_suffix <- switch(
  query_mode,
  broad = paste0(
    "[Organism] AND (COI[Gene] OR CO1[Gene] OR COX1[Gene] OR COXI[Gene])",
    " AND 1200:1800[Sequence Length]"
  ),
  complete_cds = paste0(
    "[Organism] AND (COI[Gene] OR CO1[Gene] OR COX1[Gene] OR COXI[Gene])",
    " AND 1400:1800[Sequence Length] AND complete cds[Title]",
    " NOT COII[Title]"
  ),
  folmer_barcode = paste0(
    "[Organism] AND (COI[Gene] OR CO1[Gene] OR COX1[Gene] OR COXI[Gene])",
    " AND 600:800[Sequence Length]"
  ),
  stop("QUERY_MODE must be 'broad', 'complete_cds', or 'folmer_barcode'.")
)

for (i in seq_len(nrow(orders))) {
  order_name <- orders$order[i]
  message("\n=== ", order_name, " ===")
  raw_dir <- file.path(raw_root, order_name)
  cluster_dir <- file.path(cluster_root, order_name)
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)

  raw_path <- file.path(raw_dir, paste0(order_name, "_GB.fasta"))
  cluster_path <- file.path(
    cluster_dir, paste0(order_name, "_GB_cons_cluster_Majority.fasta")
  )
  aligned_path <- file.path(alignment_root, paste0(order_name, "_COX1_reference_aligned.fasta"))
  mafft_log <- file.path(alignment_root, paste0(order_name, "_mafft.log"))

  status <- "complete"
  note <- ""
  tryCatch({
    if (!file.exists(raw_path) || force) {
      unlink(file.path(raw_dir, c("log.txt", paste0(order_name, "_GB.fasta"))))
      PrimerMiner::Download_GB(
        taxon = order_name,
        folder = NULL,
        custom_query = query_suffix,
        GB_subset = subset_size,
        setwd = raw_dir
      )
      downloaded_at_root <- file.path(project_root, paste0(order_name, "_GB.fasta"))
      if (file.exists(downloaded_at_root) && !file.exists(raw_path)) {
        file.rename(downloaded_at_root, raw_path)
      }
    }
    if (!file.exists(raw_path) || count_fasta(raw_path) < 1L) {
      stop("No GenBank sequences were downloaded.")
    }

    if (!file.exists(cluster_path) || force) {
      unlink(file.path(cluster_dir, "Vsearch"), recursive = TRUE)
      unlink(file.path(cluster_dir, "log.txt"))
      run_primerminer_clustering(raw_path, cluster_dir)
      generated <- file.path(
        cluster_dir,
        paste0(sub("\\.fasta$", "", basename(raw_path)), "_cons_cluster_Majority.fasta")
      )
      if (!identical(generated, cluster_path) && file.exists(generated)) {
        file.rename(generated, cluster_path)
      }
    }
    if (!file.exists(cluster_path) || count_fasta(cluster_path) < 1L) {
      stop("PrimerMiner clustering did not produce consensus sequences.")
    }

    if (!file.exists(aligned_path) || count_fasta(aligned_path) < 1L || force) {
      run_mafft(reference_path, cluster_path, aligned_path, mafft_log)
    }
  }, error = function(e) {
    status <<- "failed"
    note <<- conditionMessage(e)
    message("FAILED: ", note)
  })

  manifest_rows[[i]] <- data.frame(
    order = order_name,
    environment_group = orders$environment_group[i],
    query = paste0(order_name, query_suffix),
    requested_subset = subset_size,
    downloaded_sequences = count_fasta(raw_path),
    consensus_clusters = count_fasta(cluster_path),
    aligned_sequences_including_reference = count_fasta(aligned_path),
    reference_accession = reference_accession,
    reference_feature_coordinates = paste0(reference_start, "-", reference_stop),
    reference_length_bp = 1536L,
    clustering_identity = 0.97,
    status = status,
    note = note,
    raw_fasta = sub(paste0("^", project_root, "/"), "", raw_path),
    cluster_consensus_fasta = sub(paste0("^", project_root, "/"), "", cluster_path),
    reference_aligned_fasta = sub(paste0("^", project_root, "/"), "", aligned_path),
    built_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
}

manifest_path <- file.path(provenance_dir, "order_alignment_manifest.csv")
new_manifest <- do.call(rbind, manifest_rows)
if (file.exists(manifest_path) && nzchar(requested_raw)) {
  old_manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
  old_manifest <- old_manifest[!old_manifest$order %in% new_manifest$order, , drop = FALSE]
  new_manifest <- rbind(old_manifest, new_manifest)
}
new_manifest <- new_manifest[order(new_manifest$order), , drop = FALSE]
write.csv(new_manifest, manifest_path, row.names = FALSE, na = "")
message("\nWrote ", manifest_path)
print(new_manifest[, c(
  "order", "downloaded_sequences", "consensus_clusters",
  "aligned_sequences_including_reference", "status", "note"
)], row.names = FALSE)
