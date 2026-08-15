#!/usr/bin/env Rscript

# Retrieve complete mitochondrial records, extract their annotated COX1 coding
# feature, and create one site-complete FASTA per priority clade. Sampling is
# reproducible across the complete NCBI search result rather than taking only
# the newest records. Raw EFetch batches and an explicit manifest are retained.

project_root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(project_root, ".Rlib"), .libPaths()))
suppressPackageStartupMessages(library(rentrez))
source(file.path(project_root, "R", "functions.R"))

external_dir <- file.path(project_root, "data", "external", "targeted_zbj")
raw_dir <- file.path(external_dir, "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

groups <- c(
  "Hymenoptera", "Lepidoptera", "Acari", "Collembola", "Coleoptera"
)
n_requested <- as.integer(Sys.getenv("N_PER_GROUP", "300"))
if (is.na(n_requested) || n_requested < 1L) {
  stop("N_PER_GROUP must be a positive integer.")
}
force <- identical(tolower(Sys.getenv("FORCE", "false")), "true")
set.seed(20260728)

split_batches <- function(values, size = 20L) {
  split(values, ceiling(seq_along(values) / size))
}

write_fasta <- function(sequences, path) {
  lines <- unlist(lapply(seq_along(sequences), function(i) {
    sequence <- sequences[[i]]
    wrapped <- substring(
      sequence,
      seq(1L, nchar(sequence), by = 80L),
      pmin(seq(80L, nchar(sequence) + 79L, by = 80L), nchar(sequence))
    )
    c(paste0(">", names(sequences)[i]), wrapped)
  }), use.names = FALSE)
  writeLines(lines, path, useBytes = TRUE)
}

manifest_rows <- list()
for (group_name in groups) {
  message("Searching complete mitochondrial records for ", group_name)
  query <- paste0(
    group_name,
    "[Organism] AND mitochondrion[filter] AND ",
    "(\"complete genome\"[Title] OR \"mitochondrial genome\"[Title])"
  )
  search <- entrez_search(
    db = "nuccore",
    term = query,
    retmax = 10000L
  )
  if (!length(search$ids)) {
    stop("No complete mitochondrial records found for ", group_name)
  }
  sampled_ids <- sample(
    search$ids,
    min(n_requested, length(search$ids)),
    replace = FALSE
  )
  batches <- split_batches(sampled_ids)
  extracted <- list()
  extracted_index <- 1L

  for (i in seq_along(batches)) {
    batch_path <- file.path(
      external_dir,
      sprintf("%s_cds_batch_%02d.fasta", tolower(group_name), i)
    )
    if (!file.exists(batch_path) || force) {
      fetched <- entrez_fetch(
        db = "nuccore",
        id = batches[[i]],
        rettype = "fasta_cds_na",
        retmode = "text"
      )
      writeLines(fetched, batch_path, useBytes = TRUE)
      Sys.sleep(0.4)
    }
    cds <- parse_fasta(batch_path)
    cox1 <- cds[
      grepl(
        "\\[gene=(COX1|COI|CO1|COXI)\\]",
        names(cds),
        ignore.case = TRUE
      ) |
        grepl(
          "\\[protein=cytochrome (c )?oxidase subunit (I|1)\\]",
          names(cds),
          ignore.case = TRUE
        )
    ]
    if (!length(cox1)) next
    for (j in seq_along(cox1)) {
      header <- names(cox1)[j]
      accession <- sub(
        "^.*lcl\\|([^_ ]+)_cds_.*$",
        "\\1",
        header
      )
      if (!grepl("^[A-Z]{1,4}[0-9]+\\.[0-9]+$", accession)) next
      sequence <- clean_sequence(cox1[[j]])
      if (nchar(sequence) < 1400L || nchar(sequence) > 1700L) next
      extracted[[extracted_index]] <- sequence
      names(extracted)[extracted_index] <- paste(
        accession,
        group_name,
        "complete mitochondrial COX1 CDS"
      )
      extracted_index <- extracted_index + 1L
    }
  }

  extracted <- unlist(extracted, use.names = TRUE)
  accessions <- sub(" .*", "", names(extracted))
  extracted <- extracted[!duplicated(accessions)]
  output_path <- file.path(raw_dir, paste0(group_name, "_complete_COX1.fasta"))
  write_fasta(extracted, output_path)

  manifest_rows[[group_name]] <- data.frame(
    source_group = group_name,
    query = query,
    ncbi_records_found = as.integer(search$count),
    mitochondrial_records_sampled = length(sampled_ids),
    complete_cox1_features_extracted = length(extracted),
    minimum_cox1_length = if (length(extracted)) min(nchar(extracted)) else NA,
    maximum_cox1_length = if (length(extracted)) max(nchar(extracted)) else NA,
    output_fasta = sub(paste0("^", project_root, "/"), "", output_path),
    sampling_seed = 20260728L,
    retrieved_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
  message(
    group_name, ": extracted ", length(extracted),
    " complete COX1 features from ", length(sampled_ids), " sampled genomes."
  )
}

manifest <- do.call(rbind, manifest_rows)
rownames(manifest) <- NULL
write.csv(
  manifest,
  file.path(external_dir, "targeted_complete_cox1_manifest.csv"),
  row.names = FALSE,
  na = ""
)
message("Wrote targeted complete-COX1 manifest.")
