#!/usr/bin/env Rscript

project_root <- normalizePath(getwd(), mustWork = TRUE)
reference_dir <- file.path(project_root, "data", "reference")
dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)

download_ncbi_fasta <- function(accession, output_name) {
  url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
    "db=nuccore&id=", accession, "&rettype=fasta&retmode=text"
  )
  destination <- file.path(reference_dir, output_name)
  temporary <- tempfile(fileext = ".fasta")
  on.exit(unlink(temporary), add = TRUE)
  download.file(url, temporary, mode = "wb", quiet = TRUE)
  lines <- readLines(temporary, warn = FALSE)
  if (!length(lines) || !startsWith(lines[1], ">")) {
    stop("NCBI returned no FASTA for ", accession)
  }
  writeLines(lines, destination, useBytes = TRUE)
  destination
}

path <- download_ncbi_fasta(
  "FN812768.2", "its_fungal_reference_FN812768.2.fasta"
)
message("Wrote fungal ITS pilot reference: ", path)
