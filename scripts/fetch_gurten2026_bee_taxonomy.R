#!/usr/bin/env Rscript

# Fetch accession-level NCBI taxonomy for the bee centroid labels in Gurten et
# al. (2026) Supplement 5. Raw ESummary JSON and taxonomy XML batches are kept
# so every family/subfamily/tribe assignment can be audited back to NCBI.

project_root <- normalizePath(getwd(), mustWork = TRUE)
external_dir <- file.path(project_root, "data", "external", "gurten2026")
.libPaths(c(file.path(project_root, ".Rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(jsonlite)
  library(readxl)
  library(xml2)
})
source(file.path(project_root, "R", "functions.R"))

alignment_path <- file.path(external_dir, "ClusteredReferences.fasta")
bee_lookup_path <- file.path(external_dir, "supplement_7_bees.xlsx")
if (!file.exists(alignment_path) || !file.exists(bee_lookup_path)) {
  stop("Run scripts/download_gurten2026_reference.R first.")
}

force_fetch <- identical(tolower(Sys.getenv("FORCE", "false")), "true")
ncbi_email <- Sys.getenv("NCBI_EMAIL", "")
ncbi_key <- Sys.getenv("NCBI_API_KEY", "")

query_suffix <- paste0(
  "&tool=coi_primer_atlas",
  if (nzchar(ncbi_email)) paste0("&email=", URLencode(ncbi_email)) else "",
  if (nzchar(ncbi_key)) paste0("&api_key=", URLencode(ncbi_key)) else ""
)

sequences <- parse_fasta(alignment_path)
headers <- names(sequences)
accessions <- sub(" .*", "", headers)
descriptions <- sub("^[^ ]+ +", "", headers)
genera <- gsub("^\\[|\\]$", "", sub(" .*", "", descriptions))
bee_genera <- unique(read_excel(
  bee_lookup_path,
  sheet = "BeesInReference"
)$genus)
bee_accessions <- unique(accessions[genera %in% bee_genera])
if (length(bee_accessions) != 590L) {
  warning("Expected 590 bee centroid accessions; found ", length(bee_accessions))
}

batch_values <- function(values, size) {
  split(values, ceiling(seq_along(values) / size))
}

download_checked <- function(url, destination, validator) {
  if (file.exists(destination) && !force_fetch) {
    validator(destination)
    message("Verified existing ", basename(destination))
    return(invisible(destination))
  }
  temporary <- tempfile(pattern = "ncbi_taxonomy_")
  on.exit(unlink(temporary), add = TRUE)
  download.file(url, temporary, mode = "wb", quiet = FALSE)
  validator(temporary)
  if (!file.copy(temporary, destination, overwrite = force_fetch)) {
    stop("Could not write ", destination)
  }
  unlink(temporary)
  Sys.sleep(if (nzchar(ncbi_key)) 0.12 else 0.4)
  invisible(destination)
}

summary_batches <- batch_values(bee_accessions, 75L)
summary_paths <- character(length(summary_batches))
for (i in seq_along(summary_batches)) {
  ids <- paste(summary_batches[[i]], collapse = ",")
  destination <- file.path(
    external_dir,
    sprintf("bee_centroid_esummary_batch_%02d.json", i)
  )
  url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
    "?db=nuccore&retmode=json&id=", URLencode(ids, reserved = TRUE),
    query_suffix
  )
  validator <- function(path) {
    payload <- fromJSON(path, simplifyVector = FALSE)
    if (is.null(payload$result$uids) || !length(payload$result$uids)) {
      stop("NCBI ESummary returned no records for ", basename(destination))
    }
  }
  download_checked(url, destination, validator)
  summary_paths[i] <- destination
}

taxids <- character()
for (path in summary_paths) {
  payload <- fromJSON(path, simplifyVector = FALSE)
  for (uid in payload$result$uids) {
    taxids <- c(taxids, as.character(payload$result[[uid]]$taxid))
  }
}
taxids <- unique(taxids[nzchar(taxids)])

taxonomy_batches <- batch_values(taxids, 75L)
for (i in seq_along(taxonomy_batches)) {
  ids <- paste(taxonomy_batches[[i]], collapse = ",")
  destination <- file.path(
    external_dir,
    sprintf("bee_taxonomy_batch_%02d.xml", i)
  )
  url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
    "?db=taxonomy&retmode=xml&id=", URLencode(ids, reserved = TRUE),
    query_suffix
  )
  validator <- function(path) {
    document <- read_xml(path)
    if (!length(xml_find_all(document, "/TaxaSet/Taxon"))) {
      stop("NCBI taxonomy returned no records for ", basename(destination))
    }
  }
  download_checked(url, destination, validator)
}

message(
  "Verified NCBI taxonomy for ", length(bee_accessions),
  " bee centroid accessions across ", length(taxids), " taxids."
)
