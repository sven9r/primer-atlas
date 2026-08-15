#!/usr/bin/env Rscript

# Download the ODbL supplementary inputs used for the expanded BeePrime panel.
# Files are written only after their published-source download passes the
# recorded MD5 checksum.

project_root <- normalizePath(getwd(), mustWork = TRUE)
external_dir <- file.path(project_root, "data", "external", "gurten2026")
dir.create(external_dir, recursive = TRUE, showWarnings = FALSE)

files <- data.frame(
  destination = c(
    "supplement_1_empirical_amplification.docx",
    "supplement_4_penalty_figures.zip",
    "supplement_5_clustered_references.zip",
    "supplement_6_austrian_bees.xls",
    "supplement_7_bees.xlsx",
    "supplement_8_arthropods.xlsx"
  ),
  url = paste0(
    "https://binary.pensoft.net/file/",
    c("1649676", "1649679", "1649680", "1649681", "1649682", "1649683")
  ),
  md5 = c(
    "784ac388f942f62799b03ed60cbf07b5",
    "aa254c8bc9d78e0148b5dc1824fcc40d",
    "ed54684d51df88591f9704db8442e27b",
    "b62838726c483ab8e0f468345116b5d7",
    "d77f99914682c7a36d5e420f49bf5017",
    "b1cc17908299a50f7d872e0990303541"
  ),
  stringsAsFactors = FALSE
)

force_download <- identical(tolower(Sys.getenv("FORCE", "false")), "true")

for (i in seq_len(nrow(files))) {
  destination <- file.path(external_dir, files$destination[i])
  expected_md5 <- files$md5[i]

  if (file.exists(destination)) {
    current_md5 <- unname(tools::md5sum(destination))
    if (identical(current_md5, expected_md5)) {
      message("Verified existing ", basename(destination))
      next
    }
    if (!force_download) {
      stop(
        basename(destination), " exists but its checksum differs. ",
        "Inspect it or rerun with FORCE=true."
      )
    }
  }

  temporary <- tempfile(pattern = "gurten2026_")
  on.exit(unlink(temporary), add = TRUE)
  message("Downloading ", files$url[i])
  download.file(files$url[i], temporary, mode = "wb", quiet = FALSE)
  observed_md5 <- unname(tools::md5sum(temporary))
  if (!identical(observed_md5, expected_md5)) {
    stop(
      "Checksum mismatch for ", basename(destination),
      ": expected ", expected_md5, ", observed ", observed_md5
    )
  }
  if (!file.copy(temporary, destination, overwrite = force_download)) {
    stop("Could not write ", destination)
  }
  unlink(temporary)
}

alignment_path <- file.path(external_dir, "ClusteredReferences.fasta")
if (!file.exists(alignment_path) || force_download) {
  archive <- file.path(external_dir, "supplement_5_clustered_references.zip")
  members <- unzip(archive, list = TRUE)
  fasta_member <- members$Name[grepl("ClusteredReferences\\.fasta$", members$Name)]
  if (length(fasta_member) != 1L) {
    stop("Could not identify one ClusteredReferences.fasta member in Supplement 5.")
  }
  extracted <- unzip(
    archive,
    files = fasta_member,
    exdir = external_dir,
    junkpaths = TRUE,
    overwrite = force_download
  )
  if (!length(extracted) || !file.exists(alignment_path)) {
    stop("Supplement 5 downloaded but the reference FASTA was not extracted.")
  }
}

message(
  "Gurten et al. (2026) supplementary reference inputs are verified in ",
  external_dir
)
