#!/usr/bin/env Rscript

.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
suppressPackageStartupMessages(library(rsconnect))

app_files <- c(
  "app.R", "R", "www", "vendor/PrimerMiner", "data/catalog", "data/pinned",
  "data/releases/COI/latest.json", "data/releases/ITS_FUNGAL/latest.json",
  "data/reference/coi_reference_NC_001322.1.fasta",
  "data/reference/its_fungal_reference_FN812768.2.fasta",
  "data/primers.csv", "data/citations.csv", "data/primer_target_claims.csv",
  "data/derived/pair_geometry.csv", "data/derived/marker_pair_geometry.csv",
  "data/derived/primer_bindings.csv", "data/derived/folmer_region.csv",
  "data/derived/order_pair_summary.csv",
  "data/derived/beeprime_reference_manifest.csv",
  "data/derived/claimed_primer_overall_summary.csv",
  "data/derived/claimed_primer_order_summary.csv",
  "data/derived/claimed_primer_phylogeny_group_summary.csv",
  "data/derived/beeprime_bee_family_summary.csv",
  "data/derived/beeprime_bee_subfamily_summary.csv",
  "data/derived/beeprime_bee_genus_summary.csv",
  "data/derived/beeprime_hymenoptera_family_summary.csv",
  "data/derived/beeprime_hymenoptera_lineage_summary.csv",
  "data/derived/beeprime_hymenoptera_centroid_scores.csv",
  "data/derived/beeprime_empirical_validation.csv",
  "data/derived/beeprime_empirical_genus_summary.csv",
  "data/derived/targeted_zbj_overall_summary.csv",
  "data/derived/targeted_zbj_family_summary.csv",
  "data/derived/targeted_zbj_genus_summary.csv",
  "data/derived/targeted_zbj_subfamily_summary.csv",
  "data/derived/targeted_zbj_study_group_summary.csv",
  "data/derived/targeted_zbj_sequence_scores.csv.gz",
  "data/external/targeted_zbj/targeted_complete_cox1_manifest.csv",
  "data/primerminer/aligned_reference"
)
missing <- app_files[!file.exists(app_files)]
if (length(missing)) stop("Missing Connect file(s): ", paste(missing, collapse = ", "))
writeManifest(appDir = ".", appFiles = app_files, python = NULL, quarto = FALSE)
manifest <- jsonlite::read_json("manifest.json", simplifyVector = FALSE)
paths <- names(manifest$files)
bytes <- sum(file.info(paths)$size, na.rm = TRUE)
if (bytes >= 100 * 1024^2) stop("Connect bundle is too large: ", round(bytes / 1024^2, 1), " MB")
message("Connect manifest contains ", length(paths), " files / ", round(bytes / 1024^2, 1), " MB.")
