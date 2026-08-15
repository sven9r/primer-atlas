suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))

scores <- read_csv(
  "data/derived/beeprime_centroid_scores_gurten2026.csv",
  show_col_types = FALSE
)
hymenoptera <- read_csv(
  "data/derived/beeprime_hymenoptera_centroid_scores.csv",
  show_col_types = FALSE
)
bee_families <- read_csv(
  "data/derived/beeprime_bee_family_summary.csv",
  show_col_types = FALSE
)
bee_subfamilies <- read_csv(
  "data/derived/beeprime_bee_subfamily_summary.csv",
  show_col_types = FALSE
)
bee_genera <- read_csv(
  "data/derived/beeprime_bee_genus_summary.csv",
  show_col_types = FALSE
)
empirical <- read_csv(
  "data/derived/beeprime_empirical_validation.csv",
  show_col_types = FALSE
)
empirical_manifest <- read_csv(
  "data/derived/beeprime_empirical_manifest.csv",
  show_col_types = FALSE
)

stopifnot(
  nrow(scores) == 67352L,
  n_distinct(scores$accession) == 67352L,
  nrow(hymenoptera) == 6797L,
  all(hymenoptera$order == "Hymenoptera"),
  !any(is.na(hymenoptera$accession))
)

bees <- filter(scores, is_bee)
stopifnot(
  nrow(bees) == 590L,
  sum(bees$pair_scorable) == 541L,
  abs(mean(bees$pair_scorable) - 0.9169492) < 1e-6,
  abs(median(bees$pair_penalty, na.rm = TRUE) - 51.69) < 1e-6,
  setequal(
    bee_families$family,
    c(
      "Andrenidae", "Apidae", "Colletidae",
      "Halictidae", "Megachilidae", "Melittidae"
    )
  ),
  sum(bee_families$n_centroids) == 590L,
  sum(bee_subfamilies$n_centroids) == 590L,
  all(!is.na(bee_subfamilies$subfamily)),
  nrow(bee_genera) == 42L,
  sum(bee_genera$n_centroids == 0L) == 12L,
  sum(bee_genera$n_centroids) == 590L
)

apis <- filter(bee_genera, genus == "Apis")
bombus <- filter(bee_genera, genus == "Bombus")
eucera <- filter(bee_genera, genus == "Eucera")
melittidae <- filter(bee_families, family == "Melittidae")
stopifnot(
  nrow(apis) == 1L,
  apis$pair_median_penalty < 10,
  nrow(bombus) == 1L,
  bombus$pair_median_penalty < 10,
  nrow(eucera) == 1L,
  eucera$pair_median_penalty > 100,
  melittidae$n_centroids == 1L
)

target_rows <- filter(empirical, target == "yes")
listed_failures <- filter(target_rows, detected == "no")
manifest_value <- function(field) {
  empirical_manifest$value[match(field, empirical_manifest$field)]
}
stopifnot(
  nrow(target_rows) == 37L,
  sum(target_rows$detected == "yes") == 32L,
  nrow(listed_failures) == 5L,
  manifest_value("article_text_target_bee_extracts") == "38",
  manifest_value("article_text_target_bee_detected") == "33",
  setequal(
    listed_failures$genus,
    c("Ceratina", "Eucera", "Xylocopa")
  )
)

message("Expanded BeePrime reference and empirical-source checks passed.")
