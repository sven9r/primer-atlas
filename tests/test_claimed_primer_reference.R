suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))

overall <- read_csv(
  "data/derived/claimed_primer_overall_summary.csv",
  show_col_types = FALSE
)
orders <- read_csv(
  "data/derived/claimed_primer_order_summary.csv",
  show_col_types = FALSE
)
groups <- read_csv(
  "data/derived/claimed_primer_phylogeny_group_summary.csv",
  show_col_types = FALSE
)
bindings <- read_csv(
  "data/derived/primer_bindings.csv",
  show_col_types = FALSE
)
order_manifest <- read_csv(
  "data/provenance/order_alignment_manifest.csv",
  show_col_types = FALSE
)
targeted_manifest <- read_csv(
  "data/external/targeted_zbj/targeted_complete_cox1_manifest.csv",
  show_col_types = FALSE
)
targeted_overall <- read_csv(
  "data/derived/targeted_zbj_overall_summary.csv",
  show_col_types = FALSE
)
targeted_families <- read_csv(
  "data/derived/targeted_zbj_family_summary.csv",
  show_col_types = FALSE
)
targeted_groups <- read_csv(
  "data/derived/targeted_zbj_study_group_summary.csv",
  show_col_types = FALSE
)

zbj <- overall |> filter(pair_id == "ZBJ_ART")
zbj_deg <- overall |> filter(pair_id == "ZBJ_ART_DEG")
stopifnot(
  nrow(zbj) == 1L,
  nrow(zbj_deg) == 1L,
  !zbj$reference_suitable,
  !zbj_deg$reference_suitable,
  zbj$reference_suitability == "binding-region incomplete",
  zbj$pair_scorable_fraction < 0.01,
  zbj_deg$pair_scorable_fraction < 0.01,
  grepl("penalties cannot be interpreted", zbj$reference_limitation_reason)
)

zbj_bindings <- bindings |>
  filter(pair_id %in% c("ZBJ_ART", "ZBJ_ART_DEG")) |>
  arrange(direction, pair_id)
stopifnot(
  nrow(zbj_bindings) == 4L,
  setequal(zbj_bindings$alignment_start, c(33L, 220L)),
  all(zbj_bindings$alignment_start[zbj_bindings$direction == "forward"] == 33L),
  all(zbj_bindings$alignment_start[zbj_bindings$direction == "reverse"] == 220L),
  bindings$sequence[
    bindings$pair_id == "ZBJ_ART_DEG" &
      bindings$direction == "reverse"
  ] == "WAYTARTCARTTWCCRAAHCCHCC"
)

spidprey <- overall |> filter(pair_id == "NOSPID")
spider <- orders |> filter(pair_id == "NOSPID", order == "Araneae")
diptera <- orders |> filter(pair_id == "NOSPID", order == "Diptera")
stopifnot(
  spidprey$reference_suitable,
  spidprey$pair_scorable_fraction > 0.95,
  spider$n_centroids > 2500L,
  spider$pair_scorable_fraction > 0.95,
  spider$pair_median_penalty > 500,
  spider$terminal_3_mismatch_fraction_among_scorable > 0.95,
  diptera$pair_median_penalty < 50
)

nospid_groups <- groups |> filter(pair_id == "NOSPID")
stopifnot(
  setequal(
    c(
      "Acari", "Collembola", "Butterflies (Papilionoidea)",
      "Macroheterocera (macro-moth core)",
      "Microlepidoptera (operational grade)",
      "Boundary / convention-sensitive moths"
    ),
    nospid_groups$display_group
  ),
  sum(
    nospid_groups$n_centroids[
      nospid_groups$display_group %in% c(
        "Butterflies (Papilionoidea)",
        "Macroheterocera (macro-moth core)",
        "Microlepidoptera (operational grade)",
        "Boundary / convention-sensitive moths"
      )
    ]
  ) == sum(orders$n_centroids[
    orders$pair_id == "NOSPID" & orders$order == "Lepidoptera"
  ])
)

new_full_length_groups <- order_manifest |>
  filter(order %in% c("Acari", "Collembola"))
stopifnot(
  nrow(new_full_length_groups) == 2L,
  all(new_full_length_groups$status == "complete"),
  all(new_full_length_groups$downloaded_sequences >= 100L),
  all(new_full_length_groups$aligned_sequences_including_reference > 40L)
)

stopifnot(
  nrow(targeted_manifest) == 5L,
  all(targeted_manifest$complete_cox1_features_extracted >= 60L),
  all(targeted_manifest$minimum_cox1_length >= 1400L),
  all(targeted_manifest$maximum_cox1_length <= 1700L),
  all(targeted_overall$pair_scorable_fraction > 0.90)
)

targeted_hymenoptera <- targeted_overall |>
  filter(pair_id == "ZBJ_ART", source_group == "Hymenoptera")
targeted_lepidoptera <- targeted_overall |>
  filter(pair_id == "ZBJ_ART", source_group == "Lepidoptera")
targeted_coleoptera <- targeted_overall |>
  filter(pair_id == "ZBJ_ART", source_group == "Coleoptera")
targeted_hymenoptera_deg <- targeted_overall |>
  filter(pair_id == "ZBJ_ART_DEG", source_group == "Hymenoptera")
apis <- targeted_families |>
  filter(
    pair_id == "ZBJ_ART",
    source_group == "Hymenoptera",
    family == "Apidae"
  )
stopifnot(
  targeted_hymenoptera$n_sequences > 150L,
  targeted_hymenoptera$n_species > 100L,
  targeted_hymenoptera$pair_median_penalty > 500,
  targeted_lepidoptera$pair_median_penalty < 50,
  targeted_coleoptera$n_sequences > 100L,
  targeted_coleoptera$n_species > 75L,
  targeted_coleoptera$pair_scorable_fraction > 0.90,
  targeted_hymenoptera_deg$pair_median_penalty <
    targeted_hymenoptera$pair_median_penalty,
  apis$n_sequences > 50L,
  apis$pair_median_penalty > 800
)

lepidoptera_groups <- targeted_groups |>
  filter(pair_id == "ZBJ_ART", source_group == "Lepidoptera")
stopifnot(
  setequal(
    lepidoptera_groups$study_group,
    c(
      "Butterflies (Papilionoidea)",
      "Macroheterocera (macro-moth core)",
      "Microlepidoptera (operational grade)",
      "Boundary / convention-sensitive moths"
    )
  ),
  sum(lepidoptera_groups$n_sequences) ==
    targeted_lepidoptera$n_sequences
)

message("Target-claim reference, ZBJ orientation, and phylogeny-group checks passed.")
