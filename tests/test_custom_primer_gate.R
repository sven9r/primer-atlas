suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))

source("R/functions.R")

reference <- clean_sequence(
  parse_fasta("data/reference/coi_reference_NC_001322.1.fasta")[[1]]
)
primer_library <- read_csv("data/primers.csv", show_col_types = FALSE)
alignment_paths <- list.files(
  "data/primerminer/aligned_reference",
  pattern = "_COX1_reference_aligned\\.fasta$",
  full.names = TRUE
)
alignment_orders <- sub(
  "_COX1_reference_aligned\\.fasta$",
  "",
  basename(alignment_paths)
)
alignment_templates <- setNames(
  lapply(alignment_paths, function(path) {
    sequences <- parse_fasta(path)
    sequences[!grepl("^NC_001322\\.1_COX1_reference", names(sequences))]
  }),
  alignment_orders
)

stopifnot(clean_sequence("acgt-n") == "ACGT-N")

for (pair in unique(primer_library$pair_id)) {
  rows <- filter(primer_library, pair_id == pair)
  result <- locate_primer_pair(
    reference_sequence = reference,
    forward_sequence = rows$sequence[rows$direction == "forward"][1],
    reverse_sequence = rows$sequence[rows$direction == "reverse"][1],
    expected_amplicon_bp = rows$reported_amplicon_bp[1],
    alignment_templates = alignment_templates
  )
  stopifnot(
    result$forward_assessment$credible ||
      result$best$evidence_basis_forward == "installed COI alignments",
    result$reverse_assessment$credible ||
      result$best$evidence_basis_reverse == "installed COI alignments"
  )
}

beeprime <- locate_primer_pair(
  reference_sequence = reference,
  forward_sequence = "ATGAATTAATAATGATCANATYTATAAYWC",
  reverse_sequence = "GGATAWACWGTTCAWCCWGTWCC",
  alignment_templates = alignment_templates
)
stopifnot(
  beeprime$alignment_rescue_used,
  beeprime$best$start_forward == 132L,
  beeprime$best$end_forward == 161L,
  beeprime$best$start_reverse == 361L,
  beeprime$best$end_reverse == 383L,
  beeprime$best$targeted_region_bp == 199L,
  beeprime$best$aligned_amplicon_bp == 252L,
  any(
    beeprime$pair_order_support$order == "Hymenoptera" &
      beeprime$pair_order_support$n_pair_supported >= 2L
  )
)

mollcoi253 <- locate_primer_pair(
  reference_sequence = reference,
  forward_sequence = "GGAGTAGGAACTGGTTGGAC",
  reverse_sequence = "CAGCTGCTAACACAGGCA",
  expected_amplicon_bp = 253,
  alignment_templates = alignment_templates
)
stopifnot(
  mollcoi253$alignment_rescue_used,
  mollcoi253$best$start_forward == 355L,
  mollcoi253$best$start_reverse == 590L,
  mollcoi253$best$aligned_amplicon_bp == 253L,
  mollcoi253$alignment_minimum_identity == 0.70,
  any(
    mollcoi253$pair_order_support$order == "Bivalvia" &
      mollcoi253$pair_order_support$n_pair_supported >= 2L
  )
)

unrelated_28s <- tryCatch(
  locate_primer_pair(
    reference_sequence = reference,
    forward_sequence = "ACAAGTACCGTGAGGGAAAGTTG",
    reverse_sequence = "TCGGAAGGAACCAGCTACTA",
    alignment_templates = alignment_templates
  ),
  error = function(error) error
)
stopifnot(
  inherits(unrelated_28s, "error"),
  grepl("No credible COI binding pair", conditionMessage(unrelated_28s)),
  grepl("chance-match expectation", conditionMessage(unrelated_28s))
)

message("Custom-primer binding gate tests passed.")
