# Choose a complete scoring panel and verify that it reflects the retained
# sequence inputs. A partial full-panel directory must never hide missing groups.
coi_scoring_inputs <- function(root = getwd(), mode = "auto", groups = NULL,
                               minimum_per_group = 100L, minimum_total = 10000L) {
  if (!mode %in% c("auto", "full", "centroid_preview")) stop("Unknown COI scoring mode")
  if (is.null(groups)) {
    orders <- read.csv(file.path(root, "data/orders_22.csv"), stringsAsFactors = FALSE)
    groups <- orders$order[orders$order != "Bivalvia"]
  }
  full_dir <- file.path(root, "data/primerminer/full_aligned_reference")
  full_files <- list.files(full_dir, pattern = "_COX1_full_reference_aligned[.]fasta$", full.names = TRUE)
  if (mode == "auto") mode <- if (length(full_files)) "full" else "centroid_preview"
  if (mode == "centroid_preview") {
    files <- list.files(file.path(root, "data/primerminer/aligned_reference"),
                        pattern = "_COX1_reference_aligned[.]fasta$", full.names = TRUE)
    if (!length(files)) stop("No reference-coordinate preview alignments are installed")
    return(list(files = files, evaluation_mode = mode, evidence_class = "legacy_centroid_preview",
                output_prefix = "preview_order_", manifest = NULL))
  }
  files <- file.path(full_dir, paste0(groups, "_COX1_full_reference_aligned.fasta"))
  if (any(!file.exists(files))) stop("The full COI scoring panel is incomplete: ", paste(groups[!file.exists(files)], collapse = ", "))
  coi_validate_retained_state(root, groups)
  reference_path <- file.path(root, "data/reference/coi_reference_NC_001322.1.fasta")
  blocks <- lapply(seq_along(groups), function(i) {
    group <- groups[i]
    raw_path <- file.path(root, "data/primerminer/raw", group, paste0(group, "_GB.fasta"))
    sidecar <- paste0(files[i], ".inputs.json")
    if (!file.exists(sidecar)) stop("Full COI alignment lacks verified input provenance: ", group)
    inputs <- jsonlite::read_json(sidecar, simplifyVector = TRUE)
    if (!identical(inputs$evaluation_mode, "full")) stop("A preview alignment was supplied as full COI evidence")
    checked <- coi_validate_full_alignment(files[i], raw_path, reference_path, inputs, minimum_per_group)
    data.frame(order = group, evaluation_mode = "full", evidence_class = "full_order_reference_alignment",
               n_templates = checked$aligned_templates, raw_sha256 = atlas_sha256(raw_path),
               reference_sha256 = atlas_sha256(reference_path), alignment_sha256 = atlas_sha256(files[i]),
               alignment_input_sha256 = inputs$input_sha256, stringsAsFactors = FALSE)
  })
  manifest <- do.call(rbind, blocks)
  if (sum(manifest$n_templates) < minimum_total) stop("The full COI panel is below its total sequence target")
  list(files = files, evaluation_mode = "full", evidence_class = "full_order_reference_alignment",
       output_prefix = "full_order_", manifest = manifest)
}
