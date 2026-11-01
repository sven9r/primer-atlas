# BeePrime reference identity is independent of the general arthropod panel.
beeprime_reference_choices <- c(
  "Original publication alignment" = "publication",
  "Expanded alignment with new sequences" = "expanded"
)

beeprime_empty_reference <- function(reason, version = NA_character_) {
  list(available = FALSE, mode = "expanded", version = version,
       label = "Expanded alignment with new sequences", reason = reason,
       scores = data.frame(), tables = list(), geography = data.frame(),
       alignment_artifact = NULL, added_sequences = 0L)
}

beeprime_read_reference_artifact <- function(artifact, root = getwd()) {
  if (is.null(artifact) || is.null(artifact$sha256)) stop("Missing reference artifact checksum.")
  path <- artifact$local_path
  if (!is.null(path) && !grepl("^/", path)) path <- file.path(root, path)
  if (is.null(path) || !file.exists(path)) {
    if (is.null(artifact$url)) stop("Reference artifact is unavailable.")
    path <- atlas_cached_artifact(artifact$url, artifact$sha256)
  }
  if (!identical(digest::digest(file = path, algo = "sha256", serialize = FALSE), artifact$sha256)) {
    stop("Reference artifact checksum failed.")
  }
  artifact$local_path <- path
  as.data.frame(atlas_read_artifact(artifact, root), stringsAsFactors = FALSE)
}

beeprime_expanded_reference <- function(descriptor, publication_scores, root = getwd()) {
  if (is.null(descriptor)) return(beeprime_empty_reference(
    "An expanded BeePrime alignment has not been published yet. Choose the original publication alignment to inspect the available evidence."
  ))
  tryCatch({
    stopifnot(
      identical(descriptor$pair_id, "BEEPRIME"),
      identical(descriptor$reference_mode, "expanded"),
      identical(descriptor$taxonomic_scope, "bees"),
      identical(descriptor$publication_baseline, "gurten2026_suppl5"),
      length(descriptor$version) == 1L, nzchar(descriptor$version)
    )
    scores <- beeprime_read_reference_artifact(descriptor$sequence_scores, root)
    required <- c(
      "accession", "pair_id", "order", "family", "subfamily", "genus",
      "is_bee", "cluster_size", "reference_origin", "organism_label",
      "scientific_name_ncbi", "forward_penalty", "reverse_penalty",
      "forward_scorable", "reverse_scorable", "pair_scorable", "pair_penalty",
      "perfect_match_pair", "terminal_3_mismatch_pair",
      "forward_binding_sequence_primer_oriented", "reverse_binding_sequence_primer_oriented",
      "forward_mismatch_count", "reverse_mismatch_count",
      "forward_terminal_3_mismatches", "reverse_terminal_3_mismatches"
    )
    stopifnot(all(required %in% names(scores)), nrow(scores) > 0L,
              !anyDuplicated(scores$accession), all(nzchar(scores$accession)),
              all(scores$pair_id == "BEEPRIME"), all(scores$is_bee %in% TRUE),
              all(scores$order == "Hymenoptera"),
              all(scores$reference_origin %in% c("publication", "added")))
    baseline <- publication_scores$accession[publication_scores$is_bee %in% TRUE]
    stopifnot(setequal(scores$accession[scores$reference_origin == "publication"], baseline))
    added <- scores$accession[scores$reference_origin == "added"]
    stopifnot(length(added) > 0L, !any(added %in% baseline),
              descriptor$added_sequences == length(added),
              all(scores$cluster_size[scores$reference_origin == "added"] == 1))
    scores$pair_label <- "BeePrime"
    scores$lineage_group <- "Bee"
    scores$reference_mode <- "expanded"
    scores$reference_version <- descriptor$version
    summarize <- function(columns) {
      as.data.frame(summarize_expanded_primer_scores(scores, columns))
    }
    tables <- list(
      bee_family = summarize("family"),
      bee_subfamily = summarize(c("family", "subfamily")),
      bee_genus = summarize(c("family", "subfamily", "genus")),
      hymenoptera_family = summarize(c("family", "lineage_group")),
      hymenoptera_lineage = summarize("lineage_group")
    )
    geography <- data.frame(
      accession = scores$accession, country_or_territory = NA_character_,
      locality = NA_character_, normalized_location = NA_character_,
      geographic_resolution = "unresolved", latitude = NA_real_, longitude = NA_real_
    )
    if (!is.null(descriptor$geography)) {
      supplied <- beeprime_read_reference_artifact(descriptor$geography, root)
      stopifnot(all(names(geography) %in% names(supplied)),
                !anyDuplicated(supplied$accession), setequal(supplied$accession, scores$accession))
      geography <- supplied
      geography$country_or_territory[geography$country_or_territory == ""] <- NA_character_
      geography$locality[geography$locality == ""] <- NA_character_
    }
    list(available = TRUE, mode = "expanded", version = descriptor$version,
         label = "Expanded alignment with new sequences", reason = NULL,
         scores = scores, tables = tables, geography = geography,
         alignment_artifact = descriptor$alignment, added_sequences = length(added))
  }, error = function(e) beeprime_empty_reference(
    "The expanded BeePrime reference could not be verified. Choose the original publication alignment to inspect the available evidence.",
    descriptor$version
  ))
}

beeprime_reference_table <- function(reference, name) {
  if (!isTRUE(reference$available)) return(data.frame())
  x <- reference$tables[[name]]
  x$reference_mode <- reference$mode
  x$reference_version <- reference$version
  x
}

beeprime_reference_alignment_path <- function(reference, root = getwd()) {
  if (!isTRUE(reference$available)) stop("The selected reference is unavailable.")
  if (reference$mode == "publication") return(file.path(root, "data", "external", "gurten2026", "ClusteredReferences.fasta"))
  artifact <- reference$alignment_artifact
  if (is.null(artifact) || is.null(artifact$sha256)) stop("The selected alignment download is unavailable.")
  path <- artifact$local_path
  if (!is.null(path) && !grepl("^/", path)) path <- file.path(root, path)
  if (is.null(path) || !file.exists(path)) path <- atlas_cached_artifact(artifact$url, artifact$sha256)
  if (!identical(digest::digest(file = path, algo = "sha256", serialize = FALSE), artifact$sha256)) stop("Reference alignment checksum failed.")
  path
}
