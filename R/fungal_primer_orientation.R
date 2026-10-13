# Explicit opposite-role oligos derived from a frozen catalog component.
# These are computational candidates, not corrections to the cited source or
# published assays. Reversing both role and 5'-to-3' oligo preserves the raw
# exact-search query under the repository's established nucleotide model.

.fungal_orientation_sequence <- function(sequence) {
  sequence <- as.character(sequence)
  if (length(sequence) != 1L || is.na(sequence) || !nzchar(sequence) ||
      sequence != toupper(gsub("[[:space:]]", "", sequence)) ||
      !grepl("^[ACGTURYSWKMBDHVNI]+$", sequence)) {
    stop("Orientation requires a normalized ungapped nucleotide oligo.")
  }
  sequence
}

.fungal_orientation_rc <- function(sequence) {
  sequence <- .fungal_orientation_sequence(sequence)
  # Same symbols and semantics as R/functions.R::reverse_complement():
  # U is treated as T; inosine I is treated as unrestricted nucleotide N.
  complements <- c(A="T", T="A", U="A", C="G", G="C", R="Y", Y="R",
    S="S", W="W", K="M", M="K", B="V", V="B", D="H", H="D", N="N", I="N")
  paste(rev(unname(complements[strsplit(sequence, "", fixed=TRUE)[[1L]]])), collapse="")
}

.fungal_orientation_base <- function(snapshot) {
  if (!is.list(snapshot) || !is.list(snapshot$manifest) ||
      length(snapshot$manifest$id) != 1L || is.na(snapshot$manifest$id) ||
      !is.character(snapshot$manifest$id) || !nzchar(snapshot$manifest$id)) {
    stop("A selected immutable source snapshot identity is required.")
  }
  catalog <- if (is.data.frame(snapshot$frozen_component_catalog)) {
    snapshot$frozen_component_catalog
  } else snapshot$catalog
  if (!is.data.frame(catalog)) stop("A frozen component catalog is required.")
  catalog <- as.data.frame(catalog, stringsAsFactors=FALSE)
  required <- c("catalog_id", "primer_name", "sequence", "direction", "gene_locus")
  if (!all(required %in% names(catalog)) || anyNA(catalog$catalog_id) ||
      any(!nzchar(as.character(catalog$catalog_id))) || anyDuplicated(catalog$catalog_id) ||
      anyNA(catalog$primer_name) || any(!nzchar(as.character(catalog$primer_name)))) {
    stop("Frozen orientation-parent identities are incomplete or repeated.")
  }
  if (("orientation_derived" %in% names(catalog) && any(catalog$orientation_derived %in% TRUE)) ||
      ("primary_reference_status" %in% names(catalog) &&
       any(catalog$primary_reference_status %in% "derived_reverse_complement"))) {
    stop("Use the original frozen catalog as the orientation parent, not another derived variant.")
  }
  if (anyNA(catalog$direction) || any(!catalog$direction %in% c("forward", "reverse"))) {
    stop("Each frozen orientation parent requires an effective forward or reverse role.")
  }
  for (sequence in catalog$sequence) .fungal_orientation_sequence(sequence)
  if (anyDuplicated(paste(catalog$primer_name, catalog$sequence, catalog$direction, sep="\036"))) {
    stop("Orientation parents require canonical name, sequence and role identities.")
  }
  rownames(catalog) <- NULL
  catalog
}

.fungal_orientation_field <- function(catalog, field) {
  if (!field %in% names(catalog)) return(rep("", nrow(catalog)))
  value <- as.character(catalog[[field]])
  value[is.na(value)] <- ""
  value
}

# Returns ONLY explicit supplemental variants. Callers keep snapshot$catalog
# unchanged and supply these rows to the combination builder separately.
fungal_orientation_catalog <- function(snapshot) {
  original <- .fungal_orientation_base(snapshot)
  variants <- original
  if (!nrow(original)) {
    fields <- c("parent_catalog_id", "parent_primer_name", "parent_sequence", "parent_direction",
      "parent_source_direction", "parent_primary_reference_status", "parent_primary_reference_key",
      "parent_primary_reference_url", "parent_reference", "parent_source_url", "parent_reported_position",
      "parent_source_version", "primary_reference_status", "source_direction", "source_version",
      "reported_position", "orientation_origin", "derivation_rule", "parent_noncanonical_symbols",
      "sequence_derivation_note", "direction_note", "reference", "remarks")
    for (field in fields) variants[[field]] <- character()
    variants$orientation_derived <- variants$raw_query_equivalence <- logical()
    return(variants)
  }
  variants$parent_catalog_id <- as.character(original$catalog_id)
  variants$parent_primer_name <- as.character(original$primer_name)
  variants$parent_sequence <- as.character(original$sequence)
  variants$parent_direction <- as.character(original$direction)
  variants$parent_source_direction <- .fungal_orientation_field(original, "source_direction")
  variants$parent_primary_reference_status <- .fungal_orientation_field(original, "primary_reference_status")
  variants$parent_primary_reference_key <- .fungal_orientation_field(original, "primary_reference_key")
  variants$parent_primary_reference_url <- .fungal_orientation_field(original, "primary_reference_url")
  variants$parent_reference <- .fungal_orientation_field(original, "reference")
  variants$parent_source_url <- .fungal_orientation_field(original, "source_url")
  variants$parent_reported_position <- .fungal_orientation_field(original, "reported_position")
  variants$parent_source_version <- rep(snapshot$manifest$id, nrow(original))
  variants$direction <- ifelse(original$direction == "forward", "reverse", "forward")
  variants$sequence <- vapply(original$sequence, .fungal_orientation_rc, character(1), USE.NAMES=FALSE)
  variants$primer_name <- paste0(original$primer_name, " [RC]")
  identity <- vapply(seq_len(nrow(original)), function(i) paste(
    "reverse_complement_opposite_role_v1", original$catalog_id[i], original$primer_name[i],
    original$sequence[i], original$direction[i], variants$direction[i], variants$sequence[i], sep="\036"), character(1))
  variants$catalog_id <- paste0("FUNGAL_RC_", substr(vapply(identity, digest::digest,
    character(1), algo="sha256", serialize=FALSE, USE.NAMES=FALSE), 1L, 24L))
  if (anyDuplicated(variants$catalog_id) || any(variants$catalog_id %in% original$catalog_id) ||
      any(paste(variants$primer_name, variants$sequence, variants$direction, sep="\036") %in%
          paste(original$primer_name, original$sequence, original$direction, sep="\036"))) {
    stop("A derived orientation identity conflicts with its frozen source catalog.")
  }
  variants$primary_reference_status <- rep("derived_reverse_complement", nrow(original))
  variants$source_direction <- rep("derived", nrow(original))
  variants$source_version <- rep(snapshot$manifest$id, nrow(original))
  # A source's reported position describes the original component; raw sites
  # are reused only after a separate exact parent-and-query identity check.
  variants$reported_position <- rep(NA_character_, nrow(original))
  variants$orientation_derived <- rep(TRUE, nrow(original))
  variants$orientation_origin <- rep("derived_reverse_complement", nrow(original))
  variants$derivation_rule <- rep("reverse_complement_opposite_role_v1", nrow(original))
  variants$raw_query_equivalence <- rep(TRUE, nrow(original))
  variants$parent_noncanonical_symbols <- vapply(original$sequence, function(sequence) {
    paste(c("I", "U")[vapply(c("I", "U"), grepl, logical(1), x=sequence, fixed=TRUE)], collapse=",")
  }, character(1), USE.NAMES=FALSE)
  variants$sequence_derivation_note <- vapply(seq_len(nrow(original)), function(i) {
    note <- "Computational 5′→3′ reverse complement with opposite effective role; original citation describes the parent oligo only. No publication, synthesis or PCR validation is asserted for this derived candidate."
    if (grepl("I", original$sequence[i], fixed=TRUE)) note <- paste(note,
      "The repository nucleotide model normalizes inosine I to degenerate N; N is not chemically identical to inosine. Exact-search query equivalence is nucleotide-set equivalence, not demonstrated biochemical equivalence.")
    if (grepl("U", original$sequence[i], fixed=TRUE)) note <- paste(note,
      "The repository nucleotide model treats U as T; this normalization is disclosed rather than a synthesis recommendation.")
    note
  }, character(1))
  variants$direction_note <- paste0("Derived opposite role from ", original$primer_name, " (",
    original$direction, ") in ", snapshot$manifest$id, ". ", variants$sequence_derivation_note)
  variants$reference <- paste0("Derived reverse complement of ", original$primer_name, " (",
    original$direction, "); parent source: ", ifelse(nzchar(variants$parent_reference),
      variants$parent_reference, "citation unavailable"), ". Derived candidate is not a published assay.")
  prior_remarks <- .fungal_orientation_field(original, "remarks")
  variants$remarks <- paste0(variants$sequence_derivation_note,
    ifelse(nzchar(prior_remarks), paste0(" Parent catalog remarks: ", prior_remarks), ""))
  rownames(variants) <- NULL
  variants
}

# Resolve either a frozen ID or its canonical derived ID to the requested
# effective role. The original-role result is the original row, unchanged.
fungal_orientation_resolve <- function(snapshot, parent_id, direction) {
  if (length(parent_id) != 1L || is.na(parent_id) || !is.character(parent_id) || !nzchar(parent_id) ||
      length(direction) != 1L || is.na(direction) || !direction %in% c("forward", "reverse")) {
    stop("Choose one canonical parent identity and a forward or reverse role.")
  }
  original <- .fungal_orientation_base(snapshot)
  variants <- fungal_orientation_catalog(snapshot)
  index <- match(parent_id, original$catalog_id)
  if (is.na(index)) index <- match(parent_id, variants$catalog_id)
  if (is.na(index)) stop("The orientation identity does not belong to this selected frozen catalog.")
  if (original$direction[index] == direction) original[index, , drop=FALSE] else variants[index, , drop=FALSE]
}

# Shared guard for algebraic site-cache reuse. An opposite-direction relabel
# with the unchanged oligo cannot pass. Cross-collection reuse may retain the
# originating source version, while binding evidence belongs to the receiving
# collection; its parent component still has to match exactly.
fungal_orientation_validate <- function(snapshot, rows, require_source_version=TRUE) {
  if (!is.data.frame(rows) || length(require_source_version) != 1L ||
      is.na(require_source_version) || !is.logical(require_source_version)) {
    stop("Validate explicit orientation rows with a scalar source-version policy.")
  }
  if (!nrow(rows)) return(invisible(TRUE))
  expected <- fungal_orientation_catalog(snapshot)
  identity_fields <- c("catalog_id", "primer_name", "sequence", "direction", "parent_catalog_id",
    "parent_primer_name", "parent_sequence", "parent_direction", "parent_source_direction",
    "parent_primary_reference_status", "parent_primary_reference_key", "parent_primary_reference_url",
    "parent_reference", "parent_source_url", "orientation_origin", "derivation_rule",
    "primary_reference_status", "parent_noncanonical_symbols", "sequence_derivation_note")
  required <- c(identity_fields, "orientation_derived", "raw_query_equivalence", "source_version", "parent_source_version")
  if (!all(required %in% names(rows)) || anyNA(rows$catalog_id) || anyDuplicated(rows$catalog_id)) {
    stop("Derived orientation identity and provenance fields are incomplete or repeated.")
  }
  indices <- match(rows$catalog_id, expected$catalog_id)
  if (anyNA(indices)) stop("Derived orientation identity is not canonical for this frozen parent catalog.")
  for (field in identity_fields) {
    actual <- as.character(rows[[field]])
    wanted <- as.character(expected[[field]][indices])
    if (anyNA(actual) || !identical(actual, wanted)) {
      stop("Derived orientation field differs from its exact frozen parent: ", field)
    }
  }
  for (field in c("orientation_derived", "raw_query_equivalence")) {
    values <- toupper(as.character(rows[[field]]))
    if (anyNA(values) || any(values != "TRUE")) stop("Derived orientation proof flag is not true: ", field)
  }
  source <- as.character(rows$source_version)
  parent_source <- as.character(rows$parent_source_version)
  if (anyNA(source) || anyNA(parent_source) || any(!nzchar(source)) ||
      !identical(source, parent_source) ||
      (require_source_version && any(source != snapshot$manifest$id))) {
    stop("Derived orientation source provenance does not match the required originating snapshot.")
  }
  invisible(TRUE)
}
