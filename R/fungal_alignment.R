# Project exact raw-template sites onto a frozen multiple-sequence alignment.
# Alignment-column distances include gaps; physical products retain the engine's
# raw inclusive amplicon_bp. Missing or inconsistent evidence is never placed.

.fungal_alignment_accession <- function(name) {
  sub("^_R_", "", sub("[[:space:]].*$", "", name))
}

fungal_alignment_index <- function(alignment, sequences = NULL) {
  if (!is.character(alignment) || (length(alignment) &&
      (is.null(names(alignment)) || anyNA(names(alignment)) || any(!nzchar(names(alignment)))))) {
    stop("Frozen alignment must be a named character vector.")
  }
  if (anyNA(alignment)) stop("Frozen alignment contains a missing sequence.")
  alignment <- toupper(gsub("[[:space:]]", "", alignment))
  if (!length(alignment)) names(alignment) <- character()
  accessions <- .fungal_alignment_accession(names(alignment))
  if (anyDuplicated(accessions)) stop("Frozen alignment repeats an accession identity.")
  widths <- nchar(alignment)
  if (length(unique(widths)) > 1L) stop("Frozen alignment rows have different column counts.")
  n_columns <- if (length(widths)) as.integer(widths[1]) else 0L
  metadata <- data.frame(
    accession = accessions, alignment_name = names(alignment),
    orientation = ifelse(startsWith(names(alignment), "_R_"), "reverse", "forward"),
    strand = ifelse(startsWith(names(alignment), "_R_"), -1L, 1L),
    n_template_bases = integer(length(alignment)), n_gap_columns = integer(length(alignment)),
    n_ambiguous_bases = integer(length(alignment)),
    status = rep("raw_template_unavailable", length(alignment)), stringsAsFactors = FALSE
  )
  if (!is.null(sequences)) {
    if (!is.character(sequences) || anyNA(sequences) || (length(sequences) && is.null(names(sequences))) ||
        anyNA(names(sequences)) || any(!nzchar(names(sequences))) ||
        anyDuplicated(.fungal_alignment_accession(names(sequences)))) {
      stop("Raw templates must be a named character vector with unique accessions.")
    }
    names(sequences) <- .fungal_alignment_accession(names(sequences))
    sequences <- toupper(gsub("[-[:space:]]", "", sequences))
  }
  maps <- setNames(vector("list", length(alignment)), accessions)
  occupancy_blocks <- list()
  for (i in seq_along(alignment)) {
    bases <- strsplit(alignment[i], "", fixed = TRUE)[[1]]
    present <- bases != "-"
    columns <- which(present)
    template <- paste0(bases[present], collapse = "")
    metadata$n_template_bases[i] <- length(columns)
    metadata$n_gap_columns[i] <- sum(!present)
    metadata$n_ambiguous_bases[i] <- sum(!bases[present] %in% c("A", "C", "G", "T"))
    if (!grepl("^[ACGTURYSWKMBDHVNI-]*$", alignment[i])) {
      metadata$status[i] <- "invalid_alignment_bases"
      maps[i] <- list(integer())
      next
    }
    if (!is.null(sequences) && accessions[i] %in% names(sequences)) {
      raw <- unname(sequences[[accessions[i]]])
      expected <- if (metadata$strand[i] == -1L) {
        if (!exists("reverse_complement", mode = "function")) stop("Source R/functions.R for reverse-strand alignment validation.")
        reverse_complement(raw)
      } else raw
      metadata$status[i] <- if (identical(template, expected)) "verified" else "alignment_template_mismatch"
    }
    maps[[i]] <- if (metadata$strand[i] == -1L) rev(columns) else columns
    names(maps[[i]]) <- as.character(seq_along(maps[[i]]))
    runs <- rle(present)
    ends <- cumsum(runs$lengths)
    if (any(runs$values)) {
      occupancy_blocks[[length(occupancy_blocks) + 1L]] <- data.frame(
        accession = accessions[i], start_col = ends[runs$values] - runs$lengths[runs$values] + 1L,
        end_col = ends[runs$values], n_bases = runs$lengths[runs$values], stringsAsFactors = FALSE
      )
    }
  }
  occupancy <- if (length(occupancy_blocks)) do.call(rbind, occupancy_blocks) else {
    data.frame(accession = character(), start_col = integer(), end_col = integer(), n_bases = integer())
  }
  names(alignment) <- accessions
  structure(list(alignment = alignment, maps = maps, n_columns = n_columns,
                 metadata = metadata, occupancy = occupancy), class = "fungal_alignment_index")
}

# Frozen landmark intervals are raw-template coordinates, just like primer
# sites. An annotation belongs only to its own accession, including legacy
# unscoped pilot landmarks whose manifest identifies a single reference.
fungal_alignment_landmark_rows <- function(index, landmarks, reference_accession = NULL) {
  if (!inherits(index, "fungal_alignment_index")) stop("index must be produced by fungal_alignment_index().")
  if (is.null(landmarks)) landmarks <- data.frame(landmark_id = character(), label = character(),
                                                start = integer(), end = integer())
  required <- c("landmark_id", "label", "start", "end")
  if (!all(required %in% names(landmarks))) stop("Landmarks are missing raw coordinates or labels.")
  rows <- as.data.frame(landmarks, stringsAsFactors = FALSE)
  if (!"accession" %in% names(rows)) {
    reference <- if (length(reference_accession) == 1L && !is.na(reference_accession)) as.character(reference_accession) else NA_character_
    rows$accession <- rep(reference, nrow(rows))
  }
  if (!"source_note" %in% names(rows)) rows$source_note <- rep("Frozen snapshot landmark", nrow(rows))
  rows$evidence_status <- ifelse(grepl("curated", rows$source_note, ignore.case = TRUE),
                                "curated_reference", "annotated_reference")
  for (column in c("start_col", "end_col", "left_col", "right_col", "direction", "raw_bp")) {
    rows[[column]] <- rep(NA_integer_, nrow(rows))
  }
  rows$alignment_status <- rep("template_not_in_alignment", nrow(rows))
  for (i in seq_len(nrow(rows))) {
    accession <- as.character(rows$accession[i])
    match_row <- match(accession, index$metadata$accession)
    if (is.na(match_row)) next
    status <- index$metadata$status[match_row]
    if (!identical(status, "verified")) {
      rows$alignment_status[i] <- status
      next
    }
    bounds <- suppressWarnings(as.numeric(c(rows$start[i], rows$end[i])))
    map <- index$maps[[accession]]
    if (length(bounds) != 2L || any(!is.finite(bounds)) || any(bounds != floor(bounds)) ||
        bounds[1] < 1L || bounds[2] < bounds[1] || bounds[2] > length(map)) {
      rows$alignment_status[i] <- "invalid_raw_coordinates"
      next
    }
    columns <- unname(map[as.integer(bounds)])
    rows$start_col[i] <- columns[1]
    rows$end_col[i] <- columns[2]
    rows$left_col[i] <- min(columns)
    rows$right_col[i] <- max(columns)
    rows$direction[i] <- index$metadata$strand[match_row]
    rows$raw_bp[i] <- as.integer(diff(bounds) + 1L)
    rows$alignment_status[i] <- "mapped"
  }
  rows
}

fungal_alignment_pair_rows <- function(alignment = NULL, pair_sequences,
                                       pair_ids = NULL, accessions = NULL,
                                       sequences = NULL, index = NULL) {
  if (is.null(index)) {
    if (inherits(alignment, "fungal_alignment_index")) index <- alignment else {
      index <- fungal_alignment_index(alignment, sequences)
    }
  }
  if (!inherits(index, "fungal_alignment_index")) stop("index must be produced by fungal_alignment_index().")
  required <- c("pair_id", "accession", "forward_start", "forward_end",
                "reverse_start", "reverse_end", "amplicon_bp", "status")
  if (!all(required %in% names(pair_sequences))) stop("Pair evidence is missing raw site coordinates or status.")
  rows <- as.data.frame(pair_sequences, stringsAsFactors = FALSE)
  if (!is.null(pair_ids)) rows <- rows[rows$pair_id %in% pair_ids, , drop = FALSE]
  if (!is.null(accessions)) rows <- rows[rows$accession %in% accessions, , drop = FALSE]
  coordinate_columns <- c("forward_start_col", "forward_end_col", "reverse_start_col", "reverse_end_col",
                          "forward_left_col", "forward_right_col", "reverse_left_col", "reverse_right_col",
                          "product_start_col", "product_end_col", "product_alignment_span", "n_gap_columns_in_product")
  for (column in coordinate_columns) rows[[column]] <- rep(NA_integer_, nrow(rows))
  rows$alignment_status <- rep("unresolved_binding_evidence", nrow(rows))
  rows$alignment_template_status <- rep("template_not_in_alignment", nrow(rows))
  rows$forward_alignment_status <- rows$reverse_alignment_status <- rows$alignment_status
  rows$alignment_orientation <- rep(NA_character_, nrow(rows))
  rows$forward_direction <- rows$reverse_direction <- rep(NA_integer_, nrow(rows))
  if (!nrow(rows)) return(rows)
  metadata_index <- match(rows$accession, index$metadata$accession)
  for (i in seq_len(nrow(rows))) {
    matched <- metadata_index[i]
    if (is.na(matched)) {
      rows$alignment_status[i] <- "template_not_in_alignment"
      next
    }
    metadata <- index$metadata[matched, , drop = FALSE]
    rows$alignment_template_status[i] <- metadata$status
    rows$alignment_orientation[i] <- metadata$orientation
    rows$forward_direction[i] <- metadata$strand
    rows$reverse_direction[i] <- -metadata$strand
    if (metadata$status != "verified") {
      rows$alignment_status[i] <- metadata$status
      rows$forward_alignment_status[i] <- rows$reverse_alignment_status[i] <- metadata$status
      next
    }
    coordinates <- unlist(rows[i, c("forward_start", "forward_end", "reverse_start", "reverse_end")], use.names = FALSE)
    map <- index$maps[[rows$accession[i]]]
    scorable <- !is.na(rows$status[i]) && rows$status[i] == "scorable"
    # An exact component can remain visible when the opposite component is
    # missing. Low-specificity, ambiguous or undocumented coordinates stay
    # unplaced; a partial pair never receives a product span or a length tier.
    for (j in seq_len(2L)) {
      region <- c("forward", "reverse")[j]
      evidence_column <- paste0(region, "_status")
      exact <- if (evidence_column %in% names(rows)) {
        !is.na(rows[[evidence_column]][i]) && rows[[evidence_column]][i] == "exact_compatible_site"
      } else scorable
      if (!exact) next
      endpoints <- coordinates[(j - 1L) * 2L + seq_len(2L)]
      if (any(!is.finite(endpoints)) || any(endpoints != floor(endpoints)) ||
          any(endpoints < 1L) || any(endpoints > length(map)) || endpoints[1] > endpoints[2]) {
        rows[[paste0(region, "_alignment_status")]][i] <- "invalid_raw_coordinates"
        next
      }
      mapped_site <- unname(map[endpoints])
      rows[[paste0(region, "_start_col")]][i] <- mapped_site[1]
      rows[[paste0(region, "_end_col")]][i] <- mapped_site[2]
      rows[[paste0(region, "_left_col")]][i] <- min(mapped_site)
      rows[[paste0(region, "_right_col")]][i] <- max(mapped_site)
      rows[[paste0(region, "_alignment_status")]][i] <- "mapped"
    }
    if (!scorable) {
      if (rows$forward_alignment_status[i] == "mapped" || rows$reverse_alignment_status[i] == "mapped") {
        rows$alignment_status[i] <- "partial_exact_sites"
      }
      next
    }
    if (rows$forward_alignment_status[i] != "mapped" || rows$reverse_alignment_status[i] != "mapped") {
      rows$alignment_status[i] <- "invalid_binding_site_evidence"
      next
    }
    if (any(!is.finite(coordinates)) || any(coordinates != floor(coordinates)) ||
        any(coordinates < 1L) || any(coordinates > length(map)) ||
        coordinates[1] > coordinates[2] || coordinates[3] > coordinates[4] || coordinates[2] >= coordinates[3]) {
      rows$alignment_status[i] <- "invalid_raw_coordinates"
      next
    }
    product_bp <- coordinates[4] - coordinates[1] + 1L
    if (!is.finite(rows$amplicon_bp[i]) || rows$amplicon_bp[i] != product_bp) {
      rows$alignment_status[i] <- "raw_product_length_mismatch"
      next
    }
    mapped <- unname(map[coordinates])
    for (j in seq_len(4L)) rows[[coordinate_columns[j]]][i] <- mapped[j]
    rows$forward_left_col[i] <- min(mapped[1:2])
    rows$forward_right_col[i] <- max(mapped[1:2])
    rows$reverse_left_col[i] <- min(mapped[3:4])
    rows$reverse_right_col[i] <- max(mapped[3:4])
    rows$product_start_col[i] <- min(mapped[c(1, 4)])
    rows$product_end_col[i] <- max(mapped[c(1, 4)])
    rows$product_alignment_span[i] <- rows$product_end_col[i] - rows$product_start_col[i] + 1L
    rows$n_gap_columns_in_product[i] <- rows$product_alignment_span[i] - product_bp
    rows$alignment_status[i] <- "mapped"
  }
  rownames(rows) <- NULL
  rows
}

fungal_alignment_window <- function(index, mapped_rows = NULL, accessions = NULL,
                                    flank = 20L, max_columns = 120L, start = NULL) {
  if (!inherits(index, "fungal_alignment_index")) stop("index must be produced by fungal_alignment_index().")
  if (length(max_columns) != 1L || !is.finite(max_columns) || max_columns < 1L ||
      length(flank) != 1L || !is.finite(flank) || flank < 0L) stop("Window width and flank must be finite nonnegative lengths, with a positive width.")
  max_columns <- as.integer(max_columns)
  if (is.null(accessions)) accessions <- index$metadata$accession
  accessions <- unique(accessions[accessions %in% index$metadata$accession])
  if (is.null(start)) {
    support <- if (!is.null(mapped_rows) && nrow(mapped_rows)) {
      mapped_rows$product_start_col[mapped_rows$alignment_status == "mapped" & mapped_rows$accession %in% accessions]
    } else numeric()
    start <- if (length(support)) min(support, na.rm = TRUE) - flank else 1L
  }
  if (length(start) != 1L || !is.finite(start)) stop("Window start must be a finite alignment column.")
  last_start <- max(1L, index$n_columns - max_columns + 1L)
  start <- as.integer(max(1L, min(start, last_start)))
  end <- min(index$n_columns, start + max_columns - 1L)
  columns <- if (end >= start) seq.int(start, end) else integer()
  rows <- index$metadata[match(accessions, index$metadata$accession), , drop = FALSE]
  rows$sequence <- if (length(columns)) substring(index$alignment[accessions], start, end) else rep("", length(accessions))
  bases <- matrix("", nrow = length(accessions), ncol = length(columns),
                  dimnames = list(accessions, as.character(columns)))
  if (length(columns) && length(accessions)) {
    bases[] <- do.call(rbind, strsplit(rows$sequence, "", fixed = TRUE))
  }
  overlays <- if (is.null(mapped_rows)) data.frame() else {
    mapped_rows[mapped_rows$accession %in% accessions, , drop = FALSE]
  }
  if (nrow(overlays)) {
    for (region in c("product", "forward", "reverse")) {
      left <- if (region == "product") overlays$product_start_col else overlays[[paste0(region, "_left_col")]]
      right <- if (region == "product") overlays$product_end_col else overlays[[paste0(region, "_right_col")]]
      supported <- if (region == "product") overlays$alignment_status == "mapped" else {
        overlays[[paste0(region, "_alignment_status")]] == "mapped"
      }
      visible <- supported & is.finite(left) & is.finite(right) & right >= start & left <= end
      overlays[[paste0(region, "_visible")]] <- visible
      overlays[[paste0(region, "_window_start_col")]] <- ifelse(visible, pmax(left, start), NA_integer_)
      overlays[[paste0(region, "_window_end_col")]] <- ifelse(visible, pmin(right, end), NA_integer_)
    }
  }
  list(columns = columns, start_col = start, end_col = end, n_columns = index$n_columns,
       rows = rows, bases = bases, overlays = overlays)
}
