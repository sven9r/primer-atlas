# Fungal rDNA combination discovery. Coordinates and lengths belong to each
# ungapped template; no coordinate-reference placement is used to fill a site.
# Source R/functions.R before this module for the established IUPAC helpers.

fungal_catalog <- function(path = "data/catalog/unite_primers.csv") {
  catalog <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("primer_name", "sequence", "direction", "gene_locus")
  if (!all(required %in% names(catalog))) {
    stop("Fungal catalog is missing: ", paste(setdiff(required, names(catalog)), collapse = ", "))
  }
  catalog$source_direction <- catalog$direction
  catalog$direction_note <- "Source direction retained."
  # The sequence below is the reverse complement of ITS3, and anneals as the
  # reverse ITS1 primer. Preserve the imported UNITE direction for inspection.
  corrected <- catalog$primer_name == "ITS2" &
    toupper(catalog$sequence) == "GCTGCGTTCTTCATCGATGC"
  catalog$direction[corrected] <- "reverse"
  catalog$direction_note[corrected] <- paste(
    "Derived orientation correction: ITS2 is reverse; its sequence is the",
    "reverse complement of ITS3. Imported source direction is preserved."
  )
  catalog <- catalog[
    catalog$gene_locus %in% c("SSU", "SSU/ITS", "ITS", "ITS/LSU", "LSU") &
      catalog$direction %in% c("forward", "reverse"), , drop = FALSE
  ]
  catalog$sequence <- toupper(gsub("[[:space:]]", "", catalog$sequence))
  catalog$catalog_id <- paste(
    catalog$primer_name, catalog$direction, catalog$sequence, sep = "|"
  )
  rownames(catalog) <- NULL
  catalog
}

.fungal_documented_pairs <- function() {
  data.frame(
    pair_id = c("ITS1_ITS4", "ITS1F_ITS4", "ITS1F_ITS4B", "ITS3_ITS4",
                "FITS7_ITS4", "GITS7_ITS4", "ITS86F_ITS4", "ITS1F_LR21"),
    forward_primer = c("ITS1", "ITS1F", "ITS1F", "ITS3", "fITS7", "gITS7", "ITS86F", "ITS1F"),
    reverse_primer = c("ITS4", "ITS4", "ITS4B", "ITS4", "ITS4", "ITS4", "ITS4", "LR21"),
    # Freeze exact component identities as well as names: an updated oligo
    # cannot inherit the documented pair label merely by keeping its name.
    forward_sequence = c("TCCGTAGGTGAACCTGCGG", "CTTGGTCATTTAGAGGAAGTAA",
                         "CTTGGTCATTTAGAGGAAGTAA", "GCATCGATGAAGAACGCAGC",
                         "GTGARTCATCGAATCTTTG", "GTGARTCATCGARTCTTTG",
                         "GTGAATCATCGAATCTTTGAA", "CTTGGTCATTTAGAGGAAGTAA"),
    reverse_sequence = c("TCCTCCGCTTATTGATATGC", "TCCTCCGCTTATTGATATGC",
                         "CAGGAGACTTGTACACGGTCCAG", "TCCTCCGCTTATTGATATGC",
                         "TCCTCCGCTTATTGATATGC", "TCCTCCGCTTATTGATATGC",
                         "TCCTCCGCTTATTGATATGC", "ACTTCAAGCGTTTCCCTTT"),
    stringsAsFactors = FALSE
  )
}

.fungal_empty_sites <- function() {
  data.frame(catalog_id = character(), primer_name = character(),
             direction = character(), accession = character(), start = integer(),
             end = integer(), n_matches = integer(), status = character(),
             stringsAsFactors = FALSE)
}

.fungal_empty_pair_sequences <- function() {
  data.frame(pair_id = character(), accession = character(),
             forward_start = integer(), forward_end = integer(),
             reverse_start = integer(), reverse_end = integer(),
             amplicon_bp = integer(), status = character(),
             forward_status = character(), reverse_status = character(),
             stringsAsFactors = FALSE)
}

.fungal_empty_pairs <- function() {
  data.frame(pair_id = character(), pair_label = character(),
             forward_primer = character(), reverse_primer = character(),
             forward_sequence = character(), reverse_sequence = character(),
             target = character(), evidence_status = character(),
             source_citation = character(), source_url = character(),
             region = character(), n_panel = integer(), n_scorable = integer(),
             amplicon_min_bp = numeric(), amplicon_median_bp = numeric(),
             amplicon_max_bp = numeric(), stringsAsFactors = FALSE)
}

.fungal_component_text <- function(catalog, i, j, column, separator = " + ") {
  if (!column %in% names(catalog)) return("")
  values <- as.character(catalog[[column]][c(i, j)])
  paste(unique(values[!is.na(values) & nzchar(values)]), collapse = separator)
}

.fungal_exact_site <- function(reference_bases, primer_sequence, direction) {
  empty <- list(start = NA_integer_, end = NA_integer_, n_matches = 0L,
                status = "no_exact_binding_site")
  if (!valid_primer_input(primer_sequence)) {
    empty$status <- "invalid_oligo"
    return(empty)
  }
  query <- split_bases(if (direction == "reverse") {
    reverse_complement(primer_sequence)
  } else primer_sequence)
  width <- length(query)
  if (length(reference_bases) < width) {
    empty$status <- "template_shorter_than_oligo"
    return(empty)
  }
  starts <- seq_len(length(reference_bases) - width + 1L)
  compatible <- rep(TRUE, length(starts))
  for (position in seq_along(query)) {
    observed <- reference_bases[starts + position - 1L]
    # Template ambiguity never establishes an exact binding site. Primer
    # degeneracy and inosine retain the established repository interpretation.
    compatible <- compatible & observed %in% iupac[[query[position]]]
    if (!any(compatible)) break
  }
  hits <- starts[compatible]
  if (!length(hits)) return(empty)
  if (length(hits) > 1L) {
    empty$n_matches <- length(hits)
    empty$status <- "ambiguous_binding_site"
    return(empty)
  }
  candidate <- data.frame(
    compatible_identity = 1, terminal_3_mismatches = 0L,
    longest_compatible_run = width,
    random_match_evalue = random_match_evalue(query, width, length(reference_bases))
  )
  if (!isTRUE(assess_primer_binding(primer_sequence, candidate)$credible)) {
    empty$n_matches <- 1L
    empty$status <- "low_specificity_site"
    return(empty)
  }
  list(start = as.integer(hits), end = as.integer(hits + width - 1L),
       n_matches = 1L, status = "exact_compatible_site")
}

.fungal_add_taxonomy <- function(rows, taxonomy) {
  columns <- c("class", "phylum", "organism")
  if (is.null(taxonomy) || !nrow(taxonomy) || !"accession" %in% names(taxonomy)) {
    for (column in columns) rows[[column]] <- rep(NA_character_, nrow(rows))
    return(rows)
  }
  if (anyDuplicated(taxonomy$accession)) {
    stop("Taxonomy accession identifiers must be unique.")
  }
  matched <- match(rows$accession, taxonomy$accession)
  for (column in columns) {
    rows[[column]] <- if (column %in% names(taxonomy)) {
      as.character(taxonomy[[column]][matched])
    } else rep(NA_character_, nrow(rows))
  }
  rows
}

fungal_evaluate_panel <- function(catalog, sequences, taxonomy = NULL,
                                min_amplicon_bp = 100L) {
  required_helpers <- c("valid_primer_input", "split_bases", "reverse_complement",
                        "random_match_evalue", "assess_primer_binding")
  if (!all(vapply(required_helpers, exists, logical(1), mode = "function"))) {
    stop("Source R/functions.R before evaluating the fungal panel.")
  }
  if (!is.character(sequences) || is.null(names(sequences)) ||
      anyNA(names(sequences)) || any(!nzchar(names(sequences))) ||
      anyDuplicated(names(sequences)) || anyNA(sequences)) {
    stop("Sequences must be a named character vector with unique accession identifiers and no missing sequences.")
  }
  if (length(min_amplicon_bp) != 1L || !is.finite(min_amplicon_bp) || min_amplicon_bp < 1) {
    stop("min_amplicon_bp must be a positive length.")
  }
  required <- c("catalog_id", "primer_name", "sequence", "direction", "gene_locus")
  if (!all(required %in% names(catalog)) || anyDuplicated(catalog$catalog_id)) {
    stop("Catalog requires unique catalog_id values and primer name, sequence, direction, and locus columns.")
  }
  catalog <- catalog[catalog$direction %in% c("forward", "reverse") &
                       catalog$gene_locus %in% c("SSU", "SSU/ITS", "ITS", "ITS/LSU", "LSU"), , drop = FALSE]
  rownames(catalog) <- NULL
  # Alignment gaps are not product bases. Ambiguous template bases are retained
  # in the product length but cannot substantiate either complete binding site.
  sequences <- toupper(gsub("[-[:space:]]", "", sequences))
  reference_bases <- lapply(sequences, function(sequence) strsplit(sequence, "", fixed = TRUE)[[1]])
  site_blocks <- lapply(seq_along(sequences), function(sequence_index) {
    sites <- lapply(seq_len(nrow(catalog)), function(primer_index) {
      .fungal_exact_site(reference_bases[[sequence_index]],
                         catalog$sequence[primer_index], catalog$direction[primer_index])
    })
    if (!length(sites)) return(.fungal_empty_sites())
    data.frame(
      catalog_id = catalog$catalog_id, primer_name = catalog$primer_name,
      direction = catalog$direction, accession = names(sequences)[sequence_index],
      start = vapply(sites, `[[`, integer(1), "start"),
      end = vapply(sites, `[[`, integer(1), "end"),
      n_matches = vapply(sites, `[[`, integer(1), "n_matches"),
      status = vapply(sites, `[[`, character(1), "status"),
      stringsAsFactors = FALSE
    )
  })
  sites <- if (length(site_blocks)) do.call(rbind, site_blocks) else .fungal_empty_sites()
  forward_indices <- which(catalog$direction == "forward")
  reverse_indices <- which(catalog$direction == "reverse")
  combinations <- expand.grid(forward = forward_indices, reverse = reverse_indices,
                              KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  if (!nrow(combinations)) {
    return(list(sites = .fungal_add_taxonomy(sites, taxonomy),
                pair_sequences = .fungal_add_taxonomy(.fungal_empty_pair_sequences(), taxonomy),
                pairs = .fungal_empty_pairs()))
  }
  forward <- catalog[combinations$forward, , drop = FALSE]
  reverse <- catalog[combinations$reverse, , drop = FALSE]
  pair_ids <- paste0("FUNGAL_", gsub("[^A-Z0-9]+", "_", toupper(forward$primer_name)),
                     "_", gsub("[^A-Z0-9]+", "_", toupper(reverse$primer_name)))
  duplicated_ids <- duplicated(pair_ids) | duplicated(pair_ids, fromLast = TRUE)
  pair_ids[duplicated_ids] <- paste0(pair_ids[duplicated_ids], "_", forward$sequence[duplicated_ids],
                                   "_", reverse$sequence[duplicated_ids])
  documented <- .fungal_documented_pairs()
  keys <- paste(forward$primer_name, reverse$primer_name,
                forward$sequence, reverse$sequence, sep = "|")
  documented_keys <- paste(documented$forward_primer, documented$reverse_primer,
                           documented$forward_sequence, documented$reverse_sequence, sep = "|")
  documented_index <- match(keys, documented_keys)
  documented_match <- !is.na(documented_index)
  pair_ids[documented_match] <- documented$pair_id[documented_index[documented_match]]

  sequence_blocks <- lapply(seq_along(sequences), function(sequence_index) {
    template_sites <- site_blocks[[sequence_index]]
    f <- template_sites[combinations$forward, , drop = FALSE]
    r <- template_sites[combinations$reverse, , drop = FALSE]
    f_available <- f$status == "exact_compatible_site"
    r_available <- r$status == "exact_compatible_site"
    status <- rep("both_sites_missing", nrow(combinations))
    status[f_available & !r_available] <- "reverse_site_missing"
    status[!f_available & r_available] <- "forward_site_missing"
    status[f$status == "low_specificity_site" | r$status == "low_specificity_site"] <- "low_specificity_site"
    status[f$status == "ambiguous_binding_site" | r$status == "ambiguous_binding_site"] <- "ambiguous_binding_site"
    both <- f_available & r_available
    status[both] <- "outward_or_overlapping_sites"
    inward <- both & !is.na(f$end) & !is.na(r$start) & f$end < r$start
    product_length <- r$end - f$start + 1L
    status[inward] <- "below_minimum_length"
    status[inward & product_length >= min_amplicon_bp] <- "scorable"
    product_length[status != "scorable"] <- NA_integer_
    data.frame(pair_id = pair_ids, accession = names(sequences)[sequence_index],
               forward_start = f$start, forward_end = f$end,
               reverse_start = r$start, reverse_end = r$end,
               amplicon_bp = product_length, status = status,
               forward_status = f$status, reverse_status = r$status,
               stringsAsFactors = FALSE)
  })
  sequence_rows <- if (length(sequence_blocks)) do.call(rbind, sequence_blocks) else .fungal_empty_pair_sequences()
  supported_ids <- unique(sequence_rows$pair_id[sequence_rows$status == "scorable"])
  retained <- documented_match | pair_ids %in% supported_ids
  sequence_rows <- sequence_rows[sequence_rows$pair_id %in% pair_ids[retained], , drop = FALSE]
  pair_rows <- lapply(which(retained), function(i) {
    fi <- combinations$forward[i]
    ri <- combinations$reverse[i]
    lengths <- sequence_rows$amplicon_bp[
      sequence_rows$pair_id == pair_ids[i] & sequence_rows$status == "scorable"
    ]
    component_status <- if ("primary_reference_status" %in% names(catalog)) {
      catalog$primary_reference_status[c(fi, ri)]
    } else character()
    evidence <- if (documented_match[i]) "documented_pair" else if ("unpublished" %in% component_status) {
      "candidate_unpublished_component"
    } else if (length(component_status) == 2L && all(component_status == "publication_cited")) {
      "candidate_published_components"
    } else "candidate_incomplete_provenance"
    data.frame(
      pair_id = pair_ids[i], pair_label = paste(catalog$primer_name[fi], catalog$primer_name[ri], sep = " + "),
      forward_primer = catalog$primer_name[fi], reverse_primer = catalog$primer_name[ri],
      forward_sequence = catalog$sequence[fi], reverse_sequence = catalog$sequence[ri],
      target = .fungal_component_text(catalog, fi, ri, "target"),
      evidence_status = evidence,
      source_citation = .fungal_component_text(catalog, fi, ri, "reference"),
      source_url = .fungal_component_text(catalog, fi, ri, "primary_reference_url", " | "),
      region = .fungal_component_text(catalog, fi, ri, "gene_locus"),
      n_panel = length(sequences), n_scorable = length(lengths),
      amplicon_min_bp = if (length(lengths)) min(lengths) else NA_real_,
      amplicon_median_bp = if (length(lengths)) median(lengths) else NA_real_,
      amplicon_max_bp = if (length(lengths)) max(lengths) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  pairs <- if (length(pair_rows)) do.call(rbind, pair_rows) else .fungal_empty_pairs()
  if (nrow(pairs)) pairs <- pairs[order(is.na(pairs$amplicon_median_bp), pairs$amplicon_median_bp, pairs$pair_label), , drop = FALSE]
  rownames(pairs) <- rownames(sequence_rows) <- rownames(sites) <- NULL
  list(sites = .fungal_add_taxonomy(sites, taxonomy),
       pair_sequences = .fungal_add_taxonomy(sequence_rows, taxonomy), pairs = pairs)
}

fungal_subset_summary <- function(pair_sequences, pairs, accessions) {
  accessions <- unique(as.character(accessions))
  accessions <- accessions[!is.na(accessions) & nzchar(accessions)]
  result <- pairs
  result$n_panel <- rep(length(accessions), nrow(result))
  result$n_scorable <- integer(nrow(result))
  result$amplicon_min_bp <- result$amplicon_median_bp <- result$amplicon_max_bp <- rep(NA_real_, nrow(result))
  selected <- pair_sequences[pair_sequences$accession %in% accessions, , drop = FALSE]
  for (i in seq_len(nrow(result))) {
    rows <- selected[selected$pair_id == result$pair_id[i], , drop = FALSE]
    lengths <- rows$amplicon_bp[rows$status == "scorable" & is.finite(rows$amplicon_bp) & rows$amplicon_bp > 0]
    result$n_scorable[i] <- length(lengths)
    if (!length(lengths)) next
    result$amplicon_min_bp[i] <- min(lengths)
    result$amplicon_median_bp[i] <- median(lengths)
    result$amplicon_max_bp[i] <- max(lengths)
  }
  result
}

fungal_platform_summary <- function(pair_sequences, pairs, platform = "ALL",
                                    illumina_max = 500, pacbio_max = 2000,
                                    nanopore_max = 3500) {
  thresholds <- c(illumina_max, pacbio_max, nanopore_max)
  if (length(thresholds) != 3L || any(!is.finite(thresholds)) ||
      any(thresholds <= 0) || any(diff(thresholds) < 0)) {
    stop("Platform length ceilings must be positive and ordered Illumina <= PacBio <= Nanopore.")
  }
  platform <- toupper(platform)
  if (length(platform) != 1L || !platform %in% c("ALL", "ILLUMINA", "PACBIO", "NANOPORE")) {
    stop("platform must be ALL, Illumina, PacBio, or Nanopore.")
  }
  result <- pairs
  result$length_tier <- rep("Unresolved", nrow(result))
  result$tier_color <- rep("#94a3b8", nrow(result))
  result$platforms <- rep("Length evidence unavailable", nrow(result))
  result$length_range_crosses_tiers <- rep(FALSE, nrow(result))
  result$length_tier_range <- rep("Unresolved", nrow(result))
  result$length_eligible_fraction <- rep(NA_real_, nrow(result))
  result$panel_eligible_fraction <- rep(NA_real_, nrow(result))
  result$n_length_eligible <- result$n_illumina_eligible <- result$n_pacbio_eligible <-
    result$n_nanopore_eligible <- result$n_length_observed <- integer(nrow(result))
  result$platform_selected <- rep(platform, nrow(result))
  for (i in seq_len(nrow(result))) {
    rows <- pair_sequences[pair_sequences$pair_id == result$pair_id[i], , drop = FALSE]
    lengths <- rows$amplicon_bp[rows$status == "scorable" & is.finite(rows$amplicon_bp) & rows$amplicon_bp > 0]
    observed <- length(lengths)
    result$n_length_observed[i] <- observed
    if (!observed) next
    counts <- c(Illumina = sum(lengths <= illumina_max), PacBio = sum(lengths <= pacbio_max),
                Nanopore = sum(lengths <= nanopore_max))
    result$n_illumina_eligible[i] <- counts[1]
    result$n_pacbio_eligible[i] <- counts[2]
    result$n_nanopore_eligible[i] <- counts[3]
    result$n_length_eligible[i] <- counts[match(platform, c("ILLUMINA", "PACBIO", "NANOPORE"), nomatch = 3L)]
    result$length_eligible_fraction[i] <- result$n_length_eligible[i] / observed
    if ("n_panel" %in% names(result) && is.finite(result$n_panel[i]) && result$n_panel[i] > 0) {
      result$panel_eligible_fraction[i] <- result$n_length_eligible[i] / result$n_panel[i]
    }
    assign_tier <- function(bp) {
      if (bp <= illumina_max) 1L else if (bp <= pacbio_max) 2L else if (bp <= nanopore_max) 3L else 4L
    }
    tier <- assign_tier(median(lengths))
    tier_names <- c("Illumina", "PacBio", "Nanopore", "Beyond configured limits")
    minimum_tier <- assign_tier(min(lengths))
    maximum_tier <- assign_tier(max(lengths))
    # Color describes the observed median only. Range crossings and numerator /
    # denominator counts make heterogeneous lengths and incomplete panels visible.
    result$length_range_crosses_tiers[i] <- minimum_tier != maximum_tier
    result$length_tier_range[i] <- paste(unique(tier_names[c(minimum_tier, maximum_tier)]), collapse = " → ")
    result$length_tier[i] <- tier_names[tier]
    result$tier_color[i] <- c("#2563eb", "#7c3aed", "#0f766e", "#64748b")[tier]
    eligible <- counts > 0
    result$platforms[i] <- if (any(eligible)) {
      paste0(names(counts)[eligible], " (", counts[eligible], "/", observed, ")", collapse = " · ")
    } else "Exceeds all configured length ceilings"
  }
  result
}
