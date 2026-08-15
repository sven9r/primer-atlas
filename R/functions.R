iupac <- list(
  A = "A", C = "C", G = "G", T = "T", U = "T",
  R = c("A", "G"), Y = c("C", "T"), S = c("G", "C"),
  W = c("A", "T"), K = c("G", "T"), M = c("A", "C"),
  B = c("C", "G", "T"), D = c("A", "G", "T"),
  H = c("A", "C", "T"), V = c("A", "C", "G"),
  N = c("A", "C", "G", "T"), I = c("A", "C", "G", "T"),
  `-` = character()
)

iupac_codes <- names(iupac)
iupac_match_matrix <- outer(
  iupac_codes,
  iupac_codes,
  Vectorize(function(a, b) length(intersect(iupac[[a]], iupac[[b]])) > 0)
)
dimnames(iupac_match_matrix) <- list(iupac_codes, iupac_codes)

mismatch_multiplier_matrix <- outer(
  iupac_codes,
  iupac_codes,
  Vectorize(function(primer_base, template_base) {
    if (iupac_match_matrix[primer_base, template_base]) return(0)
    p <- iupac[[primer_base]]
    t <- iupac[[template_base]]
    if (!length(p) || !length(t)) return(2)
    purine <- c("A", "G")
    mean(outer(
      p, t,
      Vectorize(function(x, y) {
        if (x == y) return(0)
        if ((x %in% purine) == (y %in% purine)) return(2)
        0.5
      })
    ))
  })
)
dimnames(mismatch_multiplier_matrix) <- list(iupac_codes, iupac_codes)

clean_sequence <- function(x) {
  gsub("[^ACGTURYSWKMBDHVNI-]", "", toupper(x))
}

normalize_geo_loc_name <- function(value) {
  value <- trimws(ifelse(is.na(value), "", as.character(value)))
  value[value %in% c("", "missing", "not applicable", "not provided")] <- ""
  pieces <- strsplit(value, ":", fixed = TRUE)
  country <- vapply(pieces, function(x) {
    if (length(x) && nzchar(trimws(x[1]))) trimws(x[1]) else NA_character_
  }, character(1))
  locality <- vapply(pieces, function(x) {
    if (length(x) > 1L) trimws(paste(x[-1], collapse = ":")) else NA_character_
  }, character(1))
  locality[!nzchar(ifelse(is.na(locality), "", locality))] <- NA_character_
  normalized_location <- ifelse(
    is.na(country), NA_character_,
    ifelse(is.na(locality), country, paste(country, locality, sep = " · "))
  )
  data.frame(
    geo_loc_name_raw = ifelse(nzchar(value), value, NA_character_),
    country_or_territory = country,
    locality = locality,
    normalized_location = normalized_location,
    geographic_resolution = ifelse(
      is.na(country), "unresolved", ifelse(is.na(locality), "country", "locality")
    ),
    missing_location_reason = ifelse(
      is.na(country), "No usable geo_loc_name or legacy country qualifier", ""
    ),
    stringsAsFactors = FALSE
  )
}

parse_lat_lon <- function(value) {
  value <- trimws(ifelse(is.na(value), "", as.character(value)))
  pattern <- paste0(
    "^([0-9]+(?:\\.[0-9]+)?)\\s*([NS])\\s+",
    "([0-9]+(?:\\.[0-9]+)?)\\s*([EW])$"
  )
  matched <- regexec(pattern, toupper(value), perl = TRUE)
  pieces <- regmatches(toupper(value), matched)
  latitude <- longitude <- rep(NA_real_, length(value))
  for (i in seq_along(pieces)) {
    if (length(pieces[[i]]) != 5L) next
    latitude[i] <- as.numeric(pieces[[i]][2]) * ifelse(pieces[[i]][3] == "S", -1, 1)
    longitude[i] <- as.numeric(pieces[[i]][4]) * ifelse(pieces[[i]][5] == "W", -1, 1)
  }
  data.frame(latitude = latitude, longitude = longitude)
}

filter_pair_facets <- function(pair_ids, facets, selections) {
  result <- unique(pair_ids)
  for (facet_type in names(selections)) {
    selected <- selections[[facet_type]]
    if (is.null(selected) || !length(selected)) next
    allowed <- unique(facets$pair_id[
      facets$facet_type == facet_type & facets$facet_value %in% selected
    ])
    result <- intersect(result, allowed)
  }
  result
}

split_bases <- function(x) strsplit(clean_sequence(x), "", fixed = TRUE)[[1]]

reverse_complement <- function(x) {
  comp <- c(
    A = "T", T = "A", U = "A", C = "G", G = "C",
    R = "Y", Y = "R", S = "S", W = "W", K = "M", M = "K",
    B = "V", V = "B", D = "H", H = "D", N = "N", I = "N",
    `-` = "-"
  )
  paste(rev(unname(comp[split_bases(x)])), collapse = "")
}

degeneracy <- function(x) {
  sizes <- vapply(split_bases(x), function(b) length(iupac[[b]]), numeric(1))
  if (any(sizes == 0)) return(NA_real_)
  prod(sizes)
}

gc_bounds <- function(x) {
  chars <- split_bases(x)
  possible <- lapply(chars, function(b) iupac[[b]])
  min_gc <- sum(vapply(possible, function(z) all(z %in% c("G", "C")), logical(1)))
  max_gc <- sum(vapply(possible, function(z) any(z %in% c("G", "C")), logical(1)))
  c(min = min_gc, max = max_gc)
}

tm_range <- function(x, sodium_mM = 50) {
  x <- clean_sequence(x)
  n <- nchar(x)
  bounds <- gc_bounds(x)
  sodium_M <- max(sodium_mM, 0.1) / 1000
  estimate <- function(gc) {
    81.5 + 16.6 * log10(sodium_M) + 0.41 * (100 * gc / n) - 600 / n
  }
  sort(unname(vapply(bounds, estimate, numeric(1))))
}

bases_match <- function(a, b) {
  isTRUE(iupac_match_matrix[a, b])
}

position_penalty <- function(distance_from_3prime) {
  reference <- c(
    242.4, 202.8, 169.8, 142.4, 119.5, 100.4, 84.5, 71.2,
    60.2, 51.0, 43.3, 36.9, 31.6, 27.2, 23.5, 20.4, 17.8,
    15.7, 13.9, 12.4, 11.2, 10.2, 9.3, 8.6, 8.0, 7.5, 7.1,
    6.7, 6.4, 6.2
  )
  reference[pmin(pmax(distance_from_3prime, 1), length(reference))]
}

mismatch_type_multiplier <- function(primer_base, template_base) {
  unname(mismatch_multiplier_matrix[primer_base, template_base])
}

normalize_primer_input <- function(x) {
  toupper(gsub("[[:space:]]", "", x))
}

valid_primer_input <- function(x, min_length = 15L, max_length = 80L) {
  normalized <- normalize_primer_input(x)
  nchar(normalized) >= min_length &&
    nchar(normalized) <= max_length &&
    grepl("^[ACGTURYSWKMBDHVNI]+$", normalized)
}

longest_true_run <- function(x) {
  if (!length(x) || !any(x)) return(0L)
  max(rle(x)$lengths[rle(x)$values])
}

random_match_evalue <- function(primer_bases, compatible_matches, search_width) {
  per_base_probability <- vapply(
    primer_bases,
    function(base) length(iupac[[base]]) / 4,
    numeric(1)
  )
  distribution <- 1
  for (probability in per_base_probability) {
    distribution <- c(distribution * (1 - probability), 0) +
      c(0, distribution * probability)
  }
  probability_at_least_observed <- sum(
    distribution[(compatible_matches + 1L):length(distribution)]
  )
  max(1L, search_width - length(primer_bases) + 1L) *
    probability_at_least_observed
}

assess_primer_binding <- function(
  primer_sequence,
  candidate,
  minimum_identity = 0.80,
  maximum_terminal_3_mismatches = 1L,
  minimum_compatible_run = 7L,
  maximum_random_match_evalue = 0.05
) {
  reasons <- character()
  if (candidate$compatible_identity < minimum_identity) {
    reasons <- c(
      reasons,
      sprintf(
        "compatible identity %.1f%% is below %.0f%%",
        candidate$compatible_identity * 100,
        minimum_identity * 100
      )
    )
  }
  if (candidate$terminal_3_mismatches > maximum_terminal_3_mismatches) {
    reasons <- c(
      reasons,
      paste0(
        candidate$terminal_3_mismatches,
        " mismatch(es) occur in the terminal 3 bases"
      )
    )
  }
  if (candidate$longest_compatible_run < minimum_compatible_run) {
    reasons <- c(
      reasons,
      paste0(
        "longest compatible run is ", candidate$longest_compatible_run,
        " bp (minimum ", minimum_compatible_run, ")"
      )
    )
  }
  if (candidate$random_match_evalue > maximum_random_match_evalue) {
    reasons <- c(
      reasons,
      paste0(
        "chance-match expectation is ",
        format(candidate$random_match_evalue, digits = 2, scientific = TRUE),
        " across the COI reference (maximum ",
        maximum_random_match_evalue, ")"
      )
    )
  }
  list(
    credible = !length(reasons),
    reasons = reasons,
    summary = paste(reasons, collapse = "; ")
  )
}

alignment_binding_support <- function(
  alignment_templates,
  primer_sequence,
  start,
  end,
  direction,
  minimum_identity = 0.80,
  maximum_terminal_3_mismatches = 1L,
  minimum_compatible_run = 7L
) {
  if (!length(alignment_templates)) {
    return(list(summary = data.frame(), templates = data.frame()))
  }
  primer_sequence <- normalize_primer_input(primer_sequence)
  query <- if (direction == "reverse") {
    reverse_complement(primer_sequence)
  } else {
    primer_sequence
  }
  query_bases <- split_bases(query)
  width <- length(query_bases)

  order_rows <- lapply(names(alignment_templates), function(order_name) {
    sequences <- alignment_templates[[order_name]]
    if (!length(sequences)) return(NULL)
    do.call(rbind, lapply(seq_along(sequences), function(i) {
      template_bases <- split_bases(substr(sequences[i], start, end))
      available <- length(template_bases) == width &&
        all(template_bases %in% c("A", "C", "G", "T"))
      compatible <- if (length(template_bases) == width) {
        mapply(bases_match, query_bases, template_bases)
      } else {
        rep(FALSE, width)
      }
      mismatch_indices <- which(!compatible)
      distances <- if (direction == "forward") {
        width - mismatch_indices + 1L
      } else {
        mismatch_indices
      }
      identity <- mean(compatible)
      terminal_3 <- sum(distances <= 3L)
      compatible_run <- longest_true_run(compatible)
      data.frame(
        order = order_name,
        template = names(sequences)[i],
        available = available,
        compatible_identity = identity,
        mismatch_count = sum(!compatible),
        terminal_3_mismatches = terminal_3,
        longest_compatible_run = compatible_run,
        supported = available &&
          identity >= minimum_identity &&
          terminal_3 <= maximum_terminal_3_mismatches &&
          compatible_run >= minimum_compatible_run,
        stringsAsFactors = FALSE
      )
    }))
  })
  template_metrics <- do.call(rbind, order_rows)
  if (is.null(template_metrics) || !nrow(template_metrics)) {
    return(list(summary = data.frame(), templates = data.frame()))
  }
  grouped <- split(template_metrics, template_metrics$order)
  summary <- do.call(rbind, lapply(grouped, function(rows) {
    available <- rows$available %in% TRUE
    supported <- rows$supported %in% TRUE
    data.frame(
      order = rows$order[1],
      n_templates = nrow(rows),
      n_available = sum(available),
      n_supported = sum(supported),
      support_fraction = if (sum(available)) {
        sum(supported) / sum(available)
      } else {
        NA_real_
      },
      median_compatible_identity = if (sum(available)) {
        median(rows$compatible_identity[available])
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL
  summary <- summary[
    order(-summary$n_supported, -summary$support_fraction, summary$order),
    ,
    drop = FALSE
  ]
  list(summary = summary, templates = template_metrics)
}

alignment_support_passes <- function(
  support_summary,
  minimum_supported_templates = 2L,
  minimum_support_fraction = 0.10
) {
  if (!nrow(support_summary)) return(FALSE)
  any(
    support_summary$n_supported >= minimum_supported_templates &
      support_summary$support_fraction >= minimum_support_fraction,
    na.rm = TRUE
  )
}

combine_pair_alignment_support <- function(forward_support, reverse_support) {
  if (!nrow(forward_support$templates) || !nrow(reverse_support$templates)) {
    return(data.frame())
  }
  paired <- merge(
    forward_support$templates,
    reverse_support$templates,
    by = c("order", "template"),
    suffixes = c("_forward", "_reverse")
  )
  grouped <- split(paired, paired$order)
  summary <- do.call(rbind, lapply(grouped, function(rows) {
    available <- rows$available_forward & rows$available_reverse
    supported <- rows$supported_forward & rows$supported_reverse
    data.frame(
      order = rows$order[1],
      n_templates = nrow(rows),
      n_pair_available = sum(available),
      n_pair_supported = sum(supported),
      pair_support_fraction = if (sum(available)) {
        sum(supported) / sum(available)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL
  summary[
    order(
      -summary$n_pair_supported,
      -summary$pair_support_fraction,
      summary$order
    ),
    ,
    drop = FALSE
  ]
}

summarize_order_preference <- function(
  pair_order_support,
  expected_order = NA_character_,
  minimum_available_templates = 3L
) {
  usable <- pair_order_support[
    pair_order_support$n_pair_available >= minimum_available_templates &
      is.finite(pair_order_support$pair_support_fraction),
    ,
    drop = FALSE
  ]
  if (!nrow(usable)) {
    return(list(
      headline = "Order preference could not be assessed.",
      detail = paste(
        "Fewer than", minimum_available_templates,
        "templates had both binding sites in every installed order."
      ),
      expected = "",
      strongest_orders = character(),
      low_orders = character()
    ))
  }
  usable <- usable[
    order(
      -usable$pair_support_fraction,
      -usable$n_pair_supported,
      usable$order
    ),
    ,
    drop = FALSE
  ]
  format_order <- function(row) {
    paste0(
      row$order, " ",
      round(row$pair_support_fraction * 100, 1), "% (",
      row$n_pair_supported, "/", row$n_pair_available, ")"
    )
  }
  strong <- usable[usable$pair_support_fraction >= 0.50, , drop = FALSE]
  strongest <- if (nrow(strong)) {
    head(strong, 4L)
  } else {
    head(usable, 3L)
  }
  strongest_labels <- vapply(
    seq_len(nrow(strongest)),
    function(i) format_order(strongest[i, , drop = FALSE]),
    character(1)
  )
  top_gap <- if (nrow(usable) >= 2L) {
    usable$pair_support_fraction[1] - usable$pair_support_fraction[2]
  } else {
    usable$pair_support_fraction[1]
  }
  headline <- if (
    usable$pair_support_fraction[1] >= 0.50 &&
      top_gap >= 0.20
  ) {
    paste0(
      "Appears to prefer ", usable$order[1],
      " among the installed order alignments."
    )
  } else if (nrow(strong) == 1L) {
    paste0(
      "Strongest compatibility is with ", strong$order[1],
      "; no exclusive order preference is demonstrated."
    )
  } else if (nrow(strong) > 1L) {
    paste0(
      "Strong compatibility is shared across ",
      paste(head(strong$order, 4L), collapse = ", "),
      if (nrow(strong) > 4L) ", and additional orders" else "",
      "."
    )
  } else {
    "No installed order reaches 50% pair compatibility."
  }

  low <- usable[
    usable$pair_support_fraction <= 0.05,
    ,
    drop = FALSE
  ]
  low_labels <- head(low$order, 8L)
  discrimination <- if (length(low_labels)) {
    paste0(
      "Low compatibility suggests discrimination against ",
      paste(low_labels, collapse = ", "),
      if (nrow(low) > length(low_labels)) {
        paste0(", and ", nrow(low) - length(low_labels), " additional order(s)")
      } else {
        ""
      },
      " in this reference set."
    )
  } else {
    "No order shows ≤5% pair compatibility in the available reference set."
  }

  expected_note <- ""
  if (
    length(expected_order) == 1L &&
      !is.na(expected_order) &&
      nzchar(expected_order)
  ) {
    expected_row <- usable[usable$order == expected_order, , drop = FALSE]
    expected_note <- if (!nrow(expected_row)) {
      paste0(
        "Expected target ", expected_order,
        " cannot be assessed because fewer than ",
        minimum_available_templates, " templates contain both sites."
      )
    } else if (expected_row$pair_support_fraction >= 0.50) {
      paste0(
        "Expected target ", format_order(expected_row),
        " is strongly supported",
        if (expected_row$order == usable$order[1]) {
          " and ranks first."
        } else {
          paste0(
            ", but ", usable$order[1],
            " has higher support; the pair is not order-specific here."
          )
        }
      )
    } else if (expected_row$pair_support_fraction >= 0.10) {
      paste0(
        "Expected target ", format_order(expected_row),
        " has only partial support."
      )
    } else {
      paste0(
        "Warning: expected target ", format_order(expected_row),
        " has low support in the installed alignment."
      )
    }
  }

  list(
    headline = headline,
    detail = paste0(
      "Strongest observed compatibility: ",
      paste(strongest_labels, collapse = "; "),
      ". ", discrimination
    ),
    expected = expected_note,
    strongest_orders = strongest$order,
    low_orders = low$order
  )
}

primer_candidate_table <- function(reference_sequence, primer_sequence, direction) {
  if (!direction %in% c("forward", "reverse")) {
    stop("direction must be 'forward' or 'reverse'.")
  }
  reference_sequence <- clean_sequence(reference_sequence)
  primer_sequence <- normalize_primer_input(primer_sequence)
  if (!valid_primer_input(primer_sequence)) {
    stop("Primer must contain 15–80 IUPAC DNA bases.")
  }

  query <- if (direction == "reverse") {
    reverse_complement(primer_sequence)
  } else {
    primer_sequence
  }
  query_bases <- split_bases(query)
  width <- length(query_bases)
  if (width > nchar(reference_sequence)) {
    stop("Primer is longer than the reference sequence.")
  }

  starts <- seq_len(nchar(reference_sequence) - width + 1L)
  rows <- lapply(starts, function(start) {
    template_bases <- split_bases(substr(
      reference_sequence, start, start + width - 1L
    ))
    mismatch <- !mapply(bases_match, query_bases, template_bases)
    compatible <- !mismatch
    indices <- which(mismatch)
    distances <- if (direction == "forward") {
      width - indices + 1L
    } else {
      indices
    }
    penalty <- if (length(indices)) {
      sum(
        position_penalty(distances) *
          mapply(
            mismatch_type_multiplier,
            query_bases[indices],
            template_bases[indices]
          )
      )
    } else {
      0
    }
    data.frame(
      start = start,
      end = start + width - 1L,
      mismatch_count = length(indices),
      mismatch_fraction = length(indices) / width,
      compatible_identity = mean(compatible),
      terminal_3_mismatches = sum(distances <= 3L),
      terminal_5_mismatches = sum(distances <= 5L),
      longest_compatible_run = longest_true_run(compatible),
      random_match_evalue = random_match_evalue(
        query_bases, sum(compatible), nchar(reference_sequence)
      ),
      localization_penalty = round(penalty, 3),
      reference_window = paste(template_bases, collapse = ""),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  result[order(
    result$mismatch_count,
    result$terminal_5_mismatches,
    result$localization_penalty,
    result$start
  ), , drop = FALSE]
}

locate_primer_pair <- function(
  reference_sequence,
  forward_sequence,
  reverse_sequence,
  expected_amplicon_bp = NA_real_,
  candidate_mismatch_slack = 2L,
  max_candidates_per_primer = 100L,
  alignment_templates = NULL,
  max_alignment_rescue_candidates = 10L,
  minimum_alignment_supported_templates = 2L,
  minimum_alignment_support_fraction = 0.10
) {
  expected_supplied <- is.finite(expected_amplicon_bp) && expected_amplicon_bp > 0
  alignment_minimum_identity <- if (expected_supplied) 0.70 else 0.80
  forward <- primer_candidate_table(
    reference_sequence, forward_sequence, "forward"
  )
  reverse <- primer_candidate_table(
    reference_sequence, reverse_sequence, "reverse"
  )

  credible_sites <- function(candidate_table, primer_sequence) {
    vapply(
      seq_len(nrow(candidate_table)),
      function(i) {
        assess_primer_binding(
          primer_sequence, candidate_table[i, , drop = FALSE]
        )$credible
      },
      logical(1)
    )
  }
  forward$credible_site <- credible_sites(forward, forward_sequence)
  reverse$credible_site <- credible_sites(reverse, reverse_sequence)
  forward$evidence_basis <- ifelse(
    forward$credible_site, "coordinate reference", "rejected"
  )
  reverse$evidence_basis <- ifelse(
    reverse$credible_site, "coordinate reference", "rejected"
  )

  rescue_from_alignments <- function(
    candidate_table,
    primer_sequence,
    direction
  ) {
    if (any(candidate_table$credible_site) || !length(alignment_templates)) {
      return(candidate_table)
    }
    candidates_to_check <- seq_len(min(
      nrow(candidate_table), max_alignment_rescue_candidates
    ))
    for (i in candidates_to_check) {
      support <- alignment_binding_support(
        alignment_templates = alignment_templates,
        primer_sequence = primer_sequence,
        start = candidate_table$start[i],
        end = candidate_table$end[i],
        direction = direction,
        minimum_identity = alignment_minimum_identity
      )
      if (alignment_support_passes(
        support$summary,
        minimum_supported_templates = minimum_alignment_supported_templates,
        minimum_support_fraction = minimum_alignment_support_fraction
      )) {
        candidate_table$credible_site[i] <- TRUE
        candidate_table$evidence_basis[i] <- "installed COI alignments"
        break
      }
    }
    candidate_table
  }
  forward <- rescue_from_alignments(
    forward, forward_sequence, "forward"
  )
  reverse <- rescue_from_alignments(
    reverse, reverse_sequence, "reverse"
  )

  describe_candidate <- function(direction, candidate, sequence, assessment) {
    paste0(
      tools::toTitleCase(direction), " best accidental site ",
      candidate$start, "–", candidate$end, ": ",
      candidate$mismatch_count, "/", nchar(normalize_primer_input(sequence)),
      " incompatible bases (", round(candidate$compatible_identity * 100, 1),
      "% compatible identity), terminal-3 mismatches ",
      candidate$terminal_3_mismatches,
      ", longest compatible run ", candidate$longest_compatible_run,
      " bp, chance-match expectation ",
      format(candidate$random_match_evalue, digits = 2, scientific = TRUE),
      ". ",
      if (assessment$credible) {
        "This individual primer passed the sequence-specificity gate."
      } else {
        paste0("Rejected because ", assessment$summary, ".")
      }
    )
  }
  if (!any(forward$credible_site) || !any(reverse$credible_site)) {
    forward_assessment <- assess_primer_binding(
      forward_sequence, forward[1, , drop = FALSE]
    )
    reverse_assessment <- assess_primer_binding(
      reverse_sequence, reverse[1, , drop = FALSE]
    )
    alignment_note <- if (length(alignment_templates)) {
      paste0(
        " None of the first ",
        min(max_alignment_rescue_candidates, nrow(forward), nrow(reverse)),
        " reference-ranked candidate sites was supported by at least ",
        minimum_alignment_supported_templates,
        " templates and ",
        round(minimum_alignment_support_fraction * 100),
        "% of available templates in an installed order alignment."
      )
    } else {
      ""
    }
    stop(
      paste(
        "No credible COI binding pair was found.",
        describe_candidate(
          "forward", forward[1, , drop = FALSE],
          forward_sequence, forward_assessment
        ),
        describe_candidate(
          "reverse", reverse[1, , drop = FALSE],
          reverse_sequence, reverse_assessment
        ),
        paste(
          "The oligos may target another marker, be entered in the wrong",
          "5′→3′ orientation, or require a different COI reference."
        ),
        alignment_note
      ),
      call. = FALSE
    )
  }

  credible_forward <- forward[forward$credible_site, , drop = FALSE]
  credible_reverse <- reverse[reverse$credible_site, , drop = FALSE]
  forward_pool <- credible_forward[
    credible_forward$mismatch_count <=
      min(credible_forward$mismatch_count) + candidate_mismatch_slack,
    ,
    drop = FALSE
  ]
  reverse_pool <- credible_reverse[
    credible_reverse$mismatch_count <=
      min(credible_reverse$mismatch_count) + candidate_mismatch_slack,
    ,
    drop = FALSE
  ]
  forward_pool <- head(forward_pool, max_candidates_per_primer)
  reverse_pool <- head(reverse_pool, max_candidates_per_primer)

  combinations <- merge(
    transform(forward_pool, candidate_forward = seq_len(nrow(forward_pool))),
    transform(reverse_pool, candidate_reverse = seq_len(nrow(reverse_pool))),
    by = NULL,
    suffixes = c("_forward", "_reverse")
  )
  combinations <- combinations[
    combinations$start_forward < combinations$start_reverse,
    ,
    drop = FALSE
  ]
  if (!nrow(combinations)) {
    stop("No correctly oriented forward/reverse site combination was found.")
  }

  combinations$aligned_amplicon_bp <-
    combinations$end_reverse - combinations$start_forward + 1L
  combinations$targeted_region_bp <- pmax(
    0L,
    combinations$start_reverse - combinations$end_forward - 1L
  )
  combinations$total_mismatches <-
    combinations$mismatch_count_forward + combinations$mismatch_count_reverse
  combinations$total_terminal_5_mismatches <-
    combinations$terminal_5_mismatches_forward +
    combinations$terminal_5_mismatches_reverse
  combinations$total_localization_penalty <-
    combinations$localization_penalty_forward +
    combinations$localization_penalty_reverse
  combinations$expected_length_difference <- if (expected_supplied) {
    abs(combinations$aligned_amplicon_bp - expected_amplicon_bp)
  } else {
    NA_real_
  }

  ranking <- if (expected_supplied) {
    order(
      combinations$expected_length_difference,
      combinations$total_mismatches,
      combinations$total_terminal_5_mismatches,
      combinations$total_localization_penalty,
      combinations$aligned_amplicon_bp
    )
  } else {
    order(
      combinations$total_mismatches,
      combinations$total_terminal_5_mismatches,
      combinations$total_localization_penalty,
      combinations$aligned_amplicon_bp
    )
  }
  combinations <- combinations[ranking, , drop = FALSE]
  combinations$pair_candidate_rank <- seq_len(nrow(combinations))

  best <- combinations[1, , drop = FALSE]
  selected_candidate <- function(direction) {
    data.frame(
      compatible_identity = best[[paste0("compatible_identity_", direction)]],
      terminal_3_mismatches =
        best[[paste0("terminal_3_mismatches_", direction)]],
      longest_compatible_run =
        best[[paste0("longest_compatible_run_", direction)]],
      random_match_evalue =
        best[[paste0("random_match_evalue_", direction)]]
    )
  }
  forward_assessment <- assess_primer_binding(
    forward_sequence, selected_candidate("forward")
  )
  reverse_assessment <- assess_primer_binding(
    reverse_sequence, selected_candidate("reverse")
  )
  forward_order_support <- alignment_binding_support(
    alignment_templates,
    forward_sequence,
    best$start_forward,
    best$end_forward,
    "forward",
    minimum_identity = alignment_minimum_identity
  )
  reverse_order_support <- alignment_binding_support(
    alignment_templates,
    reverse_sequence,
    best$start_reverse,
    best$end_reverse,
    "reverse",
    minimum_identity = alignment_minimum_identity
  )
  pair_order_support <- combine_pair_alignment_support(
    forward_order_support, reverse_order_support
  )
  alignment_rescue_used <-
    best$evidence_basis_forward == "installed COI alignments" ||
    best$evidence_basis_reverse == "installed COI alignments"
  if (alignment_rescue_used && !alignment_support_passes(
    transform(
      pair_order_support,
      n_supported = n_pair_supported,
      support_fraction = pair_support_fraction
    ),
    minimum_supported_templates = minimum_alignment_supported_templates,
    minimum_support_fraction = minimum_alignment_support_fraction
  )) {
    stop(
      paste(
        "The individual oligos had plausible COI sites, but no installed order",
        "alignment supported them together as a primer pair at those coordinates."
      ),
      call. = FALSE
    )
  }
  equal_best <- if (expected_supplied) {
    combinations[
      combinations$expected_length_difference == best$expected_length_difference &
        combinations$total_mismatches == best$total_mismatches &
        combinations$total_terminal_5_mismatches == best$total_terminal_5_mismatches &
        combinations$total_localization_penalty == best$total_localization_penalty,
      ,
      drop = FALSE
    ]
  } else {
    combinations[
      combinations$total_mismatches == best$total_mismatches &
        combinations$total_terminal_5_mismatches == best$total_terminal_5_mismatches &
        combinations$total_localization_penalty == best$total_localization_penalty,
      ,
      drop = FALSE
    ]
  }

  list(
    best = best,
    forward_assessment = forward_assessment,
    reverse_assessment = reverse_assessment,
    forward_order_support = forward_order_support,
    reverse_order_support = reverse_order_support,
    pair_order_support = pair_order_support,
    alignment_rescue_used = alignment_rescue_used,
    forward_candidates = head(forward, 5L),
    reverse_candidates = head(reverse, 5L),
    pair_candidates = head(combinations, 5L),
    forward_equally_best_count = sum(
      forward$mismatch_count == forward$mismatch_count[1] &
        forward$terminal_5_mismatches == forward$terminal_5_mismatches[1] &
        forward$localization_penalty == forward$localization_penalty[1]
    ),
    reverse_equally_best_count = sum(
      reverse$mismatch_count == reverse$mismatch_count[1] &
        reverse$terminal_5_mismatches == reverse$terminal_5_mismatches[1] &
        reverse$localization_penalty == reverse$localization_penalty[1]
    ),
    equally_best_pair_count = nrow(equal_best),
    expected_amplicon_supplied = expected_supplied,
    alignment_minimum_identity = alignment_minimum_identity
  )
}

parse_fasta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  headers <- which(startsWith(lines, ">"))
  ends <- c(headers[-1] - 1, length(lines))
  sequences <- vapply(
    seq_along(headers),
    function(i) paste(lines[(headers[i] + 1):ends[i]], collapse = ""),
    character(1)
  )
  names(sequences) <- sub("^>", "", lines[headers])
  sequences
}

primer_map_domain <- function(
  pairs,
  folmer_start,
  folmer_end,
  tick_bp = 100,
  reference_length = 1536
) {
  stopifnot(
    is.data.frame(pairs),
    is.numeric(folmer_start),
    is.numeric(folmer_end),
    length(folmer_start) == 1L,
    length(folmer_end) == 1L,
    is.finite(tick_bp),
    tick_bp > 0,
    is.finite(reference_length),
    reference_length > 0
  )
  coordinate_columns <- intersect(
    c("forward_start", "forward_end", "reverse_start", "reverse_end"),
    names(pairs)
  )
  selected_coordinates <- if (length(coordinate_columns) && nrow(pairs)) {
    unlist(pairs[, coordinate_columns, drop = FALSE], use.names = FALSE)
  } else {
    numeric()
  }
  selected_coordinates <- selected_coordinates[is.finite(selected_coordinates)]
  requested_min <- min(c(folmer_start, selected_coordinates), na.rm = TRUE)
  requested_max <- max(c(folmer_end, selected_coordinates), na.rm = TRUE)
  domain_min <- max(0, floor(requested_min / tick_bp) * tick_bp)
  domain_max <- min(
    ceiling(reference_length / tick_bp) * tick_bp,
    ceiling(requested_max / tick_bp) * tick_bp
  )
  if (domain_max <= domain_min) {
    domain_max <- min(
      ceiling(reference_length / tick_bp) * tick_bp,
      domain_min + tick_bp
    )
  }
  c(min = domain_min, max = domain_max)
}

sort_primer_map_rows <- function(pairs, mode = "folmer_overlap") {
  stopifnot(is.data.frame(pairs), length(mode) == 1L)
  if (!nrow(pairs)) return(pairs)
  required <- c(
    "pair_label", "pair_start", "pair_end",
    "folmer_overlap_bp", "aligned_amplicon_bp"
  )
  missing_columns <- setdiff(required, names(pairs))
  if (length(missing_columns)) {
    stop(
      "Primer-map sorting requires columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  order_index <- switch(
    mode,
    binding_site = order(
      pairs$pair_start,
      pairs$pair_end,
      tolower(pairs$pair_label),
      na.last = TRUE
    ),
    alphabetical = order(tolower(pairs$pair_label), na.last = TRUE),
    folmer_overlap = order(
      -pairs$folmer_overlap_bp,
      tolower(pairs$pair_label),
      na.last = TRUE
    ),
    amplicon_length = order(
      pairs$aligned_amplicon_bp,
      tolower(pairs$pair_label),
      na.last = TRUE
    ),
    stop("Unknown primer-map sort mode: ", mode, call. = FALSE)
  )
  pairs[order_index, , drop = FALSE]
}

number_primer_choice_labels <- function(pair_metadata, displayed_pairs) {
  stopifnot(
    is.data.frame(pair_metadata),
    is.data.frame(displayed_pairs),
    all(c("pair_id", "pair_label") %in% names(pair_metadata)),
    "pair_id" %in% names(displayed_pairs)
  )
  displayed_index <- match(
    as.character(pair_metadata$pair_id),
    as.character(displayed_pairs$pair_id)
  )
  labels <- as.character(pair_metadata$pair_label)
  numbered <- !is.na(displayed_index)
  labels[numbered] <- paste0(
    displayed_index[numbered],
    " · ",
    labels[numbered]
  )
  stats::setNames(as.character(pair_metadata$pair_id), labels)
}

primer_use_tags <- function(use_case, target = NA_character_) {
  use_case <- ifelse(is.na(use_case), "", trimws(as.character(use_case)))
  target <- ifelse(is.na(target), "", trimws(as.character(target)))
  search_text <- tolower(paste(use_case, target))
  tags <- character()
  add_tag <- function(condition, label) {
    if (isTRUE(condition)) tags <<- c(tags, label)
  }
  add_tag(grepl("edna", search_text, fixed = TRUE), "eDNA metabarcoding")
  add_tag(
    grepl("diet", search_text, fixed = TRUE) ||
      grepl("gut content", search_text, fixed = TRUE),
    "diet metabarcoding"
  )
  add_tag(grepl("bulk", search_text, fixed = TRUE), "bulk metabarcoding")
  add_tag(
    grepl("freshwater", search_text, fixed = TRUE),
    "freshwater metabarcoding"
  )
  add_tag(
    grepl("broad metabarcoding", search_text, fixed = TRUE),
    "broad metabarcoding"
  )
  add_tag(
    grepl("regular barcoding", search_text, fixed = TRUE) ||
      grepl("fish barcoding", search_text, fixed = TRUE),
    "DNA barcoding"
  )
  add_tag(
    grepl("degraded dna", search_text, fixed = TRUE),
    "degraded-DNA assay"
  )
  add_tag(
    grepl("plant-associated", search_text, fixed = TRUE),
    "plant-associated arthropods"
  )
  add_tag(
    grepl("custom input", search_text, fixed = TRUE),
    "custom primer"
  )
  if (!length(tags) && nzchar(use_case)) {
    tags <- use_case
  }
  if (nzchar(target)) {
    tags <- c(tags, paste0("target: ", target))
  }
  unique(tags)
}

score_primer_pair_primerminer <- function(
  alignment_files,
  primer_rows,
  geometry
) {
  if (!requireNamespace("PrimerMiner", quietly = TRUE)) {
    stop(
      "PrimerMiner is required to score custom pairs. Run scripts/install_primerminer.R.",
      call. = FALSE
    )
  }
  suppressPackageStartupMessages(
    require("PrimerMiner", character.only = TRUE, quietly = TRUE)
  )
  stopifnot(
    length(alignment_files) > 0L,
    nrow(primer_rows) == 2L,
    nrow(geometry) == 1L,
    all(c("forward", "reverse") %in% primer_rows$direction)
  )

  median_safe <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) median(x) else NA_real_
  }
  quantile_safe <- function(x, probability) {
    x <- x[is.finite(x)]
    if (length(x)) {
      unname(quantile(x, probability, na.rm = TRUE))
    } else {
      NA_real_
    }
  }
  order_from_path <- function(path) {
    sub("_COX1_reference_aligned\\.fasta$", "", basename(path))
  }

  primer_template_rows <- list()
  position_rows <- list()
  template_index <- 1L
  position_index <- 1L

  for (alignment_path in alignment_files) {
    order_name <- order_from_path(alignment_path)
    for (i in seq_len(nrow(primer_rows))) {
      primer <- primer_rows[i, , drop = FALSE]
      direction <- primer$direction[1]
      alignment_start <- if (direction == "forward") {
        geometry$forward_start[1]
      } else {
        geometry$reverse_start[1]
      }
      alignment_end <- if (direction == "forward") {
        geometry$forward_end[1]
      } else {
        geometry$reverse_end[1]
      }

      evaluated <- suppressMessages(PrimerMiner::evaluate_primer(
        alignment_imp = alignment_path,
        primer_sequ = clean_sequence(primer$sequence[1]),
        start = alignment_start,
        stop = alignment_end,
        forward = direction == "forward",
        gap_NA = TRUE,
        N_NA = TRUE,
        mm_position = "Position_v1",
        mm_type = "Type_v1",
        adjacent = 2,
        sequ_names = TRUE
      ))
      evaluated <- evaluated[
        !grepl("^NC_001322\\.1_COX1_reference", evaluated$Template),
        ,
        drop = FALSE
      ]

      score_columns <- grep("^V[0-9]+$", names(evaluated), value = TRUE)
      score_matrix <- as.matrix(evaluated[, score_columns, drop = FALSE])
      mismatch_matrix <- score_matrix > 0
      mismatch_matrix[is.na(mismatch_matrix)] <- FALSE
      missing_matrix <- is.na(score_matrix)
      adjacent_count <- apply(mismatch_matrix, 1, function(x) {
        if (length(x) < 2L) return(0L)
        sum(x[-length(x)] & x[-1L])
      })
      distance_values <- as.integer(sub("^V", "", score_columns))
      terminal_count <- function(max_distance) {
        columns <- which(distance_values <= max_distance)
        if (!length(columns)) return(rep(0L, nrow(evaluated)))
        rowSums(mismatch_matrix[, columns, drop = FALSE])
      }

      primer_template_rows[[template_index]] <- data.frame(
        order = order_name,
        pair_id = primer$pair_id[1],
        pair_label = primer$pair_label[1],
        primer_name = primer$primer_name[1],
        direction = direction,
        template = as.character(evaluated$Template),
        binding_sequence_primer_oriented = evaluated$sequ,
        alignment_start = alignment_start,
        alignment_end = alignment_end,
        primerminer_penalty = evaluated$sum,
        scorable = !is.na(evaluated$sum),
        missing_binding_positions = rowSums(missing_matrix),
        mismatch_count = rowSums(mismatch_matrix),
        terminal_1_mismatches = terminal_count(1L),
        terminal_3_mismatches = terminal_count(3L),
        terminal_5_mismatches = terminal_count(5L),
        adjacent_mismatch_pairs = adjacent_count,
        stringsAsFactors = FALSE
      )
      template_index <- template_index + 1L

      primer_bases <- split_bases(primer$sequence[1])
      template_bases <- lapply(evaluated$sequ, split_bases)
      for (column_index in seq_along(score_columns)) {
        distance <- distance_values[column_index]
        primer_index <- length(primer_bases) - distance + 1L
        base_values <- vapply(
          template_bases,
          function(x) {
            if (length(x) >= primer_index) x[primer_index] else NA_character_
          },
          character(1)
        )
        position_rows[[position_index]] <- data.frame(
          order = order_name,
          pair_id = primer$pair_id[1],
          pair_label = primer$pair_label[1],
          primer_name = primer$primer_name[1],
          direction = direction,
          template = as.character(evaluated$Template),
          alignment_position = if (direction == "forward") {
            alignment_start + primer_index - 1L
          } else {
            alignment_end - primer_index + 1L
          },
          primer_position_5_to_3 = primer_index,
          distance_from_3prime = distance,
          primer_base = primer_bases[primer_index],
          template_base = base_values,
          position_penalty = score_matrix[, column_index],
          mismatch = mismatch_matrix[, column_index],
          unavailable_reason = ifelse(
            !missing_matrix[, column_index],
            "",
            ifelse(
              base_values == "-",
              "alignment gap",
              ifelse(base_values == "N", "ambiguous N", "unscorable")
            )
          ),
          stringsAsFactors = FALSE
        )
        position_index <- position_index + 1L
      }
    }
  }

  primer_template_scores <- do.call(rbind, primer_template_rows)
  primer_position_scores <- do.call(rbind, position_rows)
  forward <- primer_template_scores[
    primer_template_scores$direction == "forward",
    ,
    drop = FALSE
  ]
  reverse <- primer_template_scores[
    primer_template_scores$direction == "reverse",
    ,
    drop = FALSE
  ]
  join_keys <- c("order", "pair_id", "pair_label", "template")
  names(forward)[!names(forward) %in% join_keys] <- paste0(
    "forward_",
    names(forward)[!names(forward) %in% join_keys]
  )
  names(reverse)[!names(reverse) %in% join_keys] <- paste0(
    "reverse_",
    names(reverse)[!names(reverse) %in% join_keys]
  )
  pair_template_scores <- merge(
    forward,
    reverse,
    by = join_keys,
    all = TRUE,
    sort = FALSE
  )
  pair_template_scores$pair_scorable <-
    pair_template_scores$forward_scorable %in% TRUE &
      pair_template_scores$reverse_scorable %in% TRUE
  pair_template_scores$pair_penalty <- ifelse(
    pair_template_scores$pair_scorable,
    pair_template_scores$forward_primerminer_penalty +
      pair_template_scores$reverse_primerminer_penalty,
    NA_real_
  )
  pair_template_scores$perfect_match_pair <-
    pair_template_scores$pair_scorable &
      pair_template_scores$forward_mismatch_count == 0L &
      pair_template_scores$reverse_mismatch_count == 0L

  groups <- split(
    pair_template_scores,
    interaction(
      pair_template_scores$order,
      pair_template_scores$pair_id,
      drop = TRUE
    )
  )
  order_pair_summary <- do.call(rbind, lapply(groups, function(x) {
    scorable <- x$pair_scorable %in% TRUE
    data.frame(
      order = x$order[1],
      pair_id = x$pair_id[1],
      pair_label = x$pair_label[1],
      n_templates = nrow(x),
      n_pair_scorable = sum(scorable),
      pair_scorable_fraction = mean(scorable),
      unavailable_fraction = mean(!scorable),
      perfect_match_fraction_among_scorable = if (sum(scorable)) {
        mean(x$perfect_match_pair[scorable])
      } else {
        NA_real_
      },
      forward_median_penalty =
        median_safe(x$forward_primerminer_penalty),
      forward_p90_penalty =
        quantile_safe(x$forward_primerminer_penalty, 0.9),
      reverse_median_penalty =
        median_safe(x$reverse_primerminer_penalty),
      reverse_p90_penalty =
        quantile_safe(x$reverse_primerminer_penalty, 0.9),
      pair_median_penalty = median_safe(x$pair_penalty),
      pair_p90_penalty = quantile_safe(x$pair_penalty, 0.9),
      terminal_1_mismatch_fraction = if (sum(scorable)) {
        mean(
          x$forward_terminal_1_mismatches[scorable] +
            x$reverse_terminal_1_mismatches[scorable] > 0
        )
      } else {
        NA_real_
      },
      terminal_3_mismatch_fraction = if (sum(scorable)) {
        mean(
          x$forward_terminal_3_mismatches[scorable] +
            x$reverse_terminal_3_mismatches[scorable] > 0
        )
      } else {
        NA_real_
      },
      adjacent_mismatch_fraction = if (sum(scorable)) {
        mean(
          x$forward_adjacent_mismatch_pairs[scorable] +
            x$reverse_adjacent_mismatch_pairs[scorable] > 0
        )
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }))
  rownames(order_pair_summary) <- NULL

  list(
    order_scores = order_pair_summary[
      order(
        order_pair_summary$order,
        order_pair_summary$pair_label
      ),
      ,
      drop = FALSE
    ],
    pair_template_scores = pair_template_scores,
    primer_position_scores = primer_position_scores
  )
}

summarize_expanded_primer_scores <- function(data, group_columns) {
  median_safe <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) median(x) else NA_real_
  }
  mean_safe <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) mean(x) else NA_real_
  }
  quantile_safe <- function(x, probability) {
    x <- x[is.finite(x)]
    if (length(x)) {
      unname(quantile(x, probability, na.rm = TRUE))
    } else {
      NA_real_
    }
  }

  usable <- data
  for (column in group_columns) {
    usable[[column]][
      is.na(usable[[column]]) | !nzchar(trimws(usable[[column]]))
    ] <- paste0(column, " unresolved")
  }
  if (!nrow(usable)) return(data.frame())
  group_key <- interaction(
    usable[, group_columns, drop = FALSE],
    drop = TRUE,
    lex.order = TRUE
  )
  groups <- split(usable, group_key)
  result <- do.call(rbind, lapply(groups, function(x) {
    scorable <- x$pair_scorable %in% TRUE
    group_values <- x[1, group_columns, drop = FALSE]
    metrics <- data.frame(
      n_centroids = nrow(x),
      source_records_represented = sum(x$cluster_size, na.rm = TRUE),
      n_forward_scorable = sum(x$forward_scorable %in% TRUE),
      forward_scorable_fraction = mean(x$forward_scorable %in% TRUE),
      n_reverse_scorable = sum(x$reverse_scorable %in% TRUE),
      reverse_scorable_fraction = mean(x$reverse_scorable %in% TRUE),
      n_pair_scorable = sum(scorable),
      pair_scorable_fraction = mean(scorable),
      forward_median_penalty = median_safe(x$forward_penalty),
      reverse_median_penalty = median_safe(x$reverse_penalty),
      pair_median_penalty = median_safe(x$pair_penalty),
      pair_mean_penalty = mean_safe(x$pair_penalty),
      pair_p90_penalty = quantile_safe(x$pair_penalty, 0.9),
      perfect_match_fraction_among_scorable = if (sum(scorable)) {
        mean(x$perfect_match_pair[scorable])
      } else {
        NA_real_
      },
      terminal_3_mismatch_fraction_among_scorable = if (sum(scorable)) {
        mean(x$terminal_3_mismatch_pair[scorable])
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
    cbind(group_values, metrics)
  }))
  rownames(result) <- NULL
  result
}

score_primer_pair_expanded_reference <- function(
  alignment_path,
  taxonomy,
  primer_rows,
  geometry,
  alignment_offset
) {
  if (!requireNamespace("PrimerMiner", quietly = TRUE)) {
    stop(
      "PrimerMiner is required for expanded custom-primer scoring.",
      call. = FALSE
    )
  }
  suppressPackageStartupMessages(
    require("PrimerMiner", character.only = TRUE, quietly = TRUE)
  )
  stopifnot(
    file.exists(alignment_path),
    nrow(primer_rows) == 2L,
    nrow(geometry) == 1L,
    all(c("forward", "reverse") %in% primer_rows$direction),
    "accession" %in% names(taxonomy)
  )

  first_fasta_width <- function(path) {
    connection <- file(path, open = "r")
    on.exit(close(connection), add = TRUE)
    found_header <- FALSE
    width <- 0L
    repeat {
      line <- readLines(connection, n = 1L, warn = FALSE)
      if (!length(line)) break
      if (startsWith(line, ">")) {
        if (found_header) break
        found_header <- TRUE
      } else if (found_header) {
        width <- width + nchar(gsub("[[:space:]]", "", line))
      }
    }
    if (!found_header || width < 1L) {
      stop("Expanded-reference FASTA has no readable aligned sequence.")
    }
    width
  }
  alignment_width <- first_fasta_width(alignment_path)
  direction_limitations <- character()

  unavailable_direction <- function(primer_length) {
    data.frame(
      accession = as.character(taxonomy$accession),
      binding_sequence_primer_oriented = NA_character_,
      penalty = NA_real_,
      scorable = FALSE,
      missing_binding_positions = as.integer(primer_length),
      mismatch_count = NA_integer_,
      terminal_1_mismatches = NA_integer_,
      terminal_3_mismatches = NA_integer_,
      terminal_5_mismatches = NA_integer_,
      stringsAsFactors = FALSE
    )
  }

  evaluate_direction <- function(direction) {
    primer <- primer_rows[
      primer_rows$direction == direction,
      ,
      drop = FALSE
    ]
    start <- if (direction == "forward") {
      geometry$forward_start[1]
    } else {
      geometry$reverse_start[1]
    }
    end <- if (direction == "forward") {
      geometry$forward_end[1]
    } else {
      geometry$reverse_end[1]
    }
    start <- as.integer(start + alignment_offset)
    end <- as.integer(end + alignment_offset)
    if (end - start + 1L != nchar(clean_sequence(primer$sequence[1]))) {
      stop("Expanded-reference coordinates do not match the primer length.")
    }
    if (start < 1L || end > alignment_width) {
      direction_limitations[direction] <<- paste0(
        tools::toTitleCase(direction), " primer requests alignment columns ",
        start, "–", end, ", outside the expanded reference span 1–",
        alignment_width, "."
      )
      return(unavailable_direction(nchar(clean_sequence(primer$sequence[1]))))
    }

    evaluated <- suppressMessages(PrimerMiner::evaluate_primer(
      alignment_imp = alignment_path,
      primer_sequ = clean_sequence(primer$sequence[1]),
      start = start,
      stop = end,
      forward = direction == "forward",
      gap_NA = TRUE,
      N_NA = TRUE,
      mm_position = "Position_v1",
      mm_type = "Type_v1",
      adjacent = 2,
      sequ_names = TRUE
    ))
    score_columns <- grep("^V[0-9]+$", names(evaluated), value = TRUE)
    score_matrix <- as.matrix(evaluated[, score_columns, drop = FALSE])
    mismatch_matrix <- score_matrix > 0
    mismatch_matrix[is.na(mismatch_matrix)] <- FALSE
    missing_matrix <- is.na(score_matrix)
    distances <- as.integer(sub("^V", "", score_columns))
    terminal_count <- function(max_distance) {
      columns <- which(distances <= max_distance)
      if (!length(columns)) return(rep(0L, nrow(evaluated)))
      rowSums(mismatch_matrix[, columns, drop = FALSE])
    }

    data.frame(
      accession = as.character(evaluated$Template),
      binding_sequence_primer_oriented = evaluated$sequ,
      penalty = evaluated$sum,
      scorable = !is.na(evaluated$sum),
      missing_binding_positions = rowSums(missing_matrix),
      mismatch_count = rowSums(mismatch_matrix),
      terminal_1_mismatches = terminal_count(1L),
      terminal_3_mismatches = terminal_count(3L),
      terminal_5_mismatches = terminal_count(5L),
      stringsAsFactors = FALSE
    )
  }

  forward <- evaluate_direction("forward")
  reverse <- evaluate_direction("reverse")
  names(forward)[-1] <- paste0("forward_", names(forward)[-1])
  names(reverse)[-1] <- paste0("reverse_", names(reverse)[-1])
  scores <- merge(forward, reverse, by = "accession", sort = FALSE)
  taxonomy_match <- match(scores$accession, taxonomy$accession)
  if (anyNA(taxonomy_match)) {
    stop("Expanded-reference scores could not be mapped to taxonomy.")
  }
  scores <- cbind(
    data.frame(
      pair_id = primer_rows$pair_id[1],
      pair_label = primer_rows$pair_label[1],
      stringsAsFactors = FALSE
    ),
    scores,
    taxonomy[
      taxonomy_match,
      setdiff(names(taxonomy), "accession"),
      drop = FALSE
    ]
  )
  scores$pair_scorable <-
    scores$forward_scorable %in% TRUE &
    scores$reverse_scorable %in% TRUE
  scores$pair_penalty <- ifelse(
    scores$pair_scorable,
    scores$forward_penalty + scores$reverse_penalty,
    NA_real_
  )
  scores$perfect_match_pair <-
    scores$pair_scorable &
    scores$forward_mismatch_count == 0L &
    scores$reverse_mismatch_count == 0L
  scores$terminal_3_mismatch_pair <-
    scores$pair_scorable &
    (
      scores$forward_terminal_3_mismatches +
        scores$reverse_terminal_3_mismatches
    ) > 0L

  scorable_fraction <- mean(scores$pair_scorable)
  suitable <-
    sum(scores$pair_scorable) >= 100L &&
    scorable_fraction >= 0.50
  scores$reference_suitable <- suitable
  scores$reference_suitability <- if (suitable) {
    "usable"
  } else {
    "binding-region incomplete"
  }
  scores$reference_limitation_reason <- if (suitable) {
    ""
  } else if (length(direction_limitations)) {
    paste(
      paste(unname(direction_limitations), collapse = " "),
      paste0(
        "The expanded centroid alignment does not span both primer sites, so ",
        "pair penalties and taxon-coverage interpretations are unavailable."
      )
    )
  } else {
    paste0(
      "Only ", sum(scores$pair_scorable), "/", nrow(scores),
      " centroids (", round(scorable_fraction * 100, 1),
      "%) span both primer sites; penalties cannot be interpreted as taxon coverage."
    )
  }

  overall <- summarize_expanded_primer_scores(
    scores,
    c("pair_id", "pair_label")
  )
  overall$reference_suitable <- suitable
  overall$reference_suitability <- unique(scores$reference_suitability)
  overall$reference_limitation_reason <- unique(
    scores$reference_limitation_reason
  )
  list(scores = scores, overall = overall)
}

alignment_consensus <- function(sequences) {
  chars <- strsplit(unname(sequences), "", fixed = TRUE)
  width <- max(lengths(chars))
  matrix <- do.call(rbind, lapply(chars, function(x) {
    length(x) <- width
    x
  }))
  consensus <- apply(matrix, 2, function(column) {
    usable <- column[column %in% c("A", "C", "G", "T")]
    if (!length(usable)) return("N")
    names(which.max(table(usable)))
  })
  paste(consensus, collapse = "")
}

pair_thermal_summary <- function(primers, sodium_mM = 50, offset = 4) {
  rows <- lapply(seq_len(nrow(primers)), function(i) {
    tm <- tm_range(primers$sequence[i], sodium_mM)
    data.frame(
      pair_id = primers$pair_id[i],
      pair_label = primers$pair_label[i],
      primer_name = primers$primer_name[i],
      direction = primers$direction[i],
      sequence = primers$sequence[i],
      length = nchar(clean_sequence(primers$sequence[i])),
      degeneracy = degeneracy(primers$sequence[i]),
      tm_min = tm[1],
      tm_max = tm[2],
      reported_ta_c = primers$reported_ta_c[i],
      source_key = primers$source_key[i],
      stringsAsFactors = FALSE
    )
  })
  table <- do.call(rbind, rows)
  lower_tm <- aggregate(tm_min ~ pair_id, table, min)
  lower_tm$ta_low <- lower_tm$tm_min - offset - 1
  lower_tm$ta_high <- lower_tm$tm_min - offset + 1
  merge(table, lower_tm[, c("pair_id", "ta_low", "ta_high")], by = "pair_id")
}
