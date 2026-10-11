# Pure views of a frozen alignment. Mapping and validation live in
# R/fungal_alignment.R; alignment columns never substitute for raw product bp.

.fungal_alignment_ui_tier <- function(bp, ceilings) {
  limits <- as.numeric(ceilings)
  if (length(limits) != 3L || any(!is.finite(limits)) || any(diff(limits) < 0)) {
    stop("Alignment views require three ordered platform length ceilings")
  }
  if (!length(bp) || !is.finite(bp)) return(list(label = "Unresolved", color = "#64748b"))
  tier <- if (bp <= limits[1]) 1L else if (bp <= limits[2]) 2L else if (bp <= limits[3]) 3L else 4L
  list(label = c("Illumina", "PacBio", "Nanopore", "Beyond configured limits")[tier],
       color = c("#2563eb", "#7c3aed", "#0f766e", "#64748b")[tier])
}

.fungal_alignment_ui_status <- function(row) {
  if (!nrow(row)) return("No pair evidence on this template")
  value <- as.character(row$alignment_status[1])
  if (is.na(value) || !nzchar(value)) value <- as.character(row$status[1])
  gsub("_", " ", value, fixed = TRUE)
}

.fungal_alignment_ui_row <- function(rows, accession, pair_id) {
  if (!nrow(rows)) return(rows)
  rows[rows$accession == accession & rows$pair_id == pair_id, , drop = FALSE]
}

.fungal_alignment_ui_taxon <- function(taxonomy, accession) {
  row <- taxonomy[match(accession, taxonomy$accession), , drop = FALSE]
  values <- vapply(intersect(c("organism", "class"), names(row)), function(name) {
    value <- as.character(row[[name]][1])
    if (is.na(value) || !nzchar(value) || value == "unresolved") "" else value
  }, character(1))
  paste(values[nzchar(values)], collapse = " · ")
}

.fungal_alignment_ui_domain <- function(index, domain) {
  if (is.null(domain) || length(domain) != 2L || any(!is.finite(domain))) {
    domain <- c(1L, index$n_columns)
  }
  domain <- sort(as.integer(domain))
  domain <- pmax(1L, pmin(index$n_columns, domain))
  if (domain[1] == domain[2] && index$n_columns > 1L) {
    domain <- c(max(1L, domain[1] - 1L), min(index$n_columns, domain[2] + 1L))
  }
  domain
}

.fungal_alignment_ui_svg_text <- function(x, y, text, size = 12, color = "#475569", ...) {
  htmltools::tags$text(x = x, y = y, fill = color, `font-size` = size, text, ...)
}

.fungal_alignment_ui_axis <- function(domain, scale, bottom, top = 24) {
  ticks <- pretty(domain, n = 7)
  ticks <- ticks[ticks > domain[1] + diff(domain) * .055 & ticks < domain[2] - diff(domain) * .055]
  ticks <- unique(c(domain, ticks))
  lapply(ticks, function(column) htmltools::tagList(
    htmltools::tags$line(x1 = scale(column), x2 = scale(column), y1 = top + 8, y2 = bottom,
                        stroke = "#e2e8f0", `stroke-width` = 1),
    .fungal_alignment_ui_svg_text(scale(column), top, format(column, big.mark = ","),
                                 size = 11, `text-anchor` = "middle")
  ))
}

.fungal_alignment_ui_occupancy <- function(index, accession, domain, scale, y) {
  runs <- index$occupancy[index$occupancy$accession == accession, , drop = FALSE]
  runs <- runs[runs$end_col >= domain[1] & runs$start_col <= domain[2], , drop = FALSE]
  htmltools::tagList(
    htmltools::tags$line(x1 = scale(domain[1]), x2 = scale(domain[2]), y1 = y, y2 = y,
                        stroke = "#f1f5f9", `stroke-width` = 6),
    lapply(seq_len(nrow(runs)), function(i) htmltools::tags$line(
      x1 = scale(max(domain[1], runs$start_col[i])),
      x2 = scale(min(domain[2], runs$end_col[i])), y1 = y, y2 = y,
      stroke = "#94a3b8", `stroke-width` = 5,
      htmltools::tags$title(paste0("Non-gap alignment columns ", runs$start_col[i], "–", runs$end_col[i]))
    ))
  )
}

.fungal_alignment_ui_region_label <- function(row) {
  label <- as.character(row$label[1])
  if (grepl("ANCHOR", row$landmark_id[1], fixed = TRUE)) paste0(label, " anchor (partial)") else label
}

.fungal_alignment_ui_region_color <- function(row) {
  if (grepl("^ITS[12]$", row$label[1])) "#a16207" else "#475569"
}

.fungal_alignment_ui_region_title <- function(row) {
  paste0(.fungal_alignment_ui_region_label(row), " · ", row$accession,
         "\nRaw bases ", row$start, "–", row$end, " (", row$raw_bp, " bp)",
         "\nAlignment columns ", row$left_col, "–", row$right_col,
         "\n", if (row$evidence_status == "curated_reference") "Curated pilot boundaries" else "Frozen NCBI annotation",
         " · ", row$source_note)
}

.fungal_alignment_ui_visible_regions <- function(rows, accession, domain) {
  if (is.null(rows) || !nrow(rows)) return(NULL)
  rows[rows$accession %in% accession & rows$alignment_status == "mapped" &
       rows$right_col >= domain[1] & rows$left_col <= domain[2], , drop = FALSE]
}

# Gene labels sit above the track and spacer labels below it, so short genes
# remain legible without changing the actual width of their annotation span.
.fungal_alignment_ui_regions <- function(rows, accession, domain, scale, y) {
  rows <- .fungal_alignment_ui_visible_regions(rows, accession, domain)
  if (!length(rows) || !nrow(rows)) return(.fungal_alignment_ui_svg_text(
    scale(domain[1]), y + 4, "No frozen gene / region boundary in this window", size = 11))
  center <- (scale(pmax(domain[1], rows$left_col)) + scale(pmin(domain[2], rows$right_col))) / 2
  spacer <- grepl("^ITS[12]$", rows$label)
  label_x <- center
  for (lane in c(FALSE, TRUE)) {
    ids <- which(spacer == lane)
    ids <- ids[order(center[ids])]
    if (length(ids)) {
      for (k in seq_along(ids)) {
        label_x[ids[k]] <- max(scale(domain[1]) + 22, center[ids[k]],
                               if (k > 1L) label_x[ids[k - 1L]] + 88 else -Inf)
      }
      overflow <- max(0, tail(label_x[ids], 1) - (scale(domain[2]) - 22))
      label_x[ids] <- label_x[ids] - overflow
    }
  }
  htmltools::tagList(lapply(seq_len(nrow(rows)), function(i) {
    row <- rows[i, , drop = FALSE]
    left <- scale(max(domain[1], row$left_col))
    right <- scale(min(domain[2], row$right_col))
    color <- .fungal_alignment_ui_region_color(row)
    label_y <- y + if (spacer[i]) 29 else -15
    htmltools::tags$g(class = "fungal-alignment-region", `data-accession` = accession,
      `data-landmark-id` = row$landmark_id, `data-start-col` = row$left_col,
      `data-end-col` = row$right_col, `data-evidence` = row$evidence_status,
      htmltools::tags$title(.fungal_alignment_ui_region_title(row)),
      htmltools::tags$rect(x = left, y = y - 7, width = max(1, right - left), height = 14,
        rx = 2, fill = color, `fill-opacity` = .16, stroke = color,
        `stroke-dasharray` = if (row$evidence_status == "curated_reference") "3 2" else NULL),
      htmltools::tags$line(x1 = center[i], x2 = label_x[i], y1 = y + if (spacer[i]) 8 else -8,
                          y2 = label_y + if (spacer[i]) -10 else 3, stroke = color, `stroke-width` = .8),
      .fungal_alignment_ui_svg_text(label_x[i], label_y, .fungal_alignment_ui_region_label(row),
                                   size = 11, color = color, `text-anchor` = "middle", `font-weight` = 700)
    )
  }))
}

.fungal_alignment_ui_region_cells <- function(rows, accession, columns) {
  rows <- .fungal_alignment_ui_visible_regions(rows, accession, range(columns))
  labels <- rep("", length(columns))
  for (i in seq_along(columns)) {
    if (length(rows) && nrow(rows)) {
      at <- rows[rows$left_col <= columns[i] & rows$right_col >= columns[i], , drop = FALSE]
      if (nrow(at)) labels[i] <- paste(at$landmark_id, collapse = "|")
    }
  }
  runs <- rle(labels)
  ends <- cumsum(runs$lengths)
  lapply(seq_along(runs$lengths), function(i) {
    first <- ends[i] - runs$lengths[i] + 1L
    ids <- strsplit(runs$values[i], "|", fixed = TRUE)[[1]]
    at <- if (length(rows)) rows[rows$landmark_id %in% ids, , drop = FALSE] else NULL
    htmltools::tags$td(colspan = runs$lengths[i], class = "fungal-alignment-region-cell",
      `data-start-col` = columns[first], `data-end-col` = columns[ends[i]],
      style = if (length(at) && nrow(at)) paste0("background:", .fungal_alignment_ui_region_color(at[1, ]), "20;") else NULL,
      title = if (length(at) && nrow(at)) paste(vapply(seq_len(nrow(at)), function(j) .fungal_alignment_ui_region_title(at[j, ]), character(1)), collapse = "\n\n") else "No frozen boundary annotation",
      if (length(at) && nrow(at)) paste(vapply(seq_len(nrow(at)), function(j) .fungal_alignment_ui_region_label(at[j, ]), character(1)), collapse = " / ") else ""
    )
  })
}

.fungal_alignment_ui_pair_button <- function(pair, x, y, ns, prefix = "") {
  htmltools::tags$foreignObject(x = x, y = y - 15, width = 260, height = 32,
    htmltools::tags$div(xmlns = "http://www.w3.org/1999/xhtml",
      htmltools::tags$button(type = "button", class = "fungal-pair-choice",
        `data-input` = ns("pair_click"), `data-pair` = pair$pair_id,
        title = pair$pair_label,
        style = "white-space:nowrap;max-width:100%;overflow:hidden;text-overflow:ellipsis;line-height:1.3;font-size:13px;",
        paste0(prefix, pair$pair_label)
      )
    )
  )
}

.fungal_alignment_ui_span <- function(row, domain, scale, y, ceilings, pair_label = "") {
  if (!nrow(row)) {
    return(.fungal_alignment_ui_svg_text(scale(domain[1]) + 8, y + 4,
                                        .fungal_alignment_ui_status(row), size = 11, color = "#64748b"))
  }
  row <- row[1, , drop = FALSE]
  mapped <- identical(as.character(row$alignment_status[1]), "mapped")
  site_mapped <- function(region) {
    field <- paste0(region, "_alignment_status")
    if (field %in% names(row)) identical(as.character(row[[field]][1]), "mapped") else mapped
  }
  forward_mapped <- site_mapped("forward")
  reverse_mapped <- site_mapped("reverse")
  if (!mapped && !forward_mapped && !reverse_mapped) {
    return(.fungal_alignment_ui_svg_text(scale(domain[1]) + 8, y + 4,
      .fungal_alignment_ui_status(row), size = 11, color = "#64748b"))
  }
  tier <- .fungal_alignment_ui_tier(if (mapped) row$amplicon_bp else NA_real_, ceilings)
  if (mapped && (row$product_end_col < domain[1] || row$product_start_col > domain[2])) {
    return(.fungal_alignment_ui_svg_text(scale(domain[1]) + 8, y + 4,
      paste0("Outside this column window · ", row$amplicon_bp, " raw bp"), size = 11))
  }
  bounded <- function(value) pmax(domain[1], pmin(domain[2], value))
  arrow <- function(column, direction, region) {
    if (!is.finite(column) || !is.finite(direction) || column < domain[1] || column > domain[2]) return(NULL)
    tip <- scale(column)
    htmltools::tags$polygon(
      points = paste(paste0(tip + direction * 6, ",", y), paste0(tip - direction * 4, ",", y - 5),
                     paste0(tip - direction * 4, ",", y + 5)), fill = "#172033",
      class = "fungal-alignment-arrow", `data-primer-role` = region, `data-direction` = direction)
  }
  binding <- function(left, right, region) {
    if (!is.finite(left) || !is.finite(right) || right < domain[1] || left > domain[2]) return(NULL)
    htmltools::tags$line(x1 = scale(bounded(left)), x2 = scale(bounded(right)),
                        y1 = y, y2 = y, stroke = "#172033", `stroke-width` = 7,
                        class = "fungal-alignment-primer", `data-primer-role` = region,
                        `data-start-col` = left, `data-end-col` = right)
  }
  title <- if (mapped) paste0(pair_label, " · ", row$accession, "\nRaw product: ", row$amplicon_bp,
                  " bp including primers · ", tier$label,
                  "\nAlignment columns: ", row$product_start_col, "–", row$product_end_col,
                  "\nRaw sites: F ", row$forward_start, "–", row$forward_end,
                  "; R ", row$reverse_start, "–", row$reverse_end) else paste0(
    pair_label, " · ", row$accession, "\nProduct unresolved: ", .fungal_alignment_ui_status(row),
    "\nOnly verified exact individual binding sites are drawn. No product length or platform tier is assigned.")
  htmltools::tags$g(
    `data-accession` = row$accession, `data-pair-id` = row$pair_id,
    `data-alignment-status` = row$alignment_status,
    `data-product-bp` = if (mapped) row$amplicon_bp else NULL,
    `data-start-col` = if (mapped) row$product_start_col else NULL,
    `data-end-col` = if (mapped) row$product_end_col else NULL,
    `data-coordinate-system` = "alignment_columns",
    htmltools::tags$title(title),
    if (mapped) htmltools::tags$line(class = "fungal-alignment-product", x1 = scale(bounded(row$product_start_col)), x2 = scale(bounded(row$product_end_col)),
                        y1 = y, y2 = y, stroke = tier$color, `stroke-width` = 5, `stroke-linecap` = "round"),
    if (forward_mapped) binding(row$forward_left_col, row$forward_right_col, "forward"),
    if (reverse_mapped) binding(row$reverse_left_col, row$reverse_right_col, "reverse"),
    if (forward_mapped) arrow(row$forward_end_col, row$forward_direction, "forward"),
    if (reverse_mapped) arrow(row$reverse_start_col, row$reverse_direction, "reverse"),
    .fungal_alignment_ui_svg_text(1030, y + 4,
      if (mapped) paste0(row$amplicon_bp, " raw bp") else "Unresolved product", size = if (mapped) 12 else 10,
      color = tier$color)
  )
}

fungal_alignment_overview <- function(index, mapped_rows, pairs, taxonomy, reference_accession,
                                     selected_accessions, compare_ids, active_id, domain, ceilings, ns,
                                     landmarks = NULL, show_regions = TRUE) {
  tags <- htmltools::tags
  if (!length(index$alignment) || !is.finite(index$n_columns) || index$n_columns < 1L) {
    return(tags$div(class = "callout", "No frozen alignment is available in this snapshot."))
  }
  domain <- .fungal_alignment_ui_domain(index, domain)
  selected_accessions <- intersect(as.character(selected_accessions), names(index$alignment))
  compare_ids <- head(unique(intersect(as.character(compare_ids), pairs$pair_id)), 5L)
  if (!length(compare_ids) && length(active_id) && active_id %in% pairs$pair_id) compare_ids <- active_id
  priorities <- unique(intersect(c(as.character(active_id), compare_ids), pairs$pair_id))
  shown <- if (nrow(pairs) > 20L) {
    head(pairs[match(unique(c(priorities, pairs$pair_id)), pairs$pair_id), , drop = FALSE], 20L)
  } else pairs
  comparisons <- pairs[match(compare_ids, pairs$pair_id), , drop = FALSE]
  scale <- function(column) 290 + (column - domain[1]) / max(1L, diff(domain)) * 710
  visible_regions <- if (isTRUE(show_regions))
    .fungal_alignment_ui_visible_regions(landmarks, reference_accession, domain) else NULL
  region_height <- if (isTRUE(show_regions)) 72L else 0L
  reference_height <- 90 + region_height + nrow(shown) * 32
  reference_svg <- tags$svg(viewBox = paste("0 0 1160", reference_height), role = "img",
    `aria-label` = paste("Frozen alignment columns", paste(domain, collapse = "through"),
                        "with", nrow(shown), "candidate primer combinations",
                        if (isTRUE(show_regions)) "and template-specific genes and regions", "on", reference_accession),
    if (isTRUE(show_regions)) {
      lapply(seq_len(if (length(visible_regions)) nrow(visible_regions) else 0L), function(i) {
        row <- visible_regions[i, , drop = FALSE]
        left <- scale(max(domain[1], row$left_col))
        right <- scale(min(domain[2], row$right_col))
        tags$rect(class = "fungal-alignment-region-band", `data-accession` = reference_accession,
          `data-landmark-id` = row$landmark_id, x = left, y = 32,
          width = max(1, right - left), height = reference_height - 42,
          fill = .fungal_alignment_ui_region_color(row), `fill-opacity` = .055,
          `pointer-events` = "none")
      })
    },
    .fungal_alignment_ui_axis(domain, scale, reference_height - 10),
    if (isTRUE(show_regions)) htmltools::tagList(
      .fungal_alignment_ui_svg_text(10, 58, "Genes / regions", size = 12, `font-weight` = 700),
      .fungal_alignment_ui_svg_text(10, 76, "On this reference template", size = 10),
      .fungal_alignment_ui_regions(landmarks, reference_accession, domain, scale, 64)),
    .fungal_alignment_ui_svg_text(10, 52 + region_height, reference_accession, size = 14, color = "#172033", `font-weight` = 800),
    .fungal_alignment_ui_occupancy(index, reference_accession, domain, scale, 50 + region_height),
    lapply(seq_len(nrow(shown)), function(i) {
      pair <- shown[i, , drop = FALSE]
      row <- .fungal_alignment_ui_row(mapped_rows, reference_accession, pair$pair_id)
      y <- 82 + region_height + (i - 1L) * 32
      tags$g(
        if (length(active_id) && pair$pair_id == active_id) tags$rect(x = 0, y = y - 17, width = 1150, height = 30,
                                                                  rx = 5, fill = "#eef2ff"),
        .fungal_alignment_ui_pair_button(pair, 6, y, ns),
        .fungal_alignment_ui_span(row, domain, scale, y, ceilings, pair$pair_label)
      )
    })
  )
  block_height <- 45 + region_height + max(1L, nrow(comparisons)) * 19
  template_height <- 45 + length(selected_accessions) * block_height
  template_svg <- tags$svg(viewBox = paste("0 0 1160", template_height), role = "img",
    `aria-label` = paste("Frozen alignment non-gap occupancy and selected primer products for", length(selected_accessions), "templates"),
    .fungal_alignment_ui_axis(domain, scale, template_height - 5),
    lapply(seq_along(selected_accessions), function(i) {
      accession <- selected_accessions[i]
      y <- 49 + (i - 1L) * block_height
      detail <- .fungal_alignment_ui_taxon(taxonomy, accession)
      tags$g(
        `data-accession` = accession,
        .fungal_alignment_ui_svg_text(10, y + 4, accession, size = 13, color = "#172033", `font-weight` = 800),
        tags$text(x = 10, y = y + 19, fill = "#475569", `font-size` = 10,
                  if (nchar(detail) > 43) paste0(substr(detail, 1, 40), "…") else detail, tags$title(detail)),
        .fungal_alignment_ui_occupancy(index, accession, domain, scale, y),
        if (isTRUE(show_regions)) .fungal_alignment_ui_regions(landmarks, accession, domain, scale, y + 38),
        lapply(seq_len(nrow(comparisons)), function(j) {
          pair <- comparisons[j, , drop = FALSE]
          lane_y <- y + 26 + region_height + (j - 1L) * 19
          row <- .fungal_alignment_ui_row(mapped_rows, accession, pair$pair_id)
          tags$g(
            .fungal_alignment_ui_svg_text(260, lane_y + 4, paste0("P", j), size = 10),
            .fungal_alignment_ui_span(row, domain, scale, lane_y, ceilings, pair$pair_label)
          )
        })
      )
    })
  )
  tags$div(class = "fungal-alignment-overview", `data-coordinate-system` = "alignment_columns",
    `data-alignment-columns` = index$n_columns, `data-template-count` = length(selected_accessions),
    `data-reference-accession` = reference_accession, `data-start-col` = domain[1], `data-end-col` = domain[2],
    tags$h4("Primers on the frozen alignment"),
    tags$p(class = "tiny", paste0("Alignment columns ", format(domain[1], big.mark = ","), "–",
      format(domain[2], big.mark = ","), " of ", format(index$n_columns, big.mark = ","),
      ". Grey marks show non-gap sequence; empty stretches are alignment gaps. Product labels count raw template bases including primers.")),
    tags$h5(paste0("Candidate combinations on ", reference_accession)),
    if (isTRUE(show_regions)) tags$div(class = "fungal-alignment-window-regions",
      tags$p(class = "tiny", "Regions in this window · ", reference_accession),
      if (length(visible_regions) && nrow(visible_regions)) lapply(seq_len(nrow(visible_regions)), function(i) {
        row <- visible_regions[i, , drop = FALSE]
        clipped <- row$left_col < domain[1] || row$right_col > domain[2]
        tags$span(class = paste("fungal-alignment-window-region", if (grepl("^ITS[12]$", row$label)) "spacer" else "gene"),
          `data-accession` = row$accession, `data-landmark-id` = row$landmark_id,
          `data-start-col` = row$left_col, `data-end-col` = row$right_col,
          `data-clipped` = tolower(as.character(clipped)), title = .fungal_alignment_ui_region_title(row),
          tags$strong(.fungal_alignment_ui_region_label(row)),
          tags$span(class = "tiny", paste0(" · columns ", format(row$left_col, big.mark = ","), "–",
            format(row$right_col, big.mark = ","), if (clipped) " · partly in view" else "")))
      }) else tags$p(class = "tiny", "No frozen gene or region boundary overlaps this window.")),
    tags$p(class = "tiny fungal-alignment-pan-hint", "Scroll sideways to inspect the full alignment track and product lengths."),
    tags$div(class = "fungal-alignment-scroll", reference_svg),
    tags$h5(paste0("Selected templates on the same columns (", length(selected_accessions), ")")),
    if (nrow(comparisons)) tags$p(class = "tiny", paste0("P", seq_len(nrow(comparisons)), " = ", comparisons$pair_label, collapse = " · ")),
    if (length(selected_accessions)) tags$div(class = "fungal-alignment-scroll", template_svg)
      else tags$p("No alignment templates match the current lineage filters."),
    tags$p(class = "tiny", "Colors describe each template's raw product length tier. Dark primer blocks and arrowheads show binding spans and extension direction.",
      if (isTRUE(show_regions)) " Slate tracks mark rRNA genes; amber tracks mark ITS spacers. Each annotation belongs to its labeled template. Dashed borders identify curated pilot boundaries; partial anchors are labeled.",
      " Across-clade variable-region alignment does not establish positional homology or experimental amplification performance.")
  )
}

fungal_alignment_bases <- function(index, mapped_rows, pairs, taxonomy, accessions,
                                  start, width = 100L, active_id, ceilings,
                                  landmarks = NULL, show_regions = TRUE) {
  tags <- htmltools::tags
  accessions <- intersect(as.character(accessions), names(index$alignment))
  if (!length(accessions)) return(tags$p("No aligned templates match the current filters."))
  width <- max(1L, min(120L, as.integer(width)))
  if (!length(start) || !is.finite(start)) start <- 1L
  start <- max(1L, min(index$n_columns, as.integer(start)))
  window <- fungal_alignment_window(index, mapped_rows = mapped_rows, accessions = accessions,
                                    start = start, max_columns = width)
  columns <- window$columns
  active <- if (length(active_id)) pairs[pairs$pair_id == active_id, , drop = FALSE] else pairs[0, , drop = FALSE]
  active_label <- if (nrow(active)) active$pair_label[1] else "No active combination"
  body <- lapply(accessions, function(accession) {
    bases <- strsplit(substr(index$alignment[[accession]], min(columns), max(columns)), "", fixed = TRUE)[[1]]
    raw_positions <- match(columns, index$maps[[accession]])
    row <- .fungal_alignment_ui_row(mapped_rows, accession, active_id)
    if (nrow(row)) row <- row[1, , drop = FALSE]
    mapped <- nrow(row) > 0L && identical(as.character(row$alignment_status[1]), "mapped")
    site_mapped <- function(region) {
      field <- paste0(region, "_alignment_status")
      if (!nrow(row)) FALSE else if (field %in% names(row)) identical(as.character(row[[field]][1]), "mapped") else mapped
    }
    forward_mapped <- site_mapped("forward")
    reverse_mapped <- site_mapped("reverse")
    tier <- .fungal_alignment_ui_tier(if (mapped) row$amplicon_bp else NA_real_, ceilings)
    detail <- .fungal_alignment_ui_taxon(taxonomy, accession)
    annotation <- .fungal_alignment_ui_visible_regions(landmarks, accession, range(columns))
    htmltools::tagList(
    if (isTRUE(show_regions) && length(annotation) && nrow(annotation)) tags$tr(
      class = "fungal-alignment-region-row", `data-accession` = accession,
      tags$th(scope = "row", class = "fungal-alignment-label", accession, " · regions"),
      .fungal_alignment_ui_region_cells(landmarks, accession, columns)),
    tags$tr(
      `data-accession` = accession,
      tags$th(scope = "row", class = "fungal-alignment-label", accession,
              tags$div(class = "tiny", detail),
              tags$div(class = "tiny", if (mapped) paste0(row$amplicon_bp, " raw bp · ", tier$label)
                       else .fungal_alignment_ui_status(row))),
      lapply(seq_along(columns), function(i) {
        column <- columns[i]
        gap <- bases[i] == "-"
        in_product <- mapped && column >= row$product_start_col && column <= row$product_end_col
        in_forward <- forward_mapped && column >= row$forward_left_col && column <= row$forward_right_col
        in_reverse <- reverse_mapped && column >= row$reverse_left_col && column <= row$reverse_right_col
        classes <- c("fungal-alignment-base", if (gap) "gap", if (!gap && in_product) "product",
                     if (!gap && in_forward) "forward", if (!gap && in_reverse) "reverse")
        style <- if (!gap && (in_forward || in_reverse)) paste0("background:", tier$color, "26;border-bottom:3px solid ", tier$color, ";")
          else if (!gap && in_product) paste0("background:", tier$color, "12;") else NULL
        region <- if (isTRUE(show_regions) && !gap && length(annotation) && nrow(annotation))
          annotation[annotation$start <= raw_positions[i] & annotation$end >= raw_positions[i], , drop = FALSE] else NULL
        title <- paste0(accession, " · alignment column ", column, " · ",
                        if (gap) "gap; no raw base" else paste0("raw position ", raw_positions[i]),
                        if (length(region) && nrow(region)) paste0(" · ", paste(region$label, collapse = " / ")) else "",
                        if (!gap && in_forward) " · forward primer site" else if (!gap && in_reverse) " · reverse primer site" else "")
        tags$td(class = paste(classes, collapse = " "), style = style, title = title,
                `data-alignment-col` = column, `data-raw-position` = if (!gap) raw_positions[i] else NULL, bases[i])
      })
    ))
  })
  tags$div(class = "fungal-alignment-bases", `data-coordinate-system` = "alignment_columns",
    `data-start-col` = min(columns), `data-end-col` = max(columns), `data-template-count` = length(accessions),
    tags$h5(paste0("Actual aligned bases · columns ", min(columns), "–", max(columns))),
    tags$p(class = "tiny", paste0(active_label, ". The letters and gaps come directly from this frozen alignment. Primer sites receive dark borders; product shading excludes gap cells. Hover a base for its raw position.")),
    tags$div(class = "fungal-alignment-base-scroll",
      tags$table(class = "fungal-alignment-base-table",
        tags$thead(tags$tr(tags$th("Template / raw product"), lapply(columns, function(column) {
          tags$th(title = paste("Alignment column", column), if (column == min(columns) || column == max(columns) || column %% 10L == 0L) column else "")
        }))), tags$tbody(body)
      )
    )
  )
}
