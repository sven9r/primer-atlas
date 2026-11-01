.primer_overview_domain <- function(
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


.primer_overview_use_tags <- function(use_case, target = NA_character_) {
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


# Shared COI-style primer overview. All cohort, geometry, thermal and source
# data are explicit arguments; alignment columns must never be passed as bp.
primer_overview_map <- function(rows,primers,marker,landmarks=data.frame(),focus=NULL,
  thermal_rows=data.frame(),sodium_mM=50,ta_offset_c=4,citations=data.frame(),
  primer_facets=data.frame(),diet_uses=data.frame(),arthropod_off_target_ids=character(),
  folmer_region=NULL,domain=NULL,active_pair_id=NULL,pair_click_input=NULL,viewport_width=NULL) {
  tags <- htmltools::tags
  div <- tags$div
  tagList <- htmltools::tagList
  if(!is.null(viewport_width) && (length(viewport_width)!=1L || !is.numeric(viewport_width) ||
      !is.finite(viewport_width) || viewport_width<=0)) stop("The map viewport width must be one positive pixel measurement.")
  compact <- !is.null(viewport_width) && viewport_width<1600
  x <- rows
  if(!is.data.frame(x) || !nrow(x)) return(div(class="callout","No primer pairs match the current filters. "))
  stopifnot(is.data.frame(marker),nrow(marker)==1L)
  required <- c("pair_id","pair_label","forward_start","forward_end","reverse_start","reverse_end",
    "pair_start","pair_end","aligned_amplicon_bp","targeted_region_bp")
  if(!all(required %in% names(x))) stop("Primer overview geometry is missing: ",paste(setdiff(required,names(x)),collapse=", "))
  complete <- stats::complete.cases(x[,required[3:10],drop=FALSE])
  if(!all(complete)) stop("Primer overview rows require complete measured intervals; unresolved products cannot be drawn as complete lanes")
  for(name in c("use_case","target","placement_status")) if(!name %in% names(x)) x[[name]] <- rep(NA_character_,nrow(x))
  if(!"folmer_overlap_bp" %in% names(x)) x$folmer_overlap_bp <- rep(NA_real_,nrow(x))
  ensure_columns <- function(data,columns) {
    stopifnot(is.data.frame(data))
    for(column in columns) if(!column %in% names(data)) data[[column]] <- rep(NA_character_,nrow(data))
    data
  }
  primers <- ensure_columns(primers,c("pair_id","direction","primer_name","source_key"))
  citations <- ensure_columns(citations,c("key","url","short_citation","title"))
  primer_facets <- ensure_columns(primer_facets,c("pair_id","facet_type","facet_value"))
  diet_uses <- ensure_columns(diet_uses,c("pair_id","study_use","predator","sample_material","amplified_role","evidence_status","source_key","note"))
  thermal_rows <- ensure_columns(thermal_rows,c("pair_id","direction","primer_name","tm_min","tm_max","ta_low","ta_high","reported_ta_c"))
  landmarks <- ensure_columns(landmarks,c("landmark_id","label","start","end","landmark_type","source_note"))
  marker_id <- as.character(marker$marker_id[1L])
  reference_length <- as.numeric(marker$reference_length[1L])
  if(length(reference_length)!=1L || !is.finite(reference_length) || reference_length<=0L) stop("A positive reference length is required")
  map_primers <- primers
  if(!"source_key" %in% names(map_primers)) map_primers$source_key <- rep(NA_character_,nrow(map_primers))
  if(is.null(focus)) {
    focus <- if(marker_id=="COI" && !is.null(folmer_region)) folmer_region else list(start=1,end=reference_length)
  }
  range_values <- function(value) {
    if(is.list(value)) {
      if(!is.null(value$start) && !is.null(value$end)) value <- c(value$start,value$end)
      else if(!is.null(value$min) && !is.null(value$max)) value <- c(value$min,value$max)
    }
    value <- as.numeric(value)
    if(length(value)!=2L || any(!is.finite(value)) || value[2L]<=value[1L]) stop("A bp view requires finite increasing start/end coordinates")
    value
  }
  focus_bounds <- range_values(focus)
  focus <- list(start=focus_bounds[1L],end=focus_bounds[2L])
    width <- 1680
    left <- 112
    right <- 1000
    if(compact) {
      width <- max(240L,as.integer(floor(viewport_width)))
      left <- 44
      right <- width-24
    }
    axis_domain <- if(is.null(domain)) .primer_overview_domain(
      x,
      folmer_start = focus$start,
      folmer_end = focus$end,
      tick_bp = 100,
      reference_length = reference_length
    ) else setNames(range_values(domain),c("min","max"))
    if(!is.null(domain) && (axis_domain[1L]<0L || axis_domain[2L]>reference_length)) stop("The raw bp view lies outside the supplied reference")
    axis_min <- unname(axis_domain["min"])
    axis_max <- unname(axis_domain["max"])
    scale_x <- function(bp) {
      bounded <- pmin(pmax(bp, axis_min), axis_max)
      left + (bounded - axis_min) / (axis_max - axis_min) * (right - left)
    }
    row_h <- if(compact) 82 else 58
    first_row_y <- 88
    last_row_y <- first_row_y + (nrow(x) - 1) * row_h
    folmer_y <- last_row_y + if(compact) 70 else 54
    landmark_coordinates <- suppressWarnings(cbind(as.numeric(landmarks$start),as.numeric(landmarks$end)))
    landmark_valid <- if(nrow(landmarks)) stats::complete.cases(landmark_coordinates) &
      landmark_coordinates[,1L]>=1L & landmark_coordinates[,2L]>=landmark_coordinates[,1L] &
      (is.null(domain) | landmark_coordinates[,2L]<=reference_length) else logical()
    if(nrow(landmarks)) {
      landmarks$start <- landmark_coordinates[,1L]
      landmarks$end <- landmark_coordinates[,2L]
    }
    domain_landmarks <- landmarks[landmark_valid & landmarks$landmark_type %in% c("domain","divergent_domain"),,drop=FALSE]
    gene_landmarks <- landmarks[landmark_valid & !landmarks$landmark_type %in% c("domain","divergent_domain"),,drop=FALSE]
    domain_y <- folmer_y+30
    coi_y <- folmer_y + if(compact && marker_id=="COI" && !is.null(folmer_region)) 70L else 40L
    coi_y <- coi_y + if(nrow(domain_landmarks)) 32L else 0L
    height <- coi_y + 42
    start_col <- 1090
    end_col <- 1160
    inclusive_col <- 1250
    informative_col <- 1350
    citation_col <- 1440
    number_x <- 16
    measurements <- list()
    forward_color <- "#2563eb"
    informative_color <- "#9aa6b2"
    reverse_color <- "#dc2626"
    shared_source_color <- "#7c3aed"
    source_link <- function(source_row, x, y, color, prefix = "") {
      if (!nrow(source_row) || is.na(source_row$url[1]) ||
        !nzchar(source_row$url[1])) {
        return(tags$text(
          x = x, y = y, fill = "#6b746f", `font-size` = 10.5,
          paste0(prefix, "source unavailable")
        ))
      }
      tags$a(
        href = source_row$url[1],
        target = "_blank",
        rel = "noopener noreferrer",
        style = "cursor:pointer",
        tags$text(
          x = x, y = y, fill = color, `font-size` = 10.5,
          `text-decoration` = "underline",
          paste0(prefix, source_row$short_citation[1]),
          tags$title(source_row$title[1])
        )
      )
    }
    children <- list(
      tags$text(
        x = start_col, y = 20, `text-anchor` = "middle",
        fill = "#48554e", `font-size` = 11, `font-weight` = 700, "Start bp"
      ),
      tags$text(
        x = end_col, y = 20, `text-anchor` = "middle",
        fill = "#48554e", `font-size` = 11, `font-weight` = 700, "End bp"
      ),
      tags$text(
        x = inclusive_col, y = 20, `text-anchor` = "middle",
        fill = "#48554e", `font-size` = 11, `font-weight` = 700, "+ primers"
      ),
      tags$text(
        x = informative_col, y = 20, `text-anchor` = "middle",
        fill = "#48554e", `font-size` = 11, `font-weight` = 700, "− primers"
      ),
      tags$text(
        x = citation_col, y = 20,
        fill = "#48554e", `font-size` = 11, `font-weight` = 700, "Citation(s)"
      ),
      tags$text(
        x = right, y = 20, `text-anchor` = "end",
        fill = "#48554e", `font-size` = if(compact) 13 else 11, `font-weight` = 700,
        paste0("View ", axis_min, "–", axis_max, " bp")
      ),
      tags$line(
        x1 = left, y1 = 52, x2 = right, y2 = 52,
        stroke = "#6d756f", `stroke-width` = 1.5
      )
    )
    if(compact) children <- children[6:7]
    tick_bp <- if(is.null(domain) || axis_max-axis_min<=1500L) 100L else if(axis_max-axis_min<=5000L) 500L else 1000L
    if(compact) {
      target_step <- (axis_max-axis_min)/max(1,(right-left)/62)
      steps <- as.vector(outer(c(1,2,5),10^seq(0,8)))
      tick_bp <- min(steps[steps>=target_step])
    }
    ticks <- if(is.null(domain) && !compact) seq(axis_min,axis_max,100L) else {
      limits <- c(ceiling(axis_min/tick_bp)*tick_bp,floor(axis_max/tick_bp)*tick_bp)
      if(limits[2L]>=limits[1L]) seq(limits[1L],limits[2L],tick_bp) else numeric()
    }
    for (tick in ticks) {
      tx <- scale_x(tick)
      children <- append(children, list(
        tags$line(
          x1 = tx, y1 = 46, x2 = tx, y2 = last_row_y + 25,
          stroke = "#dedbd2", `stroke-width` = 1
        ),
        tags$text(
          x = tx, y = 40, `text-anchor` = "middle",
          fill = "#68736c", `font-size` = if(compact) 12 else 11, tick
        )
      ))
    }
    children <- append(children, list(
      tags$rect(
        x = scale_x(focus$start), y = 55,
        width = scale_x(focus$end) - scale_x(focus$start),
        height = last_row_y - 30, fill = "#f5df73", opacity = .15
      )
    ))
    for (i in seq_len(nrow(x))) {
      row <- x[i, ]
      product_color <- informative_color
      if("platform_tier" %in% names(row) && !is.na(row$platform_tier)) product_color <- switch(as.character(row$platform_tier),
        ILLUMINA="#2563eb",PACBIO="#7c3aed",NANOPORE="#0f766e",BEYOND="#64748b",informative_color)
      y <- first_row_y + (i - 1) * row_h
      f <- map_primers |> dplyr::filter(pair_id == row$pair_id, direction == "forward") |> dplyr::slice(1)
      r <- map_primers |> dplyr::filter(pair_id == row$pair_id, direction == "reverse") |> dplyr::slice(1)
      if(!nrow(f)) f <- data.frame(primer_name=if("forward_primer" %in% names(row)) row$forward_primer else "Forward primer",source_key=NA_character_)
      if(!nrow(r)) r <- data.frame(primer_name=if("reverse_primer" %in% names(row)) row$reverse_primer else "Reverse primer",source_key=NA_character_)
      forward_visible <- row$forward_end>=axis_min && row$forward_start<=axis_max
      reverse_visible <- row$reverse_end>=axis_min && row$reverse_start<=axis_max
      forward_tip_visible <- row$forward_end>=axis_min && row$forward_end<=axis_max
      reverse_tip_visible <- row$reverse_start>=axis_min && row$reverse_start<=axis_max
      product_visible <- row$pair_end>=axis_min && row$pair_start<=axis_max
      product_clipped <- row$pair_start<axis_min || row$pair_end>axis_max
      f1 <- scale_x(row$forward_start)
      f2 <- scale_x(row$forward_end)
      r1 <- scale_x(row$reverse_start)
      r2 <- scale_x(row$reverse_end)
      forward_source <- citations |> dplyr::filter(key == f$source_key[1]) |> dplyr::slice(1)
      reverse_source <- citations |> dplyr::filter(key == r$source_key[1]) |> dplyr::slice(1)
      application_values <- primer_facets$facet_value[
        primer_facets$pair_id == row$pair_id &
          primer_facets$facet_type == "application"
      ]
      application_labels <- c(
        barcoding = "DNA barcoding", bulk_community = "bulk community",
        edna = "environmental DNA", diet = "dietary metabarcoding"
      )
      use_tags <- if (length(application_values)) {
        unname(application_labels[application_values])
      } else {
        .primer_overview_use_tags(row$use_case, row$target)
      }
      diet <- diet_uses |> dplyr::filter(pair_id == row$pair_id)
      diet_detail <- if (nrow(diet)) paste(paste0(
        "\nDocumented study use: ", diet$study_use,
        "; predator: ", diet$predator,
        "; sample: ", diet$sample_material,
        "; amplified target: ", diet$amplified_role,
        "\nDiet evidence status: ", diet$evidence_status,
        "; source: ", diet$source_key,
        "\nUse note: ", diet$note
      ), collapse = "\n") else ""
      pair_temperatures <- thermal_rows |>
        dplyr::filter(pair_id == row$pair_id) |>
        dplyr::arrange(match(direction, c("forward", "reverse")))
      tm_detail <- if (nrow(pair_temperatures)) {
        paste0(
          "\nEstimated primer Tm (", sodium_mM, " mM salt): ",
          paste(sprintf(
            "%s %.1f–%.1f °C",
            pair_temperatures$primer_name,
            pair_temperatures$tm_min,
            pair_temperatures$tm_max
          ), collapse = "; "),
          "\nSuggested pair Ta (starting estimate): ",
          sprintf("%.1f–%.1f °C", pair_temperatures$ta_low[1], pair_temperatures$ta_high[1]),
          "\nReported Ta: ",
          if (all(is.na(pair_temperatures$reported_ta_c))) {
            "not reported"
          } else {
            paste0(
              paste(unique(na.omit(pair_temperatures$reported_ta_c)), collapse = ", "),
              " °C"
            )
          }
        )
      } else "\nTemperatures: unavailable for this pair"
      tooltip <- paste0(
        "Map row ", i, ": ", row$pair_label,
        "\nTags: ", paste(use_tags, collapse = " · "),
        "\nIntended target: ", row$target,
        diet_detail,
        if (row$pair_id %in% arthropod_off_target_ids) {
          "\nArthropod reference: off-target binding check only; target coverage not established here"
        } else "",
        "\nPrimers: ", f$primer_name, " + ", r$primer_name,
        tm_detail,
        "\nGene coordinates: ", row$pair_start, "–", row$pair_end,
        "\nAmplicon + primers: ", row$aligned_amplicon_bp, " bp",
        "\nAmplicon − primers (informative): ", row$targeted_region_bp, " bp",
        if (marker_id == "COI") paste0("\nInformative bp inside Folmer: ", row$folmer_overlap_bp) else "",
        "\nCoordinates derived from ", marker$coordinate_reference,
        if(!is.null(domain) && !product_visible) "\nOutside this bp window; full measured coordinates and lengths are retained." else
          if(!is.null(domain) && product_clipped) "\nPartly in this bp window; arrowheads appear only at real primer 3′ endpoints within the view." else "",
        if (!is.na(row$placement_status) && row$placement_status != "sequence_aligned") {
          paste0(
            if (marker_id == "18S") "\nPlacement provenance: " else "\nPlacement warning: ",
            row$placement_status
          )
        } else ""
      )
      number_tooltip <- paste0(
        "Map row ", i, "\nPair: ", row$pair_label,
        "\nOligos: ", f$primer_name, " + ", r$primer_name,
        "\nSorting and filtering recalculate this number."
      )
      row_children <- list(
        tags$circle(
          cx = number_x, cy = y, r = 13,
          fill = "#17241d",
          class = "primer-tooltip-target",
          `data-tooltip` = tooltip,
          `aria-label` = paste0("Show map-row details for ", row$pair_label),
          tabindex = "0",
          tags$title(number_tooltip)
        ),
        tags$text(
          x = number_x, y = y + 4,
          `text-anchor` = "middle",
          fill = "#fffdf8", `font-size` = 11, `font-weight` = 800,
          class = "primer-tooltip-target",
          `data-tooltip` = tooltip,
          `aria-label` = paste0("Show map-row details for ", row$pair_label),
          tabindex = "0",
          i,
          tags$title(number_tooltip)
        ),
        tags$line(
          x1 = left, y1 = y + row_h / 2 - 2,
          x2 = width - 12, y2 = y + row_h / 2 - 2,
          stroke = "#e5e2d9", `stroke-width` = 1
        ),
        tags$line(
          x1 = f2, y1 = y, x2 = r1, y2 = y,
          class = "primer-tooltip-target",
          `data-tooltip` = tooltip,
          `aria-label` = paste0("Show map-row details for ", row$pair_label),
          tabindex = "0",
          stroke = product_color, `stroke-width` = 5, `stroke-linecap` = "round",
          `pointer-events` = "stroke"
        ),
        tags$line(
          x1 = f1, y1 = y, x2 = f2, y2 = y,
          class = "primer-tooltip-target",
          `data-tooltip` = tooltip,
          `aria-label` = paste0("Show map-row details for ", row$pair_label),
          tabindex = "0",
          stroke = forward_color, `stroke-width` = 7, `stroke-linecap` = "round"
        ),
        tags$polygon(
          points = paste(f2, y, f2 - 8, y - 6, f2 - 8, y + 6),
          fill = forward_color,
          class = "primer-tooltip-target",
          `data-tooltip` = tooltip,
          `aria-label` = paste0("Show map-row details for ", row$pair_label),
          tabindex = "0"
        ),
        tags$text(
          x = if(compact) max(left,min(f1,right-nchar(f$primer_name)*8.4)) else f1 - 7,
          y = if(compact) y-16 else y + 4, `text-anchor` = if(compact) "start" else "end",
          class = "primer-tooltip-target",
          `data-tooltip` = tooltip,
          `aria-label` = paste0("Show map-row details for ", row$pair_label),
          tabindex = "0",
          fill = forward_color, `font-size` = if(compact) 14 else 12, `font-weight` = 700,
          f$primer_name
        ),
        tags$line(
          x1 = r1, y1 = y, x2 = r2, y2 = y,
          class = "primer-tooltip-target",
          `data-tooltip` = tooltip,
          `aria-label` = paste0("Show map-row details for ", row$pair_label),
          tabindex = "0",
          stroke = reverse_color, `stroke-width` = 7, `stroke-linecap` = "round"
        ),
        tags$polygon(
          points = paste(r1, y, r1 + 8, y - 6, r1 + 8, y + 6),
          fill = reverse_color,
          class = "primer-tooltip-target",
          `data-tooltip` = tooltip,
          `aria-label` = paste0("Show map-row details for ", row$pair_label),
          tabindex = "0"
        ),
        tags$text(
          x = if(compact) min(right,max(r2,left+nchar(r$primer_name)*8.4)) else r2 + 7,
          y = if(compact) y+22 else y + 4, `text-anchor` = if(compact) "end" else "start",
          class = "primer-tooltip-target",
          `data-tooltip` = tooltip,
          `aria-label` = paste0("Show map-row details for ", row$pair_label),
          tabindex = "0",
          fill = reverse_color, `font-size` = if(compact) 14 else 12, `font-weight` = 700,
          r$primer_name
        ),
        tags$text(
          x = start_col, y = y + 4, `text-anchor` = "middle",
          fill = "#34423a", `font-size` = 11, row$pair_start
        ),
        tags$text(
          x = end_col, y = y + 4, `text-anchor` = "middle",
          fill = "#34423a", `font-size` = 11, row$pair_end
        ),
        tags$text(
          x = inclusive_col, y = y + 4, `text-anchor` = "middle",
          fill = "#34423a", `font-size` = 11,
          paste0(row$aligned_amplicon_bp, " bp")
        ),
        tags$text(
          x = informative_col, y = y + 4, `text-anchor` = "middle",
          fill = "#34423a", `font-size` = 11,
          paste0(row$targeted_region_bp, " bp")
        )
      )
      if(compact) row_children <- row_children[1:10]
      if(!is.null(domain)) {
        remove <- c(if(!forward_visible) 5L,if(!forward_tip_visible) 6L,
          if(!reverse_visible) 8L,if(!reverse_tip_visible) 9L,
          if(!product_visible || row$reverse_start-1L<axis_min || row$forward_end+1L>axis_max) 4L)
        if(length(remove)) row_children <- row_children[-unique(remove)]
        if(!product_visible || product_clipped) row_children <- append(row_children,list(tags$text(
          x=(left+right)/2,y=y+if(compact) 42 else 20,fill="#6b746f",`font-size`=if(compact) 12 else 10,
          if(compact) paste0(if(!product_visible) "Outside this bp window" else "Partly in view", " · ",row$aligned_amplicon_bp," bp full length") else
            if(!product_visible) "Outside this bp window" else "Partly in view · full bp length retained")))
      }
      if(compact && product_visible && !product_clipped) row_children <- append(row_children,list(tags$text(
        class="primer-overview-length",x=(left+right)/2,y=y+42,`text-anchor`="middle",
        fill="#526176",`font-size`=12,
        paste0(row$aligned_amplicon_bp," bp + primers · ",row$targeted_region_bp," bp internal"))))
      if(compact) {
        html_source <- function(source_row,prefix="") {
          if(!nrow(source_row)) return(tags$span(paste0(prefix,"source unavailable")))
          label <- as.character(source_row$short_citation[1L])
          if(is.na(label) || !nzchar(label)) label <- "source unavailable"
          url <- as.character(source_row$url[1L])
          if(is.na(url) || !nzchar(url)) tags$span(paste0(prefix,label)) else
            tags$a(href=url,target="_blank",rel="noopener noreferrer",paste0(prefix,label))
        }
        sources <- if(identical(f$source_key[1],"custom") || identical(r$source_key[1],"custom")) tags$span("Session input · not saved") else
          if(identical(f$source_key[1],r$source_key[1])) html_source(forward_source) else
            tagList(html_source(forward_source,"F: "),tags$br(),html_source(reverse_source,"R: "))
        label <- paste0(i," · ",row$pair_label)
        measurements[[i]] <- tags$tr(`data-pair-id`=row$pair_id,`data-product-bp`=row$aligned_amplicon_bp,
          tags$th(scope="row",if(is.null(pair_click_input)) label else tags$button(type="button",class="pair-inspect primer-overview-measurement-button",
            `data-pair-id`=row$pair_id,`data-pair-input`=pair_click_input,`aria-pressed`=if(row$pair_id %in% active_pair_id) "true" else "false",label)),
          tags$td(row$pair_start),tags$td(row$pair_end),tags$td(paste0(row$aligned_amplicon_bp," bp")),
          tags$td(paste0(row$targeted_region_bp," bp")),tags$td(sources))
      } else if (identical(f$source_key[1], "custom") ||
        identical(r$source_key[1], "custom")) {
        row_children <- append(row_children, list(tags$text(
          x = citation_col, y = y + 4, fill = "#6b746f", `font-size` = 10.5,
          "Session input · not saved"
        )))
      } else if (identical(f$source_key[1], r$source_key[1])) {
        row_children <- append(row_children, list(
          source_link(forward_source, citation_col, y + 4, shared_source_color)
        ))
      } else {
        row_children <- append(row_children, list(
          source_link(forward_source, citation_col, y - 5, forward_color, "F: "),
          source_link(reverse_source, citation_col, y + 11, reverse_color, "R: ")
        ))
      }
      if(!is.null(pair_click_input) || !is.null(active_pair_id) || !is.null(domain)) {
        row_group <- tags$g(class="primer-overview-row",`data-pair-id`=row$pair_id,
          `data-start-bp`=row$pair_start,`data-end-bp`=row$pair_end,`data-product-bp`=row$aligned_amplicon_bp,
          `data-forward-start-bp`=row$forward_start,`data-forward-end-bp`=row$forward_end,
          `data-reverse-start-bp`=row$reverse_start,`data-reverse-end-bp`=row$reverse_end,
          `data-window-state`=if(!product_visible) "outside" else if(product_clipped) "clipped" else "full",
          `data-pair-input`=pair_click_input,role=if(!is.null(pair_click_input)) "button" else NULL,
          tabindex=if(!is.null(pair_click_input)) "0" else NULL,
          `aria-label`=if(!is.null(pair_click_input)) paste("Inspect",row$pair_label) else NULL,
          `aria-pressed`=if(!is.null(pair_click_input)) if(row$pair_id %in% active_pair_id) "true" else "false" else NULL,
          tags$rect(x=0,y=y-25,width=width-12,height=if(compact) 72 else 52,rx=5,
            fill=if(row$pair_id %in% active_pair_id) "#eef2ff" else "transparent",
            `pointer-events`="all",style=if(!is.null(pair_click_input)) "cursor:pointer" else NULL),row_children)
        children <- append(children,list(row_group))
      } else children <- append(children, row_children)
    }
    if (marker_id == "COI" && !is.null(folmer_region)) {
      children <- append(children, list(
        tags$text(x = if(compact) left else left - 14, y = if(compact) folmer_y-13 else folmer_y + 4, `text-anchor` = if(compact) "start" else "end", fill = "#6b5720", `font-size` = 12, `font-weight` = 700, "Folmer region"),
        tags$line(x1 = scale_x(folmer_region$start), y1 = folmer_y, x2 = scale_x(folmer_region$end), y2 = folmer_y, stroke = "#d7a62b", `stroke-width` = 10, `stroke-linecap` = "round"),
        tags$text(x = if(compact) left else scale_x(folmer_region$end) + 10, y = if(compact) folmer_y+23 else folmer_y + 4, fill = "#6b5720", `font-size` = if(compact) 12 else 11, paste0(folmer_region$start, "–", folmer_region$end, " · ", folmer_region$length_bp, " bp"))
      ))
    } else {
      for (j in seq_len(nrow(gene_landmarks))) {
        lm <- gene_landmarks[j, ]
        if(lm$end<axis_min || lm$start>axis_max) next
        landmark_children <- list(
          tags$rect(x = scale_x(lm$start), y = folmer_y - 9, width = max(2, scale_x(lm$end) - scale_x(lm$start)), height = 18, fill = if (j %% 2) "#9aa6a0" else "#d7dad7"),
          tags$text(x = (scale_x(lm$start) + scale_x(lm$end)) / 2, y = folmer_y + 4, `text-anchor` = "middle", fill = "#17241d", `font-size` = if(compact) 12 else 10, `font-weight` = 700, lm$label)
        )
        if(is.null(domain) && !nrow(domain_landmarks)) children <- append(children,landmark_children)
        else children <- append(children,list(tags$g(class="primer-overview-landmark",
          `data-landmark-id`=lm$landmark_id,`data-landmark-type`=lm$landmark_type,
          `data-start-bp`=lm$start,`data-end-bp`=lm$end,
          tags$title(paste(lm$label,"·",lm$start,"–",lm$end,"bp ·",lm$source_note)),landmark_children)))
      }
    }
    if(nrow(domain_landmarks)) {
      children <- append(children,list(tags$text(x=if(compact) left else left-14,y=if(compact) domain_y-13 else domain_y+4,`text-anchor`=if(compact) "start" else "end",fill="#6b746f",`font-size`=if(compact) 12 else 11,"Domains")))
      for(j in seq_len(nrow(domain_landmarks))) {
        lm <- domain_landmarks[j,,drop=FALSE]
        if(lm$end<axis_min || lm$start>axis_max) next
        partial <- "boundary_status" %in% names(lm) && !is.na(lm$boundary_status) && grepl("partial",lm$boundary_status,ignore.case=TRUE)
        label <- paste0(lm$label,if(partial) " (partial)" else "")
        children <- append(children,list(tags$g(class="primer-overview-domain",
          `data-landmark-id`=lm$landmark_id,`data-parent-landmark`=if("parent_landmark" %in% names(lm)) lm$parent_landmark else NULL,
          `data-start-bp`=lm$start,`data-end-bp`=lm$end,
          tags$title(paste(label,"·",lm$start,"–",lm$end,"bp ·",lm$source_note)),
          tags$rect(x=scale_x(lm$start),y=domain_y-8,width=max(2,scale_x(lm$end)-scale_x(lm$start)),height=16,
            fill="#d7dad7",stroke="#9aa6a0",`stroke-dasharray`=if(partial) "3 2" else NULL),
          tags$text(x=(scale_x(lm$start)+scale_x(lm$end))/2,y=domain_y+4,`text-anchor`="middle",fill="#17241d",`font-size`=if(compact) 12 else 10,`font-weight`=700,label))))
      }
    }
    children <- append(children, list(
      tags$text(x = if(compact) left else left - 14, y = if(compact) coi_y-13 else coi_y + 4, `text-anchor` = if(compact) "start" else "end", fill = "#17241d", `font-size` = 12, `font-weight` = 700, marker$display_name),
      tags$line(x1 = left, y1 = coi_y, x2 = right, y2 = coi_y, stroke = "#17241d", `stroke-width` = 8, `stroke-linecap` = "round"),
      if(!compact) tags$text(x = right + 10, y = coi_y + 4, fill = "#48554e", `font-size` = 11, paste0(axis_min, "–", min(axis_max, reference_length), " shown · ", marker$reference_version, " is ", format(reference_length, big.mark = ","), " bp"))
    ))
    div(
      class = "primer-map-wrap",
      `data-layout`=if(compact) "compact" else if(!is.null(viewport_width)) "wide" else NULL,
      if(compact) tags$style(htmltools::HTML(".primer-map-wrap[data-layout='compact']{overflow-x:visible}.primer-map-wrap[data-layout='compact']>.primer-map{width:100%!important;min-width:0!important;max-width:100%;height:auto}.primer-overview-reference{font-size:12px;color:#526176;overflow-wrap:anywhere;margin:4px 0 12px}.primer-overview-measurements{border:1px solid #dbe3ee;border-radius:8px;padding:10px 12px}.primer-overview-measurements>summary{font-size:13px;font-weight:650;cursor:pointer}.primer-overview-measurement-scroll{max-width:100%;overflow:auto;margin-top:10px}.primer-overview-measurements table{border-collapse:collapse;white-space:nowrap;font-size:13px;width:100%}.primer-overview-measurements th,.primer-overview-measurements td{padding:9px 12px;border-bottom:1px solid #e2e8f0;text-align:left;vertical-align:top}.primer-overview-measurement-button{font:inherit;color:#4338ca;background:transparent;border:0;padding:0;text-align:left;cursor:pointer}.primer-overview-measurement-button[aria-pressed='true']{font-weight:750;text-decoration:underline}.primer-overview-measurement-button:focus-visible{outline:2px solid #4f46e5;outline-offset:3px}")),
      if(!is.null(pair_click_input)) tags$script(htmltools::HTML("if(!window.atlasPrimerOverviewClick){window.atlasPrimerOverviewClick=true;function atlasPrimerOverviewSelect(e){var row=e.target.closest('.primer-overview-row[data-pair-input]');if(!row||e.target.closest('a')||!window.Shiny)return;if(e.type==='keydown'&&e.key!=='Enter'&&e.key!==' ')return;if(e.type==='keydown')e.preventDefault();Shiny.setInputValue(row.dataset.pairInput,{pair:row.dataset.pairId,nonce:Date.now()},{priority:'event'});}document.addEventListener('click',atlasPrimerOverviewSelect);document.addEventListener('keydown',atlasPrimerOverviewSelect);}")),
      tags$svg(
        class = "primer-map", viewBox = paste(0, 0, width, height),
        role = "img",
        `aria-label` = paste0(
          "Interactive ", marker$display_name, " primer binding map showing base pairs ",
          axis_min, " through ", axis_max,
          ". Rows are numbered in their current displayed order."
        ),
        children
      ),
      if(compact) tags$p(class="primer-overview-reference",paste0(axis_min,"–",min(axis_max,reference_length)," bp shown · ",marker$coordinate_reference," · ",
        marker$reference_version," is ",format(reference_length,big.mark=",")," bp. Product lengths include both primers.")),
      if(compact) tags$details(class="primer-overview-measurements",tags$summary("Measurements and sources"),
        div(class="primer-overview-measurement-scroll",tags$table(
          tags$thead(tags$tr(lapply(c("Combination","Start bp","End bp","+ primers","− primers","Citation(s)"),function(label) tags$th(scope="col",label)))),
          tags$tbody(measurements))))
    )
}
