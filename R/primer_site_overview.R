# An individual-primer overview for any raw, ungapped reference sequence.
# The caller supplies verified/calculated positions and their evidence labels.
# This renderer neither searches sequences nor infers gene/variable boundaries.

.primer_site_pack <- function(rows, scale, left, right, label_field = "primer_name",
                              font_size = 13, gap = 10) {
  if (!nrow(rows)) return(data.frame(label = character(), label_x = numeric(),
    label_left = numeric(), label_right = numeric(), text_width = numeric(), lane = integer()))
  max_chars <- max(3L, floor((right - left - 22) / (font_size * .58)))
  label <- as.character(rows[[label_field]])
  too_long <- nchar(label) > max_chars
  label[too_long] <- paste0(substr(label[too_long], 1L, max_chars - 1L), "…")
  text_width <- pmax(font_size * .7, nchar(label) * font_size * .58)
  widths <- pmin(right - left, pmax(38, text_width + 18))
  centers <- scale((rows$start + rows$end - 1) / 2)
  xmin <- pmax(left, pmin(right - widths, centers - widths / 2))
  xmax <- xmin + widths
  lane <- integer(nrow(rows)); occupied <- numeric()
  for (i in order(xmin, xmax, label, seq_len(nrow(rows)))) {
    fits <- which(occupied + gap <= xmin[i])
    chosen <- if (length(fits)) fits[[1L]] else length(occupied) + 1L
    occupied[chosen] <- xmax[i]; lane[i] <- chosen
  }
  data.frame(label = label, label_x = (xmin + xmax) / 2,
    label_left = xmin, label_right = xmax, text_width = text_width, lane = lane,
    stringsAsFactors = FALSE)
}

primer_site_overview <- function(sites, template_id, template_length,
  marker_label = "Reference", click_input = "primer_site_click",
  selected_forward = character(), selected_reverse = character(),
  viewport_width = NULL, annotations = data.frame(), domain = NULL) {
  scalar_text <- function(x) length(x) == 1L && !is.na(x) && nzchar(x)
  if (!scalar_text(template_id) || !scalar_text(marker_label) || !scalar_text(click_input))
    stop("One template identity, marker label and click input are required.")
  if (!is.numeric(template_length) || length(template_length) != 1L ||
      !is.finite(template_length) || template_length < 1 || template_length != floor(template_length))
    stop("Template length must be the positive integer length of the raw sequence.")
  if (!is.null(viewport_width) && (!is.numeric(viewport_width) || length(viewport_width) != 1L ||
      !is.finite(viewport_width) || viewport_width < 320))
    stop("Overview width must be one finite pixel width of at least 320.")
  if (is.null(domain)) domain <- c(1, template_length)
  if (!is.numeric(domain) || length(domain) != 2L || any(!is.finite(domain)) ||
      any(domain != floor(domain)) || domain[1L] < 1 || domain[2L] > template_length ||
      domain[2L] < domain[1L]) stop("Choose an inclusive raw-bp window within the template.")
  if (length(intersect(selected_forward, selected_reverse)))
    stop("A selected primer must belong to one assigned direction at a time.")
  if (!is.data.frame(sites)) stop("Primer positions must be a data frame.")
  required <- c("primer_id", "primer_name", "direction", "start", "end")
  if (!nrow(sites) && !all(required %in% names(sites)))
    sites <- data.frame(primer_id = character(), primer_name = character(),
      direction = character(), start = numeric(), end = numeric())
  if (!all(required %in% names(sites))) stop("Primer identity, direction and raw bounds are required.")
  sites <- as.data.frame(sites, stringsAsFactors = FALSE)
  for (field in c("primer_id", "primer_name", "direction"))
    if (anyNA(sites[[field]]) || any(!nzchar(as.character(sites[[field]]))))
      stop("Primer identities, names and directions cannot be missing.")
  check_owner <- function(rows) {
    if ("template_id" %in% names(rows) &&
        any(is.na(rows$template_id) | rows$template_id != template_id))
      stop("Positions and annotations must belong to the chosen raw template.")
  }
  check_owner(sites)
  direction <- tolower(as.character(sites$direction))
  if (any(!direction %in% c("forward", "reverse"))) stop("Use explicit forward/reverse primer directions.")
  sites$source_direction <- direction
  if ("role" %in% names(sites)) {
    assigned <- tolower(as.character(sites$role))
    take <- !is.na(assigned) & nzchar(assigned)
    if (any(!assigned[take] %in% c("forward", "reverse"))) stop("Assigned roles must be forward or reverse.")
    direction[take] <- assigned[take]
  }
  sites$direction <- direction
  sites$start <- suppressWarnings(as.numeric(sites$start))
  sites$end <- suppressWarnings(as.numeric(sites$end))
  missing <- is.na(sites$start) | is.na(sites$end)
  invalid <- !missing & (!is.finite(sites$start) | !is.finite(sites$end) |
    sites$start != floor(sites$start) | sites$end != floor(sites$end) |
    sites$start < 1 | sites$end < sites$start | sites$end > template_length)
  if (any(invalid)) stop("Primer bounds must be inclusive integer positions within the raw template.")
  unavailable_n <- sum(missing)
  rows <- sites[!missing, , drop = FALSE]
  # One oligo can have several supported sites. Exact duplicate site rows do
  # not create duplicate click targets or artificially increase the count.
  rows <- rows[!duplicated(rows[, required, drop = FALSE]), , drop = FALSE]
  outside <- rows$end < domain[1L] | rows$start > domain[2L]
  outside_n <- sum(outside)
  rows <- rows[!outside, , drop = FALSE]
  rows$clipped <- rows$start < domain[1L] | rows$end > domain[2L]
  rows$match_type <- if ("match_type" %in% names(rows)) as.character(rows$match_type) else rep("Site evidence supplied", nrow(rows))
  rows$match_type[is.na(rows$match_type) | !nzchar(rows$match_type)] <- "Site evidence supplied"
  given_selected <- if ("selected" %in% names(rows))
    !is.na(rows$selected) & tolower(as.character(rows$selected)) %in% c("true", "1", "yes") else rep(FALSE, nrow(rows))
  rows$selected <- given_selected |
    (rows$direction == "forward" & rows$primer_id %in% selected_forward) |
    (rows$direction == "reverse" & rows$primer_id %in% selected_reverse)

  if (!is.data.frame(annotations)) stop("Annotations must be explicit raw intervals.")
  genes <- as.data.frame(annotations, stringsAsFactors = FALSE)
  if (nrow(genes)) {
    if (!all(c("label", "start", "end") %in% names(genes)))
      stop("Annotations need a label and explicit raw start/end.")
    check_owner(genes)
    genes$start <- suppressWarnings(as.numeric(genes$start))
    genes$end <- suppressWarnings(as.numeric(genes$end))
    if (anyNA(genes[, c("label", "start", "end")]) || any(!nzchar(genes$label)) ||
        any(!is.finite(genes$start) | !is.finite(genes$end) |
          genes$start != floor(genes$start) | genes$end != floor(genes$end) |
          genes$start < 1 | genes$end < genes$start | genes$end > template_length))
      stop("Annotations must have known inclusive bounds on this raw template.")
    withheld <- if ("boundary_status" %in% names(genes))
      genes$boundary_status %in% c("unresolved", "not_detected", "invalid_coordinates", "withheld", "extraction_mismatch") else rep(FALSE, nrow(genes))
    genes <- genes[!withheld & genes$end >= domain[1L] & genes$start <= domain[2L], , drop = FALSE]
  } else genes <- data.frame(label = character(), start = numeric(), end = numeric())

  width <- if (is.null(viewport_width)) 960 else viewport_width
  left <- 30; right <- width - 24
  # Interbase boundaries give every inclusive site its actual base width:
  # start..end occupies end-start+1 units. Tick labels mark base centres.
  lower_edge <- domain[1L] - 1; upper_edge <- domain[2L]
  scale <- function(boundary) left +
    (pmax(lower_edge, pmin(upper_edge, boundary)) - lower_edge) /
    (upper_edge - lower_edge) * (right - left)
  forwards <- rows[rows$direction == "forward", , drop = FALSE]
  reverses <- rows[rows$direction == "reverse", , drop = FALSE]
  fp <- .primer_site_pack(forwards, scale, left, right)
  rp <- .primer_site_pack(reverses, scale, left, right)
  gp <- .primer_site_pack(genes, scale, left, right, label_field = "label", font_size = 12)
  n_forward <- if (nrow(fp)) max(fp$lane) else 0L
  n_reverse <- if (nrow(rp)) max(rp$lane) else 0L
  n_gene <- if (nrow(gp)) max(gp$lane) else 0L
  gene_y <- 36 + n_forward * 46 + n_gene * 25 + 20
  axis_y <- gene_y + 28 + max(1L, n_reverse) * 46 + 34
  height <- axis_y + 44
  tags <- htmltools::tags
  fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  row_tags <- function(part, packing, role) lapply(seq_len(nrow(part)), function(i) {
    row <- part[i, , drop = FALSE]; packed <- packing[i, , drop = FALSE]
    selected <- isTRUE(row$selected)
    colour <- if (role == "forward") "#087e75" else "#6e4ec2"
    pale <- if (selected) if (role == "forward") "#def4ee" else "#eee7fb" else "#f5f7fa"
    arrow_y <- if (role == "forward") gene_y - n_gene * 25 - 24 - (packed$lane - 1L) * 46 else gene_y + 28 + 24 + (packed$lane - 1L) * 46
    label_y <- if (role == "forward") arrow_y - 17 else arrow_y + 24
    x1 <- scale(row$start - 1); x2 <- scale(row$end)
    endpoint <- if (role == "forward") row$end else row$start
    head_visible <- endpoint >= domain[1L] && endpoint <= domain[2L]
    tip <- if (role == "forward") x2 else x1
    head_width <- min(7, (x2 - x1) * .45)
    base <- if (role == "forward") tip - head_width else tip + head_width
    full_title <- paste0(row$primer_name, " · ", role, " · ", template_id, ":",
      fmt(row$start), "–", fmt(row$end), " (", fmt(row$end - row$start + 1),
      " bp) · ", row$match_type,
      if (row$clipped) " · Site extends outside this window" else "",
      if (selected) " · Selected" else " · Click to select")
    tags$g(class = paste("primer-site-choice", if (selected) "is-selected" else "",
        if (row$clipped) "is-clipped" else ""),
      role = "button", tabindex = "0", `aria-label` = full_title,
      `aria-pressed` = if (selected) "true" else "false",
      `data-primer-id` = row$primer_id, `data-primer-input` = click_input,
      `data-direction` = role, `data-source-direction` = row$source_direction,
      `data-template-id` = template_id, `data-start-bp` = row$start,
      `data-end-bp` = row$end, `data-three-prime-bp` = endpoint,
      `data-label-lane` = packed$lane, `data-label-left` = packed$label_left,
      `data-label-right` = packed$label_right,
      `data-clipped` = if (row$clipped) "true" else "false",
      `data-match-type` = row$match_type,
      style = "cursor:pointer;outline-offset:3px;",
      tags$title(full_title),
      tags$line(class = "primer-site-leader", x1 = packed$label_x,
        y1 = if (role == "forward") label_y + 6 else label_y - 13,
        x2 = (x1 + x2) / 2, y2 = arrow_y,
        stroke = if (selected) colour else "#b2beca", `stroke-width` = 1),
      tags$rect(class = "primer-site-label-box", x = packed$label_left,
        y = label_y - 14, width = packed$label_right - packed$label_left,
        height = 22, rx = 6, fill = pale, stroke = if (selected) colour else "#dce3eb",
        `stroke-width` = if (selected) 1.5 else 1),
      tags$text(class = "primer-site-label", x = packed$label_x, y = label_y + 1,
        `text-anchor` = "middle", `font-size` = 13, `font-weight` = if (selected) 650 else 500,
        fill = if (selected) colour else "#344459", textLength = packed$text_width,
        lengthAdjust = "spacingAndGlyphs", packed$label),
      tags$line(class = "primer-site-shaft", x1 = x1, x2 = x2, y1 = arrow_y,
        y2 = arrow_y, stroke = colour, `stroke-width` = if (selected) 4 else 3,
        `stroke-linecap` = "butt"),
      if (head_visible) tags$polygon(class = "primer-site-arrowhead", fill = colour,
        points = paste(c(tip, arrow_y, base, arrow_y - 4, base, arrow_y + 4), collapse = " ")))
  })
  gene_tags <- lapply(seq_len(nrow(genes)), function(i) {
    row <- genes[i, , drop = FALSE]; packed <- gp[i, , drop = FALSE]
    label_y <- gene_y - 12 - (packed$lane - 1L) * 25
    tags$g(class = "primer-site-annotation", `data-label` = row$label,
      `data-start-bp` = row$start, `data-end-bp` = row$end,
      `data-label-lane` = packed$lane, `data-label-left` = packed$label_left,
      `data-label-right` = packed$label_right,
      tags$title(paste0(row$label, " · ", template_id, ":", fmt(row$start), "–", fmt(row$end))),
      tags$rect(x = scale(row$start - 1), y = gene_y,
        width = scale(row$end) - scale(row$start - 1), height = 28,
        fill = c("#cddde9", "#d9d7ee", "#cde4de")[(i - 1L) %% 3L + 1L],
        stroke = "#7b8fa2", `stroke-width` = .8),
      tags$line(x1 = packed$label_x, x2 = scale((row$start + row$end - 1) / 2),
        y1 = label_y + 4, y2 = gene_y, stroke = "#8293a4", `stroke-width` = .8),
      tags$text(x = packed$label_x, y = label_y, `text-anchor` = "middle",
        `font-size` = 12, `font-weight` = 600, fill = "#334559",
        textLength = packed$text_width, lengthAdjust = "spacingAndGlyphs", packed$label))
  })
  ticks <- unique(c(domain[1L], pretty(domain, n = max(2, floor((right - left) / 85))), domain[2L]))
  ticks <- sort(ticks[ticks >= domain[1L] & ticks <= domain[2L] & ticks == floor(ticks)])
  # Keep exact edge coordinates without colliding labels near either edge.
  if (length(ticks) > 2L) {
    mid <- ticks[-c(1L, length(ticks))]
    mid <- mid[scale(mid - .5) - scale(ticks[1L] - .5) > 34 &
      scale(utils::tail(ticks, 1L) - .5) - scale(mid - .5) > 34]
    ticks <- c(ticks[1L], mid, utils::tail(ticks, 1L))
  }
  if (length(ticks) == 2L && diff(scale(ticks - .5)) < 34) ticks <- ticks[1L]
  axis_tags <- lapply(ticks, function(bp) tags$g(class = "primer-site-tick", `data-bp` = bp,
    tags$line(x1 = scale(bp - .5), x2 = scale(bp - .5), y1 = axis_y, y2 = axis_y + 5, stroke = "#b4c1ce"),
    tags$text(x = scale(bp - .5), y = axis_y + 21, `text-anchor` = "middle",
      `font-size` = 11, fill = "#607083", fmt(bp))))
  tags$div(class = "primer-site-overview", `data-template-id` = template_id,
    `data-template-length` = template_length, `data-domain-start` = domain[1L],
    `data-domain-end` = domain[2L], `data-plotted-sites` = nrow(rows),
    `data-unavailable-sites` = unavailable_n, `data-outside-sites` = outside_n,
    tags$div(class = "primer-site-overview-heading",
      tags$strong(marker_label), tags$span(paste0(" · ", template_id, " · ",
        fmt(domain[1L]), "–", fmt(domain[2L]), " bp"))),
    tags$svg(class = "primer-site-overview-svg", xmlns = "http://www.w3.org/2000/svg",
      width = "100%", viewBox = paste(0, 0, width, height), role = "group",
      `aria-label` = paste(marker_label, "individual primer positions on", template_id),
      style = "display:block;width:100%;max-width:100%;height:auto;font-family:Arial,sans-serif;",
      tags$rect(x = left, y = gene_y, width = right - left, height = 28,
        rx = 4, fill = "#edf1f5", stroke = "#a6b5c4", `stroke-width` = 1),
      gene_tags, row_tags(forwards, fp, "forward"), row_tags(reverses, rp, "reverse"),
      tags$line(x1 = left, x2 = right, y1 = axis_y, y2 = axis_y, stroke = "#c3ced9"),
      axis_tags,
      tags$text(x = (left + right) / 2, y = height - 1, `text-anchor` = "middle",
        `font-size` = 11, fill = "#607083", "Reference position (bp)")),
    if (!nrow(rows)) tags$p(class = "tiny primer-site-empty",
      "No positioned primers in this window. Choose another template or add primers below."),
    if (unavailable_n || outside_n) tags$p(class = "tiny primer-site-omissions",
      paste0(if (unavailable_n) paste0(unavailable_n, " site", if (unavailable_n != 1L) "s" else "", " without usable coordinates. ") else "",
        if (outside_n) paste0(outside_n, " site", if (outside_n != 1L) "s" else "", " outside this window.") else "")))
}
