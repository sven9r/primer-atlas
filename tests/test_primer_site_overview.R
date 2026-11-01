#!/usr/bin/env Rscript

# Geometry/interaction checks run with no app, alignment or mutable cohort.
renderer <- new.env(parent = baseenv())
source("R/primer_site_overview.R", local = renderer)
collect <- function(x, name = NULL) {
  if (inherits(x, "shiny.tag")) return(c(
    if (is.null(name) || identical(x$name, name)) list(x) else list(),
    unlist(lapply(x$children, collect, name = name), recursive = FALSE)))
  if (is.list(x)) return(unlist(lapply(x, collect, name = name), recursive = FALSE))
  list()
}
by_class <- function(tags, class) Filter(function(tag)
  !is.null(tag$attribs$class) && class %in% strsplit(tag$attribs$class, " ", fixed = TRUE)[[1L]], tags)
num <- function(tag, field) as.numeric(tag$attribs[[field]])
expect_error <- function(expr) inherits(try(force(expr), silent = TRUE), "try-error")

sites <- data.frame(primer_id = c("f1", "f2", "r1", "r2"),
  primer_name = c("Forward A", "Forward B", "Reverse A", "Reverse B"),
  direction = c("forward", "forward", "reverse", "reverse"),
  start = c(101, 111, 701, 711), end = c(120, 130, 720, 730),
  match_type = c("exact_site", "mismatch_site", "exact_site", "exact_site"),
  template_id = "raw-example", stringsAsFactors = FALSE)
genes <- data.frame(label = c("Annotated interval", "Narrow subregion"),
  start = c(1, 390), end = c(1000, 410))
saved <- serialize(sites, NULL)
for (width in c(320, 390, 565, 960, 1280)) {
  map <- renderer$primer_site_overview(sites, "raw-example", 1000,
    "COI", selected_forward = "f1", selected_reverse = "r1", viewport_width = width,
    annotations = genes)
  svg <- collect(map, "svg")[[1L]]
  groups <- by_class(collect(svg, "g"), "primer-site-choice")
  stopifnot(length(groups) == 4L,
    identical(as.numeric(strsplit(svg$attribs$viewBox, " ", fixed = TRUE)[[1L]])[3L], width),
    setequal(vapply(groups, function(g) g$attribs[["data-primer-id"]], character(1)), sites$primer_id),
    all(vapply(groups, function(g) g$attribs[["data-primer-input"]] == "primer_site_click", logical(1))))
  f1 <- Filter(function(g) g$attribs[["data-primer-id"]] == "f1", groups)[[1L]]
  f2 <- Filter(function(g) g$attribs[["data-primer-id"]] == "f2", groups)[[1L]]
  stopifnot(f1$attribs[["aria-pressed"]] == "true", f2$attribs[["aria-pressed"]] == "false",
    all(vapply(groups, function(g) g$attribs$role == "button" && g$attribs$tabindex == "0", logical(1))))
  for (g in groups) {
    row <- sites[sites$primer_id == g$attribs[["data-primer-id"]], ]
    shaft <- by_class(collect(g, "line"), "primer-site-shaft")[[1L]]
    tip <- by_class(collect(g, "polygon"), "primer-site-arrowhead")[[1L]]
    # Each inclusive 20-bp site occupies 20/1000 of the reference bar. Primer
    # width does not grow when the label is displaced by collision packing.
    stopifnot(abs(num(shaft, "x1") - (30 + (row$start - 1) / 1000 * (width - 54))) < 1e-8,
      abs(num(shaft, "x2") - (30 + row$end / 1000 * (width - 54))) < 1e-8,
      abs(num(shaft, "x2") - num(shaft, "x1") - 20 / 1000 * (width - 54)) < 1e-8,
      num(g, "data-start-bp") == row$start, num(g, "data-end-bp") == row$end,
      abs(as.numeric(strsplit(tip$attribs$points, " ", fixed = TRUE)[[1L]])[1L] -
        if (row$direction == "forward") num(shaft, "x2") else num(shaft, "x1")) < 1e-8)
  }
  # In each direction and lane, button rectangles have actual nonoverlapping
  # extents. SVG textLength fixes the rendered label width to those extents.
  packed <- split(groups, vapply(groups, function(g) paste(
    g$attribs[["data-direction"]], g$attribs[["data-label-lane"]]), character(1)))
  for (lane in packed) {
    lo <- vapply(lane, num, numeric(1), field = "data-label-left")
    hi <- vapply(lane, num, numeric(1), field = "data-label-right")
    o <- order(lo)
    stopifnot(all(lo >= 30), all(hi <= width - 24))
    if (length(o) > 1L) stopifnot(all(tail(lo[o], -1L) > head(hi[o], -1L)))
  }
}
stopifnot(identical(serialize(sites, NULL), saved))

# Dense adjacent sites, long names and duplicate chemical components get
# separate clickable positions without horizontal overflow or label overlap.
dense <- sites[rep(1L, 24L), ]
dense$primer_id <- paste0("dense-", seq_len(nrow(dense)))
dense$primer_name <- paste0("Very long forward oligo with a distinct label ", seq_len(nrow(dense)))
dense$start <- 495 + seq_len(nrow(dense)); dense$end <- dense$start + 19
map <- renderer$primer_site_overview(dense, "raw-example", 1000, viewport_width = 320)
groups <- by_class(collect(map, "g"), "primer-site-choice")
stopifnot(length(groups) == 24L,
  length(unique(vapply(groups, function(g) g$attribs[["data-label-lane"]], numeric(1)))) == 24L,
  all(vapply(groups, num, numeric(1), field = "data-label-left") >= 30),
  all(vapply(groups, num, numeric(1), field = "data-label-right") <= 296))

# Cropping leaves global bp coordinates intact. It does not move an actual
# extension endpoint to the window edge and invent a primer arrow there.
clipped <- sites
clipped$start <- c(90, 490, 100, 500); clipped$end <- c(110, 510, 120, 520)
crop <- renderer$primer_site_overview(clipped, "raw-example", 1000,
  viewport_width = 565, domain = c(100, 500), annotations = genes)
groups <- by_class(collect(crop, "g"), "primer-site-choice")
stopifnot(length(groups) == 4L, crop$attribs[["data-template-length"]] == 1000,
  crop$attribs[["data-domain-start"]] == 100, crop$attribs[["data-domain-end"]] == 500)
for (g in groups) {
  id <- g$attribs[["data-primer-id"]]
  row <- clipped[clipped$primer_id == id, ]
  endpoint <- if (row$direction == "forward") row$end else row$start
  visible <- endpoint >= 100 && endpoint <= 500
  stopifnot(length(by_class(collect(g, "polygon"), "primer-site-arrowhead")) == as.integer(visible),
    num(g, "data-start-bp") == row$start, num(g, "data-end-bp") == row$end)
}
ticks <- by_class(collect(crop, "g"), "primer-site-tick")
stopifnot(all(c(100, 500) %in% vapply(ticks, num, numeric(1), field = "data-bp")),
  !any(c(1, 1000) %in% vapply(ticks, num, numeric(1), field = "data-bp")))

missing <- sites[1L, ]; missing$start <- missing$end <- NA_real_
empty <- renderer$primer_site_overview(missing, "raw-example", 1000)
stopifnot(empty$attribs[["data-unavailable-sites"]] == 1,
  length(by_class(collect(empty, "g"), "primer-site-choice")) == 0L)
other_window <- renderer$primer_site_overview(sites, "raw-example", 1000, domain = c(800, 900))
stopifnot(other_window$attribs[["data-outside-sites"]] == 4L)
single <- renderer$primer_site_overview(data.frame(), "one-base", 1)
stopifnot(length(collect(single, "svg")) == 1L)

bad <- sites; bad$template_id[1L] <- "another-template"
stopifnot(expect_error(renderer$primer_site_overview(bad, "raw-example", 1000)))
bad <- sites; bad$end[1L] <- 1001
stopifnot(expect_error(renderer$primer_site_overview(bad, "raw-example", 1000)),
  expect_error(renderer$primer_site_overview(sites, "raw-example", 1000, domain = c(1, 1001))),
  expect_error(renderer$primer_site_overview(sites, "raw-example", 1000, viewport_width = 100)),
  expect_error(renderer$primer_site_overview(sites, "raw-example", 1000,
    selected_forward = "f1", selected_reverse = "f1")))

# Annotation-free templates do not acquire invented V or D region intervals.
stopifnot(length(by_class(collect(map, "g"), "primer-site-annotation")) == 0L,
  length(by_class(collect(crop, "g"), "primer-site-annotation")) == 2L)
message("Shared individual-primer overview: identities, raw-bp scale, packing and truthful cropping passed.")
