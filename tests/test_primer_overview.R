#!/usr/bin/env Rscript

# The shared overview consumes explicit data; do not source the app or a cohort.
source("R/primer_overview.R")
source("R/functions.R")

collect_tags <- function(x, name = NULL) {
  if (inherits(x, "shiny.tag")) {
    own <- if (is.null(name) || identical(x$name, name)) list(x) else list()
    return(c(own, unlist(lapply(x$children, collect_tags, name = name), recursive = FALSE)))
  }
  if (is.list(x)) return(unlist(lapply(x, collect_tags, name = name), recursive = FALSE))
  list()
}
has_attr <- function(tags, name, value) Filter(function(x) identical(as.character(x$attribs[[name]]), as.character(value)), tags)
number_attr <- function(tag, name) as.numeric(tag$attribs[[name]])
html <- function(tag) paste(as.character(tag), collapse = "\n")

geometry <- read.csv("data/derived/marker_pair_geometry.csv", stringsAsFactors = FALSE)
primers <- read.csv("data/primers.csv", stringsAsFactors = FALSE)
markers <- read.csv("data/catalog/markers.csv", stringsAsFactors = FALSE)
citations <- read.csv("data/citations.csv", stringsAsFactors = FALSE)
facets <- read.csv("data/catalog/primer_pair_facets.csv", stringsAsFactors = FALSE)
diet <- read.csv("data/catalog/primer_pair_diet_uses.csv", stringsAsFactors = FALSE)
folmer <- read.csv("data/derived/folmer_region.csv", stringsAsFactors = FALSE)
fish <- geometry[geometry$pair_id == "FISH_F2_R1", , drop = FALSE]
components <- primers[primers$pair_id == fish$pair_id, , drop = FALSE]
coi <- primer_overview_map(fish, components, markers[markers$marker_id == "COI", ],
  thermal_rows = pair_thermal_summary(components), citations = citations,
  primer_facets = facets, diet_uses = diet, arthropod_off_target_ids = fish$pair_id,
  folmer_region = folmer)
coi_html <- html(coi)
stopifnot(
  grepl("<circle cx=\"16\"", coi_html, fixed = TRUE),
  grepl("Map row 1: FishF2 + FishR1", coi_html, fixed = TRUE),
  grepl("Estimated primer Tm (50 mM salt): FishF2 52.6", coi_html, fixed = TRUE),
  grepl("Suggested pair Ta (starting estimate): 47.6–49.6 °C", coi_html, fixed = TRUE),
  grepl("Reported Ta: 52 °C", coi_html, fixed = TRUE),
  grepl("Arthropod reference: off-target binding check only", coi_html, fixed = TRUE),
  grepl("707 bp", coi_html, fixed = TRUE), grepl("655 bp", coi_html, fixed = TRUE),
  grepl("View 0–800 bp", coi_html, fixed = TRUE),
  grepl("+ primers", coi_html, fixed = TRUE), grepl("− primers", coi_html, fixed = TRUE),
  grepl("Folmer region", coi_html, fixed = TRUE),
  length(has_attr(collect_tags(coi, "a"), "href", "https://doi.org/10.1098/rstb.2005.1716")) == 1L
)
# Expected pixel endpoints are independently calculated from the frozen 0..800
# bp scale and the real FishF2/R1 reference coordinates (19..44 and 700..725).
blue <- has_attr(collect_tags(coi, "line"), "stroke", "#2563eb")[[1L]]
red <- has_attr(collect_tags(coi, "line"), "stroke", "#dc2626")[[1L]]
gray <- has_attr(collect_tags(coi, "line"), "stroke", "#9aa6b2")[[1L]]
stopifnot(
  abs(number_attr(blue, "x1") - (112 + 19 / 800 * 888)) < 1e-8,
  abs(number_attr(blue, "x2") - (112 + 44 / 800 * 888)) < 1e-8,
  abs(number_attr(red, "x1") - (112 + 700 / 800 * 888)) < 1e-8,
  abs(number_attr(red, "x2") - (112 + 725 / 800 * 888)) < 1e-8,
  identical(number_attr(gray, "x1"), number_attr(blue, "x2")),
  identical(number_attr(gray, "x2"), number_attr(red, "x1")),
  length(collect_tags(coi, "polygon")) == 2L,
  length(has_attr(collect_tags(coi, "text"), "x", "1090")) == 2L,
  length(has_attr(collect_tags(coi, "text"), "x", "1440")) == 2L
)

# 18S keeps its source-supplied landmark track and legacy coordinate extent.
ssu_row <- fish
ssu_row$pair_id <- "SSU_FIXTURE"
ssu_row$pair_label <- "SSU forward + reverse"
ssu_row$forward_start <- 100L
ssu_row$forward_end <- 120L
ssu_row$reverse_start <- 410L
ssu_row$reverse_end <- 430L
ssu_row$pair_start <- 100L
ssu_row$pair_end <- 430L
ssu_row$aligned_amplicon_bp <- 331L
ssu_row$targeted_region_bp <- 289L
ssu_row$placement_status <- "source_mapped_PR2_coordinates"
ssu_components <- components
ssu_components$pair_id <- ssu_row$pair_id
ssu_landmarks <- read.csv("data/catalog/marker_landmarks.csv", stringsAsFactors = FALSE)
ssu_landmarks <- ssu_landmarks[ssu_landmarks$marker_id == "18S", ]
ssu <- primer_overview_map(ssu_row, ssu_components, markers[markers$marker_id == "18S", ],
  landmarks = ssu_landmarks, focus = c(1, 1800), citations = citations)
stopifnot(
  grepl("View 0–1800 bp", html(ssu), fixed = TRUE),
  grepl("Placement provenance: source_mapped_PR2_coordinates", html(ssu), fixed = TRUE),
  all(vapply(ssu_landmarks$label, function(label) any(vapply(collect_tags(ssu, "text"), function(tag) identical(tag$children[[1L]], label), logical(1))), logical(1))),
  length(collect_tags(ssu, "polygon")) == 2L,
  !grepl("Folmer region", html(ssu), fixed = TRUE)
)

product <- function(id, fs, fe, rs, re) {
  data.frame(pair_id = id, pair_label = id, forward_start = fs, forward_end = fe,
    reverse_start = rs, reverse_end = re, pair_start = fs, pair_end = re,
    aligned_amplicon_bp = re - fs + 1L, targeted_region_bp = rs - fe - 1L,
    forward_primer = paste(id, "F"), reverse_primer = paste(id, "R"))
}
raw_rows <- rbind(product("full",150,175,375,400), product("left",50,125,475,500),
  product("tip_out",425,475,700,725), product("outside",50,60,80,90),
  product("body_only",50,75,475,500))
raw_marker <- data.frame(marker_id = "ITS_FUNGAL", display_name = "Fungal rDNA",
  reference_length = 1000L, coordinate_reference = "Fixture.raw.1", reference_version = "Fixture.raw.1")
raw_map <- primer_overview_map(raw_rows, data.frame(), raw_marker,
  focus = c(150,400), domain = c(100,450), active_pair_id = "tip_out", pair_click_input = "fungal-pair_click")
groups <- Filter(function(tag) !is.null(tag$attribs[["data-pair-id"]]), collect_tags(raw_map,"g"))
get_group <- function(id) has_attr(groups,"data-pair-id",id)[[1L]]
stopifnot(
  length(groups) == 5L,
  length(collect_tags(get_group("full"),"polygon")) == 2L,
  length(collect_tags(get_group("left"),"polygon")) == 1L,
  length(collect_tags(get_group("tip_out"),"polygon")) == 0L,
  length(collect_tags(get_group("outside"),"polygon")) == 0L,
  length(has_attr(collect_tags(get_group("tip_out"),"line"),"stroke","#2563eb")) == 1L,
  length(has_attr(collect_tags(get_group("outside"),"line"),"stroke","#2563eb")) == 0L,
  length(has_attr(collect_tags(get_group("outside"),"line"),"stroke","#dc2626")) == 0L,
  length(has_attr(collect_tags(get_group("outside"),"line"),"stroke","#9aa6b2")) == 0L,
  length(has_attr(collect_tags(get_group("body_only"),"line"),"stroke","#9aa6b2")) == 1L,
  get_group("left")$attribs[["data-product-bp"]] == 451L,
  get_group("left")$attribs[["data-start-bp"]] == 50L,
  get_group("left")$attribs[["data-end-bp"]] == 500L,
  grepl("451 bp",html(get_group("left")),fixed=TRUE),
  grepl("Partly in view",html(get_group("left")),fixed=TRUE),
  grepl("Outside this bp window",html(get_group("outside")),fixed=TRUE),
  get_group("tip_out")$attribs[["aria-pressed"]] == "true",
  get_group("tip_out")$attribs[["data-pair-input"]] == "fungal-pair_click",
  grepl("Shiny.setInputValue",html(raw_map),fixed=TRUE)
)

# Partial parent copies and individually mapped domains remain distinct. A
# domains-only reference must not acquire an invented parent-gene interval.
regions <- data.frame(landmark_id=c("28S#1","28S#2","D1","D2","missing"),
  label=c("28S copy 1","28S copy 2","D1","D2","D3"), start=c(100,300,150,350,NA),
  end=c(250,450,200,400,NA), landmark_type=c("gene","gene","divergent_domain","divergent_domain","divergent_domain"),
  boundary_status=c("partial","partial","full","full","missing"),
  parent_landmark=c(NA,NA,"28S#1","28S#2","28S#1"),
  source_note=c("NCBI partial gene","NCBI partial gene","Projected structural reference","Projected structural reference","Boundary unavailable"))
regions_map <- primer_overview_map(raw_rows[1L,],data.frame(),raw_marker,landmarks=regions,
  focus=c(150,400),domain=c(100,450))
gene_groups <- has_attr(collect_tags(regions_map,"g"),"class","primer-overview-landmark")
domain_groups <- has_attr(collect_tags(regions_map,"g"),"class","primer-overview-domain")
stopifnot(
  length(gene_groups)==2L, length(domain_groups)==2L,
  gene_groups[[1L]]$attribs[["data-end-bp"]]==250,
  gene_groups[[2L]]$attribs[["data-start-bp"]]==300,
  domain_groups[[1L]]$attribs[["data-parent-landmark"]]=="28S#1",
  domain_groups[[2L]]$attribs[["data-start-bp"]]==350,
  number_attr(collect_tags(gene_groups[[1L]],"rect")[[1L]],"y") != number_attr(collect_tags(domain_groups[[1L]],"rect")[[1L]],"y"),
  grepl("Projected structural reference",html(regions_map),fixed=TRUE)
)
domains_only <- primer_overview_map(raw_rows[1L,],data.frame(),raw_marker,
  landmarks=regions[3:5,],focus=c(150,400),domain=c(100,450))
stopifnot(length(has_attr(collect_tags(domains_only,"g"),"class","primer-overview-landmark"))==0L)

long_marker <- raw_marker
long_marker$reference_length <- 10000L
long_view <- primer_overview_map(product("335bp",1775,1796,2090,2109),data.frame(),long_marker,
  domain=c(1,9867),focus=c(1775,2109))
long_group <- has_attr(collect_tags(long_view,"g"),"data-pair-id","335bp")[[1L]]
stopifnot(long_group$attribs[["data-product-bp"]]==335L,
  length(has_attr(collect_tags(long_view,"text"),"y","40"))==9L,
  grepl("335 bp",html(long_group),fixed=TRUE),
  grepl("No primer pairs match",html(primer_overview_map(data.frame(),data.frame(),raw_marker)),fixed=TRUE))
unresolved <- raw_rows[1L,]
unresolved$reverse_start <- NA_real_
stopifnot(inherits(try(primer_overview_map(unresolved,data.frame(),raw_marker),silent=TRUE),"try-error"))
cat("Shared primer overview geometry, sources, landmarks, clipping and click identity passed.\n")
