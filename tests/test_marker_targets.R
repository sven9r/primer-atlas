#!/usr/bin/env Rscript

source("R/marker_targets.R")
markers <- read.csv("data/catalog/markers.csv", stringsAsFactors = FALSE)
links <- read.csv("data/catalog/marker_organism_groups.csv", stringsAsFactors = FALSE)
fish <- atlas_target_marker_options("FISH", markers, links)
stopifnot(
  identical(fish$marker_id, c("12S_MT", "COI", "18S", "28S", "16S_MT")),
  identical(atlas_target_marker_selection(fish), "12S_MT"),
  identical(atlas_target_marker_selection(fish, current = "COI"), "COI"),
  identical(atlas_target_marker_selection(fish, current = "COI", group_changed = TRUE), "12S_MT"),
  identical(atlas_target_marker_selection(fish, requested = "18S", group_changed = TRUE), "18S"),
  identical(atlas_target_marker_selection(fish, requested = "ITS_FUNGAL"), "12S_MT"),
  fish$status[fish$marker_id == "12S_MT"] == "planned",
  fish$readiness[fish$marker_id == "12S_MT"] == "planned_release",
  fish$readiness[fish$marker_id == "18S"] == "catalog_pilot",
  !any(c("fit_ready", "scoring_ready") %in% names(fish))
)
stopifnot(
  identical(unname(atlas_target_marker_choices(fish)), fish$marker_id),
  any(grepl("Catalog / custom templates", names(atlas_target_marker_choices(fish, TRUE)), fixed = TRUE)),
  identical(atlas_target_marker_options("FISH", markers, links, FALSE)$marker_id,
    c("COI", "18S"))
)

arthropods <- atlas_target_marker_options("ARTHROPODS", markers, links)
fungi <- atlas_target_marker_options("FUNGI", markers, links)
all_markers <- atlas_target_marker_options("ALL", markers, links)
stopifnot(
  identical(atlas_target_marker_selection(arthropods), "COI"),
  identical(atlas_target_marker_selection(fungi), "ITS_FUNGAL"),
  identical(fungi$marker_id, c("ITS_FUNGAL", "18S", "28S")),
  setequal(all_markers$marker_id, markers$marker_id),
  identical(atlas_target_marker_selection(all_markers), "COI"),
  identical(atlas_target_marker_default("PROTISTS"), "18S"),
  identical(atlas_target_marker_default("BACTERIA"), "16S_PROK"),
  identical(atlas_target_marker_default("FISH", c("COI", "18S")), "COI")
)
for (group_id in c("METAZOA", "ARTHROPODS", "FISH", "BIRDS", "ZOOPLANKTON")) {
  rows <- atlas_target_marker_options(group_id, markers, links)
  stopifnot(all(c("18S", "28S") %in% rows$marker_id),
    rows$relationship[rows$marker_id == "28S"] == "secondary")
}
birds <- atlas_target_marker_options("BIRDS", markers, links)
stopifnot(all(c("12S_MT","COI","18S","28S") %in% birds$marker_id))
stopifnot(all(c("18S", "28S", "16S_MT") %in% arthropods$marker_id),
  !anyDuplicated(paste(links$marker_id, links$group_id)))

unknown <- atlas_target_marker_options("UNKNOWN", markers, links)
stopifnot(nrow(unknown) == 0L, is.null(atlas_target_marker_selection(unknown)),
  length(atlas_target_marker_choices(unknown)) == 0L,
  is.null(atlas_target_marker_default("FISH", character())),
  identical(atlas_target_marker_default(NULL), "COI"))

# A target switch changes only options/defaults.  The frozen registry and its
# existing marker identities, coordinate references and release statuses survive.
original_markers <- read.csv("data/catalog/markers.csv", stringsAsFactors = FALSE)
stopifnot(identical(markers, original_markers),
  identical(fish$coordinate_reference[fish$marker_id == "18S"],
    markers$coordinate_reference[markers$marker_id == "18S"]))
message("Target-group marker choices, explicit readiness and defaults passed.")
