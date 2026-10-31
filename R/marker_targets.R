# Marker choice and intended organism group are separate dimensions.  Registry
# status is a release-level property; it is not proof that a marker has a scored
# reference for the selected group (for example, COI arthropods versus fish).

atlas_target_marker_default <- function(group_id = "ALL", available_ids = NULL) {
  group_id <- if (length(group_id) == 1L && !is.na(group_id) && nzchar(group_id))
    as.character(group_id) else "ALL"
  preferred <- switch(group_id,
    ARTHROPODS = "COI", FISH = "12S_MT", BIRDS = "12S_MT",
    FUNGI = "ITS_FUNGAL", BACTERIA = "16S_PROK", ARCHAEA = "16S_PROK",
    PLANT_MICROBIOME = "16S_PROK", PROTISTS = "18S",
    PHYTOPLANKTON = "18S", ZOOPLANKTON = "18S", "COI")
  if (is.null(available_ids)) return(preferred)
  available_ids <- unique(as.character(available_ids))
  available_ids <- available_ids[!is.na(available_ids) & nzchar(available_ids)]
  if (preferred %in% available_ids) return(preferred)
  if (length(available_ids)) available_ids[[1L]] else NULL
}

atlas_target_marker_options <- function(group_id = "ALL", markers, links,
                                        include_planned = TRUE) {
  stopifnot(is.data.frame(markers), is.data.frame(links))
  if (!all(c("marker_id", "display_name", "status") %in% names(markers)))
    stop("Marker registry needs marker_id, display_name and status.")
  if (!all(c("marker_id", "group_id", "relationship") %in% names(links)))
    stop("Marker target links need marker_id, group_id and relationship.")
  if (anyDuplicated(markers$marker_id)) stop("Marker IDs must be unique.")
  group_id <- if (length(group_id) == 1L && !is.na(group_id) && nzchar(group_id))
    as.character(group_id) else "ALL"
  rows <- as.data.frame(markers, stringsAsFactors = FALSE)
  if (identical(group_id, "ALL")) {
    rows$relationship <- "available"
  } else {
    selected_links <- links[!is.na(links$group_id) & links$group_id == group_id, , drop = FALSE]
    selected_links <- selected_links[!duplicated(selected_links$marker_id), , drop = FALSE]
    rows <- rows[rows$marker_id %in% selected_links$marker_id, , drop = FALSE]
    rows$relationship <- as.character(selected_links$relationship[
      match(rows$marker_id, selected_links$marker_id)])
  }
  if (!isTRUE(include_planned))
    rows <- rows[rows$status %in% c("active", "pilot"), , drop = FALSE]

  # These are interface priorities, not a ranking of assay performance.
  priority <- switch(group_id,
    FISH = c("12S_MT", "COI", "18S", "28S", "16S_MT"),
    ARTHROPODS = c("COI", "18S", "28S", "16S_MT"),
    FUNGI = c("ITS_FUNGAL", "18S", "28S"),
    c(atlas_target_marker_default(group_id), markers$marker_id))
  rows <- rows[order(match(rows$marker_id, priority), rows$display_name,
                     na.last = TRUE), , drop = FALSE]
  rownames(rows) <- NULL
  labels <- c(COI = "COI", `12S_MT` = "12S · mitochondrial",
    `18S` = "18S · nuclear SSU", ITS_FUNGAL = "Fungal rDNA · ITS / SSU / LSU",
    `16S_PROK` = "16S · bacterial / archaeal", `16S_MT` = "16S · mitochondrial",
    `28S` = "28S · nuclear LSU")
  rows$choice_label <- unname(labels[rows$marker_id])
  missing_label <- is.na(rows$choice_label)
  rows$choice_label[missing_label] <- rows$display_name[missing_label]
  rows$readiness <- ifelse(rows$status == "active", "reference_release",
    ifelse(rows$status == "pilot", "catalog_pilot", "planned_release"))
  rows$readiness_label <- ifelse(rows$status == "active", "Published marker release",
    ifelse(rows$status == "pilot", "Catalog pilot", "Catalog / custom templates"))
  rows$group_id <- rep(group_id, nrow(rows))
  rows$is_default <- rows$marker_id %in%
    atlas_target_marker_default(group_id, rows$marker_id)
  rows
}

atlas_target_marker_choices <- function(options, show_readiness = FALSE) {
  if (!all(c("marker_id", "choice_label") %in% names(options)))
    stop("Use atlas_target_marker_options() to build marker choices.")
  labels <- options$choice_label
  if (isTRUE(show_readiness)) labels <- paste(labels, options$readiness_label, sep = " · ")
  stats::setNames(as.character(options$marker_id), labels)
}

atlas_target_marker_selection <- function(options, current = NULL,
                                           requested = NULL,
                                           group_changed = FALSE) {
  ids <- as.character(options$marker_id)
  valid <- function(x) length(x) == 1L && !is.na(x) && x %in% ids
  if (valid(requested)) return(as.character(requested))
  if (!isTRUE(group_changed) && valid(current)) return(as.character(current))
  defaults <- ids[!is.na(options$is_default) & options$is_default]
  if (length(defaults)) defaults[[1L]] else if (length(ids)) ids[[1L]] else NULL
}
