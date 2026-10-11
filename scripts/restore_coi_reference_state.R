#!/usr/bin/env Rscript
root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(root, ".Rlib"), .libPaths()))
source("R/data_layer.R")
source("R/coi_retention.R")
state <- atlas_read_manifest("COI", root, release_id = {
  id <- Sys.getenv("COI_PREVIOUS_RELEASE_ID", "")
  if (nzchar(id)) id else NULL
})
result <- coi_restore_retained_state(root, state,
  allow_bootstrap = identical(tolower(Sys.getenv("COI_ALLOW_BOOTSTRAP", "false")), "true"))
if (isTRUE(result$restored)) {
  message("Restored ", result$retained_records,
    " retained COI sequences from ", result$previous_release_id)
} else {
  membership <- coi_recover_published_membership(state$manifest, root)
  message("Explicit COI bootstrap: recovered ", nrow(membership),
    " published arthropod accession identities; original raw sequence bytes remain unavailable")
}
