#!/usr/bin/env Rscript
.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
source("R/data_layer.R")
source("R/coi_retention.R")
source("R/coi_promotion.R")
release_id <- Sys.getenv("RELEASE_ID", "")
if (!nzchar(release_id)) stop("RELEASE_ID is required")
candidate <- jsonlite::read_json(atlas_manifest_path("COI", getwd(), release_id), simplifyVector = FALSE)
proof <- coi_verify_promotion(candidate)
dir.create("logs", showWarnings = FALSE)
jsonlite::write_json(proof, file.path("logs", paste0("coi_pre_promotion_", release_id, ".json")), auto_unbox = TRUE, pretty = TRUE)
message("COI promotion retention passed against current hosted release ", proof$current_release_id,
        ": ", proof$retained_from_current, " retained records; ", proof$evidence)
