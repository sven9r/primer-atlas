source("renv/activate.R")
local_atlas_library <- file.path(getwd(), ".Rlib")
if (dir.exists(local_atlas_library)) {
  .libPaths(c(local_atlas_library, .libPaths(), file.path(R.home(), "library")))
}
