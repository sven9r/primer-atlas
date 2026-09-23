#!/usr/bin/env Rscript

project_root <- normalizePath(".", mustWork = TRUE)
library_dir <- file.path(project_root, ".Rlib")
dir.create(library_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library_dir, .libPaths()))

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  Ncpus = max(1, parallel::detectCores() - 1)
)

# Several COI enrichment scripts parse taxonomy and reference XML.  These are
# installed explicitly because they run after renv setup in the release CI.
required <- c("XML", "rentrez", "seqinr", "xml2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(
    missing,
    lib = library_dir,
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
}

commit <- "9c5e4fb7a4934f590f8df3f6300550a94b04e7f3"
package_dir <- file.path(project_root, "vendor", "PrimerMiner")
if (!file.exists(file.path(package_dir, "DESCRIPTION"))) {
  stop("The pinned vendored PrimerMiner source is missing.")
}

status <- system2(
  file.path(R.home("bin"), "R"),
  c("CMD", "INSTALL", "-l", shQuote(library_dir), shQuote(package_dir))
)
if (status != 0) stop("PrimerMiner installation failed.")

tool_version <- function(command, args = "--version") {
  executable <- Sys.which(command)
  if (!nzchar(executable)) return("not installed")

  output <- tryCatch(
    system2(executable, args, stdout = TRUE, stderr = TRUE),
    error = function(error) character()
  )
  if (length(output)) output[1] else "version unavailable"
}

tool_path <- function(command) {
  executable <- Sys.which(command)
  if (nzchar(executable)) executable else "not installed"
}

dir.create("data/provenance", recursive = TRUE, showWarnings = FALSE)
manifest <- data.frame(
  component = c("PrimerMiner", "PrimerMiner_source_commit", "vsearch", "mafft"),
  version = c(
    as.character(packageVersion("PrimerMiner")),
    commit,
    tool_version("vsearch"),
    tool_version("mafft")
  ),
  installation = c(
    "project-local .Rlib; pinned vendored source; BOLDconnectR removed for NCBI-only compatibility",
    "https://github.com/VascoElbrecht/PrimerMiner",
    tool_path("vsearch"),
    tool_path("mafft")
  ),
  stringsAsFactors = FALSE
)
write.table(
  manifest,
  "data/provenance/software_versions.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Installed PrimerMiner into", library_dir, "\n")
print(manifest, row.names = FALSE)
