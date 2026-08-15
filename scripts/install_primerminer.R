#!/usr/bin/env Rscript

project_root <- normalizePath(".", mustWork = TRUE)
library_dir <- file.path(project_root, ".Rlib")
dir.create(library_dir, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library_dir, .libPaths()))

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  Ncpus = max(1, parallel::detectCores() - 1)
)

required <- c("XML", "rentrez", "seqinr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(
    missing,
    lib = library_dir,
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
}

commit <- "9c5e4fb7a4934f590f8df3f6300550a94b04e7f3"
work_dir <- tempfile("primerminer_install_")
dir.create(work_dir)
on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)

repo_dir <- file.path(work_dir, "PrimerMiner")
status <- system2(
  "git",
  c(
    "clone", "--quiet", "--no-checkout",
    "https://github.com/VascoElbrecht/PrimerMiner.git",
    shQuote(repo_dir)
  )
)
if (status != 0) stop("Could not clone PrimerMiner.")

status <- system2("git", c("-C", shQuote(repo_dir), "checkout", "--quiet", commit))
if (status != 0) stop("Could not check out the pinned PrimerMiner commit.")

package_dir <- file.path(repo_dir, "PrimerMiner")
description_path <- file.path(package_dir, "DESCRIPTION")
description <- readLines(description_path)

# BOLDconnectR 1.0 now pulls a large spatial-analysis dependency tree and is
# not needed for the NCBI-only workflow used here. PrimerMiner's NCBI download,
# clustering, plotting, and evaluate_primer functions do not call it. The
# compatibility patch is isolated to the temporary installation source.
description <- sub(
  "^Depends: BOLDconnectR, ",
  "Depends: ",
  description
)
writeLines(description, description_path)

status <- system2(
  file.path(R.home("bin"), "R"),
  c("CMD", "INSTALL", "-l", shQuote(library_dir), shQuote(package_dir))
)
if (status != 0) stop("PrimerMiner installation failed.")

dir.create("data/provenance", recursive = TRUE, showWarnings = FALSE)
manifest <- data.frame(
  component = c("PrimerMiner", "PrimerMiner_source_commit", "vsearch", "mafft"),
  version = c(
    as.character(packageVersion("PrimerMiner")),
    commit,
    system2("vsearch", "--version", stdout = TRUE, stderr = TRUE)[1],
    system2("mafft", "--version", stdout = TRUE, stderr = TRUE)[1]
  ),
  installation = c(
    "project-local .Rlib; BOLDconnectR dependency removed for NCBI-only compatibility",
    "https://github.com/VascoElbrecht/PrimerMiner",
    Sys.which("vsearch"),
    Sys.which("mafft")
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
