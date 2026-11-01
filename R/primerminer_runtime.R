# Runtime adapter for the unchanged, pinned GPL-3 PrimerMiner evaluator.
# Attribution, source identities and the GPL-3 license accompany the component
# in R/primerminer_runtime/. No installed PrimerMiner/seqinr package is needed.
.atlas_pm_sources <- Filter(function(path) is.character(path) && length(path) == 1L &&
  basename(path) == "primerminer_runtime.R", lapply(sys.frames(), function(frame) get0("ofile", envir = frame, inherits = FALSE)))
.atlas_pm_source <- if (length(.atlas_pm_sources)) utils::tail(.atlas_pm_sources, 1L)[[1]] else NULL
atlas_evaluate_primer <- local({
  runtime_directory <- file.path(if (is.null(.atlas_pm_source)) file.path(getwd(), "R") else
    dirname(normalizePath(.atlas_pm_source, mustWork = TRUE)), "primerminer_runtime")
  expected <- c(evaluate_primer.R = "531662139dd16833c26bd1a9ee6aaf69",
    Position_v1.csv = "53297a2f513a4037972ce9488ff90360", Type_v1.csv = "16f380370c754cb5712a08227746196c")
  evaluator <- NULL
  fasta_cache <- new.env(parent = emptyenv())

  read_fasta_adapter <- function(file, seqonly = FALSE) {
    signature <- c(normalizePath(file, mustWork = TRUE), unname(tools::md5sum(file)))
    if (!identical(fasta_cache$signature, signature)) {
      lines <- readLines(file, warn = FALSE)
      starts <- which(startsWith(lines, ">"))
      if (!length(starts)) stop("Primer scoring alignment has no FASTA records")
      ends <- c(starts[-1] - 1L, length(lines))
      sequences <- vapply(seq_along(starts), function(index) {
        if (ends[index] <= starts[index]) return("")
        paste0(gsub("[[:space:]]", "", lines[seq.int(starts[index] + 1L, ends[index])]), collapse = "")
      }, character(1))
      identifiers <- sub("[[:space:]].*$", "", substring(lines[starts], 2L))
      if (any(!nzchar(identifiers)) || anyDuplicated(identifiers)) stop("Primer scoring alignment needs unique nonempty FASTA identifiers")
      if (any(!nzchar(sequences)) || length(unique(nchar(sequences))) != 1L) stop("Primer scoring alignment is empty or has unequal row widths")
      names(sequences) <- identifiers
      fasta_cache$signature <- signature
      fasta_cache$sequences <- sequences
    }
    sequences <- fasta_cache$sequences
    # PrimerMiner uses strings in seqonly mode and only names in the other
    # mode. First-token identifiers match seqinr's default FASTA-name handling.
    if (seqonly) as.list(unname(sequences)) else as.list(sequences)
  }

  function(...) {
    if (is.null(evaluator)) {
      paths <- file.path(runtime_directory, names(expected))
      if (any(!file.exists(paths)) || !identical(unname(tools::md5sum(paths)), unname(expected))) {
        stop("Pinned PrimerMiner runtime evaluator/tables are missing or changed", call. = FALSE)
      }
      isolated <- new.env(parent = baseenv())
      isolated$read.csv <- utils::read.csv
      isolated$write.csv <- utils::write.csv
      isolated$read.fasta <- read_fasta_adapter
      isolated$system.file <- function(..., package) {
        if (!identical(package, "PrimerMiner")) stop("Unexpected runtime table provider")
        runtime_directory
      }
      sys.source(paths[1], envir = isolated, keep.source = FALSE)
      evaluator <<- isolated$evaluate_primer
    }
    evaluator(...)
  }
})
rm(.atlas_pm_sources, .atlas_pm_source)
