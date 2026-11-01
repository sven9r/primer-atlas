#!/usr/bin/env Rscript
.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
source("R/functions.R")
suppressPackageStartupMessages(library(PrimerMiner))
stopifnot(as.character(packageVersion("PrimerMiner")) == "0.22")
fixture <- tempfile("primerminer_runtime_")
dir.create(fixture)
on.exit(unlink(fixture, recursive = TRUE), add = TRUE)

# The source and default weighting tables are unmodified pinned copies.
for (name in c("evaluate_primer.R", "Position_v1.csv", "Type_v1.csv")) {
  source_path <- if (name == "evaluate_primer.R") file.path("vendor/PrimerMiner/R", name) else
    file.path("vendor/PrimerMiner/inst", name)
  stopifnot(identical(readBin(source_path, "raw", file.info(source_path)$size),
    readBin(file.path("R/primerminer_runtime", name), "raw", file.info(source_path)$size)))
}

write_alignment <- function(path, sequences) {
  lines <- unlist(lapply(seq_along(sequences), function(index) {
    c(paste0(">", names(sequences)[index], " full descriptive header"),
      substr(sequences[index], 1L, 7L), substr(sequences[index], 8L, nchar(sequences[index])))
  }), use.names = FALSE)
  writeLines(lines, path)
}
evaluate <- function(engine, path, primer, forward = TRUE, adjacent = 2, gap_NA = TRUE, N_NA = TRUE) {
  suppressMessages(engine(alignment_imp = path, primer_sequ = primer, start = 1L, stop = nchar(primer),
    forward = forward, adjacent = adjacent, gap_NA = gap_NA, N_NA = N_NA, sequ_names = TRUE))
}
assert_parity <- function(path, primer, ...) {
  expected <- evaluate(PrimerMiner::evaluate_primer, path, primer, ...)
  observed <- evaluate(atlas_evaluate_primer, path, primer, ...)
  stopifnot(isTRUE(all.equal(observed, expected, tolerance = 0, check.attributes = TRUE)))
  observed
}

# Exhaust all upstream primer/template IUPAC symbols, including inosine,
# ambiguous mixtures, gaps/N and both orientations; headers and line wrapping
# exercise the native FASTA adapter rather than only the numerical function.
symbols <- strsplit("ACGTRYSWKMBDHVNI-", "", fixed = TRUE)[[1]]
path <- file.path(fixture, "iupac.fasta")
write_alignment(path, setNames(vapply(symbols, strrep, character(1), 15L), paste0("template", seq_along(symbols))))
for (primer_base in setdiff(symbols, "-")) {
  for (forward in c(TRUE, FALSE)) assert_parity(path, strrep(primer_base, 15L), forward = forward)
}
for (adjacent in c(1, 2, 3)) {
  assert_parity(path, "ACGTRYSWKMBDHVI", adjacent = adjacent)
  assert_parity(path, "ACGTRYSWKMBDHVI", forward = FALSE, adjacent = adjacent, gap_NA = FALSE, N_NA = FALSE)
}

# Preserve the table's unavailable distances beyond 30 rather than inventing
# extrapolated values, including the upstream perfect-match NA behavior.
long_path <- file.path(fixture, "long.fasta")
write_alignment(long_path, c(match = strrep("A", 31L), mismatch = strrep("C", 31L)))
long <- assert_parity(long_path, strrep("A", 31L))
stopifnot(all(is.na(long$V31)), all(is.na(long$sum)))

# A cached FASTA is invalidated by its actual content, even when path and width
# are unchanged. Invalid alignments are withheld before matrix recycling.
cached <- file.path(fixture, "cached.fasta")
write_alignment(cached, c(one = strrep("A", 15L), two = strrep("A", 15L)))
stopifnot(all(assert_parity(cached, strrep("A", 15L))$sum == 0))
write_alignment(cached, c(one = strrep("A", 15L), two = paste0(strrep("A", 14L), "C")))
stopifnot(assert_parity(cached, strrep("A", 15L))$sum[2] > 0)
bad <- file.path(fixture, "bad.fasta")
for (lines in list(c(">one", "AAAA", ">two", "AAA"), c(">one first", "AAAA", ">one second", "AAAA"),
                   c(">one", ">two", "AAAA"))) {
  writeLines(lines, bad)
  stopifnot(inherits(tryCatch(evaluate(atlas_evaluate_primer, bad, "AAAA"), error = identity), "error"))
}

# The application bundle can score both order-wide and expanded custom pairs
# in a fresh process that cannot import PrimerMiner or seqinr.
bundle <- file.path(fixture, "bundle")
dir.create(bundle)
stopifnot(file.copy("R", bundle, recursive = TRUE))
child <- file.path(fixture, "without-package.R")
writeLines(r"(
args <- commandArgs(trailingOnly = TRUE)
setwd(dirname(args[1]))
.libPaths(.Library, include.site = FALSE)
stopifnot(!requireNamespace('PrimerMiner', quietly = TRUE), !requireNamespace('seqinr', quietly = TRUE))
entry <- file.path(args[1], 'entry.R')
writeLines(paste0('source(', encodeString(file.path(args[1], 'R/functions.R'), quote = '"'), ')'), entry)
source(entry)
path <- 'OrderA_COX1_full_reference_aligned.fasta'
sequences <- c(strrep('A', 60L), strrep('A', 60L),
  paste0(strrep('A', 14L), 'C', strrep('A', 45L)), paste0('N', strrep('A', 59L)))
identifiers <- c('NC_001322.1_COX1_reference', 'AB123456.1', 'AB123457.1', 'AB123458.1')
lines <- unlist(lapply(seq_along(sequences), function(i) c(paste0('>', identifiers[i], ' header'), sequences[i])), use.names = FALSE)
writeLines(lines, path)
primers <- data.frame(pair_id = 'CUSTOM', pair_label = 'Runtime fixture', primer_name = c('F','R'),
  direction = c('forward','reverse'), sequence = c(strrep('A',15L),strrep('T',15L)))
geometry <- data.frame(forward_start = 1L, forward_end = 15L, reverse_start = 31L, reverse_end = 45L)
result <- score_primer_pair_primerminer(path, primers, geometry)
stopifnot(nrow(result$pair_template_scores) == 3L, result$order_scores$n_templates == 3L,
  result$order_scores$n_pair_scorable == 2L, sum(result$pair_template_scores$perfect_match_pair) == 1L,
  all(result$pair_template_scores$order == 'OrderA'),
  any(result$primer_position_scores$unavailable_reason == 'ambiguous N'))
expanded <- 'expanded.fasta'
writeLines(lines[-c(1L,2L)], expanded)
taxonomy <- data.frame(accession = identifiers[-1])
expanded_result <- score_primer_pair_expanded_reference(expanded, taxonomy, primers, geometry, 0L)
stopifnot(nrow(expanded_result$scores) == 3L, sum(expanded_result$scores$pair_scorable) == 2L)
for (accession in taxonomy$accession) {
  expected <- result$pair_template_scores$pair_penalty[result$pair_template_scores$template == accession]
  actual <- expanded_result$scores$pair_penalty[expanded_result$scores$accession == accession]
  stopifnot(isTRUE(all.equal(actual, expected)))
}
stopifnot(!'PrimerMiner' %in% loadedNamespaces(), !'seqinr' %in% loadedNamespaces())
cat('Fresh runtime without PrimerMiner/seqinr scored both custom pair paths.\n')
)", child)
status <- system2(file.path(R.home("bin"), "Rscript"),
  c("--vanilla", shQuote(child), shQuote(normalizePath(bundle))))
stopifnot(status == 0L)
message("Pinned PrimerMiner runtime matches installed 0.22 across IUPAC, orientation, adjacency, missingness and long-primer fixtures.")
