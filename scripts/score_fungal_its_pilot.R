#!/usr/bin/env Rscript

# Local, release-isolated ITS pilot. Requires build_fungal_its_pilot.py and
# MAFFT --addfragments --keeplength against FN812768.2 before this script.
root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(root, ".Rlib"), .libPaths()))
suppressPackageStartupMessages(library(PrimerMiner))
source(file.path(root, "R/functions.R"))

out <- file.path(root, "data/derived/its_fungal_pilot")
manifest <- read.csv(file.path(out, "panel_manifest.csv"), stringsAsFactors = FALSE)
aligned <- parse_fasta(file.path(out, "reference_aligned.fasta"))
if (length(aligned) != nrow(manifest) + 1L ||
    !setequal(names(aligned)[-1], manifest$panel_id) ||
    any(nchar(aligned) != 739L)) {
  stop("MAFFT alignment does not match the 48-record manifest and 739-bp reference")
}
aligned <- toupper(aligned)
reference <- clean_sequence(parse_fasta(file.path(
  root, "data/reference/its_fungal_reference_FN812768.2.fasta"
))[[1]])
if (aligned[[1]] != reference) stop("Reference row changed in MAFFT alignment")

pairs <- subset(read.csv(file.path(root, "data/catalog/primer_pairs.csv")), marker_id == "ITS_FUNGAL")
links <- read.csv(file.path(root, "data/catalog/pair_oligos.csv"))
oligos <- read.csv(file.path(root, "data/catalog/oligos.csv"))
geometry <- subset(read.csv(file.path(root, "data/derived/marker_pair_geometry.csv")),
                   marker_id == "ITS_FUNGAL")
pair_links <- merge(merge(subset(links, pair_id %in% pairs$pair_id), oligos,
                          by = "oligo_id", suffixes = c("", "_oligo")),
                    geometry[, c("pair_id", "forward_start", "forward_end",
                                 "reverse_start", "reverse_end", "placement_status")],
                    by = "pair_id")
pair_links$site_start <- ifelse(pair_links$direction == "forward",
                                pair_links$forward_start, pair_links$reverse_start)
pair_links$site_end <- ifelse(pair_links$direction == "forward",
                              pair_links$forward_end, pair_links$reverse_end)
pair_links$site_resolved <- pair_links$placement_status == "sequence_aligned"
sites <- pair_links[!duplicated(pair_links$oligo_id),
                    c("oligo_id", "sequence", "direction", "site_start", "site_end")]
sites$site_resolved <- vapply(sites$oligo_id, function(id) {
  any(pair_links$oligo_id == id & pair_links$site_resolved)
}, logical(1))
site_loci <- c(ITS_ITS1 = "18S", ITS_ITS1F = "18S",
               ITS_ITS3 = "5.8S", ITS_FITS7 = "5.8S", ITS_GITS7 = "5.8S",
               ITS_ITS86F = "5.8S", ITS_ITS4 = "28S",
               ITS_ITS4B = "28S", ITS_LR21 = "28S")
if (!setequal(sites$oligo_id, names(site_loci))) stop("Update site locus annotations")
sites$site_locus <- unname(site_loci[sites$oligo_id])
for (i in seq_len(nrow(sites))) {
  resolved <- pair_links[pair_links$oligo_id == sites$oligo_id[i] &
                           pair_links$site_resolved, , drop = FALSE]
  if (nrow(resolved)) {
    starts <- unique(resolved$site_start)
    ends <- unique(resolved$site_end)
    if (length(starts) != 1L || length(ends) != 1L) stop("Conflicting site coordinates")
    sites$site_start[i] <- starts
    sites$site_end[i] <- ends
  }
  if (sites$site_resolved[i] &&
      sites$site_end[i] - sites$site_start[i] + 1L != nchar(sites$sequence[i])) {
    stop("Primer length and reference window disagree: ", sites$oligo_id[i])
  }
}

site_rows <- list()
for (i in seq_len(nrow(sites))) {
  site <- sites[i, ]
  rows <- manifest[, c("panel_id", "accession", "species_hypothesis", "phylum",
                       "class_name", "marker_region", "length_bp", "source_doi",
                       "source_version")]
  rows$oligo_id <- site$oligo_id
  rows$direction <- site$direction
  rows$site_locus <- site$site_locus
  rows$reference_start <- if (site$site_resolved) site$site_start else NA_integer_
  rows$reference_end <- if (site$site_resolved) site$site_end else NA_integer_
  rows$binding_sequence_primer_oriented <- NA_character_
  rows$site_status <- "unresolved_reference_site"
  rows$reference_site_identity <- NA_real_
  rows$penalty <- NA_real_
  rows$mismatch_count <- NA_integer_
  rows$terminal_3_mismatches <- NA_integer_
  if (site$site_resolved) {
    for (j in seq_len(nrow(rows))) {
      sequence <- aligned[[rows$panel_id[j]]]
      observed <- strsplit(sequence, "", fixed = TRUE)[[1]]
      covered <- which(observed != "-")
      start <- site$site_start
      end <- site$site_end
      window <- observed[start:end]
      rows$site_status[j] <- if (!length(covered)) {
        "unaligned_fragment"
      } else if (min(covered) > start && max(covered) < end) {
        "both_ends_truncated"
      } else if (min(covered) > start) {
        "missing_5prime_site"
      } else if (max(covered) < end) {
        "missing_3prime_site"
      } else if (any(window == "-")) {
        "alignment_gap_at_site"
      } else if (any(!window %in% c("A", "C", "G", "T"))) {
        "ambiguous_base_at_site"
      } else {
        "scorable"
      }
      if (rows$site_status[j] == "scorable") {
        # MAFFT can force a fragment's terminal unrelated bases onto a short
        # reference anchor. Treat weakly homologous windows as unverified,
        # not as evidence for an extreme primer mismatch.
        reference_window <- strsplit(substr(reference, start, end), "", fixed = TRUE)[[1]]
        rows$reference_site_identity[j] <- mean(window == reference_window)
        if (rows$reference_site_identity[j] < 0.70) {
          rows$site_status[j] <- "homology_unverified_at_site"
          next
        }
        binding <- paste(window, collapse = "")
        rows$binding_sequence_primer_oriented[j] <- if (site$direction == "reverse") {
          reverse_complement(binding)
        } else binding
      }
    }
    available <- which(rows$site_status == "scorable")
    if (length(available)) {
      temporary <- tempfile(fileext = ".fasta")
      on.exit(unlink(temporary), add = TRUE)
      lines <- unlist(lapply(available, function(j) c(
        paste0(">", rows$panel_id[j]),
        substr(aligned[[rows$panel_id[j]]], site$site_start, site$site_end)
      )))
      if (length(available) == 1L) lines <- c(lines, ">duplicate_for_primerminer", lines[2])
      writeLines(lines, temporary)
      scored <- suppressMessages(PrimerMiner::evaluate_primer(
        alignment_imp = temporary, primer_sequ = site$sequence,
        start = 1L, stop = nchar(site$sequence),
        forward = site$direction == "forward", gap_NA = TRUE, N_NA = TRUE,
        mm_position = "Position_v1", mm_type = "Type_v1", adjacent = 2,
        sequ_names = TRUE
      ))
      scored <- scored[scored$Template != "duplicate_for_primerminer", , drop = FALSE]
      positions <- grep("^V[0-9]+$", names(scored), value = TRUE)
      score_matrix <- as.matrix(scored[, positions, drop = FALSE])
      lookup <- match(rows$panel_id[available], scored$Template)
      if (anyNA(lookup) || anyNA(scored$sum[lookup])) stop("PrimerMiner did not score a covered site")
      rows$penalty[available] <- scored$sum[lookup]
      rows$mismatch_count[available] <- rowSums(score_matrix[lookup, , drop = FALSE] > 0)
      terminal <- positions[as.integer(sub("V", "", positions)) <= 3L]
      rows$terminal_3_mismatches[available] <- rowSums(
        as.matrix(scored[lookup, terminal, drop = FALSE]) > 0
      )
    }
  }
  site_rows[[i]] <- rows
}
site_detail <- do.call(rbind, site_rows)
write.csv(site_detail, file.path(out, "site_scores.csv"), row.names = FALSE, na = "")

pair_rows <- list()
for (i in seq_len(nrow(pairs))) {
  pair <- pairs[i, ]
  link <- pair_links[pair_links$pair_id == pair$pair_id, ]
  forward_id <- link$oligo_id[link$direction == "forward"]
  reverse_id <- link$oligo_id[link$direction == "reverse"]
  if (length(forward_id) != 1L || length(reverse_id) != 1L) stop("Pair catalog malformed")
  f <- site_detail[site_detail$oligo_id == forward_id, ]
  r <- site_detail[site_detail$oligo_id == reverse_id, ]
  r <- r[match(f$panel_id, r$panel_id), ]
  scorable <- f$site_status == "scorable" & r$site_status == "scorable" &
    link$site_start[link$direction == "forward"] < link$site_start[link$direction == "reverse"]
  pair_rows[[i]] <- data.frame(
    f[, c("panel_id", "accession", "species_hypothesis", "phylum", "class_name",
          "marker_region", "source_doi", "source_version")],
    pair_id = pair$pair_id, forward_oligo_id = forward_id, reverse_oligo_id = reverse_id,
    pair_target_region = if (pair$pair_id %in% c("ITS1_ITS4", "ITS1F_ITS4", "ITS1F_ITS4B")) {
      "ITS1+5.8S+ITS2"
    } else if (pair$pair_id == "ITS1F_LR21") {
      "ITS1+5.8S+ITS2+LSU"
    } else {
      "ITS2"
    },
    forward_status = f$site_status, reverse_status = r$site_status,
    pair_status = ifelse(scorable, "scorable", paste0(
      "forward:", f$site_status, ";reverse:", r$site_status)),
    forward_penalty = f$penalty, reverse_penalty = r$penalty,
    pair_penalty = ifelse(scorable, f$penalty + r$penalty, NA_real_),
    pair_mismatch_count = ifelse(scorable, f$mismatch_count + r$mismatch_count, NA_integer_),
    stringsAsFactors = FALSE
  )
}
pair_detail <- do.call(rbind, pair_rows)
write.csv(pair_detail, file.path(out, "pair_scores.csv"), row.names = FALSE, na = "")

summarize <- function(data, group, status_column, penalty_column) {
  groups <- split(data, interaction(data[group], drop = TRUE, lex.order = TRUE))
  do.call(rbind, lapply(groups, function(x) {
    value <- x[[penalty_column]][x[[status_column]] == "scorable"]
    data.frame(x[1, group, drop = FALSE], n_panel = nrow(x),
               n_scorable = length(value),
               median_penalty_scorable = if (length(value)) median(value) else NA_real_,
               missing_reason_counts = paste(
                 paste(names(table(x[[status_column]][x[[status_column]] != "scorable"])),
                       as.integer(table(x[[status_column]][x[[status_column]] != "scorable"])),
                       sep = "="), collapse = " | "),
               stringsAsFactors = FALSE)
  }))
}
write.csv(summarize(site_detail, "oligo_id", "site_status", "penalty"),
          file.path(out, "site_summary.csv"), row.names = FALSE, na = "")
write.csv(summarize(pair_detail, "pair_id", "pair_status", "pair_penalty"),
          file.path(out, "pair_summary.csv"), row.names = FALSE, na = "")
write.csv(summarize(pair_detail, c("phylum", "pair_id"), "pair_status", "pair_penalty"),
          file.path(out, "clade_pair_summary.csv"), row.names = FALSE, na = "")
message("Wrote local ITS site and pair scores for ", nrow(manifest), " accessions")
