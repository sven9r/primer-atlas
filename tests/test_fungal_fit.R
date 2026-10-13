#!/usr/bin/env Rscript
# Synthetic frozen-evidence fixtures establish denominator and state invariants.
# Real V1/V3 checks independently reconcile the summaries to archived artifacts.
.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
source("R/fungal_fit.R")

expect_error <- function(expression) inherits(tryCatch(force(expression), error = identity), "error")
forward <- "AACCGGTTGACCTGACTGAC"
reverse_binding <- "CAGTCGATCGTAGCCTACGA"
sequence_for <- function(bp) paste0(forward, strrep("A", bp - 40L), reverse_binding,
                                   strrep("A", max(0L, 1000L - bp)))
exact_row <- function(accession, bp, pair_id = "P1") {
  data.frame(pair_id = pair_id, accession = accession,
    forward_start = 1, forward_end = 20, reverse_start = bp - 19, reverse_end = bp,
    amplicon_bp = bp, status = "scorable", forward_status = "exact_compatible_site",
    reverse_status = "exact_compatible_site", stringsAsFactors = FALSE)
}
accessions <- c(paste0("A.", 1:6), "B.1", "U.1")
sequences <- setNames(rep(sequence_for(500L), length(accessions)), accessions)
sequences[["A.1"]] <- paste0(substring(sequences[["A.1"]], 1, 99), "N", substring(sequences[["A.1"]], 101))
sequences[["A.2"]] <- sequence_for(501L)
sequences[["B.1"]] <- sequence_for(2000L)
taxonomy <- data.frame(accession = rev(accessions), organism = paste("Synthetic", rev(accessions)),
  phylum = c("unresolved", "Beta", rep("Alpha", 6)),
  class = c("unresolved", rep("SharedClass", 7)),
  order = c(NA, "OrderB", rep("OrderA", 6)), family = c(NA, rep("SharedFamily", 7)),
  taxid = seq_along(accessions), taxonomy_source_sha256 = rep(strrep("a", 64), length(accessions)),
  source_url = paste0("https://example.invalid/", rev(accessions)), stringsAsFactors = FALSE)
rows <- rbind(exact_row("A.1", 500), exact_row("A.2", 501), exact_row("A.3", 500),
              exact_row("A.4", 500), exact_row("A.5", 500), exact_row("B.1", 2000), exact_row("U.1", 500))
rows$status[rows$accession == "A.3"] <- "reverse_site_missing"
rows$reverse_status[rows$accession == "A.3"] <- "no_exact_binding_site"
rows$reverse_start[rows$accession == "A.3"] <- rows$reverse_end[rows$accession == "A.3"] <- NA
rows$status[rows$accession %in% c("A.4", "U.1")] <- "both_sites_missing"
rows$forward_status[rows$accession %in% c("A.4", "U.1")] <- "no_exact_binding_site"
rows$reverse_status[rows$accession %in% c("A.4", "U.1")] <- "no_exact_binding_site"
rows[rows$accession %in% c("A.4", "U.1"), c("forward_start", "forward_end", "reverse_start", "reverse_end")] <- NA
rows$status[rows$accession == "A.5"] <- "ambiguous_binding_site"
rows$forward_status[rows$accession == "A.5"] <- "ambiguous_binding_site"
rows[rows$accession == "A.5", c("forward_start", "forward_end")] <- NA
rows$amplicon_bp[rows$status != "scorable"] <- NA
pairs <- data.frame(pair_id = c("P1", "P2"), pair_label = c("Synthetic F + R", "Unevaluated fixture pair"),
  evidence_status = c("candidate_published_components", "documented_pair"), n_scorable = c(999L, 999L),
  source_citation = "Synthetic fixture", source_url = "https://example.invalid/pair", stringsAsFactors = FALSE)
snapshot <- list(manifest = list(id = "synthetic-v3"), sequences = sequences,
                 pairs = pairs, pair_sequences = rows, taxonomy = taxonomy)
fit <- fungal_fit_summary(snapshot, pair_ids = "P1", platform = "ILLUMINA")
stopifnot(identical(fit$ranks, c("phylum", "class", "order", "family")),
          fit$metadata$n_snapshot_templates == 8L, fit$metadata$n_selected_templates == 8L,
          nrow(fit$details) == 8L, all(fit$details$version_id == "synthetic-v3"))
alpha <- fit$summary[fit$summary$parent_path == "Alpha", , drop = FALSE]
beta <- fit$summary[fit$summary$parent_path == "Beta", , drop = FALSE]
unknown <- fit$summary[!fit$summary$taxon_resolved, , drop = FALSE]
stopifnot(nrow(alpha) == 1L, nrow(beta) == 1L, nrow(unknown) == 1L,
          alpha$taxon == beta$taxon, alpha$group_id != beta$group_id,
          alpha$n_panel == 6L, alpha$n_exact_products == 2L,
          alpha$n_partial_sites == 1L, alpha$n_unresolved_sites == 1L,
          alpha$n_no_exact_sites == 1L, alpha$n_unavailable == 1L,
          alpha$n_product_withheld == 0L, alpha$n_evidence_available == 5L,
          alpha$exact_product_fraction == 2 / 6,
          alpha$n_length_eligible == 1L, alpha$eligible_product_fraction == 1 / 2,
          alpha$eligible_panel_fraction == 1 / 6,
          beta$n_exact_products == 1L, beta$n_length_eligible == 0L,
          unknown$n_panel == 1L, unknown$n_taxonomy_unresolved == 1L,
          unknown$n_no_exact_sites == 1L, is.na(unknown$eligible_product_fraction))
partition <- c("n_exact_products", "n_partial_sites", "n_unresolved_sites", "n_no_exact_sites",
               "n_product_withheld", "n_unavailable")
stopifnot(all(rowSums(fit$summary[, partition]) == fit$summary$n_panel),
          all(fit$summary$taxon_key == fit$summary$group_id))
unavailable <- fit$details[fit$details$accession == "A.6", , drop = FALSE]
stopifnot(unavailable$fit_state == "evidence_unavailable", grepl("unavailable", unavailable$reason),
          is.na(unavailable$amplicon_bp), !unavailable$length_eligible)
ambiguous <- fit$details[fit$details$accession == "A.5", , drop = FALSE]
stopifnot(ambiguous$fit_state == "unresolved_sites", ambiguous$n_exact_sites == 1L,
          !ambiguous$length_eligible, grepl("Multiple", ambiguous$reason))
encoded <- fit$details[fit$details$accession == "A.1", , drop = FALSE]
stopifnot(encoded$fit_state == "exact_product", encoded$n_ambiguous_product_bases == 1L,
          encoded$has_ambiguous_product_bases, encoded$illumina_eligible,
          grepl("encoded", encoded$length_basis), grepl("unknown assembly gap", fit$metadata$length_note),
          encoded$taxonomy_source_sha256 == strrep("a", 64),
          encoded$source_url == "https://example.invalid/A.1")

# A platform changes planning eligibility, not the exact-evidence denominator.
pacbio <- fungal_fit_summary(snapshot, pair_ids = "P1", platform = "PACBIO")
stopifnot(identical(pacbio$summary$n_exact_products, fit$summary$n_exact_products),
          identical(pacbio$summary$n_panel, fit$summary$n_panel),
          pacbio$summary$n_length_eligible[pacbio$summary$parent_path == "Alpha"] == 2L,
          pacbio$summary$n_length_eligible[pacbio$summary$parent_path == "Beta"] == 1L)
subset <- fungal_fit_summary(snapshot, pair_ids = "P1", accessions = c("A.1", "A.6"))
stopifnot(nrow(subset$summary) == 1L, subset$summary$n_panel == 2L,
          subset$summary$n_exact_products == 1L, subset$summary$n_unavailable == 1L,
          subset$summary$exact_product_fraction == .5)
missing_pair <- fungal_fit_summary(snapshot, pair_ids = "P2")
stopifnot(sum(missing_pair$summary$n_panel) == 8L, sum(missing_pair$summary$n_unavailable) == 8L,
          sum(missing_pair$summary$n_no_exact_sites) == 0L,
          all(is.na(missing_pair$summary$eligible_product_fraction)))
stopifnot(nrow(fungal_fit_summary(snapshot, pair_ids = character())$details) == 0L,
          nrow(fungal_fit_summary(snapshot, accessions = character())$summary) == 0L,
          expect_error(fungal_fit_summary(snapshot, pair_ids = "not-a-frozen-pair")),
          expect_error(fungal_fit_summary(snapshot, accessions = "not-a-frozen-accession")),
          expect_error(fungal_fit_summary(snapshot, rank = "genus")),
          expect_error(fungal_fit_summary(snapshot, illumina_max = 2100)))

# Rank labels alone cannot collapse unrelated class/family paths.
family <- fungal_fit_summary(snapshot, pair_ids = "P1", rank = "family")
same_family <- family$summary[family$summary$taxon == "SharedFamily", , drop = FALSE]
stopifnot(nrow(same_family) == 2L, !anyDuplicated(same_family$group_id),
          any(grepl("Alpha › SharedClass › OrderA", same_family$parent_path, fixed = TRUE)),
          any(grepl("Beta › SharedClass › OrderB", same_family$parent_path, fixed = TRUE)))
no_family <- snapshot
no_family$taxonomy$family <- "unresolved"
no_family$taxonomy$order <- NA_character_
stopifnot(identical(fungal_fit_available_ranks(no_family), c("phylum", "class")),
          expect_error(fungal_fit_summary(no_family, rank = "family")))
missing_taxon <- snapshot
missing_taxon$taxonomy <- taxonomy[taxonomy$accession != "A.1", , drop = FALSE]
missing_taxon_fit <- fungal_fit_summary(missing_taxon, pair_ids = "P1", accessions = "A.1")
stopifnot(missing_taxon_fit$summary$n_panel == 1L,
          missing_taxon_fit$summary$n_exact_products == 1L,
          missing_taxon_fit$summary$n_taxonomy_unresolved == 1L)

# Known coordinates never upgrade a withheld or contradictory stored status.
withheld <- snapshot
withheld$pair_sequences$status[1] <- "outward_or_overlapping_sites"
withheld$pair_sequences$forward_start[1] <- 600
withheld$pair_sequences$forward_end[1] <- 619
withheld_fit <- fungal_fit_summary(withheld, pair_ids = "P1", accessions = "A.1")
stopifnot(withheld_fit$summary$n_product_withheld == 1L,
          withheld_fit$summary$n_exact_products == 0L, withheld_fit$summary$n_length_eligible == 0L)
corrupt <- snapshot
corrupt$pair_sequences$amplicon_bp[1] <- 499
corrupt_fit <- fungal_fit_summary(corrupt, pair_ids = "P1", accessions = "A.1")
stopifnot(corrupt_fit$summary$n_unavailable == 1L, corrupt_fit$summary$n_exact_products == 0L,
          grepl("inconsistent", corrupt_fit$details$reason))
old_geometry <- snapshot
old_geometry$pair_sequences$status[1] <- "reference_specific_low_confidence"
stopifnot(fungal_fit_summary(old_geometry, pair_ids = "P1", accessions = "A.1")$summary$n_unavailable == 1L)
duplicate <- snapshot
duplicate$pair_sequences <- rbind(rows, rows[1, ])
stopifnot(expect_error(fungal_fit_summary(duplicate)))

# Inclusive nested planning ceilings, including exact-supported long products
# beyond every selected sequencing ceiling.
bp <- c(452, 500, 501, 2000, 2001, 3500, 3501)
boundary <- snapshot
boundary$sequences <- setNames(vapply(bp, sequence_for, character(1)), paste0("T", bp))
boundary$taxonomy <- data.frame(accession = names(boundary$sequences), phylum = "P", class = "C")
boundary$pairs <- pairs[1, ]
boundary$pair_sequences <- do.call(rbind, lapply(seq_along(bp), function(i) exact_row(names(boundary$sequences)[i], bp[i])))
boundary_fit <- fungal_fit_summary(boundary, illumina_max = 452)
stopifnot(boundary_fit$summary$n_exact_products == 7L,
          boundary_fit$summary$n_illumina_eligible == 1L,
          boundary_fit$summary$n_pacbio_eligible == 4L,
          boundary_fit$summary$n_nanopore_eligible == 6L,
          boundary_fit$summary$n_length_eligible == 6L,
          boundary_fit$summary$eligible_product_fraction == 6 / 7,
          boundary_fit$summary$amplicon_max_bp == 3501)

# Published frozen data: preserve all21-template denominators for every pair,
# including the uncultured original reference and all nonmatching templates.
source("R/functions.R")
source("R/data_layer.R")
source("R/fungal_versions.R")
v3 <- fungal_read_snapshot("v3-c6fa87f256e1", remote = FALSE)
v3_fit <- fungal_fit_summary(v3)
stopifnot(nrow(v3_fit$details) == 18564L,
          sum(v3_fit$summary$n_exact_products) == 6858L,
          sum(v3_fit$summary$n_partial_sites) == 8636L,
          sum(v3_fit$summary$n_unresolved_sites) == 205L,
          sum(v3_fit$summary$n_no_exact_sites) == 2865L,
          sum(v3_fit$summary$n_unavailable) == 0L,
          all(rowSums(v3_fit$summary[, partition]) == v3_fit$summary$n_panel),
          all(tapply(v3_fit$summary$n_panel, v3_fit$summary$pair_id, sum) == 21L),
          all(v3_fit$summary$n_illumina_eligible <= v3_fit$summary$n_pacbio_eligible),
          all(v3_fit$summary$n_pacbio_eligible <= v3_fit$summary$n_nanopore_eligible),
          any(v3_fit$details$has_ambiguous_product_bases),
          all(v3_fit$details$version_id == v3$manifest$id))
v1 <- fungal_read_snapshot("v1-9e958f105c46", remote = FALSE)
v1_fit <- fungal_fit_summary(v1, pair_ids = "ITS3_ITS4")
stopifnot(v1_fit$summary$n_panel == 1L, v1_fit$summary$n_exact_products == 1L,
          v1_fit$summary$n_taxonomy_unresolved == 1L,
          !v1_fit$summary$taxon_resolved,
          identical(v1_fit$ranks, c("phylum", "class")),
          v1_fit$details$accession == "FN812768.2", v1_fit$details$amplicon_bp == 391)
message("Fungal-fit partitions, taxonomy paths, frozen denominators, planning thresholds and V1/V3 evidence passed.")
