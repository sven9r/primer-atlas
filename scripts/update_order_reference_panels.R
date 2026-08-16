#!/usr/bin/env Rscript

# Append-only COI reference acquisition. The initial target is 1,000 retained
# sequences per arthropod target group. Each later successful monthly release
# requests ceil(previous_n * 1.01) and adds only previously unseen accessions.

root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(root, ".Rlib"), .libPaths()))
suppressPackageStartupMessages(library(rentrez))
source(file.path(root, "R", "functions.R"))

orders <- read.csv("data/orders_22.csv", stringsAsFactors = FALSE)
orders <- orders[orders$order != "Bivalvia", , drop = FALSE]
baseline <- as.integer(Sys.getenv("COI_ORDER_BASELINE", "1000"))
growth <- as.numeric(Sys.getenv("COI_MONTHLY_GROWTH", "0.01"))
release_id <- Sys.getenv("RELEASE_ID", format(Sys.Date(), "%Y-%m-%d"))
dry_run <- identical(tolower(Sys.getenv("DRY_RUN", "false")), "true")
if (!is.finite(baseline) || baseline < 1000L) stop("COI_ORDER_BASELINE must be at least 1000.")
if (!is.finite(growth) || growth < 0.01) stop("COI_MONTHLY_GROWTH must be at least 0.01.")

ledger_path <- "data/provenance/coi_reference_accession_ledger.csv"
ledger <- if (file.exists(ledger_path)) read.csv(ledger_path, stringsAsFactors = FALSE) else data.frame(
  marker_id = character(), target_group = character(), accession = character(),
  entrez_uid = character(), first_release_id = character(), added_at_utc = character()
)
dir.create(dirname(ledger_path), recursive = TRUE, showWarnings = FALSE)

accession_from_header <- function(x) {
  token <- sub(" .*", "", sub("^>", "", x))
  token <- sub("^_R_", "", token)
  sub("[|].*$", "", token)
}

query_suffix <- paste0(
  "[Organism] AND (COI[Gene] OR CO1[Gene] OR COX1[Gene] OR COXI[Gene])",
  " AND 1200:1800[Sequence Length] NOT COII[Title]"
)
qa <- vector("list", nrow(orders))

for (i in seq_len(nrow(orders))) {
  group <- orders$order[i]
  existing <- ledger[ledger$target_group == group, , drop = FALSE]
  previous_n <- length(unique(existing$accession[nzchar(existing$accession)]))
  target_n <- if (previous_n < baseline) baseline else ceiling(previous_n * (1 + growth))
  requested_new <- max(0L, target_n - previous_n)
  query <- paste0(group, query_suffix)
  message(group, ": retained=", previous_n, ", target=", target_n, ", append=", requested_new)

  search <- rentrez::entrez_search(
    db = "nuccore", term = query,
    retmax = max(target_n * 3L, baseline * 3L), sort = "pub date"
  )
  candidate_uids <- setdiff(as.character(search$ids), existing$entrez_uid)
  chosen_uids <- head(candidate_uids, requested_new)
  added_accessions <- character()
  added_uids <- character()

  raw_dir <- file.path("data", "primerminer", "raw", group)
  raw_path <- file.path(raw_dir, paste0(group, "_GB.fasta"))
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  retained_sequences <- if (file.exists(raw_path)) parse_fasta(raw_path) else character()
  retained_accessions <- if (length(retained_sequences)) accession_from_header(names(retained_sequences)) else character()

  if (length(chosen_uids) && !dry_run) {
    for (batch_start in seq(1L, length(chosen_uids), by = 100L)) {
      batch <- chosen_uids[batch_start:min(batch_start + 99L, length(chosen_uids))]
      fetched <- rentrez::entrez_fetch(db = "nuccore", id = batch, rettype = "fasta", retmode = "text")
      scratch <- tempfile(fileext = ".fasta")
      writeLines(fetched, scratch, useBytes = TRUE)
      sequences <- parse_fasta(scratch)
      unlink(scratch)
      accessions <- accession_from_header(names(sequences))
      keep <- nzchar(accessions) & !accessions %in% c(retained_accessions, added_accessions)
      if (any(keep)) {
        retained_sequences <- c(retained_sequences, sequences[keep])
        added_accessions <- c(added_accessions, accessions[keep])
        added_uids <- c(
          added_uids,
          if (length(sequences) == length(batch)) batch[keep] else rep(NA_character_, sum(keep))
        )
      }
    }
    write_fasta <- function(sequences, path) {
      lines <- unlist(lapply(seq_along(sequences), function(j) {
        starts <- seq(1L, nchar(sequences[[j]]), by = 80L)
        c(paste0(">", names(sequences)[j]), substring(sequences[[j]], starts, pmin(starts + 79L, nchar(sequences[[j]]))))
      }), use.names = FALSE)
      writeLines(lines, path, useBytes = TRUE)
    }
    write_fasta(retained_sequences, raw_path)

    if (length(added_accessions)) {
      ledger <- rbind(ledger, data.frame(
        marker_id = "COI", target_group = group,
        accession = added_accessions, entrez_uid = added_uids,
        first_release_id = release_id,
        added_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
        stringsAsFactors = FALSE
      ))
    }
  }

  retained_after <- previous_n + length(added_accessions)
  target_met <- retained_after >= target_n
  source_exhausted <- as.numeric(search$count) <= retained_after
  qa[[i]] <- data.frame(
    marker_id = "COI", target_group = group, available_in_ncbi = search$count,
    previous_retained = previous_n, target_retained = target_n,
    requested_new = requested_new, added_new = length(added_accessions),
    target_met = target_met || source_exhausted,
    target_shortfall_reason = if (target_met) "" else if (source_exhausted) "eligible NCBI query exhausted" else "new eligible accessions were not retained",
    dry_run = dry_run, release_id = release_id, stringsAsFactors = FALSE
  )
}

ledger <- ledger[!duplicated(paste(ledger$target_group, ledger$accession)), , drop = FALSE]
write.csv(ledger, ledger_path, row.names = FALSE, na = "")
qa <- do.call(rbind, qa)
write.csv(qa, "data/provenance/coi_reference_growth_qa.csv", row.names = FALSE, na = "")
print(qa, row.names = FALSE)
if (!dry_run && any(!qa$target_met)) stop("One or more order panels did not meet their append-only target and did not exhaust the eligible query; see growth QA.")
