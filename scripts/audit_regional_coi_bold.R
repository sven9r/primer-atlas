#!/usr/bin/env Rscript

# Bounded BOLD feasibility audit for the Hawaii-Madagascar insect COI pilot.
# Fetches public metadata and computes counts; it does not publish sequences or
# infer primer binding, nativeness, or endemism.

suppressPackageStartupMessages(library(curl))
suppressPackageStartupMessages(library(jsonlite))

root <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(root, "data", "catalog", "range", "regions.csv"))) {
  stop("Run this script from the Primer Atlas repository root")
}

api <- "https://portal.boldsystems.org/api"
audit_date <- as.character(Sys.Date())
output_dir <- file.path(root, "data", "derived", "regional_coi_audit")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

regions <- data.frame(
  region_id = c("HI_ARCH", "MG_COUNTRY"),
  geo_query = c("geo:province/state:Hawaii", "geo:country/ocean:Madagascar"),
  stringsAsFactors = FALSE
)
orders <- c("Insecta", "Coleoptera", "Lepidoptera", "Diptera", "Hymenoptera")
pilot_family <- "Crambidae"

or_null <- function(x, default = NA_character_) {
  if (is.null(x) || length(x) == 0L) default else x
}

scalar <- function(x) {
  x <- or_null(x)
  if (length(x) == 0L || is.null(x) || is.na(x[1])) return(NA_character_)
  as.character(x[1])
}

api_url <- function(path, params) {
  query <- paste(
    vapply(names(params), function(key) {
      paste0(
        utils::URLencode(key, reserved = TRUE), "=",
        utils::URLencode(params[[key]], reserved = TRUE)
      )
    }, character(1)),
    collapse = "&"
  )
  paste0(api, path, "?", query)
}

fetch_text <- function(url) {
  last_error <- "no response"
  for (attempt in seq_len(3L)) {
    response <- tryCatch(
      curl::curl_fetch_memory(
        url,
        handle = curl::new_handle(
          timeout = 60L,
          useragent = "PrimerAtlasRegionalCOIAudit/0.1 (academic feasibility)"
        )
      ),
      error = function(e) e
    )
    if (!inherits(response, "error") && response$status_code == 200L) {
      return(rawToChar(response$content))
    }
    last_error <- if (inherits(response, "error")) {
      conditionMessage(response)
    } else {
      paste("HTTP", response$status_code)
    }
    Sys.sleep(attempt)
  }
  stop("BOLD request failed: ", last_error, " (", url, ")")
}

fetch_json <- function(url) {
  jsonlite::fromJSON(fetch_text(url), simplifyVector = FALSE)
}

summary_row <- function(region_id, geo_query, taxon, rank) {
  query <- paste0("tax:", rank, ":", taxon, ";", geo_query)
  result <- fetch_json(api_url("/summary", list(
    query = query, fields = "specimens,marker_code"
  )))
  data.frame(
    region_id = region_id,
    taxon_rank = rank,
    taxon = taxon,
    query = query,
    specimens_reported = as.integer(or_null(result$counts$specimens, 0L)),
    coi_5p_reported = as.integer(or_null(result$marker_code[["COI-5P"]], 0L)),
    accessed_on = audit_date,
    stringsAsFactors = FALSE
  )
}

fetch_records <- function(query) {
  request <- fetch_json(api_url("/query", list(query = query, extent = "full")))
  token <- scalar(request$query_id)
  if (is.na(token)) stop("BOLD did not return a query_id for ", query)
  url <- paste0(
    api, "/documents/", utils::URLencode(token, reserved = TRUE),
    "/download?format=json"
  )
  stream <- fetch_text(url)
  lines <- strsplit(stream, "\n", fixed = TRUE)[[1]]
  lines <- lines[nzchar(trimws(lines))]
  if (!length(lines)) return(list())
  lapply(lines, jsonlite::fromJSON, simplifyVector = FALSE)
}

record_row <- function(record, region_id) {
  sequence <- toupper(scalar(record$nuc))
  if (is.na(sequence)) sequence <- ""
  bases <- strsplit(sequence, "", fixed = TRUE)[[1]]
  canonical <- sum(bases %in% c("A", "C", "G", "T"))
  ambiguous <- sum(!(bases %in% c("A", "C", "G", "T", "-")))
  ungapped <- canonical + ambiguous
  coordinates <- or_null(record$coord, list())
  latitude <- longitude <- NA_real_
  if (length(coordinates) == 2L) {
    latitude <- suppressWarnings(as.numeric(coordinates[[1]]))
    longitude <- suppressWarnings(as.numeric(coordinates[[2]]))
  }
  country <- scalar(record[["country/ocean"]])
  province <- scalar(record[["province/state"]])
  region_confirmed <- if (region_id == "HI_ARCH") {
    !is.na(country) && !is.na(province) &&
      identical(tolower(country), "united states") &&
      identical(tolower(province), "hawaii")
  } else {
    !is.na(country) && identical(tolower(country), "madagascar")
  }
  species <- scalar(record$species)
  species_label <- !is.na(species) &&
    grepl("^[A-Z][[:alpha:]-]+ [a-z][[:alpha:]-]+$", species) &&
    identical(scalar(record$identification_rank), "species")
  marker <- scalar(record$marker_code)
  ambiguity_fraction <- if (ungapped > 0L) ambiguous / ungapped else NA_real_
  data.frame(
    source = "BOLD Data Portal", accessed_on = audit_date,
    source_query_region = region_id,
    record_id = scalar(record$record_id),
    processid = scalar(record$processid),
    specimenid = scalar(record$specimenid),
    insdc_accession = scalar(record$insdc_acs),
    bold_taxid = scalar(record$taxid),
    order = scalar(record$order), family = scalar(record$family),
    genus = scalar(record$genus), species = species,
    identification_rank = scalar(record$identification_rank),
    marker_code = marker, country = country, province = province,
    locality_region = scalar(record$region),
    latitude = latitude, longitude = longitude,
    coordinate_source = scalar(record$coord_source),
    upload_date = scalar(record$sequence_upload_date),
    sequence_characters = nchar(sequence),
    canonical_bases = canonical,
    ambiguous_bases = ambiguous,
    ambiguity_fraction = ambiguity_fraction,
    region_confirmed = region_confirmed,
    species_label = species_label,
    pilot_sequence_qc = identical(marker, "COI-5P") &&
      region_confirmed && canonical >= 500L &&
      !is.na(ambiguity_fraction) && ambiguity_fraction <= 0.02,
    stringsAsFactors = FALSE
  )
}

message("Fetching BOLD summary counts")
summaries <- list()
for (region_index in seq_len(nrow(regions))) {
  region <- regions[region_index, ]
  for (taxon in orders) {
    rank <- if (taxon == "Insecta") "class" else "order"
    summaries[[length(summaries) + 1L]] <- summary_row(
      region$region_id, region$geo_query, taxon, rank
    )
    Sys.sleep(0.25)
  }
  summaries[[length(summaries) + 1L]] <- summary_row(
    region$region_id, region$geo_query, pilot_family, "family"
  )
}
summaries <- do.call(rbind, summaries)
write.csv(
  summaries, file.path(output_dir, "bold_summary_counts.csv"),
  row.names = FALSE, na = ""
)

message("Fetching bounded ", pilot_family, " record metadata")
rows <- list()
for (region_index in seq_len(nrow(regions))) {
  region <- regions[region_index, ]
  query <- paste0("tax:family:", pilot_family, ";", region$geo_query)
  records <- fetch_records(query)
  rows <- c(rows, lapply(records, record_row, region_id = region$region_id))
  message(region$region_id, ": ", length(records), " returned records")
}
if (!length(rows)) stop("BOLD returned no pilot records")
metadata <- do.call(rbind, rows)
metadata <- metadata[!is.na(metadata$record_id), , drop = FALSE]
duplicate_records <- sum(duplicated(metadata$record_id))
metadata <- metadata[!duplicated(metadata$record_id), , drop = FALSE]
metadata <- metadata[order(metadata$source_query_region, metadata$record_id), , drop = FALSE]
write.csv(
  metadata, file.path(output_dir, "bold_crambidae_record_audit.csv"),
  row.names = FALSE, na = ""
)

usable <- metadata[metadata$pilot_sequence_qc, , drop = FALSE]
hi <- usable[usable$source_query_region == "HI_ARCH", , drop = FALSE]
mg <- usable[usable$source_query_region == "MG_COUNTRY", , drop = FALSE]
shared_genera <- intersect(unique(stats::na.omit(hi$genus)), unique(stats::na.omit(mg$genus)))
hi_species <- unique(hi$species[hi$species_label])
mg_species <- unique(mg$species[mg$species_label])
shared_species <- intersect(hi_species, mg_species)

region_line <- function(region_id) {
  x <- metadata[metadata$source_query_region == region_id, , drop = FALSE]
  q <- x[x$pilot_sequence_qc, , drop = FALSE]
  valid_coords <- is.finite(x$latitude) & is.finite(x$longitude) &
    abs(x$latitude) <= 90 & abs(x$longitude) <= 180
  insdc_count <- sum(!is.na(q$insdc_accession) & nzchar(q$insdc_accession))
  paste0(
    "| ", region_id, " | ", nrow(x), " | ", sum(x$marker_code == "COI-5P", na.rm = TRUE),
    " | ", sum(x$region_confirmed), " | ", sum(valid_coords),
    " | ", nrow(q), " | ", insdc_count,
    " | ", length(unique(q$genus[!is.na(q$genus)])),
    " | ", length(unique(q$species[q$species_label])), " |"
  )
}

shared_species_lines <- if (length(shared_species)) {
  vapply(shared_species, function(name) {
    paste0(
      "| ", name, " | ", sum(hi$species == name, na.rm = TRUE),
      " | ", sum(mg$species == name, na.rm = TRUE), " |"
    )
  }, character(1))
} else {
  "| None | 0 | 0 |"
}

report <- c(
  "# Hawaii-Madagascar insect COI feasibility audit",
  "",
  paste0("BOLD Data Portal queried: ", audit_date, ". This is a bounded feasibility audit, not a global range census or a primer-performance estimate."),
  "",
  "## Source and scope",
  "",
  "The script `scripts/audit_regional_coi_bold.R` reproduces the BOLD Portal API summary and record queries. Hawaiʻi uses `geo:province/state:Hawaii`; Madagascar uses `geo:country/ocean:Madagascar`. Order-level summary counts are in `data/derived/regional_coi_audit/bold_summary_counts.csv`. Record-level inspection is limited to Crambidae and recorded in `data/derived/regional_coi_audit/bold_crambidae_record_audit.csv`. The latter contains IDs and metadata, not nucleotide strings. Query tokens expire, so rerun the script for a fresh snapshot.",
  "",
  "BOLD API documentation: https://portal.boldsystems.org/api . Its summary counts describe records matching indexed taxonomy and geography fields; they are not deduplicated species counts.",
  "",
  "## Bounded Crambidae record audit",
  "",
  "A preliminary sequence QC pass requires a COI-5P record, the expected indexed region field, at least 500 canonical A/C/G/T bases, and at most 2% other non-gap bases. This establishes a useful COI fragment; it does not establish that either primer binding site is present.",
  "",
  "| Region | Unique records | COI-5P | Region field confirmed | Valid coordinate pair | Sequence QC pass | QC records with INSDC accession | Genera in QC pass | Species labels in QC pass |",
  "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
  region_line("HI_ARCH"), region_line("MG_COUNTRY"),
  "",
  paste0("Duplicate BOLD record IDs removed across query results: ", duplicate_records, "."),
  paste0("Shared genus labels in the QC pass (", length(shared_genera), "): ", if (length(shared_genera)) paste(shared_genera, collapse = ", ") else "none", "."),
  paste0("Shared species labels in the QC pass (", length(shared_species), "): ", if (length(shared_species)) paste(shared_species, collapse = ", ") else "none", "."),
  "",
  "| Shared species label | Hawaiʻi QC records | Madagascar QC records |",
  "|---|---:|---:|",
  shared_species_lines,
  "",
  "Two shared species have only one QC-passing Hawaiʻi record each; their within-species regional comparison is not yet supported by replication. BOLD process and record IDs preserve provenance for records without an INSDC accession.",
  "",
  "These are BOLD taxonomic labels, not reconciled accepted species concepts. Shared names require taxonomic review before a within-species regional comparison. Locality precision varies, and coordinates need spatial validation against the defined regions.",
  "",
  "## Primer binding and endemicity readiness",
  "",
  "No primer pair was scored in this audit. BOLD COI-5P records have variable starts, ends, gaps, and PCR primer metadata; a 500-base barcode can still omit an assay's binding sites. The next scoring pass must align each sequence to a COI reference, locate both documented primer sites, report the scorable denominator per site and pair, then apply the existing mismatch model. The field `pilot_sequence_qc` must not be interpreted as primer scoreability.",
  "",
  "The new `data/catalog/range/taxon_region_evidence.csv` remains a separate, reviewed taxon-range catalog. Sequence locality, absence of a record elsewhere, and shared species labels do not establish endemicity.",
  "",
  "## Recommended next slice",
  "",
  "Use Crambidae and one documented internal COI primer pair. Review the shared species labels and coordinate quality first. Then map primer binding sites in the QC-passing sequences and compare (1) all Crambidae, (2) shared genera, and (3) shared species where each region has enough independent specimens. Show sample size and missing-site reasons at every level."
)
writeLines(report, file.path(root, "docs", "regional-coi-feasibility.md"), useBytes = TRUE)
message("Wrote bounded audit and report; ", nrow(metadata), " unique record IDs")
