#!/usr/bin/env Rscript
source("R/functions.R")
source("R/data_layer.R")
source("R/unite_reference.R")
source("R/unite_candidate_pairs.R")

# This combination was withheld from the 21-template pair catalog. Its two
# component oligos were nevertheless exhaustively searched in the broad panel.
data <- unite_reference_read(version_id = "u1-afd349e6cdcd")
stopifnot(!is.null(data), nrow(data$metadata) == 15402L, nrow(data$catalog) == 112L)
candidate <- data.frame(pair_id = "FUNGAL_COMB_TEST_ITS1_ITS4B", pair_label = "ITS1 + ITS4B",
  forward_primer = "ITS1", reverse_primer = "ITS4B",
  forward_sequence = "TCCGTAGGTGAACCTGCGG", reverse_sequence = "CAGGAGACTTGTACACGGTCCAG",
  evidence_status = "documented_pair", parent_version_id = "v3-c6fa87f256e1",
  n_panel = 21L, n_scorable = 0L, amplicon_min_bp = 42, amplicon_median_bp = 42, amplicon_max_bp = 42,
  stringsAsFactors = FALSE)
stopifnot(!any(data$pairs$forward_primer == "ITS1" & data$pairs$reverse_primer == "ITS4B"))
before <- serialize(data, NULL)
result <- unite_reference_with_candidates(data, candidate)
stopifnot(nrow(result$available_pairs) == 1L, nrow(result$unavailable_pairs) == 0L,
  result$available_pairs$candidate_origin == "on_demand_component_combination",
  result$available_pairs$evidence_status == "candidate_published_components",
  nrow(result$data$pairs) == nrow(data$pairs) + 1L,
  identical(result$data$manifest, data$manifest), identical(result$data$search_manifest, data$search_manifest),
  identical(result$data$sites, data$sites), identical(result$data$site_keys, data$site_keys),
  identical(result$data$catalog, data$catalog), identical(result$data$metadata, data$metadata),
  identical(result$data$sequences, data$sequences), identical(serialize(data, NULL), before),
  identical(result$data$pairs[seq_len(nrow(data$pairs)), names(data$pairs)], data$pairs),
  # Data-frame slicing can re-encode automatic row names; compare every actual
  # frozen column's serialized bytes, as well as the complete row values/types.
  all(vapply(names(data$pairs), function(column) identical(
    serialize(result$data$pairs[[column]][seq_len(nrow(data$pairs))], NULL),
    serialize(data$pairs[[column]], NULL)), logical(1))),
  result$data$on_demand_candidate_pair_origins$candidate_origin == "on_demand_component_combination")
added <- result$data$pairs[result$data$pairs$pair_id == candidate$pair_id, ]
stopifnot(added$n_panel == 15402L, is.na(added$n_scorable),
  is.na(added$amplicon_min_bp), is.na(added$amplicon_median_bp), is.na(added$amplicon_max_bp),
  added$evidence_status != "documented_pair")
fit <- unite_reference_pair_fit(result$data, candidate$pair_id)
states <- c("n_exact_products", "n_single_sites", "n_no_exact_sites", "n_unresolved", "n_unavailable", "n_withheld")
stopifnot(nrow(fit$details) == 15402L, sum(fit$summary$n_selected) == 15402L,
  all(rowSums(fit$summary[states]) == fit$summary$n_selected),
  identical(sort(fit$details$record_id), sort(data$metadata$record_id)),
  all(fit$details$version_id == data$manifest$id),
  sum(fit$details$fit_state == "exact_product") == 12L,
  sum(fit$details$fit_state == "partial_sites") == 3647L,
  sum(fit$details$fit_state == "no_exact_sites") == 11740L,
  sum(fit$details$fit_state == "unresolved_sites") == 3L,
  all(fit$details$amplicon_bp[fit$details$fit_state == "exact_product"] ==
    fit$details$reverse_end[fit$details$fit_state == "exact_product"] -
    fit$details$forward_start[fit$details$fit_state == "exact_product"] + 1L),
  all(is.na(fit$details$amplicon_bp[fit$details$fit_state != "exact_product"])))
again <- unite_reference_with_candidates(result$data, candidate)
stopifnot(identical(again$data$pairs, result$data$pairs),
  again$available_pairs$candidate_origin == "on_demand_component_combination")

# Existing pair rows retain their full frozen definition; a caller cannot change
# one simply by reusing its identifier, or resolve a primer by name alone.
existing <- data$pairs[1L, ]
same <- unite_reference_with_candidates(data, existing)
stopifnot(identical(same$data$pairs, data$pairs), same$available_pairs$candidate_origin == "frozen_pair")
conflict <- candidate
conflict$pair_id <- existing$pair_id
invalid <- candidate
invalid$pair_id <- "INVALID_SEQUENCE"
invalid$reverse_sequence <- paste0("A", substring(invalid$reverse_sequence, 2L))
wrong_id <- candidate
wrong_id$pair_id <- "INVALID_COMPONENT_ID"
wrong_id$forward_catalog_id <- "A_DIFFERENT_FROZEN_OLIGO"
bad_sequence <- unite_reference_with_candidates(data, invalid)
bad_id <- unite_reference_with_candidates(data, wrong_id)
bad_pair_id <- unite_reference_with_candidates(data, conflict)
stopifnot(nrow(bad_sequence$available_pairs) == 0L, nrow(bad_sequence$unavailable_pairs) == 1L,
  grepl("exact name, sequence and direction", bad_sequence$unavailable_pairs$reason, fixed = TRUE),
  identical(bad_sequence$data$pairs, data$pairs),
  nrow(bad_id$available_pairs) == 0L, grepl("catalog identifier conflicts", bad_id$unavailable_pairs$reason, fixed = TRUE),
  nrow(bad_pair_id$available_pairs) == 0L,
  grepl("pair identifier conflicts", bad_pair_id$unavailable_pairs$reason, fixed = TRUE),
  identical(bad_pair_id$data$pairs, data$pairs))
duplicates <- rbind(candidate, invalid)
duplicates$pair_id <- "DUPLICATE_CONFLICT"
stopifnot(inherits(try(unite_reference_with_candidates(data, duplicates), silent = TRUE), "try-error"))
not_loaded <- unite_reference_with_candidates(NULL, candidate)
stopifnot(is.null(not_loaded$data), nrow(not_loaded$available_pairs) == 0L,
  nrow(not_loaded$unavailable_pairs) == 1L, grepl("unavailable", not_loaded$unavailable_pairs$reason, fixed = TRUE))
stopifnot(identical(unite_reference_with_candidates(data, data.frame())$data, data))

# Two independently completed synthetic histories contain identically named
# primers with different reverse sequences. No current-version component/cache
# may silently repair an unavailable historical request.
fixture <- function(version, reverse_sequence) {
  catalog <- data.frame(catalog_id = c("F|forward|ACGTAC", paste0("R|reverse|", reverse_sequence)),
    primer_name = c("F", "R"), sequence = c("ACGTAC", reverse_sequence),
    direction = c("forward", "reverse"), primary_reference_status = "publication_cited", stringsAsFactors = FALSE)
  meta <- data.frame(record_id = c("fixture1", "fixture2"), accession = c("fixture1.1", "fixture2.1"),
    organism = "Synthetic fixture", class_group_id = "fixture-class", stringsAsFactors = FALSE)
  sites <- data.frame(record_id = "fixture1", catalog_id = catalog$catalog_id,
    start = c(1,103), end = c(6,108), status = "exact_compatible_site", stringsAsFactors = FALSE)
  pairs <- data.frame(pair_id = character(), pair_label = character(), forward_primer = character(),
    reverse_primer = character(), forward_sequence = character(), reverse_sequence = character(),
    evidence_status = character(), n_panel = integer(), n_scorable = integer(),
    amplicon_min_bp = numeric(), amplicon_median_bp = numeric(), amplicon_max_bp = numeric())
  hashes <- vapply(c("sequences", "catalog", "pairs", "sites"), function(x) digest::digest(paste(version,x), algo="sha256", serialize=FALSE), character(1))
  artifacts <- lapply(hashes, function(hash) list(sha256=hash))
  coverage <- data.frame(rank="class", taxon="Synthetic class", parent_path="Synthetic", group_id="fixture-class",
    n_source=2L, n_eligible_unique=2L, n_selected=2L, shortfall_to_100=98L, taxon_resolved=TRUE, status="sparse")
  list(manifest=list(id=version,n_selected=2L,n_source=2L,artifacts=artifacts),metadata=meta,
    catalog=catalog,pairs=pairs,sites=sites,site_keys=paste(sites$record_id,sites$catalog_id,sep="\n"),coverage=coverage,
    sequences=c(fixture1=paste0("ACGTAC",strrep("A",96),reverse_complement(reverse_sequence)),fixture2=strrep("A",108)),
    search_manifest=list(complete=TRUE,scan_rule_version="exact_iupac_unique_sites_v1",
      sequence_sha256=hashes[["sequences"]],catalog_artifact_sha256=hashes[["catalog"]],
      pair_artifact_sha256=hashes[["pairs"]],site_artifact_sha256=hashes[["sites"]],
      record_ids_sha256=digest::digest(paste(sort(meta$record_id),collapse="\n"),algo="sha256",serialize=FALSE),
      n_records=2L,n_catalog=2L,n_pairs=0L,n_site_searches=4L,n_sparse_rows=2L,n_no_exact_hits=2L))
}
v1 <- fixture("synthetic-history-1", "TTGGCC")
v2 <- fixture("synthetic-history-2", "TTGGCT")
request <- data.frame(pair_id="FUNGAL_COMB_HISTORY_TEST",pair_label="F + R",forward_primer="F",reverse_primer="R",
  forward_sequence="ACGTAC",reverse_sequence="TTGGCC",stringsAsFactors=FALSE)
accepted_v1 <- unite_reference_with_candidates(v1,request)
rejected_v2 <- unite_reference_with_candidates(v2,request)
support_v1 <- unite_reference_pair_fit(accepted_v1$data, request$pair_id)
stopifnot(nrow(accepted_v1$available_pairs)==1L,nrow(rejected_v2$available_pairs)==0L,
  nrow(rejected_v2$unavailable_pairs)==1L,identical(rejected_v2$data$pairs,v2$pairs),
  sum(support_v1$summary$n_exact_products)==1L,sum(support_v1$summary$n_selected)==2L,
  all(support_v1$details$version_id=="synthetic-history-1"),
  identical(accepted_v1$data$search_manifest,v1$search_manifest),
  nrow(unite_reference_with_candidates(v1,request)$available_pairs)==1L)
incomplete <- v1
incomplete$search_manifest$complete <- FALSE
wrong_denominator <- v1
wrong_denominator$search_manifest$n_site_searches <- 3L
missing_checksum <- v1
missing_checksum$search_manifest$site_artifact_sha256 <- NULL
stopifnot(inherits(try(unite_reference_with_candidates(incomplete,request),silent=TRUE),"try-error"),
  inherits(try(unite_reference_with_candidates(wrong_denominator,request),silent=TRUE),"try-error"),
  inherits(try(unite_reference_with_candidates(missing_checksum,request),silent=TRUE),"try-error"))

# Exercise the actual Shiny evidence module across two frozen catalog versions.
# The newer version has an available default pair, but the requested historical
# oligo identity is unavailable; it must not silently display that default.
suppressPackageStartupMessages(library(shiny))
source("R/unite_reference_ui.R")
new_request <- request
new_request$pair_id <- "AVAILABLE_ONLY_IN_HISTORY_2"
new_request$reverse_sequence <- "TTGGCT"
new_working <- unite_reference_with_candidates(v2,new_request)$data
module_versions <- data.frame(id=c(v1$manifest$id,v2$manifest$id),label=c("Synthetic history 1","Synthetic history 2"),
  manifest=c("history-1.json","history-2.json"),sha256=rep(strrep("a",64),2L),n_selected=2L,stringsAsFactors=FALSE)
module_data <- function(value) {
  value$versions <- module_versions
  value$alignment_index <- data.frame(rank=character(),group_id=character(),region=character(),file=character())
  value$regions <- data.frame(record_id=character(),region=character(),boundary_status=character(),start=integer(),end=integer())
  value
}
module_v1 <- module_data(accepted_v1$data)
module_v2 <- module_data(new_working)
read_module_fixture <- function(version_id=NULL) {
  if(is.null(version_id) || identical(version_id,module_v1$manifest$id)) module_v1 else
    if(identical(version_id,module_v2$manifest$id)) module_v2 else stop("Unknown synthetic version")
}
provided_pair <- shiny::reactiveVal(request$pair_id)
shiny::testServer(uniteReferenceServer,args=list(id="candidate-history",dataset=read_module_fixture,
    pair_provider=provided_pair,initial_mode="primers"),{
  session$setInputs(mode="primers",sites_rank="class",version="synthetic-history-1")
  session$flushReact()
  stopifnot(pair()==request$pair_id,nrow(support()$details)==2L,
    sum(support()$summary$n_selected)==2L,sum(support()$summary$n_exact_products)==1L)
  session$setInputs(version="synthetic-history-2")
  session$flushReact()
  stopifnot(data()$manifest$id=="synthetic-history-2",is.null(pair()),
    nrow(support()$details)==0L,nrow(support()$summary)==0L,
    grepl("unavailable in this selected UNITE version",output$site_tally$html,fixed=TRUE),
    !grepl("0 of 2",output$site_tally$html,fixed=TRUE))
  provided_pair("NOT_A_CATALOG_PAIR")
  session$flushReact()
  stopifnot(is.null(pair()),nrow(support()$details)==0L)
  provided_pair(NULL)
  session$flushReact()
  stopifnot(is.null(pair()),nrow(support()$details)==0L,
    grepl("Choose a forward and a reverse primer",output$site_tally$html,fixed=TRUE))
  # A deliberate selection of an available current-version pair releases the
  # guard; its full denominator and exact source version then follow that pair.
  session$setInputs(pair=new_request$pair_id)
  session$flushReact()
  stopifnot(pair()==new_request$pair_id,nrow(support()$details)==2L,
    sum(support()$summary$n_exact_products)==1L,
    all(support()$details$version_id=="synthetic-history-2"))
  provided_pair(request$pair_id)
  session$flushReact()
  stopifnot(is.null(pair()),nrow(support()$details)==0L)
  session$setInputs(version="synthetic-history-1")
  session$flushReact()
  stopifnot(pair()==request$pair_id,nrow(support()$details)==2L,
    all(support()$details$version_id=="synthetic-history-1"))
})
cat("UNITE on-demand exact-component combinations, real 15,402-reference evidence, frozen immutability, historical ownership and completed-search guards passed.\n")
