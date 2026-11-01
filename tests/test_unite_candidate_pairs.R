#!/usr/bin/env Rscript
source("R/functions.R")
source("R/data_layer.R")
source("R/unite_reference.R")
source("R/fungal_primer_orientation.R")
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
    direction = c("forward", "reverse"), gene_locus="ITS", primary_reference_status = "publication_cited", stringsAsFactors = FALSE)
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

# A computational reverse-role ITS3 oligo has the parent's exact genomic query.
# Reuse all sparse parent states on this actual 15,402-record cohort; missing
# sparse rows remain no-hits under the ORIGINAL112 completion proof.
source("R/fungal_explorer.R")
source("R/fungal_versions.R")
ncbi <- fungal_read_snapshot("v3-c6fa87f256e1",remote=FALSE)
supplements <- fungal_orientation_catalog(ncbi)
rc_its3 <- supplements[supplements$parent_primer_name=="ITS3",,drop=FALSE]
rc_request <- candidate
rc_request$pair_id <- "FUNGAL_COMB_TEST_ITS1_ITS3_RC"
rc_request$pair_label <- "ITS1 + ITS3 [RC]"
rc_request$reverse_primer <- rc_its3$primer_name
rc_request$reverse_sequence <- rc_its3$sequence
rc_request$reverse_catalog_id <- rc_its3$catalog_id
saved_search <- .fungal_exact_site
.fungal_exact_site <- function(...) stop("UNITE derived components must never rescan templates")
rc_result <- unite_reference_with_candidates(data,rc_request,catalog=supplements)
rc_fit <- unite_reference_pair_fit(rc_result$data,rc_request$pair_id)
.fungal_exact_site <- saved_search
stopifnot(nrow(rc_result$available_pairs)==1L,nrow(rc_result$unavailable_pairs)==0L,
  rc_result$available_pairs$candidate_origin=="on_demand_reverse_complement_combination",
  rc_result$available_pairs$has_derived_reverse_complement,
  rc_result$available_pairs$evidence_status=="candidate_incomplete_provenance",
  rc_result$available_pairs$reverse_parent_catalog_id==rc_its3$parent_catalog_id,
  nrow(rc_result$data$catalog)==113L,length(rc_result$data$on_demand_derived_component_ids)==1L,
  identical(rc_result$data$frozen_component_catalog,data$catalog),
  identical(rc_result$data$frozen_component_sites,data$sites),
  identical(serialize(rc_result$data$manifest,NULL),serialize(data$manifest,NULL)),
  identical(serialize(rc_result$data$search_manifest,NULL),serialize(data$search_manifest,NULL)),
  rc_result$data$search_manifest$n_catalog==112L,
  rc_result$data$search_manifest$n_sparse_rows==nrow(data$sites),
  rc_result$data$catalog$parent_source_version[rc_result$data$catalog$catalog_id==rc_its3$catalog_id]==ncbi$manifest$id,
  all(rc_fit$details$version_id==data$manifest$id),nrow(rc_fit$details)==15402L,
  sum(rc_fit$summary$n_selected)==15402L,all(rowSums(rc_fit$summary[states])==rc_fit$summary$n_selected),
  sum(rc_fit$summary$n_exact_products)>0L,identical(serialize(data,NULL),before))
for (name in c("catalog","sites","pairs")) {
  for (column in names(data[[name]])) stopifnot(identical(
    serialize(rc_result$data[[name]][[column]][seq_len(nrow(data[[name]]))],NULL),serialize(data[[name]][[column]],NULL)))
}
parent <- data$sites[data$sites$catalog_id==rc_its3$parent_catalog_id,,drop=FALSE]
mirrored <- rc_result$data$sites[rc_result$data$sites$catalog_id==rc_its3$catalog_id,,drop=FALSE]
stopifnot(nrow(mirrored)==nrow(parent),nrow(rc_result$data$sites)==nrow(data$sites)+nrow(parent))
for (field in c("record_id","start","end","status","n_matches")) if (field %in% names(parent)) stopifnot(identical(mirrored[[field]],parent[[field]]))
parent_positions <- match(paste(data$metadata$record_id,rc_its3$parent_catalog_id,sep="\n"),data$site_keys)
expected_reverse_status <- data$sites$status[parent_positions]
expected_reverse_status[is.na(parent_positions)] <- "no_exact_binding_site"
stopifnot(identical(rc_fit$details$reverse_status,expected_reverse_status),
  all(is.na(rc_fit$details$amplicon_bp[rc_fit$details$fit_state!="exact_product"])))
again_rc <- unite_reference_with_candidates(rc_result$data,rc_request,catalog=supplements)
stopifnot(identical(again_rc$data$catalog,rc_result$data$catalog),identical(again_rc$data$sites,rc_result$data$sites),
  identical(again_rc$data$pairs,rc_result$data$pairs),
  again_rc$available_pairs$candidate_origin=="on_demand_reverse_complement_combination")

# Caller-supplied provenance and changed-role sequences are verified, never
# trusted merely because an ID uses the RC prefix or names the parent.
for (field in c("parent_catalog_id","parent_sequence","parent_direction","parent_reference","sequence","source_version","derivation_rule")) {
  forged <- rc_its3;forged[[field]] <- paste0(forged[[field]],"_forged")
  withheld <- unite_reference_with_candidates(data,rc_request,catalog=forged)
  stopifnot(nrow(withheld$available_pairs)==0L,nrow(withheld$unavailable_pairs)==1L,
    identical(withheld$data$catalog,data$catalog),identical(withheld$data$sites,data$sites),
    identical(withheld$data$pairs,data$pairs))
}
same_oligo <- rc_its3;same_oligo$sequence <- same_oligo$parent_sequence
stopifnot(nrow(unite_reference_with_candidates(data,rc_request,catalog=same_oligo)$available_pairs)==0L)
bad_proof <- data;bad_proof$search_manifest$complete <- FALSE
stopifnot(inherits(try(unite_reference_with_candidates(bad_proof,rc_request,catalog=supplements),silent=TRUE),"try-error"))
bad_site <- data
hit <- which(bad_site$sites$catalog_id==rc_its3$parent_catalog_id & bad_site$sites$status=="exact_compatible_site")[1L]
bad_site$sites$start[hit] <- bad_site$sites$start[hit]+1L
stopifnot(inherits(try(unite_reference_with_candidates(bad_site,rc_request,catalog=supplements),silent=TRUE),"try-error"))
truncated_working <- rc_result$data
truncated_working$sites <- truncated_working$sites[-which(truncated_working$sites$catalog_id==rc_its3$catalog_id)[1L],]
truncated_working$site_keys <- paste(truncated_working$sites$record_id,truncated_working$sites$catalog_id,sep="\n")
stopifnot(inherits(try(unite_reference_with_candidates(truncated_working,rc_request),silent=TRUE),"try-error"))
changed_working_bases <- rc_result$data
changed_working_bases$sequences[[parent$record_id[which(parent$status=="exact_compatible_site")[1L]]]] <- strrep("N",nchar(changed_working_bases$sequences[[parent$record_id[which(parent$status=="exact_compatible_site")[1L]]]]))
stopifnot(inherits(try(unite_reference_with_candidates(changed_working_bases,rc_request),silent=TRUE),"try-error"))

# The same RC source variant cannot borrow parent sites from a later history
# where the original name now denotes a different oligo. Opposite-role sites
# may also be overlapping: cache reuse establishes sites, not a valid product.
history_rc <- fungal_orientation_catalog(v1)
history_rc <- history_rc[history_rc$parent_primer_name=="F",,drop=FALSE]
history_request <- request
history_request$pair_id <- "FUNGAL_COMB_HISTORY_RC_TEST"
history_request$pair_label <- "F + F [RC]"
history_request$reverse_primer <- history_rc$primer_name
history_request$reverse_sequence <- history_rc$sequence
history_request$reverse_catalog_id <- history_rc$catalog_id
history_accepted <- unite_reference_with_candidates(v1,history_request,catalog=history_rc)
history_fit <- unite_reference_pair_fit(history_accepted$data,history_request$pair_id)
stopifnot(sum(history_fit$summary$n_withheld)==1L,sum(history_fit$summary$n_no_exact_sites)==1L,
  sum(history_fit$summary$n_exact_products)==0L,all(is.na(history_fit$details$amplicon_bp)))
different_parent <- v2
different_parent$catalog$sequence[1L] <- "ACGTAT"
different_parent$catalog$catalog_id[1L] <- "F|forward|ACGTAT"
different_parent$sites$catalog_id[different_parent$sites$catalog_id==v1$catalog$catalog_id[1L]] <- "F|forward|ACGTAT"
different_parent$site_keys <- paste(different_parent$sites$record_id,different_parent$sites$catalog_id,sep="\n")
history_rejected <- unite_reference_with_candidates(different_parent,history_request,catalog=history_rc)
stopifnot(nrow(history_rejected$available_pairs)==0L,nrow(history_rejected$unavailable_pairs)==1L,
  identical(history_rejected$data$catalog,different_parent$catalog),identical(history_rejected$data$sites,different_parent$sites))
existing_with_supplements <- unite_reference_with_candidates(data,existing,catalog=supplements)
stopifnot(identical(existing_with_supplements$data$catalog,data$catalog),identical(existing_with_supplements$data$sites,data$sites),
  identical(existing_with_supplements$data$pairs,data$pairs),existing_with_supplements$available_pairs$candidate_origin=="frozen_pair")
cat("UNITE on-demand exact-component combinations, real 15,402-reference evidence, frozen immutability, historical ownership and completed-search guards passed.\n")
