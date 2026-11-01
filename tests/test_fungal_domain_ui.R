#!/usr/bin/env Rscript
.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
suppressPackageStartupMessages({library(shiny);library(bslib);library(DT)})
source("R/functions.R")
source("R/data_layer.R")
source("R/fungal_alignment.R")
source("R/fungal_alignment_ui.R")
source("R/fungal_explorer.R")
source("R/fungal_versions.R")

# Synthetic structural-reference domains deliberately arrive out of numerical
# order. The second template is reverse-aligned; neither may borrow boundaries.
raw <- c(A=paste(rep("ACGT",1250L),collapse=""),B=paste(rep("TGCA",1250L),collapse=""))
gapped <- function(sequence) paste0(substr(sequence,1L,1000L),paste(rep("-",200L),collapse=""),substr(sequence,1001L,5000L))
index <- fungal_alignment_index(c(A=gapped(raw[["A"]]),`_R_B`=gapped(reverse_complement(raw[["B"]]))),raw)
annotations <- data.frame(accession=c(rep("A",5L),"B"),landmark_id=c("28S","D10","D2","D1","D3","D1"),
  label=c("28S","D10","D2","D1","D3","D1"),start=c(500L,3200L,1800L,1500L,NA_integer_,1500L),
  end=c(4500L,3300L,2000L,1600L,NA_integer_,1600L),
  landmark_type=c("gene",rep("divergent_domain",5L)),parent_landmark=c("",rep("28S",5L)),
  evidence_status=c("",rep("mapped_reference_domain",5L)),
  source_note=c("NCBI annotated complete forward feature",rep("Synthetic structural-reference mapping fixture",5L)),
  domain_scheme=c("",rep("Fixture scheme",5L)),source_reference=c("",rep("Fixture reference A",5L)),
  boundary_status=c("complete","full","partial","full","missing","full"),stringsAsFactors=FALSE)
mapped <- fungal_alignment_landmark_rows(index,annotations)
stopifnot(mapped$evidence_status[1L]=="annotated_reference",
  all(mapped$evidence_status[-1L]=="mapped_reference_domain"),
  mapped$alignment_status[5L]=="invalid_raw_coordinates",is.na(mapped$left_col[5L]))
available <- mapped[mapped$accession=="A" & mapped$alignment_status=="mapped",,drop=FALSE]
stopifnot(identical(.fungal_alignment_ui_domain_rows(available)$label,c("D1","D2","D10")))
domain_title <- .fungal_alignment_ui_region_title(available[available$label=="D2",,drop=FALSE])
stopifnot(grepl("Projected structural-reference domain",domain_title,fixed=TRUE),
  !grepl("Frozen NCBI annotation",domain_title,fixed=TRUE),
  grepl("D2 (partial)",domain_title,fixed=TRUE),grepl("Fixture scheme",domain_title,fixed=TRUE),
  grepl("Fixture reference A",domain_title,fixed=TRUE))

row_a <- mapped[mapped$accession=="A" & mapped$label=="D1",,drop=FALSE]
row_b <- mapped[mapped$accession=="B" & mapped$label=="D1",,drop=FALSE]
window_a <- .fungal_alignment_ui_region_window(index,row_a)
window_b <- .fungal_alignment_ui_region_window(index,row_b)
stopifnot(identical(window_a$raw,c(500L,3600L)),identical(window_a$columns,c(500L,3800L)),
  identical(window_b$raw,c(500L,3600L)),identical(window_b$columns,c(1601L,4701L)),
  identical(range(which(index$maps$B>=window_b$columns[1L] & index$maps$B<=window_b$columns[2L])),c(500L,3600L)),
  is.null(.fungal_alignment_ui_region_window(index,mapped[5L,,drop=FALSE])))
whole_parent <- .fungal_alignment_ui_region_window(index,mapped[1L,,drop=FALSE])
stopifnot(identical(whole_parent$raw,c(1L,5000L)),identical(whole_parent$columns,c(1L,5200L)))

product <- data.frame(accession="A",forward_start=1450L,reverse_end=1950L,alignment_status="mapped")
coverage <- .fungal_alignment_ui_domain_coverage(product,mapped,"A")
stopifnot(coverage$available,identical(coverage$overlaps$label,c("D1","D2")),
  identical(coverage$overlaps$coverage,c("full","partial")),identical(coverage$overlaps$overlap_bp,c(101L,151L)))
product$forward_start <- 1750L
product$reverse_end <- 2100L
partial_source <- .fungal_alignment_ui_domain_coverage(product,mapped,"A")
stopifnot(identical(partial_source$overlaps$label,"D2"),partial_source$overlaps$coverage=="partial",
  partial_source$overlaps$overlap_bp==201L)
product$forward_start <- 1L
product$reverse_end <- 100L
stopifnot(nrow(.fungal_alignment_ui_domain_coverage(product,mapped,"A")$overlaps)==0L)
product$alignment_status <- "partial_exact_sites"
stopifnot(!.fungal_alignment_ui_domain_coverage(product,mapped,"A")$available)

empty_pairs <- data.frame(pair_id=character(),pair_label=character())
empty_binding <- data.frame(pair_id=character(),accession=character())
taxonomy <- data.frame(accession=c("A","B"),organism=c("Fixture A","Fixture B"),class="Fixture class")
bp <- as.character(fungal_alignment_template_bp_view(index,empty_binding,empty_pairs,taxonomy,
  reference_accession="A",active_id=character(),domain=c(1L,5000L),ceilings=c(500,2000,3500),ns=identity,landmarks=mapped))
alignment <- as.character(fungal_alignment_overview(index,empty_binding,empty_pairs,taxonomy,"A",c("A","B"),
  character(),character(),c(1L,5200L),c(500,2000,3500),identity,landmarks=mapped))
stopifnot(grepl('class="fungal-template-bp-domain"',bp,fixed=TRUE),
  grepl('data-parent-landmark="28S"',bp,fixed=TRUE),
  grepl('data-start-base="1500"',bp,fixed=TRUE),
  !grepl('data-landmark-id="D3"',bp,fixed=TRUE),
  grepl('class="fungal-alignment-domain"',alignment,fixed=TRUE),
  grepl('data-start-col="3601"',alignment,fixed=TRUE),
  grepl('data-end-col="3701"',alignment,fixed=TRUE))
# Child rectangles are on a separate track, not painted over the 28S band.
bp_parent <- regmatches(bp,regexpr('<g class="fungal-template-bp-region"[^>]*data-landmark-id="28S"[^>]*>.*?</g>',bp))
bp_domain <- regmatches(bp,regexpr('<g class="fungal-template-bp-domain"[^>]*data-landmark-id="D2"[^>]*>.*?</g>',bp))
stopifnot(grepl('y="57"',bp_parent,fixed=TRUE),grepl('y="103"',bp_domain,fixed=TRUE),
  grepl('stroke-dasharray="3 2"',bp_domain,fixed=TRUE))
hidden <- as.character(fungal_alignment_template_bp_view(index,empty_binding,empty_pairs,taxonomy,
  "A",character(),c(1L,5000L),c(500,2000,3500),identity,landmarks=mapped,show_regions=FALSE))
stopifnot(!grepl('class="fungal-template-bp-domain"',hidden,fixed=TRUE))

# Run the real module with an isolated in-memory augmentation of frozen V2.
# No snapshot files, manifests, index entries or previous versions are changed.
versions <- fungal_read_index(remote=FALSE)
base_version <- versions$id[grepl("^v2-",versions$id)][1L]
fixture <- fungal_read_snapshot(base_version,remote=FALSE)
fixture_domains <- annotations[annotations$accession=="A" & annotations$label!="28S",,drop=FALSE]
fixture_domains$accession <- "PZ962702.1"
fixture_domains$start <- c(3500L,3000L,2600L,NA_integer_)
fixture_domains$end <- c(3600L,3150L,2700L,NA_integer_)
copy_domain <- fixture_domains[fixture_domains$label=="D2",,drop=FALSE]
copy_domain$landmark_id <- "D2_copy"
copy_domain$start <- 3300L
copy_domain$end <- 3400L
fixture_domains <- rbind(fixture_domains,copy_domain)
columns <- union(names(fixture$landmarks),names(fixture_domains))
for(name in setdiff(columns,names(fixture$landmarks))) fixture$landmarks[[name]] <- ""
fixture$landmarks <- rbind(fixture$landmarks[,columns],fixture_domains[,columns])
copy_gene <- fixture$landmarks[fixture$landmarks$accession=="PZ962702.1" & fixture$landmarks$label=="28S",,drop=FALSE]
copy_gene$landmark_id <- "28S_copy"
copy_gene$start <- 7000L
copy_gene$end <- 7600L
copy_gene$boundary_status <- "partial"
fixture$landmarks <- rbind(fixture$landmarks,copy_gene)
env <- new.env(parent=globalenv())
env$fungal_read_index <- function(...) versions
env$fungal_read_snapshot <- function(version_id,...) {
  if(identical(version_id,base_version)) fixture else fungal_read_snapshot(version_id,remote=FALSE)
}
source("R/fungal_ui.R",local=env)
messages <- new.env(parent=emptyenv())
mock <- shiny::MockShinySession$new()
mock$sendInputMessage <- function(inputId,message) messages[[inputId]] <- message
shiny::testServer(env$fungalExplorerServer,args=list(id="domains"),session=mock,{
  session$setInputs(version=base_version,platform="ALL",phylum="ALL",class="ALL",locus="ALL",search="",
    examples_only=TRUE,show_unresolved=TRUE,alignment_reference="PZ962702.1",show_regions=TRUE)
  session$flushReact()
  initial <- alignment_focus_row()
  initial_bounds <- c(max(1L,initial$forward_start-1000L),
    min(length(alignment_index()$maps[["PZ962702.1"]]),initial$reverse_end+2000L))
  initial_columns <- range(unname(alignment_index()$maps[["PZ962702.1"]][initial_bounds]))
  stopifnot(identical(as.integer(messages[["alignment_range"]]$value),as.integer(initial_columns)))
  calculation <- output$alignment_bp_accounting$html
  concise <- strsplit(calculation,'<details class="fungal-length-details"',fixed=TRUE)[[1L]][1L]
  concise_text <- gsub("[[:space:]]+"," ",gsub("<[^>]+>"," ",concise))
  stopifnot(grepl('class="fungal-selected-product"',concise,fixed=TRUE),
    grepl(paste0(format(as.integer(initial$amplicon_bp),big.mark=",",trim=TRUE)," bp"),concise_text,fixed=TRUE),
    !grepl("gap columns",concise,fixed=TRUE),!grepl("Where the bp",concise,fixed=TRUE),
    grepl('<details class="fungal-length-details">',calculation,fixed=TRUE),
    grepl("How length is calculated",calculation,fixed=TRUE),grepl("gap columns",calculation,fixed=TRUE))
  controls <- output$alignment_regions$html
  at <- function(name) regexpr(paste0('data-landmark="',name,'"'),controls,fixed=TRUE)[1L]
  stopifnot(at("D1")>0L,at("D1")<at("D2"),at("D2")<at("D10"),at("D3")==-1L,
    grepl("structural reference",controls,fixed=TRUE),
    grepl('class="fungal-annotation-evidence"',controls,fixed=TRUE),
    !grepl('class="fungal-annotation-evidence" open',controls,fixed=TRUE),
    grepl("28S · bases 2,423–5,821",controls,fixed=TRUE),
    grepl("28S (partial) · bases 7,000–7,600",controls,fixed=TRUE),at("28S_copy")>0L,
    grepl("D2 (partial) · bases 3,000–3,150",controls,fixed=TRUE),
    grepl("D2 (partial) · bases 3,300–3,400",controls,fixed=TRUE),at("D2_copy")>0L)
  stopifnot(grepl("No mapped D domain overlap",output$alignment_status$html,fixed=TRUE))
  prior_pair <- active_pair()$pair_id
  session$setInputs(region_click=list(accession="PZ962702.1",landmark="D1",nonce=1L))
  session$flushReact()
  sent <- messages[["alignment_range"]]
  if(is.null(sent)) stop("Missing domain window update; captured inputs: ",paste(ls(messages),collapse=", "))
  expected <- range(unname(alignment_index()$maps[["PZ962702.1"]][c(1600L,4700L)]))
  stopifnot(!is.null(sent),identical(as.integer(sent$value),as.integer(expected)),
    identical(active_pair()$pair_id,prior_pair),alignment_reference()=="PZ962702.1")
  session$setInputs(alignment_range=sent$value,product_coordinates="template_bases")
  session$flushReact()
  stopifnot(identical(template_base_domain(),c(1600L,4700L)),
    grepl('class="fungal-template-bp-domain"',output$alignment_overview$html,fixed=TRUE))
  session$setInputs(product_coordinates="alignment_columns")
  session$flushReact()
  stopifnot(identical(active_pair()$pair_id,prior_pair),identical(template_base_domain(),c(1600L,4700L)),
    grepl('class="fungal-alignment-domain"',output$alignment_overview$html,fixed=TRUE))
  # Candidate count, comparison scope and selected-pair changes must not
  # overwrite an explicit domain focus with the extent of all primer sites.
  messages[["alignment_range"]] <- NULL
  session$setInputs(examples_only=FALSE)
  session$flushReact()
  own_scorable <- sequence_rows()$pair_id[sequence_rows()$accession=="PZ962702.1" & sequence_rows()$status=="scorable"]
  next_pair <- head(setdiff(intersect(candidates()$pair_id,own_scorable),prior_pair),1L)
  stopifnot(length(next_pair)==1L,nrow(candidates())>20L)
  session$setInputs(pair_click=list(pair=next_pair,nonce=3L),alignment_scope="active",compare_pairs=next_pair)
  session$flushReact()
  stopifnot(is.null(messages[["alignment_range"]]),identical(active_pair()$pair_id,next_pair),
    identical(template_base_domain(),c(1600L,4700L)),alignment_focus_state()$kind=="region")
  # A different reference cannot inherit the original template's D1 bounds.
  # The fallback retains the intent so returning to that reference restores it.
  session$setInputs(alignment_reference="FJ430740.1")
  session$flushReact()
  stopifnot(nrow(focused_alignment_region())==0L,
    grepl("focused region has no verified boundaries",output$alignment_status$html,fixed=TRUE))
  fallback <- alignment_focus_row()
  if(nrow(fallback)==1L && fallback$alignment_status=="mapped") {
    interval <- data.frame(accession=fallback$accession,start=fallback$forward_start,end=fallback$reverse_end,alignment_status="mapped")
    fallback_columns <- .fungal_alignment_ui_region_window(alignment_index(),interval)$columns
  } else fallback_columns <- c(1L,alignment_index()$n_columns)
  if(!identical(as.integer(messages[["alignment_range"]]$value),as.integer(fallback_columns)))
    stop("Reference fallback differs: got ",paste(messages[["alignment_range"]]$value,collapse=","),
      "; expected ",paste(fallback_columns,collapse=","),"; reference ",alignment_reference(),
      "; product ",paste(fallback$pair_id,fallback$alignment_status,collapse=","),
      "; focus reference ",paste(focus_inputs()$reference,collapse=","),"; mode ",alignment_focus_state()$kind,
      "; captured ",paste(ls(messages),collapse=","))
  session$setInputs(alignment_reference="PZ962702.1")
  session$flushReact()
  stopifnot(identical(as.integer(messages[["alignment_range"]]$value),as.integer(expected)))
  messages[["alignment_range"]] <- NULL
  session$setInputs(region_click=list(accession="FJ430740.1",landmark="D1",nonce=2L))
  session$flushReact()
  stopifnot(is.null(messages[["alignment_range"]]))
  session$setInputs(alignment_focus=1L)
  session$flushReact()
  product <- alignment_focus_row()
  product_raw <- c(max(1L,product$forward_start-1000L),min(length(alignment_index()$maps[["PZ962702.1"]]),product$reverse_end+2000L))
  product_columns <- range(unname(alignment_index()$maps[["PZ962702.1"]][product_raw]))
  stopifnot(alignment_focus_state()$kind=="product",
    identical(as.integer(messages[["alignment_range"]]$value),as.integer(product_columns)))
  session$setInputs(alignment_full=1L)
  session$flushReact()
  stopifnot(alignment_focus_state()$kind=="full",
    identical(as.integer(messages[["alignment_range"]]$value),c(1L,alignment_index()$n_columns)))
  messages[["alignment_range"]] <- NULL
  session$setInputs(examples_only=TRUE,alignment_scope="shown",compare_pairs=character())
  session$flushReact()
  stopifnot(is.null(messages[["alignment_range"]]),alignment_focus_state()$kind=="full")
  lr_pair <- candidates()$pair_id[candidates()$pair_label=="LR0R + LR3"]
  session$setInputs(pair_click=list(pair=lr_pair,nonce=4L))
  session$flushReact()
  coverage_text <- alignment_domain_coverage()
  stopifnot(grepl("fully covered: D1",coverage_text,fixed=TRUE),
    grepl("partial overlap: D2",coverage_text,fixed=TRUE))
  session$setInputs(product_coordinates="template_bases")
  session$flushReact()
  stopifnot(identical(alignment_domain_coverage(),coverage_text))
  session$setInputs(region_click=list(accession="PZ962702.1",landmark="28S_copy",nonce=5L))
  session$flushReact()
  copy_columns <- range(unname(alignment_index()$maps[["PZ962702.1"]][c(6000L,9600L)]))
  stopifnot(identical(as.integer(messages[["alignment_range"]]$value),as.integer(copy_columns)),
    focused_alignment_region()$landmark_id=="28S_copy",active_pair()$pair_id==lr_pair)
  detail_copy <- output$evidence$html
  stopifnot(grepl("including labeled partial features",detail_copy,fixed=TRUE),
    grepl("projected 28S domains",detail_copy,fixed=TRUE),
    grepl("independently of those annotations",detail_copy,fixed=TRUE),
    !grepl("only where NCBI supplies an unambiguous complete boundary",detail_copy,fixed=TRUE))
})
message("Nested source-backed domain rendering, provenance, raw-context projection and real module selection checks passed.")
