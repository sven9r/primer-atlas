#!/usr/bin/env Rscript
.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
suppressPackageStartupMessages({library(shiny); library(bslib); library(DT)})
source("R/functions.R")
source("R/data_layer.R")
source("R/fungal_explorer.R")
source("R/fungal_versions.R")
source("R/fungal_alignment.R")
source("R/fungal_alignment_ui.R")
Sys.unsetenv("ATLAS_MANIFEST_BASE_URL")
versions <- fungal_read_index(remote=FALSE)
old_id <- versions$id[startsWith(versions$id,"v1-")][1L]
current_id <- tail(versions$id,1L)
frozen <- setNames(lapply(c(old_id,current_id),function(id) fungal_read_snapshot(id,remote=FALSE)),c(old_id,current_id))
# Add an unavailable-evidence case only to this in-memory fixture. Current
# frozen full-locus templates happen to support every candidate somewhere;
# handoff must also remain honest when a future frozen candidate has no sites.
withheld <- frozen[[current_id]]$pairs$pair_id[frozen[[current_id]]$pairs$evidence_status=="candidate_published_components"][1L]
rows <- frozen[[current_id]]$pair_sequences$pair_id==withheld
frozen[[current_id]]$pair_sequences$status[rows] <- "fixture_withheld"
for(column in c("forward_status","reverse_status")) frozen[[current_id]]$pair_sequences[[column]][rows] <- "unavailable"
for(column in c("forward_start","forward_end","reverse_start","reverse_end","amplicon_bp"))
  frozen[[current_id]]$pair_sequences[[column]][rows] <- NA_real_
env <- new.env(parent=globalenv())
env$fungal_read_index <- function(...) versions
env$fungal_read_snapshot <- function(version_id,...) {
  if(!version_id %in% names(frozen)) stop("Unavailable fixture version")
  frozen[[version_id]]
}
source("R/fungal_ui.R",local=env)

# Closed secondary views must retain their controls and complete exports. The
# chart and cross-page action stay in the primary workspace, without a second
# stylesheet load from this namespaced module.
ui <- as.character(env$fungalExplorerUI("inspect"))
details <- regmatches(ui,gregexpr("<details[^>]*>",ui))[[1L]]
stopifnot(length(details)>=8L,!any(grepl(" open",details,fixed=TRUE)),
  grepl('id="inspect-open_fit"',ui,fixed=TRUE),
  grepl('id="inspect-active_combination"',ui,fixed=TRUE),
  !grepl("atlas-workspace.css",ui,fixed=TRUE))
for(id in c("download_manifest","download_pairs","download_sequences","download_alignment_sites",
    "download_alignment_regions","download_alignment")) {
  stopifnot(length(regmatches(ui,gregexpr(paste0('id="inspect-',id,'"'),ui,fixed=TRUE))[[1L]])==1L)
}

context <- reactiveVal(NULL)
messages <- new.env(parent=emptyenv())
mock <- shiny::MockShinySession$new()
mock$sendInputMessage <- function(inputId,message) messages[[inputId]] <- message
shiny::testServer(env$fungalExplorerServer,args=list(id="inspect",selection=context),session=mock,{
  session$setInputs(version=old_id,platform="ALL",phylum="ALL",class="ALL",locus="ALL",search="",
    examples_only=TRUE,show_unresolved=TRUE,compare_pairs=character())
  session$flushReact()
  original <- active_pair()$pair_id
  target <- frozen[[current_id]]$pairs$pair_id[frozen[[current_id]]$pairs$pair_label=="LR0R + LR3"]
  stopifnot(length(target)==1L,original!=target)
  context(list(version="missing-version",pair=target,nonce=1L))
  session$flushReact()
  stopifnot(snapshot()$manifest$id==old_id,active_pair()$pair_id==original,
    grepl("unavailable",output$selection_status$html,fixed=TRUE))
  context(list(version=current_id,pair="missing-pair",nonce=2L))
  session$flushReact()
  stopifnot(snapshot()$manifest$id==old_id,active_pair()$pair_id==original)

  # A test session receives updates but does not synthesize the asynchronous
  # browser response. Neither old version nor old filters may win the handoff.
  context(list(version=current_id,pair=target,nonce=3L))
  session$flushReact()
  stopifnot(messages$version$value==current_id,!is.null(pending_selection()),
    snapshot()$manifest$id==old_id,active_pair()$pair_id==original)
  session$setInputs(version=current_id)
  session$flushReact()
  stopifnot(!is.null(pending_selection()))
  session$setInputs(examples_only=FALSE)
  session$flushReact()
  stopifnot(snapshot()$manifest$id==current_id,active_pair()$pair_id==target,
    is.null(pending_selection()),is.null(requested_pair()),
    alignment_focus_state()$kind=="product",length(selected_accessions())==21L,
    nrow(candidates())>20L,nrow(alignment_pairs())==20L,
    identical(display_alignment_pairs()$pair_id,target))
  page_one <- output$quick_candidates$html
  card_ids <- function(html) regmatches(html,gregexpr('data-pair="[^"]+"',html))[[1L]]
  stopifnot(length(card_ids(page_one))==5L)
  session$setInputs(next_candidates=1L)
  session$flushReact()
  stopifnot(length(card_ids(output$quick_candidates$html))==5L,
    !any(card_ids(page_one) %in% card_ids(output$quick_candidates$html)))

  # A newer card click cancels the queued context dropdown acknowledgement.
  # After its own acknowledgement, direct dropdown choices remain usable.
  alternatives <- head(setdiff(candidates()$pair_id,target),2L)
  session$setInputs(pair_click=list(pair=alternatives[1L],nonce=1L))
  session$flushReact()
  session$setInputs(active_combination=target)
  session$flushReact()
  stopifnot(active_pair()$pair_id==alternatives[1L])
  session$setInputs(active_combination=alternatives[1L])
  session$flushReact()
  stopifnot(is.null(combination_ack()))
  session$setInputs(active_combination=alternatives[2L],compare_pairs=target)
  session$flushReact()
  stopifnot(active_pair()$pair_id==alternatives[2L],
    setequal(display_alignment_pairs()$pair_id,c(alternatives[2L],target)))
  session$setInputs(product_display="all")
  session$flushReact()
  stopifnot(nrow(display_alignment_pairs())==20L)

  # Closing tables never truncates the science exported from the current
  # candidate set or changes which immutable alignment is downloaded.
  pairs_path <- output$download_pairs
  exported <- read.csv(pairs_path,stringsAsFactors=FALSE)
  stopifnot(nrow(exported)==nrow(candidates()),setequal(exported$pair_id,candidates()$pair_id),
    all(exported$snapshot_id==current_id))
  aligned <- atlas_read_artifact(list(local_path=output$download_alignment))
  stopifnot(identical(aligned,snapshot()$alignment))

  # An explicit unsupported frozen pair is inspectable as unresolved evidence,
  # then ordinary user filters regain control instead of forcing that pin.
  unscored <- all_pairs()[all_pairs()$n_scorable==0L & all_pairs()$evidence_status!="documented_pair",,drop=FALSE]
  stopifnot(nrow(unscored)>0L)
  unresolved <- unscored$pair_id[1L]
  context(list(version=current_id,pair=unresolved,nonce=4L))
  session$flushReact()
  stopifnot(active_pair()$pair_id==unresolved,active_pair()$n_scorable==0L,
    active_pair()$n_panel==21L,nrow(alignment_focus_row())==1L,
    alignment_focus_row()$alignment_status!="mapped")
  session$setInputs(search="no_such_primer_987654321")
  session$flushReact()
  stopifnot(nrow(candidates())==0L,nrow(active_pair())==0L,is.null(requested_pair()))
})
message("Compact fungal explorer, five-card paging, focused comparisons, complete exports, frozen async handoff and newer-selection precedence passed.")
