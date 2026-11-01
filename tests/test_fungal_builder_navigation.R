#!/usr/bin/env Rscript
# Exercise the user's workflow in the real application, including an unmapped
# pair with genuine broad evidence and actual publication-download callbacks.
.libPaths(c(file.path(getwd(),".Rlib"),.libPaths()))
env <- new.env()
source("app.R",local=env)
downloads <- list()
env$downloadHandler <- function(filename,content,...) {
  downloads[[length(downloads)+1L]] <<- list(filename=filename,content=content)
  shiny::downloadHandler(filename=filename,content=content,...)
}
versions <- fungal_read_index(remote=FALSE)
v3 <- tail(versions$id,1L)
snapshot <- fungal_read_snapshot(v3,remote=FALSE)
component <- function(name,direction) {
  id <- snapshot$catalog$catalog_id[snapshot$catalog$primer_name==name & snapshot$catalog$direction==direction]
  stopifnot(length(id)==1L)
  id
}
shiny::testServer(env$server,{
  session$setInputs(map_marker="ITS_FUNGAL",map_organism="FUNGI",map_sort="binding_site",
    primer_overview_reference="PZ962702.1",primer_overview_region="ALL",primer_overview_platform="ALL",
    fungal_builder_min_bp=100,fungal_builder_max_bp=3500,fungal_builder_evidence="ALL",
    "fungal-version"=v3,"fungal-platform"="ALL","fungal-illumina_max"=500,
    "fungal-pacbio_max"=2000,"fungal-nanopore_max"=3500)
  session$flushReact()
  initial <- fungal_builder_result()
  stopifnot(length(fungal_builder_forward())==0L,length(fungal_builder_reverse())==0L,
    nrow(initial$selected_pairs)==0L,nrow(overview_inspector_rows())==0L,
    is.null(fungal_overview_evidence()))

  # The first view is an individual-primer atlas with an empty selection.
  # Every glyph keeps its original catalog ID and original role for subsequent
  # events; effective assignment and state change independently when reversed.
  atlas_ui <- as.character(fungalPrimerAtlasNavigatorUI())
  atlas <- fungal_atlas_data()
  atlas_html <- output$fungal_primer_atlas$html
  stopifnot(grepl('value="atlas"[^>]*checked',atlas_ui),
    is.data.frame(atlas$rows),nrow(atlas$rows)>0L,
    all(c('catalog_id','direction') %in% names(atlas$rows)),
    grepl('class="fungal-atlas-primer',atlas_html,fixed=TRUE),
    grepl('data-catalog-direction="forward"',atlas_html,fixed=TRUE),
    grepl('data-catalog-direction="reverse"',atlas_html,fixed=TRUE),
    grepl('data-selection-state="unselected"',atlas_html,fixed=TRUE),
    grepl('data-primer-input="fungal_atlas_primer_click"',atlas_html,fixed=TRUE),
    grepl(paste0('data-version-id="',v3,'"'),atlas_html,fixed=TRUE),
    grepl('data-accession="PZ962702.1"',atlas_html,fixed=TRUE),
    grepl('role="button"',atlas_html,fixed=TRUE),grepl('aria-pressed=',atlas_html,fixed=TRUE))
  atlas_nonce <- 100L
  builder_echo <- function(direction,selected=NULL) {
    if(is.null(selected)) selected <- if(direction=="forward") fungal_builder_forward() else fungal_builder_reverse()
    do.call(session$setInputs,setNames(list(selected),paste0("fungal_builder_",direction)))
    session$flushReact()
  }
  echo_pending <- function() {
    for(role in c("forward","reverse")) {
      expected <- isolate(fungal_builder_ack[[role]])
      if(!is.null(expected)) builder_echo(role,expected)
    }
  }
  click_primer <- function(id,direction,version=v3,accession="PZ962702.1",echo=TRUE) {
    atlas_nonce <<- atlas_nonce+1L
    session$setInputs(fungal_atlas_primer_click=list(catalog_id=id,parent_catalog_id=id,
      direction=direction,intent="cycle",version=version,accession=accession,nonce=atlas_nonce))
    session$flushReact()
    # A direction switch updates BOTH widgets. Echo every changed role,
    # including an empty selected vector, rather than only the clicked role.
    if(echo) echo_pending()
  }
  clear_selection <- function(nonce) {
    session$setInputs(fungal_builder_clear=nonce);session$flushReact()
    builder_echo("forward");builder_echo("reverse")
  }
  check_glyph <- function(parent,state,assigned,direction) {
    nodes <- xml2::xml_find_all(xml2::read_html(output$fungal_primer_atlas$html),
      paste0('.//*[@data-catalog-id="',parent,'"]'))
    source <- snapshot$catalog[snapshot$catalog$catalog_id==parent,,drop=FALSE]
    stopifnot(length(nodes)>0L,nrow(source)==1L,
      all(xml2::xml_attr(nodes,"data-assigned-catalog-id")==assigned),
      all(xml2::xml_attr(nodes,"data-direction")==direction),
      all(xml2::xml_attr(nodes,"data-catalog-direction")==source$direction),
      all(xml2::xml_attr(nodes,"data-selection-state")==state),
      all(xml2::xml_attr(nodes,"aria-pressed")==if(state=="unselected") "false" else "true"))
  }
  chosen_forward <- component("ITS1","forward")
  chosen_reverse <- component("ITS2","reverse")
  forward_variant <- fungal_orientation_resolve(snapshot,chosen_forward,"reverse")
  reverse_variant <- fungal_orientation_resolve(snapshot,chosen_reverse,"forward")
  stopifnot(chosen_forward %in% atlas$rows$catalog_id,chosen_reverse %in% atlas$rows$catalog_id)
  for(parent in list(list(id=chosen_forward,role="forward",variant=forward_variant),
                     list(id=chosen_reverse,role="reverse",variant=reverse_variant))) {
    check_glyph(parent$id,"unselected",parent$id,parent$role)
    click_primer(parent$id,parent$role)
    selected <- if(parent$role=="forward") fungal_builder_forward() else fungal_builder_reverse()
    stopifnot(identical(selected,parent$id),nrow(fungal_builder_result()$selected_pairs)==0L)
    check_glyph(parent$id,"catalog",parent$id,parent$role)
    click_primer(parent$id,parent$role)
    opposite <- if(parent$variant$direction=="forward") fungal_builder_forward() else fungal_builder_reverse()
    stopifnot(identical(opposite,parent$variant$catalog_id),!parent$id %in% c(fungal_builder_forward(),fungal_builder_reverse()),
      nrow(fungal_builder_result()$selected_pairs)==0L,
      parent$variant$sequence==reverse_complement(snapshot$catalog$sequence[match(parent$id,snapshot$catalog$catalog_id)]))
    check_glyph(parent$id,"reversed",parent$variant$catalog_id,parent$variant$direction)
    click_primer(parent$id,parent$role)
    stopifnot(length(fungal_builder_forward())==0L,length(fungal_builder_reverse())==0L,
      nrow(fungal_builder_result()$selected_pairs)==0L,
      is.null(fungal_builder_ack$forward),is.null(fungal_builder_ack$reverse))
    check_glyph(parent$id,"unselected",parent$id,parent$role)
  }

  # A real opposite-role oligo maps at the same raw site but makes a different
  # inward-facing product. ITS2 [RC] forward plus ITS4 reverse is 391 bp; it must
  # never inherit the documented ITS3+ITS4 name or publication support.
  its4 <- component("ITS4","reverse")
  click_primer(chosen_reverse,"reverse")
  click_primer(chosen_reverse,"reverse")
  click_primer(its4,"reverse")
  derived <- fungal_builder_result()
  stopifnot(nrow(derived$selected_pairs)==1L,derived$selected_pairs$has_derived_reverse_complement,
    derived$selected_pairs$evidence_status=="candidate_incomplete_provenance",
    derived$selected_pairs$forward_catalog_id==reverse_variant$catalog_id,
    derived$selected_pairs$forward_parent_catalog_id==chosen_reverse,
    derived$selected_pairs$forward_sequence==reverse_complement(reverse_variant$parent_sequence),
    nrow(overview_rows())==1L,overview_rows()$aligned_amplicon_bp==391L,
    overview_rows()$forward_start==2090L,overview_rows()$reverse_end==2480L,
    all(derived$pair_sequences$forward_site_origin=="derived_reverse_complement_parent_frozen_sites"))
  session$setInputs(primer_overview_pair_click=list(pair=derived$selected_pairs$pair_id,nonce=901L),
    fungal_publication_rank="class",fungal_publication_taxon="ALL",fungal_publication_title="Reversed primer candidate")
  session$flushReact()
  derived_support <- unite_explorer_workspace$support()
  stopifnot(fungal_overview_broad_selection()$available,
    unite_explorer_workspace$pair()==derived$selected_pairs$pair_id,
    sum(derived_support$summary$n_selected)==15402L,nrow(derived_support$details)==15402L,
    all(derived_support$details$version_id=="u1-afd349e6cdcd"))
  derived_export <- fungal_publication_build()
  stopifnot(nrow(derived_export$table)==1L,derived_export$table$forward_orientation_derived,
    derived_export$table$forward_parent_catalog_id==chosen_reverse,
    derived_export$table$forward_parent_sequence==reverse_variant$parent_sequence,
    derived_export$table$forward_parent_direction=="reverse",
    derived_export$table$forward_effective_direction=="forward",
    derived_export$table$forward_oligo_5to3==reverse_variant$sequence,
    derived_export$table$forward_reference_status=="derived_reverse_complement",
    derived_export$table$forward_parent_source_version==v3,
    derived_export$table$product_including_primers_bp==391L,
    derived_export$table$unite_selected_references==15402L,
    derived_export$table$combination_status=="candidate_incomplete_provenance")
  clear_selection(110L)

  # Three clicks can precede either widget acknowledgment. Stale original-F
  # and reversed-R echoes must not restore either assignment after removal.
  click_primer(chosen_forward,"forward",echo=FALSE)
  old_forward_echo <- fungal_builder_forward()
  click_primer(chosen_forward,"forward",echo=FALSE)
  old_reverse_echo <- fungal_builder_reverse()
  click_primer(chosen_forward,"forward",echo=FALSE)
  latest_forward <- fungal_builder_forward();latest_reverse <- fungal_builder_reverse()
  stopifnot(length(old_forward_echo)==1L,length(old_reverse_echo)==1L,
    length(latest_forward)==0L,length(latest_reverse)==0L,
    setequal(fungal_builder_ack$forward,latest_forward),setequal(fungal_builder_ack$reverse,latest_reverse))
  builder_echo("forward",old_forward_echo);builder_echo("reverse",old_reverse_echo)
  stopifnot(identical(fungal_builder_forward(),latest_forward),identical(fungal_builder_reverse(),latest_reverse),
    setequal(fungal_builder_ack$forward,latest_forward),setequal(fungal_builder_ack$reverse,latest_reverse))
  builder_echo("forward",latest_forward);builder_echo("reverse",latest_reverse)
  stopifnot(is.null(fungal_builder_ack$forward),is.null(fungal_builder_ack$reverse))
  # Ordinary dropdown edits become authoritative after the latest echoes. A
  # direct assignment to the other role removes the same parent's prior role.
  builder_echo("forward",old_forward_echo)
  stopifnot(identical(fungal_builder_forward(),old_forward_echo),length(fungal_builder_reverse())==0L)
  builder_echo("reverse",old_reverse_echo);echo_pending()
  stopifnot(length(fungal_builder_forward())==0L,identical(fungal_builder_reverse(),old_reverse_echo),
    is.null(fungal_builder_ack$forward),is.null(fungal_builder_ack$reverse))
  builder_echo("forward",old_forward_echo);echo_pending()
  stopifnot(identical(fungal_builder_forward(),old_forward_echo),length(fungal_builder_reverse())==0L)

  # The teaching panel is now an explicit action, not the initial selection.
  session$setInputs(fungal_builder_examples=101L);session$flushReact()
  builder_echo("forward");builder_echo("reverse")
  defaults <- fungal_builder_result()
  stopifnot(length(fungal_builder_forward())==5L,length(fungal_builder_reverse())==5L,
    nrow(defaults$selected_pairs)==25L,any(!defaults$selected_pairs$pair_id %in% snapshot$pairs$pair_id),
    nrow(overview_inspector_rows())==25L,
    grepl("25[[:space:]]+requested combinations",output$primer_overview_context$html),
    is.null(fungal_builder_ack$forward),is.null(fungal_builder_ack$reverse))

  # Stale/foreign events and incorrect original roles cannot alter selections.
  retained_forward <- fungal_builder_forward();retained_reverse <- fungal_builder_reverse()
  other_version <- versions$id[versions$id!=v3][1L]
  click_primer(chosen_forward,"forward",version=other_version)
  click_primer(chosen_forward,"forward",accession="FN812768.2")
  click_primer(chosen_forward,"reverse")
  click_primer("NOT_A_CATALOG_ID","forward")
  stopifnot(identical(fungal_builder_forward(),retained_forward),identical(fungal_builder_reverse(),retained_reverse))
  extra_forward <- setdiff(atlas$rows$catalog_id[atlas$rows$direction=="forward"],retained_forward)[1L]
  extra_reverse <- setdiff(atlas$rows$catalog_id[atlas$rows$direction=="reverse"],retained_reverse)[1L]
  stopifnot(!is.na(extra_forward),!is.na(extra_reverse))
  click_primer(extra_forward,"forward");click_primer(extra_reverse,"reverse")
  # Reversing a selected oligo into a full destination is rejected atomically;
  # its original assignment must survive the attempted direction switch.
  click_primer(chosen_forward,"forward");click_primer(chosen_reverse,"reverse")
  stopifnot(identical(fungal_builder_forward(),retained_forward),identical(fungal_builder_reverse(),retained_reverse))

  pair_ids <- fungal_builder_result()$selected_pairs$pair_id
  session$setInputs(fungal_overview_mode="amplicons");session$flushReact()
  stopifnot(identical(fungal_builder_forward(),retained_forward),identical(fungal_builder_reverse(),retained_reverse),
    identical(fungal_builder_result()$selected_pairs$pair_id,pair_ids))
  session$setInputs(fungal_overview_mode="atlas");session$flushReact()
  stopifnot(identical(fungal_builder_forward(),retained_forward),identical(fungal_builder_reverse(),retained_reverse),
    identical(fungal_builder_result()$selected_pairs$pair_id,pair_ids))

  clear_selection(100L)
  stopifnot(nrow(fungal_builder_result()$selected_pairs)==0L,nrow(fungal_atlas_data()$rows)==nrow(atlas$rows),
    grepl('data-primer-input="fungal_atlas_primer_click"',output$fungal_primer_atlas$html,fixed=TRUE))
  click_primer(chosen_forward,"forward");click_primer(chosen_reverse,"reverse")
  stopifnot(nrow(fungal_builder_result()$selected_pairs)==1L,
    fungal_builder_result()$selected_pairs$pair_label=="ITS1 + ITS2")
  click_primer(chosen_forward,"forward")
  stopifnot(length(fungal_builder_forward())==0L,length(fungal_builder_reverse())==2L,
    nrow(fungal_builder_result()$selected_pairs)==0L)
  click_primer(chosen_forward,"forward")
  session$setInputs(fungal_overview_mode="amplicons");session$flushReact()
  session$setInputs(fungal_overview_mode="atlas");session$flushReact()
  stopifnot(length(fungal_builder_forward())==0L,length(fungal_builder_reverse())==1L,
    nrow(fungal_builder_result()$selected_pairs)==0L)
  session$setInputs(fungal_builder_examples=100L,fungal_overview_mode="amplicons",primer_overview_region="ITS")
  session$flushReact()
  builder_echo("forward");builder_echo("reverse")
  stopifnot(nrow(fungal_builder_result()$selected_pairs)==25L)

  # F1 × R1/R2/R3 produces exactly three choices, not a documented-pair subset.
  f <- component("ITS1F","forward")
  r <- vapply(c("ITS2","ITS4","LR21"),component,character(1L),direction="reverse")
  session$setInputs(fungal_builder_forward=f,fungal_builder_reverse=r)
  session$flushReact()
  result <- fungal_builder_result()
  stopifnot(nrow(result$selected_pairs)==3L,nrow(overview_inspector_rows())==3L,
    all(result$selected_pairs$forward_primer=="ITS1F"),
    setequal(result$selected_pairs$reverse_primer,c("ITS2","ITS4","LR21")))
  observed <- result$pair_sequences[result$pair_sequences$accession=="PZ962702.1",,drop=FALSE]
  before <- overview_rows()
  stopifnot(nrow(before)==3L,all(before$aligned_amplicon_bp==before$reverse_end-before$forward_start+1L),
    all(before$targeted_region_bp==before$reverse_start-before$forward_end-1L))
  session$setInputs(fungal_builder_max_bp=500)
  session$flushReact()
  stopifnot(nrow(fungal_builder_result()$selected_pairs)==3L,nrow(overview_inspector_rows())==3L,
    nrow(overview_rows())<3L,all(overview_rows()$aligned_amplicon_bp<=500L),
    sum(grepl("Outside current filters",output$fungal_builder_combination_table$html,fixed=TRUE))==1L)
  session$setInputs(fungal_builder_max_bp=3500,primer_overview_platform="ILLUMINA")
  session$flushReact()
  stopifnot(all(overview_rows()$aligned_amplicon_bp<=500L))
  session$setInputs(primer_overview_platform="ALL")
  session$flushReact()

  # All chosen primers survive a clicked lane and its advanced-view handoff.
  pair <- result$selected_pairs$pair_id[1L]
  session$setInputs(primer_overview_pair_click=list(pair=pair,nonce=1L))
  session$flushReact()
  stopifnot(length(fungal_builder_forward())==1L,length(fungal_builder_reverse())==3L,
    overview_selected_pair()==pair,unite_explorer_workspace$pair()==pair)

  # This pair is absent from the old pair list and has no product on this map
  # reference, yet precomputed UNITE sites support products in other taxa.
  session$setInputs(fungal_builder_forward=component("ITS1","forward"),
    fungal_builder_reverse=component("ITS4B","reverse"))
  session$flushReact()
  missing <- fungal_builder_result()$selected_pairs
  stopifnot(nrow(missing)==1L,!missing$pair_id %in% snapshot$pairs$pair_id,nrow(overview_rows())==0L,
    nrow(overview_inspector_rows())==1L)
  session$setInputs(primer_overview_pair_click=list(pair=missing$pair_id,nonce=2L))
  session$flushReact()
  support <- unite_explorer_workspace$support()
  stopifnot(overview_selected_pair()==missing$pair_id,unite_explorer_workspace$pair()==missing$pair_id,
    sum(support$summary$n_selected)==15402L,sum(support$summary$n_exact_products)==12L,
    nrow(support$details)==15402L)

  # Explicit Clear must remove the retained inspector context, as well as the
  # map choices. The previous missing-map pair must not survive in visible fit.
  session$setInputs(fungal_builder_clear=1L)
  session$flushReact()
  stopifnot(length(fungal_builder_forward())==0L,length(fungal_builder_reverse())==0L,
    nrow(fungal_builder_result()$selected_pairs)==0L,is.null(overview_clicked_pair()),
    is.null(fungal_overview_evidence()),!fungal_overview_broad_selection()$available,
    is.null(fungal_overview_broad_selection()$pair))
  session$setInputs(fungal_builder_forward=NULL,fungal_builder_reverse=NULL)
  session$flushReact()
  session$setInputs(fungal_builder_forward=component("ITS1","forward"),
    fungal_builder_reverse=component("ITS4B","reverse"))
  session$flushReact()
  session$setInputs(primer_overview_pair_click=list(pair=missing$pair_id,nonce=3L))
  session$flushReact()
  stopifnot(overview_selected_pair()==missing$pair_id,fungal_overview_broad_selection()$available)

  # Clearing the final oligo must stay empty instead of restoring defaults.
  session$setInputs(fungal_builder_forward=NULL)
  session$flushReact()
  stopifnot(length(fungal_builder_forward())==0L,nrow(fungal_builder_result()$selected_pairs)==0L,
    nrow(overview_rows())==0L,is.null(fungal_overview_evidence()),!fungal_overview_broad_selection()$available)
  session$setInputs(fungal_builder_examples=1L)
  session$flushReact()
  stopifnot(nrow(fungal_builder_result()$selected_pairs)==25L)

  # Export callbacks use the same three requested pairs and filter state.
  session$setInputs(fungal_builder_forward=f,fungal_builder_reverse=r,fungal_builder_max_bp=500,
    fungal_publication_rank="class",fungal_publication_taxon="ALL",fungal_publication_title="ITS primer comparison")
  session$flushReact()
  exported <- fungal_publication_build()
  stopifnot(nrow(exported$table)==3L,length(unique(exported$table$pair_id))==3L,
    grepl("ITS primer comparison",exported$svg,fixed=TRUE))
  output$fungal_download_table
  output$fungal_download_package
  # Probe only this workspace's registered download callbacks. Filename
  # functions for unrelated COI/BeePrime inspectors require their own inputs.
  publication_downloads <- Filter(function(download) is.function(download$filename) && grepl("fungal_publication_filename",
    paste(deparse(body(download$filename)),collapse="\n"),fixed=TRUE),downloads)
  filename <- paste0(v3,"-primer-combinations.csv")
  matches <- which(vapply(publication_downloads,function(download) identical(download$filename(),filename),logical(1L)))
  stopifnot(length(matches)==1L)
  file <- tempfile(fileext=".csv")
  publication_downloads[[matches]]$content(file)
  table <- read.csv(file,stringsAsFactors=FALSE)
  stopifnot(nrow(table)==3L,setequal(table$pair_id,result$selected_pairs$pair_id))
  unlink(file)
  zip_matches <- which(vapply(publication_downloads,function(download) identical(download$filename(),paste0(v3,"-primer-combinations.zip")),logical(1L)))
  stopifnot(length(zip_matches)==1L)
  archive <- tempfile(fileext=".zip")
  publication_downloads[[zip_matches]]$content(archive)
  members <- utils::unzip(archive,list=TRUE)$Name
  stopifnot(all(c("primer-combinations.svg","primer-combinations.csv","primer-combinations.html",
    "caption.txt","methods.txt","citations.txt","requested-pairs.csv","long-rDNA-template-evidence.csv") %in% members))
  unlink(archive)

  # The fit request arrives before its historical snapshot input is delivered
  # by the browser. Choose a pair outside that version's ordinary defaults, so
  # a default reset cannot accidentally satisfy the handoff assertion.
  historical_id <- versions$id[startsWith(versions$id,"v1-")][1L]
  historical <- fungal_read_snapshot(historical_id,remote=FALSE)
  old_defaults <- fungal_combination_defaults(historical,accession="FN812768.2",region="ITS")
  old_component <- function(pair,direction) {
    ids <- historical$catalog$catalog_id[historical$catalog$primer_name==pair[[paste0(direction,"_primer")]] &
      historical$catalog$sequence==pair[[paste0(direction,"_sequence")]] & historical$catalog$direction==direction]
    stopifnot(length(ids)==1L)
    ids
  }
  forwards <- vapply(seq_len(nrow(historical$pairs)),function(i) old_component(historical$pairs[i,],"forward"),character(1L))
  reverses <- vapply(seq_len(nrow(historical$pairs)),function(i) old_component(historical$pairs[i,],"reverse"),character(1L))
  outside_defaults <- which(!forwards %in% old_defaults$forward_ids | !reverses %in% old_defaults$reverse_ids)
  stopifnot(length(outside_defaults)>0L)
  chosen <- outside_defaults[1L]
  old_pair <- historical$pairs$pair_id[chosen]
  session$setInputs("fungal_fit-version"=historical_id,"fungal_fit-pair"=old_pair,
    "fungal_fit-rank"="class","fungal_fit-platform"="ALL")
  session$flushReact()
  stopifnot(fungal_fit_workspace$snapshot()$manifest$id==historical_id,
    fungal_workspace$snapshot()$manifest$id==v3)
  session$setInputs("fungal_fit-open_explorer"=1L)
  session$flushReact()
  stopifnot(fungal_explorer_context()$version==historical_id,fungal_explorer_context()$pair==old_pair,
    fungal_workspace$snapshot()$manifest$id==v3)
  session$setInputs("fungal-version"=historical_id,primer_overview_version=historical_id,
    primer_overview_reference="FN812768.2")
  session$flushReact()
  stopifnot(fungal_builder_state$version==historical_id,
    length(fungal_builder_forward())==1L,length(fungal_builder_reverse())==1L,
    forwards[chosen] %in% fungal_builder_forward(),reverses[chosen] %in% fungal_builder_reverse(),
    nrow(fungal_builder_result()$selected_pairs)==1L,
    old_pair %in% fungal_builder_result()$selected_pairs$pair_id,
    old_pair %in% overview_inspector_rows()$pair_id)
})
message("Empty initial selection, three-click primer orientation cycles, guarded selector acknowledgments, Cartesian counts, raw length filters, missing-map taxon handoff, clear/history reset and publication downloads passed.")
