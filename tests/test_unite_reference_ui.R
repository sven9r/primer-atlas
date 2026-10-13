#!/usr/bin/env Rscript
.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
suppressPackageStartupMessages({library(shiny);library(DT)})
source("R/functions.R")
source("R/data_layer.R")
source("R/fungal_explorer.R")
source("R/unite_reference.R")
source("R/unite_reference_ui.R")
plain_text <- function(html) trimws(gsub("[[:space:]]+"," ",gsub("<[^>]*>"," ",html)))

# Synthetic source records exercise real indexed-site and regional-alignment
# helpers. No installed snapshot or scientific data is modified.
scratch <- tempfile("unite-ui-")
dir.create(scratch)
ids <- sprintf("SH_fixture_%02d",seq_len(15L))
set.seed(1313)
f <- "GCGATCGTACGATGCCTAGC"
r <- "AAGCATCCTCGGCTATCGGA"
absent_r <- "TTAGGACCGATGTACCGGTA"
sequences <- setNames(vapply(ids,function(id) paste0(f,paste(sample(c("A","C","G","T"),140L,replace=TRUE),collapse=""),reverse_complement(r)),character(1)),ids)
substr(sequences[2L],161L,180L) <- paste(rep("A",20L),collapse="")
substr(sequences[3L],1L,20L) <- paste(rep("A",20L),collapse="")
substr(sequences[3L],161L,180L) <- paste(rep("A",20L),collapse="")
substr(sequences[4L],61L,80L) <- f
sequences[5L] <- "ACGT"
substr(sequences[6L],61L,61L) <- "N"
metadata <- data.frame(record_id=ids,accession=paste0("fixture-accession-",seq_along(ids)),
  organism=paste("Fixture organism",seq_along(ids)),species_hypothesis=paste0("SH-fixture.",seq_along(ids)),phylum="Fixture phylum",
  class_name=c(rep("Fixture class A",14L),"Unresolved class"),
  order_name=c(rep("Fixture order A",7L),rep("Fixture order B",7L),"Unresolved order"),
  family="",genus="",species="",class_status=c(rep("resolved",14L),"incertae_sedis"),
  class_group_id=c(rep("c-A",14L),"c-U"),order_group_id=c(rep("o-A",7L),rep("o-B",7L),"o-U"),
  phylum_group_id="p-A",raw_length_bp=nchar(sequences),stringsAsFactors=FALSE)
coverage <- data.frame(rank=c("class","class","order","order","order","phylum"),
  taxon=c("Fixture class A","Unresolved class","Fixture order A","Fixture order B","Unresolved order","Fixture phylum"),
  parent_path=c("Fixture phylum","Fixture phylum","Fixture phylum; Fixture class A","Fixture phylum; Fixture class A","Fixture phylum; Unresolved class",""),
  group_id=c("c-A","c-U","o-A","o-B","o-U","p-A"),
  n_source=c(102134L,3L,50000L,52134L,3L,102137L),
  n_eligible_unique=c(400L,1L,200L,200L,1L,401L),n_selected=c(14L,1L,7L,7L,1L,15L),
  shortfall_to_100=c(NA_integer_,NA_integer_,93L,93L,99L,NA_integer_),taxon_resolved=c(TRUE,FALSE,TRUE,TRUE,FALSE,TRUE),
  taxon_status=c("resolved","incertae_sedis","resolved","resolved","unresolved","resolved"),
  status=c("bounded_subset","sparse_source","bounded_subset","bounded_subset","sparse_source","bounded_subset"),stringsAsFactors=FALSE)
catalog <- data.frame(primer_name=c("Fixture F","Fixture R","Fixture absent R"),
  sequence=c(f,r,absent_r),direction=c("forward","reverse","reverse"),
  catalog_id=c("f","r","absent-r"),stringsAsFactors=FALSE)
pairs <- data.frame(pair_id=c("F_R","F_ABSENT"),pair_label=c("Fixture F + Fixture R","Fixture F + Fixture absent R"),
  forward_primer="Fixture F",reverse_primer=c("Fixture R","Fixture absent R"),forward_sequence=f,
  reverse_sequence=c(r,absent_r),evidence_status=c("documented_pair","candidate_published_components"),stringsAsFactors=FALSE)
sites <- do.call(rbind,lapply(seq_along(ids),function(i) do.call(rbind,lapply(seq_len(nrow(catalog)),function(j) {
  found <- .fungal_exact_site(strsplit(sequences[i],"",fixed=TRUE)[[1L]],catalog$sequence[j],catalog$direction[j])
  data.frame(record_id=ids[i],catalog_id=catalog$catalog_id[j],start=found$start,end=found$end,status=found$status,stringsAsFactors=FALSE)
}))))
complete <- ids[c(1:4,6:15)]
regions <- do.call(rbind,lapply(c("ITS1","5.8S","ITS2"),function(region) {
  bounds <- switch(region,ITS1=c(1L,140L),`5.8S`=c(60L,69L),ITS2=c(41L,180L))
  data.frame(record_id=complete,region=region,start=bounds[1L],end=bounds[2L],boundary_status="complete")
}))
index <- do.call(rbind,lapply(c("ITS1","5.8S","ITS2"),function(region) {
  rows <- regions[regions$region==region & regions$record_id %in% ids[1:14],,drop=FALSE]
  raw <- setNames(vapply(seq_len(nrow(rows)),function(i) substr(sequences[[rows$record_id[i]]],rows$start[i],rows$end[i]),character(1)),rows$record_id)
  aligned <- if(region=="5.8S") paste0(substr(raw,1L,3L),"-",substr(raw,4L,10L)) else
    paste0(substr(raw,1L,10L),"---",substr(raw,11L,40L),"----------",substr(raw,41L,140L))
  names(aligned) <- names(raw)
  file <- paste0(region,".fasta")
  writeLines(unlist(lapply(seq_along(aligned),function(i) c(paste0(">",names(aligned)[i]),aligned[i])),use.names=FALSE),file.path(scratch,file))
  data.frame(rank="class",group_id="c-A",region=region,file=file,stringsAsFactors=FALSE)
}))
order_index <- do.call(rbind,lapply(c("o-A","o-B"),function(group) {
  region <- if(group=="o-A") "ITS1" else "ITS2"
  owners <- metadata$record_id[metadata$order_group_id==group]
  rows <- regions[regions$region==region & regions$record_id %in% owners,,drop=FALSE]
  raw <- setNames(vapply(seq_len(nrow(rows)),function(i) substr(sequences[[rows$record_id[i]]],rows$start[i],rows$end[i]),character(1)),rows$record_id)
  # These files intentionally have different gap placement from the class MSA.
  # A wrong class owner or cached class file cannot satisfy the assertions.
  gap <- if(group=="o-A") "-----" else "-------"
  aligned <- setNames(paste0(substr(raw,1L,25L),gap,substr(raw,26L,140L)),names(raw))
  file <- paste0(group,"-",region,".fasta")
  writeLines(unlist(lapply(seq_along(aligned),function(i) c(paste0(">",names(aligned)[i]),aligned[i])),use.names=FALSE),file.path(scratch,file))
  data.frame(rank="order",group_id=group,region=region,file=file,stringsAsFactors=FALSE)
}))
index <- rbind(index,order_index)
fixture <- list(manifest=list(id="unite-ui-frozen-fixture",n_source=102137L,n_selected=15L,
  source_version="Synthetic UI fixture; no biological boundary claim",source_doi="10.0000/fixture",
  created_at="2026-10-13",source_taxonomy="Fixture taxonomy",region_annotation="Synthetic source-derived regional intervals"),
  directory=scratch,metadata=metadata,coverage=coverage,regions=regions,alignment_index=index,
  sites=sites,catalog=catalog,pairs=pairs,sequences=sequences,
  site_keys=paste(sites$record_id,sites$catalog_id,sep="\n"))

ui <- as.character(uniteReferenceUI("unite"))
for(mode in c("coverage","alignment","primers")) stopifnot(grepl(paste0('value="',mode,'"'),ui,fixed=TRUE))
stopifnot(grepl("output.available",ui,fixed=TRUE),
  grepl('id="unite-mode"',ui,fixed=TRUE))
inputs <- regmatches(ui,gregexpr("<input[^>]*>",ui))[[1L]]
stopifnot(any(grepl('value="alignment"',inputs,fixed=TRUE) & grepl('checked',inputs,fixed=TRUE)))
primer_ui <- as.character(uniteReferenceUI("unite",initial_mode="primers"))
primer_inputs <- regmatches(primer_ui,gregexpr("<input[^>]*>",primer_ui))[[1L]]
stopifnot(any(grepl('value="primers"',primer_inputs,fixed=TRUE) & grepl('checked',primer_inputs,fixed=TRUE)))

messages <- new.env(parent=emptyenv())
registered_choices <- new.env(parent=emptyenv())
registered_filters <- new.env(parent=emptyenv())
mock <- shiny::MockShinySession$new()
mock$sendInputMessage <- function(inputId,message) messages[[inputId]] <- message
mock$registerDataObj <- function(name,data,filterFunc) {
  registered_choices[[name]] <- data
  registered_filters[[name]] <- filterFunc
  paste0("/fixture-data/",name)
}
provider <- reactiveVal(NULL)
shiny::testServer(uniteReferenceServer,args=list(id="unite",dataset=function() fixture,pair_provider=provider),session=mock,{
  session$setInputs(breadth_rank="class",breadth_search="",alignment_rank="class",alignment_group="c-A",
    region="ITS1",pair="F_R",sites_rank="class",mode="alignment")
  session$flushReact()
  stopifnot(available(),mode()=="alignment",output$available=="true",is.null(output$availability),
    grepl("102,137",output$source_context$html,fixed=TRUE),
    grepl("15 selected references",gsub("[[:space:]]+"," ",gsub("<[^>]*>"," ",output$source_context$html)),fixed=TRUE),
    !grepl("short of 100",output$breadth$html,fixed=TRUE),
    grepl("Incertae sedis",output$breadth$html,fixed=TRUE))
  session$setInputs(breadth_rank="order")
  session$flushReact()
  stopifnot(grepl("93 short of 100 selected",output$breadth$html,fixed=TRUE),
    grepl("Taxonomy unresolved",output$breadth$html,fixed=TRUE),
    grepl("Click an order",output$breadth_instruction$html,fixed=TRUE))
  session$setInputs(breadth_rank="class")
  session$flushReact()
  stopifnot(length(alignment()$alignment)==12L,alignment()$n_aligned==13L,
    alignment()$n_total==14L,alignment()$n_regions==13L,alignment()$n_missing_boundaries==1L,
    nchar(window()$alignment[1L])==80L,!ids[5L] %in% names(window()$alignment),
    grepl("12 displayed of 13",plain_text(output$alignment_context$html),fixed=TRUE),
    grepl('data-state="gap"',output$alignment_plot$html,fixed=TRUE),
    grepl('data-base="N" data-state="ambiguous"',output$alignment_plot$html,fixed=TRUE),
    !grepl(paste0('class="unite-alignment-row" data-record-id="',ids[5L],'"'),output$alignment_plot$html,fixed=TRUE),
    grepl(ids[5L],output$missing_regions$html,fixed=TRUE))
  full <- alignment()$alignment
  session$setInputs(window_width=120L)
  session$flushReact()
  stopifnot(window()$n_window==120L)
  session$setInputs(window_start=999L)
  session$flushReact()
  stopifnot(window()$start==34L,window()$end==153L,
    identical(unname(window()$alignment),unname(substr(full,34L,153L))),
    grepl('data-start-column="34"',output$alignment_plot$html,fixed=TRUE))
  session$setInputs(region="5.8S")
  session$flushReact()
  stopifnot(window()$n_window==11L,window()$start==1L,
    grepl("All 11 alignment columns",output$window_control$html,fixed=TRUE),
    identical(gsub("-","",window()$alignment,fixed=TRUE),setNames(substr(sequences[names(window()$alignment)],60L,69L),names(window()$alignment))))
  session$setInputs(region="ITS2")
  session$flushReact()
  exported_alignment <- atlas_read_artifact(list(local_path=output$download_alignment))
  stopifnot(identical(exported_alignment,alignment()$alignment),length(exported_alignment)==12L)
  selected_counts <- read.csv(output$download_breadth,stringsAsFactors=FALSE)
  stopifnot(nrow(selected_counts)==2L,sum(selected_counts$n_selected)==15L,
    sum(selected_counts$n_source)==102137L,all(is.na(selected_counts$shortfall_to_100)))
  session$setInputs(breadth_search="Unresolved")
  session$flushReact()
  filtered <- read.csv(output$download_breadth,stringsAsFactors=FALSE)
  stopifnot(nrow(groups())==1L,nrow(filtered)==1L,filtered$group_id=="c-U",
    grepl("102,137",output$source_context$html,fixed=TRUE))
  session$setInputs(group_click=list(group="c-U",rank="class",nonce=1L))
  session$flushReact()
  stopifnot(alignment_group()=="c-U",messages$mode$value=="alignment",
    length(alignment()$alignment)==0L,grepl("No regional alignment",output$alignment_plot$html,fixed=TRUE),
    !grepl("complete region boundaries of",output$alignment_context$html,fixed=TRUE))
  session$setInputs(alignment_group="c-A")
  session$flushReact()
  stopifnot(alignment_group()=="c-U")
  session$setInputs(alignment_group="c-U")
  session$flushReact()
  session$setInputs(alignment_group="c-A")
  session$flushReact()
  stopifnot(alignment_group()=="c-A",length(alignment()$alignment)==12L)
  # A clicked order must load its own genuine regional MSA, not a class file.
  session$setInputs(breadth_search="",breadth_rank="order",region="ITS2")
  session$flushReact()
  session$setInputs(group_click=list(group="o-B",rank="order",nonce=2L))
  session$flushReact()
  stopifnot(alignment_rank()=="order",alignment_group()=="o-B",messages$alignment_rank$value=="order",
    length(alignment()$alignment)==7L,alignment()$n_total==7L,alignment()$n_columns==147L,
    setequal(names(alignment()$alignment),metadata$record_id[metadata$order_group_id=="o-B"]),
    all(alignment()$metadata$order_group_id=="o-B"),!ids[1L] %in% names(alignment()$alignment))
  order_alignment <- atlas_read_artifact(list(local_path=output$download_alignment))
  stopifnot(identical(order_alignment,alignment()$alignment),all(nchar(order_alignment)==147L))
  session$setInputs(alignment_rank="order",alignment_group="o-B")
  session$flushReact()
  session$setInputs(alignment_group="o-U")
  session$flushReact()
  stopifnot(alignment_group()=="o-U",alignment()$n_total==1L,length(alignment()$alignment)==0L)

  # Summary state counts must partition the retained records. A full exact bar
  # does not change source breadth or turn missing evidence into PCR failure.
  actual <- support()
  stopifnot(grepl("Documented pair",output$site_tally$html,fixed=TRUE),
    !grepl("Candidate combination",output$site_tally$html,fixed=TRUE))
  stopifnot(sum(actual$summary$n_selected)==15L,sum(actual$summary$n_exact_products)==11L,
    sum(actual$summary$n_single_sites)==1L,sum(actual$summary$n_no_exact_sites)==1L,
    sum(actual$summary$n_unresolved)==1L,sum(actual$summary$n_unavailable)==1L,
    all(rowSums(actual$summary[c("n_exact_products","n_single_sites","n_no_exact_sites","n_unresolved","n_unavailable","n_withheld")])==actual$summary$n_selected))
  session$setInputs(site_summary_rows_selected=1L)
  session$flushReact()
  stopifnot(nrow(site_details())==14L,all(site_details()$group_id=="c-A"))
  all_evidence <- read.csv(output$download_sites,stringsAsFactors=FALSE)
  stopifnot(nrow(all_evidence)==15L,setequal(all_evidence$record_id,ids),
    all(all_evidence$version_id==fixture$manifest$id),all(all_evidence$pair_id=="F_R"))
  provider("F_ABSENT")
  session$flushReact()
  stopifnot(pair()=="F_ABSENT",nrow(site_details())==0L,sum(support()$summary$n_exact_products)==0L,
    grepl("Candidate combination",output$site_tally$html,fixed=TRUE),
    grepl("independent validation",output$site_tally$html,fixed=TRUE))
  # Selection-only server updates must not overwrite the AJAX catalog with
  # a zero-column frame; later searches still need the original labels/IDs.
  keys <- grep("(^|-)pair$",ls(registered_choices),value=TRUE)
  stopifnot(length(keys)==1L)
  choices <- registered_choices[[keys]]
  stopifnot(is.data.frame(choices),all(c("label","value") %in% names(choices)),
    nrow(choices)==nrow(pairs),setequal(choices$value,pairs$pair_id))
  response <- registered_filters[[keys]](choices,list(QUERY_STRING=
    "query=&field=%5B%22label%22%5D&value=value&conju=and&maxop=2"))
  returned <- jsonlite::fromJSON(response$content)
  stopifnot(response$status==200L,setequal(returned$value,pairs$pair_id))
  session$setInputs(pair="F_ABSENT")
  session$flushReact()
  session$setInputs(pair="F_R",sites_rank="order")
  session$flushReact()
  stopifnot(pair()=="F_R",nrow(support()$summary)==3L,nrow(site_details())==0L,
    all(support()$details$group_id==metadata$order_group_id))
})

# Without an explicit group, Orders is the default and the chosen region
# prefers a group that has an actual indexed alignment file.
shiny::testServer(uniteReferenceServer,args=list(id="default",dataset=function() fixture),{
  session$flushReact()
  stopifnot(alignment_rank()=="order",alignment_group()=="o-A",alignment()$n_columns==145L,
    window()$n_window==80L,all(alignment()$metadata$order_group_id=="o-A"))
  session$setInputs(region="ITS2")
  session$flushReact()
  stopifnot(alignment_group()=="o-B",alignment()$n_columns==147L)
  session$setInputs(window_width=60L)
  session$flushReact()
  stopifnot(window()$n_window==60L,nchar(window()$alignment[1L])==60L)
})

# Two separately named frozen cohorts are selectable inside source details.
# A real identity change resets old clade/pair context atomically; mode/rank
# edits and rereading the same identity must preserve a user's selections.
version_one <- fixture
version_two <- fixture
version_two$manifest$id <- "unite-ui-second-frozen-fixture"
version_two$metadata$organism <- paste(version_two$metadata$organism,"cohort two")
history <- data.frame(id=c(version_one$manifest$id,version_two$manifest$id),label=c("Fixture U1","Fixture U2"),
  manifest=c("u1/manifest.json","u2/manifest.json"),sha256=c(paste(rep("a",64L),collapse=""),paste(rep("b",64L),collapse="")),n_selected=15L)
version_one$versions <- version_two$versions <- history
version_reader <- function(version_id=NULL) {
  if(is.null(version_id) || identical(version_id,version_one$manifest$id)) return(version_one)
  if(identical(version_id,version_two$manifest$id)) return(version_two)
  stop("Unknown frozen fixture")
}
shiny::testServer(uniteReferenceServer,args=list(id="history",dataset=version_reader),{
  session$setInputs(version=version_one$manifest$id,alignment_rank="class",alignment_group="c-U",pair="F_ABSENT",region="ITS1")
  session$flushReact()
  stopifnot(data()$manifest$id==version_one$manifest$id,pair()=="F_ABSENT",alignment_group()=="c-U",
    grepl('id="history-version"',output$source_evidence$html,fixed=TRUE),
    !grepl('id="unite-version"',ui,fixed=TRUE))
  session$setInputs(mode="primers",sites_rank="order")
  session$flushReact()
  stopifnot(pair()=="F_ABSENT",alignment_group()=="c-U")
  session$setInputs(version=version_two$manifest$id)
  session$flushReact()
  stopifnot(data()$manifest$id==version_two$manifest$id,pair()=="F_R",alignment_rank()=="order",alignment_group()=="o-A",
    length(alignment()$alignment)==6L,all(grepl("cohort two",alignment()$metadata$organism,fixed=TRUE)),nrow(site_details())==0L)
  evidence <- read.csv(output$download_sites,stringsAsFactors=FALSE)
  stopifnot(all(evidence$version_id==version_two$manifest$id),all(evidence$pair_id=="F_R"),
    all(grepl("cohort two",evidence$organism,fixed=TRUE)))
  session$setInputs(alignment_rank="order",alignment_group="o-A",pair="F_R")
  session$flushReact()
  session$setInputs(alignment_rank="class",alignment_group="c-U",pair="F_ABSENT")
  session$flushReact()
  session$setInputs(version=version_two$manifest$id)
  session$flushReact()
  stopifnot(pair()=="F_ABSENT",alignment_group()=="c-U")
  session$setInputs(version=version_one$manifest$id)
  session$flushReact()
  stopifnot(data()$manifest$id==version_one$manifest$id,pair()=="F_R",alignment_group()=="o-A",
    !any(grepl("cohort two",alignment()$metadata$organism,fixed=TRUE)))
})

shiny::testServer(uniteReferenceServer,args=list(id="empty",dataset=function() NULL,initial_mode="primers"),{
  session$flushReact()
  stopifnot(!available(),mode()=="primers",output$available=="false",
    grepl("not available yet",output$availability$html,fixed=TRUE))
})
shiny::testServer(uniteReferenceServer,args=list(id="bad",dataset=function() stop("Checksum mismatch")),{
  session$flushReact()
  stopifnot(!available(),output$available=="false",grepl("could not be verified",output$availability$html,fixed=TRUE))
})
bad <- list(alignment=c(A="ACG",B="AC"),metadata=data.frame(record_id=c("A","B")))
stopifnot(inherits(try(.unite_reference_ui_alignment_window(bad),silent=TRUE),"try-error"))
unlink(scratch,recursive=TRUE)
message("UNITE UI modes, honest breadth/sparsity, true regional-window bases, sample denominators, missing/ambiguous evidence, exports and unavailable-only states passed.")
