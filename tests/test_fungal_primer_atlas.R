#!/usr/bin/env Rscript
source("R/functions.R")
source("R/data_layer.R")
source("R/fungal_explorer.R")
source("R/fungal_versions.R")
# The new atlas is self-contained: even its raw IUPAC validation does not rely
# on mutable application globals or the selected/preset combination list.
atlas <- new.env(parent=baseenv())
source("R/fungal_primer_atlas.R",local=atlas)
collect <- function(x,name=NULL) {
  if(inherits(x,"shiny.tag")) return(c(if(is.null(name) || x$name==name) list(x) else list(),
    unlist(lapply(x$children,collect,name=name),recursive=FALSE)))
  if(is.list(x)) return(unlist(lapply(x,collect,name=name),recursive=FALSE))
  list()
}
filter_attr <- function(rows,name,value) Filter(function(row) identical(as.character(row$attribs[[name]]),as.character(value)),rows)
has_class <- function(rows,class) Filter(function(row) !is.null(row$attribs$class) && class %in% strsplit(row$attribs$class," ",fixed=TRUE)[[1L]],rows)

snapshot <- fungal_read_snapshot("v3-c6fa87f256e1",remote=FALSE)
original <- serialize(snapshot,NULL)
expected <- list(PZ962702.1=c(60L,29L,31L),PZ276773.1=c(62L,29L,33L),FN812768.2=c(15L,10L,5L))
for(accession in names(expected)) {
  value <- atlas$fungal_primer_atlas_data(snapshot,accession,"ALL")
  frozen <- as.data.frame(snapshot$sites[snapshot$sites$accession==accession & snapshot$sites$status=="exact_compatible_site",])
  stopifnot(value$n_eligible==112L,nrow(value$rows)==expected[[accession]][1L],
    sum(value$rows$direction=="forward")==expected[[accession]][2L],
    sum(value$rows$direction=="reverse")==expected[[accession]][3L],
    nrow(value$rows)+nrow(value$missing)==112L,
    setequal(value$rows$catalog_id,frozen$catalog_id),
    all(value$rows$start==frozen$start[match(value$rows$catalog_id,frozen$catalog_id)]),
    all(value$rows$end==frozen$end[match(value$rows$catalog_id,frozen$catalog_id)]),
    all(value$rows$end-value$rows$start+1L==nchar(value$rows$sequence)),
    all(value$rows$accession==accession),all(value$rows$version_id==snapshot$manifest$id))
  selected_f <- value$rows$catalog_id[value$rows$direction=="forward"][1L]
  selected_r <- value$rows$catalog_id[value$rows$direction=="reverse"][1L]
  rendered <- atlas$fungal_primer_atlas(value,selected_forward=selected_f,selected_reverse=selected_r)
  groups <- has_class(collect(rendered,"g"),"fungal-atlas-primer")
  stopifnot(setequal(vapply(groups,function(row) row$attribs[["data-catalog-id"]],character(1L)),value$rows$catalog_id),
    all(vapply(groups,function(row) row$attribs[["data-accession"]]==accession,logical(1L))),
    all(vapply(groups,function(row) row$attribs[["data-version-id"]]==snapshot$manifest$id,logical(1L))),
    all(vapply(groups,function(row) row$attribs[["data-primer-input"]]=="fungal_atlas_primer_click",logical(1L))))
  for(row in groups) {
    site <- value$rows[match(row$attribs[["data-catalog-id"]],value$rows$catalog_id),]
    stopifnot(row$attribs[["data-start-bp"]]==site$start,row$attribs[["data-end-bp"]]==site$end,
      row$attribs[["data-three-prime-bp"]]==if(site$direction=="forward") site$end else site$start,
      row$attribs[["data-direction"]]==site$direction,row$attribs$role=="button",
      row$attribs[["aria-pressed"]]==if(site$catalog_id %in% c(selected_f,selected_r)) "true" else "false")
  }
  # Labels in one role/panel/lane never overlap; source coordinates stay at the
  # arrow endpoint while collision-packed labels are independently positioned.
  packing_groups <- split(groups,vapply(groups,function(row) paste(row$attribs[["data-panel-id"]],
    row$attribs[["data-direction"]],row$attribs[["data-label-lane"]]),character(1L)))
  for(group in packing_groups) {
    starts <- vapply(group,function(row) as.numeric(row$attribs[["data-label-left"]]),numeric(1L))
    ends <- vapply(group,function(row) as.numeric(row$attribs[["data-label-right"]]),numeric(1L))
    order <- order(starts)
    if(length(order)>1L) stopifnot(all(starts[order][-1L]>head(ends[order],-1L)))
  }
  text <- paste(as.character(rendered),collapse="\n")
  stopifnot(grepl("positions on NCBI",text,fixed=TRUE),grepl("UNITE catalog",text,fixed=TRUE),
    !grepl("own_template_",text,fixed=TRUE),!grepl("V1|V2|V3",paste(value$landmarks$label,collapse=" ")))
}

# Do not pool PZ276's early and late partial LSU copies into an empty 10-kb
# ITS panel. Actual offsets and copy boundaries are preserved separately.
pz <- atlas$fungal_primer_atlas_data(snapshot,"PZ276773.1","ALL")
stopifnot(nrow(pz$panels)==2L,pz$panels$start[pz$panels$panel_id=="SSU"]==5767L,
  pz$panels$start[pz$panels$panel_id=="ITS_LSU"]==7270L,
  pz$panels$end[pz$panels$panel_id=="ITS_LSU"]==10646L,
  all(pz$panels$end-pz$panels$start<4000L))
lsu <- atlas$fungal_primer_atlas_data(snapshot,"PZ276773.1","LSU")
stopifnot(nrow(lsu$panels)==2L,identical(lsu$panels$start,c(1,8127)),identical(lsu$panels$end,c(1182,10646)))
domain <- atlas$fungal_primer_atlas_data(snapshot,"PZ962702.1","landmark:28S_D1")
stopifnot(domain$panels$start==1543L,domain$panels$end==4690L,
  all(domain$landmarks$evidence_status[domain$landmarks$is_domain]=="mapped_reference_domain"),
  !any(domain$landmarks$label %in% c("D6","D7b","D12")))
missing_region <- atlas$fungal_primer_atlas_data(snapshot,"FN812768.2","landmark:28S_D1")
stopifnot(nrow(missing_region$panels)==0L,nrow(missing_region$rows)==0L,nrow(missing_region$outside_region)==15L)

# Different oligos sharing one exact interval remain separately selectable.
same <- pz$rows[pz$rows$primer_name %in% c("fITS7","gITS7"),]
stopifnot(nrow(same)==2L,length(unique(same$catalog_id))==2L,length(unique(same$start))==1L,length(unique(same$end))==1L)
without_pairs <- snapshot
without_pairs$pairs <- data.frame()
without_pairs$pair_sequences <- data.frame()
stopifnot(identical(atlas$fungal_primer_atlas_data(without_pairs,"PZ962702.1"),
  atlas$fungal_primer_atlas_data(snapshot,"PZ962702.1")))

# Broken raw ownership/coordinates and multi-hit results are withheld from
# glyphs. No reported_position or another accession can repair missing sites.
bad <- snapshot
target <- which(bad$sites$accession=="PZ962702.1" & bad$sites$primer_name=="ITS1")
bad$sites$start[target] <- bad$sites$start[target]+1L
invalid <- atlas$fungal_primer_atlas_data(bad,"PZ962702.1")
stopifnot(!"ITS1" %in% invalid$rows$primer_name,
  invalid$missing$site_status[invalid$missing$primer_name=="ITS1"]=="invalid_raw_site_evidence")
multiple <- snapshot
multiple$sites$n_matches[target] <- 3L
ambiguous <- atlas$fungal_primer_atlas_data(multiple,"PZ962702.1")
stopifnot(!"ITS1" %in% ambiguous$rows$primer_name,"ITS1" %in% ambiguous$multiple$primer_name,
  grepl("3 compatible matches",as.character(atlas$fungal_primer_atlas(ambiguous)),fixed=TRUE))
absent <- snapshot
absent$sites <- absent$sites[-target,]
unavailable <- atlas$fungal_primer_atlas_data(absent,"PZ962702.1")
stopifnot(!"ITS1" %in% unavailable$rows$primer_name,
  unavailable$missing$site_status[unavailable$missing$primer_name=="ITS1"]=="binding_evidence_unavailable")

# Synthetic tightly clustered domains exercise every packing lane, rather than
# checking only the first lane and letting the second lane collide.
clustered <- atlas$fungal_primer_atlas_data(snapshot,"PZ962702.1")
extra <- clustered$landmarks[clustered$landmarks$is_domain,,drop=FALSE][rep(1L,6L),]
extra$landmark_id<-paste0("fixture_domain_",1:6);extra$label<-paste0("D",1:6)
extra$start<-2500;extra$end<-2505
clustered$landmarks<-rbind(clustered$landmarks[!clustered$landmarks$is_domain,],extra)
rendered<-atlas$fungal_primer_atlas(clustered)
domains<-has_class(collect(rendered,"g"),"fungal-atlas-domain")
stopifnot(length(unique(vapply(domains,function(row) as.character(row$attribs[["data-label-lane"]]),character(1L))))==6L)
for(svg in collect(rendered,"svg")) {
  height<-as.numeric(strsplit(svg$attribs$viewBox," ",fixed=TRUE)[[1L]])[4L]
  texts<-collect(svg,"text")
  stopifnot(all(vapply(texts,function(row) as.numeric(row$attribs$y)>0 && as.numeric(row$attribs$y)<height,logical(1L))))
}
historical<-fungal_read_snapshot("v1-9e958f105c46",remote=FALSE)
old<-atlas$fungal_primer_atlas_data(historical,"FN812768.2")
stopifnot(all(old$landmarks$accession=="FN812768.2"),
  all(old$landmarks$boundary_status[old$landmarks$landmark_type=="conserved_anchor"]=="partial"),
  nchar(as.character(atlas$fungal_primer_atlas(old)))>0L,
  identical(serialize(snapshot,NULL),original),
  inherits(try(atlas$fungal_primer_atlas_data(snapshot,"FOREIGN"),silent=TRUE),"try-error"))
cat("Individual fungal atlas raw sites, all-catalog membership, role/selection metadata, split gene copies, honest missing states and readable label packing passed.\n")
