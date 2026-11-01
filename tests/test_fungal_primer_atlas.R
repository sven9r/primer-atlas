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
source("R/fungal_primer_orientation.R",local=atlas)
variants <- atlas$fungal_orientation_catalog(snapshot)
parent <- snapshot$catalog[snapshot$catalog$primer_name=="ITS1F" & snapshot$catalog$direction=="forward",]
reversed <- variants[variants$parent_catalog_id==parent$catalog_id,]
context <- atlas$fungal_primer_atlas_data(snapshot,"PZ276773.1","ALL")
reversed_map <- atlas$fungal_primer_atlas(context,selected_reverse=reversed$catalog_id,variant_catalog=variants)
reversed_groups <- filter_attr(has_class(collect(reversed_map,"g"),"fungal-atlas-primer"),"data-catalog-id",parent$catalog_id)
site <- context$rows[context$rows$catalog_id==parent$catalog_id,]
stopifnot(length(reversed_groups)>0L,all(vapply(reversed_groups,function(row)
  row$attribs[["data-direction"]]=="reverse" && row$attribs[["data-catalog-direction"]]=="forward" &&
  row$attribs[["data-assigned-catalog-id"]]==reversed$catalog_id && row$attribs[["data-selection-state"]]=="reversed" &&
  row$attribs[["data-three-prime-bp"]]==site$start && row$attribs[["data-start-bp"]]==site$start &&
  row$attribs[["data-end-bp"]]==site$end && row$attribs[["aria-pressed"]]=="true",logical(1L))),
  grepl("ITS1F [RC]",as.character(reversed_map),fixed=TRUE))
forged <- variants;forged$sequence[forged$catalog_id==reversed$catalog_id] <- parent$sequence
stopifnot(inherits(try(atlas$fungal_primer_atlas(context,selected_reverse=reversed$catalog_id,variant_catalog=forged),silent=TRUE),"try-error"),
  inherits(try(atlas$fungal_primer_atlas(context,selected_forward=parent$catalog_id,selected_reverse=reversed$catalog_id,variant_catalog=variants),silent=TRUE),"try-error"))
old<-atlas$fungal_primer_atlas_data(historical,"FN812768.2")
stopifnot(all(old$landmarks$accession=="FN812768.2"),
  all(old$landmarks$boundary_status[old$landmarks$landmark_type=="conserved_anchor"]=="partial"),
  nchar(as.character(atlas$fungal_primer_atlas(old)))>0L,
  identical(serialize(snapshot,NULL),original),
  inherits(try(atlas$fungal_primer_atlas_data(snapshot,"FOREIGN"),silent=TRUE),"try-error"))

# Responsive screen layouts are rebuilt in physical pixels rather than
# shrinking a 1000-pixel image and its text, or cropping its reverse arrows.
# The export/default layout continues to use the original 65..965 scale.
legacy <- atlas$fungal_primer_atlas(context)
for(svg in collect(legacy,"svg")) {
  stopifnot(as.numeric(strsplit(svg$attribs$viewBox," ",fixed=TRUE)[[1L]])[3L]==1000,
    is.null(svg$attribs$width),is.null(svg$attribs$style),
    svg$attribs[["data-axis-left"]]==65,svg$attribs[["data-axis-right"]]==965)
}
points_matrix <- function(points) {
  values <- as.numeric(strsplit(gsub(","," ",points,fixed=TRUE),"[[:space:]]+")[[1L]])
  matrix(values,ncol=2,byrow=TRUE)
}
responsive_contexts <- list(
  all=context,
  its=atlas$fungal_primer_atlas_data(snapshot,"PZ962702.1","ITS"),
  historical=old)
found_packed_gene_band <- FALSE
for(viewport in c(600,390,320)) for(context_name in names(responsive_contexts)) {
  value <- responsive_contexts[[context_name]]
  rendered <- atlas$fungal_primer_atlas(value,viewport_width=viewport)
  stopifnot("is-responsive" %in% strsplit(rendered$attribs$class," ",fixed=TRUE)[[1L]],
    !length(has_class(collect(rendered,"p"),"fungal-atlas-scroll-cue")))
  svgs <- collect(rendered,"svg")
  stopifnot(length(svgs)==nrow(value$panels))
  for(i in seq_along(svgs)) {
    svg <- svgs[[i]];panel <- value$panels[i,,drop=FALSE]
    box <- as.numeric(strsplit(svg$attribs$viewBox," ",fixed=TRUE)[[1L]])
    stopifnot(box[3L]==viewport,svg$attribs$width==viewport,svg$attribs$height==box[4L],
      grepl(paste0("width:",viewport,"px"),svg$attribs$style,fixed=TRUE),
      grepl("min-width:0",svg$attribs$style,fixed=TRUE),
      !grepl("transform",svg$attribs$style,fixed=TRUE),
      svg$attribs[["data-domain-start"]]==panel$start,
      svg$attribs[["data-domain-end"]]==panel$end,
      svg$attribs[["data-axis-left"]]>0,svg$attribs[["data-axis-right"]]<viewport)
    left <- as.numeric(svg$attribs[["data-axis-left"]]);right <- as.numeric(svg$attribs[["data-axis-right"]])
    raw_to_pixel <- function(bp) left+(pmax(panel$start,pmin(panel$end,bp))-panel$start)/(panel$end-panel$start)*(right-left)
    baseline <- filter_attr(collect(svg,"line"),"stroke","#cbd5e1")
    stopifnot(length(baseline)==1L,baseline[[1L]]$attribs$x1==raw_to_pixel(panel$start),
      baseline[[1L]]$attribs$x2==raw_to_pixel(panel$end))
    groups <- has_class(collect(svg,"g"),"fungal-atlas-primer")
    expected_sites <- value$rows[value$rows$end>=panel$start & value$rows$start<=panel$end,,drop=FALSE]
    stopifnot(setequal(vapply(groups,function(row) row$attribs[["data-catalog-id"]],character(1L)),expected_sites$catalog_id))
    for(group in groups) {
      site <- expected_sites[expected_sites$catalog_id==group$attribs[["data-catalog-id"]],,drop=FALSE]
      label <- collect(group,"text")
      binding <- Filter(function(line) line$attribs[["stroke-width"]] %in% c(4,7),collect(group,"line"))
      stopifnot(length(label)==1L,label[[1L]]$attribs[["font-size"]]==14,
        group$attribs[["data-label-left"]]>=left,group$attribs[["data-label-right"]]<=right,
        group$attribs[["data-start-bp"]]==site$start,group$attribs[["data-end-bp"]]==site$end,
        group$attribs[["data-selection-state"]]=="unselected",
        group$attribs[["data-assigned-catalog-id"]]==site$catalog_id,
        group$attribs[["data-catalog-direction"]]==site$direction,
        length(binding)==1L,abs(binding[[1L]]$attribs$x1-raw_to_pixel(site$start))<1e-8,
        abs(binding[[1L]]$attribs$x2-raw_to_pixel(site$end))<1e-8)
      for(arrow in collect(group,"polygon")) {
        points <- points_matrix(arrow$attribs$points)
        stopifnot(all(points[,1L]>=0 & points[,1L]<=viewport),
          abs(points[1L,1L]-raw_to_pixel(group$attribs[["data-three-prime-bp"]]))<1e-8)
      }
    }
    # Every role is repacked for the actual available width; same-lane labels
    # cannot overlap even where many primers share a physical binding site.
    lanes <- split(groups,vapply(groups,function(row) paste(row$attribs[["data-direction"]],row$attribs[["data-label-lane"]]),character(1L)))
    for(lane in lanes) {
      bounds <- do.call(rbind,lapply(lane,function(row) c(row$attribs[["data-label-left"]],row$attribs[["data-label-right"]])))
      bounds <- bounds[order(bounds[,1L]),,drop=FALSE]
      if(nrow(bounds)>1L) stopifnot(all(bounds[-1L,1L]>head(bounds[,2L],-1L)))
    }
    gene_groups <- has_class(collect(svg,"g"),"fungal-atlas-landmark")
    packed_genes <- filter_attr(gene_groups,"data-label-packed","true")
    for(gene in gene_groups) stopifnot(collect(gene,"text")[[1L]]$attribs[["font-size"]]==14)
    for(gene in packed_genes) stopifnot(gene$attribs[["data-label-left"]]>=left,gene$attribs[["data-label-right"]]<=right)
    forward_groups <- filter_attr(groups,"data-direction","forward")
    if(length(packed_genes) && length(forward_groups)) {
      # Full 14px text bounds have a gap from the closest primer's arrow/body;
      # checking only x packing would miss the narrow ITS/5.8S overlap.
      gene_top <- min(vapply(packed_genes,function(gene) collect(gene,"text")[[1L]]$attribs$y-14,numeric(1L)))
      forward_bottom <- max(vapply(forward_groups,function(group) {
        binding <- Filter(function(line) line$attribs[["stroke-width"]] %in% c(4,7),collect(group,"line"))[[1L]]
        binding$attribs$y1+6
      },numeric(1L)))
      stopifnot(gene_top-forward_bottom>=12)
      if(any(vapply(packed_genes,function(gene) gene$attribs[["data-label-lane"]]>1L,logical(1L)))) found_packed_gene_band <- TRUE
    }
    for(domain in has_class(collect(svg,"g"),"fungal-atlas-domain")) {
      stopifnot(domain$attribs[["data-label-left"]]>=left,domain$attribs[["data-label-right"]]<=right)
    }
    # Tick labels and every rendered annotation remain inside the actual SVG.
    for(text in collect(svg,"text")) stopifnot(text$attribs$x>=0,text$attribs$x<=viewport,
      text$attribs$y>0,text$attribs$y<box[4L])
  }
}
stopifnot(found_packed_gene_band)

# At narrow screen widths, a reversed primer retains the exact parent click
# identity and physical interval while moving to the opposite-role lane.
reversed_site <- context$rows[context$rows$catalog_id==parent$catalog_id,,drop=FALSE]
for(viewport in c(600,390)) {
  rendered <- atlas$fungal_primer_atlas(context,selected_reverse=reversed$catalog_id,
    variant_catalog=variants,viewport_width=viewport)
  groups <- filter_attr(has_class(collect(rendered,"g"),"fungal-atlas-primer"),"data-catalog-id",parent$catalog_id)
  stopifnot(length(groups)>0L,all(vapply(groups,function(group)
    group$attribs[["data-assigned-catalog-id"]]==reversed$catalog_id &&
    group$attribs[["data-selection-state"]]=="reversed" &&
    group$attribs[["data-catalog-direction"]]=="forward" && group$attribs[["data-direction"]]=="reverse" &&
    group$attribs[["data-start-bp"]]==reversed_site$start && group$attribs[["data-end-bp"]]==reversed_site$end &&
    group$attribs[["data-three-prime-bp"]]==reversed_site$start &&
    collect(group,"text")[[1L]]$attribs[["font-size"]]==14,logical(1L))))
}
for(width in list(319,0,Inf,NA_real_,numeric(),c(390,600),"600")) {
  stopifnot(inherits(try(atlas$fungal_primer_atlas(context,viewport_width=width),silent=TRUE),"try-error"))
}
stopifnot(identical(serialize(snapshot,NULL),original))
cat("Individual fungal atlas raw sites, memberships, roles, copies, missing states and responsive pixel-width label/arrow packing passed.\n")
