# Individual frozen primer placements, independent of the preset pair list.
# Every plotted endpoint belongs to the selected raw template, never an MSA.
.fungal_atlas_iupac <- list(A="A",C="C",G="G",T="T",U="T",R=c("A","G"),Y=c("C","T"),
  S=c("C","G"),W=c("A","T"),K=c("G","T"),M=c("A","C"),B=c("C","G","T"),
  D=c("A","G","T"),H=c("A","C","T"),V=c("A","C","G"),N=c("A","C","G","T"),I=c("A","C","G","T"))
.fungal_atlas_rc <- function(sequence) paste(rev(strsplit(chartr("ACGTURYSWKMBDHVNI",
  "TGCAAYRSWMKVHDBNI",toupper(sequence)),"",fixed=TRUE)[[1L]]),collapse="")
.fungal_atlas_text <- function(row,field,default="") {
  if(!field %in% names(row) || !length(row[[field]]) || is.na(row[[field]][1L])) return(default)
  as.character(row[[field]][1L])
}
.fungal_atlas_empty_landmarks <- function() data.frame(landmark_id=character(),label=character(),start=numeric(),end=numeric(),
  landmark_type=character(),boundary_status=character(),evidence_status=character(),source_note=character(),accession=character(),is_domain=logical())

.fungal_atlas_clusters <- function(rows) {
  if(!nrow(rows)) return(list())
  rows <- rows[order(rows$start,rows$end),,drop=FALSE]
  group <- integer(nrow(rows));current<-0L;last_end<- -Inf
  for(i in seq_len(nrow(rows))) {
    if(rows$start[i]>last_end+301) current<-current+1L
    group[i]<-current;last_end<-max(last_end,rows$end[i])
  }
  split(rows,group)
}

.fungal_atlas_landmarks <- function(snapshot,accession,size) {
  rows <- snapshot$landmarks
  if(is.null(rows) || !nrow(rows)) return(.fungal_atlas_empty_landmarks())
  rows <- as.data.frame(rows,stringsAsFactors=FALSE)
  if(!all(c("landmark_id","label","start","end") %in% names(rows))) stop("Atlas landmark identity fields are incomplete")
  if(!"accession" %in% names(rows)) {
    owner <- snapshot$manifest$reference_accession
    if(length(owner)!=1L || is.na(owner) || !nzchar(owner)) owner <- NA_character_
    rows$accession <- rep(owner,nrow(rows))
  }
  rows <- rows[!is.na(rows$accession) & rows$accession==accession,,drop=FALSE]
  for(field in c("landmark_type","boundary_status","evidence_status","source_note")) {
    if(!field %in% names(rows)) rows[[field]] <- rep("",nrow(rows))
    rows[[field]][is.na(rows[[field]])] <- ""
  }
  rows$start <- suppressWarnings(as.numeric(rows$start));rows$end <- suppressWarnings(as.numeric(rows$end))
  valid <- is.finite(rows$start) & is.finite(rows$end) & rows$start==floor(rows$start) & rows$end==floor(rows$end) &
    rows$start>=1L & rows$end>=rows$start & rows$end<=size
  withheld <- rows$boundary_status %in% c("unresolved","not_detected","invalid_coordinates","withheld","extraction_mismatch")
  rows <- rows[valid & !withheld,,drop=FALSE]
  if(anyDuplicated(rows$landmark_id)) stop("Own-template atlas landmark identifiers repeat")
  domain <- rows$landmark_type %in% c("domain","divergent_domain") | rows$evidence_status=="mapped_reference_domain"
  # Unknown structural labels never become inferred gene architecture.
  visible <- domain | rows$label %in% c("18S","SSU","ITS1","5.8S","ITS2","28S","LSU")
  rows <- rows[visible,,drop=FALSE];domain<-domain[visible]
  rows$landmark_type[!nzchar(rows$landmark_type)] <- ifelse(domain[!nzchar(rows$landmark_type)],"divergent_domain",
    ifelse(rows$label[!nzchar(rows$landmark_type)] %in% c("ITS1","ITS2"),"spacer","gene"))
  empty <- !nzchar(rows$boundary_status)
  rows$boundary_status[empty & (rows$landmark_type=="conserved_anchor" | grepl("partial|conserved.anchor",rows$source_note,ignore.case=TRUE))] <- "partial"
  rows$boundary_status[!nzchar(rows$boundary_status)] <- "source_interval"
  rows$is_domain <- domain
  rows$raw_bp <- rows$end-rows$start+1L
  rows[order(rows$start,rows$end,rows$landmark_id),,drop=FALSE]
}

fungal_primer_atlas_data <- function(snapshot,accession,region="ALL") {
  if(!is.list(snapshot) || !all(c("manifest","catalog","sequences","sites") %in% names(snapshot))) stop("A frozen catalog, raw templates and individual site table are required")
  if(length(accession)!=1L || is.na(accession) || !accession %in% names(snapshot$sequences)) stop("Choose an accession in this frozen snapshot")
  sequence <- toupper(unname(snapshot$sequences[[accession]]))
  if(length(sequence)!=1L || is.na(sequence) || !nzchar(sequence) || !grepl("^[ACGTURYSWKMBDHVNI]+$",sequence)) stop("The atlas reference must be an ungapped raw nucleotide template")
  size <- nchar(sequence)
  catalog <- as.data.frame(snapshot$catalog,stringsAsFactors=FALSE)
  required <- c("catalog_id","primer_name","direction","sequence","gene_locus")
  if(!all(required %in% names(catalog)) || anyNA(catalog$catalog_id) || anyDuplicated(catalog$catalog_id)) stop("Atlas catalog identities are incomplete or repeated")
  eligible <- !is.na(catalog$direction) & catalog$direction %in% c("forward","reverse") &
    !is.na(catalog$gene_locus) & catalog$gene_locus %in% c("SSU","SSU/ITS","ITS","ITS/LSU","LSU")
  catalog <- catalog[eligible,,drop=FALSE]
  cache <- as.data.frame(snapshot$sites,stringsAsFactors=FALSE)
  if(!all(c("catalog_id","accession","start","end","n_matches","status") %in% names(cache))) stop("Atlas individual-site evidence is incomplete")
  cache <- cache[!is.na(cache$accession) & cache$accession==accession,,drop=FALSE]
  rows <- catalog
  rows$start <- rows$end <- rep(NA_integer_,nrow(rows))
  rows$n_matches <- rep(NA_integer_,nrow(rows))
  rows$site_status <- rep("binding_evidence_unavailable",nrow(rows))
  rows$placement_reason <- rep("No individual site result is frozen for this component on this template",nrow(rows))
  for(i in seq_len(nrow(rows))) {
    matches <- which(!is.na(cache$catalog_id) & cache$catalog_id==rows$catalog_id[i])
    if(!length(matches)) next
    if(length(matches)!=1L) {
      rows$site_status[i] <- "duplicate_frozen_site_rows"
      rows$placement_reason[i] <- "Repeated frozen site identities prevent a unique placement"
      next
    }
    site <- cache[matches,,drop=FALSE]
    rows$site_status[i] <- .fungal_atlas_text(site,"status","binding_evidence_unavailable")
    rows$n_matches[i] <- suppressWarnings(as.integer(site$n_matches))
    rows$placement_reason[i] <- gsub("_"," ",rows$site_status[i],fixed=TRUE)
    identity_ok <- (!"primer_name" %in% names(site) || identical(as.character(site$primer_name),as.character(rows$primer_name[i]))) &&
      (!"direction" %in% names(site) || identical(as.character(site$direction),as.character(rows$direction[i])))
    if(!identity_ok) {
      rows$site_status[i] <- "frozen_component_identity_mismatch"
      rows$placement_reason[i] <- "Frozen site role/name differs from this exact catalog component"
      next
    }
    if(rows$site_status[i]!="exact_compatible_site") next
    if(is.na(rows$n_matches[i]) || rows$n_matches[i]!=1L) {
      rows$site_status[i] <- "ambiguous_binding_site"
      rows$placement_reason[i] <- "The frozen result does not establish exactly one compatible site"
      next
    }
    bounds <- suppressWarnings(as.numeric(c(site$start,site$end)))
    oligo <- toupper(as.character(rows$sequence[i]))
    valid <- length(bounds)==2L && all(is.finite(bounds)) && all(bounds==floor(bounds)) && bounds[1L]>=1L &&
      bounds[2L]>=bounds[1L] && bounds[2L]<=size && !is.na(oligo) && nzchar(oligo) && bounds[2L]-bounds[1L]+1L==nchar(oligo)
    if(valid) {
      query <- strsplit(if(rows$direction[i]=="reverse") .fungal_atlas_rc(oligo) else oligo,"",fixed=TRUE)[[1L]]
      observed <- strsplit(substr(sequence,bounds[1L],bounds[2L]),"",fixed=TRUE)[[1L]]
      valid <- all(vapply(seq_along(query),function(j) observed[j] %in% .fungal_atlas_iupac[[query[j]]],logical(1)))
    }
    if(!valid) {
      rows$site_status[i] <- "invalid_raw_site_evidence"
      rows$placement_reason[i] <- "Frozen site bounds or observed raw bases do not match the exact oligo"
      next
    }
    rows$start[i] <- as.integer(bounds[1L]);rows$end[i] <- as.integer(bounds[2L])
    rows$placement_reason[i] <- "Unique exact IUPAC-compatible site on this raw template"
  }
  rows$accession <- rep(accession,nrow(rows))
  rows$version_id <- rep(as.character(snapshot$manifest$id),nrow(rows))
  rows$three_prime_bp <- ifelse(rows$direction=="forward",rows$end,rows$start)
  verified <- rows[rows$site_status=="exact_compatible_site" & !is.na(rows$start),,drop=FALSE]
  missing <- rows[!rows$catalog_id %in% verified$catalog_id,,drop=FALSE]
  landmarks <- .fungal_atlas_landmarks(snapshot,accession,size)
  genes <- landmarks[!landmarks$is_domain,,drop=FALSE]
  if(length(region)!=1L || is.na(region)) stop("Choose one atlas region")
  original_region <- region
  region <- switch(toupper(region),"18S"="SSU","28S"="LSU","5_8S"="5.8S",toupper(region))
  panels <- list()
  add_panel <- function(id,title,start,end,basis) {
    bounds <- pmax(1,pmin(size,c(start,end)))
    if(bounds[1L]>=bounds[2L] && size>1L) bounds <- c(max(1,bounds[1L]-1),min(size,bounds[2L]+1))
    panels[[length(panels)+1L]] <<- data.frame(panel_id=id,title=title,start=bounds[1L],end=bounds[2L],basis=basis,stringsAsFactors=FALSE)
  }
  ssu <- genes[genes$label %in% c("18S","SSU"),,drop=FALSE]
  downstream <- genes[genes$label %in% c("ITS1","5.8S","ITS2","28S","LSU"),,drop=FALSE]
  downstream_groups <- .fungal_atlas_clusters(downstream)
  its_groups <- Filter(function(x) any(x$label %in% c("ITS1","5.8S","ITS2")),downstream_groups)
  focused_downstream <- if(length(its_groups)) do.call(rbind,its_groups) else downstream
  if(startsWith(original_region,"landmark:")) {
    own <- landmarks[landmarks$landmark_id==substring(original_region,10L),,drop=FALSE]
    if(nrow(own)==1L) add_panel("focus",paste0(own$label," · mapped region"),own$start-1000,own$end+2000,"own_template_landmark_and_raw_context")
  } else if(region=="ALL" && nrow(ssu) && nrow(downstream)) {
    groups <- if(length(its_groups)) its_groups else downstream_groups
    for(j in seq_along(groups)) {
      own<-groups[[j]]
      add_panel(paste0("ITS_LSU",if(length(groups)>1L) paste0("_",j) else ""),
        if(any(own$label %in% c("ITS1","5.8S","ITS2"))) "ITS / 5.8S / LSU" else "28S / LSU",
        min(own$start)-300,max(own$end)+100,"own_template_regions_with_adjoining_SSU_tail")
    }
    add_panel("SSU","18S / SSU",min(ssu$start),max(ssu$end),"own_template_gene_intervals")
    # Any verified site beyond the annotation-focused panels gets its own
    # bounded raw-site panel; it is never silently dropped or assigned a gene.
    covered <- vapply(seq_len(nrow(verified)),function(i) any(vapply(panels,function(panel)
      verified$end[i]>=panel$start && verified$start[i]<=panel$end,logical(1))),logical(1))
    outside <- verified[!covered,,drop=FALSE]
    if(nrow(outside)) add_panel("OTHER","Additional verified raw sites",min(outside$start)-100,max(outside$end)+100,"verified_raw_sites_gene_bounds_unavailable")
  } else if(region %in% c("ALL","SSU","ITS","ITS1","5.8S","ITS2","LSU")) {
    labels <- switch(region,SSU=c("18S","SSU"),LSU=c("28S","LSU"),ITS=c("ITS1","5.8S","ITS2"),ALL=unique(genes$label),region)
    own <- genes[genes$label %in% labels,,drop=FALSE]
    if(nrow(own)) {
      if(region %in% c("ALL","ITS") && nrow(focused_downstream)) {
        groups <- .fungal_atlas_clusters(focused_downstream)
        for(j in seq_along(groups)) {
          row<-groups[[j]]
          add_panel(paste0("ITS_LSU",if(length(groups)>1L) paste0("_",j) else ""),"ITS / 5.8S / LSU",
            min(row$start)-300,max(row$end)+100,"own_template_region_intervals")
        }
      } else if(region=="LSU") {
        copies <- .fungal_atlas_clusters(own)
        for(j in seq_along(copies)) add_panel(paste0("LSU_",j),paste0("28S / LSU",if(length(copies)>1L) paste0(" · interval ",j) else ""),
          min(copies[[j]]$start),max(copies[[j]]$end),"own_template_region_intervals")
      }
      else add_panel(region,if(region=="SSU") "18S / SSU" else if(region=="LSU") "28S / LSU" else region,
        min(own$start),max(own$end),"own_template_region_intervals")
    } else {
      locus <- switch(region,SSU="SSU|18S",LSU="LSU|28S",ITS="ITS|5[.]8S",ITS1="ITS",ITS2="ITS",`5.8S`="ITS",ALL=".*")
      relevant <- verified[grepl(locus,verified$gene_locus),,drop=FALSE]
      if(nrow(relevant)) add_panel("SITES",paste0(if(region=="ALL") "Verified primer sites" else paste(region,"associated primer sites")," · gene bounds unavailable"),
        min(relevant$start)-100,max(relevant$end)+100,"reported_component_locus_without_gene_boundaries")
    }
  } else stop("Unknown atlas region")
  if(!length(panels)) panels <- data.frame(panel_id=character(),title=character(),start=numeric(),end=numeric(),basis=character()) else panels <- do.call(rbind,panels)
  shown <- if(nrow(panels)) vapply(seq_len(nrow(verified)),function(i) any(verified$end[i]>=panels$start & verified$start[i]<=panels$end),logical(1)) else rep(FALSE,nrow(verified))
  list(rows=verified[shown,,drop=FALSE],landmarks=landmarks,missing=missing,
    multiple=missing[missing$site_status %in% c("ambiguous_binding_site","duplicate_frozen_site_rows"),,drop=FALSE],
    outside_region=verified[!shown,,drop=FALSE],panels=panels,catalog=catalog,
    accession=accession,version_id=as.character(snapshot$manifest$id),reference_length=size,region=original_region,
    n_verified=nrow(verified),n_eligible=nrow(catalog))
}

.fungal_atlas_pack <- function(rows,scale,left=65,right=965,font_size=14) {
  if(!nrow(rows)) return(data.frame(catalog_id=character(),label_x=numeric(),label_left=numeric(),label_right=numeric(),lane=integer()))
  centers <- scale((rows$start+rows$end)/2)
  widths <- pmax(40,nchar(rows$primer_name)*font_size*.59+12)
  widths <- pmin(widths,right-left)
  xmin <- pmax(left,pmin(right-widths,centers-widths/2));xmax<-xmin+widths
  ids <- if("catalog_id" %in% names(rows)) rows$catalog_id else if("landmark_id" %in% names(rows)) rows$landmark_id else as.character(seq_len(nrow(rows)))
  order <- order(xmin,xmax,rows$primer_name,ids)
  lanes <- integer(nrow(rows));occupied<-numeric()
  for(i in order) {
    fits <- which(occupied+14<xmin[i])
    lane <- if(length(fits)) fits[1L] else length(occupied)+1L
    occupied[lane]<-xmax[i];lanes[i]<-lane
  }
  data.frame(catalog_id=ids,label_x=(xmin+xmax)/2,label_left=xmin,label_right=xmax,lane=lanes)
}

fungal_primer_atlas <- function(data,selected_forward=character(),selected_reverse=character(),click_input="fungal_atlas_primer_click") {
  if(!is.list(data) || !all(c("rows","landmarks","panels","missing","catalog","accession","version_id") %in% names(data))) stop("Build an explicit frozen primer-atlas context first")
  if(any(!selected_forward %in% data$catalog$catalog_id[data$catalog$direction=="forward"]) ||
      any(!selected_reverse %in% data$catalog$catalog_id[data$catalog$direction=="reverse"])) stop("Atlas selections do not belong to this frozen catalog and direction")
  tags <- htmltools::tags
  if(!nrow(data$panels)) return(tags$div(class="fungal-primer-atlas",tags$p("No individually verified primer placements are available in this view. Missing or ambiguous sites do not establish locus absence.")))
  panel_tags <- lapply(seq_len(nrow(data$panels)),function(p) {
    panel <- data$panels[p,,drop=FALSE];bounds<-c(panel$start,panel$end);left<-65;right<-965;width<-1000
    scale <- function(bp) left+(pmax(bounds[1L],pmin(bounds[2L],bp))-bounds[1L])/diff(bounds)*(right-left)
    sites <- data$rows[data$rows$end>=bounds[1L] & data$rows$start<=bounds[2L],,drop=FALSE]
    forwards<-sites[sites$direction=="forward",,drop=FALSE];reverses<-sites[sites$direction=="reverse",,drop=FALSE]
    fp <- .fungal_atlas_pack(forwards,scale);rp <- .fungal_atlas_pack(reverses,scale)
    domains<-data$landmarks[data$landmarks$is_domain & data$landmarks$end>=bounds[1L] & data$landmarks$start<=bounds[2L],,drop=FALSE]
    domain_labels<-domains;domain_labels$primer_name<-domain_labels$label
    dp<-.fungal_atlas_pack(domain_labels,scale,font_size=13)
    domain_height<-if(nrow(domains)) 12+max(dp$lane)*18 else 0
    nf<-max(c(0L,fp$lane));nr<-max(c(0L,rp$lane));step<-42
    gene_y<-90+nf*step;height<-gene_y+110+domain_height+max(1,nr)*step
    children<-list(tags$text(x=left,y=22,fill="#475569",`font-size`=14,"Raw template bases (bp)"))
    ticks<-pretty(bounds,n=6);ticks<-ticks[ticks>=bounds[1L] & ticks<=bounds[2L]]
    for(tick in ticks) children<-append(children,list(tags$line(x1=scale(tick),x2=scale(tick),y1=32,y2=height-25,stroke="#e2e8f0",`stroke-width`=1),
      tags$text(x=scale(tick),y=height-8,`text-anchor`="middle",fill="#64748b",`font-size`=13,format(tick,big.mark=",",trim=TRUE))))
    children<-append(children,list(tags$line(x1=left,x2=right,y1=gene_y,y2=gene_y,stroke="#cbd5e1",`stroke-width`=3)))
    draw_primers <- function(rows,packed,role) {
      lapply(seq_len(nrow(rows)),function(i) {
        row<-rows[i,,drop=FALSE];layout<-packed[i,,drop=FALSE]
        selected<-row$catalog_id %in% if(role=="forward") selected_forward else selected_reverse
        color<-if(role=="forward") "#2563eb" else "#dc2626"
        y<-if(role=="forward") gene_y-28-layout$lane*step else gene_y+40+domain_height+layout$lane*step
        tip<-if(role=="forward") row$end else row$start
        visible_tip<-tip>=bounds[1L] && tip<=bounds[2L]
        source<-.fungal_atlas_text(row,"reference","Original citation unavailable")
        note<-paste0(row$primer_name," · ",role," · ",row$start,"–",row$end," raw bp on ",data$accession,
          "\nUnique exact-compatible individual site; no PCR or product support inferred.",
          "\nCatalog compilation: UNITE; original source: ",source,
          "\n",.fungal_atlas_text(row,"direction_note"),
          "\n",.fungal_atlas_text(row,"remarks"))
        label_y<-if(role=="forward") y-11 else y+23
        label_height<-23
        tags$g(class=paste("fungal-atlas-primer",if(selected) "is-selected" else ""),role="button",tabindex="0",
          `aria-label`=paste(if(selected) "Remove" else "Select",role,"primer",row$primer_name),`aria-pressed`=tolower(as.character(selected)),
          `data-catalog-id`=row$catalog_id,`data-direction`=role,`data-primer-name`=row$primer_name,
          `data-start-bp`=row$start,`data-end-bp`=row$end,`data-three-prime-bp`=tip,
          `data-primer-input`=click_input,`data-version-id`=data$version_id,`data-accession`=data$accession,
          `data-panel-id`=panel$panel_id,`data-selected`=tolower(as.character(selected)),
          `data-label-lane`=layout$lane,`data-label-left`=layout$label_left,`data-label-right`=layout$label_right,style="cursor:pointer",
          tags$title(note),
          tags$line(x1=scale((row$start+row$end)/2),x2=layout$label_x,y1=gene_y+if(role=="forward") -15 else 20,
            y2=y,stroke=color,opacity=.2,`stroke-width`=1),
          tags$rect(x=layout$label_left,y=label_y-17,width=layout$label_right-layout$label_left,height=label_height,rx=4,
            fill=if(selected) if(role=="forward") "#dbeafe" else "#fee2e2" else "transparent",
            stroke=if(selected) color else "none",`pointer-events`="all"),
          tags$line(x1=scale(row$start),x2=scale(row$end),y1=y,y2=y,stroke=color,`stroke-width`=if(selected) 7 else 4,opacity=if(selected) 1 else .75),
          if(visible_tip) tags$polygon(points=if(role=="forward") paste(scale(tip),y,scale(tip)-9,y-6,scale(tip)-9,y+6)
            else paste(scale(tip),y,scale(tip)+9,y-6,scale(tip)+9,y+6),fill=color),
          tags$text(x=layout$label_x,y=label_y,`text-anchor`="middle",fill=color,`font-size`=14,`font-weight`=if(selected) 800 else 600,row$primer_name))
      })
    }
    children<-append(children,draw_primers(forwards,fp,"forward"));children<-append(children,draw_primers(reverses,rp,"reverse"))
    genes<-data$landmarks[!data$landmarks$is_domain & data$landmarks$end>=bounds[1L] & data$landmarks$start<=bounds[2L],,drop=FALSE]
    gene_labels<-genes;gene_labels$primer_name<-genes$label
    gp<-.fungal_atlas_pack(gene_labels,scale)
    for(i in seq_len(nrow(genes))) {
      row<-genes[i,,drop=FALSE];partial<-grepl("partial|anchor",paste(row$boundary_status,row$landmark_type),ignore.case=TRUE)
      label<-row$label
      color<-if(row$label %in% c("ITS1","ITS2")) "#fef3c7" else "#dbe4ed"
      x1<-scale(row$start);x2<-scale(row$end)
      narrow<-x2-x1<nchar(label)*8+8
      label_x<-if(narrow) gp$label_x[i] else (x1+x2)/2
      label_y<-if(narrow) gene_y-18-gp$lane[i]*17 else gene_y+5
      children<-append(children,list(tags$g(class="fungal-atlas-landmark",`data-landmark-id`=row$landmark_id,`data-start-bp`=row$start,
        `data-end-bp`=row$end,`data-boundary-status`=row$boundary_status,
        tags$title(paste(row$label,row$start,"–",row$end,"raw bp ·",row$boundary_status,"·",row$evidence_status,"·",row$source_note)),
        tags$rect(x=x1,y=gene_y-14,width=max(1,x2-x1),height=28,rx=2,fill=color,stroke="#64748b",`stroke-dasharray`=if(partial) "4 3" else NULL),
        if(narrow) tags$line(x1=(x1+x2)/2,x2=label_x,y1=gene_y-14,y2=label_y+4,stroke="#64748b",`stroke-width`=1),
        tags$text(x=label_x,y=label_y,`text-anchor`="middle",fill="#334155",`font-size`=14,`font-weight`=700,label))))
    }
    if(nrow(domains)) {
      for(i in seq_len(nrow(domains))) {
        row<-domains[i,,drop=FALSE]
        children<-append(children,list(tags$g(class="fungal-atlas-domain",`data-landmark-id`=row$landmark_id,
          `data-start-bp`=row$start,`data-end-bp`=row$end,`data-evidence-status`=row$evidence_status,
          `data-label-lane`=dp$lane[i],`data-label-left`=dp$label_left[i],`data-label-right`=dp$label_right[i],
          tags$title(paste(row$label,row$start,"–",row$end,"raw bp ·",row$evidence_status,"·",row$source_note)),
          tags$rect(x=scale(row$start),y=gene_y+17,width=max(1,scale(row$end)-scale(row$start)),height=7,fill="#94a3b8"),
          tags$line(x1=scale((row$start+row$end)/2),x2=dp$label_x[i],y1=gene_y+24,y2=gene_y+27+dp$lane[i]*18,stroke="#94a3b8",`stroke-width`=1),
          tags$text(x=dp$label_x[i],y=gene_y+34+dp$lane[i]*18,`text-anchor`="middle",fill="#64748b",`font-size`=13,row$label))))
      }
    }
    tags$section(class="fungal-atlas-panel",`data-panel-id`=panel$panel_id,
      tags$h4(panel$title),tags$p(class="fungal-atlas-axis-note",paste0(format(panel$start,big.mark=","),"–",format(panel$end,big.mark=",")," raw bp · ",
        if(grepl("gene_intervals|region_intervals",panel$basis)) "Source gene and region intervals" else
          if(grepl("adjoining_SSU",panel$basis)) "Source regions with adjoining SSU context" else
            if(grepl("landmark_and",panel$basis)) "Mapped region with surrounding raw bases" else "Verified sites; gene boundaries unavailable")),
      tags$div(class="fungal-atlas-scroll",tags$svg(class="fungal-atlas-svg",viewBox=paste(0,0,width,height),role="img",
        `aria-label`=paste(panel$title,"individual verified primer sites on",data$accession),children)))
  })
  tags$div(class="fungal-primer-atlas",`data-version-id`=data$version_id,`data-accession`=data$accession,
    tags$style(htmltools::HTML(".fungal-atlas-panel{margin:0 0 20px}.fungal-atlas-panel h4{font-size:1rem;margin:10px 0 2px}.fungal-atlas-axis-note,.fungal-atlas-context{font-size:.8rem;color:#64748b;margin:0 0 8px}.fungal-atlas-scroll{overflow-x:auto;max-width:100%;border:1px solid #e2e8f0;border-radius:10px;background:#fff}.fungal-atlas-svg{display:block;width:100%;min-width:1000px;height:auto}.fungal-atlas-primer:focus-visible{outline:2px solid #4f46e5}.fungal-atlas-scroll-cue{font-size:.8rem;color:#64748b;margin:4px 0}@media(min-width:1100px){.fungal-atlas-scroll-cue{display:none}}")),
    tags$p(class="fungal-atlas-context",paste0(nrow(data$rows)," uniquely placed primers in this view · UNITE catalog; positions on NCBI ",data$accession," · raw bp. Blue → forward; red ← reverse. Click to select up to five each.")),
    tags$p(class="fungal-atlas-scroll-cue","Scroll sideways to see the complete gene bar and primer labels."),panel_tags,
    tags$details(tags$summary(paste0(nrow(data$missing)," catalog primers have missing, multiple or unavailable placements on this template")),
      tags$p("These site states do not establish amplification failure or absence of the locus. Boundary primers can appear in both panels with the same selection identity."),
      tags$div(class="fungal-atlas-scroll",tags$table(class="table table-sm",tags$thead(tags$tr(tags$th("Primer"),tags$th("Direction"),tags$th("Site evidence"))),
      tags$tbody(lapply(seq_len(nrow(data$missing)),function(i) tags$tr(tags$td(data$missing$primer_name[i]),tags$td(data$missing$direction[i]),
        tags$td(paste0(data$missing$placement_reason[i],if(is.finite(data$missing$n_matches[i]) && data$missing$n_matches[i]>1L)
          paste0(" · ",data$missing$n_matches[i]," compatible matches") else "")))))))))
}
