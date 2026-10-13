# One-template fungal adapter for the shared primer overview. Product widths
# are physical encoded bases, inclusive of primers, and never alignment-column
# distances or pooled/median coordinates from other references.

.fungal_overview_text <- function(row, column, default = "") {
  if (!column %in% names(row) || !length(row[[column]]) || is.na(row[[column]][1L])) return(default)
  as.character(row[[column]][1L])
}

.fungal_overview_number <- function(row, column) {
  suppressWarnings(as.numeric(.fungal_overview_text(row, column, NA_character_)))
}

.fungal_overview_empty_pairs <- function() {
  data.frame(pair_id=character(), pair_label=character(), label=character(), start=integer(), end=integer(),
    pair_start=integer(), pair_end=integer(), forward_start=integer(), forward_end=integer(),
    reverse_start=integer(), reverse_end=integer(), aligned_amplicon_bp=integer(), targeted_region_bp=integer(),
    marker_id=character(), target=character(), use_case=character(), placement_status=character(),
    folmer_overlap_bp=numeric(), n_panel=numeric(), n_scorable=numeric(), evidence_status=character(),
    source_citation=character(), source_url=character(), accession=character(), version_id=character(),
    length_basis=character(), region_match_basis=character(), stringsAsFactors=FALSE)
}

.fungal_overview_empty_primers <- function() {
  data.frame(pair_id=character(), pair_label=character(), primer_name=character(), direction=character(),
    sequence=character(), source_key=character(), catalog_id=character(), start=integer(), end=integer(),
    reported_ta_c=numeric(), reference=character(), primary_reference_key=character(),
    primary_reference_url=character(), primary_reference_status=character(), source_url=character(),
    gene_locus=character(), remarks=character(), site_status=character(), stringsAsFactors=FALSE)
}

.fungal_overview_empty_citations <- function() {
  data.frame(key=character(), short_citation=character(), title=character(), url=character(),
    catalog_id=character(), primer_name=character(), primary_reference_key=character(),
    primary_reference_status=character(), primary_reference_url=character(), source_url=character(),
    version_id=character(), stringsAsFactors=FALSE)
}

.fungal_overview_landmarks <- function(snapshot, accession, length_bp) {
  rows <- snapshot$landmarks
  if (is.null(rows)) rows <- data.frame(landmark_id=character(),label=character(),start=integer(),end=integer())
  rows <- as.data.frame(rows, stringsAsFactors=FALSE)
  if (!all(c("landmark_id","label","start","end") %in% names(rows))) stop("Frozen landmark fields are incomplete.")
  if (!"accession" %in% names(rows)) {
    reference <- snapshot$manifest$reference_accession
    if (length(reference)!=1L || is.na(reference)) reference <- NA_character_
    rows$accession <- rep(reference,nrow(rows))
  }
  rows <- rows[!is.na(rows$accession) & rows$accession==accession,,drop=FALSE]
  for (field in c("source_note","landmark_type","boundary_status","evidence_status")) {
    if (!field %in% names(rows)) rows[[field]] <- rep("",nrow(rows))
    rows[[field]][is.na(rows[[field]])] <- ""
  }
  rows$start <- suppressWarnings(as.numeric(rows$start)); rows$end <- suppressWarnings(as.numeric(rows$end))
  valid <- is.finite(rows$start) & is.finite(rows$end) & rows$start==floor(rows$start) & rows$end==floor(rows$end) &
    rows$start>=1L & rows$end>=rows$start & rows$end<=length_bp
  withheld <- rows$boundary_status %in% c("unresolved","not_detected","invalid_coordinates","withheld","extraction_mismatch")
  rows$placement_reason <- ifelse(!valid,"Invalid raw landmark coordinates",ifelse(withheld,"Frozen landmark evidence withheld",""))
  unresolved <- rows[!valid | withheld,,drop=FALSE]
  rows <- rows[valid & !withheld,,drop=FALSE]
  if (anyDuplicated(rows$landmark_id)) stop("Own-template landmark identities repeat.")
  domain <- rows$landmark_type=="divergent_domain" | rows$evidence_status=="mapped_reference_domain"
  genes <- !domain & rows$label %in% c("18S","SSU","ITS1","5.8S","ITS2","28S","LSU")
  empty_type <- !nzchar(rows$landmark_type)
  rows$landmark_type[empty_type & !domain] <- ifelse(rows$label[empty_type & !domain] %in% c("ITS1","ITS2"),"spacer","gene")
  empty_boundary <- !nzchar(rows$boundary_status)
  rows$boundary_status[empty_boundary & (rows$landmark_type=="conserved_anchor" |
    grepl("conserved.anchor|partial",rows$source_note,ignore.case=TRUE))] <- "partial"
  rows$boundary_status[empty_boundary & grepl("complete",rows$source_note,ignore.case=TRUE)] <- "complete"
  rows$boundary_status[!nzchar(rows$boundary_status)] <- "curated_approximate"
  rows$evidence_status[!nzchar(rows$evidence_status)] <- ifelse(
    grepl("curated|conserved.anchor",rows$source_note[!nzchar(rows$evidence_status)],ignore.case=TRUE),
    "curated_reference","annotated_reference")
  rows$marker_id <- rep("ITS_FUNGAL",nrow(rows))
  rows$raw_bp <- rows$end-rows$start+1L
  rows$display_order <- seq_len(nrow(rows))
  list(landmarks=rows[genes,,drop=FALSE], domains=rows[domain,,drop=FALSE], other=rows[!genes & !domain,,drop=FALSE], unresolved=unresolved)
}

.fungal_overview_region <- function(region, landmarks) {
  if (length(region)!=1L || is.na(region) || !nzchar(region)) stop("Choose one fungal overview region.")
  if (startsWith(region,"landmark:")) {
    id <- substring(region,10L)
    rows <- landmarks[landmarks$landmark_id==id,,drop=FALSE]
    return(list(id=region,rows=rows,specific=TRUE))
  }
  region <- switch(toupper(region),"18S"="SSU","28S"="LSU","5_8S"="5.8S",toupper(region))
  if (!region %in% c("ALL","SSU","ITS","ITS1","5.8S","ITS2","LSU")) stop("Unknown fungal overview region.")
  labels <- switch(region,SSU=c("18S","SSU"),LSU=c("28S","LSU"),ITS=c("ITS1","5.8S","ITS2"),ALL=character(),region)
  list(id=region,rows=landmarks[landmarks$label %in% labels,,drop=FALSE],specific=region %in% c("ITS1","5.8S","ITS2"))
}

.fungal_overview_overlap <- function(start,end,intervals) {
  if (!nrow(intervals)) return(0L)
  bounds <- intervals[,c("start","end"),drop=FALSE]
  bounds$start <- pmax(start,bounds$start); bounds$end <- pmin(end,bounds$end)
  bounds <- bounds[bounds$end>=bounds$start,,drop=FALSE]
  if (!nrow(bounds)) return(0L)
  bounds <- bounds[order(bounds$start,bounds$end),,drop=FALSE]
  total <- 0L; left <- bounds$start[1L]; right <- bounds$end[1L]
  for (i in seq_len(nrow(bounds))) {
    if (bounds$start[i]>right+1L) { total <- total+right-left+1L; left <- bounds$start[i]; right <- bounds$end[i] }
    else right <- max(right,bounds$end[i])
  }
  as.integer(total+right-left+1L)
}

.fungal_overview_component <- function(snapshot,pair,direction) {
  name <- .fungal_overview_text(pair,paste0(direction,"_primer"))
  sequence <- .fungal_overview_text(pair,paste0(direction,"_sequence"))
  catalog <- snapshot$catalog
  index <- which(catalog$primer_name==name & catalog$sequence==sequence & catalog$direction==direction)
  if (length(index)!=1L) return(NULL)
  catalog[index,,drop=FALSE]
}

.fungal_overview_site <- function(evidence,component,direction,sequence) {
  result <- list(start=NA_integer_,end=NA_integer_,status="component_identity_unresolved")
  if (is.null(component)) return(result)
  result$status <- .fungal_overview_text(evidence,paste0(direction,"_status"),"binding_evidence_unavailable")
  if (result$status!="exact_compatible_site") return(result)
  bounds <- c(.fungal_overview_number(evidence,paste0(direction,"_start")),.fungal_overview_number(evidence,paste0(direction,"_end")))
  oligo <- .fungal_overview_text(component,"sequence")
  if (any(!is.finite(bounds)) || any(bounds!=floor(bounds)) || bounds[1]<1L || bounds[2]<bounds[1] ||
      bounds[2]>nchar(sequence) || bounds[2]-bounds[1]+1L!=nchar(oligo)) {
    result$status <- "invalid_raw_binding_coordinates"; return(result)
  }
  if (!exists("iupac") || !exists("reverse_complement",mode="function")) stop("Source R/functions.R before the fungal overview adapter.")
  query <- strsplit(if (direction=="reverse") reverse_complement(oligo) else oligo,"",fixed=TRUE)[[1L]]
  observed <- strsplit(substr(sequence,bounds[1],bounds[2]),"",fixed=TRUE)[[1L]]
  compatible <- vapply(seq_along(query),function(i) observed[i] %in% iupac[[query[i]]],logical(1))
  if (!all(compatible)) { result$status <- "raw_binding_sequence_mismatch"; return(result) }
  result$start <- as.integer(bounds[1]); result$end <- as.integer(bounds[2]); result
}

.fungal_overview_region_match <- function(selection,start,end,components,genes) {
  if (selection$id=="ALL") return(list(matches=TRUE,basis="all_requested_candidates"))
  overlap <- if (is.finite(start) && is.finite(end)) .fungal_overview_overlap(start,end,selection$rows) else 0L
  if (overlap>0L) return(list(matches=TRUE,basis="own_template_landmark_overlap"))
  if (selection$specific) return(list(matches=FALSE,basis="no_own_template_landmark_overlap"))
  # A product fully inside a different annotated region cannot acquire this
  # region from a generic catalog label. Outside annotated spans, catalog locus
  # evidence stays useful but is explicitly less specific than a raw boundary.
  if (is.finite(start) && is.finite(end) && .fungal_overview_overlap(start,end,genes)>=end-start+1L) {
    return(list(matches=FALSE,basis="own_template_landmark_nonoverlap"))
  }
  loci <- paste(vapply(components,function(x) if (is.null(x)) "" else .fungal_overview_text(x,"gene_locus"),character(1)),collapse=" ")
  matches <- switch(selection$id,SSU=grepl("SSU|18S",loci),ITS=grepl("ITS|5[.]8S",loci),LSU=grepl("LSU|28S",loci),FALSE)
  list(matches=matches,basis=if (matches) "reported_component_locus_without_local_region_boundaries" else "no_region_evidence")
}

fungal_primer_overview_default_accession <- function(snapshot,pair_rows=NULL,region="ITS") {
  if (is.null(pair_rows)) pair_rows <- snapshot$pairs
  accessions <- names(snapshot$sequences)
  if (!length(accessions)) return(NULL)
  choices <- lapply(accessions,function(accession) {
    data <- fungal_primer_overview_data(snapshot,pair_rows,accession,region,max_rows=max(1L,nrow(pair_rows)))
    data.frame(accession=accession,n_mapped=data$summary$n_mapped,
      n_complete_regions=sum(data$landmarks$boundary_status=="complete"),
      n_regions=nrow(data$landmarks),length_bp=data$reference_length,stringsAsFactors=FALSE)
  })
  rows <- do.call(rbind,choices)
  rows$accession[order(-rows$n_mapped,-rows$n_complete_regions,-rows$n_regions,rows$length_bp,rows$accession)][1L]
}

fungal_primer_overview_data <- function(snapshot,pair_rows=NULL,accession=NULL,region="ALL",max_rows=20L) {
  required <- c("manifest","sequences","catalog","pairs","pair_sequences","landmarks")
  if (!all(required %in% names(snapshot))) stop("Fungal overview requires one complete frozen snapshot.")
  if (is.null(pair_rows)) pair_rows <- snapshot$pairs
  if (!is.data.frame(pair_rows) || !"pair_id" %in% names(pair_rows) || anyNA(pair_rows$pair_id) || anyDuplicated(pair_rows$pair_id)) stop("Overview candidates need unique frozen pair identities.")
  if (length(max_rows)!=1L || !is.finite(max_rows) || max_rows!=floor(max_rows) || max_rows<0L) stop("Overview row limit must be a nonnegative integer.")
  if (is.null(accession)) accession <- fungal_primer_overview_default_accession(snapshot,pair_rows,region)
  if (length(accession)!=1L || is.na(accession) || !accession %in% names(snapshot$sequences)) stop("Choose a reference retained in this fungal snapshot.")
  sequence <- toupper(unname(snapshot$sequences[[accession]]))
  if (!nzchar(sequence) || !grepl("^[ACGTURYSWKMBDHVNI]+$",sequence)) stop("Overview coordinates require an ungapped raw template.")
  length_bp <- nchar(sequence); version <- snapshot$manifest$id
  if (length(version)!=1L || is.na(version) || !nzchar(version)) stop("Overview snapshot identity is missing.")
  annotations <- .fungal_overview_landmarks(snapshot,accession,length_bp)
  all_landmarks <- rbind(annotations$landmarks,annotations$domains,annotations$other)
  selection <- .fungal_overview_region(region,all_landmarks)
  frozen <- snapshot$pairs
  match_pair <- match(pair_rows$pair_id,frozen$pair_id)
  if (anyNA(match_pair) || anyDuplicated(frozen$pair_id)) stop("Overview pair identity is outside the selected snapshot.")
  # Counts and extra filters supplied by the module may describe its current
  # selected taxon subset. Component identities are always frozen snapshot data.
  for (field in intersect(c("forward_primer","reverse_primer","forward_sequence","reverse_sequence"),names(pair_rows))) {
    if (anyNA(pair_rows[[field]]) || any(as.character(pair_rows[[field]])!=as.character(frozen[[field]][match_pair]))) stop("Overview component identity differs from the snapshot.")
  }
  rows <- frozen[match_pair,,drop=FALSE]
  for (field in setdiff(names(pair_rows),c("forward_primer","reverse_primer","forward_sequence","reverse_sequence"))) rows[[field]] <- pair_rows[[field]]
  evidence <- snapshot$pair_sequences[snapshot$pair_sequences$accession==accession & snapshot$pair_sequences$pair_id %in% rows$pair_id,,drop=FALSE]
  if (anyDuplicated(evidence$pair_id)) stop("Frozen product identity repeats on the selected template.")
  considered <- lapply(seq_len(nrow(rows)),function(i) {
    pair <- rows[i,,drop=FALSE]
    observed <- evidence[evidence$pair_id==pair$pair_id,,drop=FALSE]
    f <- .fungal_overview_component(snapshot,pair,"forward"); r <- .fungal_overview_component(snapshot,pair,"reverse")
    fs <- .fungal_overview_site(observed,f,"forward",sequence); rs <- .fungal_overview_site(observed,r,"reverse",sequence)
    state <- .fungal_overview_text(observed,"status","binding_evidence_unavailable")
    product <- state=="scorable" && fs$status=="exact_compatible_site" && rs$status=="exact_compatible_site" &&
      fs$end<rs$start && isTRUE(.fungal_overview_number(observed,"amplicon_bp")==rs$end-fs$start+1L)
    if (state=="scorable" && !product) state <- "inconsistent_frozen_product_evidence"
    region_match <- .fungal_overview_region_match(selection,if (product) fs$start else NA_real_,if (product) rs$end else NA_real_,list(f,r),annotations$landmarks)
    list(pair=pair,evidence=observed,f=f,r=r,fs=fs,rs=rs,status=state,product=product,region=region_match)
  })
  region_candidates <- Filter(function(row) row$region$matches,considered)
  shown <- head(region_candidates,as.integer(max_rows))
  pair_meta <- .fungal_overview_empty_pairs(); primers <- .fungal_overview_empty_primers(); citations <- .fungal_overview_empty_citations()
  missing <- list()
  for (record in shown) {
    pair <- record$pair
    label <- .fungal_overview_text(pair,"pair_label",pair$pair_id)
    if (record$product) {
      row <- pair
      row$pair_label <- row$label <- label
      row$start <- row$pair_start <- row$forward_start <- record$fs$start
      row$forward_end <- record$fs$end
      row$end <- row$pair_end <- row$reverse_end <- record$rs$end
      row$reverse_start <- record$rs$start
      row$aligned_amplicon_bp <- record$rs$end-record$fs$start+1L
      row$targeted_region_bp <- record$rs$start-record$fs$end-1L
      row$marker_id <- "ITS_FUNGAL"
      row$use_case <- .fungal_overview_text(pair,"use_case","Fungal rDNA primer combination")
      row$placement_status <- "exact_template_sites"
      row$folmer_overlap_bp <- NA_real_
      row$accession <- accession; row$version_id <- version
      row$length_basis <- "raw_template_encoded_span_inclusive_primers"
      row$region_match_basis <- record$region$basis
      # Preserve every supplied metadata field without forcing empty output
      # schema columns to change type or losing candidate support counts.
      if (!nrow(pair_meta)) pair_meta <- row else pair_meta <- rbind(pair_meta,row)
    } else {
      row <- pair
      row$accession <- accession; row$version_id <- version; row$placement_status <- record$status
      row$forward_status <- record$fs$status; row$reverse_status <- record$rs$status
      row$reason <- if (record$status=="inconsistent_frozen_product_evidence") "Frozen product coordinates, sites or length disagree; product placement withheld." else
        "No unique exact-supported inward product is available on this selected raw template; missing sites can reflect fragment coverage or mismatches."
      row$region_match_basis <- record$region$basis
      missing[[length(missing)+1L]] <- row
    }
    for (direction in c("forward","reverse")) {
      component <- record[[if (direction=="forward") "f" else "r"]]
      site <- record[[if (direction=="forward") "fs" else "rs"]]
      if (is.null(component)) next
      citation <- .fungal_overview_text(component,"reference")
      primary_url <- .fungal_overview_text(component,"primary_reference_url")
      source_url <- .fungal_overview_text(component,"source_url")
      key <- paste0("fungal_",digest::digest(paste(version,component$catalog_id,citation,primary_url,source_url,sep="\n"),algo="sha256",serialize=FALSE))
      primer <- data.frame(pair_id=pair$pair_id,pair_label=label,primer_name=component$primer_name,direction=direction,
        sequence=component$sequence,source_key=key,catalog_id=component$catalog_id,start=site$start,end=site$end,
        reported_ta_c=NA_real_,reference=citation,primary_reference_key=.fungal_overview_text(component,"primary_reference_key"),
        primary_reference_url=primary_url,primary_reference_status=.fungal_overview_text(component,"primary_reference_status"),
        source_url=source_url,gene_locus=.fungal_overview_text(component,"gene_locus"),remarks=.fungal_overview_text(component,"remarks"),
        site_status=site$status,stringsAsFactors=FALSE)
      primers <- rbind(primers,primer)
      if (!key %in% citations$key) citations <- rbind(citations,data.frame(key=key,short_citation=if (nzchar(citation)) citation else "Source citation unavailable",
        title=paste(component$primer_name,"·",if (nzchar(citation)) citation else "Source citation unavailable"),url=if (nzchar(primary_url)) primary_url else source_url,
        catalog_id=component$catalog_id,primer_name=component$primer_name,primary_reference_key=.fungal_overview_text(component,"primary_reference_key"),
        primary_reference_status=.fungal_overview_text(component,"primary_reference_status"),primary_reference_url=primary_url,
        source_url=source_url,version_id=version,stringsAsFactors=FALSE))
    }
  }
  missing_pairs <- if (length(missing)) do.call(rbind,missing) else {
    empty <- rows[FALSE,,drop=FALSE]
    for (field in c("accession","version_id","placement_status","forward_status","reverse_status","reason","region_match_basis")) empty[[field]] <- character()
    empty
  }
  # Exact selected-region boundaries take precedence. In All mode the first
  # displayed candidate is the focus, so distant pairs cannot force a whole-
  # operon initial view. Other mapped products may be explicitly clipped by UI.
  if (nrow(selection$rows)) {
    focus <- list(start=min(selection$rows$start),end=max(selection$rows$end),basis="own_template_landmarks",region=selection$id,
                  landmark_ids=selection$rows$landmark_id)
  } else if (nrow(pair_meta)) {
    focus <- list(start=pair_meta$start[1L],end=pair_meta$end[1L],basis="first_displayed_exact_product",region=selection$id,landmark_ids=character())
  } else focus <- list(start=NA_integer_,end=NA_integer_,basis="region_coordinates_unavailable",region=selection$id,landmark_ids=character())
  domain <- if (is.finite(focus$start) && is.finite(focus$end)) c(min=max(1L,focus$start-1000L),max=min(length_bp,focus$end+2000L)) else numeric()
  taxonomy <- snapshot$taxonomy
  reference <- if (is.data.frame(taxonomy)) taxonomy[taxonomy$accession==accession,,drop=FALSE] else data.frame()
  source_url <- .fungal_overview_text(reference,"source_url")
  marker <- data.frame(marker_id="ITS_FUNGAL",display_name="Fungal rDNA",molecule_context="nuclear rDNA",
    coordinate_reference=accession,reference_version=accession,reference_length=length_bp,coordinate_model="raw_template_with_landmarks",
    status="exact_site_reference",source_name="Frozen NCBI template / UNITE primer catalog",source_url=source_url,
    snapshot_id=version,stringsAsFactors=FALSE)
  available <- data.frame(region_id=c("ALL","SSU","ITS","LSU"),label=c("All candidate loci","18S / SSU","ITS region","28S / LSU"),stringsAsFactors=FALSE)
  if (nrow(all_landmarks)) available <- rbind(available,data.frame(region_id=paste0("landmark:",all_landmarks$landmark_id),
    label=paste0(all_landmarks$label," · ",all_landmarks$start,"–",all_landmarks$end," bp · ",all_landmarks$boundary_status),stringsAsFactors=FALSE))
  summary <- data.frame(version_id=version,accession=accession,region=selection$id,n_panel=length(snapshot$sequences),
    n_snapshot_pairs=nrow(snapshot$pairs),n_candidates=nrow(rows),n_region_candidates=length(region_candidates),
    n_shown_candidates=length(shown),n_mapped=nrow(pair_meta),n_missing=nrow(missing_pairs),n_truncated=max(0L,length(region_candidates)-length(shown)),
    stringsAsFactors=FALSE)
  list(pair_meta=pair_meta,primers=primers,citations=citations,marker=marker,landmarks=annotations$landmarks,domains=annotations$domains,
    other_landmarks=annotations$other,unresolved_landmarks=annotations$unresolved,missing_pairs=missing_pairs,summary=summary,
    accession=accession,reference_length=length_bp,reference_metadata=reference,focus=focus,domain=domain,available_regions=available,version_id=version)
}
