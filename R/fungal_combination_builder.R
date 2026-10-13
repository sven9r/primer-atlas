# Build a selected Cartesian set without discarding unsupported combinations.
# Frozen component/pair identities and exact site results remain the inputs.
# New crosses are explicitly derived candidates, never documented assays.

.fungal_combo_bind <- function(frames) {
  frames <- Filter(function(x) is.data.frame(x) && nrow(x)>0L,frames)
  if (!length(frames)) return(data.frame())
  fields <- unique(unlist(lapply(frames,names),use.names=FALSE))
  frames <- lapply(frames,function(x) {
    x <- as.data.frame(x,stringsAsFactors=FALSE)
    for (field in setdiff(fields,names(x))) x[[field]] <- NA
    x[,fields,drop=FALSE]
  })
  rows <- do.call(rbind,frames);rownames(rows)<-NULL;rows
}

.fungal_combo_text <- function(row,field,default="") {
  if (!field %in% names(row) || !length(row[[field]]) || is.na(row[[field]][1L])) return(default)
  as.character(row[[field]][1L])
}

.fungal_combo_keys <- function(rows) paste(rows$primer_name,rows$sequence,rows$direction,sep="\036")
.fungal_combo_pair_keys <- function(rows) paste(rows$forward_primer,rows$forward_sequence,rows$reverse_primer,rows$reverse_sequence,sep="\036")

.fungal_combo_region <- function(snapshot,region) {
  if (length(region)!=1L || is.na(region)) stop("Choose one component region.")
  if (startsWith(region,"landmark:")) {
    rows <- snapshot$landmarks
    id <- substring(region,10L)
    row <- if (is.data.frame(rows)) rows[rows$landmark_id==id,,drop=FALSE] else data.frame()
    label <- if (nrow(row)) .fungal_combo_text(row,"label") else ""
    region <- if (grepl("^D[0-9]",label) || label %in% c("28S","LSU")) "LSU" else
      if (label %in% c("18S","SSU")) "SSU" else if (label %in% c("ITS1","5.8S","ITS2")) "ITS" else "ALL"
  }
  region <- switch(toupper(region),"18S"="SSU","28S"="LSU","ITS1"="ITS","ITS2"="ITS","5.8S"="ITS",toupper(region))
  if (!region %in% c("ALL","SSU","ITS","LSU")) stop("Unknown component region.")
  region
}

fungal_combination_catalog <- function(snapshot,region="ALL",catalog=NULL) {
  base <- as.data.frame(snapshot$catalog,stringsAsFactors=FALSE)
  required <- c("catalog_id","primer_name","sequence","direction","gene_locus")
  if (!all(required %in% names(base)) || anyNA(base$catalog_id) || anyDuplicated(base$catalog_id)) stop("Frozen component catalog identity is incomplete.")
  base$component_origin <- "frozen_snapshot_catalog"
  base$component_source_version <- rep(snapshot$manifest$id,nrow(base))
  if (!is.null(catalog)) {
    extra <- as.data.frame(catalog,stringsAsFactors=FALSE)
    if (!all(required %in% names(extra)) || anyNA(extra$catalog_id) || anyDuplicated(extra$catalog_id)) stop("Explicit supplemental catalog identity is incomplete.")
    repeated <- extra$catalog_id %in% base$catalog_id
    if (any(repeated)) {
      original <- base[match(extra$catalog_id[repeated],base$catalog_id),,drop=FALSE]
      if (any(.fungal_combo_keys(extra[repeated,,drop=FALSE])!=.fungal_combo_keys(original))) stop("A supplemental row changes a frozen component identity.")
    }
    # Existing identities always keep the historical source metadata. An
    # explicit supplemental catalog cannot quietly replace their citations.
    extra <- extra[!repeated,,drop=FALSE]
    if (nrow(extra)) {
      extra$component_origin <- "explicit_supplemental_catalog"
      extra$component_source_version <- if ("source_version" %in% names(extra)) as.character(extra$source_version) else
        if ("source_snapshot_utc" %in% names(extra)) as.character(extra$source_snapshot_utc) else "unversioned_explicit_catalog"
      base <- .fungal_combo_bind(list(base,extra))
    }
  }
  if (anyDuplicated(.fungal_combo_keys(base))) stop("Multiple catalog IDs share one name/sequence/role; select its canonical frozen identity.")
  if (anyNA(base$direction) || any(!base$direction %in% c("forward","reverse"))) stop("Select components with an explicit effective forward/reverse orientation.")
  if (anyNA(base$sequence) || any(!nzchar(base$sequence)) || any(base$sequence!=toupper(gsub("[[:space:]]","",base$sequence)))) stop("Component oligos must retain normalized source sequences.")
  locus <- as.character(base$gene_locus);locus[is.na(locus)]<-""
  region <- .fungal_combo_region(snapshot,region)
  eligible <- switch(region,SSU=grepl("SSU|18S",locus),ITS=grepl("ITS|5[.]8S",locus),LSU=grepl("LSU|28S",locus),rep(TRUE,nrow(base)))
  base <- base[eligible,,drop=FALSE]
  status <- if ("primary_reference_status" %in% names(base)) as.character(base$primary_reference_status) else rep("",nrow(base))
  status[is.na(status)] <- ""
  base$publication_flag <- ifelse(status=="publication_cited","cited component",ifelse(status=="unpublished","unpublished component","provenance incomplete"))
  base$region_label <- ifelse(is.na(base$gene_locus) | !nzchar(base$gene_locus),"locus unresolved",base$gene_locus)
  base$selection_label <- paste(base$primer_name,base$region_label,base$publication_flag,sep=" · ")
  duplicate_label <- duplicated(base$selection_label) | duplicated(base$selection_label,fromLast=TRUE)
  base$selection_label[duplicate_label] <- paste0(base$selection_label[duplicate_label]," · ",base$sequence[duplicate_label])
  base$parent_version_id <- rep(snapshot$manifest$id,nrow(base))
  rownames(base)<-NULL;base
}

.fungal_combo_sparse_proof <- function(snapshot) {
  proof <- snapshot$search_manifest
  base <- if (is.data.frame(snapshot$frozen_component_catalog)) snapshot$frozen_component_catalog else snapshot$catalog
  artifacts <- snapshot$manifest$artifacts
  ids <- names(snapshot$sequences)
  fields <- c(sequence_sha256="sequences",catalog_artifact_sha256="catalog",pair_artifact_sha256="pairs",site_artifact_sha256="sites")
  checks <- vapply(names(fields),function(field) {
    descriptor <- artifacts[[fields[[field]]]]
    !is.null(descriptor$sha256) && identical(proof[[field]],descriptor$sha256)
  },logical(1))
  valid <- !is.null(proof) && isTRUE(proof$complete) && identical(proof$scan_rule_version,"exact_iupac_unique_sites_v1") &&
    isTRUE(proof$n_records==length(ids)) && isTRUE(proof$n_catalog==nrow(base)) && isTRUE(proof$n_sparse_rows==nrow(snapshot$sites)) &&
    isTRUE(proof$n_site_searches==length(ids)*nrow(base)) && isTRUE(proof$n_no_exact_hits+proof$n_sparse_rows==proof$n_site_searches) &&
    identical(proof$record_ids_sha256,digest::digest(paste(sort(ids),collapse="\n"),algo="sha256",serialize=FALSE)) && all(checks)
  if (!valid) stop("Omitted UNITE hits require a complete search proof tied to this frozen cohort and catalog.")
  TRUE
}

.fungal_combo_sites <- function(snapshot,components,accessions,max_new_site_searches=5000L) {
  if(!length(accessions)) {
    empty<-.fungal_empty_sites();empty$site_evidence_origin<-character();return(empty)
  }
  sequences <- snapshot$sequences[accessions]
  cache <- snapshot$sites
  if (!is.data.frame(cache)) cache <- data.frame()
  sparse <- nrow(cache)>0L && "record_id" %in% names(cache)
  if (identical(snapshot$manifest$marker_id,"UNITE_ITS")) sparse <- TRUE
  if (sparse) .fungal_combo_sparse_proof(snapshot)
  owner <- if (sparse) "record_id" else "accession"
  keys <- if (nrow(cache)) paste(cache[[owner]],cache$catalog_id,sep="\036") else character()
  if (anyDuplicated(keys)) stop("Frozen site identity repeats.")
  if (nrow(cache) && any(!cache[[owner]] %in% names(snapshot$sequences))) stop("A frozen site belongs to another cohort.")
  frozen_catalog <- if(is.data.frame(snapshot$frozen_component_catalog)) snapshot$frozen_component_catalog else snapshot$catalog
  if(nrow(cache) && any(!cache$catalog_id %in% frozen_catalog$catalog_id)) stop("A frozen site belongs to an unknown component.")
  fresh <- components$component_origin=="explicit_supplemental_catalog"
  n_fresh <- sum(fresh)*length(accessions)
  allow_fresh <- n_fresh<=max_new_site_searches
  blocks <- lapply(seq_len(nrow(components)),function(i) {
    component <- components[i,,drop=FALSE]
    index <- match(paste(accessions,component$catalog_id,sep="\036"),keys)
    rows <- data.frame(catalog_id=rep(component$catalog_id,length(accessions)),primer_name=rep(component$primer_name,length(accessions)),
      direction=rep(component$direction,length(accessions)),accession=accessions,start=rep(NA_integer_,length(accessions)),end=rep(NA_integer_,length(accessions)),
      n_matches=integer(length(accessions)),status=rep("binding_evidence_unavailable",length(accessions)),
      site_evidence_origin=rep("frozen_site_row_unavailable",length(accessions)),stringsAsFactors=FALSE)
    found <- !is.na(index)
    for (field in c("start","end","n_matches","status")) if (any(found)) rows[[field]][found] <- cache[[field]][index[found]]
    rows$site_evidence_origin[found] <- "frozen_exact_site_results"
    omitted <- !found & !fresh[i]
    if (sparse && any(omitted)) {
      rows$status[omitted] <- "no_exact_binding_site"
      rows$site_evidence_origin[omitted] <- "verified_exhaustive_sparse_search_no_hit"
    }
    if (fresh[i] && length(accessions)) {
      if (allow_fresh) {
        computed <- lapply(sequences,function(sequence) .fungal_exact_site(strsplit(toupper(sequence),"",fixed=TRUE)[[1L]],component$sequence,component$direction))
        rows$start <- vapply(computed,`[[`,integer(1),"start");rows$end <- vapply(computed,`[[`,integer(1),"end")
        rows$n_matches <- vapply(computed,`[[`,integer(1),"n_matches");rows$status <- vapply(computed,`[[`,character(1),"status")
        rows$site_evidence_origin <- "new_exact_search_explicit_supplemental_component"
      } else rows$site_evidence_origin <- "supplemental_search_exceeds_interactive_bound"
    }
    rows$start <- suppressWarnings(as.integer(rows$start));rows$end <- suppressWarnings(as.integer(rows$end));rows$n_matches <- as.integer(rows$n_matches)
    # Independently check every cached positive site against its own raw bases.
    for (j in which(rows$status=="exact_compatible_site")) {
      bounds <- c(rows$start[j],rows$end[j]);sequence <- toupper(unname(sequences[j]))
      query <- strsplit(if (component$direction=="reverse") reverse_complement(component$sequence) else component$sequence,"",fixed=TRUE)[[1L]]
      valid <- all(is.finite(bounds)) && bounds[1]>=1L && bounds[2]>=bounds[1] && bounds[2]<=nchar(sequence) && diff(bounds)+1L==length(query) && rows$n_matches[j]==1L
      if (valid) {
        observed <- strsplit(substr(sequence,bounds[1],bounds[2]),"",fixed=TRUE)[[1L]]
        valid <- all(vapply(seq_along(query),function(k) observed[k] %in% iupac[[query[k]]],logical(1)))
      }
      if (!valid) stop("Frozen exact binding coordinates do not match their raw template and oligo.")
    }
    rows
  })
  if (!length(blocks)) {
    empty <- .fungal_empty_sites();empty$site_evidence_origin<-character();return(empty)
  }
  do.call(rbind,blocks)
}

fungal_combination_defaults <- function(snapshot,accession=NULL,region="ITS",n_forward=5L,n_reverse=5L,catalog=NULL) {
  if (any(!is.finite(c(n_forward,n_reverse))) || any(c(n_forward,n_reverse)<0) || any(c(n_forward,n_reverse)!=floor(c(n_forward,n_reverse)))) stop("Default component counts must be nonnegative integers.")
  components <- fungal_combination_catalog(snapshot,region,catalog)
  if (is.null(accession)) {
    supported <- snapshot$pair_sequences
    counts <- if (is.data.frame(supported)) table(supported$accession[supported$status=="scorable"]) else integer()
    ids <- names(snapshot$sequences)
    accession <- if (length(counts)) names(counts)[order(-as.integer(counts),names(counts))][1L] else if (length(ids)) ids[1L] else NULL
  }
  sites <- if (length(accession)==1L && accession %in% names(snapshot$sequences)) .fungal_combo_sites(snapshot,components,accession) else .fungal_empty_sites()
  exact <- sites$catalog_id[sites$status=="exact_compatible_site"]
  priorities <- list(forward=c("ITS1","ITS1F","ITS5","ITS3","fITS7","gITS7","ITS86F","NS1","LR0R"),
                     reverse=c("ITS2","ITS4","LR21","LR3","LR5","LR7","ITS4B","5.8S_Fungi"))
  pick <- function(role,n) {
    rows <- components[components$direction==role,,drop=FALSE]
    priority <- match(rows$primer_name,priorities[[role]]);priority[is.na(priority)]<-999L
    rows$catalog_id[head(order(!rows$catalog_id %in% exact,priority,rows$publication_flag!="cited component",rows$primer_name,rows$sequence),n)]
  }
  list(forward_ids=pick("forward",n_forward),reverse_ids=pick("reverse",n_reverse),accession=accession,region=.fungal_combo_region(snapshot,region))
}

fungal_combination_builder <- function(snapshot,forward_ids,reverse_ids,accession=NULL,accessions=NULL,catalog=NULL,
    platform="ALL",length_range=NULL,max_per_direction=5L,illumina_max=500,pacbio_max=2000,nanopore_max=3500,
    evidence_filter=NULL,region="ALL",max_pair_template_rows=1000000L,max_new_site_searches=5000L) {
  if (!all(c("manifest","catalog","pairs","sequences") %in% names(snapshot))) stop("A frozen component/pair/sequence snapshot is required.")
  if (!exists(".fungal_exact_site",mode="function") || !exists("fungal_platform_summary",mode="function")) stop("Source R/functions.R and R/fungal_explorer.R before building combinations.")
  forward_ids <- unique(as.character(forward_ids));reverse_ids<-unique(as.character(reverse_ids))
  if (anyNA(c(forward_ids,reverse_ids)) || any(!nzchar(c(forward_ids,reverse_ids)))) stop("Selected component identities are missing.")
  if (length(max_per_direction)!=1L || !is.finite(max_per_direction) || max_per_direction!=floor(max_per_direction) || max_per_direction<1L || max_per_direction>20L) stop("Choose an explicit per-direction bound from 1 to 20.")
  if(length(max_new_site_searches)!=1L || !is.finite(max_new_site_searches) || max_new_site_searches<0L || max_new_site_searches!=floor(max_new_site_searches)) stop("New interactive site-search bounds must be a nonnegative integer.")
  if (length(forward_ids)>max_per_direction || length(reverse_ids)>max_per_direction) stop("Selected oligos exceed the explicit per-direction bound; none were silently discarded.")
  universe <- fungal_combination_catalog(snapshot,"ALL",catalog)
  fi <- match(forward_ids,universe$catalog_id);ri<-match(reverse_ids,universe$catalog_id)
  if (anyNA(c(fi,ri)) || any(universe$direction[fi]!="forward") || any(universe$direction[ri]!="reverse")) stop("Each selected oligo must match its exact catalog identity and effective role.")
  if (is.null(accessions)) accessions <- names(snapshot$sequences)
  accessions <- unique(as.character(accessions))
  if (anyNA(accessions) || any(!accessions %in% names(snapshot$sequences))) stop("The selected reference cohort contains unknown identities.")
  if (!is.null(accession) && (length(accession)!=1L || is.na(accession) || !accession %in% accessions)) stop("The map reference must belong to the selected cohort.")
  sequences <- snapshot$sequences[accessions]
  if (anyNA(sequences) || any(!grepl("^[ACGTURYSWKMBDHVNI]+$",sequences,ignore.case=TRUE))) stop("Raw reference templates must be nonempty ungapped nucleotide strings.")
  n_combinations <- length(fi)*length(ri)
  if (!is.finite(max_pair_template_rows) || max_pair_template_rows<0 || n_combinations*length(accessions)>max_pair_template_rows) stop("Selected Cartesian set exceeds the interactive pair/template bound; reduce its scope without dropping combinations.")
  if (!is.null(length_range) && (length(length_range)!=2L || any(!is.finite(length_range)) || length_range[1]<1L || length_range[2]<length_range[1])) stop("Length bounds must be an ordered positive min/max pair.")
  categories <- c("documented_pair","candidate_published_components","candidate_unpublished_component","candidate_incomplete_provenance")
  if (identical(evidence_filter,"ALL")) evidence_filter<-NULL
  if (!is.null(evidence_filter) && any(!evidence_filter %in% categories)) stop("Unknown combination evidence category.")
  components <- universe[universe$catalog_id %in% c(forward_ids,reverse_ids),,drop=FALSE]
  sites <- .fungal_combo_sites(snapshot,components,accessions,max_new_site_searches)
  cross <- expand.grid(reverse=ri,forward=fi,KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
  known_keys <- .fungal_combo_pair_keys(snapshot$pairs)
  if (anyDuplicated(known_keys) || anyDuplicated(snapshot$pairs$pair_id)) stop("Frozen pair component identity is ambiguous.")
  definitions <- list();blocks<-list()
  site_keys <- paste(sites$catalog_id,sites$accession,sep="\036")
  for (i in seq_len(nrow(cross))) {
    f <- universe[cross$forward[i],,drop=FALSE];r<-universe[cross$reverse[i],,drop=FALSE]
    key <- paste(f$primer_name,f$sequence,r$primer_name,r$sequence,sep="\036");existing<-match(key,known_keys)
    if (!is.na(existing)) {
      pair <- as.data.frame(snapshot$pairs[existing,,drop=FALSE],stringsAsFactors=FALSE);origin<-"frozen_pair_identity"
    } else {
      status <- c(.fungal_combo_text(f,"primary_reference_status"),.fungal_combo_text(r,"primary_reference_status"))
      evidence <- if ("unpublished" %in% status) "candidate_unpublished_component" else
        if (all(status=="publication_cited")) "candidate_published_components" else "candidate_incomplete_provenance"
      combine <- function(field,sep=" + ") paste(unique(c(.fungal_combo_text(f,field),.fungal_combo_text(r,field)))[nzchar(unique(c(.fungal_combo_text(f,field),.fungal_combo_text(r,field))))],collapse=sep)
      pair <- data.frame(pair_id=paste0("FUNGAL_COMB_",substr(digest::digest(paste(f$catalog_id,r$catalog_id,sep="\036"),algo="sha256",serialize=FALSE),1L,24L)),
        pair_label=paste(f$primer_name,r$primer_name,sep=" + "),forward_primer=f$primer_name,reverse_primer=r$primer_name,
        forward_sequence=f$sequence,reverse_sequence=r$sequence,target=combine("target"),evidence_status=evidence,
        source_citation=combine("reference"),source_url=combine("primary_reference_url"," | "),region=combine("gene_locus"),stringsAsFactors=FALSE)
      origin<-"new_selected_cartesian_candidate"
    }
    pair$forward_catalog_id<-f$catalog_id;pair$reverse_catalog_id<-r$catalog_id;pair$pair_origin<-origin
    pair$parent_version_id<-snapshot$manifest$id
    pair$forward_source_citation<-.fungal_combo_text(f,"reference");pair$reverse_source_citation<-.fungal_combo_text(r,"reference")
    pair$forward_source_url<-.fungal_combo_text(f,"primary_reference_url",.fungal_combo_text(f,"source_url"));pair$reverse_source_url<-.fungal_combo_text(r,"primary_reference_url",.fungal_combo_text(r,"source_url"))
    fs<-if(length(accessions)) sites[match(paste(f$catalog_id,accessions,sep="\036"),site_keys),,drop=FALSE] else sites[FALSE,,drop=FALSE]
    rs<-if(length(accessions)) sites[match(paste(r$catalog_id,accessions,sep="\036"),site_keys),,drop=FALSE] else sites[FALSE,,drop=FALSE]
    good_f<-fs$status=="exact_compatible_site";good_r<-rs$status=="exact_compatible_site"
    state<-rep("both_sites_missing",length(accessions));state[good_f & !good_r]<-"reverse_site_missing";state[!good_f & good_r]<-"forward_site_missing"
    state[fs$status=="low_specificity_site" | rs$status=="low_specificity_site"]<-"low_specificity_site"
    state[fs$status=="ambiguous_binding_site" | rs$status=="ambiguous_binding_site"]<-"ambiguous_binding_site"
    state[fs$status %in% c("binding_evidence_unavailable","invalid_oligo","template_shorter_than_oligo") |
          rs$status %in% c("binding_evidence_unavailable","invalid_oligo","template_shorter_than_oligo")]<-"binding_evidence_unavailable"
    both<-good_f & good_r;state[both]<-"outward_or_overlapping_sites"
    bp<-rs$end-fs$start+1L;inward<-both & fs$end<rs$start
    state[inward]<-"below_minimum_length";state[inward & bp>=100L]<-"scorable";bp[state!="scorable"]<-NA_integer_
    block<-data.frame(pair_id=rep(pair$pair_id,length(accessions)),accession=accessions,forward_start=fs$start,forward_end=fs$end,
      reverse_start=rs$start,reverse_end=rs$end,amplicon_bp=bp,status=state,forward_status=fs$status,reverse_status=rs$status,
      forward_site_origin=fs$site_evidence_origin,reverse_site_origin=rs$site_evidence_origin,
      forward_catalog_id=rep(f$catalog_id,length(accessions)),reverse_catalog_id=rep(r$catalog_id,length(accessions)),
      parent_version_id=rep(snapshot$manifest$id,length(accessions)),stringsAsFactors=FALSE)
    # Existing dense frozen product rows are a second independent input. New
    # candidates cannot silently replace a conflicting result in that history.
    old<-snapshot$pair_sequences
    if (!is.na(existing) && is.data.frame(old)) {
      old<-old[old$pair_id==pair$pair_id & old$accession %in% accessions,,drop=FALSE]
      if (anyDuplicated(old$accession)) stop("Frozen product identity repeats on a template.")
      if (nrow(old)) {
        current<-block[match(old$accession,block$accession),,drop=FALSE]
        for (field in c("forward_start","forward_end","reverse_start","reverse_end","amplicon_bp","status","forward_status","reverse_status")) {
          if (!identical(as.character(old[[field]]),as.character(current[[field]]))) stop("Derived exact geometry disagrees with a frozen pair result.")
        }
      }
    }
    lengths<-bp[state=="scorable"]
    pair$n_panel<-length(accessions);pair$n_scorable<-length(lengths)
    pair$amplicon_min_bp<-if(length(lengths)) min(lengths) else NA_real_;pair$amplicon_median_bp<-if(length(lengths)) median(lengths) else NA_real_;pair$amplicon_max_bp<-if(length(lengths)) max(lengths) else NA_real_
    definitions[[i]]<-pair;blocks[[i]]<-block
  }
  selected_pairs<-if(length(definitions)) .fungal_combo_bind(definitions) else .fungal_empty_pairs()
  pair_sequences<-if(length(blocks)) do.call(rbind,blocks) else .fungal_empty_pair_sequences()
  taxonomy<-snapshot$taxonomy
  if(!is.data.frame(taxonomy) && is.data.frame(snapshot$metadata)) {
    taxonomy<-snapshot$metadata
    taxonomy$accession<-taxonomy$record_id
    if("class_name" %in% names(taxonomy)) taxonomy$class<-taxonomy$class_name
  }
  pair_sequences<-.fungal_add_taxonomy(pair_sequences,taxonomy)
  if(is.data.frame(snapshot$metadata) && "record_id" %in% names(snapshot$metadata) && "accession" %in% names(snapshot$metadata)) {
    pair_sequences$source_accession<-snapshot$metadata$accession[match(pair_sequences$accession,snapshot$metadata$record_id)]
  }
  selected_pairs<-fungal_platform_summary(pair_sequences,selected_pairs,platform,illumina_max,pacbio_max,nanopore_max)
  platform<-toupper(platform);ceiling<-switch(platform,ALL=Inf,ILLUMINA=illumina_max,PACBIO=pacbio_max,NANOPORE=nanopore_max)
  eligible<-pair_sequences$status=="scorable" & is.finite(pair_sequences$amplicon_bp) & pair_sequences$amplicon_bp<=ceiling
  if (!is.null(length_range)) eligible<-eligible & pair_sequences$amplicon_bp>=length_range[1L] & pair_sequences$amplicon_bp<=length_range[2L]
  counts<-table(pair_sequences$pair_id[eligible]);selected_pairs$n_filtered_products<-as.integer(counts[match(selected_pairs$pair_id,names(counts))]);selected_pairs$n_filtered_products[is.na(selected_pairs$n_filtered_products)]<-0L
  if (is.null(accession) && length(accessions)) {
    counts<-table(pair_sequences$accession[eligible]);accession<-if(length(counts)) names(counts)[order(-as.integer(counts),names(counts))][1L] else accessions[1L]
  }
  current<-if(length(accession)==1L) pair_sequences[pair_sequences$accession==accession,,drop=FALSE] else pair_sequences[FALSE,,drop=FALSE]
  selected_pairs$reference_status<-if(nrow(current)) current$status[match(selected_pairs$pair_id,current$pair_id)] else rep("no_template_selected",nrow(selected_pairs))
  selected_pairs$reference_amplicon_bp<-if(nrow(current)) current$amplicon_bp[match(selected_pairs$pair_id,current$pair_id)] else rep(NA_integer_,nrow(selected_pairs))
  reasons<-lapply(seq_len(nrow(selected_pairs)),function(i) {
    row<-selected_pairs[i,,drop=FALSE];codes<-character()
    if(row$reference_status!="scorable") codes<-c(codes,"no_verified_product_on_reference")
    if(!is.null(evidence_filter) && !row$evidence_status %in% evidence_filter) codes<-c(codes,"evidence_filter")
    if(row$reference_status=="scorable") {
      if(!is.null(length_range) && (row$reference_amplicon_bp<length_range[1L] || row$reference_amplicon_bp>length_range[2L])) codes<-c(codes,"length_filter")
      if(row$reference_amplicon_bp>ceiling) codes<-c(codes,"platform_filter")
    }
    paste(codes,collapse=" | ")
  })
  selected_pairs$filter_reason<-if(length(reasons)) unlist(reasons,use.names=FALSE) else character();selected_pairs$included_in_map<-!nzchar(selected_pairs$filter_reason)
  working<-snapshot;working$catalog<-universe;working$frozen_component_catalog<-snapshot$catalog
  remaining_pairs<-snapshot$pairs[!snapshot$pairs$pair_id %in% selected_pairs$pair_id,,drop=FALSE]
  working$pairs<-.fungal_combo_bind(list(remaining_pairs,selected_pairs))
  old_sequences<-snapshot$pair_sequences
  if(is.data.frame(old_sequences)) old_sequences<-old_sequences[!(old_sequences$pair_id %in% selected_pairs$pair_id & old_sequences$accession %in% accessions),,drop=FALSE]
  working$pair_sequences<-.fungal_combo_bind(list(old_sequences,pair_sequences))
  working$derived_combination_context<-list(parent_version_id=snapshot$manifest$id,forward_catalog_ids=forward_ids,reverse_catalog_ids=reverse_ids,
    accessions=accessions,site_rule="Unique exact IUPAC-compatible sites on each raw template; inward nonoverlap; >=100 encoded bp including primers",
    new_pair_status="Inferred selected candidate, never automatically a documented pair")
  if(region!="ALL" && nrow(selected_pairs) && length(accession)==1L) {
    if(!exists("fungal_primer_overview_data",mode="function")) stop("Source R/fungal_primer_overview.R for actual landmark-aware regional filtering.")
    overview<-fungal_primer_overview_data(working,selected_pairs,accession,region,max_rows=nrow(selected_pairs))
    matching<-c(overview$pair_meta$pair_id,overview$missing_pairs$pair_id)
    excluded<-!selected_pairs$pair_id %in% matching
    selected_pairs$filter_reason[excluded]<-ifelse(nzchar(selected_pairs$filter_reason[excluded]),paste(selected_pairs$filter_reason[excluded],"region_filter",sep=" | "),"region_filter")
    selected_pairs$included_in_map<-!nzchar(selected_pairs$filter_reason)
    working$pairs[match(selected_pairs$pair_id,working$pairs$pair_id),names(selected_pairs)]<-selected_pairs
  }
  mapped<-selected_pairs[selected_pairs$included_in_map,,drop=FALSE]
  unmapped<-selected_pairs[selected_pairs$reference_status!="scorable",,drop=FALSE]
  excluded<-selected_pairs[selected_pairs$reference_status=="scorable" & !selected_pairs$included_in_map,,drop=FALSE]
  summary<-data.frame(parent_version_id=snapshot$manifest$id,accession=if(length(accession)==1L) accession else NA_character_,
    n_forward=length(forward_ids),n_reverse=length(reverse_ids),n_selected_combinations=n_combinations,n_panel=length(accessions),
    n_frozen_pair_identities=sum(selected_pairs$pair_origin=="frozen_pair_identity"),n_new_candidate_identities=sum(selected_pairs$pair_origin=="new_selected_cartesian_candidate"),
    n_verified_on_reference=sum(selected_pairs$reference_status=="scorable"),n_mapped_after_filters=nrow(mapped),n_unmapped_on_reference=nrow(unmapped),
    n_verified_filtered_out=nrow(excluded),n_unavailable_site_searches=sum(sites$status=="binding_evidence_unavailable"),stringsAsFactors=FALSE)
  list(selected_pairs=selected_pairs,pairs=selected_pairs,pair_sequences=pair_sequences,sites=sites,pair_rows=mapped,map_pair_rows=mapped,
    unmapped_pairs=unmapped,filtered_out_pairs=excluded,components=components,summary=summary,working_snapshot=working,
    accession=accession,filters=list(platform=platform,length_range=length_range,evidence_filter=evidence_filter,region=region,
    figure_pair_ids=mapped$pair_id),version_id=snapshot$manifest$id)
}
