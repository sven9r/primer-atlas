# Publication artifacts share one measured context. SVG and tables never turn
# alignment columns, a panel median, or an absent site into template product bp.
.pub_text <- function(x, default="") {
  if(is.null(x) || !length(x) || is.na(x[1L])) return(default)
  as.character(x[1L])
}
.pub_field <- function(row,name,default="") .pub_text(if(name %in% names(row)) row[[name]] else NULL,default)
.pub_number <- function(row,name,default=NA_real_) suppressWarnings(as.numeric(.pub_field(row,name,as.character(default))))
.pub_frame <- function(x) if(is.null(x)) data.frame() else as.data.frame(x,stringsAsFactors=FALSE)
.pub_or <- function(x,y) if(is.null(x)) y else x
.pub_integer <- function(x) is.finite(x) & x==floor(x)
.pub_escape <- function(x) as.character(htmltools::htmlEscape(as.character(x),attribute=TRUE))
.pub_count <- function(x) if(length(x)!=1L || is.na(x)) "unavailable" else format(x,trim=TRUE,big.mark=",",scientific=FALSE)
.pub_evidence_label <- function(status) {
  labels<-c(documented_pair="Documented combination",candidate_published_components="Candidate from published components",
    candidate_unpublished_component="Candidate with an unpublished component",candidate_incomplete_provenance="Candidate with incomplete provenance")
  if(status %in% names(labels)) unname(labels[status]) else gsub("_"," ",status,fixed=TRUE)
}
.pub_reason_label <- function(reason) {
  if(!nzchar(reason)) return("")
  codes<-strsplit(reason," | ",fixed=TRUE)[[1L]]
  labels<-c(no_map="No verified product on this reference",no_verified_product_on_reference="No verified product on this reference",
    length="Outside selected length",length_filter="Outside selected length",platform="Outside platform length",platform_filter="Outside platform length",
    evidence="Outside selected evidence category",evidence_filter="Outside selected evidence category",region_filter="Outside selected region",
    region_or_focus="Outside selected region or focus",not_in_current_figure_selection="Outside current figure selection")
  paste(vapply(codes,function(code) if(code %in% names(labels)) unname(labels[code]) else gsub("_"," ",code,fixed=TRUE),character(1)),collapse="; ")
}
.pub_rc <- function(sequence) paste(rev(strsplit(chartr("ACGTURYSWKMBDHVNI","TGCAAYRSWMKVHDBNI",toupper(sequence)),"",fixed=TRUE)[[1L]]),collapse="")
.pub_iupac <- list(A="A",C="C",G="G",T="T",U="T",R=c("A","G"),Y=c("C","T"),S=c("G","C"),W=c("A","T"),
  K=c("G","T"),M=c("A","C"),B=c("C","G","T"),D=c("A","G","T"),H=c("A","C","T"),V=c("A","C","G"),N=c("A","C","G","T"),I=c("A","C","G","T"))
.pub_site_matches <- function(sequence,oligo,start,end,reverse=FALSE) {
  query <- strsplit(if(reverse) .pub_rc(oligo) else toupper(oligo),"",fixed=TRUE)[[1L]]
  observed <- strsplit(substr(sequence,start,end),"",fixed=TRUE)[[1L]]
  length(query)==length(observed) && length(query)>0L && all(vapply(seq_along(query),function(i) observed[i] %in% .pub_iupac[[query[i]]],logical(1)))
}
.pub_bind <- function(rows) {
  if(!length(rows)) return(data.frame())
  names <- unique(unlist(lapply(rows,names),use.names=FALSE))
  rows <- lapply(rows,function(row) {for(name in setdiff(names,names(row))) row[[name]] <- NA;row[,names,drop=FALSE]})
  result <- do.call(rbind,rows);rownames(result)<-NULL;result
}
.pub_orientation_fields <- c("parent_catalog_id","parent_primer_name","parent_sequence","parent_direction",
  "parent_source_direction","parent_primary_reference_status","parent_primary_reference_key","parent_primary_reference_url",
  "parent_reference","parent_source_url","parent_reported_position","parent_source_version","source_version",
  "orientation_derived","orientation_origin","derivation_rule","raw_query_equivalence","sequence_derivation_note","parent_noncanonical_symbols")
.pub_orientation_bool <- function(value) {
  text<-tolower(.pub_text(value))
  if(!nzchar(text)) return(NA)
  if(!text %in% c("true","false","1","0")) stop("Orientation provenance has an invalid boolean.")
  text %in% c("true","1")
}
.pub_component_orientation <- function(primer,pair,side,catalog) {
  id<-.pub_field(primer,"catalog_id",.pub_field(pair,paste0(side,"_catalog_id")))
  original<-if(nrow(catalog) && "catalog_id" %in% names(catalog)) catalog[catalog$catalog_id==id,,drop=FALSE] else data.frame()
  if(nrow(original)>1L) stop("Orientation component catalog identity repeats.")
  fields<-.pub_orientation_fields[vapply(.pub_orientation_fields,function(field)
    field %in% names(primer) || field %in% names(original) || paste0(side,"_",field) %in% names(pair),logical(1))]
  if(!length(setdiff(fields,"source_version")) && !startsWith(id,"FUNGAL_RC_")) return(NULL)
  values<-lapply(fields,function(field) {
    candidates<-c(.pub_field(primer,field),.pub_field(original,field),.pub_field(pair,paste0(side,"_",field)))
    candidates<-unique(candidates[nzchar(candidates)])
    if(length(candidates)>1L) stop("Orientation metadata differs between catalog, component and pair.")
    if(length(candidates)) candidates else ""
  });names(values)<-fields
  derived<-.pub_orientation_bool(values$orientation_derived)
  origin<-.pub_text(values$orientation_origin)
  if(startsWith(id,"FUNGAL_RC_") || identical(origin,"derived_reverse_complement")) {
    if(!isTRUE(derived)) stop("Derived orientation identity needs an explicit derivation flag.")
  }
  if(isTRUE(derived)) {
    required<-c("parent_catalog_id","parent_primer_name","parent_sequence","parent_direction","parent_source_version")
    if(any(!vapply(required,function(field) nzchar(.pub_text(values[[field]])),logical(1))) ||
       !identical(origin,"derived_reverse_complement") || !identical(.pub_text(values$derivation_rule),"reverse_complement_opposite_role_v1") ||
       !isTRUE(.pub_orientation_bool(values$raw_query_equivalence))) stop("Derived orientation provenance is incomplete.")
    parent<-toupper(.pub_text(values$parent_sequence));effective<-toupper(.pub_field(primer,"sequence",.pub_field(pair,paste0(side,"_sequence"))))
    expected<-if(exists("reverse_complement",mode="function")) reverse_complement(parent) else .pub_rc(gsub("I","N",parent,fixed=TRUE))
    parent_role<-.pub_text(values$parent_direction)
    if(!grepl("^[ACGTURYSWKMBDHVNI]+$",parent) || !parent_role %in% c("forward","reverse") || parent_role==side ||
       !identical(effective,expected) || !startsWith(id,"FUNGAL_RC_")) stop("Derived oligo must be the declared parent's reverse complement in the opposite role.")
    if(nrow(original)==1L && (!identical(.pub_field(original,"sequence"),effective) || !identical(.pub_field(original,"direction"),side)))
      stop("Derived effective oligo or role disagrees with its catalog row.")
    parent_row<-if(nrow(catalog)) catalog[catalog$catalog_id==.pub_text(values$parent_catalog_id),,drop=FALSE] else data.frame()
    if(nrow(parent_row)>1L || (nrow(parent_row)==1L &&
       (!identical(.pub_field(parent_row,"sequence"),parent) || !identical(.pub_field(parent_row,"direction"),parent_role) ||
        !identical(.pub_field(parent_row,"primer_name"),.pub_text(values$parent_primer_name))))) stop("Derived parent metadata disagrees with its exact catalog row.")
    status<-.pub_field(primer,"primary_reference_status",.pub_field(original,"primary_reference_status"))
    if(status!="derived_reverse_complement" || .pub_field(pair,"evidence_status") %in% c("documented_pair","candidate_published_components"))
      stop("A derived orientation cannot inherit publication support for its oligo or pair.")
    if(grepl("I",parent,fixed=TRUE) && (!nzchar(.pub_text(values$sequence_derivation_note)) || !grepl("I",.pub_text(values$parent_noncanonical_symbols),fixed=TRUE)))
      stop("Inosine normalization requires an explicit source-symbol and derivation note.")
    values$orientation_citation_scope<-"parent_primer_only"
    values$orientation_citation_caveat<-"Original references, where available, describe the parent primer only; this derived reverse-complement oligo and combination have no publication or PCR validation."
  }
  values$orientation_derived<-derived
  if("raw_query_equivalence" %in% names(values)) values$raw_query_equivalence<-.pub_orientation_bool(values$raw_query_equivalence)
  values$catalog_id<-id;values$effective_direction<-side
  values
}

primer_publication_context <- function(overview, selected_pair_ids=NULL, reference_sequence=NULL,
  broad_evidence=NULL, broad_manifest=NULL, filters=list(), domain=NULL, title="Primer combinations",citation_lookup=NULL) {
  required <- c("pair_meta","primers","citations","marker","landmarks","missing_pairs","version_id","accession","reference_length")
  if(!is.list(overview) || !all(required %in% names(overview))) stop("Publication export needs a verified shared-overview context.")
  mapped <- .pub_frame(overview$pair_meta);missing <- .pub_frame(overview$missing_pairs)
  available <- c(as.character(mapped$pair_id),as.character(missing$pair_id))
  if(anyNA(available) || anyDuplicated(available)) stop("Publication pair/template identities repeat.")
  ids <- if(is.null(selected_pair_ids)) available else as.character(selected_pair_ids)
  if(!length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) || any(!ids %in% available)) stop("Every requested publication pair must be present in the export context.")
  reference <- .pub_text(overview$accession);version <- .pub_text(overview$version_id)
  size <- as.numeric(overview$reference_length)
  if(!nzchar(reference) || !nzchar(version) || length(size)!=1L || !.pub_integer(size) || size<1L) stop("Publication reference/version/length is unavailable.")
  if(!is.null(reference_sequence)) {
    reference_sequence <- toupper(as.character(reference_sequence))
    if(length(reference_sequence)!=1L || is.na(reference_sequence) || nchar(reference_sequence)!=size ||
       !grepl("^[ACGTURYSWKMBDHVNI]+$",reference_sequence)) stop("Publication raw template is invalid or gapped.")
  }
  primers <- .pub_frame(overview$primers)
  component_catalog <- .pub_frame(overview$component_catalog)
  figure_ids <- if(is.null(filters$figure_pair_ids)) ids[ids %in% mapped$pair_id] else as.character(filters$figure_pair_ids)
  if(anyNA(figure_ids) || anyDuplicated(figure_ids) || any(!figure_ids %in% ids) || any(!figure_ids %in% mapped$pair_id)) stop("Figure pairs must be selected and have measured template products.")
  if(is.null(domain)) domain <- .pub_or(overview$domain,c(1,size))
  domain <- as.numeric(domain)
  if(length(domain)!=2L || any(!is.finite(domain)) || domain[1L]<0 || domain[2L]>size || domain[2L]<=domain[1L]) domain <- c(1,size)
  registry <- data.frame(citation_id=character(),role=character(),citation=character(),url=character(),status=character(),primary_reference_key=character(),stringsAsFactors=FALSE)
  register <- function(role,citation,url="",status="",primary_reference_key="") {
    citation <- .pub_text(citation,"Citation unavailable");url<-.pub_text(url);status<-.pub_text(status)
    if(!nzchar(citation)) citation<-"Citation unavailable"
    primary_reference_key<-.pub_text(primary_reference_key)
    same <- which(registry$role==role & registry$citation==citation & registry$url==url & registry$status==status & registry$primary_reference_key==primary_reference_key)
    if(length(same)) return(registry$citation_id[same[1L]])
    id <- paste0("C",nrow(registry)+1L)
    registry <<- rbind(registry,data.frame(citation_id=id,role=role,citation=citation,url=url,status=status,primary_reference_key=primary_reference_key,stringsAsFactors=FALSE));id
  }
  marker <- .pub_frame(overview$marker)
  reference_url <- .pub_field(marker,"source_url")
  if(!nzchar(reference_url) && grepl("^[A-Za-z]+_?[0-9]+[.][0-9]+$",reference)) reference_url <- paste0("https://www.ncbi.nlm.nih.gov/nuccore/",reference)
  ncbi_citation <- register("coordinate reference",paste("NCBI",reference,"; frozen reference",version),reference_url,"frozen accession/version")
  source_citations <- .pub_frame(overview$citations)
  primer_rows <- list();geometry_rows <- list()
  for(i in seq_along(ids)) {
    id <- ids[i];exact <- id %in% mapped$pair_id
    row <- if(exact) mapped[mapped$pair_id==id,,drop=FALSE] else missing[missing$pair_id==id,,drop=FALSE]
    f <- primers[primers$pair_id==id & primers$direction=="forward",,drop=FALSE]
    r <- primers[primers$pair_id==id & primers$direction=="reverse",,drop=FALSE]
    if(nrow(f)>1L || nrow(r)>1L) stop("Publication pair components repeat.")
    oligo <- function(primer,side) .pub_field(primer,"sequence",.pub_field(row,paste0(side,"_sequence")))
    fseq<-oligo(f,"forward");rseq<-oligo(r,"reverse")
    f_orientation<-.pub_component_orientation(f,row,"forward",component_catalog)
    r_orientation<-.pub_component_orientation(r,row,"reverse",component_catalog)
    fref<-.pub_field(f,"reference","Citation unavailable");rref<-.pub_field(r,"reference","Citation unavailable")
    furl<-.pub_field(f,"primary_reference_url",.pub_field(f,"source_url"));rurl<-.pub_field(r,"primary_reference_url",.pub_field(r,"source_url"))
    fkey<-.pub_field(f,"primary_reference_key");rkey<-.pub_field(r,"primary_reference_key")
    fcite<-register("primer component",fref,furl,.pub_field(f,"primary_reference_status","unavailable"),fkey)
    rcite<-register("primer component",rref,rurl,.pub_field(r,"primary_reference_status","unavailable"),rkey)
    start<-end<-fs<-fe<-rs<-re<-bp<-internal<-NA_real_
    flen<-if(nzchar(fseq)) nchar(fseq) else NA_real_;rlen<-if(nzchar(rseq)) nchar(rseq) else NA_real_
    if(exact) {
      fs<-.pub_number(row,"forward_start");fe<-.pub_number(row,"forward_end");rs<-.pub_number(row,"reverse_start");re<-.pub_number(row,"reverse_end")
      bounds<-c(fs,fe,rs,re)
      if(any(!.pub_integer(bounds)) || fs<1L || fe<fs || rs<=fe || re<rs || re>size ||
         !isTRUE(fe-fs+1L==flen) || !isTRUE(re-rs+1L==rlen)) stop("Publication primer spans or oligo lengths disagree.")
      start<-fs;end<-re;bp<-end-start+1L;internal<-rs-fe-1L
      if(!isTRUE(bp==.pub_number(row,"aligned_amplicon_bp")) || !isTRUE(internal==.pub_number(row,"targeted_region_bp")) || !isTRUE(bp-internal==flen+rlen)) stop("Publication product lengths disagree with inclusive template coordinates.")
      if(!is.null(reference_sequence) && (!.pub_site_matches(reference_sequence,fseq,fs,fe) || !.pub_site_matches(reference_sequence,rseq,rs,re,TRUE))) stop("Publication primer sites disagree with raw template bases.")
    } else {
      fs<-.pub_number(f,"start");fe<-.pub_number(f,"end");rs<-.pub_number(r,"start");re<-.pub_number(r,"end")
    }
    included <- id %in% figure_ids
    reason <- if(included) "" else if(!exact) "no_map" else "not_in_current_figure_selection"
    explicit <- filters$exclusion_reasons
    if(!included && nzchar(.pub_field(row,"filter_reason"))) reason<-.pub_field(row,"filter_reason")
    if(!included && !is.null(explicit) && id %in% names(explicit)) reason <- .pub_text(explicit[[id]],reason)
    if(!included && exact && reason=="not_in_current_figure_selection") {
      minbp<-suppressWarnings(as.numeric(.pub_or(filters$min_bp,.pub_or(filters$minbp,0))));maxbp<-suppressWarnings(as.numeric(.pub_or(filters$max_bp,.pub_or(filters$maxbp,Inf))))
      if(length(minbp)==1L && length(maxbp)==1L && is.finite(bp) && (bp<minbp || bp>maxbp)) reason<-"length"
      platform<-toupper(.pub_text(filters$platform,"ALL"))
      ceiling<-switch(platform,ILLUMINA=as.numeric(.pub_or(filters$illumina_max,500)),PACBIO=as.numeric(.pub_or(filters$pacbio_max,2000)),NANOPORE=as.numeric(.pub_or(filters$nanopore_max,3500)),Inf)
      if(reason=="not_in_current_figure_selection" && length(ceiling)==1L && is.finite(ceiling) && bp>ceiling) reason<-"platform"
      evidence<-tolower(.pub_text(filters$evidence,"all"))
      if(reason=="not_in_current_figure_selection" && grepl("document",evidence) && .pub_field(row,"evidence_status")!="documented_pair") reason<-"evidence"
      if(reason=="not_in_current_figure_selection" && nzchar(.pub_text(filters$region)) && toupper(.pub_text(filters$region))!="ALL") reason<-"region_or_focus"
    }
    ambiguous <- if(exact && !is.null(reference_sequence)) nchar(gsub("[ACGT]","",substr(reference_sequence,start,end))) else NA_real_
    geometry_rows[[i]]<-data.frame(selection_row=i,figure_row=if(included) match(id,figure_ids) else NA_integer_,pair_id=id,
      pair_label=.pub_field(row,"pair_label",id),forward_primer=.pub_field(f,"primer_name",.pub_field(row,"forward_primer")),reverse_primer=.pub_field(r,"primer_name",.pub_field(row,"reverse_primer")),
      forward_oligo_5to3=fseq,reverse_oligo_5to3=rseq,forward_length_bp=flen,reverse_length_bp=rlen,
      forward_start=fs,forward_end=fe,reverse_start=rs,reverse_end=re,product_start=start,product_end=end,
      product_including_primers_bp=bp,product_excluding_primers_bp=internal,n_ambiguous_product_bases=ambiguous,
      coordinate_basis="1-based inclusive encoded raw-template bases",reference_accession=reference,reference_length_bp=size,ncbi_snapshot_id=version,
      ncbi_snapshot_references=if(is.data.frame(overview$summary)) .pub_number(overview$summary,"n_panel") else NA_real_,
      ncbi_selected_references=.pub_number(row,"n_panel"),ncbi_exact_products=.pub_number(row,"n_scorable"),
      placement_status=.pub_field(row,"placement_status",if(exact) "exact_template_sites" else "unavailable"),
      combination_status=.pub_field(row,"evidence_status","candidate_incomplete_provenance"),
      combination_source_citation=.pub_field(row,"source_citation"),combination_source_url=.pub_field(row,"source_url"),
      forward_citation_id=fcite,reverse_citation_id=rcite,coordinate_citation_id=ncbi_citation,
      forward_reference=fref,reverse_reference=rref,forward_reference_url=furl,reverse_reference_url=rurl,
      forward_primary_reference_key=fkey,reverse_primary_reference_key=rkey,
      forward_reference_status=.pub_field(f,"primary_reference_status","unavailable"),reverse_reference_status=.pub_field(r,"primary_reference_status","unavailable"),
      included_in_figure=included,exclusion_reason=reason,
      window_state=if(!exact) "unmapped" else if(end<domain[1L] || start>domain[2L]) "outside" else if(start<domain[1L] || end>domain[2L]) "clipped" else "full",
      stringsAsFactors=FALSE)
    for(side in c("forward","reverse")) {
      orientation<-if(side=="forward") f_orientation else r_orientation
      if(!is.null(orientation)) for(field in names(orientation)) geometry_rows[[i]][[paste0(side,"_",field)]]<-orientation[[field]]
    }
    for(component in list(f,r)) if(nrow(component)) {
      component$source_key<-if(component$direction[1L]=="forward") fcite else rcite
      primer_rows[[length(primer_rows)+1L]]<-component
    }
  }
  geometry <- .pub_bind(geometry_rows)
  annotations <- .pub_bind(Filter(function(x) is.data.frame(x) && nrow(x),list(overview$landmarks,overview$domains)))
  if(nrow(annotations)) {
    if(any(!.pub_integer(annotations$start)) || any(!.pub_integer(annotations$end)) || any(annotations$start<1L | annotations$end<annotations$start | annotations$end>size)) stop("Publication landmark coordinates are invalid.")
    annotations$citation_id<-vapply(seq_len(nrow(annotations)),function(i) {
      row<-annotations[i,,drop=FALSE]
      if(!nzchar(.pub_field(row,"source_citation")) && !nzchar(.pub_field(row,"source_url"))) return(ncbi_citation)
      register("region/domain annotation",.pub_field(row,"source_citation",.pub_field(row,"source_note","Annotation source")),.pub_field(row,"source_url"),.pub_field(row,"evidence_status"))
    },character(1))
  }
  table <- geometry
  broad_meta <- if(is.list(broad_evidence)) .pub_or(broad_evidence$metadata,list()) else list()
  broad_summary <- .pub_frame(if(is.list(broad_evidence)) broad_evidence$summary else NULL)
  broad_version <- .pub_text(broad_meta$version_id,.pub_text(if(is.list(broad_manifest)) broad_manifest$id else NULL))
  if(!is.null(broad_manifest) && nzchar(broad_version) && !identical(broad_version,.pub_text(broad_manifest$id))) stop("Broad evidence and manifest versions differ.")
  broad_doi <- .pub_text(broad_meta$source_doi,.pub_text(if(is.list(broad_manifest)) broad_manifest$source_doi else NULL))
  broad_citation <- if(nzchar(broad_version)) register("UNITE reference dataset",paste("UNITE",broad_version,if(nzchar(broad_doi)) paste(";",broad_doi) else ""),if(nzchar(broad_doi)) paste0("https://doi.org/",sub("^https://doi.org/","",broad_doi)) else "https://unite.ut.ee/repository.php","frozen retained ITS references") else ""
  broad_rows <- list()
  for(i in seq_len(nrow(geometry))) {
    row<-geometry[i,,drop=FALSE]
    scope<-if(nrow(broad_summary) && "pair_id" %in% names(broad_summary)) broad_summary[broad_summary$pair_id==row$pair_id,,drop=FALSE] else data.frame()
    if(nrow(scope) && !nzchar(broad_version)) stop("Broad taxonomic counts require an explicit immutable source version.")
    if(!nrow(scope)) scope<-data.frame(unavailable=TRUE)
    for(j in seq_len(nrow(scope))) {
      tax<-scope[j,,drop=FALSE]
      denominator<-if("n_selected" %in% names(tax)) .pub_number(tax,"n_selected") else .pub_number(tax,"n_panel")
      exact<-.pub_number(tax,"n_exact_products")
      if(nzchar(.pub_field(tax,"version_id")) && .pub_field(tax,"version_id")!=broad_version) stop("Broad row and source versions differ.")
      if(is.finite(denominator) && (!.pub_integer(denominator) || denominator<0L)) stop("Broad taxonomic denominator is invalid.")
      if(is.finite(exact) && (!.pub_integer(exact) || exact<0L || !is.finite(denominator) || exact>denominator)) stop("Broad exact-product count exceeds its selected denominator.")
      count_fields<-c("n_exact_products","n_single_sites","n_no_exact_sites","n_unresolved","n_unavailable","n_withheld")
      if(all(count_fields %in% names(tax))) {
        counts<-vapply(count_fields,function(field) .pub_number(tax,field),numeric(1))
        if(any(!.pub_integer(counts) | counts<0L) || !is.finite(denominator) || sum(counts)!=denominator) stop("Broad outcome counts must partition the selected reference denominator.")
      }
      out<-row
      out$unite_snapshot_id<-broad_version;out$unite_source_doi<-broad_doi;out$unite_citation_id<-broad_citation
      out$taxonomic_rank<-.pub_field(tax,"rank",.pub_text(broad_meta$rank))
      out$taxon_group_id<-.pub_field(tax,"group_id",.pub_text(broad_meta$group_id))
      out$taxon<-.pub_field(tax,"taxon","Evidence unavailable")
      out$parent_path<-.pub_field(tax,"parent_path")
      out$unite_selected_references<-denominator;out$unite_exact_products<-exact
      out$unite_exact_product_fraction<-if(is.finite(denominator) && denominator>0L && is.finite(exact)) exact/denominator else NA_real_
      out$unite_length_eligible<-.pub_number(tax,"n_length_eligible")
      out$unite_one_site<-.pub_number(tax,if("n_single_sites" %in% names(tax)) "n_single_sites" else "n_partial_sites")
      out$unite_no_exact_sites<-.pub_number(tax,"n_no_exact_sites")
      out$unite_unresolved<-.pub_number(tax,if("n_unresolved" %in% names(tax)) "n_unresolved" else "n_unresolved_sites")
      out$unite_unavailable<-.pub_number(tax,"n_unavailable")
      out$unite_product_withheld<-.pub_number(tax,if("n_withheld" %in% names(tax)) "n_withheld" else "n_product_withheld")
      out$taxonomic_evidence_status<-if(is.finite(denominator) && is.finite(exact)) "retained_reference_exact_site_counts" else "unavailable"
      broad_rows[[length(broad_rows)+1L]]<-out
    }
  }
  table<-.pub_bind(broad_rows)
  if(anyDuplicated(paste(table$pair_id,table$taxon_group_id,sep="\n"))) stop("Publication pair/taxon rows repeat.")
  context<-list(overview=overview,selected_pair_ids=ids,figure_pair_ids=figure_ids,geometry_table=geometry,table=table,
    primers=.pub_bind(primer_rows),citations=registry,landmarks=annotations,marker=marker,domain=domain,title=.pub_text(title,"Primer combinations"),
    reference_accession=reference,reference_length=size,version_id=version,reference_sequence_sha256=if(is.null(reference_sequence)) "" else digest::digest(reference_sequence,algo="sha256",serialize=FALSE),
    broad_version_id=broad_version,broad_source_doi=broad_doi,filters=filters,verified_citations=citation_lookup)
  class(context)<-c("primer_publication_context","list");.pub_verified_context(context)
}

fungal_publication_context <- function(snapshot,pair_ids=NULL,accession=NULL,broad_evidence=NULL,filters=list(),
  region="ALL",domain=NULL,title="Fungal primer combinations",pair_rows=NULL,broad_manifest=NULL,citation_lookup=NULL) {
  if(!exists("fungal_primer_overview_data",mode="function")) stop("Source R/fungal_primer_overview.R before the publication adapter.")
  if(is.null(pair_ids)) pair_ids<-if(is.null(pair_rows)) snapshot$pairs$pair_id else pair_rows$pair_id
  if(is.null(pair_rows)) pair_rows<-snapshot$pairs[match(pair_ids,snapshot$pairs$pair_id),,drop=FALSE]
  # All requested rows enter the table. Figure filtering is separate and explicit.
  overview<-fungal_primer_overview_data(snapshot,pair_rows,accession,region="ALL",max_rows=length(pair_ids))
  # The shared geometry adapter retains effective primer identity; the exact
  # catalog supplies optional original/derived orientation provenance.
  overview$component_catalog<-snapshot$catalog
  if(is.null(filters$region)) filters$region<-region
  figure_ids<-if(is.null(filters$figure_pair_ids)) overview$pair_meta$pair_id else filters$figure_pair_ids
  figure_rows<-pair_rows[pair_rows$pair_id %in% figure_ids,,drop=FALSE]
  focused<-fungal_primer_overview_data(snapshot,figure_rows,overview$accession,
    region=.pub_text(filters$region,region),max_rows=nrow(figure_rows))
  overview$focus<-focused$focus
  if(is.null(domain) && length(focused$domain)==2L) domain<-focused$domain
  primer_publication_context(overview,pair_ids,snapshot$sequences[[overview$accession]],broad_evidence,broad_manifest,filters,domain,title,citation_lookup)
}

.pub_verified_context <- function(context) {
  citations<-context$citations
  if(!"catalog_citation" %in% names(citations)) citations$catalog_citation<-citations$citation
  if(!"catalog_url" %in% names(citations)) citations$catalog_url<-citations$url
  citations$citation<-citations$catalog_citation;citations$url<-citations$catalog_url
  citations$doi<-citations$verification_url<-citations$verification_note<-""
  citations$bibliography_status<-ifelse(citations$role=="primer component","short_citation_only_full_reference_unavailable","source_record")
  lookup<-.pub_frame(context$verified_citations)
  if(nrow(lookup)) {
    if(!all(c("primary_reference_key","full_citation") %in% names(lookup)) || anyNA(lookup$primary_reference_key) || anyDuplicated(lookup$primary_reference_key)) stop("Verified bibliography lookup must have unique exact primary-reference keys.")
    for(i in seq_len(nrow(citations))) {
      match<-which(lookup$primary_reference_key==citations$primary_reference_key[i] & nzchar(citations$primary_reference_key[i]))
      if(!length(match)) next
      row<-lookup[match,,drop=FALSE]
      full<-.pub_field(row,"full_citation")
      if(!nzchar(full)) next
      citations$citation[i]<-full;citations$doi[i]<-.pub_field(row,"doi")
      citations$url[i]<-.pub_field(row,"url",citations$url[i])
      citations$verification_url[i]<-.pub_field(row,"verification_url")
      citations$verification_note[i]<-.pub_field(row,"verification_note")
      citations$bibliography_status[i]<-"verified_bibliography_supplement"
    }
  }
  enrich<-function(table) {
    for(side in c("forward","reverse")) {
      source<-match(table[[paste0(side,"_citation_id")]],citations$citation_id)
      table[[paste0(side,"_full_reference")]]<-citations$citation[source]
      table[[paste0(side,"_full_reference_url")]]<-citations$url[source]
      table[[paste0(side,"_doi")]]<-citations$doi[source]
      table[[paste0(side,"_bibliography_status")]]<-citations$bibliography_status[source]
      table[[paste0(side,"_verification_url")]]<-citations$verification_url[source]
    }
    table
  }
  context$citations<-citations;context$table<-enrich(context$table);context$geometry_table<-enrich(context$geometry_table);context
}
.pub_derived_components <- function(geometry) {
  rows<-list()
  fields<-c(.pub_orientation_fields,"catalog_id","effective_direction","orientation_citation_scope","orientation_citation_caveat")
  for(side in c("forward","reverse")) {
    flag<-paste0(side,"_orientation_derived")
    if(!flag %in% names(geometry)) next
    for(i in which(geometry[[flag]] %in% TRUE)) {
      row<-geometry[i,c("selection_row","figure_row","pair_id","pair_label"),drop=FALSE]
      row$direction<-side;row$effective_oligo_5to3<-geometry[[paste0(side,"_oligo_5to3")]][i]
      for(field in fields) if(paste0(side,"_",field) %in% names(geometry)) row[[field]]<-geometry[[paste0(side,"_",field)]][i]
      rows[[length(rows)+1L]]<-row
    }
  }
  .pub_bind(rows)
}

.pub_find_svg <- function(node) {
  if(inherits(node,"shiny.tag")) {
    if(identical(node$name,"svg")) return(node)
    for(child in node$children) {value<-.pub_find_svg(child);if(!is.null(value)) return(value)}
  } else if(is.list(node)) for(child in node) {value<-.pub_find_svg(child);if(!is.null(value)) return(value)}
  NULL
}
.pub_svg_document <- function(svg,title,caption,methods,citations,metadata=list()) {
  tags<-htmltools::tags
  if(is.null(svg)) svg<-tags$svg(viewBox="0 0 1680 110",tags$text(x=16,y=70,`font-size`=16,"No exact-supported product is drawable for this template and selection."))
  bounds<-as.numeric(strsplit(svg$attribs$viewBox," ",fixed=TRUE)[[1L]])
  width<-bounds[3L];original_height<-bounds[4L]
  missing_bibliography<-if("bibliography_status" %in% names(citations)) citations$bibliography_status=="short_citation_only_full_reference_unavailable" else rep(FALSE,nrow(citations))
  refs<-if(nrow(citations)) paste0("[",citations$citation_id,"] ",citations$citation,ifelse(nzchar(citations$status),paste0(" · ",citations$status),""),ifelse(missing_bibliography," · Full bibliography unavailable",""),ifelse(nzchar(citations$url),paste0(" · ",citations$url),"")) else character()
  lines<-unlist(lapply(c("Caption",caption,"Methods",methods,"References",refs),function(text) strwrap(text,width=max(60L,floor(width/7.4)))),use.names=FALSE)
  svg$attribs$xmlns<-"http://www.w3.org/2000/svg"
  svg$attribs$width<-width;svg$attribs$height<-original_height+90L+length(lines)*18L
  svg$attribs$viewBox<-paste(0,0,width,svg$attribs$height)
  svg$attribs$style<-"font-family:Arial,Helvetica,sans-serif;background:#ffffff"
  svg$attribs$`aria-label`<-title
  children<-svg$children
  svg$children<-list(tags$title(title),tags$desc(paste(caption,methods)),tags$metadata(jsonlite::toJSON(metadata,auto_unbox=TRUE,null="null")),
    tags$rect(x=0,y=0,width=width,height=svg$attribs$height,fill="#ffffff"),
    tags$text(x=16,y=23,`font-size`=20,`font-weight`=700,title),
    tags$g(transform="translate(0 40)",children))
  for(i in seq_along(lines)) svg<-htmltools::tagAppendChild(svg,tags$text(x=16,y=original_height+66L+(i-1L)*18L,fill="#263238",`font-size`=13,lines[i]))
  # The installed SVG renderer lacks these two hyphen glyphs. Normalize only
  # display text; source strings in CSV/HTML and structured metadata are intact.
  display_text<-function(node,normalize=FALSE) {
    if(inherits(node,"shiny.tag")) {
      child_text<-node$name %in% c("text","tspan","title","desc")
      node$children<-lapply(node$children,display_text,normalize=child_text)
    } else if(is.list(node)) node<-lapply(node,display_text,normalize=normalize)
    else if(is.character(node) && normalize) node<-gsub("[\u2010\u2011]","-",node,perl=TRUE)
    node
  }
  svg<-display_text(svg)
  paste0('<?xml version="1.0" encoding="UTF-8"?>\n',as.character(svg),'\n')
}

.pub_filters_text <- function(filters) {
  values<-filters[setdiff(names(filters),c("figure_pair_ids","exclusion_reasons"))]
  if(!length(values)) return("No additional selection filters were supplied.")
  parts<-vapply(names(values),function(name) paste0(name,"=",paste(as.character(values[[name]]),collapse=",")),character(1))
  paste("Selection filters:",paste(parts,collapse="; "))
}

primer_publication_export <- function(context) {
  if(!inherits(context,"primer_publication_context")) stop("Create one verified publication context before exporting.")
  context<-.pub_verified_context(context)
  mapped<-context$overview$pair_meta[match(context$figure_pair_ids,context$overview$pair_meta$pair_id),,drop=FALSE]
  geometry<-context$geometry_table;annotations<-context$landmarks
  orientation<-.pub_derived_components(geometry)
  annotation_text<-if(nrow(annotations)) paste(paste0(annotations$label," (",annotations$boundary_status,")"),collapse=", ") else "No gene/domain coordinates are available for this reference."
  products<-geometry[geometry$included_in_figure,,drop=FALSE]
  lengths<-if(nrow(products)) paste0("Figure product spans range from ",min(products$product_including_primers_bp)," to ",max(products$product_including_primers_bp)," encoded bp including primers.") else "No exact-supported product is displayed."
  ncbi_total<-unique(geometry$ncbi_snapshot_references[is.finite(geometry$ncbi_snapshot_references)])
  broad_total<-unique(context$table$unite_selected_references[is.finite(context$table$unite_selected_references)])
  evidence_text<-if(nrow(products)) paste(vapply(unique(products$combination_status),function(status) {
    rows<-sort(products$figure_row[products$combination_status==status])
    paste0(.pub_evidence_label(status),": figure ",if(length(rows)==1L) "row " else "rows ",paste(rows,collapse=", "),".")
  },character(1)),collapse=" ") else ""
  caption<-paste(context$title,"on",context$reference_accession,"(reference",context$version_id,").",
    nrow(geometry),"requested combinations;",length(context$figure_pair_ids),"measured product rows in the figure.",
    if(length(ncbi_total)==1L) paste("The NCBI snapshot contains",.pub_count(ncbi_total),"reference templates; geometry uses the named template.") else "",
    "The axis uses 1-based encoded raw-template bp; forward and reverse oligos are reported 5′–3′.",lengths,
    "The internal span excludes both primer-binding lengths. Gene/region annotations:",annotation_text,
    if(nzchar(context$broad_version_id)) paste("Taxonomic counts use the separate UNITE ITS panel",context$broad_version_id,
      if(length(broad_total)==1L) paste("with",.pub_count(broad_total),"selected references per combination.") else "and the selected taxon denominators shown in the table.") else "UNITE taxonomic evidence was not supplied.",
    evidence_text,"Exact-site counts describe retained reference sequences; amplification performance has not been validated.")
  methods<-paste("Product start is the forward-site start; product end is the reverse-site end. Inclusive encoded length is end − start + 1; internal length is reverse-start − forward-end − 1.",
    if(nzchar(context$reference_sequence_sha256)) "Measured coordinates and oligo spans were checked against the selected raw reference; template ambiguity cannot establish an exact primer site." else "Measured coordinates and oligo lengths were reconciled; raw template bases were not supplied for a second site check.",
    "Alignment gaps contribute no product bases.",
    "NCBI geometry and UNITE taxonomic support retain independent source versions and denominators. Missing or filtered combinations remain in the table with their inclusion flag and exclusion reason.",
    "Gene boundaries retain their source completeness; mapped LSU domains retain comparative-projection status.",
    "Full bibliography supplements are matched by exact frozen primary-reference key; missing supplements retain the original short citation and an explicit unavailable-full-reference status.",.pub_filters_text(context$filters))
  if(nrow(orientation)) {
    components<-orientation[!duplicated(orientation$catalog_id),,drop=FALSE]
    descriptions<-vapply(seq_len(nrow(components)),function(i) {
      row<-components[i,,drop=FALSE]
      paste0(row$parent_primer_name," (parent ",row$parent_direction," 5′–3′ ",row$parent_sequence,
        ") → derived ",row$effective_direction," 5′–3′ ",row$effective_oligo_5to3,
        if(nzchar(.pub_field(row,"sequence_derivation_note"))) paste0("; ",row$sequence_derivation_note) else "")
    },character(1))
    caption<-paste(caption,"Derived orientation:",paste(descriptions,collapse=". "),
      "Original citations, where available, describe parent primers only; derived oligos and their combinations are unvalidated computational candidates.")
    methods<-paste(methods,"The derived oligo is the parent's reverse complement in the opposite effective role; original IDs, roles, sequences and source versions are retained.",
      "Raw-query equivalence concerns the nucleotide-set matching model. Inosine-to-N normalization, when present, is computational and does not assert chemical identity or PCR performance.")
  }
  svg<-NULL
  if(nrow(mapped)) {
    if(!exists("primer_overview_map",mode="function")) stop("Source R/primer_overview.R before scientific figure export.")
    source<-context$citations
    plot_citations<-data.frame(key=source$citation_id,short_citation=paste0("[",source$citation_id,"]"),title=source$citation,url=source$url,stringsAsFactors=FALSE)
    map<-primer_overview_map(mapped,context$primers,context$marker,context$landmarks,
      focus=list(start=context$domain[1L],end=context$domain[2L]),domain=context$domain,citations=plot_citations)
    svg<-.pub_find_svg(map)
  }
  metadata<-list(reference_accession=context$reference_accession,ncbi_version=context$version_id,unite_version=context$broad_version_id,
    raw_reference_sha256=context$reference_sequence_sha256,coordinate_basis="raw-template encoded bp;1-based inclusive",selected_pair_ids=context$selected_pair_ids,
    figure_pair_ids=context$figure_pair_ids,filters=context$filters)
  if(nrow(orientation)) metadata$derived_component_orientation<-orientation
  result<-list(svg=.pub_svg_document(svg,context$title,caption,methods,context$citations,metadata),table=context$table,
    geometry_table=context$geometry_table,caption=caption,methods=methods,citations=context$citations,landmarks=context$landmarks,context=context)
  result$html<-primer_publication_html(result)
  class(result)<-c("primer_publication_export","list");result
}

primer_publication_html <- function(export) {
  table<-export$table;escape<-.pub_escape
  number<-function(x) if(length(x)!=1L || is.na(x)) "—" else .pub_count(x)
  rows<-vapply(seq_len(nrow(table)),function(i) {
    row<-table[i,,drop=FALSE]
    oligo_cell<-function(side) {
      content<-paste0(escape(row[[paste0(side,"_primer")]]),"<br><span class=\"oligo\">",escape(row[[paste0(side,"_oligo_5to3")]]),"</span>")
      if(isTRUE(.pub_orientation_bool(row[[paste0(side,"_orientation_derived")]]))) {
        content<-paste0(content,"<br><small>Derived reverse complement; parent ",escape(.pub_field(row,paste0(side,"_parent_primer_name"))),
          " (",escape(.pub_field(row,paste0(side,"_parent_direction"))),") 5′–3′: <span class=\"oligo\">",escape(.pub_field(row,paste0(side,"_parent_sequence"))),
          "</span>. Original citation, where available, applies to the parent only. ",escape(.pub_field(row,paste0(side,"_sequence_derivation_note"))),"</small>")
      }
      content
    }
    cells<-c(paste0(if(is.na(row$figure_row)) "—" else row$figure_row," · ",escape(row$pair_label)),
      oligo_cell("forward"),oligo_cell("reverse"),
      paste0(number(row$product_including_primers_bp)," / ",number(row$product_excluding_primers_bp)),
      paste0(escape(row$reference_accession),"<br><small>",escape(row$ncbi_snapshot_id),"</small>"),
      paste0(escape(row$taxon),"<br><small>",escape(row$taxonomic_rank)," · ",escape(row$parent_path),"</small>"),
      if(is.na(row$unite_exact_products)) "Unavailable" else paste0(number(row$unite_exact_products)," / ",number(row$unite_selected_references)),
      paste0("<span data-evidence=\"",escape(row$combination_status),"\">",escape(.pub_evidence_label(row$combination_status)),"</span><br><small>",if(row$included_in_figure) "Shown in figure" else escape(.pub_reason_label(row$exclusion_reason)),"</small>"),
      paste0("[",escape(row$forward_citation_id),"] [",escape(row$reverse_citation_id),"] [",escape(row$coordinate_citation_id),"]",if(nzchar(row$unite_citation_id)) paste0(" [",escape(row$unite_citation_id),"]") else ""))
    paste0("<tr>",paste0("<td>",cells,"</td>",collapse=""),"</tr>")
  },character(1))
  references<-vapply(seq_len(nrow(export$citations)),function(i) {
    row<-export$citations[i,,drop=FALSE]
    paste0("<li id=\"",escape(row$citation_id),"\">[",escape(row$citation_id),"] ",escape(row$citation),
      if(nzchar(row$status)) paste0(" <small>",escape(row$status),"</small>") else "",
      if("bibliography_status" %in% names(row) && row$bibliography_status=="short_citation_only_full_reference_unavailable") " <small>Full bibliography unavailable</small>" else "",
      if(nzchar(row$url)) paste0(" · <a href=\"",escape(row$url),"\">",escape(row$url),"</a>") else "","</li>")
  },character(1))
  paste0('<!doctype html><html lang="en"><head><meta charset="utf-8"><title>',escape(export$context$title),'</title><style>',
    'body{font:14px/1.45 Arial,Helvetica,sans-serif;color:#172033;margin:24px}h1{font-size:22px}table{width:100%;border-collapse:collapse;font-size:12px}th,td{border:1px solid #cbd5e1;padding:8px;vertical-align:top;text-align:left}th{background:#eef2f6}.oligo{font-family:monospace;overflow-wrap:anywhere}small{color:#475569}li{margin:6px 0;overflow-wrap:anywhere}@media print{body{margin:0;font-size:10pt}table{font-size:8pt}thead{display:table-header-group}tr{break-inside:avoid}a{color:inherit;text-decoration:none}@page{size:landscape;margin:12mm}}',
    '</style></head><body><h1>',escape(export$context$title),'</h1><p>',escape(export$caption),'</p><table><thead><tr>',
    paste0('<th>',c('Figure / combination','Forward oligo (5′–3′)','Reverse oligo (5′–3′)','Product bp incl. / excl. primers','Coordinate reference','Selected taxon','UNITE exact / selected','Evidence / figure inclusion','Citations'),'</th>',collapse=''),
    '</tr></thead><tbody>',paste(rows,collapse=''),'</tbody></table><h2>Methods</h2><p>',escape(export$methods),'</p><h2>References</h2><ol>',paste(references,collapse=''),'</ol></body></html>\n')
}

primer_publication_png_available <- function() requireNamespace("rsvg",quietly=TRUE)
primer_publication_write <- function(export,path,format=c("svg","csv","html","caption","methods","citations","landmarks","png"),width=4200L) {
  format<-match.arg(format)
  if(format=="png") {
    if(!primer_publication_png_available()) stop("PNG rendering is unavailable; the standalone SVG remains available.")
    if(length(width)!=1L || !.pub_integer(width) || width<100L || width>12000L) stop("Choose a PNG width from100 to12000 pixels.")
    rsvg::rsvg_png(charToRaw(enc2utf8(export$svg)),file=path,width=as.integer(width))
  } else if(format %in% c("csv","citations","landmarks")) {
    value<-switch(format,csv=export$table,citations=export$citations,landmarks=export$landmarks)
    utils::write.csv(value,path,row.names=FALSE,na="",fileEncoding="UTF-8")
  } else writeLines(enc2utf8(export[[format]]),path,useBytes=TRUE)
  invisible(path)
}

primer_publication_alignment_export <- function(window,region,group,version,source_doi,pair=NULL,title=NULL) {
  if(!exists(".unite_reference_ui_alignment_svg",mode="function")) stop("Source R/unite_reference_ui.R before regional figure export.")
  alignment<-window$alignment
  if(!is.character(alignment) || !length(alignment) || is.null(names(alignment)) || anyDuplicated(names(alignment)) ||
     length(unique(nchar(alignment)))!=1L || any(!grepl("^[ACGTURYSWKMBDHVNI-]+$",alignment))) stop("Regional publication alignment is invalid.")
  if(!nzchar(.pub_text(version)) || !nzchar(.pub_text(group)) || !nzchar(.pub_text(region)) ||
     !isTRUE(window$end-window$start+1L==nchar(alignment[1L]))) stop("Regional window coordinates or source identity differ.")
  title<-.pub_text(title,paste(region,"regional comparison ·",group))
  refs<-data.frame(citation_id="C1",role="UNITE reference dataset",citation=paste("UNITE",version),
    url=if(nzchar(.pub_text(source_doi))) paste0("https://doi.org/",sub("^https://doi.org/","",source_doi)) else "https://unite.ut.ee/repository.php",status="frozen computational regional alignment",stringsAsFactors=FALSE)
  caption<-paste(title,".",length(alignment),"displayed references; alignment columns",window$start,"through",window$end,
    "from frozen version",version,". Colored cells preserve actual resolved or ambiguous bases; dashes are alignment gaps.",
    if(!is.null(pair)) paste("Selected combination context:",.pub_text(pair),".") else "",
    "Columns are gapped regional alignment positions; no amplicon bp or amplification performance is inferred from this panel.")
  methods<-"The exported window is the actual retained regional alignment slice. Computational ITSx boundaries and within-clade alignment scope remain source limitations; missing-boundary records are not rendered as all-gap sequences."
  svg<-.pub_find_svg(.unite_reference_ui_alignment_svg(window,region))
  records<-data.frame(record_id=names(alignment),version_id=version,group_id=group,region=region,start_column=window$start,end_column=window$end,aligned_window_sequence=unname(alignment),stringsAsFactors=FALSE)
  list(svg=.pub_svg_document(svg,title,caption,methods,refs,list(version=version,group=group,region=region,start_column=window$start,end_column=window$end,record_ids=names(alignment))),
    table=records,geometry_table=records,caption=caption,methods=methods,citations=refs,landmarks=data.frame(),context=list(title=title))
}
