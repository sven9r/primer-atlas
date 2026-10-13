#!/usr/bin/env Rscript
.libPaths(c(file.path(getwd(),".Rlib"),.libPaths()))
source("R/functions.R");source("R/data_layer.R");source("R/fungal_explorer.R")
source("R/fungal_versions.R");source("R/fungal_primer_overview.R")
source("R/unite_reference.R");source("R/fungal_combination_builder.R")

snapshot<-fungal_read_snapshot("v3-c6fa87f256e1",remote=FALSE)
before<-serialize(snapshot,NULL)
catalog<-fungal_combination_catalog(snapshot)
stopifnot(nrow(catalog)==112L,sum(catalog$direction=="forward")==40L,sum(catalog$direction=="reverse")==72L,
  all(catalog$component_origin=="frozen_snapshot_catalog"),all(catalog$component_source_version==snapshot$manifest$id),
  catalog$direction[catalog$primer_name=="ITS2"]=="reverse",!anyDuplicated(catalog$selection_label))
f<-catalog$catalog_id[match(c("ITS1","ITS1F","ITS5","ITS3","fITS7"),catalog$primer_name)]
r<-catalog$catalog_id[match(c("ITS2","ITS4","ITS4B","LR21","LR3"),catalog$primer_name)]
result<-fungal_combination_builder(snapshot,f,r,"PZ962702.1")
stopifnot(nrow(result$selected_pairs)==25L,nrow(result$pair_sequences)==25L*21L,
  result$summary$n_selected_combinations==25L,result$summary$n_forward==5L,result$summary$n_reverse==5L,
  result$summary$n_mapped_after_filters==18L,result$summary$n_unmapped_on_reference==7L,
  result$summary$n_frozen_pair_identities==19L,result$summary$n_new_candidate_identities==6L,
  nrow(result$map_pair_rows)==18L,nrow(result$unmapped_pairs)==7L,
  !any(result$selected_pairs$pair_origin=="new_selected_cartesian_candidate" & result$selected_pairs$evidence_status=="documented_pair"))
existing<-result$selected_pairs[result$selected_pairs$pair_origin=="frozen_pair_identity",]
index<-match(.fungal_combo_pair_keys(existing),.fungal_combo_pair_keys(snapshot$pairs))
stopifnot(identical(existing$pair_id,snapshot$pairs$pair_id[index]),
  identical(existing$evidence_status,snapshot$pairs$evidence_status[index]),
  all(existing$source_citation==snapshot$pairs$source_citation[index]))
map<-fungal_primer_overview_data(result$working_snapshot,result$map_pair_rows,"PZ962702.1",max_rows=25L)
short<-map$pair_meta[map$pair_meta$pair_label=="ITS1 + ITS2",]
stopifnot(nrow(map$pair_meta)==18L,short$aligned_amplicon_bp==335L,short$targeted_region_bp==296L)

# Missing crosses still have stable exact-component identities and complete
# per-template states. They are not replaced by only the combinations that work.
one<-fungal_combination_builder(snapshot,f[1L],r[1:3],"PZ962702.1")
stopifnot(one$summary$n_selected_combinations==3L,nrow(one$selected_pairs)==3L,nrow(one$map_pair_rows)==2L,
  nrow(one$unmapped_pairs)==1L,one$unmapped_pairs$pair_label=="ITS1 + ITS4B",
  startsWith(one$unmapped_pairs$pair_id,"FUNGAL_COMB_"),one$unmapped_pairs$n_scorable==0L,
  one$unmapped_pairs$reference_status=="reverse_site_missing")
repeat_one<-fungal_combination_builder(snapshot,f[1L],r[1:3],"FN812768.2")
stopifnot(identical(one$selected_pairs$pair_id,repeat_one$selected_pairs$pair_id),
  identical(one$selected_pairs$pair_id,fungal_combination_builder(snapshot,f[1L],r[1:3],"PZ962702.1")$selected_pairs$pair_id))

# Filtering changes the figure selection while retaining ALL requested crosses
# for publication tables and an explicit reason for every omitted figure row.
filtered<-fungal_combination_builder(snapshot,f,r,"PZ962702.1",platform="Illumina",length_range=c(300,400))
stopifnot(nrow(filtered$selected_pairs)==25L,
  all(filtered$map_pair_rows$reference_amplicon_bp>=300 & filtered$map_pair_rows$reference_amplicon_bp<=400),
  all(filtered$map_pair_rows$included_in_map),all(nzchar(filtered$filtered_out_pairs$filter_reason)),
  nrow(filtered$map_pair_rows)+nrow(filtered$filtered_out_pairs)+nrow(filtered$unmapped_pairs)==25L)
doc<-fungal_combination_builder(snapshot,f,r,"PZ962702.1",evidence_filter="documented_pair")
stopifnot(nrow(doc$selected_pairs)==25L,all(doc$map_pair_rows$evidence_status=="documented_pair"),
  any(grepl("evidence_filter",doc$selected_pairs$filter_reason)))
ssu<-fungal_combination_builder(snapshot,f,r,"PZ962702.1",region="SSU")
stopifnot(!"ITS3 + ITS4" %in% ssu$map_pair_rows$pair_label,
  any(grepl("region_filter",ssu$selected_pairs$filter_reason)),nrow(ssu$selected_pairs)==25L)

# Selectors may be empty, or a cohort may contain no references. Neither state
# invents products or loses the definition of the requested Cartesian set.
empty<-fungal_combination_builder(snapshot,character(),r,"PZ962702.1")
none<-fungal_combination_builder(snapshot,f,r,accessions=character())
stopifnot(empty$summary$n_selected_combinations==0L,nrow(empty$selected_pairs)==0L,nrow(empty$map_pair_rows)==0L,
  none$summary$n_selected_combinations==25L,none$summary$n_panel==0L,nrow(none$pair_sequences)==0L,
  nrow(none$unmapped_pairs)==25L,all(none$selected_pairs$reference_status=="no_template_selected"))
stopifnot(inherits(try(fungal_combination_builder(snapshot,r[1],f[1],"PZ962702.1"),silent=TRUE),"try-error"),
  inherits(try(fungal_combination_builder(snapshot,c(f,catalog$catalog_id[catalog$primer_name=="gITS7"]),r,"PZ962702.1"),silent=TRUE),"try-error"))

# A synthetic panel supports all25 combinations, proving that the 20-example
# cap is not inherited. A second template has one missing forward and an
# ambiguous reverse site: products remain on their own raw template only.
set.seed(2718)
dna<-function(n) paste(sample(c("A","C","G","T"),n,replace=TRUE),collapse="")
fseq<-replicate(5L,dna(20L));rseq<-replicate(5L,dna(20L));template<-dna(900L)
for(i in 1:5) {
  substr(template,1L+(i-1L)*25L,20L+(i-1L)*25L)<-fseq[i]
  substr(template,600L+(i-1L)*25L,619L+(i-1L)*25L)<-reverse_complement(rseq[i])
}
second<-template;substr(second,1L,20L)<-paste(rep("N",20L),collapse="")
substr(second,800L,819L)<-reverse_complement(rseq[5L])
synthetic_catalog<-data.frame(catalog_id=paste0("SYN_",c(paste0("F",1:5),paste0("R",1:5))),
  primer_name=c(paste0("SYN_F",1:5),paste0("SYN_R",1:5)),sequence=c(fseq,rseq),direction=rep(c("forward","reverse"),each=5L),
  gene_locus="ITS",target="Synthetic test",reference="Synthetic fixture; no biological evidence",primary_reference_status="publication_cited",
  primary_reference_url="https://example.invalid/synthetic",stringsAsFactors=FALSE)
taxonomy<-data.frame(accession=c("A","B"),class="SyntheticClass",phylum="SyntheticPhylum",organism="Synthetic test",stringsAsFactors=FALSE)
synthetic<-c(list(manifest=list(id="synthetic-test-only",marker_id="ITS_FUNGAL"),catalog=synthetic_catalog,sequences=c(A=template,B=second),
  taxonomy=taxonomy,landmarks=data.frame(landmark_id=character(),label=character(),start=integer(),end=integer())),
  fungal_evaluate_panel(synthetic_catalog,c(A=template,B=second),taxonomy))
# list() concatenation above deliberately shares the frozen engine outputs,
# with a single unambiguous copy of every snapshot field.
synthetic$catalog<-synthetic_catalog
all25<-fungal_combination_builder(synthetic,synthetic_catalog$catalog_id[1:5],synthetic_catalog$catalog_id[6:10],"A")
stopifnot(nrow(all25$map_pair_rows)==25L,nrow(all25$selected_pairs)==25L,nrow(all25$pair_sequences)==50L,
  all(all25$pair_sequences$amplicon_bp[all25$pair_sequences$status=="scorable"]>=100L),
  all(all25$pair_sequences$class=="SyntheticClass"))
partial<-fungal_combination_builder(synthetic,synthetic_catalog$catalog_id[1:5],synthetic_catalog$catalog_id[6:10],"B")
stopifnot(nrow(partial$selected_pairs)==25L,nrow(partial$map_pair_rows)==16L,nrow(partial$unmapped_pairs)==9L,
  sum(partial$selected_pairs$reference_status=="ambiguous_binding_site")==5L)

# An explicitly supplied oligo receives its own source provenance and fresh
# exact search. It cannot rewrite a historical component or become documented.
extra<-synthetic_catalog[1L,];extra$catalog_id<-"EXPLICIT_OTHER_F";extra$primer_name<-"Explicit other forward"
extra$sequence<-substr(template,200L,219L);extra$source_version<-"explicit-test-v1"
new<-fungal_combination_builder(synthetic,extra$catalog_id,synthetic_catalog$catalog_id[6L],"A",catalog=extra)
stopifnot(nrow(new$selected_pairs)==1L,nrow(new$map_pair_rows)==1L,
  new$selected_pairs$pair_origin=="new_selected_cartesian_candidate",
  new$selected_pairs$evidence_status!="documented_pair",
  new$components$component_source_version[new$components$catalog_id==extra$catalog_id]=="explicit-test-v1",
  any(new$sites$site_evidence_origin=="new_exact_search_explicit_supplemental_component"))
bounded<-fungal_combination_builder(synthetic,extra$catalog_id,synthetic_catalog$catalog_id[6L],"A",catalog=extra,max_new_site_searches=0L)
stopifnot(nrow(bounded$selected_pairs)==1L,nrow(bounded$map_pair_rows)==0L,
  bounded$selected_pairs$reference_status=="binding_evidence_unavailable",bounded$summary$n_unavailable_site_searches==2L)
changed<-synthetic_catalog[1L,];changed$sequence<-dna(20L)
stopifnot(inherits(try(fungal_combination_catalog(synthetic,catalog=changed),silent=TRUE),"try-error"))
changed<-synthetic_catalog[1L,];changed$reference<-"A later citation must not overwrite frozen metadata"
preserved<-fungal_combination_catalog(synthetic,catalog=changed)
stopifnot(preserved$reference[preserved$catalog_id==changed$catalog_id]==synthetic_catalog$reference[1L])

# Existing geometry and sites are separately checked. Mismatched cached bases
# or a corrupted frozen product cannot be used as a new verified placement.
bad<-snapshot;bad$sequences[["PZ962702.1"]]<-paste(rep("N",nchar(bad$sequences[["PZ962702.1"]])),collapse="")
stopifnot(inherits(try(fungal_combination_builder(bad,f,r,"PZ962702.1"),silent=TRUE),"try-error"))
bad<-snapshot;index<-which(bad$pair_sequences$pair_id=="ITS1_ITS4" & bad$pair_sequences$accession=="PZ962702.1")
bad$pair_sequences$amplicon_bp[index]<-999L
stopifnot(inherits(try(fungal_combination_builder(bad,f,r,"PZ962702.1"),silent=TRUE),"try-error"))

# Historical V1 can generate missing crosses only from its own frozen catalog
# and sites, while retaining its existing documented identities and source id.
v1<-fungal_read_snapshot("v1-9e958f105c46",remote=FALSE)
history<-fungal_combination_builder(v1,f[1L],r[1:3],"FN812768.2")
stopifnot(nrow(history$selected_pairs)==3L,history$summary$n_panel==1L,
  all(history$selected_pairs$parent_version_id==v1$manifest$id),
  history$selected_pairs$pair_id[history$selected_pairs$pair_label=="ITS1 + ITS4"]=="ITS1_ITS4",
  identical(serialize(snapshot,NULL),before))
defaults<-fungal_combination_defaults(snapshot,"PZ962702.1")
stopifnot(length(defaults$forward_ids)==5L,length(defaults$reverse_ids)==5L,
  all(defaults$forward_ids %in% catalog$catalog_id[catalog$direction=="forward"]),
  all(defaults$reverse_ids %in% catalog$catalog_id[catalog$direction=="reverse"]))

# Sparse absence is evidence only when the exhaustive proof covers the complete
# original cohort. A selected subset never changes that proof's denominator.
broad<-unite_reference_read()
stopifnot(!is.null(broad))
subset_ids<-head(broad$metadata$record_id,12L)
broad_cross<-fungal_combination_builder(broad,f[1L],r[3L],subset_ids[1L],accessions=subset_ids)
stopifnot(nrow(broad_cross$selected_pairs)==1L,broad_cross$selected_pairs$n_panel==12L,
  nrow(broad_cross$pair_sequences)==12L,
  all(broad_cross$pair_sequences$source_accession==broad$metadata$accession[match(subset_ids,broad$metadata$record_id)]),
  any(broad_cross$sites$site_evidence_origin=="verified_exhaustive_sparse_search_no_hit"))
broken_proof<-broad;broken_proof$search_manifest$complete<-FALSE
stopifnot(inherits(try(fungal_combination_builder(broken_proof,f[1L],r[3L],subset_ids[1L],accessions=subset_ids),silent=TRUE),"try-error"))
message("Fungal Cartesian25/25, missing crosses, filters, source identity, historical and sparse-proof fixtures passed.")
