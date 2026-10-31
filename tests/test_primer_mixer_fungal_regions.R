source("R/data_layer.R")
source("R/fungal_versions.R")
source("R/primer_mixer.R")

snapshot <- fungal_read_snapshot("v3-c6fa87f256e1",remote=FALSE)
original_sequences <- snapshot$sequences
ssu <- primer_mixer_fungal_templates(snapshot,"18S")
lsu <- primer_mixer_fungal_templates(snapshot,"28S")
stopifnot(nrow(ssu)==3L,nrow(lsu)==4L,sum(ssu$partial)==2L,sum(lsu$partial)==3L,
  sum(ssu$boundary_status=="complete")==1L,sum(lsu$boundary_status=="complete")==1L,
  !any(ssu$projected),!any(lsu$projected),all(ssu$marker_id=="18S"),all(lsu$marker_id=="28S"),
  all(ssu$version==snapshot$manifest$id),all(lsu$target_group=="FUNGI"),
  identical(snapshot$sequences,original_sequences))
for(rows in list(ssu,lsu)) for(i in seq_len(nrow(rows))) {
  parent <- unname(snapshot$sequences[[rows$accession_version[i]]])
  interval <- snapshot$landmarks[snapshot$landmarks$accession==rows$accession_version[i] &
    snapshot$landmarks$landmark_id==rows$source_landmark_id[i],,drop=FALSE]
  stopifnot(nrow(interval)==1L,rows$source_start[i]==interval$start,rows$source_end[i]==interval$end,
    identical(rows$sequence[i],substr(parent,rows$source_start[i],rows$source_end[i])),
    nchar(rows$sequence[i])==rows$source_end[i]-rows$source_start[i]+1L,
    rows$parent_sequence_sha256[i]==digest::digest(parent,algo="sha256",serialize=FALSE),
    rows$sequence_sha256[i]==digest::digest(rows$sequence[i],algo="sha256",serialize=FALSE))
  tax <- snapshot$taxonomy[snapshot$taxonomy$accession==rows$accession_version[i],,drop=FALSE]
  stopifnot(rows$taxonomy_source_sha256[i]==tax$taxonomy_source_sha256,
    rows$phylum[i]==tax$phylum,rows$class[i]==tax$class,rows$species[i]=="Unresolved")
}
pieces <- lsu[lsu$accession_version=="PZ276773.1",,drop=FALSE]
stopifnot(nrow(pieces)==2L,all(pieces$partial),identical(pieces$source_start,c(1L,8127L)),
  identical(pieces$source_end,c(1182L,10646L)),!anyDuplicated(pieces$template_id))

# Invalid/withheld intervals cannot become a full-gene surrogate. Projection
# provenance is explicit when a bounded gene interval is projected.
projected <- snapshot
row <- which(projected$landmarks$label=="18S" & projected$landmarks$accession=="PZ276773.1")
projected$landmarks$evidence_status[row] <- "projected_gene_interval"
projected$landmarks$mapping_method[row] <- "source_reference_projection"
projected_ssu <- primer_mixer_fungal_templates(projected,"18S")
stopifnot(sum(projected_ssu$projected)==1L,nrow(projected_ssu)==3L)
withheld <- snapshot
withheld$landmarks$boundary_status[row] <- "withheld"
stopifnot(nrow(primer_mixer_fungal_templates(withheld,"18S"))==2L)
invalid <- snapshot
invalid$landmarks$end[row] <- nchar(snapshot$sequences[["PZ276773.1"]])+1L
stopifnot(nrow(primer_mixer_fungal_templates(invalid,"18S"))==2L)
v1 <- fungal_read_snapshot("v1-9e958f105c46",remote=FALSE)
v1_ssu <- primer_mixer_fungal_templates(v1,"18S")
v1_lsu <- primer_mixer_fungal_templates(v1,"28S")
stopifnot(nrow(v1_ssu)==0L,nrow(v1_lsu)==0L,
  grepl("pilot-display anchors are withheld",attr(v1_ssu,"unavailable_reason"),fixed=TRUE))

catalog <- primer_mixer_catalog()
for(marker in c("18S","28S")) {
  aliases <- catalog[catalog$marker_id==marker & catalog$catalog_origin=="UNITE_locus_alias",,drop=FALSE]
  stopifnot(nrow(aliases)==if(marker=="18S") 30L else 32L,all(aliases$target_group=="FUNGI"),
    all(aliases$original_marker_id=="ITS_FUNGAL"),all(grepl("newly assembled pairings remain candidates",aliases$evidence_note,fixed=TRUE)))
  parents <- catalog[match(aliases$original_primer_id,catalog$primer_id),,drop=FALSE]
  for(field in c("primer_name","sequence","source_sequence","direction","citation","source_url","gene_locus"))
    stopifnot(identical(aliases[[field]],parents[[field]]))
  stopifnot(!any(aliases$primer_id %in% parents$primer_id),all(parents$marker_id=="ITS_FUNGAL"))
  again <- primer_mixer_catalog()
  stopifnot(all(aliases$primer_id %in% again$primer_id))
  f <- head(aliases$primer_id[aliases$direction=="forward"],2L)
  r <- head(aliases$primer_id[aliases$direction=="reverse"],2L)
  templates <- if(marker=="18S") ssu else lsu
  result <- primer_mixer_build(catalog,f,r,templates,marker,max_mismatches=2L)
  stopifnot(nrow(result$products)==length(f)*length(r)*nrow(templates),all(result$summary$n_templates==nrow(templates)),
    !any(result$products$status=="marker_mismatch"))
}
message("Fungal marker extraction passed:3SSU/4LSU bounded intervals,partial/projected provenance,source taxonomy,unchanged sourceDNA and30/32 source-preserving UNITE aliases.")
