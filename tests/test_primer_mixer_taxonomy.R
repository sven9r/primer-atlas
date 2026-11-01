source("R/primer_mixer.R")
source("R/primer_mixer_exports.R")

all <- primer_mixer_templates(marker_id="COI",target_group="ALL")
fish <- primer_mixer_templates(marker_id="COI",target_group="FISH")
birds <- primer_mixer_templates(marker_id="COI",target_group="BIRDS")
arthropods <- primer_mixer_templates(marker_id="COI",target_group="ARTHROPODS")
zooplankton <- primer_mixer_templates(marker_id="COI",target_group="ZOOPLANKTON")
ranks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")
stopifnot(nrow(all)==6L,nrow(fish)==3L,nrow(birds)==1L,nrow(arthropods)==1L,nrow(zooplankton)==1L,
  all(ranks %in% names(all)),!anyDuplicated(all$template_id),
  all(all$taxonomy_status=="source_taxid_verified"),all(nzchar(all$taxonomy_source_url)),
  all(grepl("^[a-f0-9]{64}$",all$taxonomy_source_sha256)),
  all(fish$target_group=="FISH"),all(birds$target_group=="BIRDS"),all(arthropods$target_group=="ARTHROPODS"),
  all(arthropods$taxid==7245L),all(arthropods$family=="Drosophilidae"),
  all(birds$taxid==9031L),all(birds$class=="Aves"),all(fish$phylum=="Chordata"),
  all(zooplankton$target_group=="ZOOPLANKTON"),all(zooplankton$family=="Daphniidae"),
  all(zooplankton$order=="Diplostraca"),all(zooplankton$class=="Branchiopoda"),
  !any(fish$template_id %in% arthropods$template_id))

# Taxonomy is joined onto source records without changing sequence bytes,
# source coordinates, source URL or accession/version provenance.
source_rows <- utils::read.csv("data/pinned/FISH_STARTER/v1/templates.csv",stringsAsFactors=FALSE,check.names=FALSE)
source_taxa <- utils::read.csv("data/pinned/FISH_STARTER/v1/taxonomy.csv",stringsAsFactors=FALSE,check.names=FALSE)
joined <- primer_mixer_templates(target_group="FISH")
si <- match(joined$template_id,source_rows$template_id)
ti <- match(joined$template_id,source_taxa$template_id)
stopifnot(nrow(joined)==13L,nrow(source_rows)==13L,identical(joined$sequence,source_rows$sequence[si]),
  identical(joined$source_url,source_rows$source_url[si]),identical(joined$version,source_rows$version[si]),
  identical(joined$source_start,source_rows$source_start[si]),
  identical(joined$family,source_taxa$family[ti]),identical(joined$order,source_taxa$order[ti]),
  identical(joined$taxonomy_source_url,source_taxa$taxonomy_source_url[ti]))

# Denominators are the explicit selected templates, not all source taxonomy
# rows or an inferred taxon inventory. Every template remains represented.
catalog <- primer_mixer_catalog()
f <- catalog$primer_id[catalog$primer_name=="FishF1" & catalog$direction=="forward"][1L]
r <- catalog$primer_id[catalog$primer_name=="FishR1" & catalog$direction=="reverse"][1L]
res <- primer_mixer_build(catalog,f,r,all,"COI",max_mismatches=2L)
stopifnot(nrow(res$products)==6L,res$summary$n_templates==6L,
  res$summary$n_matching_marker_templates==6L,sum(table(res$templates$family))==6L,
  all(res$templates$family==all$family))
res_fish <- primer_mixer_build(catalog,f,r,fish,"COI",max_mismatches=2L)
stopifnot(nrow(res_fish$products)==3L,res_fish$summary$n_templates==3L,
  sum(table(res_fish$templates$family))==3L,!any(res_fish$templates$family %in% c("Drosophilidae","Phasianidae")))
table <- primer_mixer_results_table(res)
stopifnot(identical(table$template_family,all$family[match(table$template_id,all$template_id)]),
  all(table$n_input_templates==6L),all(table$template_taxonomy_status=="source_taxid_verified"))

# Each ecological-group marker uses only its own actual frozen source record.
# The test makes no assertion that a chosen published oligo must match it.
for(marker in c("COI","18S","28S","16S_MT")) {
  own <- primer_mixer_templates(marker_id=marker,target_group="ZOOPLANKTON")
  stopifnot(nrow(own)==1L,all(own$marker_id==marker),all(own$target_group=="ZOOPLANKTON"),
    all(own$family=="Daphniidae"),all(own$taxonomy_status=="source_taxid_verified"),
    all(nzchar(own$ecology_source_url)),!any(own$template_id %in% joined$template_id))
  marker_catalog <- catalog[catalog$marker_id==marker,,drop=FALSE]
  fi <- head(marker_catalog$primer_id[marker_catalog$direction=="forward"],1L)
  ri <- head(marker_catalog$primer_id[marker_catalog$direction=="reverse"],1L)
  stopifnot(length(fi)==1L,length(ri)==1L)
  marker_result <- primer_mixer_build(catalog,fi,ri,own,marker,max_mismatches=2L)
  stopifnot(nrow(marker_result$pairs)==1L,nrow(marker_result$products)==1L,
    marker_result$summary$n_templates==1L,marker_result$summary$n_matching_marker_templates==1L,
    marker_result$products$status!="marker_mismatch")
  if(marker=="16S_MT") {
    exported_marker <- primer_mixer_results_table(marker_result)
    stopifnot(own$source_strand=="-",grepl("complement",own$source_feature_location,fixed=TRUE),
      exported_marker$template_source_strand=="-",identical(exported_marker$template_source_feature_location,own$source_feature_location))
  }
}
fish_16s <- primer_mixer_templates(marker_id="16S_MT",target_group="FISH")
arthropod_16s <- primer_mixer_templates(marker_id="16S_MT",target_group="ARTHROPODS")
all_16s <- primer_mixer_templates(marker_id="16S_MT",target_group="ALL")
stopifnot(nrow(fish_16s)==3L,nrow(arthropod_16s)==1L,nrow(all_16s)==5L,
  all(fish_16s$target_group=="FISH"),all(arthropod_16s$target_group=="ARTHROPODS"),
  arthropod_16s$template_id=="ARTHROPODS_NC_001322.1_16S_MT",arthropod_16s$source_strand=="-",
  arthropod_16s$taxid==7245L,arthropod_16s$family=="Drosophilidae",
  !any(fish_16s$template_id %in% arthropod_16s$template_id),nrow(primer_mixer_templates(target_group="ALL"))==25L)

# User-supplied headers that resemble accessions cannot inherit source-backed
# taxonomy. They remain explicitly user-provided sequences in the engine.
fasta <- primer_mixer_parse_fasta(paste0(">NC_001322.1_COX1\n",arthropods$sequence),"COI")
stopifnot(!any(ranks %in% names(fasta)),fasta$organism=="User-provided template")
unverified <- primer_mixer_build(catalog,f,r,fasta,"COI")
stopifnot(!"family" %in% names(unverified$templates),
  primer_mixer_results_table(unverified)$template_family=="Unresolved")

# Redundant marker/group identities are checked before accepting a taxonomy
# supplement, even if its template ID matches a selected source record.
fixture_root <- tempfile("mixer-taxonomy-")
fixture_path <- file.path(fixture_root,"data/pinned/ANIMAL_STARTER/v1")
dir.create(fixture_path,recursive=TRUE)
original_panel <- "data/pinned/ANIMAL_STARTER/v1"
invisible(file.copy(file.path(original_panel,c("manifest.json","templates.csv","taxonomy.csv")),fixture_path,overwrite=TRUE))
fixture <- utils::read.csv(file.path(fixture_path,"taxonomy.csv"),stringsAsFactors=FALSE,check.names=FALSE)
fixture$marker_id[fixture$template_id==arthropods$template_id] <- "12S_MT"
utils::write.csv(fixture,file.path(fixture_path,"taxonomy.csv"),row.names=FALSE)
fixture_manifest <- jsonlite::read_json(file.path(fixture_path,"manifest.json"),simplifyVector=FALSE)
fixture_manifest$artifacts[["taxonomy.csv"]]$sha256 <- digest::digest(file=file.path(fixture_path,"taxonomy.csv"),algo="sha256",serialize=FALSE)
fixture_manifest$artifacts[["taxonomy.csv"]]$bytes <- file.info(file.path(fixture_path,"taxonomy.csv"))$size
jsonlite::write_json(fixture_manifest,file.path(fixture_path,"manifest.json"),auto_unbox=TRUE,pretty=TRUE)
error <- tryCatch({.primer_mixer_template_taxonomy(arthropods,fixture_root);NULL},error=function(e)e)
stopifnot(inherits(error,"error"),grepl("disagrees with raw template identity",conditionMessage(error)))
unlink(fixture_root,recursive=TRUE)

# A corrupted or missing frozen artifact cannot silently downgrade to an
# unresolved taxonomy table, another panel, or a cached previous file read.
tamper_root <- tempfile("mixer-starter-integrity-")
tamper_panel <- file.path(tamper_root,"data/pinned/FISH_STARTER/v1")
dir.create(tamper_panel,recursive=TRUE)
invisible(file.copy(file.path("data/pinned/FISH_STARTER/v1",c("manifest.json","templates.csv","taxonomy.csv")),tamper_panel,overwrite=TRUE))
stopifnot(nrow(primer_mixer_templates(tamper_root,"COI","FISH"))==3L)
expect_integrity_error <- function(code,pattern) {
  error <- tryCatch({force(code);NULL},error=function(e)e)
  stopifnot(inherits(error,"error"),grepl("Starter reference integrity failed",conditionMessage(error)),
    grepl(pattern,conditionMessage(error)),grepl("No substitute reference was loaded",conditionMessage(error),fixed=TRUE))
}
tampered <- utils::read.csv(file.path(tamper_panel,"templates.csv"),stringsAsFactors=FALSE,check.names=FALSE)
substr(tampered$sequence[1L],1L,1L) <- ifelse(substr(tampered$sequence[1L],1L,1L)=="A","C","A")
utils::write.csv(tampered,file.path(tamper_panel,"templates.csv"),row.names=FALSE)
expect_integrity_error(primer_mixer_templates(tamper_root,"COI","FISH"),"SHA-256 mismatch")
invisible(file.copy("data/pinned/FISH_STARTER/v1/templates.csv",file.path(tamper_panel,"templates.csv"),overwrite=TRUE))
stopifnot(nrow(primer_mixer_templates(tamper_root,"COI","FISH"))==3L)
tampered <- utils::read.csv(file.path(tamper_panel,"taxonomy.csv"),stringsAsFactors=FALSE,check.names=FALSE)
tampered$family[1L] <- "FabricatedFamily"
utils::write.csv(tampered,file.path(tamper_panel,"taxonomy.csv"),row.names=FALSE)
expect_integrity_error(primer_mixer_templates(tamper_root,"COI","FISH"),"SHA-256 mismatch")
invisible(file.copy("data/pinned/FISH_STARTER/v1/taxonomy.csv",file.path(tamper_panel,"taxonomy.csv"),overwrite=TRUE))
stopifnot(nrow(primer_mixer_templates(tamper_root,"COI","FISH"))==3L)
unlink(file.path(tamper_panel,"taxonomy.csv"))
expect_integrity_error(primer_mixer_templates(tamper_root,"COI","FISH"),"missing required artifact")
invisible(file.copy("data/pinned/FISH_STARTER/v1/taxonomy.csv",file.path(tamper_panel,"taxonomy.csv"),overwrite=TRUE))
unlink(file.path(tamper_panel,"manifest.json"))
expect_integrity_error(primer_mixer_templates(tamper_root,"COI","FISH"),"missing frozen manifest")
unlink(tamper_root,recursive=TRUE)
message("Primer mixer taxonomy tests passed:authoritative rank joins,source preservation,scoped denominators,marker isolation and unclassified user FASTA.")
