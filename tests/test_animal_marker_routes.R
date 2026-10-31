# Every exposed animal-group/marker route has real, correctly scoped examples
# and source-backed components. This is a usability/identity gate, not a claim
# of broad taxonomic sampling or experimental amplification performance.
source("R/marker_targets.R")
source("R/primer_mixer.R")
source("R/animal_fit.R")
markers <- read.csv("data/catalog/markers.csv",stringsAsFactors=FALSE)
links <- read.csv("data/catalog/marker_organism_groups.csv",stringsAsFactors=FALSE)
catalog <- primer_mixer_catalog()
routes <- 0L
for(group in c("METAZOA","ARTHROPODS","FISH","BIRDS","ZOOPLANKTON")) {
  options <- atlas_target_marker_options(group,markers,links)
  for(marker in options$marker_id) {
    components <- catalog[catalog$marker_id==marker,,drop=FALSE]
    forward <- components$primer_id[components$direction=="forward"]
    reverse <- components$primer_id[components$direction=="reverse"]
    templates <- primer_mixer_templates(marker_id=marker,target_group=group)
    stopifnot(length(forward)>0L,length(reverse)>0L,nrow(templates)>0L,
      all(templates$marker_id==marker),all(nzchar(templates$sequence)),
      all(nzchar(templates$taxonomy_version)),all(nzchar(templates$taxonomy_source_sha256)))
    if(!group %in% c("METAZOA")) stopifnot(all(templates$target_group==group))
    result <- primer_mixer_build(components,forward[1L],reverse[1L],templates,marker,max_mismatches=2L)
    fit <- animal_fit_data(result,rank="family")
    stopifnot(nrow(result$products)==nrow(templates),nrow(fit$details)==nrow(templates),
      sum(fit$summary$n_panel)==nrow(templates),
      all(fit$summary$n_exact+fit$summary$n_substitution+fit$summary$n_partial+
        fit$summary$n_ambiguous+fit$summary$n_missing==fit$summary$n_panel))
    routes <- routes+1L
  }
}
stopifnot(routes==22L)
message("All ",routes," animal-group/marker routes have scoped source templates, oligos and reconciled taxonomic fit.")
