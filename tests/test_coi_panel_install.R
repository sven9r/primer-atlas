#!/usr/bin/env Rscript
# Exercise the actual R/Python boundary installer using an already-acquired
# small fixture; no acquisition, release publication or repository data writes.
.libPaths(c(file.path(getwd(), '.Rlib'), .libPaths()))
source('R/data_layer.R')
source('R/coi_retention.R')
repo <- normalizePath(getwd())
worker <- tempfile('coi_install_fixture_')
dir.create(worker)
on.exit(unlink(worker, recursive = TRUE), add = TRUE)
for(folder in c('R','scripts','data/derived/coi_reference_panel','data/provenance')) dir.create(file.path(worker,folder),recursive=TRUE)
for(file in c('R/data_layer.R','R/coi_retention.R','scripts/update_order_reference_panels.R')) stopifnot(file.copy(file.path(repo,file),file.path(worker,file)))
if(dir.exists(file.path(repo,'.Rlib'))) stopifnot(file.symlink(file.path(repo,'.Rlib'),file.path(worker,'.Rlib')))
groups <- c('OrderA','OrderB')
write.csv(data.frame(order=groups,environment_group='fixture'),file.path(worker,'data/orders_22.csv'),row.names=FALSE)
panel <- file.path(worker,'data/derived/coi_reference_panel')
headers <- paste0(c('AA000001.1','AA000002.1','BB000001.1','BB000002.1'),' original header preserved')
sequences <- c('ACGTACGTACGT','ACGTACGTACGA','TTGCTTGCTTGC','TTGCTTGCTTGA')
# Plain character inputs to vapply acquire automatic value names by default.
# Comparing those vectors with unnamed CSV columns used to reject valid hashes.
stopifnot(!is.null(names(vapply(sequences,digest::digest,character(1),algo='sha256',serialize=FALSE))))
records <- data.frame(target_group=rep(groups,each=2),accession=coi_accession(headers),header=headers,
  sequence_sha256=vapply(sequences,digest::digest,character(1),algo='sha256',serialize=FALSE,USE.NAMES=FALSE),
  header_sha256=vapply(headers,digest::digest,character(1),algo='sha256',serialize=FALSE,USE.NAMES=FALSE),
  first_release_id='fixture-release',added_at_utc='2026-10-11T00:00:00Z',stringsAsFactors=FALSE)
write_fasta <- function(headers,sequences,path) {
 dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
 writeLines(unlist(lapply(seq_along(headers),function(i)c(paste0('>',headers[i]),sequences[i])),use.names=FALSE),path)
}
write_fasta(headers,sequences,file.path(panel,'panel.fasta'))
for(name in c('panel_manifest.csv','ledger.csv','taxonomy.csv')) write.csv(records,file.path(panel,name),row.names=FALSE)
write.csv(data.frame(target_group=groups,retained_after=2,target_met=TRUE),file.path(panel,'growth_qa.csv'),row.names=FALSE)
jsonlite::write_json(list(marker_id='COI',retained_count=4),file.path(panel,'run_manifest.json'),auto_unbox=TRUE)
# Existing group inputs have one frozen accession each; installing the panel
# must preserve their original headers and bases while adding the second row.
for(i in seq_along(groups)) {
 row <- (i-1L)*2L+1L
 write_fasta(headers[row],sequences[row],file.path(worker,'data/primerminer/raw',groups[i],paste0(groups[i],'_GB.fasta')))
}
old_dir <- getwd(); previous <- Sys.getenv(c('COI_USE_EXISTING_PANEL','COI_PANEL_DIR','COI_PREVIOUS_PANEL_DIR'),unset=NA_character_)
on.exit({setwd(old_dir);for(name in names(previous))if(is.na(previous[[name]]))Sys.unsetenv(name)else do.call(Sys.setenv,setNames(list(previous[[name]]),name))},add=TRUE)
setwd(worker)
Sys.setenv(COI_USE_EXISTING_PANEL='true',COI_PANEL_DIR='data/derived/coi_reference_panel',COI_PREVIOUS_PANEL_DIR='data/derived/coi_reference_panel')
run <- function(report_failure=TRUE){
 log <- file.path(worker,'installer.log')
 status <- system2(file.path(R.home('bin'),'Rscript'),'scripts/update_order_reference_panels.R',stdout=log,stderr=log)
 if(status!=0L && report_failure)cat(paste(readLines(log,warn=FALSE),collapse='\n'),'\n')
 status
}
stopifnot(run()==0L)
installed <- coi_validate_retained_state(worker,groups)
stopifnot(nrow(installed$records)==4L,identical(installed$records$header,headers),
          identical(installed$records$sequence_sha256,records$sequence_sha256),
          identical(atlas_sha256(file.path(worker,'data/provenance/coi_reference_accession_ledger.csv')),atlas_sha256(file.path(panel,'ledger.csv'))))
raw_paths <- file.path(worker,'data/primerminer/raw',groups,paste0(groups,'_GB.fasta'))
before <- vapply(raw_paths,atlas_sha256,character(1))
# A changed header without updated provenance fails before any group is replaced.
corrupt_headers <- headers;corrupt_headers[1] <- paste0(corrupt_headers[1],' changed')
write_fasta(corrupt_headers,sequences,file.path(panel,'panel.fasta'))
stopifnot(run(FALSE)!=0L,identical(before,vapply(raw_paths,atlas_sha256,character(1))))
message('COI wrapper installs a valid acquired panel despite automatic hash-vector names and rejects changed header provenance before mutation.')
