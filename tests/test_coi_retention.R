#!/usr/bin/env Rscript
.libPaths(c(file.path(getwd(), '.Rlib'), .libPaths()))
source('R/data_layer.R')
source('R/coi_retention.R')
repo <- normalizePath(getwd())
fixture <- tempfile('coi_two_workers_')
dir.create(fixture)
on.exit(unlink(fixture, recursive = TRUE), add = TRUE)
groups <- c('Araneae', 'Diptera')
reference <- coi_read_fasta('data/reference/coi_reference_NC_001322.1.fasta')[[1]]
write_fasta <- function(values, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(unlist(lapply(seq_along(values), function(i) c(paste0('>', names(values)[i]), values[i])), use.names = FALSE), path)
}
worker <- function(name) {
  path <- file.path(fixture, name); dir.create(path)
  for (folder in c('R', 'scripts', 'data/reference', 'data/provenance')) dir.create(file.path(path, folder), recursive = TRUE)
  for (file in c('R/functions.R', 'R/data_layer.R', 'R/coi_retention.R', 'scripts/build_22_group_alignments.R', 'data/reference/coi_reference_NC_001322.1.fasta')) {
    stopifnot(file.copy(file.path(repo, file), file.path(path, file)))
  }
  if (dir.exists(file.path(repo, '.Rlib'))) stopifnot(file.symlink(file.path(repo, '.Rlib'), file.path(path, '.Rlib')))
  write.csv(data.frame(order = groups, environment_group = 'test'), file.path(path, 'data/orders_22.csv'), row.names = FALSE)
  path
}
ledger <- function(path, release, previous = NULL) {
  records <- coi_retained_records(path, groups)
  rows <- records[, c('target_group', 'accession', 'sequence_sha256', 'header_sha256')]
  rows$entrez_uid <- paste0('uid_', rows$accession)
  rows$first_release_id <- release
  rows$added_at_utc <- '2026-10-11T00:00:00Z'
  rows$taxonomy_source_sha256 <- paste0('taxonomy_', rows$accession)
  if (!is.null(previous)) {
    matches <- match(paste(rows$target_group, rows$accession), paste(previous$target_group, previous$accession))
    for (column in c('first_release_id', 'added_at_utc', 'taxonomy_source_sha256')) rows[[column]][!is.na(matches)] <- previous[[column]][matches[!is.na(matches)]]
  }
  write.csv(rows, file.path(path, 'data/provenance/coi_reference_accession_ledger.csv'), row.names = FALSE)
  rows
}
align <- function(path) {
  previous_dir <- getwd(); previous_orders <- Sys.getenv('ORDERS', unset = NA_character_)
  previous_mode <- Sys.getenv('EVALUATION_MODE', unset = NA_character_)
  on.exit({setwd(previous_dir); if(is.na(previous_orders)) Sys.unsetenv('ORDERS') else Sys.setenv(ORDERS=previous_orders); if(is.na(previous_mode)) Sys.unsetenv('EVALUATION_MODE') else Sys.setenv(EVALUATION_MODE=previous_mode)}, add=TRUE)
  setwd(path); Sys.setenv(ORDERS=paste(groups, collapse=','), EVALUATION_MODE='full')
  log <- file.path(path, 'build.log')
  status <- system2(file.path(R.home('bin'), 'Rscript'), 'scripts/build_22_group_alignments.R', stdout=log, stderr=log)
  if(status != 0L) stop(paste(readLines(log, warn=FALSE), collapse='\n'))
}
# Worker one creates a validated immutable state using the actual MAFFT builder.
first <- worker('first')
for (g in seq_along(groups)) {
  values <- setNames(c(reference, paste0(substr(reference,1,700), 'A', substr(reference,702,nchar(reference)))),
    paste0('TEST',g, c('A.1 original frozen header','B.1 original frozen header')))
  write_fasta(values, file.path(first, 'data/primerminer/raw', groups[g], paste0(groups[g], '_GB.fasta')))
}
first_ledger <- ledger(first, 'fixture-1')
align(first)
first_records <- coi_validate_retained_state(first, groups)$records
first_inputs <- lapply(groups, function(group) jsonlite::read_json(file.path(first, 'data/primerminer/full_aligned_reference', paste0(group, '_COX1_full_reference_aligned.fasta.inputs.json'))))
release1 <- file.path(fixture, 'immutable/fixture-1'); dir.create(release1, recursive=TRUE)
state1 <- coi_stage_retained_state(first, release1, 'fixture-1', 'https://data.example.invalid', groups)
manifest1 <- list(marker_id='COI', release_version='fixture-1', coi_retained_state=state1)
# A fresh worker has no raw data or ledger until the complete release is restored.
second <- worker('second')
stopifnot(!dir.exists(file.path(second, 'data/primerminer/raw')))
restored <- coi_restore_retained_state(second, list(manifest=manifest1, stale=FALSE, source='fixture immutable release'), groups=groups)
stopifnot(restored$retained_records==4L)
coi_assert_retention(first_records, coi_validate_retained_state(second, groups)$records)
stopifnot(identical(atlas_sha256(file.path(second,'data/provenance/coi_reference_accession_ledger.csv')),
                   atlas_sha256(file.path(first,'data/provenance/coi_reference_accession_ledger.csv'))))
# Simulate two newly acquired accession versions; previously frozen rows remain.
for(g in seq_along(groups)) {
  path <- file.path(second,'data/primerminer/raw',groups[g],paste0(groups[g],'_GB.fasta'))
  values <- coi_read_fasta(path)
  values <- c(values,setNames(paste0(substr(reference,1,900),'C',substr(reference,902,nchar(reference))),paste0('TEST',g,'C.1 newly acquired header')))
  write_fasta(values,path)
}
second_ledger <- ledger(second,'fixture-2',first_ledger)
align(second)
second_records <- coi_validate_retained_state(second,groups)$records
coi_assert_retention(first_records,second_records)
for(i in seq_along(groups)) {
  alignment <- file.path(second,'data/primerminer/full_aligned_reference',paste0(groups[i],'_COX1_full_reference_aligned.fasta'))
  inputs <- jsonlite::read_json(paste0(alignment,'.inputs.json'))
  stopifnot(inputs$raw_sequences==3L,inputs$aligned_templates==3L,
    !identical(inputs$input_sha256,first_inputs[[i]]$input_sha256),
    !identical(inputs$raw_sha256,first_inputs[[i]]$raw_sha256))
  coi_validate_full_alignment(alignment,file.path(second,'data/primerminer/raw',groups[i],paste0(groups[i],'_GB.fasta')),
    file.path(second,'data/reference/coi_reference_NC_001322.1.fasta'),inputs)
}
# The second successful state also restores in a third clean worker, with lineage.
release2 <- file.path(fixture,'immutable/fixture-2');dir.create(release2,recursive=TRUE)
state2 <- coi_stage_retained_state(second,release2,'fixture-2','https://data.example.invalid',groups)
third <- worker('third')
coi_restore_retained_state(third,list(manifest=list(marker_id='COI',release_version='fixture-2',coi_retained_state=state2),stale=FALSE,source='fixture immutable release'),groups=groups)
coi_assert_retention(second_records,coi_validate_retained_state(third,groups)$records)
proof <- jsonlite::read_json(file.path(third,'data/provenance/coi_reference_restore.json'))
stopifnot(proof$previous_release_id=='fixture-2',proof$retained_records==6L)
# Hash corruption, ledger-only restoration, stale pointers and rollback truncation
# all fail before replacing any existing retained inputs.
expect_error <- function(value) stopifnot(inherits(try(force(value),silent=TRUE),'try-error'))
expect_error(coi_restore_retained_state(third,list(manifest=manifest1,stale=FALSE,source='olderfixture'),groups=groups))
expect_error(coi_restore_retained_state(worker('stale'),list(manifest=manifest1,stale=TRUE),groups=groups))
expect_error(coi_restore_retained_state(worker('legacy'),list(manifest=list(marker_id='COI'),stale=FALSE),groups=groups))
ledger_only <- manifest1
ledger_only$coi_retained_state$artifacts <- ledger_only$coi_retained_state$artifacts[c('ledger','records')]
expect_error(coi_restore_retained_state(worker('ledger_only'),list(manifest=ledger_only,stale=FALSE),groups=groups))
ledger_path <- file.path(third,'data/provenance/coi_reference_accession_ledger.csv')
ledger_lines <- readLines(ledger_path)
unlink(ledger_path)
expect_error(coi_validate_retained_state(third,groups))
writeLines(ledger_lines,ledger_path)
ledger_without_header <- read.csv(ledger_path,stringsAsFactors=FALSE)
ledger_without_header$header_sha256 <- NULL
write.csv(ledger_without_header,ledger_path,row.names=FALSE)
expect_error(coi_validate_retained_state(third,groups))
writeLines(ledger_lines,ledger_path)
no_header <- file.path(fixture,'missing_header.fasta')
writeLines(c('>',reference),no_header)
expect_error(coi_read_fasta(no_header))
# Header-only changes preserve accession and bases but still violate immutable
# record provenance; missing raw inputs cannot be substituted by ledger counts.
raw_path <- file.path(third,'data/primerminer/raw',groups[1],paste0(groups[1],'_GB.fasta'))
raw_lines <- readLines(raw_path)
changed_headers <- raw_lines;changed_headers[1] <- paste0(changed_headers[1],' changed provenance')
writeLines(changed_headers,raw_path)
expect_error(coi_validate_retained_state(third,groups))
writeLines(raw_lines,raw_path)
unlink(raw_path)
expect_error(coi_validate_retained_state(third,groups))
writeLines(raw_lines,raw_path)
corrupt <- state1$artifacts$ledger$local_path
write('corrupt',corrupt,append=TRUE)
expect_error(coi_restore_retained_state(worker('corrupt'),list(manifest=manifest1,stale=FALSE),groups=groups))
# A raw-byte change cannot hide behind the same accession count.
raw_path <- file.path(third,'data/primerminer/raw',groups[1],paste0(groups[1],'_GB.fasta'))
values <- coi_read_fasta(raw_path);changed_base <- if(substr(values[1],1,1)=='A') 'C' else 'A'
values[1] <- paste0(changed_base,substr(values[1],2,nchar(values[1])));write_fasta(values,raw_path)
expect_error(coi_validate_retained_state(third,groups))
message('Two successive COI clean-worker builds retain exact prior bytes/metadata, append records, invalidate alignments and restore complete immutable state.')
