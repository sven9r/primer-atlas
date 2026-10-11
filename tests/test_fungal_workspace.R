#!/usr/bin/env Rscript
.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
suppressPackageStartupMessages({library(shiny); library(DT); library(bslib)})
source("R/functions.R")
source("R/data_layer.R")
source("R/fungal_explorer.R")
source("R/fungal_versions.R")
source("R/fungal_alignment.R")
source("R/fungal_alignment_ui.R")
source("R/fungal_ui.R")
versions <- fungal_read_index()
stopifnot(nrow(versions) >= 2L)
shiny::testServer(fungalExplorerServer, args=list(id="test"), {
  session$setInputs(version=tail(versions$id,1),platform="ALL",illumina_max=500,pacbio_max=2000,
                    nanopore_max=3500,phylum="ALL",class="ALL",locus="ALL",search="",
                    examples_only=TRUE,show_unresolved=TRUE)
  session$flushReact()
  stopifnot(length(snapshot()$sequences)>=21L,nrow(candidates())==20L,
            all(c("Illumina","PacBio","Nanopore") %in% candidates()$length_tier))
  stopifnot(alignment_index()$n_columns==nchar(snapshot()$alignment[1]),
            nrow(mapped_alignment())==20L*21L,nrow(alignment_pairs())==20L,
            alignment_reference() %in% selected_accessions())
  overview <- output$alignment_overview
  bases <- output$alignment_bases
  stopifnot(grepl("frozen alignment",overview$html,fixed=TRUE),
            grepl("fungal-alignment-base",bases$html,fixed=TRUE))
  # Platform controls restrict the actual pair-template length evidence.
  session$setInputs(platform="ILLUMINA",examples_only=FALSE,show_unresolved=FALSE)
  session$flushReact()
  illumina_ids <- candidates()$pair_id
  stopifnot(length(illumina_ids)>0,all(candidates()$n_illumina_eligible>0))
  session$setInputs(platform="NANOPORE")
  session$flushReact()
  stopifnot(all(illumina_ids %in% candidates()$pair_id),nrow(candidates())>length(illumina_ids))
  long_id <- candidates()$pair_id[candidates()$pair_label=="FF700 + LR7"]
  session$setInputs(pair_click=list(pair=long_id,nonce=1))
  session$flushReact()
  stopifnot(active_pair()$pair_id==long_id,active_pair()$amplicon_median_bp>2000,
            active_pair()$amplicon_max_bp<=3500)
  session$setInputs(alignment_scope="active")
  session$flushReact()
  stopifnot(nrow(alignment_pairs())==1L,alignment_pairs()$pair_id==long_id,
            all(alignment_focus_row()$pair_id==long_id))
  session$setInputs(alignment_scope="compared",compare_pairs=long_id)
  session$flushReact()
  stopifnot(identical(alignment_pairs()$pair_id,long_id))
  session$setInputs(alignment_scope="shown")
  # Five explicit comparisons remain visible even when a different pair is active.
  comparisons <- head(candidates()[candidates()$pair_id!=long_id,,drop=FALSE],5L)
  session$setInputs(compare_pairs=comparisons$pair_id,alignment_reference="FN812768.2")
  session$flushReact()
  overview <- output$alignment_overview
  stopifnot(grepl(paste0("P5 = ",comparisons$pair_label[5]),overview$html,fixed=TRUE),
            alignment_focus_row()$accession=="FN812768.2",
            alignment_focus_row()$alignment_status!="mapped",
            grepl("no complete exact-site product on",output$alignment_status$html,fixed=TRUE))
  session$setInputs(compare_pairs=long_id)
  # Class filters change the true denominator and data exported by the module.
  session$setInputs(class="Agaricomycetes")
  session$flushReact()
  stopifnot(length(selected_accessions())==4L,all(all_pairs()$n_panel==4L),
            all(sequence_rows()$class=="Agaricomycetes"))
  stopifnot(all(mapped_alignment()$accession %in% selected_accessions()),
            alignment_reference() %in% selected_accessions())
  session$setInputs(platform="ALL",show_unresolved=TRUE,examples_only=TRUE)
  session$flushReact()
  stopifnot(all(candidates()$n_scorable>0 | candidates()$evidence_status=="documented_pair"))
  # Incompatible lineage selections cannot leave unsupported global candidates.
  session$setInputs(phylum="Ascomycota",class="Agaricomycetes")
  session$flushReact()
  stopifnot(length(selected_accessions())==0L,nrow(candidates())==0L,nrow(active_pair())==0L)
  stopifnot(nrow(mapped_alignment())==0L,is.null(alignment_reference()))
  # Historical selection loads the original sequence and its unresolved taxonomy.
  session$setInputs(version=versions$id[1],phylum="ALL",class="ALL",platform="ALL",examples_only=TRUE)
  session$flushReact()
  stopifnot(length(snapshot()$sequences)==1L,identical(names(snapshot()$sequences),"FN812768.2"),
            all(all_pairs()$n_panel==1L),!any(all_pairs()$length_tier=="Nanopore"))
  stopifnot(length(alignment_index()$alignment)==1L,alignment_reference()=="FN812768.2",
            all(mapped_alignment()$accession=="FN812768.2"))
  # Empty searches remain a useful empty state with no stale selected evidence.
  session$setInputs(search="no_such_primer_or_clade_98765")
  session$flushReact()
  stopifnot(nrow(candidates())==0L,nrow(active_pair())==0L)
})
message("Fungal workspace platform, class, selection, and historical switching checks passed.")
