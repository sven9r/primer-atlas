#!/usr/bin/env Rscript
.libPaths(c(file.path(getwd(),".Rlib"),.libPaths()))
source("R/functions.R")
source("R/data_layer.R")
source("R/fungal_explorer.R")
source("R/fungal_versions.R")
source("R/fungal_ui.R")
source("R/fungal_primer_overview.R")

snapshot <- fungal_read_snapshot("v3-c6fa87f256e1",remote=FALSE)
starters <- fungal_starter_examples(fungal_platform_summary(snapshot$pair_sequences,snapshot$pairs))
stopifnot(nrow(starters)==20L,nrow(snapshot$pairs)==884L)
original <- serialize(snapshot,NULL)
all <- fungal_primer_overview_data(snapshot,starters,"PZ962702.1","ALL")
stopifnot(all$summary$n_candidates==20L,all$summary$n_mapped==20L,all$summary$n_missing==0L,
          all$summary$n_snapshot_pairs==884L,all$reference_length==nchar(snapshot$sequences[["PZ962702.1"]]),
          all(all$pair_meta$accession=="PZ962702.1"),all(all$pair_meta$version_id==snapshot$manifest$id),
          all(all$pair_meta$folmer_overlap_bp %in% NA_real_),all$marker$marker_id=="ITS_FUNGAL")

# This physical product spans 1,573 frozen alignment columns but just 335 raw
# bases. The shared overview uses raw bounds and subtracts both primer widths.
short <- all$pair_meta[all$pair_meta$pair_label=="ITS1 + ITS2",]
stopifnot(nrow(short)==1L,short$start==1775L,short$end==2109L,
          short$aligned_amplicon_bp==335L,short$targeted_region_bp==296L,
          short$pair_start==short$start,short$pair_end==short$end,
          short$aligned_amplicon_bp==short$end-short$start+1L,
          short$targeted_region_bp==short$reverse_start-short$forward_end-1L,
          all(unname(all$domain)==c(775,4109)))
for (i in seq_len(nrow(all$pair_meta))) {
  row <- all$pair_meta[i,]
  widths <- row$forward_end-row$forward_start+1L+row$reverse_end-row$reverse_start+1L
  stopifnot(row$aligned_amplicon_bp-row$targeted_region_bp==widths,
            row$forward_end<row$reverse_start,row$end<=all$reference_length)
}

# The adapter ignores gap placements, reversed MAFFT rows, panel medians and
# pooled geometry. A corrupt frozen encoded length is withheld rather than used.
pair <- starters[starters$pair_label=="ITS1 + ITS2",]
reverse_alignment <- snapshot
reverse_alignment$alignment <- c("_R_PZ962702.1"="----unrelated_alignment_columns----")
same <- fungal_primer_overview_data(reverse_alignment,pair,"PZ962702.1")
stopifnot(same$pair_meta$aligned_amplicon_bp==335L)
bad <- snapshot
owned <- bad$pair_sequences$pair_id==pair$pair_id & bad$pair_sequences$accession=="PZ962702.1"
bad$pair_sequences$amplicon_bp[owned] <- 1573L
withheld <- fungal_primer_overview_data(bad,pair,"PZ962702.1")
stopifnot(nrow(withheld$pair_meta)==0L,nrow(withheld$missing_pairs)==1L,
          withheld$missing_pairs$placement_status=="inconsistent_frozen_product_evidence")

# A real pair supported elsewhere cannot be placed on FN812768 by borrowing
# those coordinates. Exact individual sites remain distinct from a product.
pilot <- fungal_primer_overview_data(snapshot,starters,"FN812768.2")
stopifnot(pilot$summary$n_mapped==6L,pilot$summary$n_missing==14L,
          pilot$summary$n_mapped+pilot$summary$n_missing==20L,
          all(!pilot$missing_pairs$pair_id %in% pilot$pair_meta$pair_id),
          !any(c("start","end","pair_start","pair_end") %in% names(pilot$missing_pairs)))
ambiguous <- snapshot
ambiguous$pair_sequences$status[owned] <- "ambiguous_binding_site"
ambiguous$pair_sequences$reverse_status[owned] <- "ambiguous_binding_site"
unresolved <- fungal_primer_overview_data(ambiguous,pair,"PZ962702.1")
stopifnot(nrow(unresolved$pair_meta)==0L,
          is.na(unresolved$primers$start[unresolved$primers$direction=="reverse"]),
          unresolved$primers$start[unresolved$primers$direction=="forward"]==1775L)

# Coarse filtering uses actual raw landmark overlap. A generic ITS catalog
# label cannot turn an SSU-contained product into an ITS product, or vice versa.
ssu <- fungal_primer_overview_data(snapshot,starters,"PZ962702.1","SSU")
its <- fungal_primer_overview_data(snapshot,starters,"PZ962702.1","ITS")
lsu <- fungal_primer_overview_data(snapshot,starters,"PZ962702.1","LSU")
ssu_only <- "FF700 + SSU1196Rngs"
stopifnot(ssu_only %in% ssu$pair_meta$pair_label,!ssu_only %in% its$pair_meta$pair_label,
          !ssu_only %in% lsu$pair_meta$pair_label,
          "ITS3 + ITS4" %in% its$pair_meta$pair_label,!"ITS3 + ITS4" %in% ssu$pair_meta$pair_label,
          all(its$pair_meta$region_match_basis=="own_template_landmark_overlap"),
          identical(unname(its$domain),c(805,4422)),
          lsu$domain["min"]==1423L,lsu$domain["max"]==min(lsu$reference_length,7821L))
catalog_only <- fungal_primer_overview_data(snapshot,starters,"FN812768.2","ITS")
stopifnot(nrow(catalog_only$landmarks)==0L,
          all(catalog_only$pair_meta$region_match_basis=="reported_component_locus_without_local_region_boundaries"))

# Domain coordinates are own-accession projected structural evidence. They
# remain separate from gene/spacer bars and focus without fabricated D labels.
d1 <- fungal_primer_overview_data(snapshot,starters,"PZ962702.1","landmark:28S_D1")
stopifnot(nrow(d1$domains)==9L,nrow(d1$landmarks)==5L,
          all(d1$domains$evidence_status=="mapped_reference_domain"),
          all(d1$domains$boundary_status=="projected"),
          !any(d1$landmarks$landmark_type=="divergent_domain"),
          d1$focus$start==2543L,d1$focus$end==2690L,identical(unname(d1$domain),c(1543,4690)),
          all(d1$pair_meta$end>=2543L & d1$pair_meta$start<=2690L),
          "landmark:28S_D1" %in% d1$available_regions$region_id,
          !any(grepl("D6|D7b|D12",d1$domains$label)))
absent_domain <- fungal_primer_overview_data(snapshot,starters,"FN812768.2","landmark:28S_D1")
stopifnot(nrow(absent_domain$pair_meta)==0L,length(absent_domain$domain)==0L,
          absent_domain$focus$basis=="region_coordinates_unavailable")

# Preserve the caller's selected row order and actual limit. Missing products
# are counted within that same limited selection, not replaced with other pairs.
limited <- fungal_primer_overview_data(snapshot,starters,"PZ962702.1",max_rows=3L)
stopifnot(identical(limited$pair_meta$pair_id,starters$pair_id[1:3]),limited$summary$n_shown_candidates==3L,
          limited$summary$n_truncated==17L,all(limited$pair_meta$n_panel==starters$n_panel[1:3]),
          all(limited$pair_meta$n_scorable==starters$n_scorable[1:3]),
          all(limited$pair_meta$evidence_status==starters$evidence_status[1:3]))
empty <- fungal_primer_overview_data(snapshot,starters[FALSE,],"PZ962702.1","ITS")
zero <- fungal_primer_overview_data(snapshot,starters,"PZ962702.1",max_rows=0L)
stopifnot(nrow(empty$pair_meta)==0L,nrow(empty$primers)==0L,nrow(empty$missing_pairs)==0L,
          empty$summary$n_candidates==0L,nrow(zero$pair_meta)==0L,zero$summary$n_truncated==20L)

# Primary citation keys can collide across components. Each primer keeps its
# own full catalog identity, citation and original source URL in this namespace.
collision <- snapshot
collision$catalog$primary_reference_key <- "same_source_key"
citations <- fungal_primer_overview_data(collision,starters,"PZ962702.1")
stopifnot(!anyDuplicated(citations$citations$key))
for (i in seq_len(nrow(citations$primers))) {
  primer <- citations$primers[i,];source <- citations$citations[citations$citations$key==primer$source_key,]
  original_component <- snapshot$catalog[snapshot$catalog$catalog_id==primer$catalog_id,]
  stopifnot(nrow(source)==1L,source$catalog_id==primer$catalog_id,
            source$short_citation==original_component$reference,
            source$primary_reference_url==original_component$primary_reference_url,
            source$source_url==original_component$source_url)
}
stopifnot(nrow(citations$primers)==40L,nrow(citations$citations)>1L)

# Historical pilot approximations apply only to their captured reference, and
# its 18S/28S conserved anchors retain partial status rather than becoming genes.
v1 <- fungal_read_snapshot("v1-9e958f105c46",remote=FALSE)
historical <- fungal_primer_overview_data(v1,v1$pairs,"FN812768.2","ITS")
stopifnot(all(historical$landmarks$accession=="FN812768.2"),nrow(historical$domains)==0L,
          all(historical$landmarks$boundary_status[historical$landmarks$landmark_type=="conserved_anchor"]=="partial"),
          historical$version_id==v1$manifest$id,
          !any(historical$primers$source_key %in% all$primers$source_key))

# Default reference selection is deterministic and based on real supported
# products, with complete annotation breaking equal mapped-product counts.
chosen <- fungal_primer_overview_default_accession(snapshot,starters,"ALL")
stopifnot(chosen %in% c("PZ962702.1","PZ276773.1"),
          fungal_primer_overview_data(snapshot,starters,chosen,"ALL")$summary$n_mapped==20L,
          identical(fungal_primer_overview_default_accession(snapshot,starters,"ALL"),chosen),
          fungal_primer_overview_data(snapshot,starters,region="ALL")$accession==chosen,
          identical(serialize(snapshot,NULL),original))
foreign <- starters[1L,];foreign$pair_id <- "NOT_IN_THIS_SNAPSHOT"
stopifnot(inherits(try(fungal_primer_overview_data(snapshot,foreign,"PZ962702.1"),silent=TRUE),"try-error"),
          inherits(try(fungal_primer_overview_data(snapshot,starters,"FOREIGN"),silent=TRUE),"try-error"))
message("Fungal shared-overview raw lengths, ownership, regions, provenance and historical fixtures passed.")
