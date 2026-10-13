#!/usr/bin/env Rscript
.libPaths(c(file.path(getwd(),'.Rlib'),.libPaths()))
source('R/functions.R');source('R/data_layer.R');source('R/fungal_explorer.R');source('R/fungal_versions.R')
source('R/fungal_alignment.R');source('R/fungal_ui.R');source('R/fungal_primer_overview.R')
source('R/primer_overview.R');source('R/unite_reference.R');source('R/unite_reference_ui.R');source('R/fungal_combination_builder.R');source('R/primer_publication_export.R')

snapshot<-fungal_read_snapshot('v3-c6fa87f256e1',remote=FALSE)
original<-serialize(snapshot,NULL)
labels<-c('ITS1 + ITS2','ITS1F + ITS2','LR0R + LR3')
pairs<-snapshot$pairs[snapshot$pairs$pair_label %in% labels,,drop=FALSE]
pairs<-pairs[match(labels,pairs$pair_label),,drop=FALSE]
stopifnot(nrow(pairs)==3L)
ids<-pairs$pair_id
unite<-unite_reference_read()
stopifnot(!is.null(unite),nrow(unite$metadata)==15402L)
broad_rows<-lapply(ids,function(id) {
  fit<-unite_reference_pair_fit(unite,id,'class')
  columns<-c('n_selected','n_exact_products','n_single_sites','n_no_exact_sites','n_unresolved','n_unavailable','n_withheld')
  counts<-as.data.frame(lapply(fit$summary[,columns,drop=FALSE],sum))
  data.frame(pair_id=id,group_id='ALL',taxon='All retained references',parent_path='',rank='class',counts,stringsAsFactors=FALSE)
})
broad<-list(summary=do.call(rbind,broad_rows),metadata=list(version_id=unite$manifest$id,rank='class',group_id='ALL',source_doi=unite$manifest$source_doi,n_selected=15402L))
lookup<-read.csv('data/catalog/primer_publication_citations.csv',stringsAsFactors=FALSE)
context<-fungal_publication_context(snapshot,ids,'PZ962702.1',broad_evidence=broad,
  filters=list(min_bp=300,max_bp=340,platform='ALL',figure_pair_ids=ids[1L],region='ITS'),citation_lookup=lookup,
  title='Primers & "controls" <RNA>')
export<-primer_publication_export(context)
stopifnot(nrow(export$geometry_table)==3L,nrow(export$table)==3L,
  identical(export$table$pair_id,ids),sum(export$table$included_in_figure)==1L,
  export$table$product_including_primers_bp[1L]==335L,export$table$product_excluding_primers_bp[1L]==296L,
  export$table$forward_length_bp[1L]+export$table$reverse_length_bp[1L]==39L,
  export$table$forward_start[1L]==1775L,export$table$reverse_end[1L]==2109L,
  all(export$table$exclusion_reason[-1L]=='length'),all(export$table$ncbi_snapshot_references==21L),
  all(export$table$unite_selected_references==15402L),all(export$table$unite_snapshot_id==unite$manifest$id),
  identical(export$table$combination_status,pairs$evidence_status),
  all(export$table$unite_exact_products==broad$summary$n_exact_products),
  grepl('335',export$svg,fixed=TRUE),grepl('1-based',export$caption,fixed=TRUE),
  grepl('amplification performance has not been validated',export$caption,fixed=TRUE),
  grepl('candidate_published_components',export$html,fixed=TRUE),
  grepl('Candidate from published components',export$html,fixed=TRUE),
  grepl(paste0(.pub_evidence_label(pairs$evidence_status[1L]),': figure row 1.'),export$caption,fixed=TRUE),
  grepl('data-product-bp="335"',export$svg,fixed=TRUE),
  !grepl('data-product-bp="1573"',export$svg,fixed=TRUE),
  grepl('xmlns="http://www.w3.org/2000/svg"',export$svg,fixed=TRUE),
  !grepl('<script',export$svg,fixed=TRUE),grepl('&lt;RNA&gt;',export$html,fixed=TRUE))
stopifnot(export$table$forward_bibliography_status[1L]=='verified_bibliography_supplement',
  export$table$forward_doi[1L]=='10.1016/b978-0-12-372180-8.50042-1',
  grepl('PCR Protocols',export$table$forward_full_reference[1L],fixed=TRUE),
  identical(context$citations$catalog_citation,export$citations$catalog_citation),
  identical(context$geometry_table$forward_full_reference,export$geometry_table$forward_full_reference))
stopifnot(!grepl('\u2010',export$svg,fixed=TRUE),grepl('\u2010',export$html,fixed=TRUE),
  any(grepl('\u2010',export$table$forward_full_reference,fixed=TRUE)))

# An explicit binding-site order controls figure numbering independently of
# requested/table order. Reversed selections and ordered subsets both retain
# that order in the SVG, table links, structured context, and caption lists.
for(figure_ids in list(rev(ids),ids[c(3L,1L)])) {
  ordered<-primer_publication_export(fungal_publication_context(snapshot,ids,'PZ962702.1',broad_evidence=broad,
    filters=list(figure_pair_ids=figure_ids,region='ALL'),citation_lookup=lookup))
  svg_rows<-regmatches(ordered$svg,gregexpr('class="primer-overview-row" data-pair-id="[^"]+"',ordered$svg))[[1L]]
  svg_ids<-sub('.*data-pair-id="([^"]+)".*','\\1',svg_rows)
  expected_rows<-match(ids,figure_ids)
  stopifnot(identical(ordered$context$figure_pair_ids,figure_ids),identical(ordered$table$pair_id,ids),
    identical(ordered$geometry_table$figure_row,expected_rows),identical(ordered$table$figure_row,expected_rows),
    identical(svg_ids,figure_ids))
  for(status in unique(ordered$geometry_table$combination_status[ordered$geometry_table$included_in_figure])) {
    rows<-sort(ordered$geometry_table$figure_row[ordered$geometry_table$included_in_figure & ordered$geometry_table$combination_status==status])
    label<-paste0(.pub_evidence_label(status),': figure ',if(length(rows)==1L) 'row ' else 'rows ',paste(rows,collapse=', '),'.')
    stopifnot(grepl(label,ordered$caption,fixed=TRUE))
  }
}

# A very different frozen gap layout cannot inflate the measured335bp span.
index<-fungal_alignment_index(snapshot$alignment,snapshot$sequences)
projected<-fungal_alignment_pair_rows(pair_sequences=snapshot$pair_sequences,index=index)
row<-projected[projected$pair_id==ids[1L] & projected$accession=='PZ962702.1',,drop=FALSE]
stopifnot(row$amplicon_bp==335L,row$reverse_end_col-row$forward_start_col+1L==1573L)
changed<-snapshot;changed$alignment<-c('_R_PZ962702.1'='----unrelated gapped coordinate scale----')
changed_export<-primer_publication_export(fungal_publication_context(changed,ids,'PZ962702.1',broad_evidence=broad,
  filters=list(min_bp=300,max_bp=340,platform='ALL',figure_pair_ids=ids[1L],region='ITS'),citation_lookup=lookup,title='Primers & "controls" <RNA>'))
stopifnot(identical(export$svg,changed_export$svg),identical(export$table,changed_export$table),identical(serialize(snapshot,NULL),original))

# The table retains requests outside the focus and those unavailable on the
# selected template. They never borrow measured coordinates from another row.
unmapped<-primer_publication_export(fungal_publication_context(snapshot,ids,'FN812768.2',broad_evidence=broad,
  filters=list(figure_pair_ids=character(),region='LSU'),citation_lookup=lookup))
stopifnot(nrow(unmapped$table)==3L,!any(unmapped$table$included_in_figure),
  any(is.na(unmapped$table$product_including_primers_bp)),
  all(unmapped$table$exclusion_reason[is.na(unmapped$table$product_including_primers_bp)]=='no_map'),
  grepl('No exact-supported product',unmapped$svg,fixed=TRUE))
no_broad<-primer_publication_export(fungal_publication_context(snapshot,ids[1L],'PZ962702.1'))
stopifnot(is.na(no_broad$table$unite_selected_references),is.na(no_broad$table$unite_exact_products),
  no_broad$table$taxonomic_evidence_status=='unavailable',
  no_broad$table$forward_bibliography_status=='short_citation_only_full_reference_unavailable',
  !nzchar(no_broad$table$forward_doi),grepl('full-reference status',no_broad$methods,fixed=TRUE))
invalid<-broad;invalid$summary$n_exact_products[1L]<-15403L
stopifnot(inherits(try(fungal_publication_context(snapshot,ids,'PZ962702.1',broad_evidence=invalid),silent=TRUE),'try-error'))
invalid<-broad;invalid$summary$n_no_exact_sites[1L]<-invalid$summary$n_no_exact_sites[1L]+1L
stopifnot(inherits(try(fungal_publication_context(snapshot,ids,'PZ962702.1',broad_evidence=invalid),silent=TRUE),'try-error'))
invalid<-broad;invalid$summary$version_id<-'another-version'
stopifnot(inherits(try(fungal_publication_context(snapshot,ids,'PZ962702.1',broad_evidence=invalid),silent=TRUE),'try-error'))
wrong<-context$overview;wrong$pair_meta$aligned_amplicon_bp[wrong$pair_meta$pair_id==ids[1L]]<-1573L
stopifnot(inherits(try(primer_publication_context(wrong,ids,snapshot$sequences[['PZ962702.1']]),silent=TRUE),'try-error'))

# All25 selected crosses survive the paper table; the figure includes only the
# actual filtered products and new candidates retain their transient identities.
catalog<-fungal_combination_catalog(snapshot)
forward<-catalog$catalog_id[match(c('ITS1','ITS1F','ITS5','ITS3','fITS7'),catalog$primer_name)]
reverse<-catalog$catalog_id[match(c('ITS2','ITS4','ITS4B','LR21','LR3'),catalog$primer_name)]
crosses<-fungal_combination_builder(snapshot,forward,reverse,'PZ962702.1',platform='Illumina',length_range=c(300,400),region='ITS')
all25<-primer_publication_export(fungal_publication_context(crosses$working_snapshot,crosses$selected_pairs$pair_id,
  'PZ962702.1',filters=list(min_bp=300,max_bp=400,platform='Illumina',region='ITS',figure_pair_ids=crosses$map_pair_rows$pair_id),citation_lookup=lookup))
stopifnot(nrow(all25$table)==25L,identical(all25$table$pair_id,crosses$selected_pairs$pair_id),
  identical(all25$table$included_in_figure,crosses$selected_pairs$included_in_map),
  identical(all25$table$exclusion_reason,crosses$selected_pairs$filter_reason),
  identical(all25$table$combination_status,crosses$selected_pairs$evidence_status),
  all(all25$table$product_including_primers_bp[all25$table$included_in_figure]>=300L),
  all(all25$table$product_including_primers_bp[all25$table$included_in_figure]<=400L),
  !any(all25$table$combination_status[startsWith(all25$table$pair_id,'FUNGAL_COMB_')]=='documented_pair'),
  identical(serialize(snapshot,NULL),original))

scratch<-tempfile('primer-publication-');dir.create(scratch)
for(format in c('svg','csv','html','caption','methods','citations','landmarks')) {
  path<-file.path(scratch,paste0('export.',format));primer_publication_write(export,path,format)
  stopifnot(file.exists(path),file.info(path)$size>0L)
}
roundtrip<-read.csv(file.path(scratch,'export.csv'),stringsAsFactors=FALSE)
stopifnot(nrow(roundtrip)==nrow(export$table),identical(roundtrip$pair_id,ids),
  identical(roundtrip$product_including_primers_bp,as.integer(export$table$product_including_primers_bp)),
  identical(roundtrip$forward_full_reference,export$table$forward_full_reference))
if(requireNamespace('xml2',quietly=TRUE)) {
  document<-xml2::read_xml(export$svg)
  metadata<-jsonlite::fromJSON(xml2::xml_text(xml2::xml_find_first(document,'.//*[local-name()="metadata"]')))
  stopifnot(xml2::xml_name(document)=='svg',length(xml2::xml_find_all(document,'.//*[local-name()="g" and @class="primer-overview-row"]'))==1L,
    metadata$reference_accession=='PZ962702.1',metadata$ncbi_version==snapshot$manifest$id,
    metadata$unite_version==unite$manifest$id)
  xml2::read_html(export$html)
}
if(primer_publication_png_available()) {
  image<-file.path(scratch,'export.png');primer_publication_write(export,image,'png',width=1680L)
  bytes<-readBin(image,'raw',8L);stopifnot(identical(as.integer(bytes),c(137L,80L,78L,71L,13L,10L,26L,10L)))
}

# Export the genuine selected regional window with its gapped axis and exact
# displayed base strings. Its caption cannot relabel columns as product bp.
owner<-unite$alignment_index[unite$alignment_index$rank=='order' & unite$alignment_index$region=='ITS2' & unite$alignment_index$n_input>1L,,drop=FALSE][1L,]
alignment<-unite_reference_alignment(unite,owner$group_id,'ITS2',limit=4L)
window<-.unite_reference_ui_alignment_window(alignment,start=2L,width=60L)
regional<-primer_publication_alignment_export(window,'ITS2',owner$group_id,unite$manifest$id,unite$manifest$source_doi,pair=labels[1L])
stopifnot(identical(regional$table$aligned_window_sequence,unname(window$alignment)),
  all(regional$table$start_column==window$start),all(regional$table$end_column==window$end),
  grepl('alignment gaps',regional$caption,fixed=TRUE),grepl('data-state="gap"',regional$svg,fixed=TRUE),
  grepl(unite$manifest$id,regional$caption,fixed=TRUE))
if(requireNamespace('xml2',quietly=TRUE)) xml2::read_xml(regional$svg)
badwindow<-window;badwindow$end<-badwindow$end+1L
stopifnot(inherits(try(primer_publication_alignment_export(badwindow,'ITS2',owner$group_id,unite$manifest$id,unite$manifest$source_doi),silent=TRUE),'try-error'))
unlink(scratch,recursive=TRUE)
message('Publication SVG/HTML/CSV lengths, selection, source versions, bibliography, true regional-window sequences and optionalPNG passed.')
