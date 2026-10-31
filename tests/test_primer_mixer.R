source("R/primer_mixer.R")

expect_error <- function(code, pattern) {
  error <- tryCatch({ force(code); NULL }, error=function(e)e)
  stopifnot(inherits(error,"error"),grepl(pattern,conditionMessage(error)))
}
set.seed(7319)
dna <- function(n) paste(sample(c("A","C","G","T"),n,replace=TRUE),collapse="")
forward <- dna(20L); reverse <- dna(21L)
spacer <- dna(83L); left <- dna(17L); right <- dna(19L)
raw <- paste0(left,forward,spacer,.primer_mixer_complement(reverse),right)
f <- primer_mixer_custom("Forward",forward,"forward","12S_MT")
r <- primer_mixer_custom("Reverse",reverse,"reverse","12S_MT")
catalog <- rbind(f,r)
templates <- data.frame(template_id="fish_1",sequence=raw,marker_id="12S_MT",organism="A test template",stringsAsFactors=FALSE)
result <- primer_mixer_build(catalog,f$primer_id,r$primer_id,templates,"12S_MT")
p <- result$products
stopifnot(nrow(p)==1L,p$status=="exact_product",p$forward_start==18L,p$forward_end==37L,
  p$reverse_start==121L,p$reverse_end==141L,p$amplicon_bp==124L,p$internal_bp==83L,
  p$forward_mismatches==0L,p$reverse_mismatches==0L,result$summary$n_templates==1L,
  result$summary$n_exact_products==1L,p$passes_length_filter)

# Every requested component contributes to the Cartesian set. The length
# filter changes eligible rows without erasing measured sites or requests.
f2 <- primer_mixer_custom("Forward alternate",paste0(ifelse(substr(forward,1L,1L)=="A","C","A"),substring(forward,2L)),"forward","12S_MT")
r2 <- primer_mixer_custom("Reverse alternate",paste0(ifelse(substr(reverse,1L,1L)=="G","T","G"),substring(reverse,2L)),"reverse","12S_MT")
cart <- primer_mixer_build(rbind(catalog,f2,r2),c(f$primer_id,f2$primer_id),c(r$primer_id,r2$primer_id),templates,"12S_MT",max_mismatches=1L,min_bp=125L)
stopifnot(nrow(cart$pairs)==4L,nrow(cart$products)==4L,sum(cart$products$status=="exact_product")==1L,
  sum(cart$products$status=="mismatch_product")==3L,!any(cart$products$passes_length_filter),all(cart$products$amplicon_bp==124L))

# Terminal mismatch counts follow the supplied oligo's 5'->3' orientation.
# Reverse 3' corresponds to the leftmost site base on the displayed template.
mutate <- function(sequence,position) {
  base <- substr(sequence,position,position)
  substr(sequence,position,position) <- ifelse(base=="A","C","A")
  sequence
}
mutated <- templates
mutated$sequence <- mutate(mutate(raw,p$forward_end),p$reverse_start)
strict <- primer_mixer_build(catalog,f$primer_id,r$primer_id,mutated,"12S_MT")
loose <- primer_mixer_build(catalog,f$primer_id,r$primer_id,mutated,"12S_MT",max_mismatches=1L)
stopifnot(strict$products$status=="no_eligible_sites",loose$products$status=="mismatch_product",
  loose$products$forward_mismatches==1L,loose$products$reverse_mismatches==1L,
  loose$products$forward_terminal_5_mismatches==1L,loose$products$reverse_terminal_5_mismatches==1L,
  all(loose$sites$terminal_base_mismatch),identical(loose$sites$mismatch_positions_5_to_3,c("20","21")))

# Equal sites and unknown reference bases are withheld, never called exact.
repeated <- templates; repeated$sequence <- paste0(raw,dna(13L),raw)
ties <- primer_mixer_build(catalog,f$primer_id,r$primer_id,repeated,"12S_MT")
stopifnot(ties$products$status=="ambiguous_sites",all(is.na(ties$sites$start)),all(ties$sites$n_best_sites==2L))
unknown <- templates; substr(unknown$sequence,p$forward_start,p$forward_start) <- "N"
ambiguous <- primer_mixer_build(catalog,f$primer_id,r$primer_id,unknown,"12S_MT")
stopifnot(ambiguous$products$status=="ambiguous_sites",ambiguous$sites$status[ambiguous$sites$direction=="forward"]=="ambiguous_template_candidate")
degenerate <- f
degenerate$sequence <- paste0("N",substring(forward,2L))
degenerate$primer_id <- "degenerate_forward"
compatible <- primer_mixer_build(rbind(degenerate,r),degenerate$primer_id,r$primer_id,templates,"12S_MT")
stopifnot(compatible$products$status=="exact_product")

# Missing sequence, too-short fragments, partial sites and wrong markers
# retain distinct states. A fish marker cannot borrow COI coordinates.
missing <- templates; missing$sequence <- NA_character_
short <- templates; short$sequence <- "ACGTACGT"
partial <- templates; partial$sequence <- substr(raw,1L,p$reverse_start-1L)
wrong <- templates; wrong$marker_id <- "COI"
stopifnot(primer_mixer_build(catalog,f$primer_id,r$primer_id,missing,"12S_MT")$products$status=="template_unavailable",
  primer_mixer_build(catalog,f$primer_id,r$primer_id,short,"12S_MT")$products$status=="template_too_short",
  primer_mixer_build(catalog,f$primer_id,r$primer_id,partial,"12S_MT")$products$status=="one_eligible_site",
  primer_mixer_build(catalog,f$primer_id,r$primer_id,wrong,"12S_MT")$products$status=="marker_mismatch")
none <- primer_mixer_build(catalog,f$primer_id,r$primer_id,NULL,"12S_MT")
stopifnot(nrow(none$pairs)==1L,nrow(none$products)==0L,none$summary$n_templates==0L,
  "sequence" %in% names(none$templates))

# Source catalog and input chemicals survive normalization. Invalid input
# fails explicitly instead of silently stripping punctuation or alignment gaps.
converted <- primer_mixer_custom("Inosine", "acgt acgtIUU", "forward", "COI")
stopifnot(converted$source_sequence=="ACGTACGTIUU",converted$sequence=="ACGTACGTNTT",
  grepl("not biochemical",converted$normalization_note),converted$source_key=="custom")
reversed <- primer_mixer_reverse(converted)
stopifnot(reversed$direction=="reverse",reversed$parent_primer_id==converted$primer_id,
  reversed$parent_sequence=="ACGTACGTIUU",reversed$sequence==.primer_mixer_complement(converted$sequence),
  reversed$catalog_origin=="derived_reverse_complement",grepl("no empirical validation",reversed$normalization_note),
  reversed$citation==converted$citation,reversed$primer_id!=converted$primer_id,
  converted$direction=="forward",identical(primer_mixer_reverse(converted)$primer_id,reversed$primer_id))
expect_error(primer_mixer_reverse(reversed),"original component")
expect_error(primer_mixer_custom("bad","ACGT----ACGT","forward","COI"),"ungapped")
expect_error(primer_mixer_custom("bad","ACGTZACGT","forward","COI"),"IUPAC")
expect_error(primer_mixer_custom("tiny","ACGT","forward","COI"),"8 to 120")
expect_error(primer_mixer_build(catalog,f$primer_id,r$primer_id,templates,"COI"),"chosen marker")
expect_error(primer_mixer_build(catalog,f$primer_id,r$primer_id,templates,"12S_MT",max_mismatches=6L),"0 to 5")
expect_error(primer_mixer_build(catalog,rep(f$primer_id,6L),r$primer_id,templates,"12S_MT",max_total_bases=3L),"total template")

fasta <- primer_mixer_parse_fasta(paste0(">One\n",raw,"\n>Two\n",raw),"12S_MT")
stopifnot(nrow(fasta)==2L,all(fasta$marker_id=="12S_MT"),identical(fasta$sequence,c(raw,raw)))
expect_error(primer_mixer_parse_fasta(">repeat\nACGT\n>repeat\nACGT","COI"),"unique")
expect_error(primer_mixer_parse_fasta("ACGT","COI"),"FASTA")
expect_error(primer_mixer_parse_fasta(">gap\nACGT-ACGT","COI"),"ungapped")

# The shared map adapter reports actual raw bp from one selected template.
map <- primer_mixer_overview_data(result,"fish_1")
stopifnot(nrow(map$rows)==1L,nrow(map$primers)==2L,map$rows$pair_start==18L,
  map$rows$pair_end==141L,map$rows$aligned_amplicon_bp==124L,
  map$rows$targeted_region_bp==83L,map$marker$coordinate_reference=="fish_1",
  map$reference_length==nchar(raw),map$marker$marker_id=="12S_MT",map$rows$platform_tier=="ILLUMINA")
stopifnot(identical(primer_mixer_platform_tier(c(500L,501L,2000L,2001L,3500L,3501L,NA_integer_)),
  c("ILLUMINA","PACBIO","PACBIO","NANOPORE","NANOPORE","BEYOND",NA_character_)))
expect_error(primer_mixer_platform_tier(500.5),"inclusive integer")
# Planning colors follow measured inclusive products on different molecules,
# including the exact tier boundaries; internal spans cannot set these tiers.
planned_bp <- c(500L,501L,2000L,2001L,3500L,3501L)
planning_templates <- data.frame(template_id=paste0("planning_",planned_bp),
  sequence=vapply(planned_bp,function(bp) paste0(left,forward,dna(bp-nchar(forward)-nchar(reverse)),.primer_mixer_complement(reverse),right),character(1)),
  marker_id="12S_MT",stringsAsFactors=FALSE)
planning <- primer_mixer_build(catalog,f$primer_id,r$primer_id,planning_templates,"12S_MT")
stopifnot(identical(planning$products$amplicon_bp,planned_bp))
actual_tiers <- vapply(planning_templates$template_id,function(id) primer_mixer_overview_data(planning,id)$rows$platform_tier,character(1))
stopifnot(identical(unname(actual_tiers),c("ILLUMINA","PACBIO","PACBIO","NANOPORE","NANOPORE","BEYOND")))
filtered_map <- primer_mixer_overview_data(cart,"fish_1")
stopifnot(nrow(filtered_map$rows)==0L,nrow(cart$pairs)==4L)
scanned <- primer_mixer_sites(catalog,templates,"12S_MT")
stopifnot(nrow(scanned)==2L,all(scanned$status=="exact_site"),identical(scanned$start,c(18L,121L)))
expect_error(primer_mixer_sites(catalog,templates,"12S_MT",max_search_bases=2L),"Narrow")

library <- primer_mixer_catalog()
stopifnot(!anyDuplicated(library$primer_id),all(c("COI","18S","ITS_FUNGAL") %in% library$marker_id),
  all(library$direction %in% c("forward","reverse")),all(nzchar(library$source_sequence)),
  nrow(primer_mixer_templates(marker_id="COI",target_group="ARTHROPODS"))==1L,
  !any(primer_mixer_templates(marker_id="12S_MT",target_group="ARTHROPODS")$marker_id=="12S_MT"))
pr2_source <- utils::read.csv("data/catalog/pr2_18s_primers.csv",stringsAsFactors=FALSE)
stopifnot(all(paste0("pr2_18s_",pr2_source$primer_id) %in% library$primer_id))
fish_templates <- primer_mixer_templates(target_group="FISH")
if(nrow(fish_templates)) stopifnot(is.logical(fish_templates$partial),
  all(fish_templates$partial[fish_templates$accession_version %in% c("FJ915075.1","FJ710886.1")]))
message("Generic primer mixer tests passed: Cartesian candidates, measured coordinates, mismatches, 3-prime orientation, ambiguity, marker isolation and source identity.")
