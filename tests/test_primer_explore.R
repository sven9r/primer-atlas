# The shared Shiny investigation flow must keep marker/group context, the
# selected chemical oligos, raw template evidence and display filters together.
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(dplyr))
`%||%` <- function(x,y) if(is.null(x)) y else x
source("R/primer_mixer.R")
source("R/primer_mixer_exports.R")
source("R/primer_site_overview.R")
source("R/primer_overview.R")
source("R/animal_fit.R")
source("R/primer_explore_server.R")

chosen_marker <- shiny::reactiveVal("12S_MT")
chosen_group <- shiny::reactiveVal("FISH")
shiny::testServer(primerExploreServer,args=list(marker=chosen_marker,target_group=chosen_group),{
  session$setInputs(platform="ALL",mismatches="0",min_bp=0,max_bp=10000,region="ALL",width=620)
  session$flushReact()
  stopifnot(active_marker()=="12S_MT",nrow(templates())==3L,
    all(templates()$marker_id=="12S_MT"),all(templates()$target_group=="FISH"),
    !length(forward()),!length(reverse()))
  acknowledge <- function() {
    session$setInputs(forward=shiny::isolate(state$forward),reverse=shiny::isolate(state$reverse))
    session$flushReact()
  }
  acknowledge()
  session$setInputs(example=1L)
  session$flushReact()
  acknowledge()
  initial <- result()
  stopifnot(length(forward())==4L,length(reverse())==4L,nrow(initial$pairs)==16L,
    nrow(initial$products)==48L,all(initial$summary$n_templates==3L),
    all(initial$summary$n_exact_products<=initial$summary$n_templates))
  # Exercise actual registered download callbacks through the public mock
  # output API. Tables/figures must carry the same input sequences and sources.
  result_file <- output$download_results
  primer_file <- output$download_primers
  template_file <- output$download_templates
  figure_file <- output$download_figure
  table_file <- output$download_table
  package_file <- output$download_package
  stopifnot(all(file.exists(c(result_file,primer_file,template_file,figure_file,table_file,package_file))))
  exported <- utils::read.csv(result_file,stringsAsFactors=FALSE,check.names=FALSE)
  oligos <- utils::read.csv(primer_file,stringsAsFactors=FALSE,check.names=FALSE)
  stopifnot(nrow(exported)==48L,nrow(oligos)==8L,
    setequal(exported$pair_id,initial$pairs$pair_id),
    all(c("forward_primer_id","reverse_primer_id","forward_citation","reverse_citation",
      "minimum_product_bp","maximum_product_bp","platform_ceiling_bp") %in% names(exported)),
    all(exported$allowed_substitutions_per_primer==0L),
    identical(exported$forward_sequence,initial$pairs$forward_sequence[match(exported$pair_id,initial$pairs$pair_id)]),
    all(nzchar(oligos$citation)),all(nzchar(oligos$source_url)))
  expected_hashes <- vapply(templates()$sequence,function(x) digest::digest(x,algo="sha256",serialize=FALSE),character(1))
  stopifnot(identical(exported$template_sha256,unname(expected_hashes[match(exported$template_id,templates()$template_id)])))
  template_lines <- readLines(template_file,warn=FALSE)
  stopifnot(sum(startsWith(template_lines,">"))==3L,
    identical(template_lines[!startsWith(template_lines,">")],templates()$sequence))
  svg <- paste(readLines(figure_file,warn=FALSE),collapse="\n")
  raw_id <- selected_reference()
  mapped_source <- initial$components$source_url[initial$components$primer_name=="12S-V5f"][1L]
  stopifnot(grepl("<svg",svg,fixed=TRUE),grepl(raw_id,svg,fixed=TRUE),
    grepl(expected_hashes[match(raw_id,templates()$template_id)],svg,fixed=TRUE),
    grepl(mapped_source,svg,fixed=TRUE),grepl("PCR",svg,fixed=TRUE),grepl("candidate",svg,ignore.case=TRUE))
  table_html <- paste(readLines(table_file,warn=FALSE),collapse="\n")
  table_xml <- xml2::read_html(table_html)
  stopifnot(length(xml2::xml_find_all(table_xml,"//table[contains(@class,'combination-table')]/tbody/tr"))==16L,
    grepl(mapped_source,table_html,fixed=TRUE),grepl("@page{size:A4 landscape",table_html,fixed=TRUE))
  package_directory <- tempfile("module-evidence-")
  dir.create(package_directory)
  files <- utils::unzip(package_file,list=TRUE)$Name
  stopifnot(all(c("template-results.csv","primer-components.csv","primer-sites.csv","source-template-metadata.csv",
    "templates.fasta","paper-table.html","amplicons.svg","methods.json","METHODS.txt","SHA256SUMS") %in% files),
    !"taxonomic-fit.csv" %in% files)
  utils::unzip(package_file,exdir=package_directory)
  checksums <- readLines(file.path(package_directory,"SHA256SUMS"),warn=FALSE)
  names <- sub("^[a-f0-9]{64}[[:space:]]+","",checksums)
  hashes <- substr(checksums,1L,64L)
  stopifnot(setequal(names,setdiff(files,"SHA256SUMS")),all(file.exists(file.path(package_directory,names))))
  observed <- vapply(file.path(package_directory,names),digest::digest,character(1),algo="sha256",file=TRUE)
  stopifnot(identical(unname(observed),hashes),
    nrow(utils::read.csv(file.path(package_directory,"template-results.csv")))==48L)
  methods <- jsonlite::read_json(file.path(package_directory,"methods.json"),simplifyVector=TRUE)
  stopifnot(methods$n_templates==3L,nrow(methods$requested_combinations)==16L)
  unlink(package_directory,recursive=TRUE)
  exact0 <- sum(initial$products$status=="exact_product")
  session$setInputs(mismatches="1")
  session$flushReact()
  stopifnot(sum(result()$products$status=="exact_product")==exact0,
    sum(result()$products$status=="mismatch_product")>0L,
    result()$settings$max_mismatches==1L)
  requested_ids <- result()$pairs$pair_id
  session$setInputs(platform="ILLUMINA",max_bp=150)
  session$flushReact()
  stopifnot(identical(result()$pairs$pair_id,requested_ids),result()$settings$max_bp==150,
    all(result()$products$amplicon_bp[result()$products$passes_length_filter]<=150L))
  session$setInputs(platform="ALL",max_bp=10000,mismatches="0")
  session$flushReact()
  # A platform ceiling below the requested minimum removes mapped products
  # without deleting candidates or pretending the underlying exact sites failed.
  session$setInputs(platform="ILLUMINA",min_bp=600)
  session$flushReact()
  stopifnot(nrow(result()$pairs)==16L,!any(result()$products$passes_length_filter),
    sum(result()$products$status=="exact_product")==exact0)
  session$setInputs(platform="ALL",min_bp=0)
  session$flushReact()

  # Custom oligos combine with curated components without acquiring citations
  # or changing any original catalog row.
  template_sequence <- templates()$sequence[1L]
  custom_sequence <- substr(template_sequence,1L,20L)
  before_catalog <- base_catalog
  session$setInputs(custom_name="My additional forward",custom_sequence=custom_sequence,
    custom_direction="forward",add_custom=1L)
  session$flushReact()
  acknowledge()
  stopifnot(length(forward())==5L,length(reverse())==4L,nrow(result()$pairs)==20L,
    any(result()$components$catalog_origin=="custom_session_input"),identical(before_catalog,base_catalog))
  added <- state$custom[state$custom$primer_name=="My additional forward",,drop=FALSE]
  stopifnot(nrow(added)==1L,added$source_key=="custom",!nzchar(added$source_url))

  # Two explicit FASTA templates replace the starter panel only in this
  # marker/group session; header names do not become inferred taxonomy.
  fasta <- paste0(">input_A\n",template_sequence,"\n>input_B\n",templates()$sequence[2L])
  session$setInputs(fasta=fasta,use_fasta=1L)
  session$flushReact()
  stopifnot(nrow(templates())==2L,all(templates()$version=="session_input"),
    nrow(result()$products)==40L,all(result()$summary$n_templates==2L),
    identical(templates()$template_id,c("input_A","input_B")))
  saved_forward <- forward(); saved_reverse <- reverse()
  chosen_marker("COI")
  session$flushReact()
  stopifnot(active_marker()=="COI",!length(forward()),!length(reverse()),
    all(templates()$marker_id=="COI"),all(templates()$target_group=="FISH"),nrow(templates())==3L)
  chosen_marker("12S_MT")
  session$flushReact()
  acknowledge()
  stopifnot(setequal(forward(),saved_forward),setequal(reverse(),saved_reverse),
    identical(templates()$template_id,c("input_A","input_B")),nrow(result()$pairs)==20L)
  session$setInputs(reset_templates=1L)
  session$flushReact()
  stopifnot(nrow(templates())==3L,all(templates()$target_group=="FISH"))

  # The visible derived arrow must be the third click in the three-state
  # cycle. It cannot create a chain of nested reverse-complement children.
  session$setInputs(clear=1L)
  session$flushReact()
  acknowledge()
  candidate <- catalog()[catalog()$marker_id=="12S_MT" & catalog()$direction=="forward",,drop=FALSE][1L,]
  parent_id <- candidate$primer_id
  session$setInputs(primer_click=list(id=parent_id,nonce=1L))
  session$flushReact()
  acknowledge()
  stopifnot(identical(forward(),parent_id),!length(reverse()))
  session$setInputs(primer_click=list(id=parent_id,nonce=2L))
  session$flushReact()
  acknowledge()
  stopifnot(!length(forward()),length(reverse())==1L)
  child_id <- reverse()
  child <- catalog()[match(child_id,catalog()$primer_id),,drop=FALSE]
  stopifnot(child$catalog_origin=="derived_reverse_complement",child$parent_primer_id==parent_id,
    child$sequence==.primer_mixer_complement(candidate$sequence),child$direction=="reverse",
    candidate$primer_id %in% base_catalog$primer_id)
  session$setInputs(primer_click=list(id=child_id,nonce=3L))
  session$flushReact()
  acknowledge()
  stopifnot(!length(forward()),!length(reverse()),!any(startsWith(catalog()$primer_id,"rc_rc_")))

  # Target-group changes cannot transfer the fish cohort into arthropods.
  chosen_group("ARTHROPODS")
  chosen_marker("COI")
  session$flushReact()
  stopifnot(all(templates()$target_group=="ARTHROPODS"),all(templates()$marker_id=="COI"),
    !any(templates()$template_id %in% primer_mixer_templates(marker_id="COI",target_group="FISH")$template_id))
  for(next_marker in c("18S","28S")) {
    chosen_marker(next_marker)
    session$flushReact()
    acknowledge()
    session$setInputs(example=if(next_marker=="18S") 10L else 11L)
    session$flushReact()
    acknowledge()
    stopifnot(nrow(templates())>=1L,all(templates()$target_group=="ARTHROPODS"),
      all(templates()$marker_id==next_marker),nrow(result()$pairs)>=1L,
      nrow(result()$products)==nrow(result()$pairs)*nrow(templates()),
      all(result()$summary$n_templates==nrow(templates())))
  }
  chosen_group("BIRDS")
  chosen_marker("COI")
  session$flushReact()
  acknowledge()
  session$setInputs(example=12L)
  session$flushReact()
  acknowledge()
  selected_names <- result()$components$primer_name
  stopifnot(all(c("BirdF1","BirdR1") %in% selected_names),all(templates()$target_group=="BIRDS"),
    all(templates()$marker_id=="COI"),nrow(templates())>=1L)
  for(next_marker in c("12S_MT","18S","28S")) {
    chosen_marker(next_marker)
    session$flushReact()
    stopifnot(nrow(templates())>=1L,all(templates()$target_group=="BIRDS"),
      all(templates()$marker_id==next_marker),!any(grepl("^FISH_",templates()$template_id)))
  }
})

# Optional navigation callbacks receive the actual selected oligo identities
# and remain quiet for an empty selection or an unrelated version control.
callback_marker <- shiny::reactiveVal("COI")
callback_group <- shiny::reactiveVal("ARTHROPODS")
callback_report <- shiny::reactiveVal(NULL)
evidence_calls <- list(); version_calls <- character()
shiny::testServer(primerExploreServer,args=list(marker=callback_marker,target_group=callback_group,
  taxonomic_report=callback_report,
  on_evidence=function(pair,components) evidence_calls[[length(evidence_calls)+1L]] <<- list(pair=pair,components=components),
  on_version=function(version) version_calls <<- c(version_calls,version)),{
  session$setInputs(platform="ALL",mismatches="0",min_bp=0,max_bp=10000,region="ALL",width=620)
  session$flushReact()
  session$setInputs(broader_evidence=1L,snapshot_version="unrelated-marker-version")
  session$flushReact()
  stopifnot(!length(evidence_calls),!length(version_calls))
  session$setInputs(example=1L)
  session$flushReact()
  before_pairs <- result()$pairs
  session$setInputs(broader_evidence=2L)
  session$flushReact()
  stopifnot(length(evidence_calls)==1L,nrow(evidence_calls[[1L]]$pair)==1L,
    evidence_calls[[1L]]$pair$pair_id %in% before_pairs$pair_id,
    all(c(evidence_calls[[1L]]$pair$forward_id,evidence_calls[[1L]]$pair$reverse_id) %in% evidence_calls[[1L]]$components$primer_id),
    identical(before_pairs,result()$pairs),!length(version_calls))
  # The common package preserves explicit user taxonomy separately from the
  # original NCBI source-template metadata and includes both taxonomic tables.
  ids <- result()$templates$template_id
  taxonomy_csv <- paste0("template_id,family\n",ids[1L],",My reviewed & supplied family\n")
  user_rows <- animal_fit_parse_taxonomy_csv(taxonomy_csv,ids)
  callback_report(animal_fit_data(result(),"family",user_taxonomy=user_rows))
  session$flushReact()
  annotated_package <- output$download_package
  annotated_html <- paste(readLines(output$download_table,warn=FALSE),collapse="\n")
  stopifnot(grepl("My reviewed &amp; supplied family",annotated_html,fixed=TRUE),
    grepl("user_provided_metadata",annotated_html,fixed=TRUE))
  annotated_directory <- tempfile("annotated-evidence-")
  dir.create(annotated_directory)
  annotated_files <- utils::unzip(annotated_package,list=TRUE)$Name
  stopifnot(all(c("taxonomic-fit.csv","taxonomic-template-evidence.csv") %in% annotated_files))
  utils::unzip(annotated_package,exdir=annotated_directory)
  original_metadata <- utils::read.csv(file.path(annotated_directory,"source-template-metadata.csv"),stringsAsFactors=FALSE)
  annotated_metadata <- utils::read.csv(file.path(annotated_directory,"taxonomic-template-evidence.csv"),stringsAsFactors=FALSE)
  stopifnot(all(original_metadata$taxonomy_status=="source_taxid_verified"),
    all(annotated_metadata$taxonomy_status=="user_provided_metadata"),
    all(annotated_metadata$taxonomy_source_sha256==user_rows$taxonomy_source_sha256),
    all(annotated_metadata$taxonomy_version=="session_metadata"))
  unlink(annotated_directory,recursive=TRUE)
})
message("Shared primer exploration module tests passed:16 combinations,template denominators,custom components,marker context,raw FASTA and visible RC click cycle.")
