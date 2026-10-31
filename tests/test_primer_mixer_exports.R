suppressPackageStartupMessages(library(dplyr))
source("R/primer_mixer.R")
source("R/primer_mixer_exports.R")
source("R/primer_overview.R")

catalog <- primer_mixer_catalog()
templates <- primer_mixer_templates(marker_id="28S",target_group="FISH")
catalog <- catalog[catalog$marker_id=="28S" & catalog$catalog_origin=="fish_assay_catalog",,drop=FALSE]
result <- primer_mixer_build(catalog,catalog$primer_id[catalog$direction=="forward"],
  catalog$primer_id[catalog$direction=="reverse"],templates,"28S")
result$settings$platform_ceiling_bp <- 2000
stopifnot(nrow(result$pairs)==1L,result$products$status=="exact_product")
data <- primer_mixer_overview_data(result,templates$template_id[1L])
map <- primer_overview_map(data$rows,data$primers,data$marker,citations=data$citations,
  focus=data$focus,domain=data$domain,viewport_width=NULL)
svg <- primer_mixer_export_svg(map,result,templates$template_id[1L],domain=data$domain)
document <- xml2::read_xml(svg)
ns <- c(s="http://www.w3.org/2000/svg")
metadata_node <- xml2::xml_find_first(document,"//s:metadata[@id='primer-atlas-mixer-metadata']",ns)
stopifnot(!inherits(metadata_node,"xml_missing"))
metadata <- jsonlite::fromJSON(xml2::xml_text(metadata_node),simplifyVector=FALSE)
sha <- digest::digest(templates$sequence[1L],algo="sha256",serialize=FALSE)
stopifnot(metadata$schema_version==1L,metadata$marker_id=="28S",
  metadata$displayed_template_sha256==sha,metadata$settings$platform_ceiling_bp==2000,
  metadata$n_requested_combinations==1L,metadata$n_input_templates==1L,
  metadata$displayed_template[[1L]]$template_id==templates$template_id[1L],
  metadata$displayed_template[[1L]]$sequence==templates$sequence[1L],
  metadata$components[[1L]]$primer_id==result$components$primer_id[1L],
  grepl("PCR probability",metadata$evidence_scope,fixed=TRUE),
  !grepl("amplification_probability",svg,fixed=TRUE))
footer <- xml2::xml_find_first(document,"//s:g[@class='primer-mixer-export-provenance']",ns)
footer_text <- xml2::xml_text(footer)
stopifnot(grepl("PREDICTED RefSeq MODEL",footer_text,fixed=TRUE),
  grepl(result$components$citation[1L],footer_text,fixed=TRUE),
  grepl(result$components$source_url[1L],footer_text,fixed=TRUE),
  grepl(sha,footer_text,fixed=TRUE),grepl("Candidate combinations",footer_text,fixed=TRUE))
vb <- as.numeric(strsplit(xml2::xml_attr(document,"viewBox")," ",fixed=TRUE)[[1L]])
ys <- as.numeric(xml2::xml_attr(xml2::xml_find_all(footer,".//s:text",ns),"y"))
stopifnot(max(ys)<vb[4L],xml2::xml_attr(document,"width")=="1680",as.numeric(xml2::xml_attr(document,"height"))==vb[4L])

table <- primer_mixer_results_table(result)
stopifnot(nrow(table)==1L,table$forward_primer_id==result$pairs$forward_id,
  table$reverse_primer_id==result$pairs$reverse_id,table$template_sha256==sha,
  table$template_source==templates$source_url,table$template_partial==templates$partial,
  table$forward_citation==result$components$citation[result$components$direction=="forward"],
  table$reverse_citation==result$components$citation[result$components$direction=="reverse"],
  table$minimum_product_bp==0,table$maximum_product_bp==Inf,table$platform_ceiling_bp==2000,
  table$n_input_templates==1L,grepl("no PCR probability",table$evidence_scope,fixed=TRUE))

# Derived/custom chemicals retain parent provenance in the export; XML
# escaping and literal backslashes in user names must survive unchanged.
parent <- primer_mixer_custom("Parent <A> \\3", "ACGTACGTTGCATGCA", "forward", "COI")
child <- primer_mixer_reverse(parent)
other <- primer_mixer_custom("Other & forward", "TGTACGACCGTATCAG", "forward", "COI")
raw <- paste0(other$sequence,"TCGACGTCAGATCGTTGCATGACATGCA",parent$sequence)
custom_templates <- data.frame(template_id="Custom <template> & 1",sequence=raw,marker_id="COI",version="session_input",source_url="",stringsAsFactors=FALSE)
custom <- primer_mixer_build(.primer_mixer_bind(list(child,other)),other$primer_id,child$primer_id,custom_templates,"COI")
custom$settings$platform_ceiling_bp <- Inf
custom_data <- primer_mixer_overview_data(custom,custom_templates$template_id)
custom_map <- primer_overview_map(custom_data$rows,custom_data$primers,custom_data$marker,
  citations=custom_data$citations,focus=custom_data$focus,domain=custom_data$domain,viewport_width=360)
custom_svg <- primer_mixer_export_svg(custom_map,custom,custom_templates$template_id)
custom_xml <- xml2::read_xml(custom_svg)
custom_metadata <- jsonlite::fromJSON(xml2::xml_text(xml2::xml_find_first(custom_xml,"//s:metadata",ns)),simplifyVector=FALSE)
stopifnot(custom_metadata$settings$platform_ceiling_bp=="unbounded",
  custom_metadata$displayed_template[[1L]]$template_id=="Custom <template> & 1",
  grepl("Parent <A> \\3",xml2::xml_text(custom_xml),fixed=TRUE))
custom_csv <- primer_mixer_results_table(custom)
stopifnot(custom_csv$reverse_parent_primer_id==parent$primer_id,custom_csv$reverse_parent_sequence==parent$source_sequence,
  grepl("no empirical validation",custom_csv$reverse_normalization_note,fixed=TRUE))

empty <- primer_mixer_build(catalog,character(),character(),templates,"28S")
empty_table <- primer_mixer_results_table(empty)
stopifnot(nrow(empty_table)==0L,all(c("forward_citation","template_sha256","platform_ceiling_bp") %in% names(empty_table)))

# The formatted table retains every requested pairing, both source-oriented
# oligos, explicit counts and source/model boundaries in a printable artifact.
html <- primer_mixer_export_html(result)
html_xml <- xml2::read_html(html)
combination_rows <- xml2::xml_find_all(html_xml,"//table[contains(@class,'combination-table')]/tbody/tr")
stopifnot(length(combination_rows)==1L,grepl("@page{size:A4 landscape",html,fixed=TRUE),
  grepl(result$components$citation[1L],xml2::xml_text(html_xml),fixed=TRUE),
  grepl(result$components$source_sequence[1L],xml2::xml_text(html_xml),fixed=TRUE),
  grepl("PREDICTED RefSeq MODEL",xml2::xml_text(html_xml),fixed=TRUE),
  grepl("no PCR probability",xml2::xml_text(html_xml),fixed=TRUE),
  !length(xml2::xml_find_all(html_xml,"//script")))
escaped_html <- primer_mixer_export_html(custom)
escaped_xml <- xml2::read_html(escaped_html)
stopifnot(grepl("Parent <A> \\3",xml2::xml_text(escaped_xml),fixed=TRUE),
  !length(xml2::xml_find_all(escaped_xml,"//a[contains(@href,'javascript:')]")))

`%||%` <- function(x,y) if(is.null(x)) y else x
source("R/animal_fit.R")
all_catalog <- primer_mixer_catalog()
fish_12s <- all_catalog[all_catalog$marker_id=="12S_MT" & all_catalog$catalog_origin=="fish_assay_catalog",,drop=FALSE]
panel_12s <- primer_mixer_templates(marker_id="12S_MT",target_group="FISH")
many <- primer_mixer_build(all_catalog,fish_12s$primer_id[fish_12s$direction=="forward"],fish_12s$primer_id[fish_12s$direction=="reverse"],panel_12s,"12S_MT",max_mismatches=2L)
user_csv <- paste0("template_id,family\n",panel_12s$template_id[1L],",User & supplied family\n")
user_taxonomy <- animal_fit_parse_taxonomy_csv(user_csv,panel_12s$template_id)
tax_report <- animal_fit_data(many,"family",user_taxonomy=user_taxonomy)
many_html <- primer_mixer_export_html(many,tax_report)
many_xml <- xml2::read_html(many_html)
stopifnot(length(xml2::xml_find_all(many_xml,"//table[contains(@class,'combination-table')]/tbody/tr"))==16L,
  length(xml2::xml_find_all(many_xml,"//table[contains(@class,'taxonomy-table')]/tbody/tr"))==nrow(tax_report$summary),
  grepl("User & supplied family",xml2::xml_text(many_xml),fixed=TRUE),
  grepl("user_provided_metadata",xml2::xml_text(many_xml),fixed=TRUE),
  grepl(user_taxonomy$taxonomy_source_sha256,xml2::xml_text(many_xml),fixed=TRUE),
  grepl("session_metadata",xml2::xml_text(many_xml),fixed=TRUE))
bad_tax <- tax_report; bad_tax$metadata$marker_id <- "COI"
error <- tryCatch({primer_mixer_export_html(many,bad_tax);NULL},error=function(e)e)
stopifnot(inherits(error,"error"),grepl("marker differs",conditionMessage(error)))
stopifnot(grepl("No forward × reverse combinations",primer_mixer_export_html(empty),fixed=TRUE))
message("Primer mixer exports passed:actual template/source metadata,visible citations,source flags,settings,raw sequence hashes,derived oligos and portable SVG bounds.")
