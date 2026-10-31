# Portable evidence exports for the session-local primer mixer. These helpers
# describe the supplied sequences and candidate combinations; they never
# inherit experimental validation from a component's parent publication.

.primer_mixer_export_field <- function(rows,field,default="") {
  if(!field %in% names(rows)) return(rep(default,nrow(rows)))
  value <- rows[[field]]
  if(is.character(value)) value[is.na(value)] <- default
  value
}

.primer_mixer_export_scope <- paste(
  "Candidate combinations evaluated on explicit input sequences.",
  "Sequence compatibility does not establish PCR success; no PCR probability",
  "or experimental validation of a new pairing is inferred.",
  "Parent citations describe the original oligos and do not validate derived reverse-complement candidates."
)

primer_mixer_results_table <- function(result) {
  required <- c("products","pairs","components","templates","summary","settings","method")
  if(!all(required %in% names(result))) stop("Export requires a complete primer mixer result.")
  rows <- result$products
  pi <- match(rows$pair_id,result$pairs$pair_id)
  ti <- match(rows$template_id,result$templates$template_id)
  pairs <- result$pairs[pi,,drop=FALSE]
  template <- result$templates[ti,,drop=FALSE]
  rows$pair_label <- pairs$pair_label
  for(role in c("forward","reverse")) {
    component <- result$components[match(pairs[[paste0(role,"_id")]],result$components$primer_id),,drop=FALSE]
    mapping <- c(primer_id="primer_id",primer_name="primer_name",sequence="sequence",source_sequence="source_sequence",
      direction="direction",citation="citation",source_url="source_url",catalog_origin="catalog_origin",target_group="target_group",
      target_scope="target_scope",parent_primer_id="parent_primer_id",parent_primer_name="parent_primer_name",parent_sequence="parent_sequence",
      parent_direction="parent_direction",original_primer_id="original_primer_id",original_marker_id="original_marker_id",
      original_gene_locus="original_gene_locus",gene_locus="gene_locus",evidence_note="evidence_note",normalization_note="normalization_note")
    for(name in names(mapping)) rows[[paste0(role,"_",name)]] <- .primer_mixer_export_field(component,mapping[[name]])
  }
  ranks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")
  for(field in c("version","source_url","organism","accession_version","source_start","source_end","source_strand","source_feature_location",
      "source_status","partial","projected","boundary_status","evidence_status","source_landmark_id","source_landmark_type",
      "source_note","mapping_method","parent_sequence_sha256","target_group","normalization_note",ranks,"taxid","taxonomy_status",
      "taxonomy_source_url","taxonomy_version","taxonomy_source_sha256")) {
    name <- if(field=="source_url") "template_source" else paste0("template_",field)
    rows[[name]] <- .primer_mixer_export_field(template,field,default=if(field %in% ranks) "Unresolved" else "")
  }
  rows$planning_platform_tier <- primer_mixer_platform_tier(rows$amplicon_bp)
  rows$template_sha256 <- vapply(template$sequence,function(sequence) {
    if(is.na(sequence) || !nzchar(sequence)) NA_character_ else digest::digest(sequence,algo="sha256",serialize=FALSE)
  },character(1))
  rows$allowed_substitutions_per_primer <- rep(result$settings$max_mismatches,nrow(rows))
  rows$minimum_product_bp <- rep(result$settings$min_bp,nrow(rows))
  rows$maximum_product_bp <- rep(result$settings$max_bp,nrow(rows))
  rows$platform_ceiling_bp <- rep(if(is.null(result$settings$platform_ceiling_bp)) Inf else result$settings$platform_ceiling_bp,nrow(rows))
  si <- match(rows$pair_id,result$summary$pair_id)
  rows$n_input_templates <- result$summary$n_templates[si]
  rows$n_matching_marker_templates <- result$summary$n_matching_marker_templates[si]
  rows$coordinate_model <- rep("One-based inclusive positions on the selected ungapped raw template; internal bp excludes both oligos",nrow(rows))
  rows$length_filter_basis <- rep("Inclusive product bp, including both oligos; filter does not alter site evidence or the template denominator",nrow(rows))
  rows$evidence_scope <- rep(.primer_mixer_export_scope,nrow(rows))
  rows$search_method <- rep(result$method,nrow(rows))
  rows
}

.primer_mixer_export_wrap <- function(text,width) {
  width <- max(16L,as.integer(width))
  words <- strsplit(trimws(gsub("[[:space:]]+"," ",text))," ",fixed=TRUE)[[1L]]
  if(!length(words)) return("")
  broken <- unlist(lapply(words,function(word) {
    if(nchar(word)<=width) return(word)
    starts <- seq.int(1L,nchar(word),by=width)
    substring(word,starts,pmin(nchar(word),starts+width-1L))
  }),use.names=FALSE)
  output <- character(); line <- ""
  for(word in broken) {
    candidate <- if(nzchar(line)) paste(line,word) else word
    if(nchar(candidate)>width && nzchar(line)) { output <- c(output,line); line <- word }
    else line <- candidate
  }
  c(output,line)
}

primer_mixer_export_svg <- function(map_tag,result,template_id,domain=NULL) {
  ti <- match(template_id,result$templates$template_id)
  if(length(ti)!=1L || is.na(ti)) stop("Export one explicit mixer template at a time.")
  template <- result$templates[ti,,drop=FALSE]
  if(is.na(template$sequence[1L]) || !nzchar(template$sequence[1L])) stop("A raw sequence is required for a map export.")
  html <- as.character(map_tag)
  svg <- regmatches(html,regexpr("<svg\\b[\\s\\S]*?</svg>",html,perl=TRUE))
  if(!length(svg) || !nzchar(svg)) stop("Choose a template with a mapped product before exporting the figure.")
  open <- regmatches(svg,regexpr("^<svg\\b[^>]*>",svg,perl=TRUE))
  parsed <- regmatches(open,regexec('viewBox="([^"]+)"',open,perl=TRUE))[[1L]]
  if(length(parsed)!=2L) stop("The map SVG needs a finite viewBox for a portable export.")
  bounds <- suppressWarnings(as.numeric(strsplit(parsed[2L],"[[:space:]]+")[[1L]]))
  if(length(bounds)!=4L || any(!is.finite(bounds)) || any(bounds[3:4]<=0)) stop("The map SVG viewBox is invalid.")
  width <- bounds[3L]; original_height <- bounds[4L]
  if(is.null(domain)) domain <- c(1L,nchar(template$sequence[1L]))
  if(length(domain)!=2L || any(!is.finite(domain)) || domain[1L]<1L || domain[2L]<domain[1L] || domain[2L]>nchar(template$sequence[1L])) stop("The export domain must be an inclusive window on this template.")
  products <- primer_mixer_results_table(result)
  displayed <- products[products$template_id==template_id & products$passes_length_filter,,drop=FALSE]
  hashes <- vapply(result$templates$sequence,function(sequence) {
    if(is.na(sequence) || !nzchar(sequence)) NA_character_ else digest::digest(sequence,algo="sha256",serialize=FALSE)
  },character(1))
  template_panel <- result$templates[,setdiff(names(result$templates),c("sequence","source_sequence")),drop=FALSE]
  template_panel$sequence_sha256 <- hashes
  settings <- lapply(result$settings,function(value) {
    if(!is.numeric(value) || all(is.finite(value))) return(value)
    if(length(value)==1L && is.infinite(value)) return("unbounded")
    lapply(value,function(one) if(is.infinite(one)) "unbounded" else one)
  })
  metadata <- list(schema_version=1L,artifact_type="primer_atlas_candidate_amplicon_map",
    marker_id=result$marker_id,coordinate_model="One-based inclusive positions on one ungapped raw input template",
    display_domain=c(start=domain[1L],end=domain[2L]),displayed_template=template,
    displayed_template_sha256=hashes[ti],template_panel=template_panel,
    components=result$components,requested_combinations=result$pairs,displayed_products=displayed,
    n_requested_combinations=nrow(result$pairs),n_input_templates=nrow(result$templates),
    settings=settings,search_method=result$method,evidence_scope=.primer_mixer_export_scope)
  metadata_json <- jsonlite::toJSON(metadata,auto_unbox=TRUE,pretty=TRUE,na="null",null="null",dataframe="rows",digits=NA)
  cite <- unique(result$components[,intersect(c("citation","source_url"),names(result$components)),drop=FALSE])
  sources <- if(nrow(cite)) vapply(seq_len(nrow(cite)),function(i) paste(
    .primer_mixer_export_field(cite,"citation")[i],.primer_mixer_export_field(cite,"source_url")[i]),character(1)) else "Custom oligos; no published validation assigned."
  record_status <- .primer_mixer_export_field(template,"source_status")[1L]
  partial <- if("partial" %in% names(template) && tolower(as.character(template$partial[1L])) %in% c("true","1","yes")) " · partial source record" else ""
  limits <- paste0("Allowed substitutions per primer: ",result$settings$max_mismatches,
    "; product bp range: ",result$settings$min_bp,"–",if(is.finite(result$settings$max_bp)) result$settings$max_bp else "unbounded",
    if(!is.null(result$settings$platform_ceiling_bp)) paste0("; platform ceiling: ",if(is.finite(result$settings$platform_ceiling_bp)) result$settings$platform_ceiling_bp else "unbounded") else "",
    ". Length includes both primers.")
  caption <- c(paste0("Primer Atlas · ",result$marker_id," · ",template_id," · raw bp ",domain[1L],"–",domain[2L]),
    paste0(nrow(displayed)," mapped products on this template; ",nrow(result$pairs)," requested combinations across ",nrow(result$templates)," input templates."),
    .primer_mixer_export_scope,limits,
    paste0("Template: ",.primer_mixer_export_field(template,"version")[1L]," · ",record_status,partial," · ",.primer_mixer_export_field(template,"source_url")[1L]),
    paste0("Template SHA-256: ",hashes[ti]),paste0("Oligo source: ",sources))
  font_size <- 12
  max_chars <- max(16L,floor((width-48)/(font_size*.58)))
  lines <- unlist(lapply(caption,.primer_mixer_export_wrap,width=max_chars),use.names=FALSE)
  footer_height <- 36+length(lines)*18
  new_height <- original_height+footer_height
  updated_open <- sub('viewBox="[^"]+"',paste0('viewBox="',paste(c(bounds[1:3],new_height),collapse=" "),'"'),open,perl=TRUE)
  if(!grepl(' xmlns="',updated_open,fixed=TRUE)) updated_open <- sub("^<svg",'<svg xmlns="http://www.w3.org/2000/svg"',updated_open)
  if(grepl(' width="',updated_open,fixed=TRUE)) updated_open <- sub(' width="[^"]+"',paste0(' width="',width,'"'),updated_open,perl=TRUE)
  else updated_open <- sub(">$",paste0(' width="',width,'">'),updated_open)
  if(grepl(' height="',updated_open,fixed=TRUE)) updated_open <- sub(' height="[^"]+"',paste0(' height="',new_height,'"'),updated_open,perl=TRUE)
  else updated_open <- sub(">$",paste0(' height="',new_height,'">'),updated_open)
  escape <- function(x) as.character(htmltools::htmlEscape(x))
  provenance <- paste0('<title>Primer Atlas candidate amplicon map · ',escape(template_id),'</title>',
    '<desc>',escape(.primer_mixer_export_scope),' ',escape(result$method),'</desc>',
    '<metadata id="primer-atlas-mixer-metadata">',escape(metadata_json),'</metadata>')
  footer <- paste0('<g class="primer-mixer-export-provenance" aria-label="Candidate scope, settings and source citations">',
    '<rect x="',bounds[1L],'" y="',bounds[2L]+original_height,'" width="',width,'" height="',footer_height,'" fill="white"/>',
    paste0('<text x="',bounds[1L]+24,'" y="',bounds[2L]+original_height+24+18*(seq_along(lines)-1L),
      '" font-family="Arial,sans-serif" font-size="',font_size,'" fill="#344459">',escape(lines),'</text>',collapse=""),'</g>')
  body <- substring(svg,nchar(open)+1L)
  if(!endsWith(body,"</svg>")) stop("The map SVG has no closing element.")
  paste0(updated_open,provenance,substring(body,1L,nchar(body)-6L),footer,"</svg>")
}

# A standalone, print-ready table and evidence record. All dynamic content is
# passed as text through htmltools; source links accept HTTP(S) URLs only.
# The optional taxonomic report must describe the exact selected cohort.
primer_mixer_export_html <- function(result,taxonomic_report=NULL) {
  required <- c("pairs","components","templates","products","summary","settings","method","marker_id")
  if(!is.list(result) || !all(required %in% names(result))) stop("Formatted export requires a complete primer mixer result.")
  tags <- htmltools::tags
  field <- .primer_mixer_export_field
  text <- function(value,empty="—") {
    if(!length(value) || is.na(value[1L]) || !nzchar(as.character(value[1L]))) empty else as.character(value[1L])
  }
  link <- function(url,label) {
    url <- text(url,""); label <- text(label)
    if(grepl("^https?://",url,ignore.case=TRUE)) tags$a(href=url,rel="noopener",label) else label
  }
  fmt <- function(value) if(!length(value) || is.na(value)) "—" else if(is.infinite(value)) "Unbounded" else format(value,big.mark=",",scientific=FALSE,trim=TRUE)
  bool <- function(value) length(value)==1L && !is.na(value) && tolower(as.character(value)) %in% c("true","1","yes")
  table <- function(headers,rows,class="") tags$table(class=paste("report-table",class),
    tags$thead(tags$tr(lapply(headers,tags$th))),tags$tbody(rows))
  cells <- function(values) tags$tr(lapply(values,tags$td))
  span <- function(values) {
    values <- values[is.finite(values)]
    if(!length(values)) "—" else if(min(values)==max(values)) fmt(min(values)) else paste0(fmt(min(values)),"–",fmt(max(values)))
  }
  pair_rows <- lapply(seq_len(nrow(result$pairs)),function(i) {
    pair <- result$pairs[i,,drop=FALSE]
    f <- result$components[match(pair$forward_id,result$components$primer_id),,drop=FALSE]
    r <- result$components[match(pair$reverse_id,result$components$primer_id),,drop=FALSE]
    summary <- result$summary[match(pair$pair_id,result$summary$pair_id),,drop=FALSE]
    products <- result$products[result$products$pair_id==pair$pair_id,,drop=FALSE]
    oligo <- function(component) {
      source_sequence <- text(field(component,"source_sequence"),"")
      sequence <- text(field(component,"sequence"),"")
      if(!nzchar(source_sequence)) source_sequence <- sequence
      htmltools::tagList(tags$strong(text(component$primer_name)),tags$div(class="sequence",source_sequence),
        if(!identical(source_sequence,sequence)) tags$div(class="model-note","Search model: ",sequence))
    }
    cells(list(tags$div(tags$strong(pair$pair_label),tags$small(class="record-id",pair$pair_id)),
      oligo(f),oligo(r),paste0(summary$n_exact_products," / ",summary$n_templates),
      paste0(summary$n_mismatch_products," / ",summary$n_templates),summary$n_length_eligible,
      htmltools::tagList(tags$div(span(products$amplicon_bp)," all"),
        tags$small(span(products$amplicon_bp[products$passes_length_filter])," length eligible"))))
  })
  pair_table <- table(c("Candidate combination","Forward oligo · 5′→3′","Reverse oligo · 5′→3′",
    "Exact / input records","Substitution / input records","Length eligible","Measured product bp + primers"),pair_rows,"combination-table")
  panel <- result$templates
  if(!is.null(taxonomic_report)) {
    if(!is.list(taxonomic_report) || !all(c("summary","details","templates","metadata") %in% names(taxonomic_report))) stop("Taxonomic report requires summary, details, templates and metadata.")
    if(!is.data.frame(taxonomic_report$templates) || !"template_id" %in% names(taxonomic_report$templates) ||
        anyDuplicated(taxonomic_report$templates$template_id) || !setequal(taxonomic_report$templates$template_id,panel$template_id)) stop("Taxonomic report templates must match the exact selected input cohort.")
    if(!is.null(taxonomic_report$metadata$marker_id) && !identical(as.character(taxonomic_report$metadata$marker_id),as.character(result$marker_id))) stop("Taxonomic report marker differs from the selected candidate result.")
    if(nrow(taxonomic_report$summary) && (!"pair_id" %in% names(taxonomic_report$summary) ||
        any(!taxonomic_report$summary$pair_id %in% result$pairs$pair_id))) stop("Taxonomic report contains a combination outside the selected candidates.")
    taxa <- taxonomic_report$templates[match(panel$template_id,taxonomic_report$templates$template_id),,drop=FALSE]
    for(identity in intersect(c("sequence","version"),intersect(names(taxa),names(panel))))
      if(!identical(as.character(taxa[[identity]]),as.character(panel[[identity]]))) stop("Taxonomic report ",identity," differs from the selected source cohort.")
    for(name in intersect(c("kingdom","phylum","class","order","family","genus","species","taxid",
        "taxonomy_status","taxonomy_source_url","taxonomy_version","taxonomy_source_sha256"),names(taxa))) panel[[name]] <- taxa[[name]]
  }
  panel_rows <- lapply(seq_len(nrow(panel)),function(i) {
    row <- panel[i,,drop=FALSE]
    status <- c(text(field(row,"source_status"),"Source status not supplied"),
      if(bool(field(row,"partial",FALSE))) "Partial source interval",
      if(bool(field(row,"projected",FALSE))) "Projected interval")
    sequence <- field(row,"sequence",NA_character_)
    sha <- if(length(sequence)==1L && !is.na(sequence) && nzchar(sequence)) digest::digest(sequence,algo="sha256",serialize=FALSE) else "Unavailable"
    ancestry <- intersect(c("class","order","family","genus","species"),names(row))
    ancestry <- vapply(ancestry,function(rank) paste0(rank,": ",text(row[[rank]],"Unresolved")),character(1))
    location <- text(field(row,"source_feature_location"),"")
    if(!nzchar(location) && "source_start" %in% names(row)) location <- paste0(text(row$source_start),"..",text(row$source_end))
    cells(list(htmltools::tagList(link(field(row,"source_url"),row$template_id),tags$div(class="model-note",text(field(row,"organism"))),
        tags$small(class="hash",paste0("SHA-256 ",sha))),
      paste(text(field(row,"marker_id")),if(length(sequence)==1L && !is.na(sequence)) paste0(nchar(sequence)," bp") else "sequence unavailable",sep=" · "),
      htmltools::tagList(tags$div(text(field(row,"version"))),tags$small(text(field(row,"accession_version"))," ",location," ",text(field(row,"source_strand"),""))),
      paste(status,collapse=" · "),htmltools::tagList(tags$div(paste(ancestry,collapse=" · ")),
        tags$small(text(field(row,"taxonomy_status"),"Unresolved taxonomy")," · ",text(field(row,"taxonomy_version"))),
        tags$small(class="hash",text(field(row,"taxonomy_source_sha256"),"No taxonomy source hash")))))
  })
  source_rows <- lapply(seq_len(nrow(result$components)),function(i) {
    row <- result$components[i,,drop=FALSE]
    cells(list(htmltools::tagList(tags$strong(row$primer_name),tags$small(class="record-id",row$primer_id)),
      htmltools::tagList(link(field(row,"source_url"),field(row,"citation")),tags$small(text(field(row,"source_url"),"No publication URL assigned"))),
      htmltools::tagList(tags$div(text(field(row,"target_scope"),"Scope not supplied")),
        tags$small(text(field(row,"evidence_note"),"")),tags$small(text(field(row,"normalization_note"),"")),
        if(text(field(row,"catalog_origin"),"")=="derived_reverse_complement") tags$small(
          "Parent: ",text(field(row,"parent_primer_id"))," · ",text(field(row,"parent_sequence"))))))
  })
  detail_rows <- lapply(seq_len(nrow(result$products)),function(i) {
    row <- result$products[i,,drop=FALSE]
    label <- result$pairs$pair_label[match(row$pair_id,result$pairs$pair_id)]
    cells(list(label,row$template_id,gsub("_"," ",row$status),fmt(row$amplicon_bp),fmt(row$internal_bp),
      paste0(fmt(row$forward_mismatches)," / ",fmt(row$reverse_mismatches)),
      paste0(fmt(row$forward_terminal_5_mismatches)," / ",fmt(row$reverse_terminal_5_mismatches)),
      if(bool(row$passes_length_filter)) "Yes" else "No"))
  })
  taxonomy_section <- NULL
  if(!is.null(taxonomic_report)) {
    summary <- taxonomic_report$summary
    columns <- intersect(c("pair_label","rank","taxon","parent_path","n_panel","n_exact","n_substitution","n_partial","n_ambiguous","n_missing","n_scorable","n_length_eligible"),names(summary))
    tax_rows <- lapply(seq_len(nrow(summary)),function(i) cells(lapply(summary[i,columns,drop=FALSE],function(value) text(value))))
    headers <- c(pair_label="Combination",rank="Rank",taxon="Taxon",parent_path="Parent hierarchy",n_panel="Input records",n_exact="Exact",
      n_substitution="Substitution",n_partial="Partial / short",n_ambiguous="Ambiguous",n_missing="Missing / withheld",n_scorable="Two unique sites",n_length_eligible="Length eligible")
    taxonomy_section <- tags$section(class="report-section",tags$h2("Taxonomic reference fit"),
      tags$p("These hierarchy-scoped denominators contain the selected input records. User-provided taxonomy remains labeled as user metadata; it does not become NCBI-verified taxonomy."),
      if(nrow(summary)) table(unname(headers[columns]),tax_rows,"taxonomy-table") else tags$p("No taxonomic summaries are available for this selection."),
      tags$details(class="audit-details",tags$summary("Taxonomic report metadata and accession-level evidence"),
        tags$pre(jsonlite::toJSON(taxonomic_report,auto_unbox=TRUE,pretty=TRUE,na="null",null="null",dataframe="rows",digits=NA))))
  }
  limits <- paste0("At most ",result$settings$max_mismatches," substitutions per primer; inclusive product range ",fmt(result$settings$min_bp),"–",fmt(result$settings$max_bp)," bp",
    if(!is.null(result$settings$platform_ceiling_bp)) paste0("; project platform ceiling ",fmt(result$settings$platform_ceiling_bp)," bp") else "", ".")
  metadata <- list(schema_version=1L,artifact_type="primer_atlas_candidate_table",marker_id=result$marker_id,
    n_requested_combinations=nrow(result$pairs),n_input_templates=nrow(result$templates),settings=result$settings,
    pairs=result$pairs,components=result$components,template_provenance=panel[,setdiff(names(panel),c("sequence","source_sequence")),drop=FALSE],
    evidence_scope=.primer_mixer_export_scope,search_method=result$method)
  css <- "@page{size:A4 landscape;margin:12mm}*{box-sizing:border-box}body{margin:0;background:#eef2f6;color:#233343;font:14px/1.5 Arial,sans-serif}main{max-width:1600px;margin:28px auto;padding:32px;background:white}h1{font-size:28px;line-height:1.2;margin:0 0 12px}h2{font-size:18px;margin:26px 0 10px}p{margin:8px 0 14px}.eyebrow{font-size:11px;text-transform:uppercase;letter-spacing:.1em;color:#627487}.scope{border-left:4px solid #557a94;background:#f1f6f9;padding:12px 16px}.report-table{width:100%;border-collapse:collapse;table-layout:auto;margin:12px 0 20px;font-size:12px}.report-table th{padding:9px 8px;text-align:left;background:#eaf0f5;color:#294157;border-bottom:2px solid #aebfcd}.report-table td{vertical-align:top;padding:9px 8px;border-bottom:1px solid #dbe3ea;overflow-wrap:anywhere}.report-table tr:nth-child(even) td{background:#f8fafc}.sequence{font:11px/1.55 monospace;overflow-wrap:anywhere;margin-top:3px;min-width:130px}.model-note,small{display:block;font-size:10px;color:#536779;margin-top:4px}.record-id,.hash{font:9px/1.5 monospace;overflow-wrap:anywhere}a{color:#225c84;text-decoration:none}details{border:1px solid #dbe3ea;padding:10px 12px;margin:16px 0}summary{cursor:pointer;font-weight:600}pre{font:10px/1.5 monospace;white-space:pre-wrap;overflow-wrap:anywhere}.footer{font-size:11px;color:#536779;margin-top:24px;border-top:1px solid #ccd7e0;padding-top:10px}@media print{body{background:white;font-size:10pt}main{max-width:none;margin:0;padding:0}h1{font-size:20pt}h2{font-size:13pt;break-after:avoid}.report-table{font-size:8pt}.report-table thead{display:table-header-group}.report-table tr{break-inside:avoid}.sequence{font-size:8pt;min-width:100px}a{color:inherit}.audit-details{display:none}.report-section{break-inside:auto}.scope{background:white}.record-id,.hash{font-size:6.5pt}}"
  html <- tags$html(lang="en",tags$head(tags$meta(charset="utf-8"),tags$meta(name="viewport",content="width=device-width,initial-scale=1"),
    tags$title(paste("Primer Atlas",result$marker_id,"candidate combinations")),tags$style(htmltools::HTML(css))),
    tags$body(tags$main(tags$div(class="eyebrow","Primer Atlas · reproducible candidate table"),
      tags$h1(paste(result$marker_id,"primer combinations")),
      tags$p(nrow(result$pairs)," requested combinations across ",nrow(result$templates)," explicit input templates."),
      tags$div(class="scope",.primer_mixer_export_scope),tags$p(limits),
      tags$section(class="report-section",tags$h2("Combination table"),
        tags$p("Oligo sequences are shown 5′→3′ in their original/source chemistry; a differing normalized search sequence is disclosed. Measured span ranges summarize actual inward products across templates; they are not a shared alignment coordinate or a predicted PCR length."),
        if(nrow(result$pairs)) pair_table else tags$p("No forward × reverse combinations were selected.")),
      taxonomy_section,
      tags$section(class="report-section",tags$h2("Coordinate templates and source provenance"),
        tags$p("Local product coordinates start at one on each supplied/extracted sequence. Parent accession coordinates and strand identify its source orientation. Partial and predicted-model sources remain explicit."),
        if(nrow(panel)) table(c("Template / organism / sequence hash","Marker / raw length","Frozen version / source interval","Source boundary evidence","Taxonomy / metadata source"),panel_rows) else tags$p("No templates were supplied.")),
      tags$section(class="report-section",tags$h2("Primer component citations and scope"),
        if(nrow(result$components)) table(c("Component identity","Source citation","Original scope / candidate notes"),source_rows) else tags$p("No components selected.")),
      tags$section(class="report-section",tags$h2("Methods and interpretation"),tags$p(result$method),
        tags$p("Inclusive product bp = reverse-site end − forward-site start + 1. Internal bp = reverse-site start − forward-site end − 1. Missing or partial sites can reflect an incomplete reference or sequence differences; a missing match is not labeled PCR failure. Length/platform filters change eligible products without changing source sequences or the input-record denominators.")),
      tags$details(class="audit-details",tags$summary("All accession-level template results"),
        table(c("Combination","Template","Evidence state","Product bp + primers","Internal bp","F / R substitutions","F / R last-five 3′ substitutions","Length eligible"),detail_rows)),
      tags$details(class="audit-details",tags$summary("Structured source and method metadata"),
        tags$pre(jsonlite::toJSON(metadata,auto_unbox=TRUE,pretty=TRUE,na="null",null="null",dataframe="rows",digits=NA))),
      tags$p(class="footer","This report records the exact candidate selection, settings and source records supplied to Primer Atlas. Source-specific empirical assay evidence must be cited and interpreted separately."))))
  rendered <- htmltools::renderTags(html)
  root_open <- regmatches(rendered$html,regexpr("^<html\\b[^>]*>",rendered$html,perl=TRUE))
  if(!length(root_open) || !nzchar(root_open)) stop("Could not render the standalone table document.")
  paste0("<!DOCTYPE html>\n",root_open,"\n<head>",rendered$head,"</head>",substring(rendered$html,nchar(root_open)+1L))
}
