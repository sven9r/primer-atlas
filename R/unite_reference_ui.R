# Read-only views of one frozen UNITE panel. Regional alignments are supplied
# by the reference reader; this module never realigns or invents missing rows.
.unite_reference_ui_value <- function(value, fallback) if(is.null(value)) fallback else value
.unite_reference_ui_count <- function(value) format(as.integer(value),big.mark=",",trim=TRUE)
.unite_reference_ui_empty <- function(message) shiny::div(class="unite-empty",message)

.unite_reference_ui_alignment_window <- function(result,start=1L,width=80L) {
  alignment <- result$alignment
  if(is.null(alignment)) alignment <- character()
  if(!is.character(alignment) || (length(alignment) && (is.null(names(alignment)) || anyDuplicated(names(alignment)))))
    stop("Regional alignment lacks unique record identities")
  lengths <- nchar(alignment)
  if(length(lengths) && (anyNA(lengths) || length(unique(lengths))!=1L)) stop("Regional alignment rows have unequal widths")
  columns <- if(length(lengths)) lengths[1L] else 0L
  width <- min(120L,max(1L,as.integer(width)))
  start <- max(1L,min(as.integer(start),max(1L,columns-width+1L)))
  end <- min(columns,start+width-1L)
  window <- if(columns) substring(alignment,start,end) else setNames(character(),character())
  metadata <- result$metadata
  if(is.null(metadata)) metadata <- data.frame(record_id=character())
  if(length(alignment) && (!"record_id" %in% names(metadata) || !all(names(alignment) %in% metadata$record_id)))
    stop("Regional alignment metadata does not match its records")
  missing <- result$missing_metadata
  if(is.null(missing)) missing <- metadata[!metadata$record_id %in% names(alignment),,drop=FALSE]
  list(alignment=window,metadata=metadata[match(names(alignment),metadata$record_id),,drop=FALSE],
       missing_metadata=missing,start=start,end=end,n_columns=columns,n_window=if(columns) end-start+1L else 0L)
}

.unite_reference_ui_alignment_svg <- function(window,region) {
  if(!length(window$alignment)) return(.unite_reference_ui_empty("No regional alignment is available for this group."))
  # Twelve pixels per column keeps actual nucleotides legible. Narrow panes
  # scroll within the figure; its axis remains alignment columns, never bp.
  cell <- 12L; label_width <- 190L; row_height <- 30L; top <- 48L
  figure_width <- label_width+cell*window$n_window+12L
  figure_height <- top+row_height*length(window$alignment)+12L
  svg <- htmltools::tags$svg(xmlns="http://www.w3.org/2000/svg",role="img",
    `aria-label`=paste("UNITE",region,"regional alignment; columns",window$start,"to",window$end),
    viewBox=paste(0,0,figure_width,figure_height),width=figure_width,height=figure_height,
    `data-region`=region,`data-start-column`=window$start,`data-end-column`=window$end,
    htmltools::tags$text(x=0,y=18,fill="#475569",`font-size`=11,"Organism / accession"),
    htmltools::tags$text(x=label_width,y=18,fill="#475569",`font-size`=11,"Alignment columns (gaps included)"))
  ticks <- unique(c(window$start,seq.int(window$start,window$end)[seq.int(window$start,window$end)%%10L==0L],window$end))
  for(column in ticks) svg <- htmltools::tagAppendChild(svg,htmltools::tags$text(
    x=label_width+(column-window$start)*cell+cell/2,y=37,`text-anchor`="middle",fill="#475569",`font-size`=10,column))
  colors <- c(A="#dcfce7",C="#dbeafe",G="#fef3c7",T="#ffe4e6")
  for(i in seq_along(window$alignment)) {
    record <- names(window$alignment)[i]
    metadata <- window$metadata[i,,drop=FALSE]
    organism <- if("organism" %in% names(metadata)) as.character(metadata$organism) else ""
    organism[is.na(organism)] <- ""
    label <- record
    for(column in c("accession","species_hypothesis")) if(column %in% names(metadata) &&
        !is.na(metadata[[column]]) && nzchar(as.character(metadata[[column]]))) label <- as.character(metadata[[column]])
    primary <- if(nzchar(organism)) organism else label
    primary_display <- if(nchar(primary)>28L) paste0(substr(primary,1L,27L),"…") else primary
    secondary <- if("accession" %in% names(metadata) && !is.na(metadata$accession)) metadata$accession else label
    y <- top+(i-1L)*row_height
    row <- htmltools::tags$g(class="unite-alignment-row",`data-record-id`=record,
      htmltools::tags$title(paste(label,organism,"· record",record)),
      htmltools::tags$text(x=0,y=y+11,fill="#172033",`font-size`=11,primary_display),
      htmltools::tags$text(x=0,y=y+24,fill="#64748b",`font-size`=9,secondary))
    bases <- strsplit(window$alignment[i],"",fixed=TRUE)[[1L]]
    for(j in seq_along(bases)) {
      base <- bases[j]
      upper <- toupper(base)
      state <- if(base=="-") "gap" else if(upper %in% names(colors)) "resolved" else "ambiguous"
      fill <- if(state=="gap") "#ffffff" else if(state=="resolved") colors[[upper]] else "#ede9fe"
      ink <- if(state=="gap") "#64748b" else if(state=="ambiguous") "#6d28d9" else "#172033"
      x <- label_width+(j-1L)*cell
      row <- htmltools::tagAppendChild(row,htmltools::tags$g(class="unite-alignment-base",
        `data-record-id`=record,`data-column`=window$start+j-1L,`data-base`=base,`data-state`=state,
        htmltools::tags$title(paste(record,"· alignment column",window$start+j-1L,"·",base,
          if(state=="gap") "alignment gap" else if(state=="ambiguous") "ambiguous / IUPAC base" else "resolved base")),
        htmltools::tags$rect(x=x,y=y,width=cell-1L,height=22,fill=fill),
        htmltools::tags$text(x=x+cell/2,y=y+16,`text-anchor`="middle",fill=ink,`font-size`=12,
          `font-family`="ui-monospace, SFMono-Regular, Consolas, monospace",base)))
    }
    svg <- htmltools::tagAppendChild(svg,row)
  }
  shiny::div(class="unite-alignment-scroll",svg)
}

uniteReferenceUI <- function(id,initial_mode="alignment") {
  initial_mode <- match.arg(initial_mode,c("coverage","alignment","primers"))
  ns <- shiny::NS(id)
  shiny::div(class="unite-reference atlas-workspace",
    shiny::tags$style(shiny::HTML(".unite-reference{display:flex;flex-direction:column;gap:14px;min-width:0}.unite-reference .unite-heading{margin:0;font-size:1.2rem}.unite-reference .unite-context{font-size:.84rem;color:#526176;line-height:1.5}.unite-reference .unite-controls{display:flex;flex-wrap:wrap;gap:12px;align-items:end}.unite-reference .unite-controls>.shiny-input-container{flex:1 1 180px;width:auto!important;min-width:0;margin-bottom:0}.unite-reference .unite-controls>.unite-region-select{flex:1 1 260px}.unite-reference .unite-region-select .shiny-input-container{width:100%!important;margin:0}.unite-reference .unite-pane{min-width:0}.unite-reference .unite-empty{padding:14px;border:1px dashed #cbd5e1;border-radius:8px;color:#526176}.unite-reference .unite-breadth-scroll{overflow-x:auto}.unite-reference .unite-breadth-table{width:100%;border-collapse:collapse;font-size:.84rem}.unite-reference .unite-breadth-table th,.unite-reference .unite-breadth-table td{padding:10px 8px;border-bottom:1px solid #e2e8f0;vertical-align:top}.unite-reference .unite-breadth-table th{text-align:left;font-size:.76rem;color:#526176;font-weight:650}.unite-reference .unite-breadth-table .number{text-align:right;font-variant-numeric:tabular-nums}.unite-reference .unite-group-button{padding:0;border:0;background:transparent;text-align:left;color:#4338ca;font:inherit;font-weight:700;cursor:pointer}.unite-reference .unite-group-button:hover{text-decoration:underline}.unite-reference .unite-group-button:focus-visible{outline:3px solid #4338ca;outline-offset:3px}.unite-reference .unite-group-meta{font-size:.73rem;color:#64748b;line-height:1.4}.unite-reference .unite-shortfall{font-size:.75rem;color:#92400e}.unite-reference .unite-count-track{height:5px;background:#e2e8f0;border-radius:3px;margin-top:5px;overflow:hidden}.unite-reference .unite-count-fill{display:block;height:100%;background:#64748b;min-width:1px}.unite-reference .unite-alignment-scroll{overflow:auto;max-height:430px;border:1px solid #dbe3ee;border-radius:8px;padding:10px;overscroll-behavior-x:contain}.unite-reference .unite-alignment-scroll svg{display:block;max-width:none;flex:none}.unite-reference .unite-alignment-legend{display:flex;flex-wrap:wrap;gap:6px 16px;font-size:.76rem;color:#526176;margin:8px 0}.unite-reference .unite-alignment-key{display:inline-block;width:12px;height:12px;border:1px solid #cbd5e1;margin-right:5px;vertical-align:-1px}.unite-reference .unite-missing{font-size:.8rem;color:#526176;margin:8px 0}.unite-reference .unite-missing summary{cursor:pointer}.unite-reference .unite-site-tally{margin:10px 0;font-size:.85rem}.unite-reference .unite-pagination{display:flex;justify-content:space-between;align-items:center;gap:10px;margin-top:10px}.unite-reference .dataTables_wrapper{max-width:100%;overflow-x:auto}.unite-reference .unite-source summary{cursor:pointer;font-size:.85rem;color:#526176}.unite-reference .unite-source p{font-size:.8rem;line-height:1.5}.unite-reference .unite-scroll-cue{font-size:.74rem;color:#64748b;margin:7px 0}@media(max-width:420px){.unite-reference .unite-controls>.shiny-input-container{flex-basis:100%}.unite-reference .unite-breadth-table th,.unite-reference .unite-breadth-table td{padding:8px 5px}.unite-reference .unite-breadth-table{font-size:.78rem}}")),
    shiny::tags$script(shiny::HTML("if(!window.atlasUniteGroupClick){window.atlasUniteGroupClick=true;document.addEventListener('click',function(e){var b=e.target.closest('[data-unite-group]');if(b&&!b.disabled&&window.Shiny){Shiny.setInputValue(b.dataset.input,{group:b.dataset.uniteGroup,rank:b.dataset.rank,nonce:Date.now()},{priority:'event'});}});}")),
    shiny::tags$span(style="display:none",shiny::textOutput(ns("available"),inline=TRUE)),
    shiny::tags$style(shiny::HTML(".unite-reference .unite-body{display:flex;flex-direction:column;gap:14px;min-width:0}.unite-reference .unite-group-button:disabled{color:#526176;cursor:default;text-decoration:none}.unite-reference .unite-window-controls{display:grid;grid-template-columns:130px minmax(0,1fr);gap:12px;align-items:end}.unite-reference .unite-window-controls .shiny-input-container{width:100%!important;min-width:0;margin-bottom:0}@media(max-width:420px){.unite-reference .unite-window-controls{grid-template-columns:1fr}}")),
    shiny::uiOutput(ns("availability")),
    shiny::conditionalPanel(condition="output.available == 'true'",ns=ns,class="unite-body",
    shiny::h3(class="unite-heading","UNITE reference evidence"),
    shiny::uiOutput(ns("source_context")),
    shiny::div(class="atlas-segmented",shiny::radioButtons(ns("mode"),NULL,
      c("Reference breadth"="coverage","Regional alignment"="alignment","Primer sites"="primers"),initial_mode,inline=TRUE)),
    shiny::conditionalPanel(condition="input.mode == 'coverage'",ns=ns,
      shiny::div(class="unite-pane",
        shiny::div(class="unite-controls",shiny::selectInput(ns("breadth_rank"),"Group references by",c("Classes"="class"),"class"),
          shiny::textInput(ns("breadth_search"),"Find a group",placeholder="Taxon or parent lineage")),
        shiny::uiOutput(ns("breadth_instruction")),
        shiny::uiOutput(ns("breadth")),
        shiny::div(class="unite-pagination",shiny::actionButton(ns("previous"),"Previous",class="btn-outline-secondary btn-sm"),
          shiny::uiOutput(ns("breadth_page"),inline=TRUE),shiny::actionButton(ns("next"),"Next",class="btn-outline-secondary btn-sm")),
        shiny::downloadButton(ns("download_breadth"),"Group counts (CSV)",class="btn-outline-secondary btn-sm"))),
    shiny::conditionalPanel(condition="input.mode == 'alignment'",ns=ns,
      shiny::div(class="unite-pane",
        shiny::div(class="unite-controls",shiny::selectInput(ns("alignment_rank"),"Taxonomic rank",c("Orders"="order","Classes"="class"),"order"),
          shiny::selectInput(ns("alignment_group"),"Reference group",choices=NULL),
          shiny::div(class="atlas-segmented unite-region-select",shiny::radioButtons(ns("region"),"Region",
            c("ITS1"="ITS1","5.8S"="5.8S","ITS2"="ITS2"),"ITS1",inline=TRUE))),
        shiny::uiOutput(ns("alignment_context")),shiny::uiOutput(ns("window_control")),
        shiny::div(class="unite-alignment-legend",
          shiny::span(shiny::span(class="unite-alignment-key",style="background:#dcfce7"),"A / C / G / T · base colors"),
          shiny::span(shiny::span(class="unite-alignment-key",style="background:#ffffff"),"− · alignment gap"),
          shiny::span(shiny::span(class="unite-alignment-key",style="background:#ede9fe"),"Ambiguous / IUPAC base"),
          shiny::span("Missing regional evidence is listed separately")),
        shiny::uiOutput(ns("alignment_plot")),
        shiny::uiOutput(ns("window_note")),
        shiny::uiOutput(ns("missing_regions")),
        shiny::downloadButton(ns("download_alignment"),"Displayed regional alignment (FASTA)",class="btn-outline-secondary btn-sm"))),
    shiny::conditionalPanel(condition="input.mode == 'primers'",ns=ns,
      shiny::div(class="unite-pane",
        shiny::div(class="unite-controls",shiny::selectizeInput(ns("pair"),"Primer combination",choices=NULL,
          options=list(placeholder="Search primer combinations")),
          shiny::selectInput(ns("sites_rank"),"Taxonomic rank",c("Classes"="class"),"class")),
        shiny::uiOutput(ns("site_tally")),
        shiny::p(class="unite-context","Exact sites describe these retained sequences, not PCR success. Click a taxon to inspect product, single-site, missing and ambiguous evidence."),
        DT::DTOutput(ns("site_summary")),shiny::uiOutput(ns("site_records_context")),DT::DTOutput(ns("site_records")),
        shiny::downloadButton(ns("download_sites"),"Complete primer-site evidence (CSV)",class="btn-outline-secondary btn-sm"))),
    shiny::tags$details(class="unite-source",shiny::tags$summary("Source, selection and annotation evidence"),shiny::uiOutput(ns("source_evidence")))))
}

uniteReferenceServer <- function(id,dataset=unite_reference_read,pair_provider=NULL,initial_mode="alignment") {
  initial_mode <- match.arg(initial_mode,c("coverage","alignment","primers"))
  shiny::moduleServer(id,function(input,output,session) {
    supports_versions <- is.function(dataset) && "version_id" %in% names(formals(dataset))
    raw_data <- shiny::reactive(tryCatch(if(supports_versions) dataset(version_id=input$version)
      else if(is.function(dataset)) dataset() else dataset,
      error=function(e) structure(list(message=conditionMessage(e)),class="unite_unavailable")))
    available <- shiny::reactive({
      value <- raw_data()
      is.list(value) && is.data.frame(value$metadata) && nrow(value$metadata)>0L && is.data.frame(value$coverage)
    })
    mode <- shiny::reactive(if(length(input$mode)==1L && input$mode %in% c("coverage","alignment","primers")) input$mode else initial_mode)
    output$available <- shiny::renderText(if(available()) "true" else "false")
    shiny::outputOptions(output,"available",suspendWhenHidden=FALSE)
    output$availability <- shiny::renderUI({
      if(available()) return(NULL)
      value <- raw_data()
      .unite_reference_ui_empty(if(inherits(value,"unite_unavailable"))
        "The UNITE reference panel could not be verified. No alignment or primer-site results are shown."
      else "The UNITE reference panel is not available yet. Regional alignments and primer-site evidence appear after a verified frozen snapshot is installed.")
    })
    data <- shiny::reactive({
      value <- raw_data()
      shiny::validate(shiny::need(available(),
        "A complete frozen UNITE panel is not installed."))
      value
    })
    ranks <- shiny::reactive({
      available <- intersect(c("class","order","phylum","family"),unique(as.character(data()$coverage$rank)))
      shiny::validate(shiny::need(length(available)>0L,"This release has no taxonomic group counts."))
      available
    })
    rank_value <- function(value) if(length(value)==1L && value %in% ranks()) value else ranks()[1L]
    breadth_rank <- shiny::reactive(rank_value(input$breadth_rank))
    selected_group <- shiny::reactiveVal(NULL)
    alignment_ranks <- shiny::reactive(intersect(c("order","class"),ranks()))
    alignment_rank <- shiny::reactive({
      clicked <- selected_group()
      if(!is.null(clicked) && clicked$rank %in% alignment_ranks()) return(clicked$rank)
      if(length(input$alignment_rank)==1L && input$alignment_rank %in% alignment_ranks()) return(input$alignment_rank)
      alignment_ranks()[1L]
    })
    sites_rank <- shiny::reactive(rank_value(input$sites_rank))
    groups <- shiny::reactive(unite_reference_groups(data(),rank=breadth_rank(),search=.unite_reference_ui_value(input$breadth_search,"")))
    alignment_groups <- shiny::reactive(unite_reference_groups(data(),rank=alignment_rank(),search=""))
    page <- shiny::reactiveVal(1L)
    pages <- shiny::reactive(max(1L,ceiling(nrow(groups())/8L)))
    current_page <- shiny::reactive(min(max(1L,page()),pages()))
    shiny::observeEvent(list(breadth_rank(),input$breadth_search),page(1L),ignoreInit=TRUE)
    shiny::observeEvent(input$previous,page(max(1L,current_page()-1L)))
    shiny::observeEvent(input[["next"]],page(min(pages(),current_page()+1L)))
    shiny::observeEvent(ranks(),{
      labels <- c(class="Classes",order="Orders",phylum="Phyla",family="Families")
      for(name in c("breadth_rank","sites_rank")) shiny::updateSelectInput(session,name,
        choices=setNames(ranks(),labels[ranks()]),selected=shiny::isolate(rank_value(input[[name]])))
      shiny::updateSelectInput(session,"alignment_rank",choices=setNames(alignment_ranks(),labels[alignment_ranks()]),selected=alignment_rank())
    })
    shiny::observeEvent(input$group_click,{
      value <- input$group_click
      if(!is.list(value) || length(value$group)!=1L || !value$group %in% groups()$group_id) return()
      if(!breadth_rank() %in% c("class","order")) return()
      selected_group(list(id=value$group,rank=breadth_rank()))
      shiny::updateSelectInput(session,"alignment_rank",selected=breadth_rank())
      shiny::updateSelectInput(session,"alignment_group",selected=value$group)
      shiny::updateRadioButtons(session,"mode",selected="alignment")
    })
    default_group <- function(rows,rank,region) {
      if(!nrow(rows)) return(NULL)
      index <- data()$alignment_index
      if("rank" %in% names(index)) index <- index[index$rank==rank,,drop=FALSE]
      regional <- rows[rows$group_id %in% index$group_id[index$region==region],,drop=FALSE]
      if(nrow(regional)) return(regional$group_id[which.max(regional$n_selected)])
      built <- rows[rows$group_id %in% index$group_id,,drop=FALSE]
      if(nrow(built)) return(built$group_id[which.max(built$n_selected)])
      rows$group_id[which.max(rows$n_selected)]
    }
    alignment_group <- shiny::reactive({
      rows <- alignment_groups()
      if(!nrow(rows)) return(NULL)
      clicked <- selected_group()
      if(!is.null(clicked) && identical(clicked$rank,alignment_rank()) && clicked$id %in% rows$group_id) return(clicked$id)
      if(length(input$alignment_group)==1L && input$alignment_group %in% rows$group_id) return(input$alignment_group)
      default_group(rows,alignment_rank(),region())
    })
    shiny::observeEvent(list(input$alignment_group,input$alignment_rank),{
      clicked <- selected_group()
      if(!is.null(clicked) && identical(input$alignment_group,clicked$id) && identical(input$alignment_rank,clicked$rank)) selected_group(NULL)
    },ignoreInit=TRUE)
    shiny::observeEvent(list(alignment_groups(),region()),{
      rows <- alignment_groups()
      shiny::updateSelectInput(session,"alignment_group",choices=setNames(rows$group_id,
        paste0(rows$taxon," · ",rows$n_selected," selected")),selected=alignment_group())
    })
    region <- shiny::reactive(if(length(input$region)==1L && input$region %in% c("ITS1","5.8S","ITS2")) input$region else "ITS1")
    alignment <- shiny::reactive({
      if(is.null(alignment_group())) return(list(alignment=character(),metadata=data.frame(record_id=character()),
        n_total=0L,n_aligned=0L,n_columns=0L,reason="No reference groups are available.",scope=""))
      unite_reference_alignment(data(),group_id=alignment_group(),region=region(),limit=12L)
    })
    window_width <- shiny::reactive({
      value <- suppressWarnings(as.integer(input$window_width))
      if(length(value)==1L && !is.na(value) && value %in% c(60L,80L,120L)) value else 80L
    })
    window <- shiny::reactive(.unite_reference_ui_alignment_window(alignment(),
      start=.unite_reference_ui_value(input$window_start,1L),width=window_width()))
    pair_catalog <- shiny::reactive(data()$pairs)
    requested_pair <- shiny::reactiveVal(NULL)
    default_pair <- shiny::reactive({
      ids <- pair_catalog()$pair_id
      if(!length(ids)) return(NULL)
      preferred <- pair_catalog()$pair_id[pair_catalog()$pair_label=="ITS1 + ITS4"]
      if(length(preferred)) preferred[1L] else ids[1L]
    })
    pair <- shiny::reactive({
      ids <- pair_catalog()$pair_id
      if(!length(ids)) return(NULL)
      requested <- requested_pair()
      if(length(requested)==1L && requested %in% ids) return(requested)
      if(length(input$pair)==1L && input$pair %in% ids) return(input$pair)
      default_pair()
    })
    shiny::observeEvent(pair_catalog(),shiny::updateSelectizeInput(session,"pair",
      choices=setNames(pair_catalog()$pair_id,pair_catalog()$pair_label),selected=pair(),server=TRUE))
    if(!is.null(pair_provider)) shiny::observeEvent(pair_provider(),{
      provided <- pair_provider()
      if(is.list(provided)) provided <- .unite_reference_ui_value(provided$pair_id,provided$pair)
      if(length(provided)==1L && provided %in% pair_catalog()$pair_id) {
        requested_pair(provided)
        shiny::updateSelectizeInput(session,"pair",selected=provided,server=TRUE)
      }
    },ignoreNULL=TRUE)
    shiny::observeEvent(input$pair,{
      if(identical(input$pair,requested_pair())) requested_pair(NULL)
    },ignoreNULL=TRUE)
    support <- shiny::reactive({
      if(is.null(pair())) return(list(summary=data.frame(),details=data.frame()))
      unite_reference_pair_fit(data(),pair_id=pair(),rank=sites_rank())
    })
    chosen_site_group <- shiny::reactiveVal(NULL)
    previous_identity <- shiny::reactiveVal(NULL)
    shiny::observeEvent(if(available()) data()$manifest$id else NULL,{
      identity <- data()$manifest$id
      previous <- previous_identity()
      previous_identity(identity)
      if(is.null(previous) || identical(identity,previous)) return()
      rank <- alignment_ranks()[1L]
      rows <- unite_reference_groups(data(),rank=rank,search="")
      group <- default_group(rows,rank,region())
      selected_group(if(is.null(group)) NULL else list(id=group,rank=rank))
      requested_pair(default_pair())
      chosen_site_group(NULL)
      shiny::updateSelectInput(session,"alignment_rank",selected=rank)
      shiny::updateSelectInput(session,"alignment_group",selected=group)
      shiny::updateSelectizeInput(session,"pair",selected=default_pair(),server=TRUE)
    },ignoreNULL=TRUE,priority=10L)
    shiny::observeEvent(list(pair(),sites_rank(),data()$manifest$id),chosen_site_group(NULL),priority=10L)
    shiny::observeEvent(input$site_summary_rows_selected,{
      row <- input$site_summary_rows_selected
      rows <- support()$summary
      if(length(row)!=1L || row<1L || row>nrow(rows)) chosen_site_group(NULL) else chosen_site_group(rows$group_id[row])
    })
    site_group <- shiny::reactive({
      group <- chosen_site_group()
      if(is.null(group) || !group %in% support()$summary$group_id) return(NULL)
      group
    })
    site_details <- shiny::reactive({
      rows <- support()$details
      group <- site_group()
      if(is.null(group)) return(rows[0,,drop=FALSE])
      rows[rows$group_id==group,,drop=FALSE]
    })
    output$source_context <- shiny::renderUI({
      manifest <- data()$manifest
      shiny::div(class="unite-context",strong("UNITE · "),
        .unite_reference_ui_count(manifest$n_source)," source records · ",
        .unite_reference_ui_count(nrow(data()$metadata))," selected references · ",
        "separate within-group regional alignments")
    })
    output$breadth_page <- shiny::renderUI(shiny::span(class="unite-context",paste("Page",current_page(),"of",pages())))
    output$breadth_instruction <- shiny::renderUI(shiny::p(class="unite-context",
      "Source counts describe the release. Selected counts describe this retained panel. ",
      if(breadth_rank() %in% c("class","order")) paste0("Click ",if(breadth_rank()=="order") "an " else "a ",breadth_rank()," to inspect its own regional alignment.")
      else "Regional alignments are available within class and order groups; this rank shows reference breadth."))
    output$breadth <- shiny::renderUI({
      rows <- groups()
      if(!nrow(rows)) return(.unite_reference_ui_empty("No taxonomic groups match this search."))
      rows <- rows[seq.int((current_page()-1L)*8L+1L,min(current_page()*8L,nrow(rows))),,drop=FALSE]
      shiny::div(class="unite-breadth-scroll",shiny::tags$table(class="unite-breadth-table",
        shiny::tags$thead(shiny::tags$tr(shiny::tags$th("Taxon"),shiny::tags$th(class="number","Source"),
          shiny::tags$th(class="number","Unique eligible"),shiny::tags$th(class="number","Selected"))),
        shiny::tags$tbody(lapply(seq_len(nrow(rows)),function(i) {
          row <- rows[i,,drop=FALSE]
          fraction <- if(row$n_source>0L) row$n_selected/row$n_source else 0
          shiny::tags$tr(`data-group-id`=row$group_id,
            shiny::tags$td(shiny::tags$button(type="button",class="unite-group-button",
              disabled=if(!breadth_rank() %in% c("class","order")) "disabled" else NULL,
              `data-input`=session$ns("group_click"),`data-unite-group`=row$group_id,`data-rank`=breadth_rank(),row$taxon),
              shiny::div(class="unite-group-meta",row$parent_path),
              if(is.finite(row$shortfall_to_100) && row$shortfall_to_100>0L)
                shiny::div(class="unite-shortfall",paste(row$shortfall_to_100,"short of 100 selected")),
              if(!isTRUE(row$taxon_resolved)) shiny::div(class="unite-group-meta",
                if("taxon_status" %in% names(row) && identical(as.character(row$taxon_status),"incertae_sedis")) "Incertae sedis"
                else "Taxonomy unresolved")),
            shiny::tags$td(class="number",.unite_reference_ui_count(row$n_source)),
            shiny::tags$td(class="number",.unite_reference_ui_count(row$n_eligible_unique)),
            shiny::tags$td(class="number",.unite_reference_ui_count(row$n_selected),
              shiny::div(class="unite-count-track",role="img",`aria-label`=paste(row$n_selected,"selected of",row$n_source,"source records"),
                if(row$n_selected>0L) shiny::span(class="unite-count-fill",style=paste0("width:",round(100*fraction,3),"%;")))))
        }))))
    })
    output$alignment_context <- shiny::renderUI({
      result <- alignment()
      shiny::div(class="unite-context",strong(region()," · "),length(result$alignment)," displayed of ",result$n_aligned,
        " regional alignment exemplars · ",.unite_reference_ui_count(result$n_columns)," alignment columns",
        if(!is.null(result$n_regions)) shiny::p(result$n_regions," selected references have complete region boundaries of ",result$n_total,
          " in this group",if(!is.null(result$n_missing_boundaries)) paste0(" · ",result$n_missing_boundaries," missing / incomplete boundaries"))
        else shiny::p(result$n_total," selected references in this group; regional boundary counts are unavailable in this alignment result."),
        if(nzchar(.unite_reference_ui_value(result$scope,""))) shiny::tags$details(class="unite-missing",
          shiny::tags$summary("Regional alignment scope"),shiny::p(result$scope)),
        if(nzchar(.unite_reference_ui_value(result$reason,""))) shiny::p(result$reason))
    })
    output$window_control <- shiny::renderUI({
      result <- alignment()
      columns <- if(length(result$alignment)) nchar(result$alignment[1L]) else 0L
      if(!columns) return(NULL)
      maximum <- max(1L,columns-window_width()+1L)
      shiny::div(class="unite-window-controls",
        shiny::selectInput(session$ns("window_width"),"Window columns",c("60"=60L,"80"=80L,"120"=120L),selected=window_width()),
        if(maximum==1L) shiny::p(class="unite-context",paste("All",columns,"alignment columns are shown."))
        else shiny::sliderInput(session$ns("window_start"),paste("First column ·",window_width(),"column window"),min=1L,max=maximum,
          value=shiny::isolate(max(1L,min(.unite_reference_ui_value(input$window_start,1L),maximum))),step=1L,width="100%"))
    })
    output$window_note <- shiny::renderUI({
      current <- window()
      if(!current$n_window) return(NULL)
      shiny::p(class="unite-scroll-cue",paste0("Columns ",current$start,"–",current$end," · ",current$n_window,
        " columns shown. Scroll sideways when needed. Alignment columns include gaps; they are not product bp."))
    })
    output$alignment_plot <- shiny::renderUI(.unite_reference_ui_alignment_svg(window(),region()))
    output$missing_regions <- shiny::renderUI({
      missing <- window()$missing_metadata
      if(!nrow(missing)) return(NULL)
      total <- alignment()$n_missing_boundaries
      shiny::tags$details(class="unite-missing",shiny::tags$summary(if(!is.null(total))
        paste(total,"references have missing / incomplete region boundaries") else paste(nrow(missing),"example references have no displayed regional alignment")),
        shiny::p("Missing regional evidence is not an all-gap sequence and does not establish that the locus is absent."),
        if(!is.null(total) && total>nrow(missing)) shiny::p("Showing",nrow(missing),"examples of",total,"missing-boundary records."),
        shiny::tags$ul(lapply(seq_len(nrow(missing)),function(i) shiny::tags$li(missing$record_id[i],
          if("status" %in% names(missing)) paste0(" · ",gsub("_"," ",missing$status[i],fixed=TRUE))))))
    })
    output$site_tally <- shiny::renderUI({
      rows <- support()$summary
      if(!nrow(rows)) return(.unite_reference_ui_empty("No primer-site evidence is available for this combination."))
      selected_pair <- pair_catalog()[pair_catalog()$pair_id==pair(),,drop=FALSE]
      status <- if("evidence_status" %in% names(selected_pair) && nrow(selected_pair)) selected_pair$evidence_status[1L] else "unknown"
      documented <- identical(status,"documented_pair")
      provenance <- if(documented) "Documented pair" else "Candidate combination"
      detail <- if(documented) "Published as a primer pair." else switch(status,
        candidate_published_components="Components have published sources; this combination needs independent validation.",
        candidate_unpublished_component="Includes an unpublished component; this combination needs independent validation.",
        "Pair provenance is incomplete; this combination needs independent validation.")
      shiny::tagList(shiny::p(class="unite-context",shiny::strong(provenance)," · ",detail),
        shiny::p(class="unite-site-tally",strong(sum(rows$n_exact_products)," of ",sum(rows$n_selected)," selected references have an exact product"),
          " · ",sum(rows$n_single_sites)," one-site records · ",sum(rows$n_no_exact_sites)," no exact sites · ",sum(rows$n_unavailable)," unavailable"))
    })
    output$site_summary <- DT::renderDT({
      rows <- support()$summary
      if(!nrow(rows)) return(NULL)
      display <- data.frame(Taxon=rows$taxon,`Exact products`=paste(rows$n_exact_products,"/",rows$n_selected),
        `One exact site`=rows$n_single_sites,`No exact sites`=rows$n_no_exact_sites,
        `Ambiguous / unresolved`=rows$n_unresolved,Unavailable=rows$n_unavailable,
        Withheld=.unite_reference_ui_value(rows$n_withheld,rep(0L,nrow(rows))),check.names=FALSE)
      DT::datatable(display,rownames=FALSE,selection="single",escape=TRUE,
        options=list(pageLength=8L,lengthChange=FALSE,scrollX=TRUE,order=list()),class="compact stripe nowrap")
    })
    output$site_records_context <- shiny::renderUI({
      group <- site_group()
      if(is.null(group)) return(NULL)
      row <- support()$summary[support()$summary$group_id==group,,drop=FALSE]
      shiny::h4(paste(row$taxon,"·",nrow(site_details()),"reference records"))
    })
    output$site_records <- DT::renderDT({
      rows <- site_details()
      if(!nrow(rows)) return(NULL)
      labels <- c(exact_product="Exact product",partial_sites="One exact site",no_exact_sites="No exact sites found",
        unresolved_sites="Ambiguous / unresolved",evidence_unavailable="Evidence unavailable",product_withheld="Product withheld")
      states <- ifelse(rows$fit_state %in% names(labels),unname(labels[rows$fit_state]),gsub("_"," ",rows$fit_state,fixed=TRUE))
      display <- data.frame(Accession=rows$accession,Organism=rows$organism,
        Evidence=states,`Encoded product (bp)`=rows$amplicon_bp,Reason=rows$reason,check.names=FALSE)
      if("n_ambiguous_product_bases" %in% names(rows)) display$`Ambiguous product bases` <- rows$n_ambiguous_product_bases
      if("species_hypothesis" %in% names(rows)) display <- data.frame(`SH hypothesis`=rows$species_hypothesis,display,check.names=FALSE)
      DT::datatable(display,rownames=FALSE,selection="none",escape=TRUE,
        options=list(pageLength=5L,lengthChange=FALSE,scrollX=TRUE,order=list()),class="compact stripe nowrap")
    })
    output$source_evidence <- shiny::renderUI({
      manifest <- data()$manifest
      doi <- manifest$source_doi
      versions <- data()$versions
      shiny::tagList(
        if(is.data.frame(versions) && nrow(versions)>1L && supports_versions) shiny::selectInput(session$ns("version"),
          "Frozen reference version",choices=setNames(versions$id,versions$label),selected=manifest$id)
        else if(is.data.frame(versions) && nrow(versions)==1L) shiny::p(versions$label[1L]),
        shiny::p(manifest$source_version),
        shiny::p("Reference identity: ",shiny::tags$code(manifest$id)," · ",.unite_reference_ui_value(manifest$created_at,manifest$built_at_utc)),
        shiny::p(.unite_reference_ui_value(manifest$source_taxonomy,manifest$taxonomy_source)," · ",
          .unite_reference_ui_value(manifest$region_annotation,"ITSx 1.1.3 computational region boundaries")),
        if(length(doi)==1L && !is.na(doi) && nzchar(doi)) shiny::p(shiny::tags$a(
          href=if(grepl("^https://",doi)) doi else paste0("https://doi.org/",doi),target="_blank",rel="noopener","UNITE source DOI ↗")),
        shiny::p("Source attribution: Abarenkov and colleagues (2025), UNITE general FASTA release for Fungi. ",
          shiny::tags$a(href="https://creativecommons.org/licenses/by/4.0/",target="_blank",rel="noopener","CC BY 4.0")),
        shiny::p("Selection counts describe this retained reference panel, not environmental abundance or universal fungal coverage. Regional alignments are built separately within a taxonomic group from profile-derived ITS1, 5.8S or ITS2 intervals. No missing record is replaced with a fabricated all-gap row."))
    })
    output$download_breadth <- shiny::downloadHandler(filename=function() paste0(data()$manifest$id,"-",breadth_rank(),"-breadth.csv"),
      content=function(file) write.csv(groups(),file,row.names=FALSE,na=""))
    output$download_alignment <- shiny::downloadHandler(filename=function() paste0(data()$manifest$id,"-",region(),"-regional-alignment.fasta"),
      content=function(file) {
        rows <- alignment()$alignment
        shiny::req(length(rows)>0L)
        writeLines(unlist(lapply(seq_along(rows),function(i) c(paste0(">",names(rows)[i]),rows[i])),use.names=FALSE),file)
      })
    output$download_sites <- shiny::downloadHandler(filename=function() paste0(data()$manifest$id,"-",pair(),"-primer-sites.csv"),
      content=function(file) write.csv(support()$details,file,row.names=FALSE,na=""))
    list(data=data,available=available,mode=mode,groups=groups,alignment=alignment,window=window,pair=pair,support=support,
      alignment_group=alignment_group,site_details=site_details)
  })
}
