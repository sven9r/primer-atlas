# A self-contained Shiny workspace; every view uses one immutable snapshot.
fungal_evidence_label <- function(value) {
  labels <- c(documented_pair = "Documented pair", candidate_published_components = "Candidate · cited components",
              candidate_unpublished_component = "Candidate · unpublished component",
              candidate_incomplete_provenance = "Candidate · incomplete provenance")
  unname(ifelse(value %in% names(labels), labels[value], value))
}

fungal_starter_examples <- function(rows, maximum = 20L) {
  if (!nrow(rows)) return(rows)
  priority <- c("ITS1 + ITS2", "ITS1F + ITS2", "ITS5 + ITS2", "ITS3 + ITS4", "fITS7 + ITS4",
                "gITS7 + ITS4", "ITS86F + ITS4", "ITS1 + ITS4", "ITS1F + ITS4", "ITS5 + ITS4",
                "ITS1 + LR21", "ITS1F + LR21", "LR0R + LR3", "ITS1 + LR3", "ITS1 + LR5",
                "NS7a + LR3", "SSU515Fngs + LR3", "NS31 + LB-z", "FF700 + LR7", "NS1 + LB-y")
  rank <- match(rows$pair_label, priority)
  selected <- which(!is.na(rank))
  selected <- selected[order(rank[selected])]
  remaining <- setdiff(seq_len(nrow(rows)), selected)
  # Fill missing named examples with represented length tiers, avoiding a
  # shortlist consisting entirely of nearby oligo variants.
  tiers <- c("Illumina", "PacBio", "Nanopore", "Unresolved", "Beyond configured limits")
  groups <- lapply(tiers, function(tier) {
    ids <- remaining[rows$length_tier[remaining] == tier]
    ids[order(-rows$n_scorable[ids], rows$pair_label[ids])]
  })
  for (i in seq_len(max(c(0L, lengths(groups))))) {
    selected <- c(selected, unlist(lapply(groups, function(ids) if (length(ids) >= i) ids[i] else integer()), use.names = FALSE))
  }
  rows[head(unique(selected), maximum), , drop = FALSE]
}

fungalExplorerUI <- function(id) {
  ns <- NS(id)
  versions <- fungal_read_index()
  choices <- setNames(versions$id, paste0(versions$label, " · ", substr(versions$id, 4, 11)))
  div(class = "fungal-workspace",
      tags$script(HTML("document.addEventListener('click', function(e) { var b=e.target.closest('.fungal-pair-choice'); if(b && window.Shiny) Shiny.setInputValue(b.dataset.input, {pair:b.dataset.pair, nonce:Date.now()}, {priority:'event'}); });")),
      div(class = "fungal-version-bar",
          selectInput(ns("version"), "Reference version", choices, selected = tail(versions$id, 1)),
          uiOutput(ns("snapshot_meta")),
          downloadButton(ns("download_manifest"), "Snapshot manifest", class = "btn-outline-secondary btn-sm")),
      div(class = "fungal-controls",
          radioButtons(ns("platform"), "Sequencing platform", c("All platforms"="ALL", "Illumina"="ILLUMINA", "PacBio"="PACBIO", "Nanopore"="NANOPORE"), "ALL", inline=TRUE),
          div(class = "fungal-platform-legend",
              span(class="fungal-platform-badge illumina", "Illumina"),
              span(class="fungal-platform-badge pacbio", "PacBio"),
              span(class="fungal-platform-badge nanopore", "Nanopore"),
              span(class="fungal-platform-badge unresolved", "Unresolved")),
          tags$details(tags$summary("Adjust project length ceilings"),
              layout_columns(col_widths=c(4,4,4),
                  sliderInput(ns("illumina_max"), "Illumina ceiling (bp)", min=452,max=500,value=500,step=1),
                  numericInput(ns("pacbio_max"), "PacBio ceiling (bp)",2000,min=500,max=20000,step=100),
                  numericInput(ns("nanopore_max"), "Nanopore ceiling (bp)",3500,min=2000,max=50000,step=100))),
          layout_columns(col_widths=c(4,4,4),
              selectInput(ns("phylum"),"Phylum / clade",c("All represented phyla"="ALL")),
              selectInput(ns("class"),"Class",c("All represented classes"="ALL")),
              selectInput(ns("locus"),"Primer loci",c("All rDNA loci"="ALL"))),
          layout_columns(col_widths=c(4,4,4),
              selectizeInput(ns("forward"),"Forward primers",choices=NULL,multiple=TRUE,options=list(placeholder="All supported forward primers")),
              selectizeInput(ns("reverse"),"Reverse primers",choices=NULL,multiple=TRUE,options=list(placeholder="All supported reverse primers")),
              textInput(ns("search"),"Search combinations",placeholder="Primer name, target, or citation")),
          div(class="d-flex flex-wrap gap-3",
              checkboxInput(ns("examples_only"),"Show 20 starter examples",TRUE),
              checkboxInput(ns("show_unresolved"),"Include unresolved documented pairs",TRUE),
              actionButton(ns("reset"),"Reset filters",class="btn-outline-secondary btn-sm"))),
      uiOutput(ns("discovery_status")),
      navset_card_tab(
          nav_panel("Combination list", DTOutput(ns("pairs"))),
          nav_panel("Pairing matrix", uiOutput(ns("matrix"))),
          nav_panel("Version changes", uiOutput(ns("version_changes"))),
          nav_panel("Frozen primer catalog", DTOutput(ns("catalog")))),
      selectizeInput(ns("compare_pairs"),"Compare up to five combinations",choices=NULL,multiple=TRUE,
          options=list(placeholder="Choose combinations to compare their product lengths",maxItems=5,plugins=list("remove_button"))),
      uiOutput(ns("comparison")),
      div(class="fungal-evidence", uiOutput(ns("evidence"))),
      div(class="d-flex flex-wrap gap-2 mt-3",
          downloadButton(ns("download_pairs"),"Current combinations (CSV)",class="btn-outline-primary btn-sm"),
          downloadButton(ns("download_sequences"),"Selected reference sequences (FASTA)",class="btn-outline-secondary btn-sm"),
          downloadButton(ns("download_alignment"),"Frozen alignment (FASTA)",class="btn-outline-secondary btn-sm")),
      p(class="tiny mt-3", "Exact IUPAC-compatible sites establish candidate geometry on these templates. Unmatched sites remain unresolved. Length eligibility uses raw template products including primers, excluding adapters; colors describe the observed median and badges give per-template counts. This selected panel does not establish broad fungal amplification performance."))
}

fungalExplorerServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    index <- fungal_read_index()
    updateSelectInput(session,"version",choices=setNames(index$id,paste0(index$label," · ",substr(sub("^v[0-9]+-","",index$id),1,8))),selected=tail(index$id,1))
    selected_pair <- reactiveVal(NULL)
    snapshot <- reactive({
      validate(need(nrow(index)>0,"No complete fungal snapshots are installed."))
      version <- input$version
      if (is.null(version)) version <- tail(index$id,1)
      validate(need(version %in% index$id,"The selected historical snapshot is unavailable."))
      fungal_read_snapshot(version)
    })
    selected_accessions <- reactive({
      rows <- snapshot()$taxonomy
      if (!is.null(input$phylum) && input$phylum != "ALL") rows <- rows[!is.na(rows$phylum) & rows$phylum == input$phylum,,drop=FALSE]
      if (!is.null(input$class) && input$class != "ALL") rows <- rows[!is.na(rows$class) & rows$class == input$class,,drop=FALSE]
      rows$accession
    })
    sequence_rows <- reactive({
      rows <- snapshot()$pair_sequences
      rows[rows$accession %in% selected_accessions(),,drop=FALSE]
    })
    all_pairs <- reactive({
      data <- snapshot()
      pairs <- fungal_subset_summary(data$pair_sequences,data$pairs,selected_accessions())
      ceilings <- c(input$illumina_max %||% 500,input$pacbio_max %||% 2000,input$nanopore_max %||% 3500)
      validate(need(all(is.finite(ceilings)) && all(diff(ceilings)>=0),"Length ceilings must be ordered Illumina ≤ PacBio ≤ Nanopore."))
      fungal_platform_summary(sequence_rows(),pairs,input$platform %||% "ALL",ceilings[1],ceilings[2],ceilings[3])
    })
    candidates <- reactive({
      rows <- all_pairs()
      if (!length(selected_accessions())) return(rows[FALSE,,drop=FALSE])
      rows <- rows[rows$n_scorable>0 | (isTRUE(input$show_unresolved) & rows$evidence_status=="documented_pair"),,drop=FALSE]
      if (!identical(input$platform %||% "ALL","ALL")) {
        rows <- rows[rows$n_length_eligible>0 | (isTRUE(input$show_unresolved) & rows$n_scorable==0 & rows$evidence_status=="documented_pair"),,drop=FALSE]
      }
      if (!isTRUE(input$show_unresolved)) rows <- rows[rows$n_scorable>0,,drop=FALSE]
      if (!is.null(input$locus) && input$locus!="ALL") rows <- rows[rows$region==input$locus,,drop=FALSE]
      if (length(input$forward)) rows <- rows[rows$forward_primer %in% input$forward,,drop=FALSE]
      if (length(input$reverse)) rows <- rows[rows$reverse_primer %in% input$reverse,,drop=FALSE]
      query <- trimws(input$search %||% "")
      if (nzchar(query)) {
        text <- paste(rows$pair_label,rows$target,rows$source_citation)
        rows <- rows[grepl(tolower(query),tolower(text),fixed=TRUE),,drop=FALSE]
      }
      if (isTRUE(input$examples_only)) rows <- fungal_starter_examples(rows)
      rows
    })
    observeEvent(snapshot(), {
      data <- snapshot()
      values <- function(column) {
        value <- unique(as.character(data$taxonomy[[column]]))
        value <- sort(value[!is.na(value) & nzchar(value)])
        c("All represented groups"="ALL",setNames(value,value))
      }
      updateSelectInput(session,"phylum",choices=values("phylum"),selected="ALL")
      updateSelectInput(session,"class",choices=values("class"),selected="ALL")
      loci <- sort(unique(data$pairs$region))
      updateSelectInput(session,"locus",choices=c("All rDNA loci"="ALL",setNames(loci,loci)),selected="ALL")
      updateSelectizeInput(session,"forward",choices=sort(unique(data$pairs$forward_primer)),selected=character())
      updateSelectizeInput(session,"reverse",choices=sort(unique(data$pairs$reverse_primer)),selected=character())
      selected_pair(NULL)
    })
    observeEvent(list(snapshot(),input$phylum), {
      taxonomy <- snapshot()$taxonomy
      if (!is.null(input$phylum) && input$phylum != "ALL") taxonomy <- taxonomy[!is.na(taxonomy$phylum) & taxonomy$phylum==input$phylum,,drop=FALSE]
      classes <- sort(unique(as.character(taxonomy$class)))
      classes <- classes[!is.na(classes) & nzchar(classes)]
      selected <- if(!is.null(input$class) && input$class %in% classes) input$class else "ALL"
      updateSelectInput(session,"class",choices=c("All represented classes"="ALL",setNames(classes,classes)),selected=selected)
    })
    observeEvent(input$reset, {
      for (name in c("phylum","class","locus")) updateSelectInput(session,name,selected="ALL")
      for (name in c("forward","reverse")) updateSelectizeInput(session,name,selected=character())
      updateTextInput(session,"search",value="")
      updateRadioButtons(session,"platform",selected="ALL")
      updateCheckboxInput(session,"examples_only",value=TRUE)
      updateCheckboxInput(session,"show_unresolved",value=TRUE)
    })
    observeEvent(input$pair_click, {
      pair <- input$pair_click$pair
      if (length(pair)==1 && pair %in% candidates()$pair_id) selected_pair(pair)
    })
    observe({
      rows <- candidates()
      choices <- setNames(rows$pair_id,rows$pair_label)
      current <- isolate(input$compare_pairs)
      if (is.null(current)) {
        targets <- c("ITS3 + ITS4","ITS1 + ITS4","FF700 + LR7")
        current <- rows$pair_id[match(targets,rows$pair_label,nomatch=0)]
      }
      updateSelectizeInput(session,"compare_pairs",choices=choices,selected=intersect(current,rows$pair_id))
    })
    active_pair <- reactive({
      rows <- candidates()
      if (!nrow(rows)) return(rows)
      current <- selected_pair()
      if (is.null(current) || !current %in% rows$pair_id) current <- rows$pair_id[1]
      rows[match(current,rows$pair_id),,drop=FALSE]
    })
    output$snapshot_meta <- renderUI({
      data <- snapshot()
      counts <- data$taxonomy
      class_count <- length(setdiff(unique(na.omit(counts$class)),c("unresolved","")))
      div(strong(length(data$sequences)," templates · ",class_count," resolved classes"),
          div(class="tiny",data$manifest$id," · ",data$manifest$created_at),
          div(class="tiny",data$manifest$transport_source %||% "Bundled pinned snapshot",
              if(isTRUE(attr(index,"stale"))) paste0(" · ",attr(index,"fallback_reason"))),
          div(class="tiny",data$manifest$description))
    })
    output$discovery_status <- renderUI({
      rows <- all_pairs()
      shown <- candidates()
      counts <- c(Illumina=sum(rows$n_illumina_eligible>0),PacBio=sum(rows$n_pacbio_eligible>0),Nanopore=sum(rows$n_nanopore_eligible>0))
      div(class="fungal-discovery-status",
          strong(nrow(shown)," combinations shown · ",length(selected_accessions())," selected templates"),
          div(class="tiny",paste0(names(counts),": ",counts," pairs with at least one length-eligible template",collapse=" · ")),
          div(class="tiny","Platform counts are nested; shorter products remain eligible on the longer-read systems. Each pair reports supported / selected templates."))
    })
    pair_button <- function(row, text=row$pair_label) {
      tags$button(type="button",class="fungal-pair-choice",`data-input`=session$ns("pair_click"),`data-pair`=row$pair_id,
                  title=paste(row$pair_label,row$platforms),text)
    }
    length_text <- function(rows) ifelse(is.finite(rows$amplicon_median_bp),
        paste0(format(rows$amplicon_median_bp,trim=TRUE)," bp (",rows$amplicon_min_bp,"–",rows$amplicon_max_bp,")"),"Unresolved")
    output$pairs <- renderDT({
      rows <- candidates()
      display <- data.frame(
          Pair=vapply(seq_len(nrow(rows)),function(i) as.character(pair_button(rows[i,,drop=FALSE])),character(1)),
          `Median (range)`=length_text(rows),
          `Median tier`=paste0('<span class="fungal-platform-badge ',tolower(rows$length_tier),'">',htmltools::htmlEscape(rows$length_tier),'</span>'),
          `Exact support`=paste0(rows$n_scorable," / ",rows$n_panel),
          `Length eligibility`=rows$platforms,
          Evidence=fungal_evidence_label(rows$evidence_status),
          check.names=FALSE)
      datatable(display,escape=FALSE,rownames=FALSE,class="compact stripe",selection="none",
                options=list(pageLength=20,lengthMenu=c(20,50,100),scrollX=TRUE,order=list()))
    },server=FALSE)
    output$matrix <- renderUI({
      rows <- candidates()
      if (!nrow(rows)) return(p("No pairs match these filters."))
      forwards <- head(unique(rows$forward_primer),24)
      reverses <- head(unique(rows$reverse_primer),24)
      div(class="fungal-matrix",
          p(class="tiny","Choose a product to inspect its template geometry. Empty cells have no supported combination in this view. For a larger search, turn off starter examples and select primers."),
          tags$table(class="table table-sm",
              tags$thead(tags$tr(tags$th("Forward ↓ / reverse →"),lapply(reverses,tags$th))),
              tags$tbody(lapply(forwards,function(forward) {
                tags$tr(tags$th(scope="row",forward),lapply(reverses,function(reverse) {
                  pair <- rows[rows$forward_primer==forward & rows$reverse_primer==reverse,,drop=FALSE]
                  if (!nrow(pair)) return(tags$td("—"))
                  tags$td(style=paste0("border-bottom:3px solid ",pair$tier_color[1]),
                      pair_button(pair[1,,drop=FALSE],if (is.finite(pair$amplicon_median_bp[1])) paste0(pair$amplicon_median_bp[1]," bp") else "Unresolved"))
                }))
              }))))
    })
    output$catalog <- renderDT({
      rows <- snapshot()$catalog
      display <- rows[,intersect(c("primer_name","sequence","direction","source_direction","direction_note","gene_locus","target","remarks","reference","primary_reference_status"),names(rows)),drop=FALSE]
      datatable(display,rownames=FALSE,options=list(pageLength=10,scrollX=TRUE),class="compact stripe")
    })
    output$version_changes <- renderUI({
      current <- snapshot()
      if (nrow(index)<2) return(p("One complete snapshot is available."))
      baseline <- fungal_read_snapshot(index$id[1])
      added <- setdiff(names(current$sequences),names(baseline$sequences))
      removed <- setdiff(names(baseline$sequences),names(current$sequences))
      new_classes <- setdiff(unique(na.omit(current$taxonomy$class)),unique(na.omit(baseline$taxonomy$class)))
      new_classes <- new_classes[nzchar(new_classes) & tolower(new_classes)!="unresolved"]
      common <- intersect(current$pairs$pair_id,baseline$pairs$pair_id)
      before <- baseline$pairs$amplicon_median_bp[match(common,baseline$pairs$pair_id)]
      after <- current$pairs$amplicon_median_bp[match(common,current$pairs$pair_id)]
      changed <- sum((is.na(before)!=is.na(after)) | (!is.na(before)&!is.na(after)&before!=after))
      div(p(strong("Compared with ",baseline$manifest$label)),
          p(length(added)," added accessions · ",length(removed)," removed accessions · ",changed," shared pairs with changed median length or resolved length evidence."),
          p("New represented classes: ",if(length(new_classes)) paste(new_classes,collapse=", ") else "none"),
          p(class="tiny","V1 preserves the original coordinate reference and eight recorded geometries. Its exact-site explorer results were derived when this baseline was captured; low-confidence legacy spans are retained as provenance and remain unresolved for length eligibility."),
          tags$details(tags$summary("Accessions added in the selected snapshot"),tags$pre(paste(added,collapse="\n"))))
    })
    output$comparison <- renderUI({
      rows <- candidates()
      rows <- rows[rows$pair_id %in% input$compare_pairs,,drop=FALSE]
      if(!nrow(rows)) return(NULL)
      limit <- max(input$nanopore_max %||% 3500,rows$amplicon_max_bp,na.rm=TRUE)
      scale <- function(bp) 230 + bp/limit*620
      height <- 80+nrow(rows)*65
      div(class="fungal-map",
          tags$svg(viewBox=paste("0 0 1000",height),role="img",`aria-label`="Comparison of selected amplicon product lengths",
              tags$text(x=10,y=20,`font-size`=13,"Product lengths: median bars and observed min–max ranges"),
              lapply(seq_len(nrow(rows)),function(i) {
                  pair <- rows[i,,drop=FALSE]
                  y <- 45+i*65
                  tagList(tags$text(x=10,y=y+4,`font-size`=12,pair$pair_label),
                      if(is.finite(pair$amplicon_median_bp)) tagList(
                          tags$line(x1=scale(0),x2=scale(pair$amplicon_median_bp),y1=y,y2=y,stroke=pair$tier_color,`stroke-width`=10),
                          tags$line(x1=scale(pair$amplicon_min_bp),x2=scale(pair$amplicon_max_bp),y1=y-12,y2=y-12,stroke=pair$tier_color,`stroke-width`=2),
                          tags$text(x=870,y=y+4,`font-size`=12,paste0(pair$amplicon_median_bp," bp")),
                          tags$text(x=230,y=y+25,`font-size`=10,paste0(pair$n_scorable," / ",pair$n_panel," exact-support templates · ",pair$length_tier)))
                      else tags$text(x=230,y=y,`font-size`=12,"Unresolved length evidence"))
              }),
              lapply(c(input$illumina_max %||% 500,input$pacbio_max %||% 2000,input$nanopore_max %||% 3500),function(bp) tagList(
                  tags$line(x1=scale(bp),x2=scale(bp),y1=40,y2=height-5,stroke="#94a3b8",`stroke-dasharray`="3 5"),
                  tags$text(x=scale(bp),y=35,`font-size`=10,paste0(bp," bp"))))))
    })
    output$evidence <- renderUI({
      pair <- active_pair()
      if (!nrow(pair)) return(div(class="callout","No combinations match these filters. Reset filters or expand the candidate list."))
      data <- snapshot()
      rows <- sequence_rows()
      rows <- rows[rows$pair_id==pair$pair_id,,drop=FALSE]
      supported <- rows[rows$status=="scorable",,drop=FALSE]
      examples <- head(supported[order(supported$amplicon_bp),,drop=FALSE],6)
      components <- data$catalog[(data$catalog$primer_name==pair$forward_primer & data$catalog$sequence==pair$forward_sequence) |
                                  (data$catalog$primer_name==pair$reverse_primer & data$catalog$sequence==pair$reverse_sequence),,drop=FALSE]
      div(h4(pair$pair_label),
          p(fungal_evidence_label(pair$evidence_status)," · ",pair$region," · component targets: ",pair$target),
          p(strong(length_text(pair))," · exact support ",pair$n_scorable," / ",pair$n_panel," selected templates"),
          p(pair$platforms),
          p(class="tiny", "Gene landmarks are shown only where NCBI supplies an unambiguous complete boundary. Missing or partial annotations do not change the product length measured between exact primer sites."),
          if (isTRUE(pair$length_range_crosses_tiers)) div(class="callout","Template lengths span multiple tiers: ",pair$length_tier_range,". Inspect the per-template counts before choosing a ceiling."),
          div(class="fungal-map", lapply(seq_len(nrow(examples)),function(i) {
              row <- examples[i,,drop=FALSE]
              sequence <- data$sequences[[row$accession]]
              start <- max(1,row$forward_start-100)
              end <- min(nchar(sequence),row$reverse_end+100)
              scale <- function(bp) 130+(bp-start)/max(1,end-start)*650
              template_tier <- if(row$amplicon_bp <= (input$illumina_max %||% 500)) "Illumina" else if(row$amplicon_bp <= (input$pacbio_max %||% 2000)) "PacBio" else if(row$amplicon_bp <= (input$nanopore_max %||% 3500)) "Nanopore" else "Beyond limits"
              color <- c(Illumina="#2563eb",PacBio="#7c3aed",Nanopore="#0f766e",`Beyond limits`="#64748b")[[template_tier]]
              tags$svg(viewBox="0 0 920 95",role="img",`aria-label`=paste(row$accession,pair$pair_label,row$amplicon_bp,"bp",template_tier),
                  tags$text(x=8,y=30,`font-size`=12,row$accession),
                  tags$line(x1=130,x2=780,y1=32,y2=32,stroke="#cbd5e1",`stroke-width`=2),
                  tags$line(x1=scale(row$forward_start),x2=scale(row$reverse_end),y1=32,y2=32,stroke=color,`stroke-width`=8),
                  tags$text(x=scale(row$forward_start),y=17,`font-size`=16,"→"),
                  tags$text(x=scale(row$reverse_end),y=17,`font-size`=16,"←"),
                  tags$text(x=800,y=35,`font-size`=12,paste0(row$amplicon_bp," bp")),
                  tags$text(x=130,y=60,`font-size`=10,paste0("Template positions ",start,"–",end," · ",template_tier)),
                  lapply(seq_len(nrow(data$landmarks)),function(j) {
                    landmark <- data$landmarks[j,,drop=FALSE]
                    if (!"accession" %in% names(landmark) || is.na(landmark$accession) || landmark$accession!=row$accession || landmark$end<start || landmark$start>end) return(NULL)
                    tags$text(x=scale(max(start,landmark$start)),y=84,`font-size`=10,landmark$label)
                  }))
          })),
          if(!nrow(supported)) div(class="callout","No unique inward exact-site product is supported by the selected templates. Legacy reference-specific lengths are withheld."),
          tags$details(tags$summary("Primer sequences, sources, and missing evidence"),
              p("Forward (5′→3′): ",tags$code(pair$forward_sequence)),
              p("Reverse (5′→3′): ",tags$code(pair$reverse_sequence)),
              p(pair$source_citation),
              if("remarks" %in% names(components)) lapply(seq_len(nrow(components)),function(i) {
                remark <- components$remarks[i]
                if(!is.na(remark) && nzchar(remark)) p(strong(components$primer_name[i],": "),remark)
              }),
              lapply(strsplit(pair$source_url," | ",fixed=TRUE)[[1]],function(url) if(grepl("^https://",url)) tags$a(href=url,target="_blank",rel="noopener","Component provenance ↗ ")),
              tags$pre(paste(names(table(rows$status)),as.integer(table(rows$status)),sep=": ",collapse="\n"))),
          DTOutput(session$ns("sequence_evidence")))
    })
    output$sequence_evidence <- renderDT({
      pair <- active_pair()
      rows <- sequence_rows()
      rows <- rows[rows$pair_id %in% pair$pair_id,,drop=FALSE]
      taxonomy <- snapshot()$taxonomy
      for (column in intersect(c("annotation_status","annotation_note"),names(taxonomy))) rows[[column]] <- taxonomy[[column]][match(rows$accession,taxonomy$accession)]
      datatable(rows[,intersect(c("accession","organism","phylum","class","amplicon_bp","status","forward_status","reverse_status","annotation_status","annotation_note"),names(rows)),drop=FALSE],
                rownames=FALSE,options=list(pageLength=5,scrollX=TRUE),class="compact stripe")
    })
    output$download_pairs <- downloadHandler(filename=function() paste0(snapshot()$manifest$id,"-combinations.csv"),content=function(file) {
      rows <- candidates()
      rows$snapshot_id <- snapshot()$manifest$id
      rows$illumina_ceiling <- input$illumina_max %||% 500
      rows$pacbio_ceiling <- input$pacbio_max %||% 2000
      rows$nanopore_ceiling <- input$nanopore_max %||% 3500
      write.csv(rows,file,row.names=FALSE,na="")
    })
    write_fasta <- function(sequences,file) writeLines(unlist(lapply(seq_along(sequences),function(i) c(paste0(">",names(sequences)[i]),sequences[i])),use.names=FALSE),file)
    output$download_sequences <- downloadHandler(filename=function() paste0(snapshot()$manifest$id,"-templates.fasta"),content=function(file) write_fasta(snapshot()$sequences[selected_accessions()],file))
    output$download_alignment <- downloadHandler(filename=function() paste0(snapshot()$manifest$id,"-alignment.fasta"),content=function(file) write_fasta(snapshot()$alignment,file))
    output$download_manifest <- downloadHandler(filename=function() paste0(snapshot()$manifest$id,"-manifest.json"),content=function(file) jsonlite::write_json(snapshot()$manifest,file,pretty=TRUE,auto_unbox=TRUE,null="null",na="null"))
    list(snapshot=snapshot,candidates=candidates,all_pairs=all_pairs,selected_accessions=selected_accessions,active_pair=active_pair)
  })
}
