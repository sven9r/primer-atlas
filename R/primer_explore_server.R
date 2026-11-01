primerExploreServer <- function(id, marker, target_group, fungal_snapshot=NULL, root=".", on_evidence=NULL, on_version=NULL, on_fit=NULL, taxonomic_report=NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    base_catalog <- primer_mixer_catalog(root)
    state <- shiny::reactiveValues(key=NULL, forward=character(), reverse=character(),
      custom=data.frame(), supplied=NULL, cache=list(), message="")
    ack <- new.env(parent=emptyenv()); ack$forward <- ack$reverse <- NULL
    context <- shiny::reactive(paste(target_group(), marker(), sep="::"))
    active_marker <- shiny::reactive({ shiny::req(marker()); as.character(marker()) })
    catalog <- shiny::reactive({
      rows <- base_catalog[base_catalog$marker_id==active_marker(),,drop=FALSE]
      .primer_mixer_bind(list(rows,state$custom))
    })
    selection <- function(role) intersect(state[[role]],catalog()$primer_id[catalog()$direction==role])
    forward <- shiny::reactive(selection("forward"))
    reverse <- shiny::reactive(selection("reverse"))
    sync <- function() {
      rows <- shiny::isolate(catalog())
      for(role in c("forward","reverse")) {
        own <- rows[rows$direction==role,,drop=FALSE]
        scope <- ifelse(nzchar(own$target_scope),paste0(" · ",own$target_scope),"")
        labels <- paste0(own$primer_name,scope," · ",own$source_sequence)
        chosen <- intersect(shiny::isolate(state[[role]]),own$primer_id)
        state[[role]] <- chosen
        ack[[role]] <- chosen
        shiny::updateSelectizeInput(session,role,choices=stats::setNames(own$primer_id,labels),selected=chosen,server=FALSE)
      }
    }
    shiny::observeEvent(context(), {
      if(!is.null(shiny::isolate(state$key))) {
        previous <- shiny::isolate(state$cache)
        previous[[shiny::isolate(state$key)]] <- shiny::isolate(list(forward=state$forward,reverse=state$reverse,custom=state$custom,supplied=state$supplied))
        state$cache <- previous
      }
      saved <- shiny::isolate(state$cache[[context()]])
      state$key <- context()
      state$forward <- if(is.null(saved)) character() else saved$forward
      state$reverse <- if(is.null(saved)) character() else saved$reverse
      state$custom <- if(is.null(saved)) data.frame() else saved$custom
      state$supplied <- if(is.null(saved)) NULL else saved$supplied
      state$message <- ""
      sync()
    },priority=100,ignoreInit=FALSE)
    for(role in c("forward","reverse")) local({
      assigned <- role
      shiny::observeEvent(input[[assigned]], {
        values <- input[[assigned]]
        if(is.null(values)) values <- character()
        expected <- ack[[assigned]]
        if(!is.null(expected)) {
          if(!setequal(values,expected)) return()
          ack[[assigned]] <- NULL
        }
        rows <- catalog()
        state[[assigned]] <- head(intersect(values,rows$primer_id[rows$direction==assigned]),5L)
      },ignoreInit=TRUE,ignoreNULL=FALSE)
    })
    snapshot <- shiny::reactive({
      if(is.null(fungal_snapshot) || !(active_marker()=="ITS_FUNGAL" ||
          (target_group()=="FUNGI" && active_marker() %in% c("18S","28S")))) return(NULL)
      fungal_snapshot()
    })
    pinned_templates <- shiny::reactive({
      if(target_group()=="FUNGI" && active_marker() %in% c("18S","28S") && !is.null(snapshot()))
        return(primer_mixer_fungal_templates(snapshot(),active_marker()))
      if(active_marker()=="ITS_FUNGAL" && !is.null(snapshot())) {
        data <- snapshot(); taxonomy <- data$taxonomy
        ids <- names(data$sequences)
        return(data.frame(template_id=ids, sequence=unname(unlist(data$sequences,use.names=FALSE)),
          marker_id="ITS_FUNGAL",organism=as.character(taxonomy$organism[match(ids,taxonomy$accession)]),
          source_url=paste0("https://www.ncbi.nlm.nih.gov/nuccore/",ids),version=data$manifest$id,
          target_group="FUNGI",source_status="Frozen NCBI rDNA snapshot",stringsAsFactors=FALSE))
      }
      group <- target_group()
      if(group %in% c("ALL","METAZOA")) return(primer_mixer_templates(root,active_marker(),"ALL"))
      primer_mixer_templates(root,active_marker(),group)
    })
    templates <- shiny::reactive({
      if(!is.null(state$supplied)) state$supplied else pinned_templates()
    })
    selected_reference <- shiny::reactive({
      rows <- templates()
      if(!nrow(rows)) return(NULL)
      value <- input$reference
      if(length(value)==1L && value %in% rows$template_id) return(value)
      if(active_marker()=="ITS_FUNGAL" && is.null(state$supplied) && !is.null(snapshot())) {
        preferred <- fungal_primer_overview_default_accession(snapshot())
        if(preferred %in% rows$template_id) return(preferred)
      }
      if(target_group()=="FUNGI" && active_marker() %in% c("18S","28S") && is.null(state$supplied) && "partial" %in% names(rows))
        return(rows$template_id[order(rows$partial %in% TRUE,-nchar(rows$sequence))][1L])
      if(target_group()=="FISH" && active_marker()=="18S" && is.null(state$supplied)) {
        predicted <- if("source_status" %in% names(rows)) grepl("model|predicted",rows$source_status,ignore.case=TRUE) else rep(FALSE,nrow(rows))
        return(rows$template_id[order(predicted,-nchar(rows$sequence))][1L])
      }
      rows$template_id[1L]
    })
    shiny::observeEvent(templates(), {
      rows <- templates()
      predicted <- if("source_status" %in% names(rows)) grepl("MODEL",rows$source_status) else rep(FALSE,nrow(rows))
      partial <- if("partial" %in% names(rows)) rows$partial %in% TRUE else rep(FALSE,nrow(rows))
      labels <- if(nrow(rows)) paste0(rows$organism," · ",rows$template_id,
        ifelse(predicted," · predicted model",""),ifelse(partial," · partial","")) else character()
      selected <- shiny::isolate(selected_reference())
      shiny::updateSelectInput(session,"reference",choices=stats::setNames(rows$template_id,labels),selected=selected)
    },ignoreInit=FALSE)
    shiny::observeEvent(input$use_fasta, {
      tryCatch({
        value <- primer_mixer_parse_fasta(input$fasta,active_marker())
        if(!nrow(value)) stop("Paste at least one FASTA template first.")
        state$supplied <- value
        shiny::showNotification(paste(nrow(value),"session templates loaded."),type="message")
      },error=function(e) shiny::showNotification(conditionMessage(e),type="error",duration=8))
    })
    shiny::observeEvent(input$reset_templates,{ state$supplied <- NULL })
    annotations <- shiny::reactive({
      id <- selected_reference(); rows <- templates()
      if(is.null(id)) return(data.frame())
      if(active_marker()=="ITS_FUNGAL" && is.null(state$supplied) && !is.null(snapshot()))
        return(.fungal_atlas_landmarks(snapshot(),id,nchar(rows$sequence[match(id,rows$template_id)])))
      if(!is.null(state$supplied)) return(data.frame(label="Supplied sequence · marker assigned by user",
        start=1L,end=nchar(rows$sequence[match(id,rows$template_id)]),stringsAsFactors=FALSE))
      if(target_group()=="FUNGI" && active_marker() %in% c("18S","28S") && !is.null(snapshot()) &&
          all(c("source_start","source_end","accession_version") %in% names(rows))) {
        row <- rows[rows$template_id==id,,drop=FALSE]
        gene <- data.frame(landmark_id="extracted_interval",label=paste(active_marker(),"source interval",if(isTRUE(row$partial)) "· partial" else ""),
          start=1L,end=nchar(row$sequence),landmark_type="gene",source_note=row$source_note,stringsAsFactors=FALSE)
        source <- snapshot(); parent <- row$accession_version
        own <- .fungal_atlas_landmarks(source,parent,nchar(source$sequences[[parent]]))
        own <- own[own$is_domain & own$end>=row$source_start & own$start<=row$source_end,,drop=FALSE]
        if(nrow(own)) {
          own$parent_start <- own$start; own$parent_end <- own$end
          own$start <- pmax(own$start,row$source_start)-row$source_start+1L
          own$end <- pmin(own$end,row$source_end)-row$source_start+1L
          own$label <- paste0(own$label,ifelse(own$parent_start<row$source_start | own$parent_end>row$source_end," · partial displayed interval",""))
        }
        return(.primer_mixer_bind(list(gene,own)))
      }
      label <- switch(active_marker(),COI="COI / COX1",`12S_MT`="12S rRNA",`18S`="18S rRNA",`28S`="28S rRNA",
        `16S_MT`="Mitochondrial 16S",`16S_PROK`="16S rRNA",active_marker())
      data.frame(label=label,start=1L,end=nchar(rows$sequence[match(id,rows$template_id)]),stringsAsFactors=FALSE)
    })
    output$region_ui <- shiny::renderUI({
      choices <- c("Full supplied marker"="ALL")
      if(active_marker()=="ITS_FUNGAL" && is.null(state$supplied)) {
        choices <- c("ITS1 / 5.8S / ITS2"="ITS","18S / SSU"="SSU","28S / LSU"="LSU","All rDNA"="ALL")
        own <- annotations()
        if(nrow(own) && "landmark_id" %in% names(own)) choices <- c(choices,
          stats::setNames(paste0("landmark:",own$landmark_id),paste0(own$label," · annotated interval")))
      }
      if(target_group()=="FUNGI" && active_marker() %in% c("18S","28S") && is.null(state$supplied)) {
        own <- annotations()
        if(nrow(own) && "landmark_id" %in% names(own)) choices <- c(choices,
          stats::setNames(paste0("landmark:",own$landmark_id),paste0(own$label," · annotated interval")))
      }
      current <- shiny::isolate(input$region)
      if(is.null(current) || !current %in% choices) current <- unname(choices[1L])
      shiny::tagList(shiny::selectInput(ns("region"),"Region",choices,selected=current),
        shiny::checkboxInput(ns("all_sites"),"Show all mapped primer labels",value=FALSE))
    })
    # Closed map settings must still initialize the region used by the overview.
    shiny::outputOptions(output,"region_ui",suspendWhenHidden=FALSE)
    domain <- shiny::reactive({
      id <- selected_reference(); rows <- templates()
      if(is.null(id)) return(NULL)
      size <- nchar(rows$sequence[match(id,rows$template_id)])
      region <- input$region %||% if(active_marker()=="ITS_FUNGAL") "ITS" else "ALL"
      if(!is.null(region) && startsWith(region,"landmark:") && is.null(state$supplied)) {
        own <- annotations(); own <- own[own$landmark_id==sub("^landmark:","",region),,drop=FALSE]
        if(nrow(own)!=1L) return(NULL)
        return(c(max(1,own$start-1000),min(size,own$end+2000)))
      }
      if(active_marker()!="ITS_FUNGAL" || is.null(region) || region=="ALL" || !is.null(state$supplied)) return(c(1,size))
      labels <- switch(region,ITS=c("ITS1","5.8S","ITS2"),SSU=c("18S","SSU"),LSU=c("28S","LSU"),character())
      own <- annotations(); own <- own[own$label %in% labels,,drop=FALSE]
      if(!nrow(own)) return(NULL)
      c(max(1,min(own$start)-1000),min(size,max(own$end)+2000))
    })
    max_mismatches <- shiny::reactive({
      value <- suppressWarnings(as.integer(input$mismatches))
      if(length(value)!=1L || is.na(value)) 0L else value
    })
    length_limits <- shiny::reactive({
      ceiling <- switch(input$platform %||% "ALL",ILLUMINA=500,PACBIO=2000,NANOPORE=3500,Inf)
      lower <- input$min_bp %||% 0; upper <- input$max_bp %||% 10000
      if(!is.finite(lower) || !is.finite(upper) || lower<0 || upper<lower) stop("Choose an increasing, nonnegative product-length range.")
      # An incompatible platform ceiling is an empty product filter, not an
      # invalid search or a change to the underlying template evidence.
      c(lower,upper,ceiling)
    })
    result <- shiny::reactive({
      limits <- length_limits()
      value <- primer_mixer_build(catalog(),forward(),reverse(),templates(),active_marker(),
        max_mismatches=max_mismatches(),min_bp=limits[1],max_bp=limits[2])
      if(nrow(value$products)) {
        value$products$passes_length_filter <- value$products$passes_length_filter &
          !is.na(value$products$amplicon_bp) & value$products$amplicon_bp<=limits[3]
        value$summary$n_length_eligible <- vapply(value$summary$pair_id,function(pair)
          sum(value$products$passes_length_filter[value$products$pair_id==pair]),integer(1))
      }
      value$settings$platform_ceiling_bp <- limits[3]
      value
    })
    output$scope <- shiny::renderUI({
      rows <- templates()
      shiny::span(class="atlas-status-chip",if(!is.null(state$supplied)) "Session templates" else "Reference examples",
        " · ",nrow(rows)," templates · ",nrow(catalog())," catalog / custom oligos")
    })
    output$version_ui <- shiny::renderUI({
      if(is.null(snapshot()) || is.null(on_version) || !is.null(state$supplied)) return(NULL)
      versions <- fungal_read_index(remote=FALSE)
      shiny::selectInput(ns("snapshot_version"),"Frozen fungal reference version",stats::setNames(versions$id,versions$label),selected=snapshot()$manifest$id)
    })
    shiny::observeEvent(input$snapshot_version, {
      if(!is.null(snapshot()) && !is.null(on_version) && input$snapshot_version!=snapshot()$manifest$id)
        on_version(input$snapshot_version)
    },ignoreInit=TRUE)
    output$template_sources <- shiny::renderUI({
      rows <- templates()
      if(!nrow(rows)) return(shiny::p(attr(rows,"unavailable_reason") %||%
        "No marker-matched starter templates are bundled for this target group. Add FASTA templates to investigate it."))
      shiny::tagList(shiny::p(class="atlas-context-line","These are explicit sequence examples, not a representative taxonomic coverage panel. Product coordinates start at base 1 of each extracted or supplied marker. Predicted and partial sources retain those labels."),
        shiny::tags$ul(lapply(seq_len(nrow(rows)),function(i) shiny::tags$li(
          if(nzchar(rows$source_url[i])) shiny::tags$a(href=rows$source_url[i],target="_blank",rel="noopener",rows$template_id[i]) else rows$template_id[i],
          " · ",rows$organism[i]," · ",nchar(rows$sequence[i])," bp · ",rows$version[i],
          if("source_status" %in% names(rows)) paste0(" · ",rows$source_status[i]),
          if("partial" %in% names(rows) && isTRUE(rows$partial[i])) " · partial source record"))))
    })
    site_data <- shiny::reactive({
      id <- selected_reference(); rows <- templates()
      if(is.null(id)) return(data.frame())
      components <- catalog()
      components <- components[nchar(components$sequence)>=8L & nchar(components$sequence)<=120L,,drop=FALSE]
      own <- rows[rows$template_id==id,,drop=FALSE]
      limit <- max(1L,min(500L,floor(2500000/nchar(own$sequence))))
      n_total <- nrow(components)
      if(n_total>limit) {
        selected <- which(components$primer_id %in% c(forward(),reverse()))
        chosen <- unique(c(selected,seq_len(n_total)))
        components <- components[head(chosen,limit),,drop=FALSE]
      }
      sites <- primer_mixer_sites(components,own,active_marker(),max_mismatches())
      attr(sites,"scan_notice") <- if(n_total>limit) paste0("Overview search is limited to ",nrow(components),
        " of ",n_total," components on this long template, including your selected primers. The full catalog remains searchable below.") else ""
      sites
    })
    output$overview <- shiny::renderUI({
      id <- selected_reference()
      if(is.null(id)) return(shiny::div(class="callout",attr(templates(),"unavailable_reason") %||%
        "Add FASTA templates below to draw this marker and test combinations."))
      if(is.null(domain())) return(shiny::div(class="callout","This template has no verified boundary for the selected region. Choose another region or template."))
      sites <- site_data(); scan_notice <- attr(sites,"scan_notice"); rows <- catalog()
      sites <- sites[sites$status %in% c("exact_site","mismatch_site"),,drop=FALSE]
      sites$primer_name <- rows$primer_name[match(sites$primer_id,rows$primer_id)]
      sites$match_type <- sites$status
      variants <- rows[rows$catalog_origin=="derived_reverse_complement",,drop=FALSE]
      active_variants <- variants[variants$primer_id %in% c(forward(),reverse()),,drop=FALSE]
      if(nrow(variants)) sites <- sites[!sites$primer_id %in% variants$primer_id | sites$primer_id %in% active_variants$primer_id,,drop=FALSE]
      if(nrow(active_variants)) sites <- sites[!sites$primer_id %in% active_variants$parent_primer_id,,drop=FALSE]
      sites <- sites[sites$end>=domain()[1] & sites$start<=domain()[2],,drop=FALSE]
      n_sites <- nrow(sites)
      per_direction <- if((input$width %||% 650)<800) 5L else 7L
      if(!isTRUE(input$all_sites) && n_sites>2L*per_direction) {
        picked <- unlist(lapply(c("forward","reverse"),function(role) {
          indices <- which(sites$direction==role)
          indices <- indices[order(sites$start[indices],sites$primer_name[indices])]
          if(length(indices)>per_direction) indices[unique(round(seq(1,length(indices),length.out=per_direction)))] else indices
        }))
        picked <- unique(c(picked,which(sites$primer_id %in% c(forward(),reverse()))))
        sites <- sites[picked,,drop=FALSE]
      }
      shiny::tagList(primer_site_overview(sites,id,nchar(templates()$sequence[match(id,templates()$template_id)]),
        active_marker(),click_input=ns("primer_click"),selected_forward=forward(),selected_reverse=reverse(),
        viewport_width=input$width %||% 650,annotations=annotations(),domain=domain()),
        shiny::p(class="atlas-context-line",nrow(sites)," of ",n_sites," uniquely placed primers shown in this region. Search the full marker catalog in Mix & match.",
          if(n_sites==0L) " No unique sites meet the current sequence-match setting; this does not establish PCR failure."),
        if(nzchar(scan_notice %||% "")) shiny::p(class="atlas-context-line",scan_notice))
    })
    output$selection_count <- shiny::renderUI(shiny::span(class="atlas-status-chip",
      length(forward())," forward × ",length(reverse())," reverse = ",length(forward())*length(reverse())," combinations"))
    shiny::observeEvent(input$clear,{ state$forward <- state$reverse <- character(); sync() })
    load_example <- function() {
      rows <- catalog(); group <- target_group()
      preferred <- switch(active_marker(),
        `12S_MT`=if(group=="BIRDS") list(f=c("MiBird-U-F","12S-V5f"),r=c("MiBird-U-R","12S-V5r")) else
          list(f=c("MiFish-U-F","MiFish-E-F","L1848 / teleo_F","12S-V5f"),r=c("MiFish-U-R","MiFish-E-R","H1913 / teleo_R","12S-V5r")),
        COI=if(group=="FISH") list(f=c("FishF1","FishF2","coi.175f"),r=c("FishR1","FishR2","coi.226r")) else if(group=="BIRDS")
          list(f=c("BirdF1"),r=c("BirdR1","BirdR2","BirdR3")) else
          list(f=c("mlCOIintF","BF2","fwhF2"),r=c("Fol-degen-rev","BR2","fwhR2")),
        `18S`=list(f=c("TAReuk454FWD1","1380F","1389F"),r=c("TAReukREV3","1510R")),
        `28S`=if(group=="ARTHROPODS") list(f=c("D2F"),r=c("D3R")) else list(f=c("C1′"),r=c("C2")),
        ITS_FUNGAL=list(f=c("ITS1","ITS1F","ITS5","ITS3","fITS7"),r=c("ITS2","ITS4","ITS4B","LR21","LR3")),list(f=character(),r=character()))
      choose <- function(role,names) {
        own <- rows[rows$direction==role,,drop=FALSE]
        ids <- unlist(lapply(names,function(name) head(own$primer_id[own$primer_name==name],1L)),use.names=FALSE)
        if(!length(ids)) ids <- head(own$primer_id,3L)
        unique(ids)
      }
      state$forward <- choose("forward",preferred$f); state$reverse <- choose("reverse",preferred$r)
      sync()
    }
    shiny::observeEvent(input$example,load_example())
    shiny::observeEvent(input$add_custom, {
      tryCatch({
        row <- primer_mixer_custom(input$custom_name,input$custom_sequence,input$custom_direction,active_marker())
        role <- row$direction
        if(length(shiny::isolate(state[[role]]))>=5L) stop("Remove one selected primer in this direction before adding another.")
        custom <- shiny::isolate(state$custom)
        if(nrow(custom)) custom <- custom[custom$primer_id!=row$primer_id,,drop=FALSE]
        state$custom <- .primer_mixer_bind(list(custom,row))
        state[[role]] <- unique(c(shiny::isolate(state[[role]]),row$primer_id))
        state$message <- paste(row$primer_name,"added as a session candidate.",row$normalization_note)
        sync()
      },error=function(e) shiny::showNotification(conditionMessage(e),type="error",duration=8))
    })
    output$custom_status <- shiny::renderUI(shiny::p(class="atlas-context-line",state$message))
    shiny::observeEvent(input$primer_click, {
      id <- input$primer_click$id
      rows <- catalog(); own <- rows[rows$primer_id %in% id,,drop=FALSE]
      if(nrow(own)!=1L) return()
      if(own$catalog_origin=="derived_reverse_complement") {
        state$forward <- setdiff(state$forward,id); state$reverse <- setdiff(state$reverse,id)
        sync(); return()
      }
      variant <- primer_mixer_reverse(own)
      variant_id <- variant$primer_id
      role <- own$direction; opposite <- if(role=="forward") "reverse" else "forward"
      if(variant_id %in% c(state$forward,state$reverse)) {
        state$forward <- setdiff(state$forward,variant_id); state$reverse <- setdiff(state$reverse,variant_id)
      } else if(id %in% state[[role]]) {
        if(length(state[[opposite]])>=5L) { shiny::showNotification("The opposite direction already has five primers. Remove one first."); return() }
        state$custom <- .primer_mixer_bind(list(state$custom[state$custom$primer_id!=variant_id,,drop=FALSE],variant))
        state[[role]] <- setdiff(state[[role]],id)
        state[[opposite]] <- c(state[[opposite]],variant_id)
      } else {
        if(length(state[[role]])>=5L) { shiny::showNotification("This direction already has five primers. Remove one first."); return() }
        state[[role]] <- c(state[[role]],id)
      }
      sync()
    })
    map_data <- shiny::reactive({
      value <- result(); id <- selected_reference()
      if(is.null(id) || !nrow(value$pairs)) return(NULL)
      data <- primer_mixer_overview_data(value,id)
      data$domain <- domain(); data
    })
    output$result_status <- shiny::renderUI({
      value <- result()
      shiny::span(class="atlas-status-chip",nrow(value$pairs)," combinations × ",nrow(value$templates)," templates")
    })
    make_map <- function(width=NULL) {
      data <- map_data()
      if(is.null(data) || !nrow(data$rows)) return(shiny::div(class="callout",
        if(!length(forward()) || !length(reverse())) "Select at least one forward and one reverse primer, or load an example." else
          "No unambiguous inward product meets these settings on the coordinate template. Inspect template results for site evidence, or change the template or filters."))
      if(is.null(data$domain)) return(shiny::div(class="callout","This region has no verified bounds on the selected template."))
      annotation <- annotations(); annotation$landmark_id <- paste0("region_",seq_len(nrow(annotation)))
      primer_overview_map(data$rows,data$primers,data$marker,landmarks=annotation,
        focus=list(start=data$domain[1],end=data$domain[2]),domain=data$domain,citations=data$citations,
        active_pair_id=input$inspect_pair,pair_click_input=ns("pair_click"),viewport_width=width)
    }
    output$map <- shiny::renderUI(shiny::tagList(make_map(input$width %||% 650),
      if(!is.null(map_data()) && nrow(map_data()$rows)) shiny::p(class="atlas-context-line",
        "Product length colors: ",shiny::span(style="color:#2563eb","Illumina ≤500 bp")," · ",
        shiny::span(style="color:#7c3aed","PacBio ≤2,000 bp")," · ",
        shiny::span(style="color:#0f766e","Nanopore ≤3,500 bp")," · gray above 3,500 bp. Each product uses its shortest eligible tier; these are project planning ceilings.")))
    shiny::observeEvent(result()$pairs, {
      pairs <- result()$pairs; selected <- shiny::isolate(input$inspect_pair)
      if(length(selected)!=1L || !selected %in% pairs$pair_id) selected <- head(pairs$pair_id,1L)
      shiny::updateSelectInput(session,"inspect_pair",choices=stats::setNames(pairs$pair_id,pairs$pair_label),selected=selected)
    })
    shiny::observeEvent(input$pair_click, {
      click <- input$pair_click
      pair <- if(is.list(click)) click$pair_id %||% click$pair else click
      if(length(pair)==1L && pair %in% result()$pairs$pair_id)
        shiny::updateSelectInput(session,"inspect_pair",selected=pair)
    })
    output$pair_summaries <- shiny::renderUI({
      value <- result(); sums <- value$summary
      if(!nrow(sums)) return(NULL)
      ordered <- sums[order(-sums$n_length_eligible,-sums$n_exact_products,sums$pair_label),,drop=FALSE]
      shown <- head(ordered,4L)
      all_rows <- data.frame(Combination=sums$pair_label,Templates=sums$n_templates,
        `Exact products`=sums$n_exact_products,`With substitutions`=sums$n_mismatch_products,
        `Length eligible`=sums$n_length_eligible,`Ambiguous`=sums$n_ambiguous,
        `One site / short record`=sums$n_partial,`No eligible inward product`=sums$n_no_eligible_sites,
        `Unavailable`=sums$n_unavailable,check.names=FALSE)
      shiny::tagList(shiny::p(class="atlas-context-line","Sequence evidence across these templates. A missing site can reflect mismatches, truncation or an absent region; these counts are not PCR success rates."),
        shiny::div(class="explore-fit-grid",lapply(seq_len(nrow(shown)),function(i) {
          row <- shown[i,,drop=FALSE]
          shiny::div(class="explore-fit-card",
            shiny::strong(row$pair_label),
            shiny::span(class="explore-fit-count",row$n_exact_products," / ",row$n_templates," exact products"),
            shiny::div(class="explore-fit-track",shiny::span(style=paste0("width:",if(row$n_templates) 100*row$n_exact_products/row$n_templates else 0,"%"))),
            shiny::span(class="atlas-context-line",row$n_mismatch_products," with substitutions · ",row$n_length_eligible," within length filter"))
        })),
        shiny::tags$details(class="atlas-detail",shiny::tags$summary(paste("All",nrow(sums),"combination results")),
          display_table(all_rows)))
    })
    display_table <- function(rows) {
      shiny::div(class="explore-table-wrap",shiny::tags$table(class="table table-sm",
        shiny::tags$thead(shiny::tags$tr(lapply(names(rows),shiny::tags$th))),
        shiny::tags$tbody(lapply(seq_len(nrow(rows)),function(i) shiny::tags$tr(lapply(rows[i,,drop=FALSE],function(value)
          shiny::tags$td(if(is.na(value)) "—" else as.character(value))))))))
    }
    inspect_pair <- shiny::reactive({
      pairs <- result()$pairs
      if((input$inspect_pair %||% "") %in% pairs$pair_id) input$inspect_pair else head(pairs$pair_id,1L)
    })
    output$broader_evidence_ui <- shiny::renderUI({
      if(is.null(on_evidence) || !(active_marker()=="ITS_FUNGAL" ||
          (active_marker()=="COI" && target_group() %in% c("ALL","ARTHROPODS","METAZOA")))) return(NULL)
      shiny::tagList(shiny::actionButton(ns("broader_evidence"),
        if(active_marker()=="ITS_FUNGAL") "Inspect broader UNITE reference support" else "Open full arthropod COI reference fit",
        class="btn-outline-primary btn-sm"),shiny::p(class="atlas-context-line",
          "The broader reference collection has its own version and denominator. Exact oligo identities must be available; custom and modified primers do not inherit another assay's scores."))
    })
    output$fit_action <- shiny::renderUI({
      if(is.null(on_fit) || !nrow(result()$pairs)) return(NULL)
      shiny::actionButton(ns("open_fit"),"Check taxonomic fit for these combinations",class="btn-primary btn-sm")
    })
    shiny::observeEvent(input$open_fit, {
      value <- result()
      if(!nrow(value$pairs) || is.null(on_fit)) return()
      if(active_marker()=="ITS_FUNGAL") {
        pair <- value$pairs[value$pairs$pair_id %in% inspect_pair(),,drop=FALSE]
        if(nrow(pair)!=1L || is.null(on_evidence) || !isTRUE(on_evidence(pair,value$components))) return()
      }
      on_fit()
    })
    shiny::observeEvent(input$broader_evidence, {
      value <- result(); pair <- value$pairs[value$pairs$pair_id %in% inspect_pair(),,drop=FALSE]
      if(!is.null(on_evidence) && nrow(pair)==1L) on_evidence(pair,value$components)
    })
    output$template_results <- shiny::renderUI({
      value <- result(); pair <- inspect_pair()
      if(!length(pair)) return(shiny::p("Choose a combination first."))
      rows <- value$products[value$products$pair_id==pair,,drop=FALSE]
      if(!nrow(rows)) return(shiny::p("No templates are loaded. Add marker-matched FASTA to evaluate this candidate."))
      out <- data.frame(Template=rows$template_id,State=gsub("_"," ",rows$status),
        `Product bp`=rows$amplicon_bp,`Internal bp`=rows$internal_bp,
        `F substitutions`=rows$forward_mismatches,`R substitutions`=rows$reverse_mismatches,
        `F 3′ last 5`=rows$forward_terminal_5_mismatches,`R 3′ last 5`=rows$reverse_terminal_5_mismatches,
        `Length eligible`=rows$passes_length_filter,check.names=FALSE)
      display_table(out)
    })
    output$binding_trace <- shiny::renderUI({
      value <- result(); id <- selected_reference(); pair <- value$pairs[value$pairs$pair_id %in% inspect_pair(),,drop=FALSE]
      if(is.null(id) || nrow(pair)!=1L) return(NULL)
      components <- value$components[value$components$primer_id %in% c(pair$forward_id,pair$reverse_id),,drop=FALSE]
      shiny::tagList(shiny::h4("Primer bases on the coordinate template"),
        lapply(seq_len(nrow(components)),function(i) {
          component <- components[i,,drop=FALSE]
          site <- value$sites[value$sites$template_id==id & value$sites$primer_id==component$primer_id,,drop=FALSE]
          if(nrow(site)!=1L) return(NULL)
          shiny::div(class="explore-binding-trace",shiny::strong(component$primer_name," · ",gsub("_"," ",site$status)),
            shiny::div(class="atlas-context-line",if(!is.na(site$start)) paste0("Raw bases ",site$start,"–",site$end,"; both rows shown in primer 5′→3′ orientation.")),
            shiny::tags$pre(paste0("Primer    5′ ",component$sequence," 3′\n",
              "Template  5′ ",if(nzchar(site$binding_sequence_primer_oriented)) site$binding_sequence_primer_oriented else "No unique eligible full-length site"," 3′")),
            if(nzchar(site$mismatch_positions_5_to_3)) shiny::p(class="atlas-context-line","Substitution positions in primer orientation: ",site$mismatch_positions_5_to_3))
        }))
    })
    output$selected_sources <- shiny::renderUI({
      rows <- result()$components
      if(!nrow(rows)) return(shiny::p("Select primers to view their source sequences and citations."))
      shiny::tagList(shiny::p(class="atlas-context-line",result()$method),
        shiny::tags$ul(lapply(seq_len(nrow(rows)),function(i) shiny::tags$li(shiny::strong(rows$primer_name[i]),
          " · ",rows$source_sequence[i]," · ",rows$direction[i]," · ",
          if(nzchar(rows$source_url[i])) shiny::tags$a(href=rows$source_url[i],target="_blank",rel="noopener",rows$citation[i]) else rows$citation[i],
          shiny::div(class="atlas-context-line",rows$evidence_note[i]," ",rows$normalization_note[i])))))
    })
    output$downloads_ui <- shiny::renderUI({
      shiny::div(class="atlas-card-actions",
        shiny::downloadButton(ns("download_results"),"Template results CSV",class="btn-sm"),
        shiny::downloadButton(ns("download_primers"),"Primer sequences CSV",class="btn-sm"),
        shiny::downloadButton(ns("download_table"),"Paper table HTML",class="btn-sm"),
        shiny::downloadButton(ns("download_package"),"Evidence package ZIP",class="btn-sm"),
        if(!is.null(map_data()) && nrow(map_data()$rows)>0L && !is.null(domain()))
          shiny::downloadButton(ns("download_figure"),"Amplicon SVG",class="btn-sm"),
        if(nrow(templates())) shiny::downloadButton(ns("download_templates"),"Templates FASTA",class="btn-sm"))
    })
    output$download_results <- shiny::downloadHandler(
      filename=function() paste0("primer-atlas-",active_marker(),"-template-results.csv"),
      content=function(file) {
        utils::write.csv(primer_mixer_results_table(result()),file,row.names=FALSE,na="")
      })
    output$download_primers <- shiny::downloadHandler(filename=function() paste0("primer-atlas-",active_marker(),"-primers.csv"),
      content=function(file) utils::write.csv(result()$components,file,row.names=FALSE,na=""))
    output$download_templates <- shiny::downloadHandler(filename=function() paste0("primer-atlas-",active_marker(),"-templates.fasta"),
      content=function(file) {
        rows <- templates()
        writeLines(unlist(lapply(seq_len(nrow(rows)),function(i) c(paste0(">",rows$template_id[i]," | ",rows$version[i]),rows$sequence[i]))),file)
      })
    output$download_figure <- shiny::downloadHandler(filename=function() paste0("primer-atlas-",active_marker(),"-amplicons.svg"),content=function(file) {
      writeLines(primer_mixer_export_svg(make_map(NULL),result(),selected_reference(),domain()),file)
    })
    current_taxonomic_report <- function() {
      if(is.null(taxonomic_report)) return(NULL)
      report <- taxonomic_report()
      if(is.null(report)) return(NULL)
      value <- result()
      if(!identical(as.character(report$metadata$marker_id),as.character(value$marker_id)) ||
          !setequal(report$templates$template_id,value$templates$template_id))
        stop("The taxonomic report does not belong to the current marker/template cohort.")
      report
    }
    write_template_fasta <- function(file) {
      rows <- templates()
      sequence <- if("source_sequence" %in% names(rows)) rows$source_sequence else rows$sequence
      writeLines(unlist(lapply(seq_len(nrow(rows)),function(i) c(paste0(">",rows$template_id[i]," | ",rows$version[i]),sequence[i]))),file)
    }
    output$download_table <- shiny::downloadHandler(filename=function() paste0("primer-atlas-",active_marker(),"-paper-table.html"),content=function(file)
      writeLines(primer_mixer_export_html(result(),current_taxonomic_report()),file))
    output$download_package <- shiny::downloadHandler(filename=function() paste0("primer-atlas-",active_marker(),"-evidence.zip"),content=function(file) {
      value <- result(); report <- current_taxonomic_report()
      directory <- tempfile("primer-atlas-evidence-"); dir.create(directory)
      on.exit(unlink(directory,recursive=TRUE),add=TRUE)
      utils::write.csv(primer_mixer_results_table(value),file.path(directory,"template-results.csv"),row.names=FALSE,na="")
      utils::write.csv(value$components,file.path(directory,"primer-components.csv"),row.names=FALSE,na="")
      utils::write.csv(value$sites,file.path(directory,"primer-sites.csv"),row.names=FALSE,na="")
      utils::write.csv(value$templates[,setdiff(names(value$templates),c("sequence","source_sequence")),drop=FALSE],
        file.path(directory,"source-template-metadata.csv"),row.names=FALSE,na="")
      write_template_fasta(file.path(directory,"templates.fasta"))
      writeLines(primer_mixer_export_html(value,report),file.path(directory,"paper-table.html"))
      if(!is.null(map_data()) && nrow(map_data()$rows)>0L && !is.null(domain()))
        writeLines(primer_mixer_export_svg(make_map(NULL),value,selected_reference(),domain()),file.path(directory,"amplicons.svg"))
      if(!is.null(report)) {
        utils::write.csv(report$summary,file.path(directory,"taxonomic-fit.csv"),row.names=FALSE,na="")
        utils::write.csv(report$details,file.path(directory,"taxonomic-template-evidence.csv"),row.names=FALSE,na="")
      }
      metadata <- list(app_version="2026.10.31",marker_id=value$marker_id,method=value$method,settings=value$settings,
        coordinate_reference=selected_reference(),display_domain=domain(),components=value$components,
        requested_combinations=value$pairs,n_templates=nrow(value$templates),
        taxonomy=if(!is.null(report)) report$metadata else NULL)
      jsonlite::write_json(metadata,file.path(directory,"methods.json"),pretty=TRUE,auto_unbox=TRUE,na="null",null="null")
      writeLines(c("Primer Atlas: session candidate evidence",value$method,.primer_mixer_export_scope,
        "Inclusive encoded product bp = reverse site end - forward site start + 1.",
        "Internal bp = reverse site start - forward site end - 1.",
        "Source-template metadata retains original taxonomy. Taxonomic evidence additionally retains any explicitly user-provided session annotation.",
        "Predicted and partial source records remain labelled. Encoded sequence spans do not establish genomic PCR products or experimental amplification.",
        "Every requested combination is retained. Length/platform filters change displayed eligibility without changing template denominators."),
        file.path(directory,"METHODS.txt"))
      artifacts <- list.files(directory,full.names=TRUE)
      hashes <- vapply(artifacts,digest::digest,character(1),algo="sha256",file=TRUE)
      writeLines(paste(hashes,basename(artifacts)),file.path(directory,"SHA256SUMS"))
      zip::zipr(file,list.files(directory,full.names=TRUE),root=directory,include_directories=FALSE)
    })
    list(catalog=catalog,templates=templates,forward=forward,reverse=reverse,result=result,
      selected_reference=selected_reference,site_data=site_data,domain=domain,state=state,load_example=load_example)
  })
}
