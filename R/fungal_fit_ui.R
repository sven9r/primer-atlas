# A focused reference-support view; no COI mismatch penalties are reused here.
fungal_fit_bar <- function(value, total, label = "Exact product support") {
  fraction <- if (total > 0) value / total else 0
  div(class = "atlas-evidence-bar",
      div(class = "atlas-evidence-bar__label", span(label),
          strong(class = "atlas-evidence-bar__value", paste(value, "of", total))),
      div(class = "atlas-evidence-bar__track", role = "img",
          `aria-label` = paste(label, value, "of", total, "retained references"),
          span(class = "atlas-evidence-bar__fill", style = paste0("width:", round(100 * fraction, 2), "%;"))))
}

fungalFitUI <- function(id) {
  ns <- shiny::NS(id)
  index <- fungal_read_index()
  versions <- setNames(index$id, sub(" ·.*$", "", index$label))
  div(class = "atlas-workspace fungal-fit-workspace",
      tags$script(HTML("(function(){ var pending=null; document.addEventListener('click', function(e) { var b=e.target.closest('.atlas-fit-taxon'); if(b && window.Shiny) { pending=b.dataset.inspector; Shiny.setInputValue(b.dataset.input, {group:b.dataset.group, nonce:Date.now()}, {priority:'event'}); } }); new MutationObserver(function(){ if(!pending) return; var p=document.getElementById(pending); if(p && p.querySelector('.atlas-fit-inspector')) { pending=null; window.requestAnimationFrame(function(){ p.scrollIntoView({behavior:'smooth',block:'start'}); }); } }).observe(document.documentElement,{childList:true,subtree:true}); })();")),
      div(class = "atlas-workspace-toolbar atlas-fit-toolbar",
          selectInput(ns("version"), "Reference version", versions, selected = tail(index$id, 1L)),
          selectizeInput(ns("pair"), "Primer combination", choices = NULL,
              options = list(placeholder = "Search primer combinations")),
          selectInput(ns("rank"), "Taxonomic rank", c("Classes" = "class", "Phyla" = "phylum"), "class")),
      div(class = "atlas-fit-context", uiOutput(ns("context")),
          actionButton(ns("open_explorer"), "Show primer product →", class = "btn-outline-primary btn-sm")),
      uiOutput(ns("metrics")),
      div(id = ns("taxon_groups"), class = "atlas-panel atlas-fit-results",
          div(class = "atlas-fit-results-heading", h3("Reference support by taxon"),
              p(class = "atlas-context-line", "Click a taxon to inspect its references.")),
          div(class = "atlas-fit-list-tools",
              textInput(ns("search"), "Find a taxon", placeholder = "Class, family or parent lineage"),
              selectInput(ns("platform"), "Encoded product ceiling",
                  c("All products" = "ALL", "Illumina · ≤500 bp" = "ILLUMINA",
                    "PacBio · ≤2,000 bp" = "PACBIO", "Nanopore · ≤3,500 bp" = "NANOPORE"), "ALL")),
          uiOutput(ns("taxa")),
          div(class = "atlas-fit-pagination", actionButton(ns("previous"), "Previous", class = "btn-outline-secondary btn-sm"),
              uiOutput(ns("page_status"), inline = TRUE),
              actionButton(ns("next"), "Next", class = "btn-outline-secondary btn-sm"))),
      uiOutput(ns("taxon_inspector")),
      tags$details(class = "atlas-detail",
          tags$summary("How fungal reference support is measured"),
          p("A supported product has two unique exact IUPAC-compatible binding sites in an inward orientation on the retained sequence. A missing exact site is a search result on that record; it does not establish PCR failure or absence of the locus."),
          p("Bars use all retained templates in the displayed taxon as their denominator. Partial, ambiguous, withheld and unavailable evidence remain separate. The ceiling selector uses encoded reference spans including primer sites; ambiguity and possible assembly gaps remain visible in sequence evidence."),
          p("This bounded reference panel is not a representative sample of all fungi. Mismatch penalties and experimental amplification probabilities have not been measured in this view.")),
      tags$details(class = "atlas-detail",
          tags$summary("Download this reference-support view"),
          div(class = "atlas-card-actions",
              downloadButton(ns("download_summary"), "Taxon counts (CSV)", class = "btn-outline-secondary btn-sm"),
              downloadButton(ns("download_details"), "Reference evidence (CSV)", class = "btn-outline-secondary btn-sm"))))
}

fungalFitServer <- function(id, selection = NULL) {
  moduleServer(id, function(input, output, session) {
    index <- fungal_read_index()
    latest <- tail(index$id, 1L)
    requested <- reactiveVal(NULL)
    chosen_group <- reactiveVal(NULL)
    page <- reactiveVal(1L)
    snapshot <- reactive({
      version <- input$version %||% latest
      validate(need(version %in% index$id, "Choose an installed fungal reference version."))
      fungal_read_snapshot(version)
    })
    pair_catalog <- reactive(snapshot()$pairs)
    selected_pair <- reactive({
      rows <- pair_catalog()
      target <- requested()
      if (!is.null(target) && identical(target$version, snapshot()$manifest$id) && target$pair %in% rows$pair_id) return(target$pair)
      if (!is.null(input$pair) && input$pair %in% rows$pair_id) return(input$pair)
      preferred <- match("ITS1 + ITS4", rows$pair_label)
      if (!is.na(preferred)) rows$pair_id[preferred] else rows$pair_id[1L]
    })
    if (!is.null(selection)) observeEvent(selection(), {
      value <- selection()
      req(is.list(value), value$version %in% index$id)
      frozen <- fungal_read_snapshot(value$version)
      req(length(value$pair) == 1L, value$pair %in% frozen$pairs$pair_id)
      requested(value)
      chosen_group(NULL)
      page(1L)
      updateSelectInput(session, "version", selected = value$version)
      updateTextInput(session, "search", value = "")
    }, ignoreNULL = TRUE)
    observeEvent(snapshot(), {
      rows <- pair_catalog()
      order <- order(rows$evidence_status != "documented_pair", rows$pair_label)
      rows <- rows[order, , drop = FALSE]
      updateSelectizeInput(session, "pair", choices = setNames(rows$pair_id, rows$pair_label),
          selected = isolate(selected_pair()), server = TRUE)
      available <- fungal_fit_available_ranks(snapshot())
      labels <- c(phylum = "Phyla", class = "Classes", order = "Orders", family = "Families")
      rank <- isolate(input$rank)
      if (is.null(rank) || !rank %in% available) rank <- if ("class" %in% available) "class" else available[1L]
      updateSelectInput(session, "rank", choices = setNames(available, labels[available]), selected = rank)
      chosen_group(NULL)
      page(1L)
    }, ignoreInit = FALSE)
    observeEvent(input$pair, {
      value <- requested()
      if (!is.null(value) && identical(input$version, value$version) && identical(input$pair, value$pair)) requested(NULL)
      chosen_group(NULL)
      page(1L)
    }, ignoreNULL = TRUE)
    fit_data <- reactive({
      available <- fungal_fit_available_ranks(snapshot())
      rank <- input$rank %||% "class"
      if (!rank %in% available) rank <- available[1L]
      fungal_fit_summary(snapshot(), rank = rank, pair_ids = selected_pair(), platform = input$platform %||% "ALL")
    })
    groups <- reactive({
      rows <- fit_data()$summary
      search <- trimws(input$search %||% "")
      if (nzchar(search)) rows <- rows[grepl(tolower(search), tolower(paste(rows$taxon, rows$parent_path)), fixed = TRUE), , drop = FALSE]
      rows[order(-rows$n_panel, rows$parent_path, rows$taxon), , drop = FALSE]
    })
    observeEvent(list(input$rank, input$search), {
      chosen_group(NULL)
      page(1L)
    }, ignoreInit = TRUE)
    observeEvent(input$platform, page(1L), ignoreInit = TRUE)
    pages <- reactive(max(1L, ceiling(nrow(groups()) / 8L)))
    current_page <- reactive(min(page(), pages()))
    observeEvent(input$previous, page(max(1L, current_page() - 1L)))
    observeEvent(input[["next"]], page(min(pages(), current_page() + 1L)))
    observeEvent(input$taxon_click, {
      value <- input$taxon_click$group
      if (length(value) == 1L && value %in% groups()$group_id) chosen_group(value)
    })
    selected_taxon <- reactive({
      value <- chosen_group()
      if (is.null(value) || !value %in% groups()$group_id) return(NULL)
      value
    })
    taxon_details <- reactive({
      rows <- fit_data()$details
      value <- selected_taxon()
      if (is.null(value)) return(rows[0, , drop = FALSE])
      rows[rows$group_id == value, , drop = FALSE]
    })
    output$context <- renderUI({
      data <- fit_data()
      pair <- pair_catalog()[match(selected_pair(), pair_catalog()$pair_id), , drop = FALSE]
      div(class = "atlas-context-line",
          span(class = "atlas-status-chip", if (pair$evidence_status == "documented_pair") "Documented pair" else "Candidate combination"),
          span(paste(data$metadata$n_snapshot_templates, "retained templates ·", sub(" ·.*$", "", snapshot()$manifest$label))),
          span("Bounded reference panel"))
    })
    output$metrics <- renderUI({
      rows <- fit_data()$summary
      total <- sum(rows$n_panel)
      supported <- sum(rows$n_exact_products)
      partial <- sum(rows$n_partial_sites)
      div(class = "atlas-fit-metrics",
          div(class = "atlas-metric", span("Exact product support"), strong(paste(supported, "of", total)),
              span("retained reference templates")),
          div(class = "atlas-metric", span("One exact site"), strong(partial), span("partial pair evidence")),
          div(class = "atlas-metric", span("Taxonomy groups"), strong(nrow(rows)), span("includes unresolved assignments")))
    })
    output$taxa <- renderUI({
      rows <- groups()
      if (!nrow(rows)) return(div(class = "atlas-context-line", "No taxa match this search. Clear the search to return to the full panel."))
      rows <- rows[seq.int((current_page() - 1L) * 8L + 1L, min(current_page() * 8L, nrow(rows))), , drop = FALSE]
      div(class = "atlas-action-grid atlas-fit-taxon-grid", lapply(seq_len(nrow(rows)), function(i) {
        row <- rows[i, , drop = FALSE]
        other <- row$n_panel - row$n_exact_products - row$n_partial_sites - row$n_unavailable
        tags$button(type = "button", class = paste("atlas-taxon-card atlas-fit-taxon", if (identical(selected_taxon(), row$group_id)) "is-selected"),
            `data-input` = session$ns("taxon_click"), `data-group` = row$group_id,
            `data-inspector` = session$ns("taxon_inspector"),
            `aria-pressed` = if (identical(selected_taxon(), row$group_id)) "true" else "false",
            div(class = "atlas-taxon-card__heading", strong(row$taxon), span(class = "atlas-status-chip", paste(row$n_panel, "refs"))),
            if (nzchar(row$parent_path)) span(class = "atlas-taxon-card__meta", row$parent_path),
            fungal_fit_bar(row$n_exact_products, row$n_panel),
            span(class = "atlas-taxon-card__meta", paste(row$n_partial_sites, "partial ·", other, "other site states ·", row$n_unavailable, "unavailable")),
            if (row$platform != "ALL") span(class = "atlas-taxon-card__meta", paste(row$n_length_eligible, "of", row$n_exact_products, "exact products within the encoded bp ceiling")),
            span(class = "atlas-taxon-card__action", "Inspect references →"))
      }))
    })
    output$page_status <- renderUI(span(class = "atlas-context-line", paste("Page", current_page(), "of", pages())))
    output$taxon_inspector <- renderUI({
      value <- selected_taxon()
      if (is.null(value)) return(NULL)
      row <- groups()[match(value, groups()$group_id), , drop = FALSE]
      div(class = "atlas-panel atlas-fit-inspector",
          div(class = "atlas-fit-results-heading", h3(row$taxon),
              p(class = "atlas-context-line", row$parent_path)),
          div(class = "atlas-card-actions",
              tags$a(href = paste0("#", session$ns("taxon_groups")), class = "btn btn-outline-secondary btn-sm", "Back to taxon groups ↑"),
              downloadButton(session$ns("download_taxon"), "These references (CSV)", class = "btn-outline-secondary btn-sm")),
          p(class = "atlas-context-line", paste(row$n_exact_products, "exact products ·", row$n_partial_sites, "partial site pairs ·",
              row$n_unresolved_sites, "ambiguous/unresolved ·", row$n_no_exact_sites, "no exact sites ·",
              row$n_product_withheld, "withheld geometry ·", row$n_unavailable, "unavailable")),
          DTOutput(session$ns("references")))
    })
    output$references <- renderDT({
      rows <- taxon_details()
      state_labels <- c(exact_product = "Exact product", partial_sites = "One exact site", unresolved_sites = "Ambiguous / unresolved sites",
          no_exact_sites = "No exact sites found", product_withheld = "Product withheld", evidence_unavailable = "Evidence unavailable")
      display <- data.frame(Accession = rows$accession, Organism = rows$organism,
          `Site evidence` = unname(state_labels[rows$fit_state]), `Encoded product (bp)` = rows$amplicon_bp,
          `Ambiguous product bases` = rows$n_ambiguous_product_bases, Reason = rows$reason, check.names = FALSE)
      datatable(display, rownames = FALSE, escape = TRUE, options = list(pageLength = 5, lengthChange = FALSE, scrollX = TRUE), class = "compact stripe")
    })
    output$download_summary <- downloadHandler(filename = function() paste0(snapshot()$manifest$id, "-fungal-fit-", fit_data()$metadata$rank, ".csv"),
        content = function(file) write.csv(groups(), file, row.names = FALSE, na = ""))
    output$download_details <- downloadHandler(filename = function() paste0(snapshot()$manifest$id, "-fungal-fit-reference-evidence.csv"),
        content = function(file) {
          rows <- fit_data()$details
          rows <- rows[rows$group_id %in% groups()$group_id, , drop = FALSE]
          write.csv(rows, file, row.names = FALSE, na = "")
        })
    output$download_taxon <- downloadHandler(filename = function() paste0(snapshot()$manifest$id, "-selected-fungal-taxon.csv"),
        content = function(file) write.csv(taxon_details(), file, row.names = FALSE, na = ""))
    list(snapshot = snapshot, selected_pair = selected_pair, fit_data = fit_data,
        selected_taxon = selected_taxon, taxon_details = taxon_details, groups = groups)
  })
}
