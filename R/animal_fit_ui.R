# Shared target/marker fit view. Its only computational input is the current
# explorer's actual candidate/template result, never a different cohort's scores.
animalFitUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::tags$style(shiny::HTML(".animal-fit-toolbar{display:grid;grid-template-columns:minmax(130px,1fr) minmax(180px,2fr);gap:14px}.animal-fit-pairs{margin-top:8px}.animal-fit-taxa{display:grid;gap:14px}.animal-fit-taxon-card{padding:16px;border:1px solid #dce5ed;border-radius:12px;background:#fff}.animal-fit-taxon-heading{display:flex;justify-content:space-between;gap:12px}.animal-fit-taxon-heading h4{margin:0;font-size:17px}.animal-fit-lineage{font-size:12px;color:#66788b;margin:5px 0 12px;overflow-wrap:anywhere}.animal-fit-cell-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(155px,1fr));gap:10px}.animal-fit-cell{appearance:none;text-align:left;background:#f8fafc;border:1px solid #e0e8f0;border-radius:9px;padding:11px;color:#263d52;cursor:pointer;min-width:0}.animal-fit-cell:hover,.animal-fit-cell:focus-visible{border-color:#087e75;background:#eef8f5}.animal-fit-cell strong{display:block;font-size:12px;overflow-wrap:anywhere}.animal-fit-cell-count{display:block;font-size:17px;font-weight:650;margin:7px 0}.animal-fit-cell-note{display:block;font-size:11px;color:#66788b;line-height:1.45}.animal-fit-stack{display:flex;height:7px;background:#e2e8f0;border-radius:5px;overflow:hidden;margin-bottom:8px}.animal-fit-pagination{display:flex;gap:12px;align-items:center;justify-content:flex-end;margin-top:14px}.animal-fit-legend{display:flex;flex-wrap:wrap;gap:10px;font-size:11px;color:#64768b;margin:12px 0}.animal-fit-legend i{display:inline-block;width:8px;height:8px;border-radius:2px;margin-right:4px}.animal-fit-table-wrap{overflow-x:auto}.animal-fit-table-wrap table{font-size:12px}.animal-fit-context{padding:12px 0;color:#5d7083}.animal-fit-inspector{margin-top:18px}.animal-fit-provenance{font-size:12px;overflow-wrap:anywhere}@media(max-width:500px){.animal-fit-toolbar{grid-template-columns:1fr}.animal-fit-cell-grid{grid-template-columns:1fr 1fr}.animal-fit-taxon-card{padding:12px}.animal-fit-cell{padding:9px}}")),
    shiny::div(class = "animal-fit-workspace atlas-workspace", id = ns("workspace"),
      shiny::uiOutput(ns("context")),
      shiny::div(class = "animal-fit-toolbar",
        shiny::selectInput(ns("rank"), "Taxonomic rank",
          c("Classes" = "class", "Orders" = "order", "Families" = "family", "Genera" = "genus", "Species" = "species", "Phyla" = "phylum"), "class"),
        shiny::textInput(ns("search"), "Find a taxon", placeholder = "Taxon name or parent lineage")),
      shiny::div(class = "animal-fit-pairs", shiny::selectizeInput(ns("pairs"), "Compare combinations",
        choices = NULL, multiple = TRUE, options = list(maxItems = 4, plugins = list("remove_button"),
          placeholder = "Choose up to four selected combinations"))),
      shiny::uiOutput(ns("matrix")),
      shiny::div(class = "animal-fit-pagination",
        shiny::actionButton(ns("previous"), "Previous", class = "btn-outline-secondary btn-sm"),
        shiny::uiOutput(ns("page_status"), inline = TRUE),
        shiny::actionButton(ns("next"), "Next", class = "btn-outline-secondary btn-sm")),
      shiny::tags$details(class = "atlas-detail animal-fit-inspector", id = ns("inspector"),
        shiny::tags$summary("Inspect templates and sequence evidence"), shiny::uiOutput(ns("inspector_content"))),
      shiny::tags$details(class = "atlas-detail",
        shiny::tags$summary("Sources, denominators and taxonomy metadata"),
        shiny::uiOutput(ns("sources")),
        shiny::p(class = "atlas-context-line", "Every cell uses all selected templates in its full parent lineage as the denominator. Exact, substitution, partial, ambiguous and no-product/unavailable states partition that set. Both-site counts are reported separately. These sequence searches do not predict PCR success or establish coverage of an entire taxon."),
        shiny::tags$details(class = "atlas-detail",
          shiny::tags$summary("Assign taxonomy to your templates"),
          shiny::p(class = "atlas-context-line", "Headers never assign taxonomy. Paste CSV with exact template_id values and any of kingdom, phylum, class, order, family, genus or species. Each supplied row replaces that template's lineage metadata; omitted ranks stay Unresolved. Sequence sources remain unchanged. Your annotation is labelled user-provided metadata."),
          shiny::textAreaInput(ns("taxonomy_csv"), "Session taxonomy CSV", rows = 5,
            placeholder = "template_id,phylum,class,order,family,genus,species\nmy_template,Chordata,Actinopterygii,,,Danio,Danio rerio"),
          shiny::div(class = "atlas-card-actions",
            shiny::actionButton(ns("apply_taxonomy"), "Apply taxonomy metadata", class = "btn-outline-primary btn-sm"),
            shiny::actionButton(ns("clear_taxonomy"), "Restore source taxonomy", class = "btn-outline-secondary btn-sm")),
          shiny::uiOutput(ns("taxonomy_status"))),
        shiny::div(class = "atlas-card-actions",
          shiny::downloadButton(ns("download_summary"), "Taxon counts CSV", class = "btn-sm"),
          shiny::downloadButton(ns("download_details"), "Template evidence CSV", class = "btn-sm"),
          shiny::downloadButton(ns("download_taxonomy"), "Current taxonomy CSV", class = "btn-sm")))),
    shiny::tags$script(shiny::HTML(paste0("(function(){var id=", jsonlite::toJSON(ns("workspace"), auto_unbox = TRUE),
      ";document.addEventListener('click',function(e){var root=document.getElementById(id);var cell=e.target.closest('.animal-fit-cell');if(!root||!cell||!root.contains(cell)||!window.Shiny)return;Shiny.setInputValue(",
      jsonlite::toJSON(ns("cell_click"), auto_unbox = TRUE),
      ",{group:cell.dataset.taxonKey,pair:cell.dataset.pairId,nonce:Date.now()},{priority:'event'});var detail=document.getElementById(",
      jsonlite::toJSON(ns("inspector"), auto_unbox = TRUE), ");if(detail)detail.open=true;});})();")))
  )
}

animalFitServer <- function(id, marker, target_group, result, root = ".", taxonomy = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    page <- shiny::reactiveVal(1L)
    selected_cell <- shiny::reactiveVal(NULL)
    user_taxonomy <- shiny::reactiveVal(NULL)
    user_note <- shiny::reactiveVal("")
    value <- shiny::reactive({
      own <- result()
      shiny::req(is.list(own), !is.null(own$pairs), !is.null(own$templates), !is.null(own$products))
      own
    })
    cohort <- shiny::reactive({
      own <- value(); templates <- own$templates
      digest::digest(list(marker = marker(), target = target_group(), ids = templates$template_id,
        sequence = templates$sequence, versions = templates$version), algo = "sha256")
    })
    shiny::observeEvent(cohort(), {
      user_taxonomy(NULL); user_note(""); selected_cell(NULL); page(1L)
    }, priority = 100)
    selected_rank <- shiny::reactive({
      rank <- input$rank
      if (length(rank) != 1L || !rank %in% .animal_fit_ranks) "class" else rank
    })
    data <- shiny::reactive(animal_fit_data(value(), selected_rank(),
      taxonomy = if (is.function(taxonomy)) taxonomy() else taxonomy,
      user_taxonomy = user_taxonomy()))
    default_pairs <- shiny::reactive({
      own <- value(); ids <- own$pairs$pair_id
      summary <- own$summary
      if (is.data.frame(summary) && nrow(summary) && "pair_id" %in% names(summary)) {
        index <- match(ids, summary$pair_id)
        count <- function(field) {
          x <- if (field %in% names(summary)) summary[[field]][index] else rep(0, length(ids))
          x[is.na(x)] <- 0; x
        }
        ids <- ids[order(-count("n_length_eligible"), -count("n_exact_products"), seq_along(ids))]
      }
      utils::head(ids, 4L)
    })
    shiny::observeEvent(value()$pairs, {
      rows <- value()$pairs; selected <- shiny::isolate(input$pairs)
      selected <- intersect(selected, rows$pair_id)
      if (!length(selected)) selected <- default_pairs()
      shiny::updateSelectizeInput(session, "pairs", choices = stats::setNames(rows$pair_id, rows$pair_label),
        selected = selected, server = FALSE)
    }, ignoreInit = FALSE)
    selected_pairs <- shiny::reactive({
      ids <- value()$pairs$pair_id
      own <- input$pairs
      if (is.null(own)) return(default_pairs())
      utils::head(intersect(own, ids), 4L)
    })
    scoped <- shiny::reactive({
      rows <- data()$summary
      if (!nrow(rows)) return(rows)
      rows <- rows[rows$pair_id %in% selected_pairs(), , drop = FALSE]
      search <- trimws(.animal_fit_or(input$search, ""))
      if (nzchar(search)) rows <- rows[grepl(tolower(search), tolower(paste(rows$taxon, rows$parent_path)), fixed = TRUE), , drop = FALSE]
      rows
    })
    groups <- shiny::reactive({
      rows <- scoped()
      if (!nrow(rows)) return(rows)
      rows[!duplicated(rows$taxon_key), , drop = FALSE]
    })
    pages <- shiny::reactive(max(1L, ceiling(nrow(groups()) / 8L)))
    current_page <- shiny::reactive(min(page(), pages()))
    shiny::observeEvent(list(selected_rank(), input$search, selected_pairs()), page(1L))
    shiny::observeEvent(input$previous, page(max(1L, current_page() - 1L)))
    shiny::observeEvent(input[["next"]], page(min(pages(), current_page() + 1L)))
    shiny::observeEvent(input$cell_click, {
      click <- input$cell_click; rows <- scoped()
      if (is.list(click) && length(click$group) == 1L && length(click$pair) == 1L &&
          nrow(rows) && any(rows$taxon_key == click$group & rows$pair_id == click$pair))
        selected_cell(list(group = click$group, pair = click$pair))
    })
    selected_details <- shiny::reactive({
      rows <- data()$details; chosen <- selected_cell()
      if (is.null(chosen) || !nrow(rows)) return(rows[FALSE, , drop = FALSE])
      rows[rows$taxon_key == chosen$group & rows$pair_id == chosen$pair, , drop = FALSE]
    })
    output$context <- shiny::renderUI({
      own <- data()
      templates <- own$templates
      n_partial <- if("partial" %in% names(templates)) sum(templates$partial %in% TRUE) else 0L
      n_predicted <- if("source_status" %in% names(templates)) sum(grepl("model|predicted",templates$source_status,ignore.case=TRUE),na.rm=TRUE) else 0L
      shiny::div(class = "animal-fit-context", shiny::strong(paste(marker(), "·", target_group())),
        " · ", own$metadata$n_templates, " selected templates · ", own$metadata$n_pairs, " combinations",
        shiny::div(class = "atlas-status-chip", if (own$metadata$n_templates <= 20L) "Small selected panel" else "Selected sequence panel",
          " · counts describe these records, not the whole taxon"),
        if(n_partial || n_predicted) shiny::div(class="atlas-context-line",n_partial," partial source records · ",n_predicted," predicted sequence models. These source flags are separate from primer-site completeness."),
        shiny::div(class = "atlas-context-line", "Fit is calculated from the primers and templates selected in Explore. Reference cohorts are never substituted when a marker or target changes."))
    })
    colour <- c(exact = "#087e75", substitution = "#7257c9", partial = "#d99a31", ambiguous = "#bd7194", missing = "#c3cdd8")
    output$matrix <- shiny::renderUI({
      own <- data()
      if (!own$metadata$n_pairs) return(shiny::div(class = "callout", "Choose forward and reverse primers in Explore, then return here to inspect taxonomic fit."))
      if (!own$metadata$n_templates) return(shiny::div(class = "callout", "Load marker-matched reference or FASTA templates in Explore first."))
      if (!length(selected_pairs())) return(shiny::div(class = "callout", "Choose a combination above."))
      taxa <- groups()
      if (!nrow(taxa)) return(shiny::div(class = "callout", "No taxon matches this search."))
      indices <- seq.int((current_page() - 1L) * 8L + 1L, min(current_page() * 8L, nrow(taxa)))
      rows <- scoped()
      shiny::tagList(shiny::div(class = "animal-fit-legend",
        lapply(names(colour), function(state) shiny::span(shiny::tags$i(style = paste0("background:", colour[[state]])),
          switch(state, exact = "Exact product", substitution = "With substitutions", partial = "Partial sites", ambiguous = "Ambiguous", missing = "No product / unavailable")))),
        shiny::div(class = "animal-fit-taxa", lapply(indices, function(i) {
          taxon <- taxa[i, , drop = FALSE]; cells <- rows[rows$taxon_key == taxon$taxon_key, , drop = FALSE]
          cells <- cells[match(selected_pairs(), cells$pair_id, nomatch = 0L), , drop = FALSE]
          shiny::div(class = "animal-fit-taxon-card", `data-taxon-key` = taxon$taxon_key,
            shiny::div(class = "animal-fit-taxon-heading", shiny::h4(taxon$taxon), shiny::span(paste0(taxon$n_panel,if(taxon$n_panel==1L) " template" else " templates"))),
            shiny::div(class = "animal-fit-lineage", taxon$parent_path),
            shiny::div(class = "animal-fit-cell-grid", lapply(seq_len(nrow(cells)), function(j) {
              row <- cells[j, , drop = FALSE]
              counts <- c(exact = row$n_exact, substitution = row$n_substitution, partial = row$n_partial,
                ambiguous = row$n_ambiguous, missing = row$n_missing)
              shiny::tags$button(type = "button", class = "animal-fit-cell",
                `data-taxon-key` = row$taxon_key, `data-pair-id` = row$pair_id,
                `aria-label` = paste(row$taxon, row$pair_label, row$n_exact, "exact products of", row$n_panel, "templates; inspect sequence evidence"),
                shiny::strong(row$pair_label), shiny::span(class = "animal-fit-cell-count", row$n_exact, " / ", row$n_panel, " exact"),
                shiny::div(class = "animal-fit-stack", lapply(names(counts), function(state) shiny::span(
                  style = paste0("width:", 100 * counts[[state]] / row$n_panel, "%;background:", colour[[state]]),
                  title = paste(state, counts[[state]], "of", row$n_panel)))),
                shiny::span(class = "animal-fit-cell-note", row$n_substitution, " with substitutions · ", row$n_scorable, " with both unique sites"),
                shiny::span(class = "animal-fit-cell-note", row$n_partial, " partial · ", row$n_ambiguous, " ambiguous · ", row$n_missing, " no product / unavailable"))
            })))
        })))
    })
    output$page_status <- shiny::renderUI(shiny::span(class = "atlas-context-line", "Page ", current_page(), " / ", pages(), " · ", nrow(groups()), " taxa"))
    table <- function(rows) shiny::div(class = "animal-fit-table-wrap", shiny::tags$table(class = "table table-sm",
      shiny::tags$thead(shiny::tags$tr(lapply(names(rows), shiny::tags$th))),
      shiny::tags$tbody(lapply(seq_len(nrow(rows)), function(i) shiny::tags$tr(lapply(rows[i, , drop = FALSE], function(x)
        shiny::tags$td(if (is.na(x)) "—" else as.character(x))))))))
    output$inspector_content <- shiny::renderUI({
      rows <- selected_details()
      if (!nrow(rows)) return(shiny::p(class = "atlas-context-line", "Click a taxon–combination cell to inspect every template in its denominator."))
      field <- function(name) if (name %in% names(rows)) rows[[name]] else rep(NA, nrow(rows))
      source_flags <- vapply(seq_len(nrow(rows)),function(i) {
        flags <- c(if(isTRUE(field("partial")[i])) "partial source",
          if(!is.na(field("source_status")[i]) && grepl("model|predicted",field("source_status")[i],ignore.case=TRUE)) "predicted model")
        if(length(flags)) paste(flags,collapse="; ") else "Not flagged"
      },character(1))
      out <- data.frame(Template = rows$template_id, Organism = field("organism"), State = gsub("_", " ", rows$status),
        `Source flags`=source_flags,
        `F substitutions` = field("forward_mismatches"), `R substitutions` = field("reverse_mismatches"),
        `Product bp` = field("amplicon_bp"), `Both sites` = rows$two_unique_sites,
        `Length eligible` = rows$length_eligible, `Sequence version` = field("version"), check.names = FALSE)
      shiny::tagList(shiny::h4(rows$taxon[1L], " · ", rows$pair_label[1L]),
        shiny::p(class = "atlas-context-line", nrow(rows), " templates in the full denominator · ", sum(rows$two_unique_sites), " with two unique eligible sites. Geometry-withheld results remain explicit."),
        table(out), shiny::tags$details(class = "atlas-detail animal-fit-provenance", shiny::tags$summary("Exact sequence and taxonomy sources"),
          shiny::tags$ul(lapply(seq_len(nrow(rows)), function(i) shiny::tags$li(shiny::strong(rows$template_id[i]),
            " · sequence: ", if (nzchar(.animal_fit_text(field("source_url")[i]))) shiny::tags$a(href = field("source_url")[i], target = "_blank", rel = "noopener", field("version")[i]) else field("version")[i],
            " · taxonomy: ", if (nzchar(.animal_fit_text(rows$taxonomy_source_url[i]))) shiny::tags$a(href = rows$taxonomy_source_url[i], target = "_blank", rel = "noopener", rows$taxonomy_version[i]) else .animal_fit_or(rows$taxonomy_version[i], "Unresolved"),
            " · ", rows$taxonomy_status[i])))))
    })
    output$sources <- shiny::renderUI({
      own <- data()
      # Component identity, rather than a name-only lookup, carries the
      # original citation into both the source list and evidence CSV.
      components <- own$components
      shiny::tagList(shiny::p(class = "atlas-context-line", "Sequence versions: ",
        if (nzchar(own$metadata$sequence_versions)) own$metadata$sequence_versions else "Unspecified",
        ". Taxonomy versions: ", if (nzchar(own$metadata$taxonomy_versions)) own$metadata$taxonomy_versions else "Unresolved / unsourced", "."),
        shiny::p(class = "atlas-context-line", "The evidence CSV retains exact source URLs and hashes, original component citations, effective and parent sequences, raw search states, substitution/length settings and full parent lineages."),
        if (nrow(components)) shiny::tags$ul(lapply(seq_len(nrow(components)), function(i) {
          citation <- if ("citation" %in% names(components)) .animal_fit_text(components$citation[i]) else "Source unavailable"
          url <- if ("source_url" %in% names(components)) .animal_fit_text(components$source_url[i]) else ""
          shiny::tags$li(shiny::strong(components$primer_name[i]), " · ",
            if (nzchar(url)) shiny::tags$a(href = url, target = "_blank", rel = "noopener", citation) else citation,
            if ("catalog_origin" %in% names(components) && components$catalog_origin[i] == "derived_reverse_complement") " · parent provenance only; derived oligo unvalidated")
        })))
    })
    shiny::observeEvent(input$apply_taxonomy, {
      tryCatch({
        parsed <- animal_fit_parse_taxonomy_csv(input$taxonomy_csv, value()$templates$template_id)
        user_taxonomy(parsed); user_note(paste(nrow(parsed), "templates assigned user-provided taxonomy metadata for this session.")); selected_cell(NULL)
      }, error = function(e) shiny::showNotification(conditionMessage(e), type = "error", duration = 8))
    })
    shiny::observeEvent(input$clear_taxonomy, { user_taxonomy(NULL); user_note("Source taxonomy restored."); selected_cell(NULL) })
    output$taxonomy_status <- shiny::renderUI(shiny::p(class = "atlas-context-line", user_note()))
    download_details <- shiny::reactive({
      own <- data()$details; summaries <- scoped()
      if (!nrow(own) || !nrow(summaries)) return(own[FALSE, , drop = FALSE])
      keys <- .animal_fit_tuple(summaries$pair_id, summaries$taxon_key)
      own <- own[.animal_fit_tuple(own$pair_id, own$taxon_key) %in% keys, , drop = FALSE]
      index <- match(.animal_fit_tuple(own$pair_id, own$taxon_key), keys)
      for (field in c("n_panel", "n_exact", "n_substitution", "n_partial", "n_ambiguous", "n_missing", "n_scorable"))
        own[[paste0("taxon_", field)]] <- summaries[[field]][index]
      own
    })
    output$download_summary <- shiny::downloadHandler(filename = function() paste0("primer-atlas-", marker(), "-taxon-counts.csv"),
      content = function(file) utils::write.csv(scoped(), file, row.names = FALSE, na = ""))
    output$download_details <- shiny::downloadHandler(filename = function() paste0("primer-atlas-", marker(), "-taxon-template-evidence.csv"),
      content = function(file) utils::write.csv(download_details(), file, row.names = FALSE, na = ""))
    output$download_taxonomy <- shiny::downloadHandler(filename = function() paste0("primer-atlas-", marker(), "-template-taxonomy.csv"),
      content = function(file) {
        rows <- data()$templates
        # This file can be edited and pasted directly into the metadata input.
        # Full source/version/hash fields remain in the template evidence CSV.
        fields <- intersect(c("template_id", .animal_fit_ranks), names(rows))
        utils::write.csv(rows[, fields, drop = FALSE], file, row.names = FALSE, na = "")
      })
    list(data = data, selected_pairs = selected_pairs, selected_details = selected_details,
      scoped = scoped, page = current_page, user_taxonomy = user_taxonomy, download_details = download_details)
  })
}
