# Compact controls for the shared primer overview; the server keeps the two
# oligo selections authoritative and computes their Cartesian product.
fungalPrimerAtlasNavigatorUI <- function() {
  shiny::tagList(
    shiny::tags$script(shiny::HTML("if(!window.fungalAtlasSelection){window.fungalAtlasSelection=true;function pick(b){if(b&&window.Shiny)Shiny.setInputValue(b.dataset.primerInput,{catalog_id:b.dataset.catalogId,direction:b.dataset.catalogDirection||b.dataset.direction,intent:'cycle',version:b.dataset.versionId,accession:b.dataset.accession,nonce:Date.now()},{priority:'event'});}document.addEventListener('click',function(e){pick(e.target.closest('.fungal-atlas-primer[data-primer-input]'));});document.addEventListener('keydown',function(e){var b=e.target.closest('.fungal-atlas-primer[data-primer-input]');if(b&&!e.repeat&&(e.key==='Enter'||e.key===' ')){e.preventDefault();pick(b);}});}")),
    shiny::tags$style(shiny::HTML(".fungal-atlas-intro{font-size:.88rem;color:#526176;margin:0 0 10px}.fungal-atlas-selection{display:flex;align-items:center;justify-content:space-between;flex-wrap:wrap;gap:12px;margin:12px 0 16px;padding:12px 16px;border:1px solid #dbe3ee;border-radius:10px;background:#fff;font-size:.85rem}.fungal-atlas-selection>span{color:#64748b}")),
    shiny::div(class="atlas-segmented fungal-atlas-mode",
      shiny::radioButtons("fungal_overview_mode",NULL,
        c("Primer atlas"="atlas","Selected amplicons"="amplicons"),"atlas",inline=TRUE)),
    shiny::uiOutput("fungal_atlas_selection"),
    shiny::conditionalPanel(condition="input.fungal_overview_mode == 'atlas'",
      shiny::p(class="fungal-atlas-intro","Click a primer to select it. Click again to reverse its direction; a third click removes it. Blue arrows are forward, red arrows reverse. Reversed variants use reverse-complemented sequences."),
      shiny::uiOutput("fungal_primer_atlas")))
}

fungalCombinationBuilderUI <- function() {
  shiny::tagList(
    shiny::tags$script(shiny::HTML("if(!window.atlasBuilderInspect){window.atlasBuilderInspect=true;document.addEventListener('click',function(e){var b=e.target.closest('.pair-inspect[data-pair-input]');if(b&&window.Shiny)Shiny.setInputValue(b.dataset.pairInput,{pair:b.dataset.pairId,nonce:Date.now()},{priority:'event'});});}")),
    shiny::tags$style(shiny::HTML(".fungal-builder{border:1px solid #dbe3ee;border-radius:14px;background:#fff;padding:16px;min-width:0}.fungal-builder-heading{display:flex;justify-content:space-between;align-items:center;gap:12px;margin-bottom:12px}.fungal-builder-heading h3{font-size:1rem;margin:0}.fungal-builder-heading p{font-size:.78rem;color:#64748b;margin:0}.fungal-builder-fields{display:grid;grid-template-columns:minmax(0,1fr) minmax(0,1fr);gap:16px}.fungal-builder-fields .shiny-input-container{width:100%!important;margin:0;min-width:0}.fungal-builder-fields .selectize-input{min-height:72px;border-radius:8px;padding:9px}.fungal-builder-fields label{font-size:.82rem;font-weight:650}.fungal-builder-filterbar{display:flex;align-items:end;flex-wrap:wrap;gap:12px;margin-top:12px}.fungal-builder-filterbar .shiny-input-container{width:140px!important;max-width:100%;margin:0}.fungal-builder-filterbar label{font-size:.78rem}.fungal-builder-actions{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin-top:12px}.fungal-builder-summary{font-size:.82rem;color:#475569;line-height:1.5}.fungal-builder-table-wrap{overflow-x:auto}.fungal-builder-table{width:100%;border-collapse:collapse;font-size:.8rem}.fungal-builder-table th,.fungal-builder-table td{padding:9px 8px;border-bottom:1px solid #e2e8f0;text-align:left;vertical-align:top}.fungal-builder-table th{color:#526176;font-weight:650}.fungal-builder-table .pair-inspect{border:0;background:transparent;color:#4338ca;font:inherit;font-weight:650;text-align:left;padding:0;cursor:pointer}.fungal-publication-actions{display:flex;gap:8px;flex-wrap:wrap}.fungal-oligo-evidence{display:grid;grid-template-columns:minmax(0,1fr) minmax(0,1fr);gap:12px}.fungal-oligo-evidence>div{padding:12px;border:1px solid #e2e8f0;border-radius:10px;font-size:.8rem;overflow-wrap:anywhere}.fungal-oligo-evidence code{font-size:.75rem}.fungal-oligo-evidence p{margin:6px 0}.fungal-publication-scope{font-size:.8rem;color:#526176;line-height:1.5}@media(max-width:600px){.fungal-builder-fields,.fungal-oligo-evidence{grid-template-columns:1fr}.fungal-builder-heading{align-items:start;flex-direction:column}.fungal-builder-filterbar .shiny-input-container{flex:1 1 120px;width:auto!important}}")),
    shiny::div(class="fungal-builder",
      shiny::div(class="fungal-builder-heading",
        shiny::h3("Your primer selection"),
        shiny::p("UNITE catalog · select on the map or search here · up to five per direction")),
      shiny::uiOutput("fungal_builder_controls"),
      shiny::conditionalPanel(condition="input.fungal_overview_mode == 'amplicons'",
      shiny::div(class="fungal-builder-filterbar",
        shiny::numericInput("fungal_builder_min_bp","Minimum product bp",100,min=100,step=25),
        shiny::numericInput("fungal_builder_max_bp","Maximum product bp",3500,min=100,step=25))),
      shiny::div(class="fungal-builder-actions",
        shiny::actionButton("fungal_builder_examples","Load 5 × 5 example",class="btn-outline-secondary btn-sm")),
      shiny::tags$details(class="atlas-detail",
        shiny::tags$summary("Refine evidence and inspect all requested combinations"),
        shiny::selectInput("fungal_builder_evidence","Publication evidence",
          c("All combinations"="ALL","Documented pairs only"="DOCUMENTED","Both primers have cited papers"="PUBLISHED"),"ALL"),
        shiny::p(class="fungal-builder-summary","Length filters use the selected reference's product, including both primers. Missing positions remain available below for investigation in other taxa."),
        shiny::uiOutput("fungal_builder_combination_table"))))
}

fungalPublicationUI <- function() {
  shiny::tags$details(id="fungal-publication-export",class="atlas-detail",
    shiny::tags$summary("3 · Export figure, table and citations"),
    shiny::p(class="fungal-publication-scope","Exports include all mapped combinations that pass the current filters, across every map page. The complete requested-combination table retains missing-site and filtered states."),
    shiny::p(class="fungal-publication-scope","To include a DNA alignment figure, open Regional alignment above and choose its group, region and window. The package includes that exact displayed alignment and its FASTA."),
    shiny::textInput("fungal_publication_title","Figure title","Fungal primer combinations"),
    shiny::div(class="fungal-builder-fields",
      shiny::selectInput("fungal_publication_rank","Summarize reference evidence by",c("Classes"="class","Orders"="order"),"class"),
      shiny::uiOutput("fungal_publication_taxon_control")),
    shiny::div(class="fungal-publication-actions",
      shiny::downloadButton("fungal_download_figure","Primer figure (SVG)",class="btn-outline-primary btn-sm"),
      shiny::uiOutput("fungal_publication_png_control"),
      shiny::downloadButton("fungal_download_table","Combination table (CSV)",class="btn-outline-primary btn-sm"),
      shiny::downloadButton("fungal_download_print_table","Publication table (HTML)",class="btn-outline-secondary btn-sm"),
      shiny::downloadButton("fungal_download_package","Figure and evidence package (ZIP)",class="btn-primary btn-sm")),
    shiny::p(class="fungal-publication-scope","SVG is editable vector artwork. The package records the exact reference, snapshot, filters, primer sources, taxonomic denominators, figure caption and methods. These computational observations do not establish PCR performance."))
}
