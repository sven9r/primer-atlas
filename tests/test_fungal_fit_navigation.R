#!/usr/bin/env Rscript
# Exercise cross-page navigation in the real app server. Evidence calculations
# have separate unit checks; this test guards routing and immutable context.
.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
env <- new.env()
source("app.R", local = env)
# Capture the app's actual download callback and still register it with Shiny.
# This avoids reaching into private MockShinySession download internals.
downloads <- list()
env$downloadHandler <- function(filename, content, ...) {
  downloads[[length(downloads) + 1L]] <<- list(filename = filename, content = content)
  shiny::downloadHandler(filename = filename, content = content, ...)
}
ui <- as.character(env$ui)
stopifnot(grepl("Taxonomic fit", ui, fixed = TRUE),
          !grepl("Check taxonomic fit · COI", ui, fixed = TRUE))
for (route in c("start", "find", "fit", "taxa", "regions", "thermal", "sources")) {
  stopifnot(grepl(paste0('data-value="', route, '"'), ui, fixed = TRUE))
}
buttons <- c("go_fit", "go_fungal_fit", "go_fungal_explore", "go_coi_explore",
             "go_fungal_library", "go_coi_library", "go_fungal_fit_library",
             "go_find", "go_18s", "go_18s_start")
for (id in buttons) {
  matches <- regmatches(ui, gregexpr(paste0('id="', id, '"'), ui, fixed = TRUE))[[1L]]
  stopifnot(length(matches) == 1L)
}
for (id in c("fit_marker", "fungal-open_fit", "fungal_fit-open_explorer")) {
  stopifnot(grepl(paste0('id="', id, '"'), ui, fixed = TRUE))
}

versions <- fungal_read_index(remote = FALSE)
v1_id <- versions$id[startsWith(versions$id, "v1-")][1]
v3_id <- tail(versions$id, 1)
snapshot <- fungal_read_snapshot(v3_id, remote = FALSE)
pair <- snapshot$pairs$pair_id[snapshot$pairs$pair_label == "ITS1 + ITS4"][1]
stopifnot(length(pair) == 1L, !is.na(pair))

shiny::testServer(env$server, {
  session$setInputs(main_nav = "start", fit_marker = "COI",
    map_organism = "ALL", map_marker = "COI", map_sort = "binding_site",
    amplicon_range = c(100, 1800), pair_select = c("MCO", "ANML"),
    application_filter = character(), environment_filter = character(),
    intent_filter = character(), show_custom = TRUE,
    order_pair_select = c("MCO", "ANML"), detail_order = "Araneae",
    detail_pair = "MCO", detail_example = "p90",
    lineage_marker = "COI", lineage_pair = "MCO", lineage_target = "Araneae",
    lineage_country = "All", lineage_locality = "All",
    "fungal-version" = v3_id, "fungal-platform" = "ALL",
    "fungal-illumina_max" = 500, "fungal-pacbio_max" = 2000,
    "fungal-nanopore_max" = 3500, "fungal-phylum" = "ALL",
    "fungal-class" = "ALL", "fungal-locus" = "ALL", "fungal-search" = "",
    "fungal-examples_only" = TRUE, "fungal-show_unresolved" = TRUE,
    "fungal_fit-version" = v3_id, "fungal_fit-pair" = pair,
    "fungal_fit-rank" = "class", "fungal_fit-platform" = "ALL")
  session$flushReact()
  sent <- list()
  session$sendInputMessage <- function(inputId, message) sent[[inputId]] <<- message

  # Each landing entry specifies the evidence mode and uses the same fit page.
  session$setInputs(go_fungal_fit = 1L)
  session$flushReact()
  stopifnot(sent$main_nav$value == "fit", sent$fit_marker$value == "ITS_FUNGAL")
  session$setInputs(fit_marker = "ITS_FUNGAL")
  session$flushReact()
  session$setInputs(go_fit = 1L)
  session$flushReact()
  stopifnot(sent$main_nav$value == "fit", sent$fit_marker$value == "COI",
    input$detail_order == "Araneae", input$detail_pair == "MCO",
    identical(input$order_pair_select, c("MCO", "ANML")))

  # Marker-specific primer entries cannot inherit another marker's old scope.
  session$setInputs(go_fungal_explore = 1L)
  session$flushReact()
  stopifnot(sent$main_nav$value == "find", sent$map_marker$value == "ITS_FUNGAL",
            sent$map_organism$value == "FUNGI")
  session$setInputs(go_coi_explore = 1L)
  session$flushReact()
  stopifnot(sent$main_nav$value == "find", sent$map_marker$value == "COI")
  session$setInputs(go_find = 1L)
  session$flushReact()
  stopifnot(sent$main_nav$value == "find")

  # Explorer→fit transfers exactly one snapshot/pair, not a current catalog
  # selection assembled independently by the destination module.
  session$setInputs("fungal-pair_click" = list(pair = pair, nonce = 1L))
  session$flushReact()
  stopifnot(fungal_workspace$snapshot()$manifest$id == v3_id,
            fungal_workspace$active_pair()$pair_id == pair)
  session$setInputs("fungal-open_fit" = 1L)
  session$flushReact()
  context <- fungal_fit_context()
  stopifnot(context$version == v3_id, context$pair == pair,
            length(context$nonce) == 1L, sent$main_nav$value == "fit",
            sent$fit_marker$value == "ITS_FUNGAL")

  # Apply browser input updates explicitly: MockShinySession records an
  # updateSelectInput message but does not synthesize the client response.
  session$setInputs(fit_marker = "ITS_FUNGAL", "fungal_fit-version" = v3_id,
                    "fungal_fit-pair" = pair)
  session$flushReact()
  stopifnot(fungal_fit_workspace$snapshot()$manifest$id == v3_id,
            fungal_fit_workspace$selected_pair() == pair,
            sum(fungal_fit_workspace$fit_data()$summary$n_panel) == length(snapshot$sequences),
            nrow(fungal_fit_workspace$fit_data()$details) == length(snapshot$sequences),
            is.null(fungal_fit_workspace$selected_taxon()))
  group <- fungal_fit_workspace$fit_data()$summary$group_id[1]
  session$setInputs("fungal_fit-taxon_click" = list(group = group, nonce = 1L))
  session$flushReact()
  details <- fungal_fit_workspace$taxon_details()
  stopifnot(fungal_fit_workspace$selected_taxon() == group, nrow(details) > 0L,
            all(details$group_id == group), all(details$pair_id == pair),
            all(details$version_id == v3_id))

  session$setInputs(go_fungal_library = 1L)
  session$flushReact()
  stopifnot(sent$lineage_organism$value == "FUNGI", sent$lineage_marker$value == "ITS_FUNGAL")
  session$setInputs(lineage_organism = "FUNGI", lineage_marker = "ITS_FUNGAL")
  session$flushReact()
  stopifnot(grepl(paste(length(snapshot$sequences), "retained templates"), output$fungal_library_context$html, fixed = TRUE))
  matches <- vapply(downloads, function(value) {
    filename <- tryCatch(value$filename(), error = function(e) "")
    if (length(filename) != 1L) return(FALSE)
    grepl("-reference-metadata[.]csv$", filename)
  }, logical(1))
  stopifnot(sum(matches) == 1L)
  library_download <- downloads[[which(matches)]]
  file <- tempfile(fileext = ".csv")
  library_download$content(file)
  records <- read.csv(file, stringsAsFactors = FALSE)
  stopifnot(nrow(records) == length(snapshot$sequences),
            setequal(records$accession, snapshot$taxonomy$accession),
            all(records$snapshot_id == v3_id),
            library_download$filename() == paste0(v3_id, "-reference-metadata.csv"))
  unlink(file)
  session$setInputs(go_fungal_fit_library = 1L)
  session$flushReact()
  stopifnot(sent$main_nav$value == "fit", sent$fit_marker$value == "ITS_FUNGAL")
  session$setInputs("fungal_fit-open_explorer" = 1L)
  session$flushReact()
  context <- fungal_explorer_context()
  stopifnot(context$version == v3_id, context$pair == pair,
            length(context$nonce) == 1L, sent$main_nav$value == "find",
            sent$map_marker$value == "ITS_FUNGAL")

  # Returning from historical evidence must preserve that version too.
  historical <- fungal_read_snapshot(v1_id, remote = FALSE)
  old_pair <- historical$pairs$pair_id[historical$pairs$n_scorable > 0][1]
  session$setInputs("fungal_fit-version" = v1_id, "fungal_fit-pair" = old_pair)
  session$flushReact()
  stopifnot(length(fungal_fit_workspace$snapshot()$sequences) == 1L,
            fungal_fit_workspace$snapshot()$manifest$id == v1_id,
            sum(fungal_fit_workspace$fit_data()$summary$n_panel) == 1L,
            nrow(fungal_fit_workspace$fit_data()$details) == 1L,
            is.null(fungal_fit_workspace$selected_taxon()))
  file <- tempfile(fileext = ".csv")
  library_download$content(file)
  records <- read.csv(file, stringsAsFactors = FALSE)
  stopifnot(nrow(records) == 1L, records$accession == "FN812768.2",
            records$snapshot_id == v1_id,
            library_download$filename() == paste0(v1_id, "-reference-metadata.csv"))
  unlink(file)
  session$setInputs(go_coi_library = 1L)
  session$flushReact()
  stopifnot(sent$lineage_organism$value == "ALL", sent$lineage_marker$value == "COI")
  session$setInputs("fungal_fit-open_explorer" = 2L)
  session$flushReact()
  stopifnot(fungal_explorer_context()$version == v1_id,
            fungal_explorer_context()$pair == old_pair)

  # The existing COI inspector still opens exact taxa/sequences with its own
  # pair/order. Fungal mode/context changes must not leak into these controls.
  session$setInputs(open_lineage = 1L)
  session$flushReact()
  stopifnot(sent$main_nav$value == "taxa", sent$lineage_marker$value == "COI",
            sent$lineage_pair$value == "MCO", sent$lineage_target$value == "Araneae",
            input$detail_pair == "MCO", input$detail_order == "Araneae",
            identical(input$order_pair_select, c("MCO", "ANML")))

  # An empty explorer selection cannot silently open an unrelated prior/default
  # fit pair. Both the exact context and current route remain unchanged.
  session$setInputs("fungal-search" = "no_such_pair_987654321")
  session$flushReact()
  stopifnot(nrow(fungal_workspace$active_pair()) == 0L)
  previous_context <- fungal_fit_context()
  previous_route <- sent$main_nav$value
  session$setInputs("fungal-open_fit" = 2L)
  session$flushReact()
  stopifnot(identical(fungal_fit_context(), previous_context),
            identical(sent$main_nav$value, previous_route))
})
message("Animal/fungal fit routes, frozen-context round trips, exact library CSVs, historical/empty states and preserved COI drill-down passed.")
