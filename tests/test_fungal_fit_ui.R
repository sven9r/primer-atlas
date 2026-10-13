#!/usr/bin/env Rscript
.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
library(shiny)
library(DT)
source("R/functions.R")
source("R/data_layer.R")
source("R/fungal_explorer.R")
source("R/fungal_versions.R")
source("R/fungal_fit.R")
source("R/fungal_fit_ui.R")
Sys.unsetenv("ATLAS_MANIFEST_BASE_URL")
versions <- fungal_read_index(remote = FALSE)
latest <- tail(versions$id, 1L)
old <- versions$id[startsWith(versions$id, "v1-")][1L]
data <- fungal_read_snapshot(latest, remote = FALSE)
pair <- data$pairs$pair_id[match("ITS1 + ITS4", data$pairs$pair_label)]
alternate <- data$pairs$pair_id[match("LR0R + LR3", data$pairs$pair_label)]
context <- reactiveVal(NULL)

shiny::testServer(fungalFitServer, args = list(id = "fit", selection = context), {
  session$setInputs(version = latest, pair = pair, rank = "class", platform = "ALL", search = "")
  session$flushReact()
  expected <- fungal_fit_summary(data, rank = "class", pair_ids = pair)
  stopifnot(nrow(fit_data()$details) == length(data$sequences),
      sum(fit_data()$summary$n_panel) == length(data$sequences),
      is.null(selected_taxon()), nrow(taxon_details()) == 0L,
      is.null(output$taxon_inspector))
  cards <- output$taxa$html
  for (i in seq_len(nrow(expected$summary))) {
    row <- expected$summary[i, , drop = FALSE]
    stopifnot(grepl(as.character(row$taxon), cards, fixed = TRUE),
      grepl(paste(row$n_exact_products, "of", row$n_panel), cards, fixed = TRUE))
  }
  chosen <- expected$summary$group_id[which.max(expected$summary$n_panel)]
  session$setInputs(taxon_click = list(group = chosen, nonce = 1L))
  session$flushReact()
  details <- taxon_details()
  stopifnot(selected_taxon() == chosen, nrow(details) > 0L,
      all(details$group_id == chosen), all(details$pair_id == pair),
      all(details$version_id == latest),
      grepl('fit-references', output$taxon_inspector$html, fixed = TRUE))
  before <- fit_data()$summary$n_exact_products
  session$setInputs(platform = "ILLUMINA")
  session$flushReact()
  stopifnot(identical(before, fit_data()$summary$n_exact_products),
      nrow(taxon_details()) == nrow(details), all(taxon_details()$platform == "ILLUMINA"),
      grepl("encoded bp ceiling", output$taxa$html, fixed = TRUE))
  session$setInputs(search = "a-taxon-not-in-this-panel")
  session$flushReact()
  stopifnot(nrow(groups()) == 0L, is.null(selected_taxon()), nrow(taxon_details()) == 0L,
      is.null(output$taxon_inspector), grepl("No taxa match", output$taxa$html, fixed = TRUE))
  session$setInputs(search = "", pair = alternate)
  session$flushReact()
  stopifnot(selected_pair() == alternate, all(fit_data()$details$pair_id == alternate), is.null(selected_taxon()))
  # A cross-page context pins version/pair until the real browser applies the
  # updates; a later direct pair choice must resume normal user control.
  context(list(version = latest, pair = pair, nonce = 1L))
  session$flushReact()
  stopifnot(selected_pair() == pair)
  session$setInputs(pair = pair)
  session$flushReact()
  session$setInputs(pair = alternate)
  session$flushReact()
  stopifnot(selected_pair() == alternate)
  historical <- fungal_read_snapshot(old, remote = FALSE)
  old_pair <- historical$pairs$pair_id[match("ITS3 + ITS4", historical$pairs$pair_label)]
  session$setInputs(version = old, pair = old_pair, rank = "class", platform = "ALL")
  session$flushReact()
  stopifnot(nrow(fit_data()$details) == 1L, sum(fit_data()$summary$n_panel) == 1L,
      all(fit_data()$details$version_id == old), is.null(selected_taxon()),
      grepl("Unresolved class", output$taxa$html, fixed = TRUE))
})
message("Rendered fungal fit denominators, card drill-down, independent platform planning, search reset, historical versions and context release passed.")
