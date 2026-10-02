#!/usr/bin/env Rscript

app <- source("app.R")$value
publication_bees <- as.data.frame(beeprime_hymenoptera_scores |> filter(is_bee %in% TRUE))
stopifnot(nrow(publication_bees) == 590L,
          !beeprime_expanded_reference(NULL, publication_bees)$available)

# A new reference must retain the publication baseline and actually add bees.
expanded <- publication_bees
expanded$pair_id <- "BEEPRIME"
expanded$reference_origin <- "publication"
new <- expanded[1, ]
new$accession <- "ATLAS_TEST_NEW_BEE.1"
new$cluster_size <- 1L
new$reference_origin <- "added"
new$pair_penalty <- 999
new$forward_penalty <- 999
new$reverse_penalty <- 0
expanded <- rbind(expanded, new)
scratch <- tempfile("reference-choice-")
dir.create(scratch)
on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
score_path <- file.path(scratch, "scores.csv")
write.csv(expanded, score_path, row.names = FALSE, na = "")
descriptor <- list(
  pair_id = "BEEPRIME", reference_mode = "expanded", taxonomic_scope = "bees",
  publication_baseline = "gurten2026_suppl5", version = "test-expanded-1",
  added_sequences = 1L,
  sequence_scores = list(local_path = score_path,
    sha256 = digest::digest(file = score_path, algo = "sha256", serialize = FALSE))
)
verified <- beeprime_expanded_reference(descriptor, publication_bees)
stopifnot(verified$available, nrow(verified$scores) == 591L,
          verified$added_sequences == 1L,
          sum(verified$tables$bee_family$n_centroids) == 591L,
          all(is.na(verified$geography$country_or_territory)),
          all(beeprime_reference_table(verified, "bee_family")$reference_version == "test-expanded-1"))
wrong_scope <- descriptor
wrong_scope$taxonomic_scope <- "all_arthropods"
stopifnot(!beeprime_expanded_reference(wrong_scope, publication_bees)$available)
wrong_checksum <- descriptor
wrong_checksum$sequence_scores$sha256 <- strrep("0", 64)
stopifnot(!beeprime_expanded_reference(wrong_checksum, publication_bees)$available)
no_growth <- descriptor
no_growth$added_sequences <- 0L
stopifnot(!beeprime_expanded_reference(no_growth, publication_bees)$available)
incomplete_baseline <- expanded[-1, ]
write.csv(incomplete_baseline, file.path(scratch, "incomplete.csv"), row.names = FALSE)
missing <- descriptor
missing$sequence_scores$local_path <- file.path(scratch, "incomplete.csv")
missing$sequence_scores$sha256 <- digest::digest(file = missing$sequence_scores$local_path, algo = "sha256", serialize = FALSE)
stopifnot(!beeprime_expanded_reference(missing, publication_bees)$available)

shiny::testServer(app$serverFuncSource(), {
  session$setInputs(map_marker = "COI", map_organism = "ALL", amplicon_range = c(110, 710),
    pair_select = "MCO", sodium = 50, ta_offset = 4,
    beeprime_reference_mode = "publication", beeprime_taxon_scope = "bee_family",
    beeprime_taxon_sort = "taxon", lineage_marker = "COI", lineage_pair = "BEEPRIME",
    lineage_target = "All", lineage_country = "All", lineage_locality = "All",
    region_pair = "BEEPRIME", region_marker = "COI", detail_pair = "BEEPRIME",
    detail_order = "Hymenoptera", claimed_taxon_rank = "family")
  stopifnot(sum(beeprime_taxon_data()$n_centroids) == 590L)
  session$setInputs(lineage_reference_mode = "expanded")
  stopifnot(beeprime_reference_mode() == "expanded", !beeprime_reference()$available,
            nrow(beeprime_taxon_data()) == 0L, nrow(lineage_raw_scores()) == 0L,
            nrow(region_pair_scores()) == 0L, nrow(claimed_taxon_data()) == 0L,
            grepl("has not been published", output$lineage_reference_status$html, fixed = TRUE))
  session$setInputs(region_reference_mode = "publication")
  stopifnot(beeprime_reference_mode() == "publication",
            sum(beeprime_taxon_data()$n_centroids) == 590L)
})

old_release <- coi_release_state
coi_release_state$manifest$reference_panels$beeprime_expanded <- descriptor
shiny::testServer(app$serverFuncSource(), {
  session$setInputs(map_marker = "COI", map_organism = "ALL", amplicon_range = c(110, 710),
    pair_select = "MCO", sodium = 50, ta_offset = 4,
    beeprime_reference_mode = "publication", beeprime_taxon_scope = "bee_family",
    beeprime_taxon_sort = "taxon", lineage_marker = "COI", lineage_pair = "BEEPRIME",
    lineage_target = "All", lineage_country = "All", lineage_locality = "All",
    region_pair = "BEEPRIME", region_marker = "COI", detail_pair = "BEEPRIME",
    detail_order = "Hymenoptera", claimed_taxon_rank = "family")
  session$setInputs(beeprime_reference_mode = "expanded")
  stopifnot(beeprime_reference()$available, sum(beeprime_taxon_data()$n_centroids) == 591L,
            nrow(lineage_raw_scores()) == 591L, nrow(region_pair_scores()) == 591L,
            all(beeprime_taxon_data()$reference_version == "test-expanded-1"),
            sum(claimed_taxon_data()$n_centroids) == 591L)
  # General primers cannot inherit a selected bee alignment.
  session$setInputs(lineage_pair = "MCO", region_pair = "MCO")
  stopifnot(!"ATLAS_TEST_NEW_BEE.1" %in% lineage_raw_scores()$accession,
            !"ATLAS_TEST_NEW_BEE.1" %in% region_pair_scores()$accession,
            pair_sequence_evidence_spec("MCO")$evidence_class == "full_order_reference_alignment")
  session$setInputs(beeprime_reference_mode = "publication")
  stopifnot(sum(beeprime_taxon_data()$n_centroids) == 590L)
})
coi_release_state <- old_release
message("Publication/expanded reference routing, provenance, integrity and no-fallback checks passed.")
