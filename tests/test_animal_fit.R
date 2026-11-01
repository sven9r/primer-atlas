#!/usr/bin/env Rscript
source("R/animal_fit.R")
source("R/animal_fit_ui.R")
source("R/primer_mixer.R")
expect_error <- function(expr) inherits(try(force(expr), silent = TRUE), "try-error")

templates <- data.frame(template_id = paste0("t", 1:7), marker_id = "COI",
  sequence = rep("ACGTACGTACGT", 7), organism = paste("Fixture", 1:7),
  kingdom = "Metazoa", phylum = "Chordata", class = "Class A",
  order = c("Order A", "Order A", "Order B", "Order B", "Order C", "Order C", "Unresolved"),
  family = c("Shared family", "Shared family", "Shared family", "Unresolved", "Family C", "Family C", "Unresolved"),
  genus = "Unresolved", species = "Unresolved", version = "fixture-v1",
  source_url = "https://example.org/fixture-sequence", taxonomy_status = "fixture_source",
  taxonomy_version = "fixture-taxonomy-v1", taxonomy_source_url = "https://example.org/fixture-taxonomy",
  taxonomy_source_sha256 = "fixture-hash", stringsAsFactors = FALSE)
templates$marker_id[6L] <- "18S"
templates[7L, c("kingdom", "phylum", "class")] <- "Unresolved"
templates$version[7L] <- "session_input"
templates$organism[7L] <- "Danio rerio in a user header"
pairs <- data.frame(pair_id = c("pair-A", "pair-B"), pair_label = c("F1 + R1", "F1 + R2"))
products <- expand.grid(template_id = templates$template_id, pair_id = pairs$pair_id,
  stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
products$status <- rep(c("exact_product", "mismatch_product", "one_eligible_site", "ambiguous_sites",
  "no_eligible_sites", "non_inward_or_overlapping_sites", "template_unavailable"), 2L)
products$forward_status <- rep(c("exact_site", "exact_site", "exact_site", "multiple_equally_good_sites",
  "no_site_within_mismatch_limit", "exact_site", "sequence_unavailable"), 2L)
products$reverse_status <- rep(c("exact_site", "mismatch_site", "no_site_within_mismatch_limit", "exact_site",
  "no_site_within_mismatch_limit", "exact_site", "sequence_unavailable"), 2L)
products$passes_length_filter <- products$status %in% c("exact_product", "mismatch_product")
products$amplicon_bp <- ifelse(products$passes_length_filter, 600L, NA_integer_)
products$forward_mismatches <- 0L; products$reverse_mismatches <- ifelse(products$status == "mismatch_product", 1L, 0L)
result <- list(pairs = pairs, templates = templates, products = products, marker_id = "COI",
  settings = list(max_mismatches = 1L, min_bp = 0, max_bp = 10000), method = "Explicit controlled fixture")
saved <- serialize(result, NULL)
class <- animal_fit_data(result)
stopifnot(identical(serialize(result, NULL), saved), sum(class$summary$n_panel[class$summary$pair_id == "pair-A"]) == 7L)
for (pair in pairs$pair_id) {
  own <- class$summary[class$summary$pair_id == pair, ]
  stopifnot(sum(own$n_exact) == 1L, sum(own$n_substitution) == 1L, sum(own$n_partial) == 1L,
    sum(own$n_ambiguous) == 1L, sum(own$n_missing) == 3L,
    sum(own$n_scorable) == 2L,
    all(rowSums(own[, c("n_exact", "n_substitution", "n_partial", "n_ambiguous", "n_missing")]) == own$n_panel))
}
wrong_marker <- class$details[class$details$template_id == "t6", ]
stopifnot(all(wrong_marker$status == "marker_mismatch"), !any(wrong_marker$two_unique_sites),
  all(class$details$taxon[class$details$template_id == "t7"] == "Unresolved class"))
family <- animal_fit_data(result, "family")
shared <- family$summary[family$summary$pair_id == "pair-A" & family$summary$taxon == "Shared family", ]
stopifnot(nrow(shared) == 2L, length(unique(shared$taxon_key)) == 2L, setequal(shared$n_panel, c(1L, 2L)),
  all(grepl("Order [AB]", shared$parent_path)))

# Missing and contradictory rows remain in the denominator as unavailable.
incomplete <- result; incomplete$products <- products[-1L, ]
missing <- animal_fit_data(incomplete)
stopifnot(sum(missing$summary$n_panel[missing$summary$pair_id == "pair-A"]) == 7L,
  missing$details$status[missing$details$template_id == "t1" & missing$details$pair_id == "pair-A"] == "evidence_unavailable")
duplicated <- result; duplicated$products <- rbind(products, products[1L, ])
duplicates <- animal_fit_data(duplicated)
stopifnot(duplicates$details$status[duplicates$details$template_id == "t1" & duplicates$details$pair_id == "pair-A"] == "inconsistent_duplicate_evidence",
  sum(duplicates$summary$n_exact[duplicates$summary$pair_id == "pair-A"]) == 0L)

# Product ceilings affect planning counts only, leaving all fit denominators
# and raw exact/substitution categories unchanged.
filtered <- result; filtered$products$passes_length_filter <- FALSE
filtered$settings$max_bp <- 500
ceiling <- animal_fit_data(filtered)
stopifnot(identical(class$summary[, c("pair_id", "taxon_key", "n_panel", "n_exact", "n_substitution", "n_scorable")],
  ceiling$summary[, c("pair_id", "taxon_key", "n_panel", "n_exact", "n_substitution", "n_scorable")]),
  sum(class$summary$n_length_eligible) == 4L, sum(ceiling$summary$n_length_eligible) == 0L)

# An accession-like user header cannot inherit supplemental source taxonomy.
fake_taxonomy <- templates[7L, c("template_id", "marker_id", "kingdom", "phylum", "class")]
fake_taxonomy$kingdom <- "Metazoa"; fake_taxonomy$phylum <- "Chordata"; fake_taxonomy$class <- "Actinopteri"
guarded <- animal_fit_data(result, taxonomy = fake_taxonomy)
stopifnot(all(guarded$details$taxon[guarded$details$template_id == "t7"] == "Unresolved class"))
csv <- "template_id,phylum,class,order,family,genus,species\nt7,Chordata,Actinopteri,Cypriniformes,Danionidae,Danio,Danio rerio"
metadata <- animal_fit_parse_taxonomy_csv(csv, templates$template_id)
assigned <- animal_fit_data(result, user_taxonomy = metadata)
own <- assigned$details[assigned$details$template_id == "t7", ]
stopifnot(all(own$taxon == "Actinopteri"), all(own$taxonomy_status == "user_provided_metadata"),
  all(own$taxonomy_version == "session_metadata"), all(own$source_url == templates$source_url[7L]),
  nchar(metadata$taxonomy_source_sha256) == 64L,
  expect_error(animal_fit_parse_taxonomy_csv(paste0(csv, "\nt7,Chordata,Actinopteri,,,,"), templates$template_id)),
  expect_error(animal_fit_parse_taxonomy_csv("template_id,class\nnot_loaded,Aves", templates$template_id)),
  expect_error(animal_fit_parse_taxonomy_csv("template_id,sequence\nt1,ACGT", templates$template_id)))

# Source-backed starter results exercise every offered fish marker. Their
# lineages arrive from the record-taxid proof, never parsing organism headers.
catalog <- primer_mixer_catalog()
for (marker in c("COI", "12S_MT", "18S", "28S")) {
  cohort <- primer_mixer_templates(marker_id = marker, target_group = "FISH")
  components <- catalog[catalog$marker_id == marker & catalog$target_group %in% "FISH", ]
  f <- utils::head(components$primer_id[components$direction == "forward"], 1L)
  r <- utils::head(components$primer_id[components$direction == "reverse"], 1L)
  actual <- primer_mixer_build(catalog, f, r, cohort, marker)
  fit <- animal_fit_data(actual, "family")
  stopifnot(sum(fit$summary$n_panel) == nrow(cohort), nrow(fit$details) == nrow(cohort),
    all(rowSums(fit$summary[, c("n_exact", "n_substitution", "n_partial", "n_ambiguous", "n_missing")]) == fit$summary$n_panel),
    all(fit$details$taxonomy_status == "source_taxid_verified"))
}

# Real zooplankton records keep their own marker, source ecology and NCBI
# lineage. This checks available evidence and denominators, not an expected hit
# rate or a transfer of fish/arthropod primer performance to Daphnia.
zoo <- primer_mixer_templates(target_group = "ZOOPLANKTON")
stopifnot(nrow(zoo) == 4L, setequal(zoo$marker_id, c("COI", "18S", "28S", "16S_MT")),
  all(zoo$family == "Daphniidae"), all(zoo$class == "Branchiopoda"),
  all(zoo$order == "Diplostraca"), all(nzchar(zoo$ecology_source_url)))
for (marker in c("COI", "18S", "28S", "16S_MT")) {
  cohort <- primer_mixer_templates(marker_id = marker, target_group = "ZOOPLANKTON")
  components <- catalog[catalog$marker_id == marker, , drop = FALSE]
  f <- utils::head(components$primer_id[components$direction == "forward"], 1L)
  r <- utils::head(components$primer_id[components$direction == "reverse"], 1L)
  own_catalog <- catalog
  if (!length(f) || !length(r)) {
    # A marker with no published oligos in this catalog is still investigable
    # through explicitly user-provided candidates. These short test inputs are
    # not represented as published assays and carry no expected exact result.
    custom_f <- primer_mixer_custom("Controlled custom forward", "ACGTACGT", "forward", marker)
    custom_r <- primer_mixer_custom("Controlled custom reverse", "TTGGAACC", "reverse", marker)
    own_catalog <- .primer_mixer_bind(list(catalog, custom_f, custom_r))
    f <- custom_f$primer_id; r <- custom_r$primer_id
  }
  actual <- primer_mixer_build(own_catalog, f, r, cohort, marker)
  fit <- animal_fit_data(actual, "family")
  stopifnot(nrow(cohort) == 1L, nrow(fit$details) == 1L,
    sum(fit$summary$n_panel) == 1L, identical(fit$summary$taxon, "Daphniidae"),
    grepl("Arthropoda", fit$summary$parent_path, fixed = TRUE),
    grepl("Branchiopoda", fit$summary$parent_path, fixed = TRUE),
    grepl("Diplostraca", fit$summary$parent_path, fixed = TRUE),
    all(fit$details$taxonomy_status == "source_taxid_verified"),
    all(fit$details$marker_id == marker),
    rowSums(fit$summary[, c("n_exact", "n_substitution", "n_partial", "n_ambiguous", "n_missing")]) == 1L)
}

# Module defaults prioritize real eligible products, with manual choices kept.
many_pairs <- data.frame(pair_id = paste0("p", 1:6), pair_label = paste("Pair", 1:6))
module_result <- result
module_result$pairs <- many_pairs
module_result$products <- do.call(rbind, lapply(many_pairs$pair_id, function(id) {
  own <- products[products$pair_id == "pair-A", ]; own$pair_id <- id; own
}))
module_result$summary <- data.frame(pair_id = many_pairs$pair_id,
  n_length_eligible = c(0L, 0L, 0L, 0L, 0L, 2L), n_exact_products = c(0L, 0L, 0L, 0L, 1L, 2L))
current <- shiny::reactiveVal(module_result)
shiny::testServer(animalFitServer, args = list(marker = function() "COI", target_group = function() "METAZOA",
  result = function() current()), {
  session$flushReact()
  stopifnot(identical(selected_pairs(), c("p6", "p5", "p1", "p2")),
    grepl("Small selected panel", output$context$html, fixed = TRUE))
  session$setInputs(rank = "family", pairs = "p3", search = "Shared family")
  session$flushReact()
  stopifnot(identical(selected_pairs(), "p3"), nrow(groups()) == 2L)
  group <- groups()$taxon_key[groups()$n_panel == 2L]
  session$setInputs(cell_click = list(group = group, pair = "p3"))
  session$flushReact()
  stopifnot(nrow(selected_details()) == 2L, all(selected_details()$template_id %in% c("t1", "t2")),
    nrow(download_details()) == 3L, all(download_details()$taxon_n_panel %in% c(1L, 2L)))
  session$setInputs(taxonomy_csv = csv, apply_taxonomy = 1L)
  session$flushReact()
  stopifnot(nrow(user_taxonomy()) == 1L, data()$templates$class[data()$templates$template_id == "t7"] == "Actinopteri")
  changed <- module_result; changed$templates$sequence[7L] <- "ACGTACGTACGA"
  current(changed); session$flushReact()
  stopifnot(is.null(user_taxonomy()), identical(selected_pairs(), "p3"))
})

# At most eight lineage cards appear on a page, regardless of collection size.
paged <- result
paged$pairs <- pairs[1L, , drop = FALSE]
paged$templates <- templates[rep(1L, 17L), , drop = FALSE]
paged$templates$template_id <- paste0("page-template-", 1:17)
paged$templates$class <- sprintf("Fixture class %02d", 1:17)
paged$products <- products[rep(1L, 17L), , drop = FALSE]
paged$products$template_id <- paged$templates$template_id
shiny::testServer(animalFitServer, args = list(marker = function() "COI", target_group = function() "METAZOA",
  result = function() paged), {
  session$flushReact()
  cards <- function() length(regmatches(output$matrix$html,
    gregexpr('class="animal-fit-taxon-card"', output$matrix$html, fixed = TRUE))[[1L]])
  stopifnot(cards() == 8L, pages() == 3L)
  session$setInputs(`next` = 1L); session$flushReact()
  stopifnot(current_page() == 2L, cards() == 8L)
  session$setInputs(`next` = 2L); session$flushReact()
  stopifnot(current_page() == 3L, cards() == 1L)
  session$setInputs(search = "Fixture class 17"); session$flushReact()
  stopifnot(current_page() == 1L, cards() == 1L)
})

html <- as.character(animalFitUI("testfit"))
stopifnot(grepl('id="testfit-taxonomy_csv"', html, fixed = TRUE),
  grepl("Headers never assign taxonomy", html, fixed = TRUE),
  grepl('id="testfit-download_details"', html, fixed = TRUE))
message("Animal fit: hierarchy-scoped partitions, missing evidence, source taxonomy, custom metadata, ceilings and module drill-down passed.")
