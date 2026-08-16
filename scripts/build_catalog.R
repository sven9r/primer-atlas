#!/usr/bin/env Rscript

project_root <- normalizePath(getwd(), mustWork = TRUE)
source(file.path(project_root, "R", "functions.R"))

catalog_dir <- file.path(project_root, "data", "catalog")
dir.create(catalog_dir, recursive = TRUE, showWarnings = FALSE)

write_catalog <- function(x, name) {
  write.csv(
    x, file.path(catalog_dir, name), row.names = FALSE, na = "",
    fileEncoding = "UTF-8"
  )
}

markers <- data.frame(
  marker_id = c("COI", "12S_MT", "ITS_FUNGAL", "18S", "16S_PROK", "16S_MT", "28S"),
  display_name = c(
    "Animal COI", "Vertebrate mitochondrial 12S", "Fungal ITS", "Eukaryotic 18S",
    "Bacterial / archaeal 16S", "Animal mitochondrial 16S", "Eukaryotic 28S"
  ),
  molecule_context = c(
    "mitochondrial", "mitochondrial", "nuclear rDNA", "nuclear rDNA",
    "prokaryotic rDNA", "mitochondrial", "nuclear rDNA"
  ),
  coordinate_reference = c(
    "NC_001322.1", "NCBI vertebrate panel release-specific", "FN812768.2", "PR2/SILVA release-specific",
    "SILVA release-specific", "NCBI taxon-panel release-specific",
    "SILVA LSU release-specific"
  ),
  reference_version = c(
    "NC_001322.1", "planned", "FN812768.2", "planned", "planned", "planned", "planned"
  ),
  reference_length = c(1536L, NA, 739L, NA, NA, NA, NA),
  coordinate_model = c(
    "sequence", "sequence_alignment", "sequence_with_landmarks", "landmark_alignment",
    "landmark_alignment", "sequence_alignment", "landmark_alignment"
  ),
  status = c("active", "planned", "pilot", "planned", "planned", "planned", "planned"),
  screening_artifact_id = c(
    "coi_order_alignments", "", "its_fn812768_screen", "", "", "", ""
  ),
  expanded_artifact_id = c(
    "coi_gurten2026", "", "its_unite_future", "", "", "", ""
  ),
  source_name = c("NCBI / PrimerMiner", "NCBI", "UNITE / NCBI", "PR2 / SILVA", "SILVA", "NCBI", "SILVA"),
  source_url = c(
    "https://www.ncbi.nlm.nih.gov/nuccore/NC_001322.1",
    "https://www.ncbi.nlm.nih.gov/",
    "https://unite.ut.ee/primers.php",
    "https://pr2-database.org/", "https://www.arb-silva.de/",
    "https://www.ncbi.nlm.nih.gov/", "https://www.arb-silva.de/"
  ),
  available_lineage_ranks = c(
    "order|family|subfamily|genus|sequence", "planned", "primer_catalog", "planned",
    "planned", "planned", "planned"
  ),
  geographic_fields = c("country|locality|lat_lon", "planned", "planned", "planned", "planned", "planned", "planned"),
  stringsAsFactors = FALSE
)
write_catalog(markers, "markers.csv")

organism_groups <- data.frame(
  group_id = c(
    "METAZOA", "ARTHROPODS", "FISH", "BIRDS", "ZOOPLANKTON",
    "FUNGI", "BACTERIA", "ARCHAEA", "PROTISTS", "PHYTOPLANKTON",
    "PLANT_MICROBIOME"
  ),
  label = c(
    "Animals", "Arthropods", "Fish", "Birds", "Zooplankton",
    "Fungi", "Bacteria", "Archaea", "Protists", "Phytoplankton",
    "Plant microbiomes"
  ),
  icon = c("🐾", "🕷️", "🐟", "🐦", "🦐", "🍄", "🦠", "◉", "🧫", "🌊", "🌿"),
  description = c(
    "Broad metazoan diversity", "Insects, arachnids and other arthropods",
    "Marine and freshwater fishes", "Avian environmental DNA",
    "Metazoan plankton communities", "Fungal communities and phylogenetic profiling",
    "Environmental and host-associated bacteria", "Archaeal community profiling",
    "Microbial eukaryotes and protists", "Algal and protist plankton",
    "Bacteria associated with roots, leaves and rhizospheres"
  ),
  sort_order = seq_len(11),
  stringsAsFactors = FALSE
)
write_catalog(organism_groups, "organism_groups.csv")

marker_group_rows <- list(
  COI = c(METAZOA = "primary", ARTHROPODS = "primary", FISH = "secondary", ZOOPLANKTON = "secondary"),
  `12S_MT` = c(METAZOA = "primary", FISH = "primary", BIRDS = "primary"),
  ITS_FUNGAL = c(FUNGI = "primary"),
  `18S` = c(METAZOA = "secondary", FUNGI = "secondary", PROTISTS = "primary", PHYTOPLANKTON = "primary", ZOOPLANKTON = "primary"),
  `16S_PROK` = c(BACTERIA = "primary", ARCHAEA = "primary", PLANT_MICROBIOME = "primary"),
  `16S_MT` = c(METAZOA = "primary", ARTHROPODS = "primary", FISH = "primary", ZOOPLANKTON = "secondary"),
  `28S` = c(FUNGI = "primary", PROTISTS = "primary", PHYTOPLANKTON = "primary", METAZOA = "secondary")
)
marker_organism_groups <- do.call(rbind, lapply(names(marker_group_rows), function(marker_id) {
  relationships <- marker_group_rows[[marker_id]]
  data.frame(
    marker_id = marker_id,
    group_id = names(relationships),
    relationship = unname(relationships),
    evidence_note = "Established or complementary use documented in the 2016–2026 primer compendium; marker membership does not imply uniform primer coverage.",
    source_url = "https://github.com/sven9r/primer-atlas/blob/main/docs/marker-targets.md",
    stringsAsFactors = FALSE
  )
}))
rownames(marker_organism_groups) <- NULL
write_catalog(marker_organism_groups, "marker_organism_groups.csv")

landmarks <- data.frame(
  marker_id = c(
    "COI", "COI",
    rep("ITS_FUNGAL", 5),
    rep("18S", 9), rep("16S_PROK", 9), rep("28S", 4)
  ),
  landmark_id = c(
    "COX1", "FOLMER",
    "18S_ANCHOR", "ITS1", "5_8S", "ITS2", "28S_ANCHOR",
    paste0("V", 1:9), paste0("V", 1:9), paste0("D", 1:4)
  ),
  label = c(
    "COX1", "Folmer region", "18S", "ITS1", "5.8S", "ITS2", "28S",
    paste0("V", 1:9), paste0("V", 1:9), paste0("D", 1:4)
  ),
  start = c(
    1L, 42L, 1L, 26L, 339L, 501L, 719L,
    seq(1L, 1601L, length.out = 9),
    seq(1L, 1401L, length.out = 9),
    c(1L, 401L, 801L, 1201L)
  ),
  end = c(
    1536L, 699L, 25L, 338L, 500L, 718L, 739L,
    seq(200L, 1800L, length.out = 9),
    seq(175L, 1575L, length.out = 9),
    c(400L, 800L, 1200L, 1600L)
  ),
  landmark_type = c(
    "gene", "focus_region", "conserved_anchor", "variable_region",
    "conserved_gene", "variable_region", "conserved_anchor",
    rep("variable_region", 18), rep("domain", 4)
  ),
  display_order = sequence(c(2, 5, 9, 9, 4)),
  source_note = c(
    "NC_001322.1 COX1 coordinates", "Folmer barcode reference region",
    rep("FN812768.2 feature architecture; boundaries curated for the pilot display", 5),
    rep("Planned release-specific landmark alignment", 22)
  ),
  stringsAsFactors = FALSE
)
landmarks$start <- as.integer(round(landmarks$start))
landmarks$end <- as.integer(round(landmarks$end))
write_catalog(landmarks, "marker_landmarks.csv")

legacy <- read.csv(
  file.path(project_root, "data", "primers.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
legacy$marker_id <- "COI"

oligo_key <- paste(legacy$marker_id, legacy$primer_name, legacy$sequence, sep = "|")
unique_oligo <- !duplicated(oligo_key)
oligos <- data.frame(
  oligo_id = paste0(
    "COI_", toupper(gsub("[^A-Za-z0-9]+", "_", legacy$primer_name[unique_oligo]))
  ),
  marker_id = "COI",
  primer_name = legacy$primer_name[unique_oligo],
  sequence = legacy$sequence[unique_oligo],
  source_key = legacy$source_key[unique_oligo],
  source_note = legacy$source_note[unique_oligo],
  stringsAsFactors = FALSE
)
oligo_id_by_key <- setNames(oligos$oligo_id, oligo_key[unique_oligo])

pair_oligos <- data.frame(
  pair_id = legacy$pair_id,
  oligo_id = unname(oligo_id_by_key[oligo_key]),
  direction = legacy$direction,
  stringsAsFactors = FALSE
)

pair_groups <- split(legacy, legacy$pair_id)
pairs <- do.call(rbind, lapply(pair_groups, function(x) {
  data.frame(
    pair_id = x$pair_id[1], pair_label = x$pair_label[1], marker_id = "COI",
    reported_amplicon_bp = suppressWarnings(max(x$reported_amplicon_bp, na.rm = TRUE)),
    original_use_case = paste(unique(x$use_case), collapse = " | "),
    target_summary = paste(unique(x$target), collapse = " / "),
    reported_ta_c = if (all(is.na(x$reported_ta_c))) NA else min(x$reported_ta_c, na.rm = TRUE),
    catalog_status = "curated",
    stringsAsFactors = FALSE
  )
}))
rownames(pairs) <- NULL

applications_for <- function(pair_id, text) {
  text <- tolower(text)
  values <- character()
  if (grepl("barcoding", text)) values <- c(values, "barcoding")
  if (grepl("edna|water", text)) values <- c(values, "edna")
  if (grepl("diet|gut content", text)) values <- c(values, "diet")
  if (grepl("bulk|metabarcoding|plant-associated", text) && !grepl("diet|gut content", text)) {
    values <- c(values, "bulk_community")
  }
  if (!length(values)) values <- "bulk_community"
  unique(values)
}
environment_for <- function(text) {
  text <- tolower(text)
  values <- character()
  if (grepl("freshwater|water edna", text)) values <- c(values, "freshwater")
  if (grepl("marine|mollusk", text)) values <- c(values, "marine")
  if (grepl("terrestrial|arthropod|bee|insect|spider|plant", text)) values <- c(values, "terrestrial")
  if (grepl("diet|gut content|plant-associated", text)) values <- c(values, "host_associated")
  unique(values)
}
broad_pairs <- c("FOLMER", "MCO", "LERAY_XT", "ANML", "ARF")
exclusion_pairs <- c("NOSPID", "NOSPI2_LAURELIN", "NOPLANT")

facets <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  row <- pairs[i, ]
  text <- paste(row$original_use_case, row$target_summary)
  applications <- applications_for(row$pair_id, text)
  environments <- environment_for(text)
  design <- if (row$pair_id %in% exclusion_pairs) {
    "exclusion_blocking"
  } else if (row$pair_id %in% broad_pairs) {
    "broad"
  } else {
    "target_enriched"
  }
  rbind(
    data.frame(pair_id = row$pair_id, facet_type = "application", facet_value = applications),
    data.frame(pair_id = row$pair_id, facet_type = "target_taxon", facet_value = row$target_summary),
    if (length(environments)) data.frame(pair_id = row$pair_id, facet_type = "environment", facet_value = environments),
    data.frame(pair_id = row$pair_id, facet_type = "design_intent", facet_value = design)
  )
}))
rownames(facets) <- NULL

# Curated fungal ITS pilot pairs from the UNITE resource. Oligo metadata is
# imported separately and coordinates are resolved against FN812768.2.
unite_path <- file.path(catalog_dir, "unite_primers.csv")
if (file.exists(unite_path)) {
  unite <- read.csv(unite_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"primary_reference_key" %in% names(unite)) {
    ref_slug <- function(x) paste0("its_primary_", gsub("(^_|_$)", "", gsub("[^a-z0-9]+", "_", tolower(ifelse(nzchar(x), x, "not_reported")))))
    unite$primary_reference_key <- ref_slug(unite$reference)
    unite$primary_reference_url <- ifelse(
      !nzchar(unite$reference) | grepl("unpublished", unite$reference, ignore.case = TRUE),
      "https://unite.ut.ee/primers.php",
      paste0("https://search.crossref.org/?q=", utils::URLencode(unite$reference, reserved = TRUE))
    )
    unite$primary_reference_status <- ifelse(!nzchar(unite$reference), "not_reported", ifelse(grepl("unpublished", unite$reference, ignore.case = TRUE), "unpublished", "publication_cited"))
  }
  its_pair_specs <- data.frame(
    pair_id = c(
      "ITS1_ITS4", "ITS1F_ITS4", "ITS1F_ITS4B", "ITS3_ITS4",
      "FITS7_ITS4", "GITS7_ITS4", "ITS86F_ITS4", "ITS1F_LR21"
    ),
    pair_label = c(
      "ITS1 + ITS4", "ITS1F + ITS4", "ITS1F + ITS4B", "ITS3 + ITS4",
      "fITS7 + ITS4", "gITS7 + ITS4", "ITS86F + ITS4", "ITS1F + LR21"
    ),
    forward = c("ITS1", "ITS1F", "ITS1F", "ITS3", "fITS7", "gITS7", "ITS86F", "ITS1F"),
    reverse = c("ITS4", "ITS4", "ITS4B", "ITS4", "ITS4", "ITS4", "ITS4", "LR21"),
    applications = c(
      "barcoding|bulk_community", "barcoding|bulk_community", "barcoding|bulk_community",
      "barcoding", "bulk_community|edna", "bulk_community|edna",
      "bulk_community|edna", "barcoding"
    ),
    design_intent = c("broad", "target_enriched", "target_enriched", "broad", "target_enriched", "target_enriched", "target_enriched", "target_enriched"),
    stringsAsFactors = FALSE
  )
  select_unite <- function(name, direction) {
    hit <- unite[unite$primer_name == name & unite$direction == direction, , drop = FALSE]
    if (!nrow(hit)) stop("UNITE primer missing from snapshot: ", name)
    hit[1, , drop = FALSE]
  }
  its_oligos <- do.call(rbind, lapply(unique(c(its_pair_specs$forward, its_pair_specs$reverse)), function(name) {
    direction <- if (name %in% its_pair_specs$forward) "forward" else "reverse"
    hit <- select_unite(name, direction)
    data.frame(
      oligo_id = paste0("ITS_", toupper(gsub("[^A-Za-z0-9]+", "_", name))),
      marker_id = "ITS_FUNGAL", primer_name = name, sequence = hit$sequence,
      source_key = hit$primary_reference_key,
      source_note = paste(hit$remarks, "Original citation:", ifelse(nzchar(hit$reference), hit$reference, "not reported"), "UNITE compilation snapshot", sep = " · "),
      stringsAsFactors = FALSE
    )
  }))
  oligos <- rbind(oligos, its_oligos)
  its_pair_oligos <- do.call(rbind, lapply(seq_len(nrow(its_pair_specs)), function(i) {
    x <- its_pair_specs[i, ]
    data.frame(
      pair_id = x$pair_id,
      oligo_id = c(
        paste0("ITS_", toupper(gsub("[^A-Za-z0-9]+", "_", x$forward))),
        paste0("ITS_", toupper(gsub("[^A-Za-z0-9]+", "_", x$reverse)))
      ),
      direction = c("forward", "reverse"), stringsAsFactors = FALSE
    )
  }))
  pair_oligos <- rbind(pair_oligos, its_pair_oligos)
  its_pairs <- data.frame(
    pair_id = its_pair_specs$pair_id, pair_label = its_pair_specs$pair_label,
    marker_id = "ITS_FUNGAL", reported_amplicon_bp = NA,
    original_use_case = "fungal ITS barcoding / metabarcoding",
    target_summary = "Fungi", reported_ta_c = 55,
    catalog_status = "pilot", stringsAsFactors = FALSE
  )
  pairs <- rbind(pairs, its_pairs)
  its_facets <- do.call(rbind, lapply(seq_len(nrow(its_pair_specs)), function(i) {
    x <- its_pair_specs[i, ]
    rbind(
      data.frame(pair_id = x$pair_id, facet_type = "application", facet_value = strsplit(x$applications, "|", fixed = TRUE)[[1]]),
      data.frame(pair_id = x$pair_id, facet_type = "target_taxon", facet_value = "Fungi"),
      data.frame(pair_id = x$pair_id, facet_type = "environment", facet_value = c("terrestrial", "host_associated")),
      data.frame(pair_id = x$pair_id, facet_type = "design_intent", facet_value = x$design_intent)
    )
  }))
  facets <- rbind(facets, its_facets)
}

claims <- read.csv(
  file.path(project_root, "data", "primer_target_claims.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
claims$marker_id <- "COI"
claims$claim_id <- paste0("CLAIM_", seq_len(nrow(claims)))
claims$evidence_type <- "publication_claim"
claims$curator_note <- claims$claim_note

if (exists("its_pair_specs")) {
  its_claims <- do.call(rbind, lapply(seq_len(nrow(its_pair_specs)), function(i) {
    spec <- its_pair_specs[i, ]
    selected <- unite[unite$primer_name %in% c(spec$forward, spec$reverse), , drop = FALSE]
    wording <- paste(unique(selected$remarks[nzchar(selected$remarks)]), collapse = " | ")
    references <- paste(unique(selected$reference[nzchar(selected$reference)]), collapse = " | ")
    data.frame(
      pair_id = spec$pair_id,
      claim_label = paste("UNITE fungal ITS catalog entry:", spec$pair_label),
      claim_kind = "source_catalog_entry", target_rank = "kingdom",
      target_taxa = "Fungi", contrast_taxa = "",
      pair_validation_status = "source_catalog_import",
      expanded_reference_role = "screening_only",
      source_keys = paste(unique(c(selected$primary_reference_key, "unite_primers")), collapse = "|"),
      claim_note = paste(wording, references, sep = " · "),
      marker_id = "ITS_FUNGAL", claim_id = paste0("ITS_CLAIM_", i),
      evidence_type = "attributed_database_snapshot",
      curator_note = "Imported from the versioned UNITE primer snapshot; recommendation and exclusion wording is preserved, not reinterpreted.",
      stringsAsFactors = FALSE
    )
  }))
  claims <- rbind(claims, its_claims)
}

write_catalog(oligos, "oligos.csv")
write_catalog(unique(pair_oligos), "pair_oligos.csv")
write_catalog(pairs, "primer_pairs.csv")
write_catalog(unique(facets), "primer_pair_facets.csv")
write_catalog(claims, "claims.csv")
sources <- read.csv(file.path(project_root, "data", "citations.csv"), stringsAsFactors = FALSE)
if (exists("unite")) {
  primary_sources <- unique(unite[, c(
    "primary_reference_key", "reference", "primary_reference_url",
    "primary_reference_status"
  )])
  primary_sources <- data.frame(
    key = primary_sources$primary_reference_key,
    category = "original primer reference",
    short_citation = ifelse(nzchar(primary_sources$reference), primary_sources$reference, "Original reference not reported"),
    title = ifelse(nzchar(primary_sources$reference), paste("Primer source:", primary_sources$reference), "Original primer reference not reported by compilation"),
    year = suppressWarnings(as.integer(sub(".*?([12][0-9]{3}).*", "\\1", primary_sources$reference))),
    doi = "", url = primary_sources$primary_reference_url,
    note = paste("Primary-reference status:", primary_sources$primary_reference_status, "· primer metadata compiled by UNITE."),
    stringsAsFactors = FALSE
  )
  primary_sources$year[is.na(primary_sources$year) | primary_sources$year < 1900L | primary_sources$year > 2100L] <- NA_integer_
  sources <- rbind(sources[!sources$key %in% primary_sources$key, , drop = FALSE], primary_sources)
}
if (!"unite_primers" %in% sources$key) {
  sources <- rbind(sources, data.frame(
    key = "unite_primers", category = "primer catalog",
    short_citation = "UNITE primer resource",
    title = "UNITE fungal ITS primer resource", year = 2026,
    doi = "", url = "https://unite.ut.ee/primers.php",
    note = "Versioned attributed snapshot imported offline; not scraped during user sessions.",
    stringsAsFactors = FALSE
  ))
}
write_catalog(sources, "sources.csv")

message(
  "Catalog build complete: ", nrow(pairs), " pairs across ",
  length(unique(pairs$marker_id)), " populated markers."
)
