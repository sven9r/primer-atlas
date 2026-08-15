local_r_library <- file.path(getwd(), ".Rlib")
if (dir.exists(local_r_library)) {
  .libPaths(c(local_r_library, .libPaths()))
}

library(shiny)
library(bslib)
library(readr)
library(dplyr)
library(stringr)
library(DT)

source("R/functions.R")
source("R/data_layer.R")

coi_release_state <- tryCatch(atlas_read_manifest("COI"), error = function(e) NULL)
release_table <- function(name, fallback, col_types = NULL) {
  artifact <- if (!is.null(coi_release_state)) coi_release_state$manifest$artifacts[[name]] else NULL
  if (!is.null(artifact)) {
    value <- tryCatch(as_tibble(atlas_read_artifact(artifact)), error = function(e) NULL)
    if (!is.null(value)) return(value)
  }
  read_csv(fallback, col_types = col_types, show_col_types = FALSE)
}

primers <- read_csv("data/primers.csv", show_col_types = FALSE)
citations <- read_csv("data/citations.csv", show_col_types = FALSE)
marker_registry <- read_csv("data/catalog/markers.csv", show_col_types = FALSE)
marker_landmarks <- read_csv("data/catalog/marker_landmarks.csv", show_col_types = FALSE)
catalog_pairs <- read_csv("data/catalog/primer_pairs.csv", show_col_types = FALSE)
catalog_oligos <- read_csv("data/catalog/oligos.csv", show_col_types = FALSE)
catalog_links <- read_csv("data/catalog/pair_oligos.csv", show_col_types = FALSE)
primer_facets <- read_csv("data/catalog/primer_pair_facets.csv", show_col_types = FALSE)
catalog_claims <- read_csv("data/catalog/claims.csv", show_col_types = FALSE)
its_primer_rows <- catalog_links |>
  inner_join(catalog_oligos, by = "oligo_id") |>
  inner_join(catalog_pairs |> select(pair_id, pair_label, reported_amplicon_bp,
                                     original_use_case, target_summary, reported_ta_c), by = "pair_id") |>
  filter(marker_id == "ITS_FUNGAL") |>
  transmute(
    pair_id, pair_label, primer_name, direction, sequence,
    reported_amplicon_bp, use_case = original_use_case, target = target_summary,
    reported_ta_c, source_key, source_note
  )
primers <- bind_rows(primers, its_primer_rows)
if (!"unite_primers" %in% citations$key) {
  citations <- bind_rows(citations, tibble(
    key = "unite_primers", category = "primer catalog",
    short_citation = "UNITE primer resource",
    title = "UNITE fungal ITS primer resource", year = 2026,
    doi = NA_character_, url = "https://unite.ut.ee/primers.php",
    note = "Versioned attributed snapshot imported offline; not scraped during user sessions."
  ))
}
primer_bindings <- read_csv("data/derived/primer_bindings.csv", show_col_types = FALSE)
pair_meta <- read_csv("data/derived/marker_pair_geometry.csv", show_col_types = FALSE) |>
  arrange(marker_id, folmer_overlap_rank)
folmer_region <- read_csv("data/derived/folmer_region.csv", show_col_types = FALSE)
coi_reference <- clean_sequence(
  parse_fasta("data/reference/coi_reference_NC_001322.1.fasta")[[1]]
)
its_reference <- clean_sequence(
  parse_fasta("data/reference/its_fungal_reference_FN812768.2.fasta")[[1]]
)

alignment_files <- list.files(
  "data/primerminer/aligned_reference",
  pattern = "_COX1_reference_aligned\\.fasta$",
  full.names = TRUE
)
alignment_orders <- sub(
  "_COX1_reference_aligned\\.fasta$",
  "",
  basename(alignment_files)
)
order_alignment_templates <- setNames(
  lapply(alignment_files, function(path) {
    sequences <- parse_fasta(path)
    sequences[!grepl("^NC_001322\\.1_COX1_reference", names(sequences))]
  }),
  alignment_orders
)
order_scores <- read_csv("data/derived/order_pair_summary.csv", show_col_types = FALSE)
pair_template_scores <- release_table(
  "pair_templates", "data/derived/pair_template_scores.csv",
  col_types = cols(template = col_character())
)
delayedAssign("primer_position_scores", {
  release_table(
    "primer_positions", "data/derived/primer_position_scores.csv",
    col_types = cols(template = col_character())
  ) |>
    mutate(unavailable_reason = coalesce(unavailable_reason, ""))
})

beeprime_reference_manifest <- read_csv(
  "data/derived/beeprime_reference_manifest.csv",
  show_col_types = FALSE
)
expanded_reference_alignment_path <- file.path(
  "data", "external", "gurten2026", "ClusteredReferences.fasta"
)
expanded_reference_taxonomy <- release_table(
  "taxonomy", "data/derived/gurten2026_centroid_taxonomy.csv",
  col_types = cols(accession = col_character())
)
manifest_value <- function(field) {
  beeprime_reference_manifest$value[
    match(field, beeprime_reference_manifest$field)
  ]
}
column_start <- function(value) as.integer(sub("-.*$", "", value))
column_end <- function(value) as.integer(sub("^.*-", "", value))
beeprime_forward_binding <- primer_bindings |>
  filter(pair_id == "BEEPRIME", direction == "forward")
beeprime_reverse_binding <- primer_bindings |>
  filter(pair_id == "BEEPRIME", direction == "reverse")
expanded_reference_offsets <- c(
  column_start(manifest_value("beeprime_forward_alignment_columns")) -
    beeprime_forward_binding$alignment_start,
  column_end(manifest_value("beeprime_forward_alignment_columns")) -
    beeprime_forward_binding$alignment_end,
  column_start(manifest_value("beeprime_reverse_alignment_columns")) -
    beeprime_reverse_binding$alignment_start,
  column_end(manifest_value("beeprime_reverse_alignment_columns")) -
    beeprime_reverse_binding$alignment_end
)
stopifnot(length(unique(expanded_reference_offsets)) == 1L)
expanded_reference_alignment_offset <- unique(expanded_reference_offsets)
beeprime_bee_family <- read_csv(
  "data/derived/beeprime_bee_family_summary.csv",
  show_col_types = FALSE
)
beeprime_bee_subfamily <- read_csv(
  "data/derived/beeprime_bee_subfamily_summary.csv",
  show_col_types = FALSE
)
beeprime_bee_genus <- read_csv(
  "data/derived/beeprime_bee_genus_summary.csv",
  show_col_types = FALSE
)
beeprime_hymenoptera_family <- read_csv(
  "data/derived/beeprime_hymenoptera_family_summary.csv",
  show_col_types = FALSE
)
beeprime_hymenoptera_lineage <- read_csv(
  "data/derived/beeprime_hymenoptera_lineage_summary.csv",
  show_col_types = FALSE
)
beeprime_hymenoptera_scores <- read_csv(
  "data/derived/beeprime_hymenoptera_centroid_scores.csv",
  col_types = cols(accession = col_character()),
  show_col_types = FALSE
)
beeprime_empirical <- read_csv(
  "data/derived/beeprime_empirical_validation.csv",
  show_col_types = FALSE
)
beeprime_empirical_genus <- read_csv(
  "data/derived/beeprime_empirical_genus_summary.csv",
  show_col_types = FALSE
)
primer_target_claims <- read_csv(
  "data/primer_target_claims.csv",
  show_col_types = FALSE
)
claimed_primer_overall <- read_csv(
  "data/derived/claimed_primer_overall_summary.csv",
  show_col_types = FALSE
)
claimed_primer_orders <- read_csv(
  "data/derived/claimed_primer_order_summary.csv",
  show_col_types = FALSE
)
claimed_primer_families <- release_table("families", "data/derived/claimed_primer_family_summary.csv")
claimed_primer_genera <- release_table("genera", "data/derived/claimed_primer_genus_summary.csv")
claimed_primer_subfamilies <- release_table("subfamilies", "data/derived/claimed_primer_subfamily_summary.csv")
claimed_primer_phylogeny_groups <- read_csv(
  "data/derived/claimed_primer_phylogeny_group_summary.csv",
  show_col_types = FALSE
)
targeted_zbj_manifest <- read_csv(
  "data/external/targeted_zbj/targeted_complete_cox1_manifest.csv",
  show_col_types = FALSE
)
targeted_zbj_overall <- read_csv(
  "data/derived/targeted_zbj_overall_summary.csv",
  show_col_types = FALSE
)
targeted_zbj_families <- read_csv(
  "data/derived/targeted_zbj_family_summary.csv",
  show_col_types = FALSE
)
targeted_zbj_genera <- read_csv(
  "data/derived/targeted_zbj_genus_summary.csv",
  show_col_types = FALSE
)
targeted_zbj_subfamilies <- read_csv(
  "data/derived/targeted_zbj_subfamily_summary.csv",
  show_col_types = FALSE
)
targeted_zbj_study_groups <- read_csv(
  "data/derived/targeted_zbj_study_group_summary.csv",
  show_col_types = FALSE
)
targeted_zbj_sequence_scores <- read_csv(
  "data/derived/targeted_zbj_sequence_scores.csv.gz",
  show_col_types = FALSE
)

citation_link <- function(row, label = NULL) {
  if (!nrow(row) || is.na(row$url[1]) || !nzchar(row$url[1])) {
    return(span(class = "tiny", "Source link pending"))
  }
  tags$a(
    href = row$url[1],
    target = "_blank",
    if (is.null(label)) row$short_citation[1] else label
  )
}

atlas_theme <- bs_theme(
  version = 5,
  bg = "#f4f2ec",
  fg = "#17241d",
  primary = "#16664b",
  secondary = "#d8e7df",
  base_font = font_google("Manrope"),
  heading_font = font_google("DM Serif Display")
)

app_css <- "
:root { --ink:#17241d; --green:#16664b; --mint:#dcebe3; --paper:#f4f2ec;
  --coral:#de6b4f; --gold:#d7a62b; --line:#d5d3ca; }
body { letter-spacing:-.01em; }
.navbar { border-bottom:1px solid var(--line); background:rgba(244,242,236,.96)!important; }
.navbar > .container-fluid { flex-wrap:nowrap; }
.navbar-nav { flex-wrap:nowrap; overflow-x:auto; }
.navbar-brand { font-family:'DM Serif Display',serif; font-size:1.08rem; white-space:nowrap; }
.institution-banner { display:flex; align-items:center; justify-content:space-between;
  gap:24px; padding:9px 24px; background:#fffdf8; border-bottom:1px solid var(--line); }
.banner-home { display:flex; align-items:center; gap:14px; min-width:250px; }
.banner-home img { width:108px; height:68px; object-fit:contain; }
.banner-home strong { display:block; font-family:'DM Serif Display',serif;
  font-size:1.2rem; line-height:1.05; }
.banner-home span { display:block; margin-top:3px; color:#65716a; font-size:.72rem; }
.banner-support { display:flex; align-items:center; justify-content:flex-end;
  gap:20px; min-width:0; }
.banner-support-label { color:#65716a; font-size:.68rem; text-transform:uppercase;
  letter-spacing:.1em; white-space:nowrap; }
.banner-logo { object-fit:contain; display:block; }
.banner-logo.fulbright { width:112px; height:48px; }
.banner-logo.walter { width:128px; height:58px; }
.banner-logo.berkeley { width:66px; height:66px; }
.hero { padding:.85rem 0 .6rem; }
.eyebrow { color:var(--green); text-transform:uppercase; letter-spacing:.13em;
  font-weight:800; font-size:.72rem; }
.hero h1 { font-size:clamp(2rem,3.2vw,3.15rem); line-height:1; max-width:none;
  white-space:nowrap; margin-bottom:.35rem; }
.hero p { max-width:820px; font-size:.96rem; color:#536159; margin-bottom:.25rem; }
.bslib-card { border:1px solid var(--line); box-shadow:none; border-radius:18px; overflow:hidden; }
.card-header { background:#fbfaf6; border-bottom:1px solid var(--line); font-weight:800; }
.sidebar { background:#eceae2!important; border-right:1px solid var(--line)!important;
  position:relative; z-index:3; }
.bslib-sidebar-layout > .main { min-width:0; overflow:hidden; }
.metric-row { display:grid; grid-template-columns:repeat(auto-fit,minmax(145px,1fr)); gap:12px; margin-bottom:16px; }
.metric { background:#fffdf8; border:1px solid var(--line); border-radius:16px; padding:14px 16px; }
.metric b { display:block; font-family:'DM Serif Display'; font-size:2rem; line-height:1; }
.metric span { font-size:.75rem; text-transform:uppercase; letter-spacing:.08em; color:#68736c; }
.primer-map-wrap { overflow-x:auto; max-width:100%; padding:8px 2px 16px;
  position:relative; z-index:1; }
.primer-map { min-width:1350px; width:100%; }
.primer-tooltip-target { cursor:pointer; }
.primer-tooltip-target:focus-visible {
  outline:3px solid rgba(22,102,75,.55); outline-offset:4px;
}
.primer-hover-tooltip {
  display:none; position:fixed; z-index:10000;
  width:min(440px,calc(100vw - 28px)); max-height:min(420px,calc(100vh - 28px));
  overflow:auto; padding:15px 17px 12px;
  border:2px solid #6f7b74; border-radius:14px;
  background:#fffdf8; color:#17241d;
  box-shadow:0 14px 38px rgba(23,36,29,.24);
  font-size:.9rem; line-height:1.48; font-weight:600;
  white-space:pre-line; pointer-events:none;
}
.primer-hover-tooltip.visible { display:block; }
.primer-hover-tooltip.pinned {
  border-color:var(--green); box-shadow:0 16px 42px rgba(22,102,75,.28);
}
.primer-hover-tooltip::after {
  content:'Hover preview · click primer name or amplicon to pin';
  display:block; margin-top:10px; padding-top:8px;
  border-top:1px solid var(--line); color:#65716a;
  font-size:.7rem; font-weight:700; white-space:normal;
}
.primer-hover-tooltip.pinned::after {
  content:'Pinned · click the same item again or press Escape to close';
  color:var(--green);
}
.map-legend { display:flex; flex-wrap:wrap; gap:14px; align-items:center;
  padding:2px 2px 10px; font-size:.75rem; color:#59655e; }
.legend-line { display:inline-block; width:28px; height:5px; border-radius:999px;
  margin-right:6px; vertical-align:middle; }
.legend-forward { background:#2563eb; }
.legend-informative { background:#9aa6b2; }
.legend-reverse { background:#dc2626; }
.legend-number {
  display:inline-flex; width:24px; height:24px; border-radius:50%;
  align-items:center; justify-content:center; margin-right:6px;
  background:#17241d; color:#fffdf8; font-size:.7rem; font-weight:800;
  vertical-align:middle;
}
.pair-card { display:grid; grid-template-columns:1.2fr 1fr auto; gap:14px; align-items:center;
  border-top:1px solid var(--line); padding:12px 4px; }
.pair-card:first-child { border-top:0; }
.pill { display:inline-block; border-radius:999px; padding:4px 9px; background:var(--mint);
  color:var(--green); font-size:.72rem; font-weight:800; margin:2px 3px 2px 0; }
.pill.coral { background:#f7dfd8; color:#9b3d29; }
.mono { font-family:ui-monospace,SFMono-Regular,Menlo,monospace; font-size:.78rem; }
.heatmap-scroll-shell { width:100%; min-width:0; position:relative; }
.heatmap-scroll-top, .heatmap-scroll-body {
  width:100%; max-width:100%; overflow-x:auto; overflow-y:hidden;
  scrollbar-gutter:stable;
}
.heatmap-scroll-top {
  height:18px; position:sticky; top:0; z-index:4;
  background:var(--paper); border-bottom:1px solid var(--line);
}
.heatmap-scroll-top-inner { height:1px; }
.heatmap { display:grid; gap:4px; align-items:stretch; min-width:100%; overflow:visible; }
.heat-cell { min-width:82px; border:0; border-radius:8px; padding:10px 6px; text-align:center;
  color:#183027; font-weight:800; cursor:pointer; font-family:inherit; }
.heat-cell:hover { box-shadow:inset 0 0 0 2px rgba(22,102,75,.55); }
.heat-cell:focus-visible { outline:3px solid rgba(22,102,75,.55); outline-offset:1px; }
.heat-cell-static { cursor:default; }
.heat-head { font-size:.68rem; writing-mode:vertical-rl; transform:rotate(180deg);
  min-height:104px; display:flex; justify-content:flex-start; align-items:center; padding:6px; }
.order-label { font-weight:800; padding:10px 8px; white-space:nowrap; }
.risk-na { background:repeating-linear-gradient(135deg,#e5e3dc,#e5e3dc 7px,#d7d4cb 7px,#d7d4cb 14px); }
.heatmap-reading-guide { display:grid; grid-template-columns:minmax(0,1.2fr) minmax(300px,.8fr);
  gap:16px; margin:0 0 18px; padding:16px; border:1px solid var(--line);
  border-radius:14px; background:#fffdf8; }
.heatmap-reading-guide h5 { margin:0 0 8px; }
.heatmap-reading-guide p { margin:0 0 7px; font-size:.82rem; color:#48554e; }
.heatmap-legend-grid { display:grid; grid-template-columns:repeat(4,minmax(115px,1fr));
  gap:9px; margin-top:10px; }
.heatmap-legend-item { display:flex; align-items:center; gap:8px; min-width:0;
  font-size:.75rem; color:#48554e; }
.heatmap-swatch { width:34px; height:25px; flex:0 0 34px; border-radius:6px;
  border:1px solid rgba(23,36,29,.12); }
.heatmap-swatch.low { background:#b9ddc9; }
.heatmap-swatch.mid { background:#efd787; }
.heatmap-swatch.high { background:#efaa90; }
.heatmap-swatch.missing {
  background:repeating-linear-gradient(135deg,#e5e3dc,#e5e3dc 7px,#d7d4cb 7px,#d7d4cb 14px);
}
.investigation-tip { border-left:4px solid var(--green); padding-left:13px; }
.investigation-tip ol { margin:7px 0 0; padding-left:20px; color:#48554e; font-size:.8rem; }
.investigation-tip li { margin-bottom:5px; }
.order-heatmap-card { overflow:visible; }
.order-heatmap-card > .card-header { border-radius:17px 17px 0 0; }
.order-heatmap-card > .card-body { overflow:visible!important; }
.order-detail-layout, .order-detail-layout > * { min-width:0; }
.order-detail-layout {
  align-items:start; overflow:visible!important; position:relative;
}
.order-detail-card, .deep-detail-card { min-width:0; }
.order-detail-card { align-self:start; height:auto!important; }
.order-inspector-card {
  overflow:visible!important; position:relative; z-index:30;
}
.order-inspector-card > .card-body {
  overflow:visible!important;
}
.order-inspector-card .selectize-control {
  position:relative; z-index:40;
}
.order-inspector-card .selectize-control.dropdown-active {
  z-index:60;
}
.order-inspector-card .selectize-dropdown {
  z-index:70; min-width:100%!important;
  box-shadow:0 10px 28px rgba(23,36,29,.18);
}
.order-inspector-card .selectize-dropdown-content {
  max-height:360px!important; overflow-y:auto!important;
  overscroll-behavior:contain;
}
.order-inspector-controls {
  display:grid; grid-template-columns:repeat(3,minmax(180px,1fr));
  gap:12px; align-items:end;
}
.order-inspector-actions {
  display:flex; flex-wrap:wrap; gap:12px; align-items:center; margin-top:2px;
}
.order-inspector-actions .tiny { margin:0; }
.order-detail-card > .card-body, .deep-detail-card > .card-body {
  overflow:visible!important; min-width:0;
}
.comparison-trace-shell {
  min-height:260px; position:relative;
}
.comparison-trace-shell > .shiny-bound-output {
  min-height:260px;
}
body.detail-switching .comparison-trace-shell::before {
  content:'Loading selected comparison…';
  position:absolute; z-index:5; top:12px; left:12px;
  border-radius:999px; padding:7px 11px;
  background:#fffdf8; border:1px solid var(--line);
  color:var(--green); font-size:.75rem; font-weight:800;
  box-shadow:0 3px 10px rgba(23,36,29,.08);
}
body.detail-switching .comparison-trace-shell > .shiny-bound-output {
  opacity:.38; transition:opacity .15s ease;
}
body.detail-switching .deep-detail-card {
  opacity:.55; pointer-events:none; cursor:progress;
  transition:opacity .15s ease;
}
.detail-table-scroll { width:100%; max-width:100%; overflow:auto; }
.deep-detail-card .datatables.html-widget-output {
  flex:0 0 auto!important; height:auto!important; min-height:0!important;
}
.deep-detail-card .dataTables_wrapper { width:100%; max-width:100%; overflow:hidden; }
.deep-detail-card .dataTables_scroll { width:100%; max-width:100%; }
.deep-detail-card .dataTables_scrollBody { overflow:auto!important; }
.sequence-strip { display:flex; flex-wrap:wrap; gap:3px; margin:10px 0 18px; }
.base { width:28px; height:42px; border-radius:6px; background:#e6e4dc; text-align:center;
  font-family:ui-monospace,monospace; font-weight:800; padding-top:4px; position:relative; }
.base.mismatch { background:#f1b09b; color:#7c2516; }
.base.unavailable { background:#d8d6d0; color:#59615c; }
.base small { display:block; font-size:.55rem; color:#65716a; margin-top:2px; }
.table-atlas { width:100%; border-collapse:separate; border-spacing:0; font-size:.84rem; }
.table-atlas th { position:sticky; top:0; background:#eae8e0; text-align:left; padding:9px; }
.table-atlas td { padding:9px; border-top:1px solid var(--line); vertical-align:top; }
.cite-card { border-top:1px solid var(--line); padding:15px 2px; }
.cite-card:first-child { border-top:0; }
.callout { border-left:4px solid var(--gold); background:#fbf4da; padding:12px 15px;
  border-radius:0 12px 12px 0; }
.taxon-info { border-left-color:#2563eb; background:#eef5ff; }
.taxon-alert { border-left-color:#dc2626; background:#fff0ed; }
.deep-dive-controls { display:grid; grid-template-columns:minmax(220px,1fr) minmax(180px,.7fr);
  gap:12px; align-items:end; margin:4px 0 12px; }
.evidence-badge { display:inline-block; border-radius:999px; padding:2px 7px;
  background:#e7e5de; color:#58635c; font-size:.68rem; font-weight:800; white-space:nowrap; }
.evidence-badge.small-n { background:#fbefd0; color:#8a5a00; }
.evidence-badge.empirical { background:#dcebe3; color:#16664b; }
.scroll-table { overflow:auto; max-height:620px; border:1px solid var(--line);
  border-radius:12px; }
.dataTables_wrapper { font-size:.8rem; }
.dataTables_wrapper .dataTables_filter input,
.dataTables_wrapper .dataTables_length select,
.dataTables_wrapper thead input {
  border:1px solid var(--line); border-radius:7px; background:#fffdf8;
  color:var(--ink); padding:5px 7px;
}
.dataTables_wrapper thead input { width:100%!important; min-width:76px; }
.dataTables_wrapper table.dataTable tbody tr.selected > * {
  box-shadow:inset 0 0 0 9999px rgba(22,102,75,.13)!important;
  color:var(--ink)!important;
}
.sequence-source-table table.dataTable tbody tr { cursor:pointer; }
.sequence-source-table .dataTables_empty { cursor:default; }
.sequence-source-table .recalculating,
.sequence-source-table .shiny-recalculating {
  opacity:.55; pointer-events:none;
}
.sequence-evidence-panel { margin-top:18px; padding-top:16px;
  border-top:1px solid var(--line); }
.sequence-evidence-actions { display:flex; flex-wrap:wrap; gap:8px;
  margin:8px 0 12px; }
.map-clear-btn { color:var(--green)!important; border-color:var(--green)!important;
  background:transparent!important; }
.map-clear-btn:hover { color:#fffdf8!important; background:var(--green)!important; }
.custom-panel { border-top:1px solid var(--line); border-bottom:1px solid var(--line);
  padding:10px 0; margin:10px 0; }
.custom-panel summary { cursor:pointer; color:var(--green); font-weight:800; margin-bottom:8px; }
.custom-panel textarea { font-family:ui-monospace,SFMono-Regular,Menlo,monospace; font-size:.78rem; }
.tiny { font-size:.75rem; color:#65716a; }
@media(max-width:900px){
  .institution-banner { align-items:flex-start; padding:8px 14px; gap:10px; }
  .banner-support { gap:8px; flex-wrap:wrap; }
  .banner-support-label { width:100%; text-align:right; }
  .banner-logo.fulbright { width:88px; }
  .banner-logo.walter { width:102px; }
  .banner-logo.berkeley { width:52px; height:52px; }
  .heatmap-reading-guide { grid-template-columns:1fr; }
  .heatmap-legend-grid { grid-template-columns:repeat(2,minmax(115px,1fr)); }
  .order-inspector-controls { grid-template-columns:1fr; }
}
@media(max-width:700px){
  .institution-banner { display:block; }
  .banner-support { justify-content:flex-start; margin-top:8px; }
  .banner-support-label { text-align:left; }
  .hero h1 { white-space:normal; }
  .metric-row{grid-template-columns:1fr}.pair-card{grid-template-columns:1fr}
}
"

app_js <- "
(function() {
  var primerTooltip = null;
  var pinnedPrimerTooltipTarget = null;

  function ensurePrimerTooltip() {
    if (primerTooltip && primerTooltip.isConnected) return primerTooltip;
    primerTooltip = document.createElement('div');
    primerTooltip.className = 'primer-hover-tooltip';
    primerTooltip.setAttribute('role', 'tooltip');
    primerTooltip.setAttribute('aria-hidden', 'true');
    document.body.appendChild(primerTooltip);
    return primerTooltip;
  }

  function placePrimerTooltip(target) {
    var tooltip = ensurePrimerTooltip();
    var targetRect = target.getBoundingClientRect();
    var tooltipRect = tooltip.getBoundingClientRect();
    var gap = 14;
    var left = targetRect.right + gap;
    if (left + tooltipRect.width > window.innerWidth - 12) {
      left = targetRect.left - tooltipRect.width - gap;
    }
    left = Math.max(12, Math.min(left, window.innerWidth - tooltipRect.width - 12));
    var top = targetRect.top + (targetRect.height - tooltipRect.height) / 2;
    top = Math.max(12, Math.min(top, window.innerHeight - tooltipRect.height - 12));
    tooltip.style.left = left + 'px';
    tooltip.style.top = top + 'px';
  }

  function showPrimerTooltip(target, pinned) {
    var message = target && target.getAttribute('data-tooltip');
    if (!message) return;
    var tooltip = ensurePrimerTooltip();
    tooltip.textContent = message;
    tooltip.classList.add('visible');
    tooltip.classList.toggle('pinned', !!pinned);
    tooltip.setAttribute('aria-hidden', 'false');
    window.requestAnimationFrame(function() {
      placePrimerTooltip(target);
    });
  }

  function hidePrimerTooltip() {
    if (!primerTooltip) return;
    primerTooltip.classList.remove('visible', 'pinned');
    primerTooltip.setAttribute('aria-hidden', 'true');
    pinnedPrimerTooltipTarget = null;
  }

  document.addEventListener('mouseover', function(event) {
    var target = event.target.closest(
      '.primer-tooltip-target[data-tooltip]'
    );
    if (!target || pinnedPrimerTooltipTarget) return;
    showPrimerTooltip(target, false);
  });

  document.addEventListener('mousemove', function(event) {
    var target = event.target.closest(
      '.primer-tooltip-target[data-tooltip]'
    );
    if (!target || pinnedPrimerTooltipTarget) return;
    placePrimerTooltip(target);
  });

  document.addEventListener('mouseout', function(event) {
    var target = event.target.closest(
      '.primer-tooltip-target[data-tooltip]'
    );
    if (!target || pinnedPrimerTooltipTarget) return;
    if (event.relatedTarget && target.contains(event.relatedTarget)) return;
    hidePrimerTooltip();
  });

  document.addEventListener('focusin', function(event) {
    var target = event.target.closest(
      '.primer-tooltip-target[data-tooltip]'
    );
    if (target && !pinnedPrimerTooltipTarget) {
      showPrimerTooltip(target, false);
    }
  });

  document.addEventListener('focusout', function(event) {
    if (!pinnedPrimerTooltipTarget &&
        event.target.closest('.primer-tooltip-target[data-tooltip]')) {
      hidePrimerTooltip();
    }
  });

  document.addEventListener('click', function(event) {
    var target = event.target.closest(
      '.primer-tooltip-target[data-tooltip]'
    );
    if (target) {
      if (pinnedPrimerTooltipTarget === target) {
        hidePrimerTooltip();
      } else {
        pinnedPrimerTooltipTarget = target;
        showPrimerTooltip(target, true);
      }
    } else if (pinnedPrimerTooltipTarget) {
      hidePrimerTooltip();
    }
  });

  document.addEventListener('keydown', function(event) {
    var target = event.target.closest(
      '.primer-tooltip-target[data-tooltip]'
    );
    if (event.key === 'Escape') {
      hidePrimerTooltip();
    } else if (target && (event.key === 'Enter' || event.key === ' ')) {
      event.preventDefault();
      if (pinnedPrimerTooltipTarget === target) {
        hidePrimerTooltip();
      } else {
        pinnedPrimerTooltipTarget = target;
        showPrimerTooltip(target, true);
      }
    }
  });

  function holdDeepDetailUntilIdle() {
    document.body.classList.add('detail-switching');
    var sawBusy = false;
    var fallback = null;

    function release() {
      document.removeEventListener('shiny:busy', onBusy);
      document.removeEventListener('shiny:idle', onIdle);
      if (fallback) window.clearTimeout(fallback);
      window.setTimeout(function() {
        document.body.classList.remove('detail-switching');
      }, 80);
    }
    function onBusy() {
      sawBusy = true;
    }
    function onIdle() {
      if (sawBusy) release();
    }

    document.addEventListener('shiny:busy', onBusy);
    document.addEventListener('shiny:idle', onIdle);
    fallback = window.setTimeout(release, 2500);
  }

  function wireHeatmapScrollers() {
    document.querySelectorAll('.heatmap-scroll-shell').forEach(function(shell) {
      var top = shell.querySelector('.heatmap-scroll-top');
      var spacer = shell.querySelector('.heatmap-scroll-top-inner');
      var body = shell.querySelector('.heatmap-scroll-body');
      var grid = shell.querySelector('.heatmap');
      if (!top || !spacer || !body || !grid) return;

      function sizeSpacer() {
        spacer.style.width = grid.scrollWidth + 'px';
        top.scrollLeft = body.scrollLeft;
      }
      sizeSpacer();

      if (shell.dataset.scrollWired === 'true') return;
      shell.dataset.scrollWired = 'true';
      var syncing = false;
      top.addEventListener('scroll', function() {
        if (syncing) return;
        syncing = true;
        body.scrollLeft = top.scrollLeft;
        syncing = false;
      });
      body.addEventListener('scroll', function() {
        if (syncing) return;
        syncing = true;
        top.scrollLeft = body.scrollLeft;
        syncing = false;
      });
      if (window.ResizeObserver) {
        var observer = new ResizeObserver(sizeSpacer);
        observer.observe(grid);
        observer.observe(body);
      }
    });
  }

  document.addEventListener('click', function(event) {
    var cell = event.target.closest('.heat-cell[data-order][data-pair]');
    if (!cell || !window.Shiny) return;
    holdDeepDetailUntilIdle();
    Shiny.setInputValue('heatmap_cell_click', {
      order: cell.dataset.order,
      pair: cell.dataset.pair,
      nonce: Date.now()
    }, {priority: 'event'});
    var detail = document.querySelector('.order-detail-layout');
    if (detail) {
      detail.scrollIntoView({behavior:'smooth', block:'start'});
      window.setTimeout(function() {
        detail.scrollIntoView({behavior:'smooth', block:'start'});
      }, 450);
    }
  });

  document.addEventListener('change', function(event) {
    if (event.target && event.target.matches(
      '#detail_order, #detail_pair, #detail_example'
    )) {
      holdDeepDetailUntilIdle();
    }
  });

  document.addEventListener('click', function(event) {
    var row = event.target.closest(
      '.sequence-source-table table.dataTable tbody tr'
    );
    if (!row || row.querySelector('.dataTables_empty')) return;
    var card = row.closest('.deep-detail-card');
    window.setTimeout(function() {
      var panel = card && card.querySelector('.sequence-evidence-panel');
      if (panel) panel.scrollIntoView({behavior:'smooth', block:'start'});
    }, 650);
  });

  document.addEventListener('shiny:connected', wireHeatmapScrollers);
  document.addEventListener('shiny:value', function() {
    window.requestAnimationFrame(wireHeatmapScrollers);
  });
  new MutationObserver(function() {
    window.requestAnimationFrame(wireHeatmapScrollers);
  }).observe(document.documentElement, {childList:true, subtree:true});

  if (window.Shiny) {
    Shiny.addCustomMessageHandler('scroll_order_detail', function(message) {
      window.setTimeout(function() {
        var detail = document.querySelector('.order-detail-layout');
        if (detail) detail.scrollIntoView({behavior:'smooth', block:'start'});
      }, 180);
    });
  } else {
    document.addEventListener('shiny:connected', function() {
      Shiny.addCustomMessageHandler('scroll_order_detail', function(message) {
        window.setTimeout(function() {
          var detail = document.querySelector('.order-detail-layout');
          if (detail) detail.scrollIntoView({behavior:'smooth', block:'start'});
        }, 180);
      });
    }, {once:true});
  }
})();
"

institution_banner <- div(
  class = "institution-banner",
  div(
    class = "banner-home",
    tags$img(
      src = "evolab-berkeley.png",
      alt = "EvoLab UC Berkeley"
    ),
    div(
      strong("EvoLab · UC Berkeley"),
      span("Open-source multi-marker primer decision support")
    )
  ),
  div(
    class = "banner-support",
    span(class = "banner-support-label", "Research home & support"),
    tags$img(
      class = "banner-logo fulbright",
      src = "fulbright.png",
      alt = "Fulbright"
    ),
    tags$img(
      class = "banner-logo walter",
      src = "walter-benjamin-dfg.png",
      alt = "DFG Walter Benjamin Programme"
    ),
    tags$img(
      class = "banner-logo berkeley",
      src = "uc-berkeley.png",
      alt = "University of California, Berkeley"
    )
  )
)

hero <- div(
  class = "hero",
  div(class = "eyebrow", "Multi-marker primer decision support"),
  h1("See the gap before it becomes bias."),
  p(
    "Compare binding geometry, thermal compatibility, and lineage-level mismatch ",
    "risk without collapsing them into a single overconfident score."
  )
)

ui <- page_navbar(
  id = "main_nav",
  title = "Primer Atlas",
  theme = atlas_theme,
  header = tagList(
    tags$style(HTML(app_css)),
    tags$script(HTML(app_js)),
    institution_banner,
    uiOutput("data_status")
  ),
  nav_panel(
    "Primer map",
    div(
      class = "container-fluid",
      hero,
      layout_sidebar(
        sidebar = sidebar(
          width = 340,
          h4("Build the view"),
          selectInput(
            "map_marker", "Marker",
            choices = setNames(
              marker_registry$marker_id[marker_registry$status %in% c("active", "pilot")],
              paste0(marker_registry$display_name[marker_registry$status %in% c("active", "pilot")],
                     ifelse(marker_registry$status[marker_registry$status %in% c("active", "pilot")] == "pilot", " · pilot", ""))
            ),
            selected = "COI"
          ),
          selectInput(
            "application_filter", "Applications",
            choices = c("barcoding", "bulk_community", "edna", "diet"),
            multiple = TRUE
          ),
          selectInput(
            "environment_filter", "Environments",
            choices = c("freshwater", "marine", "terrestrial", "host_associated"),
            multiple = TRUE
          ),
          selectInput(
            "intent_filter", "Design intents",
            choices = c("broad", "target_enriched", "exclusion_blocking"),
            multiple = TRUE
          ),
          selectInput(
            "map_sort", "Order rows by",
            choices = c(
              "Binding site along marker (5′→3′)" = "binding_site",
              "Primer-pair name (A–Z)" = "alphabetical",
              "Folmer overlap (most bp first)" = "folmer_overlap",
              "Amplicon length (shortest first)" = "amplicon_length"
            ),
            selected = "folmer_overlap"
          ),
          fluidRow(
            column(
              6,
              actionButton(
                "select_all", "Select all",
                class = "btn-primary btn-sm w-100"
              )
            ),
            column(
              6,
              actionButton(
                "clear_all", "Clear",
                class = "map-clear-btn btn-sm w-100"
              )
            )
          ),
          sliderInput(
            "amplicon_range", "Amplicon + primers (bp)",
            min = 100, max = 750, value = c(100, 750), step = 10
          ),
          checkboxGroupInput(
            "pair_select", "Primer pairs",
            choices = setNames(pair_meta$pair_id[pair_meta$marker_id == "COI"], pair_meta$pair_label[pair_meta$marker_id == "COI"]),
            selected = c(
              "MCO", "LERAY_XT", "ZBJ_ART", "FWH2",
              "BF2_BR2", "ANML", "BEEPRIME",
              "NOSPID", "NOSPI2_LAURELIN"
            )
          ),
          tags$details(
            class = "custom-panel",
            open = "open",
            tags$summary("Add a custom primer pair"),
            selectInput(
              "custom_marker", "Marker for this pair",
              choices = setNames(
                marker_registry$marker_id[marker_registry$status %in% c("active", "pilot")],
                marker_registry$display_name[marker_registry$status %in% c("active", "pilot")]
              ), selected = "COI"
            ),
            textInput("custom_pair_name", "Pair name", value = "My primer pair"),
            textAreaInput(
              "custom_forward", "Forward primer (5′→3′)",
              value = "", rows = 2,
              placeholder = "IUPAC sequence, e.g. GCHCCHGAYATRGCHTTYCC"
            ),
            textAreaInput(
              "custom_reverse", "Reverse primer (5′→3′)",
              value = "", rows = 2,
              placeholder = "IUPAC sequence, e.g. TCDGGRTGNCCRAARAAYCA"
            ),
            conditionalPanel(
              condition = "input.custom_marker == 'COI'",
              selectInput(
                "custom_target_order",
                "Expected target group (optional)",
                choices = c("Not specified" = "__NONE__", alignment_orders),
                selected = "__NONE__"
              )
            ),
            textInput(
              "custom_expected_amplicon", "Expected amplicon + primers (optional)",
              value = "", placeholder = "e.g. 461"
            ),
            actionButton(
              "locate_custom", "Add pair to map",
              class = "btn-primary btn-sm w-100"
            ),
            checkboxInput("show_custom", "Show added custom pairs", value = TRUE),
            uiOutput("custom_location_status"),
            uiOutput("custom_pair_queue"),
            actionButton(
              "clear_custom", "Clear custom pairs",
              class = "btn-outline-secondary btn-sm w-100"
            )
          ),
          uiOutput("marker_coordinate_note")
        ),
        card(
          full_screen = TRUE,
          card_header(uiOutput("map_card_title", inline = TRUE)),
          uiOutput("overview_metrics"),
          div(
            class = "map-legend",
            span(
              span(class = "legend-number", "1"),
              "Dynamic row number"
            ),
            span(span(class = "legend-line legend-forward"), "Forward primer"),
            span(span(class = "legend-line legend-informative"), "Amplicon − primers (informative)"),
            span(span(class = "legend-line legend-reverse"), "Reverse primer")
          ),
          uiOutput("primer_map"),
          uiOutput("custom_diagnostic_selector"),
          uiOutput("custom_location_detail")
        )
      )
    )
  ),
  nav_panel(
    "Order lens",
    div(
      class = "container-fluid",
      hero,
      card(
        fill = FALSE,
        class = "order-heatmap-card",
        card_header("Which orders are at risk — and why?"),
        div(
          class = "callout",
          strong("Use this heatmap for triage, not as a PCR verdict. "),
          "Mismatch fit and reference-site availability are shown separately. ",
          "Neither number is an amplification probability or a universal pass/fail score."
        ),
        br(),
        div(
          class = "heatmap-reading-guide",
          div(
            h5("How to read each cell"),
            p(
              strong("Large number — raw median PrimerMiner pair penalty. "),
              "It is calculated only among templates for which both binding sites ",
              "are scorable. Zero means no scored mismatch in the median template; ",
              "larger values indicate a poorer sequence fit in this model."
            ),
            p(
              strong("Percentage — binding-site availability. "),
              "This is the fraction of reference templates containing scorable ",
              "sequence at both primer windows. It is not taxon coverage and does ",
              "not report amplification success."
            ),
            p(
              strong("Colour — relative within the primers currently shown. "),
              "The scale is log-transformed and capped at the displayed 95th ",
              "percentile, so colours can change when the primer selection changes."
            ),
            div(
              class = "heatmap-legend-grid",
              div(
                class = "heatmap-legend-item",
                span(class = "heatmap-swatch low"),
                span("Lower displayed penalty")
              ),
              div(
                class = "heatmap-legend-item",
                span(class = "heatmap-swatch mid"),
                span("Intermediate")
              ),
              div(
                class = "heatmap-legend-item",
                span(class = "heatmap-swatch high"),
                span("Higher displayed penalty")
              ),
              div(
                class = "heatmap-legend-item",
                span(class = "heatmap-swatch missing"),
                span("No scorable pair")
              )
            )
          ),
          div(
            class = "investigation-tip",
            h5("How to investigate a warning"),
            tags$ol(
              tags$li("Hover for p90, denominator, and terminal-3′ mismatch information."),
              tags$li("Click any heatmap cell to load that order–primer comparison below."),
              tags$li(
                "Switch the template example between median, p90, worst, and most incomplete ",
                "to inspect individual mismatch components."
              ),
              tags$li(
                "When a target-specific panel appears, filter family/genus rows and click ",
                "a row to open exact organisms, accessions, and binding sequences."
              )
            )
          )
        ),
        selectizeInput(
          "order_pair_select",
          "Primer pairs shown",
          choices = setNames(pair_meta$pair_id[pair_meta$marker_id == "COI"], pair_meta$pair_label[pair_meta$marker_id == "COI"]),
          selected = pair_meta$pair_id[pair_meta$marker_id == "COI"],
          multiple = TRUE,
          width = "100%",
          options = list(
            plugins = list("remove_button"),
            placeholder = "Select one or more primer pairs"
          )
        ),
        div(
          class = "d-flex gap-2 mb-3",
          actionButton(
            "order_select_all",
            "Select all",
            class = "btn-sm btn-success"
          ),
          actionButton(
            "order_clear",
            "Clear",
            class = "btn-sm btn-outline-success map-clear-btn"
          )
        ),
        uiOutput("order_matrix")
      ),
      br(),
      layout_columns(
        col_widths = c(12, 12),
        fill = FALSE,
        fillable = FALSE,
        class = "order-detail-layout",
        card(
          fill = FALSE,
          class = "order-detail-card order-inspector-card",
          card_header("Inspect one comparison"),
          div(
            class = "order-inspector-controls",
            selectInput(
              "detail_order", "Order / study group",
              choices = if (nrow(order_scores)) {
                sort(unique(c(order_scores$order, "Acari", "Collembola")))
              } else {
                "No data"
              }
            ),
            selectInput(
              "detail_pair", "Primer pair",
              choices = setNames(pair_meta$pair_id[pair_meta$marker_id == "COI"], pair_meta$pair_label[pair_meta$marker_id == "COI"]),
              selected = "MCO"
            ),
            selectInput(
              "detail_example", "Template example",
              choices = c(
                "Median scorable template" = "median",
                "90th-percentile scorable template" = "p90",
                "Highest-penalty scorable template" = "worst",
                "Most incomplete template" = "missing"
              ),
              selected = "p90"
            )
          ),
          div(
            class = "order-inspector-actions",
            downloadButton(
              "download_detail",
              "Download raw detail",
              class = "btn-sm"
            ),
            actionButton(
              "open_lineage", "Open in Lineage Explorer",
              class = "btn-sm btn-outline-success"
            ),
            p(
              class = "tiny",
              "Penalty, terminal mismatch, adjacency, and unavailable alignment ",
              "positions remain inspectable at each primer base."
            )
          )
        ),
        card(
          fill = FALSE,
          class = "order-detail-card",
          card_header("PrimerMiner component trace"),
          div(
            class = "comparison-trace-shell",
            uiOutput("mismatch_detail")
          )
        )
      ),
      br(),
      conditionalPanel(
        condition = "input.detail_order == 'Hymenoptera' && input.detail_pair == 'BEEPRIME'",
        card(
          fill = FALSE,
          class = "deep-detail-card",
          card_header("BeePrime taxonomic validation · expanded author reference"),
          uiOutput("beeprime_deep_metrics"),
          uiOutput("beeprime_interpretation"),
          div(
            class = "callout taxon-info",
            strong("Use this panel for BeePrime interpretation. "),
            "It scores the authors’ 99.5%-clustered COI reference and retains family, ",
            "subfamily, and genus labels. The older 50-record order panel above remains ",
            "available only as a reproducibility trace."
          ),
          div(
            class = "deep-dive-controls",
            selectInput(
              "beeprime_taxon_scope",
              "Taxonomic view",
              choices = c(
                "Bee genera" = "bee_genus",
                "Bee subfamilies" = "bee_subfamily",
                "Bee families" = "bee_family",
                "All Hymenoptera families" = "hymenoptera_family",
                "Broad Hymenoptera lineages" = "hymenoptera_lineage"
              ),
              selected = "bee_genus"
            ),
            selectInput(
              "beeprime_taxon_sort",
              "Sort rows",
              choices = c(
                "Highest median penalty first" = "penalty",
                "Largest reference sample first" = "sample",
                "Taxon A–Z" = "taxon"
              ),
              selected = "penalty"
            )
          ),
          downloadButton(
            "download_beeprime_taxa",
            "Download selected taxon table",
            class = "btn-sm"
          ),
          br(), br(),
          p(
            class = "tiny mb-2",
            "Search globally, filter individual columns, or sort by any header. ",
            "Click one row to open its exact organisms, accessions, and primer-binding ",
            "sequences below."
          ),
          div(
            class = "sequence-source-table",
            DTOutput("beeprime_taxon_table", fill = FALSE)
          ),
          div(
            class = "sequence-evidence-panel beeprime-sequence-panel",
            h5("Exact sequence evidence behind the selected BeePrime row"),
            uiOutput("beeprime_sequence_intro"),
            uiOutput("beeprime_sequence_downloads"),
            DTOutput("beeprime_sequence_table", fill = FALSE)
          ),
          br(),
          layout_columns(
            col_widths = c(7, 5),
            div(
              h5("Wet-lab check from the publication"),
              uiOutput("beeprime_empirical_summary"),
              uiOutput("beeprime_empirical_failures")
            ),
            div(
              class = "callout",
              strong("Interpretation boundary"),
              p(
                "Lower PrimerMiner penalties indicate a closer sequence fit. They are ",
                "not PCR probabilities and have no universal failure threshold."
              ),
              p(
                class = "mb-0",
                "Subfamily and tribe labels come from NCBI taxonomy for the mapped bee ",
                "centroids; missing ranks remain unresolved. “Broad lineages” are ",
                "taxonomic groupings, not ecological guild assignments."
              )
            )
          )
        )
      ),
      br(),
      conditionalPanel(
        condition = "input.detail_pair && input.detail_pair != 'BEEPRIME'",
        card(
          fill = FALSE,
          class = "deep-detail-card",
          card_header("Expanded taxonomic drill-down · exact sequence evidence"),
          uiOutput("claimed_primer_metrics"),
          uiOutput("claimed_primer_interpretation"),
          uiOutput("claimed_phylogeny_group_table"),
          div(
            class = "callout taxon-info",
            strong("Every mapped primer pair opens the same evidence path. "),
            "This panel audits whether the selected primer sites are present in the ",
            "67,352-centroid alignment, then resolves the selected order through ",
            "family or genus to exact accessions and binding sequences. Incomplete ",
            "binding regions remain visible but are blocked from biological penalty ",
            "interpretation."
          ),
          div(
            class = "deep-dive-controls",
            selectInput(
              "claimed_taxon_rank",
              "Resolution within selected order",
              choices = c(
                "Families" = "family",
                "Genera" = "genus",
                "Subfamilies where available" = "subfamily"
              ),
              selected = "family"
            ),
            selectInput(
              "claimed_taxon_sort",
              "Sort rows",
              choices = c(
                "Highest median penalty first" = "penalty",
                "Lowest site coverage first" = "coverage",
                "Largest reference sample first" = "sample",
                "Taxon A–Z" = "taxon"
              ),
              selected = "penalty"
            )
          ),
          downloadButton(
            "download_claimed_taxa",
            "Download selected taxon table",
            class = "btn-sm"
          ),
          br(), br(),
          uiOutput("claimed_taxon_notice"),
          div(
            class = "sequence-source-table",
            DTOutput("claimed_taxon_table", fill = FALSE)
          ),
          div(
            class = "sequence-evidence-panel claimed-sequence-panel",
            h5("Exact sequence evidence behind the selected row"),
            uiOutput("claimed_sequence_intro"),
            uiOutput("claimed_sequence_downloads"),
            DTOutput("claimed_sequence_table", fill = FALSE)
          ),
          br(),
          layout_columns(
            col_widths = c(7, 5),
            uiOutput("claimed_primer_sources"),
            div(
              class = "callout",
              strong("Interpretation boundary"),
              p(
                "Sparse site availability and mismatch penalty are independent ",
                "problems. Both denominators remain visible at every rank."
              ),
              p(
                class = "mb-0",
                "A published target claim is treated as a hypothesis. Hybrid primer ",
                "combinations do not inherit validation from their individual components."
              )
            )
          )
        )
      )
    )
  ),
  nav_panel(
    "Lineage Explorer",
    div(
      class = "container-fluid",
      hero,
      div(
        class = "callout",
        strong("Reference evidence, not regional amplification probability. "),
        "Geography frequently proxies taxonomy and sampling bias. Every geographic view reports located/all eligible sequences, missing locations, and taxonomic composition."
      ),
      br(),
      layout_sidebar(
        sidebar = sidebar(
          width = 330,
          h4("Persistent evidence scope"),
          selectInput(
            "lineage_marker", "Marker",
            choices = setNames(
              marker_registry$marker_id[marker_registry$status %in% c("active", "pilot")],
              marker_registry$display_name[marker_registry$status %in% c("active", "pilot")]
            ), selected = "COI"
          ),
          selectizeInput("lineage_pair", "Primer pair", choices = NULL),
          selectizeInput("lineage_target", "Target group", choices = "All", selected = "All"),
          selectizeInput("lineage_country", "Reference country / territory", choices = "All", selected = "All"),
          selectizeInput("lineage_locality", "Reference locality", choices = "All", selected = "All"),
          uiOutput("lineage_denominator"),
          tags$a(
            id = "propose_primer_link", class = "btn btn-primary btn-sm w-100",
            href = "https://github.com/sven9r/primer-atlas/issues/new?template=propose-primer.yml",
            target = "_blank", "Propose this primer"
          )
        ),
        navset_card_tab(
          id = "lineage_tabs",
          nav_panel("Overview", uiOutput("lineage_overview")),
          nav_panel("Families", DTOutput("lineage_families")),
          nav_panel("Subfamilies", uiOutput("lineage_subfamilies_ui")),
          nav_panel("Genera", uiOutput("lineage_genera_ui")),
          nav_panel("Sequences", uiOutput("lineage_sequence_note"), DTOutput("lineage_sequences")),
          nav_panel("Sources and downloads", uiOutput("lineage_sources"))
        )
      )
    )
  ),
  nav_panel(
    "Thermal window",
    div(
      class = "container-fluid",
      hero,
      layout_sidebar(
        sidebar = sidebar(
          width = 290,
          numericInput("sodium", "Monovalent salt (mM)", value = 50, min = 1, max = 500),
          sliderInput("ta_offset", "Ta below limiting Tm (°C)", min = 2, max = 8, value = 4, step = 0.5),
          div(
            class = "callout tiny",
            "Degenerate bases create a Tm interval. The pair window is derived from ",
            "the lower-Tm primer and is intended as a gradient-PCR starting point."
          )
        ),
        card(
          full_screen = TRUE,
          card_header("Primer melting ranges & pair annealing windows"),
          uiOutput("thermal_metrics"),
          div(style = "overflow:auto; max-height:680px;", uiOutput("thermal_table"))
        )
      )
    )
  ),
  nav_panel(
    "Evidence & citations",
    div(
      class = "container-fluid",
      hero,
      layout_columns(
        col_widths = c(7, 5),
        card(
          card_header("Evidence registry"),
          downloadButton("download_bib", "Download BibTeX", class = "btn-primary btn-sm"),
          uiOutput("citation_list")
        ),
        card(
          card_header("What each source supports"),
          h4("Reference alignments"),
          p(
            "The current build contains ", length(alignment_files),
            " reproducible COX1 reference groups: 20 explicitly listed arthropod ",
            "orders, the Acari and Collembola composite groups, and Bivalvia for ",
            "mollusk-focused primer validation. PrimerMiner retrieval, 97% VSEARCH ",
            "clustering, and MAFFT alignment all use versioned reference NC_001322.1. ",
            "Raw, clustered, aligned, and manifest files are retained."
          ),
          p(
            "The BeePrime taxonomic drill-down instead uses Gurten et al. (2026) ",
            "Supplement 5: 67,352 author-clustered 99.5% COI centroids derived from ",
            "316,254 NCBI sequences. Author supplements supply genus/family mappings; ",
            "saved accession-level NCBI taxonomy adds subfamily and tribe where available."
          ),
          h4("Mismatch risk"),
          p("PrimerMiner combines mismatch position, mismatch type, adjacency, and ambiguity. Missing or gapped binding sites are reported as unavailable coverage rather than assigned a failure score. Its default matrices remain a hypothesis requiring empirical calibration."),
          h4("PCR evidence"),
          p("Stadhouders et al. demonstrate strong position- and type-dependent effects in controlled 5′-nuclease assays. Those results motivate the direction of risk, not a universal probability."),
          h4("Temperature"),
          p("The current estimator is deliberately transparent and fast. A production release should add a full nearest-neighbor engine with Mg²⁺, primer concentration, and chemistry-specific corrections."),
          div(class = "callout", "Open-source principle: every number should reveal its data source, assumptions, and uncertainty.")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  output$data_status <- renderUI({
    if (!is.null(coi_release_state) && isTRUE(coi_release_state$stale)) {
      div(class = "callout tiny m-2", strong("Stale-data notice: "), "R2 was unavailable, so the repository-pinned successful release is active.")
    }
  })
  observeEvent(input$select_all, {
    updateCheckboxGroupInput(
      session, "pair_select",
      selected = pair_meta$pair_id[pair_meta$marker_id == input$map_marker]
    )
  })
  observeEvent(input$clear_all, {
    updateCheckboxGroupInput(session, "pair_select", selected = character())
  })

  custom_location <- reactiveVal(NULL)
  custom_pair_collection <- reactiveVal(list())
  custom_pair_counter <- reactiveVal(0L)
  selected_custom_id <- reactiveVal(NULL)
  observeEvent(input$map_marker, {
    marker_pairs <- pair_meta |> filter(marker_id == input$map_marker)
    choices <- setNames(marker_pairs$pair_id, marker_pairs$pair_label)
    selected <- if (identical(input$map_marker, "COI")) {
      intersect(c("MCO", "LERAY_XT", "ZBJ_ART", "FWH2", "BF2_BR2", "ANML", "BEEPRIME", "NOSPID", "NOSPI2_LAURELIN"), marker_pairs$pair_id)
    } else {
      marker_pairs$pair_id
    }
    updateCheckboxGroupInput(session, "pair_select", choices = choices, selected = selected)
    if (nrow(marker_pairs)) {
      range <- range(marker_pairs$aligned_amplicon_bp, na.rm = TRUE)
      updateSliderInput(
        session, "amplicon_range",
        min = max(20, floor(range[1] / 10) * 10),
        max = ceiling(range[2] / 10) * 10,
        value = c(max(20, floor(range[1] / 10) * 10), ceiling(range[2] / 10) * 10)
      )
    }
    updateSelectInput(session, "custom_marker", selected = input$map_marker)
    if (input$map_marker != "COI") updateSelectInput(session, "map_sort", selected = "binding_site")
    session$onFlushed(function() {
      updateCheckboxGroupInput(session, "pair_select", choices = choices, selected = selected)
    }, once = TRUE)
  }, ignoreInit = FALSE)

  output$map_card_title <- renderUI({
    marker <- marker_registry |> filter(marker_id == input$map_marker) |> slice(1)
    paste0("Interactive ", marker$display_name, " binding-site overview")
  })
  output$marker_coordinate_note <- renderUI({
    marker <- marker_registry |> filter(marker_id == input$map_marker) |> slice(1)
    div(
      class = "tiny",
      "Coordinates are sequence-aligned on ", strong(marker$coordinate_reference),
      ". Low-confidence reference-specific placements are labelled and should not be treated as universal binding sites."
    )
  })
  observeEvent(input$locate_custom, {
    pair_name <- trimws(input$custom_pair_name)
    if (!nzchar(pair_name)) pair_name <- "My primer pair"
    forward_sequence <- normalize_primer_input(input$custom_forward)
    reverse_sequence <- normalize_primer_input(input$custom_reverse)
    if (!valid_primer_input(forward_sequence)) {
      custom_location(list(error = "Forward primer must contain 15–80 IUPAC DNA bases."))
      return()
    }
    if (!valid_primer_input(reverse_sequence)) {
      custom_location(list(error = "Reverse primer must contain 15–80 IUPAC DNA bases."))
      return()
    }

    expected_text <- trimws(input$custom_expected_amplicon)
    expected_amplicon <- suppressWarnings(as.numeric(expected_text))
    selected_marker <- if (is.null(input$custom_marker)) "COI" else input$custom_marker
    marker_row <- marker_registry |> filter(marker_id == selected_marker) |> slice(1)
    marker_length <- marker_row$reference_length
    if (nzchar(expected_text) && (!is.finite(expected_amplicon) ||
      expected_amplicon < 20 || expected_amplicon > marker_length)) {
      custom_location(list(
        error = paste0("Expected amplicon must be blank or a number from 20 to ", marker_length, " bp for ", marker_row$display_name, ".")
      ))
      return()
    }
    if (!nzchar(expected_text)) expected_amplicon <- NA_real_
    target_order <- trimws(input$custom_target_order)
    if (!nzchar(target_order) || identical(target_order, "__NONE__")) {
      target_order <- NA_character_
    }

    selected_reference <- if (selected_marker == "COI") coi_reference else its_reference
    selected_templates <- if (selected_marker == "COI") order_alignment_templates else NULL
    located <- tryCatch(
      locate_primer_pair(
        reference_sequence = selected_reference,
        forward_sequence = forward_sequence,
        reverse_sequence = reverse_sequence,
        expected_amplicon_bp = expected_amplicon,
        alignment_templates = selected_templates
      ),
      error = function(e) e
    )
    if (inherits(located, "error")) {
      custom_location(list(error = conditionMessage(located)))
      return()
    }

    custom_pair_counter(custom_pair_counter() + 1L)
    custom_id <- paste0("CUSTOM_", custom_pair_counter())
    best <- located$best
    targeted_start <- best$end_forward + 1L
    targeted_end <- best$start_reverse - 1L
    folmer_overlap <- if (selected_marker == "COI") max(
      0L, min(targeted_end, folmer_region$end) - max(targeted_start, folmer_region$start) + 1L
    ) else 0L
    geometry <- tibble(
      pair_id = custom_id,
      pair_label = pair_name,
      use_case = "custom input",
      target = if (is.na(target_order)) {
        "User-supplied pair"
      } else {
        paste0("Expected target: ", target_order)
      },
      forward_primer = "Custom forward",
      reverse_primer = "Custom reverse",
      forward_start = best$start_forward,
      forward_end = best$end_forward,
      reverse_start = best$start_reverse,
      reverse_end = best$end_reverse,
      pair_start = best$start_forward,
      pair_end = best$end_reverse,
      aligned_amplicon_bp = best$aligned_amplicon_bp,
      targeted_region_start = targeted_start,
      targeted_region_end = targeted_end,
      targeted_region_bp = best$targeted_region_bp,
      folmer_overlap_bp = folmer_overlap,
      folmer_coverage_fraction = if (selected_marker == "COI") folmer_overlap / folmer_region$length_bp else NA_real_,
      reference_accession = marker_row$reference_version,
      sources = "custom user input",
      folmer_overlap_rank = NA_integer_,
      marker_id = selected_marker,
      placement_status = "sequence_aligned"
    )
    custom_primers <- tibble(
      pair_id = custom_id,
      pair_label = pair_name,
      primer_name = c(paste0(pair_name, " F"), paste0(pair_name, " R")),
      direction = c("forward", "reverse"),
      sequence = c(forward_sequence, reverse_sequence),
      reported_amplicon_bp = NA_real_,
      use_case = "custom input",
      target = if (is.na(target_order)) {
        "User-supplied pair"
      } else {
        paste0("Expected target: ", target_order)
      },
      reported_ta_c = NA_real_,
      source_key = "custom",
      source_note = "Entered interactively; not saved to the primer library"
    )
    order_preference <- if (selected_marker == "COI") summarize_order_preference(
      located$pair_order_support,
      expected_order = target_order
    ) else list(headline = "Fungal ITS pilot placement", detail = "This session pair was screened only against FN812768.2.", expected = "Expanded fungal lineage evidence is not yet available.")
    order_detail <- if (selected_marker == "COI") tryCatch(
      withProgress(
        expr = score_primer_pair_primerminer(
          alignment_files = alignment_files,
          primer_rows = custom_primers,
          geometry = geometry
        ),
        message = paste0("Scoring ", pair_name, " across installed orders"),
        detail = "Calculating PrimerMiner penalties and position-level traces.",
        value = 0.6
      ),
      error = function(e) e
    ) else list(order_scores = tibble(), pair_template_scores = tibble(), primer_position_scores = tibble())
    if (inherits(order_detail, "error")) {
      custom_location(list(
        error = paste0(
          "A credible COI location was found, but the order-level PrimerMiner ",
          "analysis failed: ", conditionMessage(order_detail)
        )
      ))
      return()
    }

    expanded_detail <- if (selected_marker == "COI" && file.exists(expanded_reference_alignment_path)) tryCatch(
      withProgress(
        expr = score_primer_pair_expanded_reference(
          alignment_path = expanded_reference_alignment_path,
          taxonomy = expanded_reference_taxonomy,
          primer_rows = custom_primers,
          geometry = geometry,
          alignment_offset = expanded_reference_alignment_offset
        ),
        message = paste0(
          "Building family, genus, and exact-sequence evidence for ",
          pair_name
        ),
        detail = "Scoring both primers on 67,352 clustered COI centroids.",
        value = 0.85
      ),
      error = function(e) e
    ) else list(overall = list(
      reference_suitable = FALSE,
      reference_limitation_reason = if (selected_marker == "COI") {
        "Expanded exact-sequence scoring is not bundled on this worker; curated evidence remains available as lazy release artifacts."
      } else {
        "Expanded fungal ITS lineage evidence is not available in the pilot release."
      }
    ), scores = tibble())
    if (inherits(expanded_detail, "error")) {
      custom_location(list(
        error = paste0(
          "The pair passed the COI-location and order checks, but the expanded ",
          "sequence drill-down failed: ", conditionMessage(expanded_detail)
        )
      ))
      return()
    }
    expanded_scores_path <- if (selected_marker == "COI" && nrow(expanded_detail$scores)) tempfile(
      pattern = paste0(tolower(custom_id), "_expanded_scores_"),
      fileext = ".rds"
    ) else ""
    if (selected_marker == "COI" && nrow(expanded_detail$scores)) saveRDS(expanded_detail$scores, expanded_scores_path, compress = FALSE)

    result <- c(
      located,
      list(
        error = NULL,
        pair_name = pair_name,
        marker_id = selected_marker,
        expected_amplicon = expected_amplicon,
        target_order = target_order,
        order_preference = order_preference,
        geometry = geometry,
        primers = custom_primers,
        order_scores = order_detail$order_scores,
        pair_template_scores = order_detail$pair_template_scores,
        primer_position_scores = order_detail$primer_position_scores,
        expanded_scores_path = expanded_scores_path,
        expanded_overall = expanded_detail$overall
      )
    )
    custom_location(result)
    collection <- custom_pair_collection()
    collection[[custom_id]] <- result
    custom_pair_collection(collection)
    selected_custom_id(custom_id)
  }, ignoreInit = TRUE)

  observeEvent(input$clear_custom, {
    collection <- custom_pair_collection()
    paths <- vapply(
      collection,
      function(x) {
        if (is.null(x$expanded_scores_path)) "" else x$expanded_scores_path
      },
      character(1)
    )
    paths <- paths[nzchar(paths) & file.exists(paths)]
    if (length(paths)) unlink(paths)
    custom_pair_collection(list())
    custom_location(NULL)
    selected_custom_id(NULL)
  })

  session$onSessionEnded(function() {
    collection <- isolate(custom_pair_collection())
    paths <- vapply(
      collection,
      function(x) {
        if (is.null(x$expanded_scores_path)) "" else x$expanded_scores_path
      },
      character(1)
    )
    paths <- paths[nzchar(paths) & file.exists(paths)]
    if (length(paths)) unlink(paths)
  })

  filtered_pairs <- reactive({
    req(input$amplicon_range)
    x <- pair_meta |>
      filter(
        marker_id == input$map_marker,
        aligned_amplicon_bp >= input$amplicon_range[1],
        aligned_amplicon_bp <= input$amplicon_range[2]
      )
    facet_match <- function(rows, type, selected) {
      if (is.null(selected) || !length(selected)) return(rows)
      allowed <- primer_facets |>
        filter(facet_type == type, facet_value %in% selected) |>
        pull(pair_id) |>
        unique()
      rows |> filter(pair_id %in% allowed)
    }
    x <- facet_match(x, "application", input$application_filter)
    x <- facet_match(x, "environment", input$environment_filter)
    x <- facet_match(x, "design_intent", input$intent_filter)
    if (is.null(input$pair_select) || !length(input$pair_select)) {
      x <- x[0, ]
    } else {
      x <- x |> filter(pair_id %in% input$pair_select)
    }
    custom <- custom_pair_collection()
    if (length(custom) && isTRUE(input$show_custom)) {
      custom_rows <- bind_rows(lapply(custom, `[[`, "geometry")) |>
        filter(marker_id == input$map_marker)
      x <- bind_rows(x, custom_rows)
    }
    x
  })

  sorted_filtered_pairs <- reactive({
    x <- filtered_pairs()
    if (!nrow(x)) return(x)
    sort_mode <- if (is.null(input$map_sort) || !nzchar(input$map_sort)) {
      "folmer_overlap"
    } else {
      input$map_sort
    }
    sort_primer_map_rows(x, sort_mode)
  })

  sorted_pair_choices <- reactive({
    sort_mode <- if (is.null(input$map_sort) || !nzchar(input$map_sort)) {
      "folmer_overlap"
    } else {
      input$map_sort
    }
    sort_primer_map_rows(pair_meta |> filter(marker_id == input$map_marker), sort_mode)
  })

  pair_choice_label_signature <- reactiveVal(NULL)
  observe({
    choices <- number_primer_choice_labels(
      sorted_pair_choices(),
      sorted_filtered_pairs()
    )
    signature <- paste(names(choices), choices, sep = "=", collapse = "|")
    if (identical(signature, pair_choice_label_signature())) return()
    current <- isolate(input$pair_select)
    if (is.null(current)) return()
    pair_choice_label_signature(signature)
    freezeReactiveValue(input, "pair_select")
    updateCheckboxGroupInput(
      session,
      "pair_select",
      choices = choices,
      selected = current
    )
  })

  displayed_primers <- reactive({
    custom <- custom_pair_collection()
    if (length(custom) && isTRUE(input$show_custom)) {
      bind_rows(primers, bind_rows(lapply(custom, `[[`, "primers")))
    } else {
      primers
    }
  })

  output$custom_location_status <- renderUI({
    custom <- custom_location()
    if (is.null(custom)) {
      return(div(class = "tiny", "Paste both oligos and click Add pair to map."))
    }
    if (!is.null(custom$error)) {
      return(div(
        class = "callout tiny mt-2",
        strong("Not added: "),
        custom$error
      ))
    }
    best <- custom$best
    supported_orders <- if (nrow(custom$pair_order_support) &&
      all(c("n_pair_supported", "pair_support_fraction") %in% names(custom$pair_order_support))) {
      custom$pair_order_support |>
        filter(n_pair_supported >= 2, pair_support_fraction >= 0.10) |>
        slice_head(n = 3)
    } else tibble()
    div(
      class = "tiny mt-2",
      strong("Accepted and added: "),
      paste0(
        best$start_forward, "–", best$end_forward, " / ",
        best$start_reverse, "–", best$end_reverse,
        " · ", best$aligned_amplicon_bp, " bp + primers / ",
        best$targeted_region_bp, " bp − primers"
      ),
      br(),
      paste(length(custom_pair_collection()), "custom pair(s) currently on the map."),
      br(),
      strong("Order compatibility: "),
      custom$order_preference$headline,
      if (nzchar(custom$order_preference$expected)) {
        tagList(br(), custom$order_preference$expected)
      },
      if (!isTRUE(custom$expanded_overall$reference_suitable)) {
        tagList(
          br(),
          strong("Expanded sequence drill-down: "),
          custom$expanded_overall$reference_limitation_reason,
          " The pair remains accepted as a marker-specific screening placement, not as proof of amplification."
        )
      },
      if (custom$alignment_rescue_used && nrow(supported_orders)) {
        tagList(
          br(),
          strong("Taxon-aware evidence: "),
          paste0(
            "the single Drosophila reference failed, but the pair is supported ",
            "by installed COI order alignments (",
            paste0(
              supported_orders$order, " ",
              supported_orders$n_pair_supported, "/",
              supported_orders$n_pair_available,
              collapse = "; "
            ),
            "). This validates the COI location, not primer specificity."
          )
        )
      }
    )
  })

  output$custom_pair_queue <- renderUI({
    collection <- custom_pair_collection()
    if (!length(collection)) return(NULL)
    shown <- sorted_filtered_pairs()
    displayed_numbers <- setNames(seq_len(nrow(shown)), shown$pair_id)
    div(
      class = "tiny mb-2",
      tagList(lapply(names(collection), function(pair_id) {
        number <- unname(displayed_numbers[pair_id])
        label <- if (length(number) && !is.na(number)) {
          paste0(number, " · ", collection[[pair_id]]$pair_name)
        } else {
          collection[[pair_id]]$pair_name
        }
        span(class = "pill", label)
      }))
    )
  })

  observeEvent(input$custom_diagnostic_select, {
    if (
      !is.null(input$custom_diagnostic_select) &&
        input$custom_diagnostic_select %in% names(custom_pair_collection())
    ) {
      selected_custom_id(input$custom_diagnostic_select)
    }
  }, ignoreInit = TRUE)

  output$custom_diagnostic_selector <- renderUI({
    collection <- custom_pair_collection()
    if (!length(collection)) return(NULL)
    current <- selected_custom_id()
    if (is.null(current) || !current %in% names(collection)) {
      current <- tail(names(collection), 1)
    }
    selectInput(
      "custom_diagnostic_select",
      "Custom-pair statistics",
      choices = setNames(
        names(collection),
        vapply(collection, `[[`, character(1), "pair_name")
      ),
      selected = current,
      width = "420px"
    )
  })

  selected_custom_result <- reactive({
    collection <- custom_pair_collection()
    if (!length(collection)) return(NULL)
    selected <- selected_custom_id()
    if (is.null(selected) || !selected %in% names(collection)) {
      selected <- tail(names(collection), 1)
    }
    collection[[selected]]
  })

  order_pair_choices <- reactive({
    curated <- setNames(pair_meta$pair_id[pair_meta$marker_id == "COI"], pair_meta$pair_label[pair_meta$marker_id == "COI"])
    collection <- custom_pair_collection()
    collection <- collection[vapply(collection, function(x) identical(x$marker_id, "COI"), logical(1))]
    if (!length(collection)) return(curated)
    custom <- setNames(
      names(collection),
      paste0(
        vapply(collection, `[[`, character(1), "pair_name"),
        " · custom"
      )
    )
    c(curated, custom)
  })

  order_lens_scores <- reactive({
    collection <- custom_pair_collection()
    if (!length(collection)) return(order_scores)
    bind_rows(
      order_scores,
      bind_rows(lapply(collection, `[[`, "order_scores"))
    )
  })

  order_lens_template_scores <- reactive({
    collection <- custom_pair_collection()
    if (!length(collection)) return(pair_template_scores)
    bind_rows(
      pair_template_scores,
      bind_rows(lapply(collection, `[[`, "pair_template_scores"))
    )
  })

  order_lens_position_scores <- reactive({
    collection <- custom_pair_collection()
    if (!length(collection)) return(primer_position_scores)
    bind_rows(
      primer_position_scores,
      bind_rows(lapply(collection, `[[`, "primer_position_scores"))
    ) |>
      mutate(unavailable_reason = coalesce(unavailable_reason, ""))
  })

  order_lens_primers <- reactive({
    collection <- custom_pair_collection()
    if (!length(collection)) return(primers)
    bind_rows(
      primers,
      bind_rows(lapply(collection, `[[`, "primers"))
    )
  })

  observeEvent(custom_pair_collection(), {
    choices <- order_pair_choices()
    values <- unname(choices)
    collection <- custom_pair_collection()
    current_matrix <- isolate(input$order_pair_select)
    if (is.null(current_matrix) && !length(collection)) {
      current_matrix <- pair_meta$pair_id
    }
    selected_matrix <- intersect(current_matrix, values)
    if (length(collection)) {
      selected_matrix <- unique(c(selected_matrix, names(collection)))
    }
    updateSelectizeInput(
      session,
      "order_pair_select",
      choices = choices,
      selected = selected_matrix,
      server = TRUE
    )

    current_detail <- isolate(input$detail_pair)
    selected_detail <- if (length(collection)) {
      tail(names(collection), 1)
    } else if (!is.null(current_detail) && current_detail %in% values) {
      current_detail
    } else {
      "MCO"
    }
    updateSelectInput(
      session,
      "detail_pair",
      choices = choices,
      selected = selected_detail
    )
  }, ignoreInit = FALSE)

  observeEvent(input$order_select_all, {
    choices <- order_pair_choices()
    updateSelectizeInput(
      session,
      "order_pair_select",
      choices = choices,
      selected = unname(choices),
      server = TRUE
    )
  })

  observeEvent(input$order_clear, {
    updateSelectizeInput(
      session,
      "order_pair_select",
      choices = order_pair_choices(),
      selected = character(),
      server = TRUE
    )
  })

  observeEvent(input$heatmap_cell_click, {
    click <- input$heatmap_cell_click
    req(click$order, click$pair)
    valid_pair <- click$pair %in% unname(order_pair_choices())
    valid_order <- click$order %in% unique(order_lens_scores()$order)
    if (!valid_pair || !valid_order) return()
    freezeReactiveValue(input, "detail_order")
    freezeReactiveValue(input, "detail_pair")
    updateSelectInput(
      session,
      "detail_order",
      selected = click$order
    )
    updateSelectInput(
      session,
      "detail_pair",
      choices = order_pair_choices(),
      selected = click$pair
    )
    session$sendCustomMessage("scroll_order_detail", list())
  })

  output$overview_metrics <- renderUI({
    x <- filtered_pairs()
    source_count <- primers |>
      filter(pair_id %in% x$pair_id) |>
      pull(source_key) |>
      unique() |>
      length()
    div(
      class = "metric-row",
      div(class = "metric", tags$b(nrow(x)), span("primer pairs shown")),
      div(
        class = "metric",
        tags$b(if (nrow(x)) paste0(min(x$aligned_amplicon_bp), "–", max(x$aligned_amplicon_bp)) else "—"),
        span("amplicon + primers (bp)")
      ),
      div(
        class = "metric",
        tags$b(if (nrow(x)) paste0(min(x$targeted_region_bp), "–", max(x$targeted_region_bp)) else "—"),
        span("amplicon − primers (bp)")
      ),
      div(class = "metric", tags$b(source_count), span("linked source records"))
      ,div(
        class = "metric",
        tags$b(sum(x$placement_status != "sequence_aligned", na.rm = TRUE)),
        span("reference-specific placement warnings")
      )
    )
  })

  output$primer_map <- renderUI({
    x <- sorted_filtered_pairs()
    if (!nrow(x)) return(div(class = "callout", "No primer pairs match the current filters."))
    map_primers <- displayed_primers()
    marker <- marker_registry |> filter(marker_id == input$map_marker) |> slice(1)
    landmarks <- marker_landmarks |> filter(marker_id == input$map_marker) |> arrange(display_order)
    reference_length <- marker$reference_length
    focus <- if (input$map_marker == "COI") {
      list(start = folmer_region$start, end = folmer_region$end)
    } else {
      list(start = min(landmarks$start), end = max(landmarks$end))
    }
    width <- 1680
    left <- 112
    right <- 1000
    axis_domain <- primer_map_domain(
      x,
      folmer_start = focus$start,
      folmer_end = focus$end,
      tick_bp = 100,
      reference_length = reference_length
    )
    axis_min <- unname(axis_domain["min"])
    axis_max <- unname(axis_domain["max"])
    scale_x <- function(bp) {
      bounded <- pmin(pmax(bp, axis_min), axis_max)
      left + (bounded - axis_min) / (axis_max - axis_min) * (right - left)
    }
    row_h <- 58
    first_row_y <- 88
    last_row_y <- first_row_y + (nrow(x) - 1) * row_h
    folmer_y <- last_row_y + 54
    coi_y <- folmer_y + 40
    height <- coi_y + 42
    start_col <- 1090
    end_col <- 1160
    inclusive_col <- 1250
    informative_col <- 1350
    citation_col <- 1440
    number_x <- 16
    forward_color <- "#2563eb"
    informative_color <- "#9aa6b2"
    reverse_color <- "#dc2626"
    shared_source_color <- "#7c3aed"
    source_link <- function(source_row, x, y, color, prefix = "") {
      if (!nrow(source_row) || is.na(source_row$url[1]) ||
        !nzchar(source_row$url[1])) {
        return(tags$text(
          x = x, y = y, fill = "#6b746f", `font-size` = 10.5,
          paste0(prefix, "source unavailable")
        ))
      }
      tags$a(
        href = source_row$url[1],
        target = "_blank",
        rel = "noopener noreferrer",
        style = "cursor:pointer",
        tags$text(
          x = x, y = y, fill = color, `font-size` = 10.5,
          `text-decoration` = "underline",
          paste0(prefix, source_row$short_citation[1]),
          tags$title(source_row$title[1])
        )
      )
    }
    children <- list(
      tags$text(
        x = start_col, y = 20, `text-anchor` = "middle",
        fill = "#48554e", `font-size` = 11, `font-weight` = 700, "Start bp"
      ),
      tags$text(
        x = end_col, y = 20, `text-anchor` = "middle",
        fill = "#48554e", `font-size` = 11, `font-weight` = 700, "End bp"
      ),
      tags$text(
        x = inclusive_col, y = 20, `text-anchor` = "middle",
        fill = "#48554e", `font-size` = 11, `font-weight` = 700, "+ primers"
      ),
      tags$text(
        x = informative_col, y = 20, `text-anchor` = "middle",
        fill = "#48554e", `font-size` = 11, `font-weight` = 700, "− primers"
      ),
      tags$text(
        x = citation_col, y = 20,
        fill = "#48554e", `font-size` = 11, `font-weight` = 700, "Citation(s)"
      ),
      tags$text(
        x = right, y = 20, `text-anchor` = "end",
        fill = "#48554e", `font-size` = 11, `font-weight` = 700,
        paste0("View ", axis_min, "–", axis_max, " bp")
      ),
      tags$line(
        x1 = left, y1 = 52, x2 = right, y2 = 52,
        stroke = "#6d756f", `stroke-width` = 1.5
      )
    )
    ticks <- seq(axis_min, axis_max, 100)
    for (tick in ticks) {
      tx <- scale_x(tick)
      children <- append(children, list(
        tags$line(
          x1 = tx, y1 = 46, x2 = tx, y2 = last_row_y + 25,
          stroke = "#dedbd2", `stroke-width` = 1
        ),
        tags$text(
          x = tx, y = 40, `text-anchor` = "middle",
          fill = "#68736c", `font-size` = 11, tick
        )
      ))
    }
    children <- append(children, list(
      tags$rect(
        x = scale_x(focus$start), y = 55,
        width = scale_x(focus$end) - scale_x(focus$start),
        height = last_row_y - 30, fill = "#f5df73", opacity = .15
      )
    ))
    for (i in seq_len(nrow(x))) {
      row <- x[i, ]
      y <- first_row_y + (i - 1) * row_h
      f <- map_primers |> filter(pair_id == row$pair_id, direction == "forward") |> slice(1)
      r <- map_primers |> filter(pair_id == row$pair_id, direction == "reverse") |> slice(1)
      f1 <- scale_x(row$forward_start)
      f2 <- scale_x(row$forward_end)
      r1 <- scale_x(row$reverse_start)
      r2 <- scale_x(row$reverse_end)
      forward_source <- citations |> filter(key == f$source_key[1]) |> slice(1)
      reverse_source <- citations |> filter(key == r$source_key[1]) |> slice(1)
      use_tags <- primer_use_tags(row$use_case, row$target)
      tooltip <- paste0(
        "Map row ", i, ": ", row$pair_label,
        "\nTags: ", paste(use_tags, collapse = " · "),
        "\nPrimers: ", f$primer_name, " + ", r$primer_name,
        "\nGene coordinates: ", row$pair_start, "–", row$pair_end,
        "\nAmplicon + primers: ", row$aligned_amplicon_bp, " bp",
        "\nAmplicon − primers (informative): ", row$targeted_region_bp, " bp",
        if (input$map_marker == "COI") paste0("\nInformative bp inside Folmer: ", row$folmer_overlap_bp) else "",
        "\nCoordinates derived from ", marker$coordinate_reference,
        if (!is.na(row$placement_status) && row$placement_status != "sequence_aligned") paste0("\nPlacement warning: ", row$placement_status) else ""
      )
      number_tooltip <- paste0(
        "Map row ", i, "\nPair: ", row$pair_label,
        "\nOligos: ", f$primer_name, " + ", r$primer_name,
        "\nSorting and filtering recalculate this number."
      )
      row_children <- list(
        tags$circle(
          cx = number_x, cy = y, r = 13,
          fill = "#17241d",
          tags$title(number_tooltip)
        ),
        tags$text(
          x = number_x, y = y + 4,
          `text-anchor` = "middle",
          fill = "#fffdf8", `font-size` = 11, `font-weight` = 800,
          i,
          tags$title(number_tooltip)
        ),
        tags$line(
          x1 = left, y1 = y + row_h / 2 - 2,
          x2 = width - 12, y2 = y + row_h / 2 - 2,
          stroke = "#e5e2d9", `stroke-width` = 1
        ),
        tags$line(
          x1 = f2, y1 = y, x2 = r1, y2 = y,
          class = "primer-tooltip-target",
          `data-tooltip` = tooltip,
          `aria-label` = paste0("Show map-row details for ", row$pair_label),
          tabindex = "0",
          stroke = informative_color, `stroke-width` = 5, `stroke-linecap` = "round",
          `pointer-events` = "stroke"
        ),
        tags$line(
          x1 = f1, y1 = y, x2 = f2, y2 = y,
          class = "primer-tooltip-target",
          `data-tooltip` = tooltip,
          `aria-label` = paste0("Show map-row details for ", row$pair_label),
          tabindex = "0",
          stroke = forward_color, `stroke-width` = 7, `stroke-linecap` = "round"
        ),
        tags$polygon(
          points = paste(f2, y, f2 - 8, y - 6, f2 - 8, y + 6),
          fill = forward_color,
          class = "primer-tooltip-target",
          `data-tooltip` = tooltip,
          `aria-label` = paste0("Show map-row details for ", row$pair_label),
          tabindex = "0"
        ),
        tags$text(
          x = f1 - 7, y = y + 4, `text-anchor` = "end",
          class = "primer-tooltip-target",
          `data-tooltip` = tooltip,
          `aria-label` = paste0("Show map-row details for ", row$pair_label),
          tabindex = "0",
          fill = forward_color, `font-size` = 12, `font-weight` = 700,
          f$primer_name
        ),
        tags$line(
          x1 = r1, y1 = y, x2 = r2, y2 = y,
          class = "primer-tooltip-target",
          `data-tooltip` = tooltip,
          `aria-label` = paste0("Show map-row details for ", row$pair_label),
          tabindex = "0",
          stroke = reverse_color, `stroke-width` = 7, `stroke-linecap` = "round"
        ),
        tags$polygon(
          points = paste(r1, y, r1 + 8, y - 6, r1 + 8, y + 6),
          fill = reverse_color,
          class = "primer-tooltip-target",
          `data-tooltip` = tooltip,
          `aria-label` = paste0("Show map-row details for ", row$pair_label),
          tabindex = "0"
        ),
        tags$text(
          x = r2 + 7, y = y + 4, `text-anchor` = "start",
          class = "primer-tooltip-target",
          `data-tooltip` = tooltip,
          `aria-label` = paste0("Show map-row details for ", row$pair_label),
          tabindex = "0",
          fill = reverse_color, `font-size` = 12, `font-weight` = 700,
          r$primer_name
        ),
        tags$text(
          x = start_col, y = y + 4, `text-anchor` = "middle",
          fill = "#34423a", `font-size` = 11, row$pair_start
        ),
        tags$text(
          x = end_col, y = y + 4, `text-anchor` = "middle",
          fill = "#34423a", `font-size` = 11, row$pair_end
        ),
        tags$text(
          x = inclusive_col, y = y + 4, `text-anchor` = "middle",
          fill = "#34423a", `font-size` = 11,
          paste0(row$aligned_amplicon_bp, " bp")
        ),
        tags$text(
          x = informative_col, y = y + 4, `text-anchor` = "middle",
          fill = "#34423a", `font-size` = 11,
          paste0(row$targeted_region_bp, " bp")
        )
      )
      if (identical(f$source_key[1], "custom") ||
        identical(r$source_key[1], "custom")) {
        row_children <- append(row_children, list(tags$text(
          x = citation_col, y = y + 4, fill = "#6b746f", `font-size` = 10.5,
          "Session input · not saved"
        )))
      } else if (identical(f$source_key[1], r$source_key[1])) {
        row_children <- append(row_children, list(
          source_link(forward_source, citation_col, y + 4, shared_source_color)
        ))
      } else {
        row_children <- append(row_children, list(
          source_link(forward_source, citation_col, y - 5, forward_color, "F: "),
          source_link(reverse_source, citation_col, y + 11, reverse_color, "R: ")
        ))
      }
      children <- append(children, row_children)
    }
    if (input$map_marker == "COI") {
      children <- append(children, list(
        tags$text(x = left - 14, y = folmer_y + 4, `text-anchor` = "end", fill = "#6b5720", `font-size` = 12, `font-weight` = 700, "Folmer region"),
        tags$line(x1 = scale_x(folmer_region$start), y1 = folmer_y, x2 = scale_x(folmer_region$end), y2 = folmer_y, stroke = "#d7a62b", `stroke-width` = 10, `stroke-linecap` = "round"),
        tags$text(x = scale_x(folmer_region$end) + 10, y = folmer_y + 4, fill = "#6b5720", `font-size` = 11, paste0(folmer_region$start, "–", folmer_region$end, " · ", folmer_region$length_bp, " bp"))
      ))
    } else {
      for (j in seq_len(nrow(landmarks))) {
        lm <- landmarks[j, ]
        children <- append(children, list(
          tags$rect(x = scale_x(lm$start), y = folmer_y - 9, width = max(2, scale_x(lm$end) - scale_x(lm$start)), height = 18, fill = if (j %% 2) "#9aa6a0" else "#d7dad7"),
          tags$text(x = (scale_x(lm$start) + scale_x(lm$end)) / 2, y = folmer_y + 4, `text-anchor` = "middle", fill = "#17241d", `font-size` = 10, `font-weight` = 700, lm$label)
        ))
      }
    }
    children <- append(children, list(
      tags$text(x = left - 14, y = coi_y + 4, `text-anchor` = "end", fill = "#17241d", `font-size` = 12, `font-weight` = 700, marker$display_name),
      tags$line(x1 = left, y1 = coi_y, x2 = right, y2 = coi_y, stroke = "#17241d", `stroke-width` = 8, `stroke-linecap` = "round"),
      tags$text(x = right + 10, y = coi_y + 4, fill = "#48554e", `font-size` = 11, paste0(axis_min, "–", min(axis_max, reference_length), " shown · ", marker$reference_version, " is ", format(reference_length, big.mark = ","), " bp"))
    ))
    div(
      class = "primer-map-wrap",
      tags$svg(
        class = "primer-map", viewBox = paste(0, 0, width, height),
        role = "img",
        `aria-label` = paste0(
          "Interactive ", marker$display_name, " primer binding map showing base pairs ",
          axis_min, " through ", axis_max,
          ". Rows are numbered in their current displayed order."
        ),
        children
      )
    )
  })

  output$custom_location_detail <- renderUI({
    custom <- selected_custom_result()
    if (is.null(custom)) return(NULL)
    best <- custom$best
    geometry <- custom$geometry
    f_tm <- tm_range(custom$primers$sequence[custom$primers$direction == "forward"])
    r_tm <- tm_range(custom$primers$sequence[custom$primers$direction == "reverse"])
    pair_candidates <- custom$pair_candidates |>
      transmute(
        Rank = pair_candidate_rank,
        `Forward site` = paste0(start_forward, "–", end_forward),
        `Reverse site` = paste0(start_reverse, "–", end_reverse),
        `Amplicon bp` = aligned_amplicon_bp,
        `F + R mismatches` = paste0(
          mismatch_count_forward, " + ", mismatch_count_reverse
        ),
        `F + R identity` = paste0(
          round(compatible_identity_forward * 100, 1), "% + ",
          round(compatible_identity_reverse * 100, 1), "%"
        ),
        `Terminal-3 F + R` = paste0(
          terminal_3_mismatches_forward, " + ", terminal_3_mismatches_reverse
        ),
        `Terminal-5 mismatches` = total_terminal_5_mismatches,
        `Length difference` = ifelse(
          is.na(expected_length_difference), "—", expected_length_difference
        )
      )
    alternative_rows <- bind_rows(
      custom$forward_candidates |>
        mutate(Direction = "forward", Rank = row_number()),
      custom$reverse_candidates |>
        mutate(Direction = "reverse", Rank = row_number())
    ) |>
      transmute(
        Direction,
        Rank,
        Site = paste0(start, "–", end),
        Mismatches = mismatch_count,
        `Compatible identity` = paste0(round(compatible_identity * 100, 1), "%"),
        `Terminal-3` = terminal_3_mismatches,
        `Terminal-5` = terminal_5_mismatches,
        `Longest compatible run` = longest_compatible_run,
        `Chance-match expectation` = format(
          random_match_evalue, digits = 2, scientific = TRUE
        ),
        `Relative localization component` = localization_penalty
      )

    html_table <- function(x) {
      tags$table(
        class = "table-atlas",
        tags$thead(tags$tr(lapply(names(x), tags$th))),
        tags$tbody(lapply(seq_len(nrow(x)), function(i) {
          tags$tr(lapply(x[i, , drop = FALSE], tags$td))
        }))
      )
    }
    ambiguity <- custom$equally_best_pair_count > 1L ||
      custom$forward_equally_best_count > 1L ||
      custom$reverse_equally_best_count > 1L
    expected_note <- if (custom$expected_amplicon_supplied) {
      paste0(
        "Expected length ", custom$expected_amplicon, " bp; selected candidate differs by ",
        best$expected_length_difference, " bp."
      )
    } else {
      paste(
        "No expected length supplied; candidates were ranked by total mismatch,",
        "terminal-five mismatch, relative localization component, then amplicon length."
      )
    }
    order_support_rows <- custom$pair_order_support |>
      filter(n_pair_available > 0) |>
      transmute(
        Order = order,
        `Templates inspected` = n_templates,
        `Both sites available` = n_pair_available,
        `Pair-supported templates` = n_pair_supported,
        `Support among available` = paste0(
          round(pair_support_fraction * 100, 1), "%"
        )
      )
    preference_callout <- div(
      class = paste(
        "callout tiny mt-2",
        if (grepl("^Warning:", custom$order_preference$expected)) {
          "taxon-alert"
        } else {
          "taxon-info"
        }
      ),
      strong("Order-level compatibility warning: "),
      custom$order_preference$headline,
      br(),
      custom$order_preference$detail,
      if (nzchar(custom$order_preference$expected)) {
        tagList(br(), strong(custom$order_preference$expected))
      },
      br(),
      em(
        "This summarizes sequence compatibility in the installed order ",
        "alignments; it does not demonstrate PCR specificity or amplification probability."
      )
    )
    evidence_callout <- if (custom$alignment_rescue_used) {
      div(
        class = "callout tiny mt-2",
        strong("Taxon-aware COI gate passed: "),
        "at least one primer did not pass on the Drosophila coordinate ",
        "reference alone, but the same coordinate pair is supported by at ",
        "least two templates and at least 10% of available templates in an ",
        "installed order alignment. This validates a plausible COI location; ",
        "it does not by itself demonstrate taxonomic specificity."
      )
    } else {
      div(
        class = "callout tiny mt-2",
        strong("Coordinate-reference COI gate passed: "),
        "each selected site has at least 80% compatible identity, no more than ",
        "one mismatch in the terminal three primer bases, an uninterrupted ",
        "compatible run of at least 7 bp, and a chance-match expectation no greater than 0.05 ",
        "across NC_001322.1."
      )
    }

    card(
      class = "mt-3",
      card_header(paste0("Custom placement diagnostics · ", custom$pair_name)),
      div(
        class = "metric-row",
        div(
          class = "metric",
          tags$b(geometry$aligned_amplicon_bp),
          span("amplicon + primers (bp)")
        ),
        div(
          class = "metric",
          tags$b(geometry$targeted_region_bp),
          span("amplicon − primers (informative bp)")
        ),
        div(
          class = "metric",
          tags$b(geometry$folmer_overlap_bp),
          span("informative bp inside Folmer")
        ),
        div(
          class = "metric",
          tags$b(paste0(best$mismatch_count_forward, " + ", best$mismatch_count_reverse)),
          span("reference mismatches F + R")
        )
      ),
      div(
        class = if (ambiguity) "callout" else "tiny",
        strong(if (ambiguity) "Ambiguity warning: " else "Placement rule: "),
        expected_note,
        if (ambiguity) {
          paste0(
            " Equally best sites: forward ", custom$forward_equally_best_count,
            ", reverse ", custom$reverse_equally_best_count,
            ", paired combinations ", custom$equally_best_pair_count, "."
          )
        }
      ),
      preference_callout,
      evidence_callout,
      h5(class = "mt-3", "Selected oligos and thermal estimates"),
      tags$table(
        class = "table-atlas",
        tags$thead(tags$tr(lapply(
          c(
            "Primer", "Sequence 5′→3′", "Selected site", "Reference mismatches",
            "Compatible identity", "Terminal-3", "Longest run",
            "Chance-match expectation", "Placement evidence", "Tm at 50 mM Na⁺"
          ),
          tags$th
        ))),
        tags$tbody(
          tags$tr(
            tags$td("Forward"),
            tags$td(class = "mono", custom$primers$sequence[1]),
            tags$td(paste0(best$start_forward, "–", best$end_forward)),
            tags$td(best$mismatch_count_forward),
            tags$td(paste0(round(best$compatible_identity_forward * 100, 1), "%")),
            tags$td(best$terminal_3_mismatches_forward),
            tags$td(paste0(best$longest_compatible_run_forward, " bp")),
            tags$td(format(
              best$random_match_evalue_forward, digits = 2, scientific = TRUE
            )),
            tags$td(best$evidence_basis_forward),
            tags$td(sprintf("%.1f–%.1f °C", f_tm[1], f_tm[2]))
          ),
          tags$tr(
            tags$td("Reverse"),
            tags$td(class = "mono", custom$primers$sequence[2]),
            tags$td(paste0(best$start_reverse, "–", best$end_reverse)),
            tags$td(best$mismatch_count_reverse),
            tags$td(paste0(round(best$compatible_identity_reverse * 100, 1), "%")),
            tags$td(best$terminal_3_mismatches_reverse),
            tags$td(paste0(best$longest_compatible_run_reverse, " bp")),
            tags$td(format(
              best$random_match_evalue_reverse, digits = 2, scientific = TRUE
            )),
            tags$td(best$evidence_basis_reverse),
            tags$td(sprintf("%.1f–%.1f °C", r_tm[1], r_tm[2]))
          )
        )
      ),
      if (nrow(order_support_rows)) {
        tagList(
          h5(class = "mt-3", "Pair support in installed COI order alignments"),
          p(
            class = "tiny",
            "A template is counted only when both binding windows are present ",
            "and both primers pass the sequence-compatibility checks. These ",
            "broad order alignments are localization evidence, not a specificity assay."
          ),
          div(style = "overflow:auto;", html_table(order_support_rows))
        )
      },
      h5(class = "mt-3", "Leading paired placements"),
      div(style = "overflow:auto;", html_table(pair_candidates)),
      tags$details(
        class = "mt-3",
        tags$summary(strong("Show leading sites for each individual primer")),
        div(style = "overflow:auto;", html_table(alternative_rows))
      ),
      p(
        class = "tiny mt-3",
        "This locates the pair on NC_001322.1 and keeps it only for this session. ",
        "Its complete PrimerMiner order matrix and position-level traces are ",
        "available in the Order lens. ",
        "High mismatch, alignment-rescued, or ambiguous placements should be ",
        "checked before interpretation."
      )
    )
  })

  output$order_matrix <- renderUI({
    scores <- order_lens_scores()
    if (!nrow(scores)) {
      return(div(class = "callout", "Reference alignments are not installed. Run the alignment and scoring scripts."))
    }
    pairs <- intersect(input$order_pair_select, unname(order_pair_choices()))
    custom_ids <- names(custom_pair_collection())
    pairs <- c(
      intersect(custom_ids, pairs),
      setdiff(pairs, custom_ids)
    )
    if (!length(pairs)) {
      return(div(
        class = "callout",
        "Select at least one primer pair to display the order comparison."
      ))
    }
    scores <- scores |> filter(pair_id %in% pairs)
    orders <- sort(unique(scores$order))
    grid <- c("150px", rep("minmax(82px,1fr)", length(pairs)))
    cells <- list(div())
    finite_penalties <- scores$pair_median_penalty[
      is.finite(scores$pair_median_penalty)
    ]
    penalty_cap <- if (length(finite_penalties)) {
      max(1, unname(quantile(finite_penalties, .95, na.rm = TRUE)))
    } else {
      1
    }
    penalty_palette <- colorRampPalette(c("#b9ddc9", "#efd787", "#efaa90"))(101)
    penalty_color <- function(value) {
      if (is.na(value)) return(NULL)
      scaled <- pmin(1, log1p(value) / log1p(penalty_cap))
      penalty_palette[max(1, round(scaled * 100) + 1)]
    }
    choices <- order_pair_choices()
    for (pair in pairs) {
      label <- names(choices)[match(pair, choices)]
      cells <- append(cells, list(div(class = "heat-head", label)))
    }
    for (ord in orders) {
      cells <- append(cells, list(div(class = "order-label", ord)))
      for (pair in pairs) {
        row <- scores |> filter(order == ord, pair_id == pair)
        if (!nrow(row)) {
          cells <- append(cells, list(div(
            class = "heat-cell heat-cell-static risk-na",
            "—"
          )))
        } else if (is.na(row$pair_median_penalty)) {
          cells <- append(cells, list(
            tags$button(
              type = "button",
              class = "heat-cell risk-na",
              `data-order` = ord,
              `data-pair` = pair,
              `aria-label` = paste0(
                ord, " · ", row$pair_label,
                ": no pair-scorable templates. Open detailed comparison."
              ),
              title = paste0(
                "No pair-scorable templates; ", row$n_templates,
                " consensus templates inspected. Click for details."
              ),
              div("NA"), tags$small("0% available")
            )
          ))
        } else {
          cells <- append(cells, list(
            tags$button(
              type = "button",
              class = "heat-cell",
              `data-order` = ord,
              `data-pair` = pair,
              `aria-label` = paste0(
                ord, " · ", row$pair_label,
                ": raw median penalty ", round(row$pair_median_penalty, 1),
                ", both sites available in ",
                round(row$pair_scorable_fraction * 100),
                "% of templates. Open detailed comparison."
              ),
              style = paste0("background:", penalty_color(row$pair_median_penalty), ";"),
              title = paste0(
                "Raw pair median ", round(row$pair_median_penalty, 1),
                "; p90 ", round(row$pair_p90_penalty, 1),
                "; both binding sites scorable in ", row$n_pair_scorable,
                "/", row$n_templates, " templates",
                "; terminal-3 mismatch fraction ",
                round(row$terminal_3_mismatch_fraction * 100),
                "%. Click for details."
              ),
              div(round(row$pair_median_penalty, 1)),
              tags$small(paste0(
                round(row$pair_scorable_fraction * 100),
                "% available"
              ))
            )
          ))
        }
      }
    }
    tagList(
      div(
        class = "heatmap-scroll-shell",
        div(
          class = "heatmap-scroll-top",
          tabindex = "0",
          role = "region",
          `aria-label` = "Horizontal scrollbar for the order-primer heatmap",
          div(class = "heatmap-scroll-top-inner")
        ),
        div(
          class = "heatmap-scroll-body",
          div(
            class = "heatmap",
            style = paste0(
              "grid-template-columns:",
              paste(grid, collapse = " "),
              ";"
            ),
            cells
          )
        )
      ),
      p(
        class = "tiny mt-2",
        "The synchronized scrollbar above remains available while reading lower rows. ",
        "Click a cell to open its detailed PrimerMiner trace."
      )
    )
  })

  output$mismatch_detail <- renderUI({
    req(
      input$detail_order,
      input$detail_pair,
      input$detail_example,
      cancelOutput = TRUE
    )
    candidates <- order_lens_template_scores() |>
      filter(order == input$detail_order, pair_id == input$detail_pair)
    pair <- order_lens_primers() |> filter(pair_id == input$detail_pair)
    summary_row <- order_lens_scores() |>
      filter(order == input$detail_order, pair_id == input$detail_pair)
    if (!nrow(candidates) || nrow(pair) < 2 || !nrow(summary_row)) {
      return(div(class = "callout", "No alignment detail is available for this comparison."))
    }

    if (input$detail_example == "missing") {
      example <- candidates |>
        mutate(
          missing_total = coalesce(forward_missing_binding_positions, 0) +
            coalesce(reverse_missing_binding_positions, 0)
        ) |>
        slice_max(missing_total, n = 1, with_ties = FALSE)
    } else {
      scorable <- candidates |> filter(pair_scorable, is.finite(pair_penalty))
      if (!nrow(scorable)) {
        example <- candidates |>
          mutate(
            missing_total = coalesce(forward_missing_binding_positions, 0) +
              coalesce(reverse_missing_binding_positions, 0)
          ) |>
          slice_min(missing_total, n = 1, with_ties = FALSE)
      }
      if (nrow(scorable)) {
        target <- switch(
          input$detail_example,
          median = median(scorable$pair_penalty),
          p90 = unname(quantile(scorable$pair_penalty, .9)),
          worst = max(scorable$pair_penalty)
        )
        example <- scorable |>
          mutate(distance_to_target = abs(pair_penalty - target)) |>
          slice_min(distance_to_target, n = 1, with_ties = FALSE)
      }
    }

    selected_positions <- order_lens_position_scores() |>
      filter(
        order == input$detail_order,
        pair_id == input$detail_pair,
        template == example$template
      )

    strip <- function(primer_row, direction) {
      detail <- selected_positions |>
        filter(primer_name == primer_row$primer_name, direction == !!direction) |>
        arrange(primer_position_5_to_3)
      if (!nrow(detail)) return(div(class = "callout", "No position records."))
      div(
        h5(primer_row$primer_name, span(class = "pill", direction)),
        div(
          class = "sequence-strip",
          lapply(seq_len(nrow(detail)), function(i) {
            row <- detail[i, ]
            div(
              class = paste(
                "base",
                if (nzchar(row$unavailable_reason)) "unavailable" else if (row$mismatch) "mismatch"
              ),
              title = paste0(
                "Primer ", row$primer_base, " / template ", row$template_base,
                "; reference position ", row$alignment_position,
                "; distance from 3′ ", row$distance_from_3prime,
                "; component penalty ", ifelse(is.na(row$position_penalty), "NA", row$position_penalty)
              ),
              row$primer_base,
              tags$small(paste0(
                row$template_base, " · ",
                ifelse(is.na(row$position_penalty), "NA", round(row$position_penalty, 1))
              ))
            )
          })
        )
      )
    }

    component_rows <- selected_positions |>
      filter(mismatch | nzchar(unavailable_reason)) |>
      arrange(direction, distance_from_3prime)
    component_table <- if (!nrow(component_rows)) {
      div(class = "callout", "This template has no mismatches or unavailable bases in either binding site.")
    } else {
      tags$table(
        class = "table-atlas",
        tags$thead(tags$tr(lapply(
          c("Primer", "Dir.", "Ref. pos.", "3′ distance", "Primer/template", "Component", "Availability"),
          tags$th
        ))),
        tags$tbody(lapply(seq_len(nrow(component_rows)), function(i) {
          row <- component_rows[i, ]
          tags$tr(
            tags$td(row$primer_name),
            tags$td(row$direction),
            tags$td(row$alignment_position),
            tags$td(row$distance_from_3prime),
            tags$td(class = "mono", paste0(row$primer_base, "/", row$template_base)),
            tags$td(ifelse(is.na(row$position_penalty), "NA", round(row$position_penalty, 2))),
            tags$td(ifelse(nzchar(row$unavailable_reason), row$unavailable_reason, "scorable mismatch"))
          )
        }))
      )
    }

    div(
      if (grepl("^CUSTOM_", input$detail_pair)) {
        div(
          class = "callout taxon-info mb-3",
          strong("Session custom pair: "),
          "these PrimerMiner penalties and position-level components were ",
          "calculated on demand against the installed order alignments. The ",
          "expanded panel below independently scores the same pair on 67,352 ",
          "taxonomically annotated COI centroids."
        )
      } else if (identical(input$detail_order, "Hymenoptera") &&
          identical(input$detail_pair, "BEEPRIME")) {
        div(
          class = "callout taxon-alert mb-3",
          strong("Legacy-panel warning: "),
          "these metrics come from 50 randomly retrieved Hymenoptera records ",
          "clustered into anonymous 97% consensuses. The shown 90th-percentile ",
          "score is one example—not a bee-wide estimate. Use the expanded ",
          "taxonomic panel below."
        )
      } else {
        div(
          class = "callout taxon-alert mb-3",
          strong("Legacy-panel warning: "),
          "the component trace uses the original 50-record order build. Use the ",
          "shared expanded-reference panel below for family/genus denominators and ",
          "exact accession-level evidence."
        )
      },
      div(
        class = "metric-row",
        div(class = "metric", tags$b(summary_row$n_templates), span(
          if (identical(input$detail_order, "Hymenoptera") &&
              identical(input$detail_pair, "BEEPRIME")) {
            "legacy consensus templates"
          } else {
            "consensus templates"
          }
        )),
        div(class = "metric", tags$b(paste0(round(summary_row$pair_scorable_fraction * 100), "%")), span("both sites available")),
        div(class = "metric", tags$b(ifelse(is.na(summary_row$forward_median_penalty), "NA", round(summary_row$forward_median_penalty, 1))), span("forward median")),
        div(class = "metric", tags$b(ifelse(is.na(summary_row$reverse_median_penalty), "NA", round(summary_row$reverse_median_penalty, 1))), span("reverse median")),
        div(class = "metric", tags$b(ifelse(is.na(example$pair_penalty), "NA", round(example$pair_penalty, 1))), span("shown pair penalty"))
      ),
      p(
        class = "tiny",
        "Template: ", code(example$template),
        " · small labels show template base and final position component after adjacency adjustment."
      ),
      strip(pair |> filter(direction == "forward") |> slice(1), "forward"),
      strip(pair |> filter(direction == "reverse") |> slice(1), "reverse"),
      div(class = "detail-table-scroll", component_table),
      div(
        class = "callout mt-3",
        "Raw PrimerMiner components are relative penalties. Coral marks a scored mismatch; ",
        "gray marks a gap/N that is excluded from the score and counted in coverage."
      )
    )
  })
  outputOptions(output, "mismatch_detail", suspendWhenHidden = FALSE)

  beeprime_taxon_data <- reactive({
    genus_to_taxon <- beeprime_bee_genus |>
      distinct(family, subfamily, genus)
    empirical_mapped <- beeprime_empirical_genus |>
      left_join(genus_to_taxon, by = "genus")

    empirical_by_family <- empirical_mapped |>
      filter(!is.na(family)) |>
      group_by(family) |>
      summarise(
        n_empirical_tested = sum(n_empirical_tested),
        n_empirical_detected = sum(n_empirical_detected),
        .groups = "drop"
      )
    empirical_by_subfamily <- empirical_mapped |>
      filter(!is.na(family), !is.na(subfamily)) |>
      group_by(family, subfamily) |>
      summarise(
        n_empirical_tested = sum(n_empirical_tested),
        n_empirical_detected = sum(n_empirical_detected),
        .groups = "drop"
      )

    selected <- switch(
      input$beeprime_taxon_scope,
      bee_family = beeprime_bee_family |>
        left_join(empirical_by_family, by = "family") |>
        mutate(taxon = family, context = "Bee family"),
      bee_subfamily = beeprime_bee_subfamily |>
        left_join(empirical_by_subfamily, by = c("family", "subfamily")) |>
        mutate(taxon = paste(family, subfamily, sep = " · "), context = "Bee subfamily"),
      bee_genus = beeprime_bee_genus |>
        left_join(
          beeprime_empirical_genus |>
            select(genus, n_empirical_tested, n_empirical_detected),
          by = "genus"
        ) |>
        mutate(taxon = genus, context = paste(family, subfamily, sep = " · ")),
      hymenoptera_family = beeprime_hymenoptera_family |>
        mutate(
          taxon = family,
          context = lineage_group,
          n_empirical_tested = NA_real_,
          n_empirical_detected = NA_real_
        ),
      hymenoptera_lineage = beeprime_hymenoptera_lineage |>
        mutate(
          taxon = lineage_group,
          context = "Broad taxonomic lineage",
          n_empirical_tested = NA_real_,
          n_empirical_detected = NA_real_
        )
    )

    selected <- selected |>
      mutate(
        n_empirical_tested = coalesce(n_empirical_tested, 0),
        n_empirical_detected = coalesce(n_empirical_detected, 0)
      )

    switch(
      input$beeprime_taxon_sort,
      sample = selected |> arrange(desc(n_centroids), taxon),
      taxon = selected |> arrange(taxon),
      selected |> arrange(desc(pair_median_penalty), taxon)
    )
  })

  output$beeprime_deep_metrics <- renderUI({
    bee <- beeprime_hymenoptera_lineage |> filter(lineage_group == "Bee")
    div(
      class = "metric-row",
      div(class = "metric", tags$b("67,352"), span("author COI centroids")),
      div(class = "metric", tags$b(format(bee$n_centroids, big.mark = ",")), span("bee centroids")),
      div(
        class = "metric",
        tags$b(paste0(round(bee$pair_scorable_fraction * 100, 1), "%")),
        span("bee pairs scorable")
      ),
      div(
        class = "metric",
        tags$b(round(bee$pair_median_penalty, 1)),
        span("bee median penalty")
      ),
      div(class = "metric", tags$b("33 / 38"), span("article PCR detections"))
    )
  })

  output$beeprime_interpretation <- renderUI({
    lineage <- beeprime_hymenoptera_lineage |>
      select(lineage_group, pair_median_penalty) |>
      arrange(pair_median_penalty)
    lineage_text <- paste0(
      lineage$lineage_group,
      " ", round(lineage$pair_median_penalty, 1),
      collapse = "; "
    )
    div(
      class = "callout mb-3",
      strong("Expanded result: BeePrime looks bee-preferring, not broadly Hymenoptera-targeting. "),
      "Bees have the lowest lineage-level median penalty (", lineage_text, "). ",
      "This is an in-silico selectivity pattern, not a categorical amplification guarantee."
    )
  })

  output$beeprime_taxon_table <- DT::renderDT({
    x <- beeprime_taxon_data()
    validate(need(nrow(x), "No rows are available for this taxonomic view."))
    evidence_text <- vapply(seq_len(nrow(x)), function(i) {
      row <- x[i, ]
      if (row$n_centroids == 0) {
        paste0(
          "0 mapped centroids · ",
          if ("species_listed_in_supplement" %in% names(row)) {
            paste0(row$species_listed_in_supplement, " species listed")
          } else {
            "catalog only"
          }
        )
      } else {
        paste0(
          format(row$n_centroids, big.mark = ","), " centroids · ",
          format(row$source_records_represented, big.mark = ","), " records",
          if (row$n_centroids < 5) " · small n" else ""
        )
      }
    }, character(1))
    both_sites <- ifelse(
      x$n_centroids == 0,
      "not mapped",
      paste0(
        x$n_pair_scorable, "/", x$n_centroids, " (",
        round(x$pair_scorable_fraction * 100, 1), "%)"
      )
    )
    wet_lab <- ifelse(
      x$n_empirical_tested > 0,
      paste0(x$n_empirical_detected, "/", x$n_empirical_tested),
      "not tested"
    )
    display <- data.frame(
      Taxon = x$taxon,
      Context = x$context,
      Evidence = evidence_text,
      `Both sites` = both_sites,
      `Forward median` = round(x$forward_median_penalty, 1),
      `Reverse median` = round(x$reverse_median_penalty, 1),
      `Pair median` = round(x$pair_median_penalty, 1),
      `Pair p90` = round(x$pair_p90_penalty, 1),
      `Terminal-3′ mismatch` =
        round(x$terminal_3_mismatch_fraction_among_scorable * 100, 1),
      `Wet-lab detection` = wet_lab,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    table <- DT::datatable(
      display,
      rownames = FALSE,
      filter = "top",
      selection = list(mode = "single", target = "row"),
      class = "compact stripe hover order-column",
      options = list(
        pageLength = 20,
        lengthMenu = c(10, 20, 50, 100),
        scrollX = TRUE,
        scrollY = "520px",
        scrollCollapse = TRUE,
        autoWidth = TRUE,
        order = list(),
        search = list(regex = FALSE, caseInsensitive = TRUE),
        language = list(
          search = "Search all columns:",
          select = list(rows = "%d row selected")
        ),
        columnDefs = list(list(
          targets = 4:8,
          render = JS(
            "function(data, type, row, meta) {",
            "  if (data === null || data === '' || data === 'NA') {",
            "    return type === 'display' ? '—' : null;",
            "  }",
            "  return data;",
            "}"
          )
        ))
      )
    )
    DT::formatRound(
      table,
      columns = c(
        "Forward median", "Reverse median", "Pair median", "Pair p90",
        "Terminal-3′ mismatch"
      ),
      digits = 1
    )
  }, server = TRUE)

  beeprime_selected_taxon <- reactive({
    selected <- input$beeprime_taxon_table_rows_selected
    x <- beeprime_taxon_data()
    if (length(selected) != 1L || selected < 1L || selected > nrow(x)) {
      return(NULL)
    }
    x[selected, , drop = FALSE]
  })

  beeprime_selected_sequence_data <- reactive({
    selected <- beeprime_selected_taxon()
    if (is.null(selected)) return(tibble())
    x <- beeprime_hymenoptera_scores
    x <- switch(
      input$beeprime_taxon_scope,
      bee_family = x |>
        filter(is_bee %in% TRUE, family == selected$family),
      bee_subfamily = x |>
        filter(
          is_bee %in% TRUE,
          family == selected$family,
          subfamily == selected$subfamily
        ),
      bee_genus = x |>
        filter(
          is_bee %in% TRUE,
          family == selected$family,
          genus == selected$genus
        ),
      hymenoptera_family = x |>
        filter(family == selected$family),
      hymenoptera_lineage = x |>
        filter(lineage_group == selected$lineage_group),
      tibble()
    )
    x |>
      arrange(desc(cluster_size), scientific_name_ncbi, organism_label, accession)
  })

  beeprime_sequence_display <- reactive({
    x <- beeprime_selected_sequence_data()
    if (!nrow(x)) return(data.frame())
    organism <- coalesce(
      x$scientific_name_ncbi,
      x$organism_label,
      x$accession
    )
    accession_link <- paste0(
      "<a href=\"https://www.ncbi.nlm.nih.gov/nuccore/",
      vapply(x$accession, URLencode, character(1), reserved = TRUE),
      "\" target=\"_blank\" rel=\"noopener\">",
      x$accession,
      "</a>"
    )
    data.frame(
      Accession = accession_link,
      Organism = organism,
      Lineage = x$lineage_group,
      Family = x$family,
      Subfamily = x$subfamily,
      Genus = x$genus,
      `Records represented` = x$cluster_size,
      `Both sites` = ifelse(x$pair_scorable %in% TRUE, "yes", "no"),
      `Pair penalty` = round(x$pair_penalty, 2),
      `F mismatches` = x$forward_mismatch_count,
      `R mismatches` = x$reverse_mismatch_count,
      `Terminal-3′ mismatches` =
        coalesce(x$forward_terminal_3_mismatches, 0) +
        coalesce(x$reverse_terminal_3_mismatches, 0),
      `Forward binding sequence (primer-oriented 5′→3′)` =
        x$forward_binding_sequence_primer_oriented,
      `Reverse binding sequence (primer-oriented 5′→3′)` =
        x$reverse_binding_sequence_primer_oriented,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })

  output$beeprime_sequence_intro <- renderUI({
    selected <- beeprime_selected_taxon()
    if (is.null(selected)) {
      return(div(
        class = "callout",
        "Select a family, subfamily, genus, or broad-lineage row above to inspect ",
        "the exact scored reference sequences."
      ))
    }
    x <- beeprime_selected_sequence_data()
    if (!nrow(x)) {
      return(div(
        class = "callout taxon-alert",
        strong(paste0(selected$taxon, ". ")),
        "No scored centroid sequence is mapped to this row. This can happen for ",
        "a publication-listed genus that is absent from the clustered reference."
      ))
    }
    div(
      class = "callout taxon-info",
      strong(paste0(selected$taxon, ". ")),
      paste0(
        format(nrow(x), big.mark = ","), " exact clustered centroid",
        if (nrow(x) == 1L) "" else "s",
        " were scored, representing ",
        format(sum(x$cluster_size, na.rm = TRUE), big.mark = ","),
        " source records. Click an accession to open its NCBI record; the binding ",
        "columns show the exact primer-oriented sites used for the penalty."
      )
    )
  })

  output$beeprime_sequence_downloads <- renderUI({
    if (is.null(beeprime_selected_taxon()) ||
        !nrow(beeprime_selected_sequence_data())) {
      return(NULL)
    }
    div(
      class = "sequence-evidence-actions",
      downloadButton(
        "download_beeprime_sequence_evidence",
        "Download evidence CSV",
        class = "btn-sm"
      ),
      downloadButton(
        "download_beeprime_sequences_fasta",
        "Download exact FASTA",
        class = "btn-sm"
      )
    )
  })

  output$beeprime_sequence_table <- DT::renderDT({
    x <- beeprime_sequence_display()
    if (!nrow(x)) {
      return(DT::datatable(
        data.frame(`Exact sequence evidence` = character()),
        rownames = FALSE,
        selection = "none",
        options = list(
          dom = "t",
          language = list(
            emptyTable = "Select a BeePrime taxon row to load exact sequences."
          )
        )
      ))
    }
    table <- DT::datatable(
      x,
      rownames = FALSE,
      filter = "top",
      escape = FALSE,
      selection = "none",
      class = "compact stripe hover order-column",
      options = list(
        pageLength = 15,
        lengthMenu = c(10, 15, 25, 50, 100),
        scrollX = TRUE,
        scrollY = "380px",
        scrollCollapse = TRUE,
        autoWidth = TRUE,
        order = list(list(6, "desc"), list(1, "asc")),
        search = list(regex = FALSE, caseInsensitive = TRUE),
        language = list(search = "Search exact evidence:")
      )
    )
    DT::formatStyle(
      table,
      columns = c(
        "Forward binding sequence (primer-oriented 5′→3′)",
        "Reverse binding sequence (primer-oriented 5′→3′)"
      ),
      fontFamily = "ui-monospace, SFMono-Regular, Menlo, monospace",
      whiteSpace = "nowrap"
    )
  })

  output$beeprime_empirical_summary <- renderUI({
    tagList(
      p(
        strong("Article text: "),
        "33 of 38 bee extracts amplified (86.8%). ",
        tags$a(
          href = "https://doi.org/10.3897/mbmg.10.183708",
          target = "_blank",
          "Open the paper"
        )
      ),
      p(
        strong("Supplement 1 table: "),
        "32 of 37 listed target rows amplified. The one-row difference from the ",
        "article is retained as a source discrepancy; both sources identify five failures."
      )
    )
  })

  output$beeprime_empirical_failures <- renderUI({
    failures <- beeprime_empirical |>
      filter(target == "yes", detected == "no")
    tags$table(
      class = "table-atlas",
      tags$thead(tags$tr(tags$th("Five listed failures"), tags$th("Genus"))),
      tags$tbody(lapply(seq_len(nrow(failures)), function(i) {
        tags$tr(tags$td(failures$taxon[i]), tags$td(failures$genus[i]))
      }))
    )
  })

  output$download_beeprime_taxa <- downloadHandler(
    filename = function() {
      paste0("BeePrime_", input$beeprime_taxon_scope, "_PrimerMiner_summary.csv")
    },
    content = function(file) {
      write.csv(beeprime_taxon_data(), file, row.names = FALSE, na = "")
    }
  )

  beeprime_selected_sequence_filename <- reactive({
    selected <- beeprime_selected_taxon()
    req(selected)
    taxon_slug <- gsub(
      "(^_+|_+$)",
      "",
      gsub("[^A-Za-z0-9]+", "_", selected$taxon)
    )
    paste("BEEPRIME", input$beeprime_taxon_scope, taxon_slug, sep = "_")
  })

  output$download_beeprime_sequence_evidence <- downloadHandler(
    filename = function() {
      paste0(beeprime_selected_sequence_filename(), "_sequence_evidence.csv")
    },
    content = function(file) {
      x <- beeprime_selected_sequence_data()
      req(nrow(x))
      write.csv(x, file, row.names = FALSE, na = "")
    }
  )

  output$download_beeprime_sequences_fasta <- downloadHandler(
    filename = function() {
      paste0(beeprime_selected_sequence_filename(), "_exact_sequences.fasta")
    },
    content = function(file) {
      selected_rows <- beeprime_selected_sequence_data()
      req(nrow(selected_rows))
      source_path <- file.path(
        "data", "external", "gurten2026", "ClusteredReferences.fasta"
      )
      validate(need(file.exists(source_path), "The source FASTA is unavailable."))
      sequences <- parse_fasta(source_path)
      fasta_accessions <- sub(
        "^_R_",
        "",
        sub(" .*", "", names(sequences))
      )
      keep <- fasta_accessions %in% selected_rows$accession
      validate(need(
        any(keep),
        "None of the selected accessions could be recovered from the source FASTA."
      ))
      sequences <- sequences[keep]
      lines <- unlist(lapply(seq_along(sequences), function(i) {
        sequence <- sequences[[i]]
        starts <- seq(1L, nchar(sequence), by = 80L)
        c(
          paste0(">", names(sequences)[i]),
          substring(
            sequence,
            starts,
            pmin(starts + 79L, nchar(sequence))
          )
        )
      }), use.names = FALSE)
      writeLines(lines, file, useBytes = TRUE)
    }
  )

  custom_detail_active <- reactive({
    req(input$detail_pair)
    grepl("^CUSTOM_", input$detail_pair)
  })

  custom_detail_record <- reactive({
    req(custom_detail_active())
    collection <- custom_pair_collection()
    record <- collection[[input$detail_pair]]
    validate(need(
      !is.null(record),
      "This custom primer pair is no longer available in the current session."
    ))
    record
  })

  claimed_pair_overall <- reactive({
    req(input$detail_pair)
    if (isTRUE(custom_detail_active())) {
      return(custom_detail_record()$expanded_overall)
    }
    claimed_primer_overall |>
      filter(pair_id == input$detail_pair) |>
      slice(1)
  })

  targeted_zbj_active <- reactive({
    input$detail_pair %in% c("ZBJ_ART", "ZBJ_ART_DEG") &&
      input$detail_order %in% targeted_zbj_overall$source_group
  })

  targeted_zbj_scope <- reactive({
    req(targeted_zbj_active())
    targeted_zbj_overall |>
      filter(
        pair_id == input$detail_pair,
        source_group == input$detail_order
      ) |>
      slice(1)
  })

  claimed_reference_usable <- reactive({
    if (isTRUE(targeted_zbj_active())) {
      scope <- targeted_zbj_scope()
      return(
        nrow(scope) == 1L &&
          scope$n_pair_scorable >= 50L &&
          scope$pair_scorable_fraction >= 0.80
      )
    }
    isTRUE(claimed_pair_overall()$reference_suitable)
  })

  claimed_selected_scope <- reactive({
    req(input$detail_pair, input$detail_order)
    if (isTRUE(custom_detail_active())) {
      x <- claimed_raw_scores() |>
        filter(order %in% claimed_selected_orders())
      if (!nrow(x)) return(data.frame())
      return(
        summarize_expanded_primer_scores(
          x,
          c("pair_id", "pair_label")
        )
      )
    }
    if (isTRUE(targeted_zbj_active())) {
      return(
        targeted_zbj_scope() |>
          mutate(
            n_centroids = n_sequences,
            source_records_represented = n_sequences
          )
      )
    }
    if (input$detail_order %in% c("Acari", "Collembola")) {
      return(
        claimed_primer_phylogeny_groups |>
          filter(
            pair_id == input$detail_pair,
            display_group == input$detail_order
          ) |>
          slice(1)
      )
    }
    claimed_primer_orders |>
      filter(
        pair_id == input$detail_pair,
        order == input$detail_order
      ) |>
      slice(1)
  })

  claimed_taxon_data <- reactive({
    req(input$detail_pair, input$detail_order, input$claimed_taxon_rank)
    if (isTRUE(custom_detail_active())) {
      x <- claimed_raw_scores() |>
        filter(order %in% claimed_selected_orders())
      if (!nrow(x)) return(data.frame())
      selected <- switch(
        input$claimed_taxon_rank,
        family = summarize_expanded_primer_scores(
          x,
          c("pair_id", "pair_label", "order", "family")
        ) |>
          mutate(taxon = family, context = order),
        genus = summarize_expanded_primer_scores(
          x |>
            filter(!is.na(family), nzchar(family)),
          c("pair_id", "pair_label", "order", "family", "genus")
        ) |>
          mutate(taxon = genus, context = paste(order, family, sep = " · ")),
        subfamily = summarize_expanded_primer_scores(
          x |>
            filter(
              !is.na(family), nzchar(family),
              !is.na(subfamily), nzchar(subfamily)
            ),
          c(
            "pair_id", "pair_label", "order", "family", "subfamily"
          )
        ) |>
          mutate(
            taxon = subfamily,
            context = paste(order, family, sep = " · ")
          )
      )
      return(switch(
        input$claimed_taxon_sort,
        coverage = selected |>
          arrange(pair_scorable_fraction, desc(n_centroids), taxon),
        sample = selected |> arrange(desc(n_centroids), taxon),
        taxon = selected |> arrange(taxon),
        selected |> arrange(desc(pair_median_penalty), taxon)
      ))
    }
    if (isTRUE(targeted_zbj_active())) {
      selected <- switch(
        input$claimed_taxon_rank,
        family = targeted_zbj_families |>
          filter(
            pair_id == input$detail_pair,
            source_group == input$detail_order
          ) |>
          mutate(taxon = family, context = source_group),
        genus = targeted_zbj_genera |>
          filter(
            pair_id == input$detail_pair,
            source_group == input$detail_order
          ) |>
          mutate(taxon = genus, context = paste(source_group, family, sep = " · ")),
        subfamily = targeted_zbj_subfamilies |>
          filter(
            pair_id == input$detail_pair,
            source_group == input$detail_order
          ) |>
          mutate(
            taxon = subfamily,
            context = paste(source_group, family, sep = " · ")
          )
      )
      selected <- selected |>
        mutate(
          n_centroids = n_sequences,
          source_records_represented = n_sequences
        )
      return(switch(
        input$claimed_taxon_sort,
        coverage = selected |>
          arrange(pair_scorable_fraction, desc(n_centroids), taxon),
        sample = selected |> arrange(desc(n_centroids), taxon),
        taxon = selected |> arrange(taxon),
        selected |> arrange(desc(pair_median_penalty), taxon)
      ))
    }
    selected_orders <- switch(
      input$detail_order,
      Acari = c(
        "Ixodida", "Mesostigmata", "Sarcoptiformes", "Trombidiformes",
        "Opilioacarida", "Holothyrida"
      ),
      Collembola = c(
        "Entomobryomorpha", "Poduromorpha", "Symphypleona", "Neelipleona"
      ),
      input$detail_order
    )

    selected <- switch(
      input$claimed_taxon_rank,
      family = claimed_primer_families |>
        filter(pair_id == input$detail_pair, order %in% selected_orders) |>
        mutate(taxon = family, context = order),
      genus = claimed_primer_genera |>
        filter(pair_id == input$detail_pair, order %in% selected_orders) |>
        mutate(taxon = genus, context = paste(order, family, sep = " · ")),
      subfamily = claimed_primer_subfamilies |>
        filter(pair_id == input$detail_pair, order %in% selected_orders) |>
        mutate(
          taxon = subfamily,
          context = paste(order, family, sep = " · ")
        )
    )

    switch(
      input$claimed_taxon_sort,
      coverage = selected |>
        arrange(pair_scorable_fraction, desc(n_centroids), taxon),
      sample = selected |> arrange(desc(n_centroids), taxon),
      taxon = selected |> arrange(taxon),
      selected |> arrange(desc(pair_median_penalty), taxon)
    )
  })

  claimed_selected_orders <- reactive({
    req(input$detail_order)
    switch(
      input$detail_order,
      Acari = c(
        "Ixodida", "Mesostigmata", "Sarcoptiformes", "Trombidiformes",
        "Opilioacarida", "Holothyrida"
      ),
      Collembola = c(
        "Entomobryomorpha", "Poduromorpha", "Symphypleona", "Neelipleona"
      ),
      input$detail_order
    )
  })

  claimed_taxon_display <- reactive({
    x <- claimed_taxon_data()
    if (!nrow(x)) return(data.frame())
    usable <- claimed_reference_usable()
    targeted_active <- isTRUE(targeted_zbj_active())
    evidence <- if (targeted_active) {
      paste0(
        format(x$n_sequences, big.mark = ","), " complete COX1 · ",
        format(x$n_species, big.mark = ","), " species"
      )
    } else {
      paste0(
        format(x$n_centroids, big.mark = ","), " centroids · ",
        format(x$source_records_represented, big.mark = ","), " records"
      )
    }
    unavailable <- rep("—", nrow(x))
    data.frame(
      Taxon = x$taxon,
      Context = x$context,
      Evidence = evidence,
      `Both sites` = paste0(
        x$n_pair_scorable, "/", x$n_centroids, " (",
        round(x$pair_scorable_fraction * 100, 1), "%)"
      ),
      `Forward median` = if (usable) {
        round(x$forward_median_penalty, 1)
      } else {
        unavailable
      },
      `Reverse median` = if (usable) {
        round(x$reverse_median_penalty, 1)
      } else {
        unavailable
      },
      `Pair median` = if (usable) {
        round(x$pair_median_penalty, 1)
      } else {
        unavailable
      },
      `Pair p90` = if (usable) {
        round(x$pair_p90_penalty, 1)
      } else {
        unavailable
      },
      `Terminal-3′ mismatch` = if (usable) {
        round(x$terminal_3_mismatch_fraction_among_scorable * 100, 1)
      } else {
        unavailable
      },
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })

  claimed_selected_taxon <- reactive({
    selected <- input$claimed_taxon_table_rows_selected
    x <- claimed_taxon_data()
    if (length(selected) != 1L || selected < 1L || selected > nrow(x)) {
      return(NULL)
    }
    x[selected, , drop = FALSE]
  })

  claimed_raw_scores <- reactive({
    req(input$detail_pair)
    if (isTRUE(custom_detail_active())) {
      path <- custom_detail_record()$expanded_scores_path
      validate(need(
        !is.null(path) && file.exists(path),
        "Expanded custom-primer sequence evidence is no longer available."
      ))
      return(readRDS(path))
    }
    if (input$detail_pair %in% c("ZBJ_ART", "ZBJ_ART_DEG")) {
      return(
        targeted_zbj_sequence_scores |>
          filter(pair_id == input$detail_pair)
      )
    }
    path <- file.path(
      "data", "derived",
      paste0("claimed_", tolower(input$detail_pair), "_centroid_scores.csv.gz")
    )
    validate(need(file.exists(path), paste("Raw sequence evidence is missing:", path)))
    read_csv(path, show_col_types = FALSE)
  }) |>
    bindCache(input$detail_pair)

  claimed_selected_sequence_data <- reactive({
    selected <- claimed_selected_taxon()
    if (is.null(selected)) return(tibble())
    x <- claimed_raw_scores()
    if (isTRUE(targeted_zbj_active())) {
      x <- x |>
        filter(source_group == input$detail_order)
    } else {
      x <- x |>
        filter(order %in% claimed_selected_orders())
    }
    x <- switch(
      input$claimed_taxon_rank,
      family = x |> filter(family == selected$taxon),
      genus = x |>
        filter(
          genus == selected$taxon,
          family == selected$family
        ),
      subfamily = x |>
        filter(
          subfamily == selected$taxon,
          family == selected$family
        )
    )
    if ("cluster_size" %in% names(x)) {
      x |> arrange(desc(cluster_size), organism_label, accession)
    } else {
      x |> arrange(species, scientific_name, accession)
    }
  })

  claimed_sequence_display <- reactive({
    x <- claimed_selected_sequence_data()
    if (!nrow(x)) return(data.frame())
    pick <- function(column, default = NA) {
      if (column %in% names(x)) x[[column]] else rep(default, nrow(x))
    }
    targeted_active <- isTRUE(targeted_zbj_active())
    organism <- if (targeted_active) {
      coalesce(
        as.character(pick("species", NA_character_)),
        as.character(pick("scientific_name", NA_character_)),
        as.character(pick("ncbi_title", NA_character_)),
        as.character(pick("accession", NA_character_))
      )
    } else {
      coalesce(
        as.character(pick("scientific_name_ncbi", NA_character_)),
        as.character(pick("organism_label", NA_character_)),
        as.character(pick("accession", NA_character_))
      )
    }
    accession <- as.character(pick("accession", ""))
    accession_link <- paste0(
      "<a href=\"https://www.ncbi.nlm.nih.gov/nuccore/",
      vapply(accession, URLencode, character(1), reserved = TRUE),
      "\" target=\"_blank\" rel=\"noopener\">",
      accession,
      "</a>"
    )
    pair_penalty <- if (claimed_reference_usable()) {
      round(as.numeric(pick("pair_penalty", NA_real_)), 2)
    } else {
      rep("—", nrow(x))
    }
    data.frame(
      Accession = accession_link,
      Organism = organism,
      Family = as.character(pick("family", NA_character_)),
      Genus = as.character(pick("genus", NA_character_)),
      `Records represented` = if (targeted_active) {
        rep(1L, nrow(x))
      } else {
        as.integer(pick("cluster_size", 1L))
      },
      `Both sites` = ifelse(
        pick("pair_scorable", FALSE) %in% TRUE,
        "yes",
        "no"
      ),
      `Pair penalty` = pair_penalty,
      `F mismatches` = as.integer(pick("forward_mismatch_count", NA_integer_)),
      `R mismatches` = as.integer(pick("reverse_mismatch_count", NA_integer_)),
      `Terminal-3′ mismatches` =
        as.integer(pick("forward_terminal_3_mismatches", 0L)) +
        as.integer(pick("reverse_terminal_3_mismatches", 0L)),
      `Forward binding sequence (primer-oriented 5′→3′)` =
        as.character(pick(
          "forward_binding_sequence_primer_oriented",
          NA_character_
        )),
      `Reverse binding sequence (primer-oriented 5′→3′)` =
        as.character(pick(
          "reverse_binding_sequence_primer_oriented",
          NA_character_
        )),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })

  output$claimed_primer_metrics <- renderUI({
    overall <- claimed_pair_overall()
    scope <- claimed_selected_scope()
    if (!nrow(overall)) {
      return(div(class = "callout", "No expanded-reference summary is available."))
    }
    scope_n <- if (nrow(scope)) scope$n_centroids else 0
    scope_available <- if (nrow(scope)) {
      paste0(
        scope$n_pair_scorable, "/", scope$n_centroids, " (",
        round(scope$pair_scorable_fraction * 100, 1), "%)"
      )
    } else {
      "not represented"
    }
    usable <- claimed_reference_usable()
    penalty_text <- if (usable) {
      if (nrow(scope) && is.finite(scope$pair_median_penalty)) {
        round(scope$pair_median_penalty, 1)
      } else {
        "NA"
      }
    } else {
      "—"
    }
    reference_value <- if (isTRUE(targeted_zbj_active())) {
      format(scope$n_sequences, big.mark = ",")
    } else {
      format(overall$n_centroids, big.mark = ",")
    }
    reference_label <- if (isTRUE(targeted_zbj_active())) {
      "complete COX1 sequences"
    } else {
      "author COI centroids"
    }
    div(
      class = "metric-row",
      div(
        class = "metric",
        tags$b(reference_value),
        span(reference_label)
      ),
      div(
        class = "metric",
        tags$b(format(scope_n, big.mark = ",")),
        span(paste(input$detail_order, "centroids"))
      ),
      div(
        class = "metric",
        tags$b(scope_available),
        span("selected-group sites present")
      ),
      div(
        class = "metric",
        tags$b(penalty_text),
        span("selected-group median penalty")
      ),
      div(
        class = "metric",
        tags$b(if (usable) "usable" else "blocked"),
        span("reference suitability")
      )
    )
  })

  output$claimed_primer_interpretation <- renderUI({
    overall <- claimed_pair_overall()
    claim <- primer_target_claims |>
      filter(pair_id == input$detail_pair) |>
      slice(1)
    scope <- claimed_selected_scope()
    if (!nrow(overall)) return(NULL)

    if (isTRUE(targeted_zbj_active())) {
      targeted <- targeted_zbj_scope()
      comparison_pair <- if (input$detail_pair == "ZBJ_ART") {
        "ZBJ_ART_DEG"
      } else {
        "ZBJ_ART"
      }
      comparison <- targeted_zbj_overall |>
        filter(
          pair_id == comparison_pair,
          source_group == input$detail_order
        ) |>
        slice(1)
      manifest <- targeted_zbj_manifest |>
        filter(source_group == input$detail_order) |>
        slice(1)
      comparison_label <- if (comparison_pair == "ZBJ_ART_DEG") {
        "degenerate variant"
      } else {
        "original ZBJ"
      }
      return(div(
        class = "callout mb-3",
        strong("Site-complete ZBJ reference. "),
        input$detail_order, " is represented by ", targeted$n_sequences,
        " complete COX1 sequences from ", targeted$n_species,
        " species; ", targeted$n_pair_scorable, " (",
        round(targeted$pair_scorable_fraction * 100, 1),
        "%) contain both scorable sites. The median mismatch penalty is ",
        round(targeted$pair_median_penalty, 1), " versus ",
        round(comparison$pair_median_penalty, 1), " for the ",
        comparison_label, ". ",
        if (
          input$detail_pair == "ZBJ_ART_DEG"
        ) {
          paste0(
            "The lower in-silico mismatch score must not be called better PCR: ",
            "Elbrecht et al. (2019) observed reduced amplification efficiency ",
            "for this more-degenerate variant. "
          )
        },
        "The complete mitochondrial records were sampled reproducibly from ",
        manifest$ncbi_records_found, " available NCBI records. Family and genus ",
        "rows reveal which lineages drive the group result; raw sequence counts ",
        "are not a balanced abundance sample."
      ))
    }

    if (!isTRUE(overall$reference_suitable)) {
      return(div(
        class = "callout taxon-alert mb-3",
        strong("Reference-site failure — not primer failure. "),
        overall$reference_limitation_reason,
        " The sparse penalty subset is withheld from biological interpretation. ",
        "The exact available and unavailable sequences remain inspectable below, ",
        "but a site-complete, taxonomically balanced COI reference is required ",
        "before treating these penalties as lineage performance."
      ))
    }

    claim_text <- if (isTRUE(custom_detail_active())) {
      "Session custom primer pair"
    } else if (nrow(claim)) {
      claim$claim_label
    } else {
      "Expanded reference result"
    }
    selected_text <- if (nrow(scope) && is.finite(scope$pair_median_penalty)) {
      paste0(
        input$detail_order, " has median penalty ",
        round(scope$pair_median_penalty, 1), " across ",
        scope$n_pair_scorable, " scorable centroids; ",
        round(scope$terminal_3_mismatch_fraction_among_scorable * 100, 1),
        "% have at least one terminal-3′ mismatch."
      )
    } else {
      paste0(input$detail_order, " is not represented in this reference.")
    }

    if (
      input$detail_pair %in% c("NOSPID", "NOSPI2_LAURELIN") &&
        identical(input$detail_order, "Araneae") &&
        nrow(scope)
    ) {
      selected_text <- paste0(
        selected_text,
        " This strong host-order mismatch pattern is consistent with spider ",
        "suppression, but the family/genus table is needed to find escape lineages."
      )
    }

    div(
      class = "callout mb-3",
      strong(claim_text, ". "),
      selected_text,
      if (nrow(claim) && claim$pair_validation_status == "hybrid_unvalidated_pair") {
        paste0(
          " This exact hybrid pair is unvalidated and does not inherit full-pair ",
          "support from the two component-primer papers."
        )
      },
      " PrimerMiner penalties are relative mismatch-risk scores, not PCR probabilities."
    )
  })

  output$claimed_phylogeny_group_table <- renderUI({
    req(input$detail_pair, input$detail_order)
    group_type <- switch(
      input$detail_order,
      Lepidoptera = "Lepidoptera study lens",
      Acari = "Composite clade",
      Collembola = "Composite clade",
      NULL
    )
    if (is.null(group_type)) return(NULL)
    if (isTRUE(targeted_zbj_active())) {
      x <- targeted_zbj_study_groups |>
        filter(
          pair_id == input$detail_pair,
          source_group == input$detail_order
        ) |>
        mutate(
          display_group = study_group,
          n_centroids = n_sequences,
          group_caveat = case_when(
            study_group == "Microlepidoptera (operational grade)" ~
              "Operational and non-monophyletic; use this as a study lens.",
            study_group == "Macroheterocera (macro-moth core)" ~
              "Phylogenetic macro-moth core rather than a body-size cutoff.",
            study_group == "Butterflies (Papilionoidea)" ~
              "Shown separately instead of forcing butterflies into a moth size bin.",
            study_group == "Boundary / convention-sensitive moths" ~
              "Placement in a micro/macro bin varies by convention.",
            study_group == "Acari" ~
              "Composite arachnid clade; not a single modern order.",
            study_group == "Collembola" ~
              "Composite hexapod clade; not a single modern order.",
            TRUE ~ ""
          )
        )
    } else {
      x <- claimed_primer_phylogeny_groups |>
        filter(
          pair_id == input$detail_pair,
          display_group_type == group_type
        )
      if (input$detail_order %in% c("Acari", "Collembola")) {
        x <- x |> filter(display_group == input$detail_order)
      }
    }
    if (!nrow(x)) return(NULL)
    usable <- claimed_reference_usable()
    rows <- lapply(seq_len(nrow(x)), function(i) {
      row <- x[i, ]
      tags$tr(
        tags$td(strong(row$display_group)),
        tags$td(paste0(
          format(row$n_centroids, big.mark = ","), " centroids"
        )),
        tags$td(paste0(
          row$n_pair_scorable, "/", row$n_centroids, " (",
          round(row$pair_scorable_fraction * 100, 1), "%)"
        )),
        tags$td(if (
          usable && is.finite(row$pair_median_penalty)
        ) round(row$pair_median_penalty, 1) else "—"),
        tags$td(row$group_caveat)
      )
    })
    tagList(
      h5(if (identical(input$detail_order, "Lepidoptera")) {
        "Micro-/macrolepidoptera study lens"
      } else {
        paste(input$detail_order, "composite reference")
      }),
      div(
        class = "scroll-table mb-3",
        tags$table(
          class = "table-atlas",
          tags$thead(tags$tr(lapply(
            c("Group", "Evidence", "Both sites", "Pair median", "Definition note"),
            tags$th
          ))),
          tags$tbody(rows)
        )
      )
    )
  })

  output$claimed_taxon_notice <- renderUI({
    if (!nrow(claimed_taxon_data())) {
      return(div(
        class = "callout",
        "No rows are available at this rank. Try families or a represented group."
      ))
    }
    if (!claimed_reference_usable()) {
      return(div(
        class = "callout taxon-alert mb-2",
        strong("Why the penalty columns contain “—”. "),
        "The selected reference does not contain both primer-binding regions for ",
        "enough sequences. The app therefore shows binding-site availability but ",
        "does not turn a sparse, biased subset into a biological penalty estimate. ",
        "This is unavailable evidence, not zero penalty and not PCR failure."
      ))
    }
    p(
      class = "tiny mb-2",
      "Search globally, filter individual columns, or sort by any header. Select ",
      "one row to open the exact scored sequences and organisms below."
    )
  })

  output$claimed_taxon_table <- DT::renderDT({
    x <- claimed_taxon_display()
    if (!nrow(x)) {
      return(DT::datatable(
        data.frame(`Taxonomic summary` = character()),
        rownames = FALSE,
        selection = "none",
        options = list(
          dom = "t",
          language = list(
            emptyTable = "Choose a mapped primer pair."
          )
        )
      ))
    }
    DT::datatable(
      x,
      rownames = FALSE,
      filter = "top",
      selection = list(mode = "single", target = "row"),
      class = "compact stripe hover order-column",
      options = list(
        pageLength = 20,
        lengthMenu = c(10, 20, 50, 100),
        scrollX = TRUE,
        scrollY = "520px",
        scrollCollapse = TRUE,
        autoWidth = TRUE,
        order = list(),
        search = list(regex = FALSE, caseInsensitive = TRUE),
        language = list(
          search = "Search all columns:",
          select = list(rows = "%d row selected")
        )
      )
    )
  }, server = TRUE)

  output$claimed_sequence_intro <- renderUI({
    selected <- claimed_selected_taxon()
    if (is.null(selected)) {
      return(div(
        class = "callout",
        "Select a family, genus, or subfamily row above to inspect the exact ",
        "reference sequences used in its score."
      ))
    }
    x <- claimed_selected_sequence_data()
    if (!nrow(x)) {
      return(div(
        class = "callout taxon-alert",
        "No sequence-level rows could be recovered for this selection."
      ))
    }
    if (isTRUE(targeted_zbj_active())) {
      return(div(
        class = "callout taxon-info",
        strong(paste0(selected$taxon, ". ")),
        paste0(
          nrow(x), " versioned complete-COX1 record",
          if (nrow(x) == 1L) "" else "s",
          " were scored. Each accession opens the exact NCBI record; the FASTA ",
          "download contains the complete COX1 coding sequences used for alignment."
        )
      ))
    }
    represented <- sum(x$cluster_size, na.rm = TRUE)
    div(
      class = "callout taxon-info",
      strong(paste0(selected$taxon, ". ")),
      paste0(
        nrow(x), " exact 99.5%-clustered centroid sequence",
        if (nrow(x) == 1L) "" else "s",
        " were scored, representing ", format(represented, big.mark = ","),
        " source records. Gurten et al. Supplement 5 exposes each scored centroid ",
        "accession and cluster size, but not the identifiers of every cluster member. ",
        "Therefore the table and FASTA identify the exact scored representatives; ",
        "the complete membership of all represented records cannot be reconstructed ",
        "from that supplement."
      )
    )
  })

  output$claimed_sequence_downloads <- renderUI({
    if (is.null(claimed_selected_taxon())) return(NULL)
    div(
      class = "sequence-evidence-actions",
      downloadButton(
        "download_claimed_sequence_evidence",
        "Download evidence CSV",
        class = "btn-sm"
      ),
      downloadButton(
        "download_claimed_sequences_fasta",
        "Download exact FASTA",
        class = "btn-sm"
      )
    )
  })

  output$claimed_sequence_table <- DT::renderDT({
    x <- claimed_sequence_display()
    if (!nrow(x)) {
      return(DT::datatable(
        data.frame(`Exact sequence evidence` = character()),
        rownames = FALSE,
        selection = "none",
        options = list(
          dom = "t",
          language = list(
            emptyTable = "Select a summary row to load exact sequence evidence."
          )
        )
      ))
    }
    table <- DT::datatable(
      x,
      rownames = FALSE,
      filter = "top",
      escape = FALSE,
      selection = "none",
      class = "compact stripe hover order-column",
      options = list(
        pageLength = 15,
        lengthMenu = c(10, 15, 25, 50, 100),
        scrollX = TRUE,
        scrollY = "380px",
        scrollCollapse = TRUE,
        autoWidth = TRUE,
        order = list(list(4, "desc"), list(1, "asc")),
        search = list(regex = FALSE, caseInsensitive = TRUE),
        language = list(search = "Search exact evidence:")
      )
    )
    DT::formatStyle(
      table,
      columns = c(
        "Forward binding sequence (primer-oriented 5′→3′)",
        "Reverse binding sequence (primer-oriented 5′→3′)"
      ),
      fontFamily = "ui-monospace, SFMono-Regular, Menlo, monospace",
      whiteSpace = "nowrap"
    )
  })

  output$claimed_primer_sources <- renderUI({
    claim <- primer_target_claims |>
      filter(pair_id == input$detail_pair) |>
      slice(1)
    pair_rows <- order_lens_primers() |>
      filter(pair_id == input$detail_pair)
    keys <- if (nrow(claim)) {
      trimws(str_split(claim$source_keys, "\\|")[[1]])
    } else {
      unique(pair_rows$source_key)
    }
    keys <- unique(c(keys, "beeprime2026_suppl5"))
    if (identical(input$detail_order, "Lepidoptera")) {
      keys <- unique(c(keys, "mitter2017", "kristensen2007"))
    }
    if (isTRUE(targeted_zbj_active())) {
      keys <- unique(c(keys, "ncbi_nucleotide", "rentrez2017"))
    }
    source_rows <- citations |> filter(key %in% keys)
    div(
      h5(if (nrow(claim)) {
        "Claim and validation sources"
      } else {
        "Primer and expanded-reference sources"
      }),
      p(if (isTRUE(custom_detail_active())) {
        paste0(
          "The primer sequences were entered in this session. Their coordinates ",
          "were alignment-derived on NC_001322.1; the taxonomic drill-down was ",
          "calculated on the Gurten et al. clustered COI reference."
        )
      } else if (nrow(claim)) {
        claim$claim_note
      } else {
        paste0(
          "Primer provenance comes from the curated primer library. Family, genus, ",
          "accession, and binding-site evidence comes from the shared Gurten et al. ",
          "67,352-centroid reference."
        )
      }),
      tagList(lapply(seq_len(nrow(source_rows)), function(i) {
        div(
          class = "cite-card",
          strong(source_rows$short_citation[i]),
          div(source_rows$title[i]),
          citation_link(source_rows[i, ], "Open source ↗")
        )
      }))
    )
  })

  output$download_claimed_taxa <- downloadHandler(
    filename = function() {
      paste0(
        input$detail_pair, "_", input$detail_order, "_",
        input$claimed_taxon_rank, "_expanded_reference.csv"
      )
    },
    content = function(file) {
      write.csv(claimed_taxon_data(), file, row.names = FALSE, na = "")
    }
  )

  selected_sequence_filename <- reactive({
    selected <- claimed_selected_taxon()
    req(selected)
    taxon_slug <- gsub(
      "(^_+|_+$)",
      "",
      gsub("[^A-Za-z0-9]+", "_", selected$taxon)
    )
    paste(
      input$detail_pair,
      input$detail_order,
      input$claimed_taxon_rank,
      taxon_slug,
      sep = "_"
    )
  })

  output$download_claimed_sequence_evidence <- downloadHandler(
    filename = function() {
      paste0(selected_sequence_filename(), "_sequence_evidence.csv")
    },
    content = function(file) {
      x <- claimed_selected_sequence_data()
      req(nrow(x))
      write.csv(x, file, row.names = FALSE, na = "")
    }
  )

  output$download_claimed_sequences_fasta <- downloadHandler(
    filename = function() {
      paste0(selected_sequence_filename(), "_exact_sequences.fasta")
    },
    content = function(file) {
      selected_rows <- claimed_selected_sequence_data()
      req(nrow(selected_rows))
      source_path <- if (isTRUE(targeted_zbj_active())) {
        file.path(
          "data", "external", "targeted_zbj", "raw",
          paste0(input$detail_order, "_complete_COX1.fasta")
        )
      } else {
        file.path(
          "data", "external", "gurten2026", "ClusteredReferences.fasta"
        )
      }
      validate(need(file.exists(source_path), "The source FASTA is unavailable."))
      sequences <- parse_fasta(source_path)
      fasta_accessions <- sub(
        "^_R_",
        "",
        sub(" .*", "", names(sequences))
      )
      keep <- fasta_accessions %in% selected_rows$accession
      validate(need(
        any(keep),
        "None of the selected accessions could be recovered from the source FASTA."
      ))
      sequences <- sequences[keep]
      lines <- unlist(lapply(seq_along(sequences), function(i) {
        sequence <- sequences[[i]]
        starts <- seq(1L, nchar(sequence), by = 80L)
        c(
          paste0(">", names(sequences)[i]),
          substring(
            sequence,
            starts,
            pmin(starts + 79L, nchar(sequence))
          )
        )
      }), use.names = FALSE)
      writeLines(lines, file, useBytes = TRUE)
    }
  )

  output$download_detail <- downloadHandler(
    filename = function() {
      paste0(input$detail_order, "_", input$detail_pair, "_PrimerMiner_components.csv")
    },
    content = function(file) {
      detail <- order_lens_position_scores() |>
        filter(order == input$detail_order, pair_id == input$detail_pair)
      write.csv(detail, file, row.names = FALSE, na = "")
    }
  )

  observeEvent(input$open_lineage, {
    updateSelectInput(session, "lineage_marker", selected = "COI")
    updateSelectizeInput(session, "lineage_pair", selected = input$detail_pair)
    updateSelectizeInput(session, "lineage_target", selected = input$detail_order)
    updateNavbarPage(session, "main_nav", selected = "Lineage Explorer")
  })

  observeEvent(input$lineage_marker, {
    available <- catalog_pairs |> filter(marker_id == input$lineage_marker)
    updateSelectizeInput(
      session, "lineage_pair",
      choices = setNames(available$pair_id, available$pair_label),
      selected = available$pair_id[1], server = TRUE
    )
    targets <- if (input$lineage_marker == "COI") sort(unique(claimed_primer_orders$order)) else "Fungi"
    updateSelectizeInput(session, "lineage_target", choices = c("All", targets), selected = "All", server = TRUE)
  }, ignoreInit = FALSE)

  lineage_geography <- reactive({
    if (input$lineage_marker != "COI") return(tibble())
    release_table("geography", "data/derived/reference_geography.csv") |>
      mutate(
        country_or_territory = na_if(country_or_territory, ""),
        locality = na_if(locality, "")
      )
  }) |> bindCache(input$lineage_marker)

  observe({
    geo <- lineage_geography()
    countries <- if (nrow(geo)) sort(unique(na.omit(geo$country_or_territory))) else character()
    updateSelectizeInput(session, "lineage_country", choices = c("All", countries), selected = "All", server = TRUE)
  })
  observeEvent(input$lineage_country, {
    geo <- lineage_geography()
    if (nrow(geo) && !is.null(input$lineage_country) && input$lineage_country != "All") {
      geo <- geo |> filter(country_or_territory == input$lineage_country)
    }
    localities <- if (nrow(geo)) sort(unique(na.omit(geo$locality))) else character()
    updateSelectizeInput(session, "lineage_locality", choices = c("All", localities), selected = "All", server = TRUE)
  }, ignoreInit = FALSE)

  lineage_scope <- reactive({
    if (input$lineage_marker != "COI") return(list(taxonomy = tibble(), geography = tibble()))
    taxonomy <- expanded_reference_taxonomy
    geo <- lineage_geography()
    joined <- taxonomy |> left_join(geo, by = "accession")
    selected <- joined
    if (!is.null(input$lineage_target) && input$lineage_target != "All") selected <- selected |> filter(order == input$lineage_target)
    if (!is.null(input$lineage_country) && input$lineage_country != "All") selected <- selected |> filter(country_or_territory == input$lineage_country)
    if (!is.null(input$lineage_locality) && input$lineage_locality != "All") selected <- selected |> filter(locality == input$lineage_locality)
    list(taxonomy = selected, geography = joined)
  })

  output$lineage_denominator <- renderUI({
    if (input$lineage_marker != "COI") return(div(class = "callout tiny", "Fungal ITS geography and lineage partitions are scheduled for the expanded reference release."))
    all_rows <- lineage_scope()$geography
    selected <- lineage_scope()$taxonomy
    located <- sum(!is.na(all_rows$country_or_territory))
    div(
      class = "callout tiny",
      strong(format(located, big.mark = ","), " / ", format(nrow(all_rows), big.mark = ","), " eligible centroids have reference geography."),
      br(), format(nrow(all_rows) - located, big.mark = ","), " lack a usable location. Current filters retain ", format(nrow(selected), big.mark = ","), " centroids."
    )
  })

  lineage_summary_table <- function(source, rank) {
    if (input$lineage_marker != "COI" || is.null(input$lineage_pair)) return(tibble())
    x <- source |> filter(pair_id == input$lineage_pair)
    if (!is.null(input$lineage_target) && input$lineage_target != "All") x <- x |> filter(order == input$lineage_target)
    selected_accessions <- lineage_scope()$taxonomy$accession
    geo_filtered <- (!is.null(input$lineage_country) && input$lineage_country != "All") ||
      (!is.null(input$lineage_locality) && input$lineage_locality != "All")
    if (geo_filtered) {
      counts <- lineage_scope()$taxonomy |> count(.data[[rank]], name = "located_centroids")
      names(counts)[1] <- rank
      x <- x |> inner_join(counts, by = rank)
    }
    x
  }

  output$lineage_overview <- renderUI({
    marker <- marker_registry |> filter(marker_id == input$lineage_marker) |> slice(1)
    if (input$lineage_marker != "COI") {
      claims <- catalog_claims |> filter(marker_id == input$lineage_marker, pair_id == input$lineage_pair)
      return(tagList(
        h3(marker$display_name, " pilot"),
        div(class = "callout", strong("Screening reference only. "), "Primer coordinates are aligned to ", marker$coordinate_reference, ". Family-to-sequence claims are intentionally unavailable until a fungal expanded reference is released."),
        if (nrow(claims)) tags$ul(lapply(claims$claim_note, tags$li))
      ))
    }
    tagList(
      h3("Reference fit by target group"),
      p("Penalty is calculated only when both binding sites are scorable. Site availability and mismatch fit remain separate."),
      DTOutput("lineage_overview_table"),
      uiOutput("lineage_geography_interpretation"),
      DTOutput("lineage_geography_stratified")
    )
  })
  output$lineage_overview_table <- renderDT({
    summary <- claimed_primer_orders |> filter(pair_id == input$lineage_pair)
    if (!is.null(input$lineage_target) && input$lineage_target != "All") summary <- summary |> filter(order == input$lineage_target)
    datatable(summary, options = list(pageLength = 12, scrollX = TRUE), filter = "top", rownames = FALSE)
  }, server = TRUE)

  output$lineage_families <- renderDT({
    x <- lineage_summary_table(claimed_primer_families, "family")
    datatable(x, filter = "top", selection = "single", rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
  }, server = TRUE)
  output$lineage_subfamilies_ui <- renderUI({
    if (input$lineage_marker != "COI") return(div(class = "callout", "Subfamily evidence is not available for this marker release."))
    DTOutput("lineage_subfamilies")
  })
  output$lineage_subfamilies <- renderDT({
    datatable(lineage_summary_table(claimed_primer_subfamilies, "subfamily"), filter = "top", selection = "single", rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
  }, server = TRUE)
  output$lineage_genera_ui <- renderUI({
    if (input$lineage_marker != "COI") return(div(class = "callout", "Genus evidence is not available for this marker release."))
    DTOutput("lineage_genera")
  })
  output$lineage_genera <- renderDT({
    datatable(lineage_summary_table(claimed_primer_genera, "genus"), filter = "top", selection = "single", rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
  }, server = TRUE)

  lineage_raw_scores <- reactive({
    req(input$lineage_marker == "COI", input$lineage_pair)
    path <- file.path("data", "derived", paste0("claimed_", tolower(input$lineage_pair), "_centroid_scores.csv.gz"))
    artifact_name <- paste0("exact_", tolower(input$lineage_pair))
    artifact <- if (!is.null(coi_release_state)) coi_release_state$manifest$artifacts[[artifact_name]] else NULL
    if (!is.null(artifact)) {
      value <- tryCatch(as_tibble(atlas_read_artifact(artifact)), error = function(e) NULL)
      if (!is.null(value)) return(value)
    }
    if (!file.exists(path)) return(tibble())
    read_csv(path, show_col_types = FALSE)
  }) |> bindCache(input$lineage_pair)

  output$lineage_geography_interpretation <- renderUI({
    if (input$lineage_marker != "COI") return(NULL)
    if (is.null(input$lineage_country) || input$lineage_country == "All") {
      return(div(class = "callout tiny", "Choose a reference country or territory to compare its mismatch pattern with the global reference after stratifying by order."))
    }
    div(
      class = "callout tiny",
      strong("Taxon-stratified reference-geography check. "),
      "For each order represented at the selected location, the table compares its local median with the same order globally. Remaining differences are descriptive reference patterns—not regional amplification probabilities—and may reflect sparse or biased sampling."
    )
  })
  output$lineage_geography_stratified <- renderDT({
    req(input$lineage_marker == "COI", input$lineage_country != "All")
    x <- lineage_raw_scores() |> left_join(lineage_geography(), by = "accession")
    summarize_order <- function(rows, prefix) {
      rows |>
        filter(!is.na(order), nzchar(order)) |>
        group_by(order) |>
        summarise(
          n = n(), n_pair_scorable = sum(pair_scorable %in% TRUE),
          pair_scorable_fraction = mean(pair_scorable %in% TRUE),
          pair_median_penalty = if (any(pair_scorable %in% TRUE)) median(pair_penalty[pair_scorable %in% TRUE], na.rm = TRUE) else NA_real_,
          .groups = "drop"
        ) |>
        rename_with(~paste0(prefix, .x), -order)
    }
    global <- summarize_order(x, "global_")
    local <- x |> filter(country_or_territory == input$lineage_country)
    if (!is.null(input$lineage_locality) && input$lineage_locality != "All") local <- local |> filter(locality == input$lineage_locality)
    comparison <- summarize_order(local, "selected_") |>
      left_join(global, by = "order") |>
      mutate(median_difference_within_order = selected_pair_median_penalty - global_pair_median_penalty)
    datatable(comparison, filter = "top", rownames = FALSE, options = list(pageLength = 12, scrollX = TRUE))
  }, server = TRUE)

  lineage_selected_family <- reactive({
    row <- input$lineage_families_rows_selected
    table <- lineage_summary_table(claimed_primer_families, "family")
    if (length(row) && row <= nrow(table)) table$family[row] else NULL
  })
  lineage_sequence_data <- reactive({
    if (input$lineage_marker != "COI") return(tibble())
    x <- lineage_raw_scores()
    if (!nrow(x)) return(x)
    if (!is.null(input$lineage_target) && input$lineage_target != "All" && "order" %in% names(x)) x <- x |> filter(order == input$lineage_target)
    family <- lineage_selected_family()
    if (!is.null(family) && "family" %in% names(x)) x <- x |> filter(.data$family == family)
    x <- x |> left_join(lineage_geography(), by = "accession")
    if (!is.null(input$lineage_country) && input$lineage_country != "All") x <- x |> filter(country_or_territory == input$lineage_country)
    if (!is.null(input$lineage_locality) && input$lineage_locality != "All") x <- x |> filter(locality == input$lineage_locality)
    x
  })
  output$lineage_sequence_note <- renderUI({
    if (input$lineage_marker != "COI") return(div(class = "callout", "Exact fungal sequence evidence will appear only after the expanded-reference artifact is released."))
    family <- lineage_selected_family()
    div(class = "tiny", if (is.null(family)) "Select a family row to restrict exact sequence evidence." else paste("Showing exact evidence for", family), " Sequence artifacts load only when this tab is opened.")
  })
  output$lineage_sequences <- renderDT({
    x <- lineage_sequence_data()
    keep <- intersect(c("accession", "organism_label", "scientific_name_ncbi", "order", "family", "subfamily", "genus", "cluster_size", "pair_scorable", "pair_penalty", "forward_binding_sequence_primer_oriented", "reverse_binding_sequence_primer_oriented", "country_or_territory", "locality"), names(x))
    datatable(x[, keep, drop = FALSE], filter = "top", rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
  }, server = TRUE)
  output$lineage_sources <- renderUI({
    claims <- catalog_claims |> filter(marker_id == input$lineage_marker, pair_id == input$lineage_pair)
    marker <- marker_registry |> filter(marker_id == input$lineage_marker) |> slice(1)
    tagList(
      h4("Coordinate and reference source"), tags$a(href = marker$source_url, target = "_blank", marker$source_name),
      h4("Publication claims"),
      if (!nrow(claims)) p("No target-specific publication claim is registered; primer citations remain attached to the oligos.") else tags$ul(lapply(seq_len(nrow(claims)), function(i) tags$li(claims$claim_label[i], ": ", claims$claim_note[i]))),
      div(class = "callout tiny", "Original third-party licensing and attribution apply. Annual software and redistributable-data snapshots are prepared for Zenodo DOI deposition.")
    )
  })

  thermal <- reactive({
    pair_thermal_summary(primers, input$sodium, input$ta_offset)
  })

  output$thermal_metrics <- renderUI({
    x <- thermal()
    div(
      class = "metric-row",
      div(class = "metric", tags$b(nrow(x)), span("individual primers")),
      div(class = "metric", tags$b(sum(x$degeneracy > 1, na.rm = TRUE)), span("degenerate primers")),
      div(class = "metric", tags$b(paste0(round(min(x$ta_low)), "–", round(max(x$ta_high)), "°")), span("pair Ta starting range"))
    )
  })

  output$thermal_table <- renderUI({
    x <- thermal() |> arrange(pair_label, direction)
    header <- tags$tr(lapply(
      c("Pair", "Primer", "Dir.", "Sequence (5′→3′)", "Variants", "Tm range", "Suggested pair Ta", "Reported Ta", "Source"),
      tags$th
    ))
    rows <- lapply(seq_len(nrow(x)), function(i) {
      row <- x[i, ]
      src <- citations |> filter(key == row$source_key)
      tags$tr(
        tags$td(strong(row$pair_label)),
        tags$td(row$primer_name),
        tags$td(row$direction),
        tags$td(class = "mono", row$sequence),
        tags$td(format(row$degeneracy, scientific = FALSE, big.mark = ",")),
        tags$td(sprintf("%.1f–%.1f °C", row$tm_min, row$tm_max)),
        tags$td(sprintf("%.1f–%.1f °C", row$ta_low, row$ta_high)),
        tags$td(ifelse(is.na(row$reported_ta_c), "—", paste0(row$reported_ta_c, " °C"))),
        tags$td(citation_link(src))
      )
    })
    tags$table(class = "table-atlas", tags$thead(header), tags$tbody(rows))
  })

  output$citation_list <- renderUI({
    tagList(lapply(seq_len(nrow(citations)), function(i) {
      row <- citations[i, ]
      div(
        class = "cite-card",
        span(class = paste("pill", if (row$category == "primer") "coral"), row$category),
        strong(row$short_citation),
        div(row$title),
        div(
          class = "tiny",
          if (!is.na(row$doi) && nzchar(row$doi)) paste0("doi:", row$doi, " · ") else "",
          row$note
        ),
        citation_link(row, "Open source ↗")
      )
    }))
  })

  output$download_bib <- downloadHandler(
    filename = "coi-primer-atlas-references.bib",
    content = function(file) {
      entries <- vapply(seq_len(nrow(citations)), function(i) {
        row <- citations[i, ]
        key <- row$key
        doi_line <- if (!is.na(row$doi) && nzchar(row$doi)) paste0("  doi = {", row$doi, "},\n") else ""
        paste0(
          "@misc{", key, ",\n",
          "  title = {", row$title, "},\n",
          "  year = {", row$year, "},\n",
          doi_line,
          "  url = {", row$url, "}\n",
          "}"
        )
      }, character(1))
      writeLines(entries, file)
    }
  )
}

shinyApp(ui, server)
