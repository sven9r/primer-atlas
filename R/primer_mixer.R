# Session-local, marker-independent primer exploration on explicit raw templates.
# A sliding, full-oligo Hamming search reports sequence compatibility; it does
# not estimate amplification probability, accept indels or borrow another
# marker's coordinates. Frozen reference catalogs and scores are never edited.

.primer_mixer_text <- function(x, default = "") {
  if (is.null(x)) return(default)
  x <- as.character(x)
  x[is.na(x)] <- default
  x
}

.primer_mixer_hash <- function(x) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("The primer mixer requires digest for stable sequence identities.")
  substr(digest::digest(x, algo = "sha256", serialize = FALSE), 1L, 20L)
}

.primer_mixer_sequence <- function(sequence, custom = FALSE) {
  if (length(sequence) != 1L || is.na(sequence)) stop("Provide one nucleotide sequence.")
  original <- toupper(gsub("[[:space:]]", "", as.character(sequence)))
  if (!nzchar(original) || grepl("[^ACGTURYSWKMBDHVNI]", original)) {
    stop("Use an ungapped nucleotide sequence with IUPAC symbols only; do not include adapters, punctuation or FASTA headers in an oligo field.")
  }
  if (custom && (nchar(original) < 8L || nchar(original) > 120L)) stop("A custom oligo must contain 8 to 120 nucleotides.")
  note <- c(if (grepl("I", original, fixed = TRUE)) "Inosine is modeled as N; this is a nucleotide-model assumption, not biochemical validation.",
    if (grepl("U", original, fixed = TRUE)) "U is modeled as T.")
  list(sequence = chartr("UI", "TN", original), source_sequence = original,
    normalization_note = paste(note, collapse = " "))
}

.primer_mixer_complement <- function(sequence) {
  complement <- c(A="T", C="G", G="C", T="A", R="Y", Y="R", S="S", W="W",
    K="M", M="K", B="V", D="H", H="D", V="B", N="N")
  paste(rev(unname(complement[strsplit(sequence, "", fixed = TRUE)[[1L]]])), collapse = "")
}

.primer_mixer_bind <- function(frames) {
  frames <- Filter(function(x) is.data.frame(x) && nrow(x), frames)
  if (!length(frames)) return(data.frame())
  fields <- unique(unlist(lapply(frames, names), use.names = FALSE))
  frames <- lapply(frames, function(x) {
    for (field in setdiff(fields, names(x))) x[[field]] <- rep(NA, nrow(x))
    x[, fields, drop = FALSE]
  })
  result <- do.call(rbind, frames)
  rownames(result) <- NULL
  result
}

# The source sequence and citation are retained even when the matching model
# converts U or inosine. A custom row has no inherited assay validation.
primer_mixer_custom <- function(name, sequence, direction, marker_id) {
  if (length(name) != 1L || is.na(name) || !nzchar(trimws(name)) || nchar(name) > 100L) stop("Name the custom oligo using at most 100 characters.")
  if (length(direction) != 1L || !direction %in% c("forward", "reverse")) stop("Choose forward or reverse for the custom oligo.")
  if (length(marker_id) != 1L || is.na(marker_id) || !nzchar(marker_id)) stop("Choose a marker for the custom oligo.")
  seq <- .primer_mixer_sequence(sequence, custom = TRUE)
  data.frame(primer_id = paste0("custom_", .primer_mixer_hash(paste(marker_id, direction, trimws(name), seq$source_sequence, sep = "\036"))),
    primer_name = trimws(name), sequence = seq$sequence, direction = direction, marker_id = marker_id,
    source_sequence = seq$source_sequence, source_url = "", citation = "User-provided oligo; no published validation assigned",
    source_key = "custom", catalog_origin = "custom_session_input", target_group = "",
    target_scope = "User-defined candidate", evidence_note = "Session input. Sequence compatibility does not establish PCR performance.",
    normalization_note = seq$normalization_note, stringsAsFactors = FALSE)
}

# A second map click constructs an explicit opposite-role oligo. The original
# source citation is retained as parent provenance only, not as validation of
# the derived component or any newly assembled pair.
primer_mixer_reverse <- function(component) {
  if (!is.data.frame(component) || nrow(component)!=1L || !all(c("primer_id","primer_name","sequence","direction","marker_id") %in% names(component))) stop("Reverse one identified catalog component at a time.")
  if ("catalog_origin" %in% names(component) && component$catalog_origin %in% "derived_reverse_complement") stop("Choose the original component rather than deriving a second reverse-complement variant.")
  original <- .primer_mixer_sequence(component$sequence)
  if (!component$direction %in% c("forward","reverse")) stop("The parent component needs a forward or reverse direction.")
  child <- component
  child$parent_primer_id <- component$primer_id
  child$parent_primer_name <- component$primer_name
  child$parent_sequence <- if("source_sequence" %in% names(component)) component$source_sequence else component$sequence
  child$parent_direction <- component$direction
  child$direction <- if(component$direction=="forward") "reverse" else "forward"
  child$sequence <- .primer_mixer_complement(original$sequence)
  child$source_sequence <- child$sequence
  child$primer_name <- paste0(component$primer_name," [RC]")
  child$primer_id <- paste0("rc_",.primer_mixer_hash(paste(component$marker_id,component$primer_id,child$sequence,child$direction,sep="\036")))
  child$catalog_origin <- "derived_reverse_complement"
  parent_note <- if("normalization_note" %in% names(component)) .primer_mixer_text(component$normalization_note) else original$normalization_note
  child$normalization_note <- paste(c(parent_note,"Explicit reverse-complement candidate in the opposite role. Parent citations describe the original oligo; no empirical validation is inherited."),collapse=" ")
  child$evidence_note <- "Derived reverse-complement candidate. Retained citations establish parent provenance only; sequence compatibility does not establish PCR performance."
  child
}

# Normalize unfrozen working catalogs for a shared selector. Stable original
# catalog IDs are retained where supplied. Duplicated legacy pair membership
# does not create duplicate chemical components.
primer_mixer_catalog <- function(root = ".") {
  read <- function(path) {
    path <- file.path(root, path)
    if (!file.exists(path)) return(data.frame())
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  fields <- function(rows, field, default = "") {
    if (!field %in% names(rows)) return(rep(default, nrow(rows)))
    .primer_mixer_text(rows[[field]], default)
  }
  blocks <- list()
  legacy <- read("data/primers.csv")
  pairs <- read("data/catalog/primer_pairs.csv")
  citations <- read("data/citations.csv")
  if (nrow(legacy)) {
    marker <- if (nrow(pairs)) pairs$marker_id[match(legacy$pair_id, pairs$pair_id)] else rep("COI", nrow(legacy))
    marker[is.na(marker)] <- "COI"
    cited <- if (nrow(citations) && "key" %in% names(citations)) match(legacy$source_key, citations$key) else rep(NA_integer_, nrow(legacy))
    blocks[[length(blocks)+1L]] <- data.frame(primer_name = fields(legacy, "primer_name"), sequence = fields(legacy, "sequence"),
      direction = fields(legacy, "direction"), marker_id = marker, source_key = fields(legacy, "source_key"),
      source_url = if ("url" %in% names(citations)) fields(citations, "url")[cited] else rep("", nrow(legacy)),
      citation = if ("short_citation" %in% names(citations)) fields(citations, "short_citation")[cited] else fields(legacy, "source_key"),
      catalog_origin = "curated_atlas", target_scope = fields(legacy, "target"),
      evidence_note = "Published component provenance; a newly assembled combination is a candidate.", stringsAsFactors = FALSE)
  }
  pr2 <- read("data/catalog/pr2_18s_primers.csv")
  if (nrow(pr2)) blocks[[length(blocks)+1L]] <- data.frame(primer_id = paste0("pr2_18s_", pr2$primer_id),
    primer_name = fields(pr2, "primer_name"), sequence = fields(pr2, "sequence"), direction = fields(pr2, "direction"),
    marker_id = "18S", source_url = fields(pr2, "source_url"), citation = fields(pr2, "reference"),
    catalog_origin = "PR2_primer_2.1.1", target_scope = fields(pr2, "specificity"),
    evidence_note = paste("PR2 catalog component; applicability to the selected animal group requires template evidence.", fields(pr2, "remarks")), stringsAsFactors = FALSE)
  unite <- read("data/catalog/unite_primers.csv")
  if (nrow(unite)) blocks[[length(blocks)+1L]] <- data.frame(primer_name = fields(unite, "primer_name"),
    sequence = fields(unite, "sequence"), direction = fields(unite, "direction"), marker_id = "ITS_FUNGAL",
    source_url = fields(unite, "primary_reference_url"), citation = fields(unite, "reference"),
    catalog_origin = "UNITE_primer_catalog", target_scope = fields(unite, "target"),
    gene_locus=fields(unite,"gene_locus"),evidence_note = fields(unite, "remarks"), stringsAsFactors = FALSE)
  fish <- read("data/catalog/fish_assay_primers.csv")
  if (nrow(fish)) {
    fish$catalog_origin <- "fish_assay_catalog"
    blocks[[length(blocks)+1L]] <- fish
  }
  animal <- read("data/catalog/animal_assay_primers.csv")
  if(nrow(animal)) {
    animal$catalog_origin <- "animal_assay_catalog"
    blocks[[length(blocks)+1L]] <- animal
  }
  rows <- .primer_mixer_bind(blocks)
  if (!nrow(rows)) return(rows)
  required <- c("primer_id", "primer_name", "sequence", "direction", "marker_id", "source_sequence", "source_url", "citation", "source_key",
    "catalog_origin", "target_group", "target_scope", "evidence_note", "normalization_note")
  for (field in required) rows[[field]] <- if (field %in% names(rows)) .primer_mixer_text(rows[[field]]) else rep("", nrow(rows))
  valid <- rows$direction %in% c("forward", "reverse") & nzchar(rows$primer_name) & nzchar(rows$marker_id)
  normalized <- lapply(rows$sequence, function(x) tryCatch(.primer_mixer_sequence(x), error = function(e) NULL))
  valid <- valid & vapply(normalized, function(x) !is.null(x), logical(1))
  excluded <- rows[!valid, , drop = FALSE]
  rows <- rows[valid, , drop = FALSE]
  normalized <- normalized[valid]
  rows$source_sequence <- vapply(normalized, `[[`, character(1), "source_sequence")
  rows$sequence <- vapply(normalized, `[[`, character(1), "sequence")
  rows$normalization_note <- vapply(normalized, `[[`, character(1), "normalization_note")
  identities <- paste(rows$marker_id, rows$primer_name, rows$sequence, rows$direction, rows$source_key, sep = "\036")
  missing_id <- !nzchar(rows$primer_id)
  rows$primer_id[missing_id] <- paste0("atlas_", vapply(identities[missing_id], .primer_mixer_hash, character(1)))
  # Prefer explicitly curated fish components over duplicate legacy components.
  rows <- rows[order(rows$catalog_origin != "fish_assay_catalog", rows$marker_id, rows$primer_name), , drop = FALSE]
  identities <- paste(rows$marker_id, rows$primer_name, rows$sequence, rows$direction, rows$source_key, sep = "\036")
  # PR2 separately identifies equal chemical sequences reported by different
  # papers. Keep those identities and citations, rather than dropping records.
  identified_source <- rows$catalog_origin %in% c("PR2_primer_2.1.1","fish_assay_catalog","animal_assay_catalog")
  identities[identified_source] <- paste0(identities[identified_source],"\036",rows$primer_id[identified_source])
  rows <- rows[!duplicated(identities), , drop = FALSE]
  if (anyDuplicated(rows$primer_id)) stop("Primer catalog IDs identify more than one oligo.")
  rows$source_key[!nzchar(rows$source_key)] <- paste0("mixer_", rows$primer_id[!nzchar(rows$source_key)])
  # Locus aliases expose the same original oligos in the matching SSU/LSU
  # selector. Their chemical sequence, citation and assay context stay intact;
  # an alias does not create another published assay or validated pairing.
  aliases <- lapply(c("18S","28S"),function(marker) {
    eligible <- if(marker=="18S") c("SSU","SSU/ITS") else c("LSU","ITS/LSU")
    own <- rows[rows$catalog_origin=="UNITE_primer_catalog" & !is.na(rows$gene_locus) & rows$gene_locus %in% eligible,,drop=FALSE]
    if(!nrow(own)) return(own)
    own$original_primer_id <- own$primer_id; own$original_marker_id <- own$marker_id
    own$original_gene_locus <- own$gene_locus
    own$marker_id <- rep(marker,nrow(own)); own$target_group <- rep("FUNGI",nrow(own))
    own$primer_id <- paste0("unite_locus_",vapply(paste(marker,own$original_primer_id,sep="\036"),.primer_mixer_hash,character(1)))
    own$catalog_origin <- rep("UNITE_locus_alias",nrow(own))
    own$evidence_note <- paste(own$evidence_note,
      "Selector alias for the original UNITE SSU/LSU component. Original chemistry, citation and assay context are unchanged; newly assembled pairings remain candidates.")
    own
  })
  rows <- .primer_mixer_bind(c(list(rows),aliases))
  if(anyDuplicated(rows$primer_id)) stop("Locus alias identities repeat.")
  rownames(rows) <- NULL
  attr(rows, "excluded_catalog_rows") <- excluded
  rows
}

# FASTA headers are display identities, not inferred taxonomy. All coordinates
# subsequently reported are on these input sequences, starting at base one.
primer_mixer_parse_fasta <- function(text, marker_id, max_templates = 100L, max_total_bases = 200000L) {
  if (length(text) != 1L || is.na(text) || !nzchar(trimws(text))) return(data.frame(template_id=character(), sequence=character(), marker_id=character(), stringsAsFactors=FALSE))
  lines <- trimws(strsplit(text, "\n", fixed = TRUE)[[1L]])
  lines <- lines[nzchar(lines)]
  starts <- which(startsWith(lines, ">"))
  if (!length(starts) || starts[1L] != 1L) stop("Paste FASTA with a >header before each sequence.")
  if (length(starts) > max_templates) stop("Use at most ", max_templates, " templates per exploration.")
  ends <- c(starts[-1L] - 1L, length(lines))
  ids <- substring(lines[starts], 2L)
  if (any(!nzchar(ids)) || anyDuplicated(ids)) stop("FASTA headers must be nonempty and unique.")
  sequences <- vapply(seq_along(starts), function(i) {
    if (ends[i] <= starts[i]) stop("Every FASTA header needs a sequence.")
    .primer_mixer_sequence(paste(lines[seq.int(starts[i]+1L, ends[i])], collapse = ""))$sequence
  }, character(1))
  if (sum(nchar(sequences)) > max_total_bases) stop("Use at most ", max_total_bases, " total template bases per exploration.")
  data.frame(template_id=ids, sequence=sequences, marker_id=marker_id, organism="User-provided template",
    source_url="", version="session_input", stringsAsFactors=FALSE)
}

.primer_mixer_empty_sites <- function() data.frame(template_id=character(), primer_id=character(), direction=character(),
  status=character(), start=integer(), end=integer(), mismatches=integer(), terminal_5_mismatches=integer(),
  terminal_base_mismatch=logical(), mismatch_positions_5_to_3=character(), n_best_sites=integer(),
  best_start=integer(), best_end=integer(), best_mismatches=integer(), ambiguous_template_bases=integer(),
  binding_sequence_primer_oriented=character(), stringsAsFactors=FALSE)

.primer_mixer_find_site <- function(sequence, primer, direction, max_mismatches) {
  out <- list(status="sequence_unavailable", start=NA_integer_, end=NA_integer_, mismatches=NA_integer_,
    terminal_5_mismatches=NA_integer_, terminal_base_mismatch=NA, mismatch_positions_5_to_3="",
    n_best_sites=0L, best_start=NA_integer_, best_end=NA_integer_, best_mismatches=NA_integer_,
    ambiguous_template_bases=NA_integer_, binding_sequence_primer_oriented="")
  if (is.na(sequence) || !nzchar(sequence)) return(out)
  sequence <- .primer_mixer_sequence(sequence)$sequence
  query <- if (direction == "reverse") .primer_mixer_complement(primer) else primer
  width <- nchar(query)
  n <- nchar(sequence) - width + 1L
  if (n < 1L) { out$status <- "template_shorter_than_oligo"; return(out) }
  codes <- list(A="A", C="C", G="G", T="T", R=c("A","G"), Y=c("C","T"), S=c("G","C"), W=c("A","T"),
    K=c("G","T"), M=c("A","C"), B=c("C","G","T"), D=c("A","G","T"), H=c("A","C","T"), V=c("A","C","G"), N=c("A","C","G","T"))
  template_bases <- strsplit(sequence, "", fixed=TRUE)[[1L]]
  query_bases <- strsplit(query, "", fixed=TRUE)[[1L]]
  starts <- seq_len(n)
  mismatches <- ambiguous <- integer(n)
  for (k in seq_len(width)) {
    observed <- template_bases[starts+k-1L]
    incompatible <- !observed %in% codes[[query_bases[k]]]
    uncertain <- !observed %in% c("A", "C", "G", "T")
    # An ambiguous reference cannot be declared an exact match even if one of
    # its possible bases would match; count only definite incompatibilities.
    if (any(uncertain)) incompatible[uncertain] <- vapply(observed[uncertain], function(x) !length(intersect(codes[[x]], codes[[query_bases[k]]])), logical(1))
    mismatches <- mismatches + as.integer(incompatible)
    ambiguous <- ambiguous + as.integer(uncertain)
  }
  complete <- which(ambiguous == 0L)
  if (!length(complete)) {
    out$status <- "ambiguous_template_sequence"
    out$ambiguous_template_bases <- min(ambiguous)
    return(out)
  }
  best <- complete[mismatches[complete] == min(mismatches[complete])]
  out$best_mismatches <- as.integer(mismatches[best[1L]])
  out$n_best_sites <- length(best)
  if (length(best)==1L) {
    out$best_start <- best[1L]
    out$best_end <- best[1L]+width-1L
  }
  if (out$best_mismatches > max_mismatches) {
    out$status <- if (any(ambiguous > 0L & mismatches <= max_mismatches)) "ambiguous_template_candidate" else "no_site_within_mismatch_limit"
    return(out)
  }
  # Uncertain windows with equal/better compatibility prevent unique placement.
  if (any(ambiguous > 0L & mismatches <= out$best_mismatches)) { out$status <- "ambiguous_template_candidate"; return(out) }
  if (length(best)>1L) { out$status <- "multiple_equally_good_sites"; return(out) }
  out$start <- out$best_start
  out$end <- out$best_end
  out$mismatches <- out$best_mismatches
  observed <- template_bases[seq.int(out$start, out$end)]
  bad <- which(!vapply(seq_len(width), function(k) observed[k] %in% codes[[query_bases[k]]], logical(1)))
  primer_bad <- sort(if (direction == "reverse") width-bad+1L else bad)
  out$mismatch_positions_5_to_3 <- paste(primer_bad, collapse=";")
  out$terminal_5_mismatches <- sum(primer_bad > width-min(5L,width))
  out$terminal_base_mismatch <- width %in% primer_bad
  out$ambiguous_template_bases <- 0L
  binding <- paste(observed, collapse="")
  out$binding_sequence_primer_oriented <- if (direction == "reverse") .primer_mixer_complement(binding) else binding
  out$status <- if (out$mismatches == 0L) "exact_site" else "mismatch_site"
  out
}

# Scan individual components for an overview without requiring any pair
# selection. The explicit cap bounds the work on an interactive Shiny worker.
primer_mixer_sites <- function(components, templates, marker_id, max_mismatches = 0L,
  max_search_bases = 2500000L) {
  if (!is.data.frame(components) || !all(c("primer_id","sequence","direction","marker_id") %in% names(components))) stop("Site scanning requires identified marker-specific oligos.")
  if (!is.data.frame(templates) || !all(c("template_id","sequence","marker_id") %in% names(templates))) stop("Site scanning requires explicit raw template identities and marker IDs.")
  if (anyNA(components$primer_id) || anyDuplicated(components$primer_id) || anyNA(templates$template_id) || anyDuplicated(templates$template_id)) stop("Component and template identities must be unique.")
  if (length(max_mismatches)!=1L || !is.finite(max_mismatches) || max_mismatches!=floor(max_mismatches) || max_mismatches<0L || max_mismatches>5L) stop("Choose a maximum of 0 to 5 substitutions per oligo.")
  if (nrow(components)>500L || nrow(templates)>100L || nrow(components)*sum(nchar(templates$sequence),na.rm=TRUE)>max_search_bases) stop("Narrow the overview to fewer primers or templates before scanning.")
  if (!nrow(components) || !nrow(templates)) return(.primer_mixer_empty_sites())
  if (any(!components$direction %in% c("forward","reverse"))) stop("Assign forward or reverse directions before scanning.")
  .primer_mixer_bind(lapply(seq_len(nrow(templates)),function(ti) {
    .primer_mixer_bind(lapply(seq_len(nrow(components)),function(ci) {
      matching_marker <- !is.na(templates$marker_id[ti]) && !is.na(components$marker_id[ci]) &&
        templates$marker_id[ti]==marker_id && components$marker_id[ci]==marker_id
      value <- .primer_mixer_find_site(if (matching_marker) templates$sequence[ti] else NA_character_,
        .primer_mixer_sequence(components$sequence[ci])$sequence,components$direction[ci],max_mismatches)
      if (!matching_marker) value$status <- "marker_mismatch"
      data.frame(template_id=templates$template_id[ti],primer_id=components$primer_id[ci],direction=components$direction[ci],
        as.data.frame(value,stringsAsFactors=FALSE),stringsAsFactors=FALSE)
    }))
  }))
}

# The compact starter panel is a coordinate/sequence example, never a claim
# of broad taxonomic coverage. Fish templates cannot silently become an
# arthropod reference or the converse.
.primer_mixer_read_starter_csv <- function(path) {
  manifest_path <- file.path(dirname(path),"manifest.json")
  name <- basename(path)
  fail <- function(reason) stop("Starter reference integrity failed for ",path,": ",reason,". No substitute reference was loaded.",call.=FALSE)
  if(!file.exists(manifest_path)) fail("missing frozen manifest")
  if(!file.exists(path)) fail("missing required artifact")
  manifest <- tryCatch(jsonlite::read_json(manifest_path,simplifyVector=FALSE),error=function(e) fail("unreadable frozen manifest"))
  if(!is.list(manifest) || !is.list(manifest$artifacts)) fail("invalid frozen artifact manifest")
  artifact <- manifest$artifacts[[name]]
  if(!is.list(artifact) || !is.character(artifact$sha256) || length(artifact$sha256)!=1L || is.na(artifact$sha256) || !grepl("^[a-fA-F0-9]{64}$",artifact$sha256)) fail("missing or invalid SHA-256 descriptor")
  observed <- digest::digest(file=path,algo="sha256",serialize=FALSE)
  if(!identical(tolower(observed),tolower(artifact$sha256))) fail("artifact SHA-256 mismatch")
  if(!is.null(artifact$bytes) && !isTRUE(file.info(path)$size==as.numeric(artifact$bytes))) fail("artifact byte count mismatch")
  rows <- tryCatch(utils::read.csv(path,stringsAsFactors=FALSE,check.names=FALSE),error=function(e) fail("unreadable CSV artifact"))
  if(name=="templates.csv" && !is.null(manifest$n_templates) && nrow(rows)!=manifest$n_templates) fail("template count differs from frozen manifest")
  if(name=="taxonomy.csv") {
    supplement <- manifest$taxonomy_supplement
    if(!is.list(supplement) || length(supplement$version)!=1L || !nzchar(supplement$version)) fail("missing taxonomy supplement identity")
    if(is.null(supplement$n_rows) || nrow(rows)!=supplement$n_rows) fail("taxonomy count differs from frozen supplement")
    if(!"taxonomy_version" %in% names(rows) || anyNA(rows$taxonomy_version) || any(rows$taxonomy_version!=supplement$version)) fail("taxonomy version differs from frozen supplement")
    if(!is.list(supplement$sources) || !length(supplement$sources)) fail("missing taxonomy source descriptors")
    source_hashes <- vapply(supplement$sources,function(source) if(!is.list(source) || is.null(source$uncompressed_sha256)) "" else source$uncompressed_sha256,character(1))
    if(!"taxonomy_source_sha256" %in% names(rows) || anyNA(rows$taxonomy_source_sha256) ||
        any(!rows$taxonomy_source_sha256 %in% source_hashes)) fail("taxonomy source hash is outside the frozen supplement")
  }
  rows
}

.primer_mixer_template_taxonomy <- function(rows,root) {
  ranks <- c("superkingdom","kingdom","phylum","class","order","family","genus","species")
  for(rank in ranks) if(!rank %in% names(rows)) rows[[rank]] <- rep("Unresolved",nrow(rows))
  defaults <- c(taxid="",taxonomy_status="No source-backed taxonomy supplement for this template",
    taxonomy_source_url="",taxonomy_version="",taxonomy_source_sha256="")
  for(field in names(defaults)) if(!field %in% names(rows)) rows[[field]] <- rep(defaults[[field]],nrow(rows))
  paths <- file.path(root,c("data/pinned/FISH_STARTER/v1/taxonomy.csv","data/pinned/ANIMAL_STARTER/v1/taxonomy.csv"))
  declared <- vapply(paths,function(path) dir.exists(dirname(path)),logical(1))
  tables <- lapply(paths[declared],.primer_mixer_read_starter_csv)
  taxonomy <- .primer_mixer_bind(tables)
  if(!nrow(taxonomy) || !nrow(rows)) return(rows)
  if(!"template_id" %in% names(taxonomy) || anyNA(taxonomy$template_id) || any(!nzchar(taxonomy$template_id)) || anyDuplicated(taxonomy$template_id)) stop("Taxonomy supplements require unique, nonmissing raw template identities.")
  index <- match(rows$template_id,taxonomy$template_id)
  matched <- which(!is.na(index))
  if(!length(matched)) return(rows)
  # A template ID cannot license a taxonomy row for a different marker, group
  # or accession version. Validate redundant source identity before joining.
  for(field in intersect(c("marker_id","target_group","group_id","accession_version"),intersect(names(rows),names(taxonomy)))) {
    original <- as.character(rows[[field]][matched]); incoming <- as.character(taxonomy[[field]][index[matched]])
    known <- !is.na(original) & nzchar(original) & !is.na(incoming) & nzchar(incoming)
    if(any(original[known]!=incoming[known])) stop("Taxonomy supplement disagrees with raw template identity: ",field,".")
  }
  fields <- setdiff(names(taxonomy),c("template_id","marker_id","group_id","target_group","accession_version"))
  for(field in fields) {
    # Sequence source fields keep their original meaning. Taxonomy's own raw
    # archives and names are retained separately, never overwriting them.
    destination <- if(field %in% c(ranks,"taxid","taxonomy_status") || startsWith(field,"taxonomy_")) field else paste0("taxonomy_",field)
    if(!destination %in% names(rows)) rows[[destination]] <- rep("",nrow(rows))
    value <- taxonomy[[field]][index[matched]]
    if(field %in% ranks) value[is.na(value) | !nzchar(as.character(value))] <- "Unresolved"
    rows[[destination]][matched] <- value
  }
  rows
}

# Extract only individually bounded SSU/LSU source intervals. Disjoint partial
# landmarks remain separate molecules; unknown gene boundaries never become
# whole-rDNA templates relabeled as 18S or 28S.
primer_mixer_fungal_templates <- function(snapshot,marker_id) {
  if(length(marker_id)!=1L || !marker_id %in% c("18S","28S")) stop("Choose 18S or 28S for a bounded fungal extraction.")
  if(!is.list(snapshot) || !all(c("manifest","sequences","landmarks","taxonomy") %in% names(snapshot))) stop("Fungal extraction requires one complete immutable snapshot.")
  empty <- data.frame(template_id=character(),sequence=character(),marker_id=character(),organism=character(),source_url=character(),
    version=character(),target_group=character(),accession_version=character(),source_start=integer(),source_end=integer(),
    source_status=character(),partial=logical(),projected=logical(),boundary_status=character(),evidence_status=character(),stringsAsFactors=FALSE)
  attr(empty,"unavailable_reason") <- "This frozen version has no independently bounded source gene interval for the selected marker. Curated pilot-display anchors are withheld from gene extraction."
  genes <- if(marker_id=="18S") c("18S","SSU") else c("28S","LSU")
  landmarks <- as.data.frame(snapshot$landmarks,stringsAsFactors=FALSE)
  if(!nrow(landmarks)) return(empty)
  if(!all(c("landmark_id","label","start","end") %in% names(landmarks))) stop("Fungal extraction landmarks are incomplete.")
  if(!"accession" %in% names(landmarks)) landmarks$accession <- rep(snapshot$manifest$reference_accession,nrow(landmarks))
  for(field in c("boundary_status","evidence_status","source_note","mapping_method","landmark_type")) {
    if(!field %in% names(landmarks)) landmarks[[field]] <- rep("",nrow(landmarks))
    landmarks[[field]] <- .primer_mixer_text(landmarks[[field]])
  }
  landmarks <- landmarks[!is.na(landmarks$label) & landmarks$label %in% genes &
    !landmarks$boundary_status %in% c("unresolved","not_detected","invalid_coordinates","withheld","extraction_mismatch") &
    !landmarks$landmark_type %in% c("domain","divergent_domain","conserved_anchor") & landmarks$evidence_status!="mapped_reference_domain",,drop=FALSE]
  if(!nrow(landmarks)) return(empty)
  records <- lapply(seq_len(nrow(landmarks)),function(i) {
    feature <- landmarks[i,,drop=FALSE]; accession <- as.character(feature$accession)
    if(length(accession)!=1L || is.na(accession) || !accession %in% names(snapshot$sequences)) return(NULL)
    source <- unname(snapshot$sequences[[accession]])
    if(length(source)!=1L || is.na(source) || !nzchar(source)) return(NULL)
    bounds <- suppressWarnings(as.numeric(c(feature$start,feature$end)))
    if(any(!is.finite(bounds)) || any(bounds!=floor(bounds)) || bounds[1L]<1L || bounds[2L]<bounds[1L] || bounds[2L]>nchar(source)) return(NULL)
    status <- feature$boundary_status
    if(!nzchar(status)) status <- if(grepl("partial",feature$source_note,ignore.case=TRUE)) "partial" else
      if(grepl("complete",feature$source_note,ignore.case=TRUE)) "complete" else "source_interval"
    projected <- grepl("project|map",paste(feature$evidence_status,feature$mapping_method,status),ignore.case=TRUE)
    partial <- status!="complete" || grepl("partial",feature$source_note,ignore.case=TRUE)
    tax <- snapshot$taxonomy[snapshot$taxonomy$accession==accession,,drop=FALSE]
    if(nrow(tax)>1L) stop("Frozen fungal taxonomy repeats a source accession.")
    field <- function(name,default="") if(nrow(tax) && name %in% names(tax)) .primer_mixer_text(tax[[name]][1L],default) else default
    sequence <- substr(source,bounds[1L],bounds[2L])
    row <- data.frame(template_id=paste0("FUNGI_",accession,"_",marker_id,"_",feature$landmark_id),sequence=sequence,marker_id=marker_id,
      organism=field("organism","Unresolved fungus"),source_url=paste0("https://www.ncbi.nlm.nih.gov/nuccore/",accession),
      version=as.character(snapshot$manifest$id),target_group="FUNGI",accession_version=accession,
      source_start=as.integer(bounds[1L]),source_end=as.integer(bounds[2L]),source_feature_location=paste0(bounds[1L],"..",bounds[2L]),
      source_status=paste("Frozen NCBI fungal",if(projected) "projected" else "source",marker_id,"interval;",status),
      partial=partial,projected=projected,boundary_status=status,evidence_status=feature$evidence_status,
      source_landmark_id=feature$landmark_id,source_landmark_type=feature$landmark_type,
      source_note=feature$source_note,mapping_method=feature$mapping_method,
      parent_sequence_sha256=digest::digest(source,algo="sha256",serialize=FALSE),
      sequence_sha256=digest::digest(sequence,algo="sha256",serialize=FALSE),
      taxid=field("taxid"),taxonomy_status="frozen_NCBI_taxonomy",taxonomy_version=as.character(snapshot$manifest$id),
      taxonomy_source_url=if(nzchar(field("taxid"))) paste0("https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=",field("taxid")) else "",
      taxonomy_source_sha256=field("taxonomy_source_sha256"),stringsAsFactors=FALSE)
    for(rank in c("superkingdom","kingdom","phylum","class","order","family","genus","species")) {
      value <- field(rank,"Unresolved")
      if(!nzchar(value) || tolower(value) %in% c("unresolved","unknown","unclassified")) value <- "Unresolved"
      row[[rank]] <- value
    }
    row
  })
  rows <- .primer_mixer_bind(records)
  if(!nrow(rows)) return(empty)
  if(anyDuplicated(rows$template_id)) stop("Fungal landmark extraction identities repeat.")
  rows
}

primer_mixer_templates <- function(root = ".", marker_id = NULL, target_group = NULL) {
  result <- data.frame(template_id=character(),sequence=character(),marker_id=character(),organism=character(),
    source_url=character(),version=character(),target_group=character(),stringsAsFactors=FALSE)
  empty_contract <- result
  group <- if(is.null(target_group)) NULL else toupper(as.character(target_group))
  if(length(group)>1L || (!is.null(group) && is.na(group))) stop("Choose one target group when loading starter templates.")
  all_groups <- is.null(group) || group %in% c("ALL","METAZOA")
  fish_path <- file.path(root,"data/pinned/FISH_STARTER/v1/templates.csv")
  fish_declared <- dir.exists(dirname(fish_path))
  if (fish_declared && (all_groups || group=="FISH")) {
    fish <- .primer_mixer_read_starter_csv(fish_path)
    if (!all(c("template_id","sequence","marker_id") %in% names(fish))) stop("Fish starter template contract is incomplete.")
    fish$target_group <- rep("FISH",nrow(fish))
    result <- .primer_mixer_bind(list(result,fish))
  }
  animal_path <- file.path(root,"data/pinned/ANIMAL_STARTER/v1/templates.csv")
  animal_declared <- dir.exists(dirname(animal_path))
  if(animal_declared) {
    animal <- .primer_mixer_read_starter_csv(animal_path)
    if(!all(c("template_id","sequence","marker_id") %in% names(animal)) ||
       !any(c("target_group","group_id") %in% names(animal))) stop("Animal starter template contract is incomplete.")
    if(!"target_group" %in% names(animal)) animal$target_group <- animal$group_id
    if(!all_groups) animal <- animal[toupper(animal$target_group) %in% group,,drop=FALSE]
    result <- .primer_mixer_bind(list(result,animal))
  }
  coi_path <- file.path(root,"data/reference/coi_reference_NC_001322.1.fasta")
  if (file.exists(coi_path) && (all_groups || group %in% c("ARTHROPODS","ARTHROPODA"))) {
    lines <- readLines(coi_path,warn=FALSE)
    seq <- .primer_mixer_sequence(paste(lines[!startsWith(lines,">")],collapse=""))$sequence
    result <- .primer_mixer_bind(list(result,data.frame(template_id="NC_001322.1_COX1",sequence=seq,marker_id="COI",organism="Drosophila yakuba",
      source_url="https://www.ncbi.nlm.nih.gov/nuccore/NC_001322.1",version="NC_001322.1",target_group="ARTHROPODS",stringsAsFactors=FALSE)))
  }
  # A target changes before its asynchronous marker selector acknowledges the
  # new choice. No cohort is valid during e.g. FUNGI + COI. Frame binding may
  # drop all-empty inputs, so restore the typed zero-row template contract
  # before filtering or taxonomy enrichment; NULL IDs must never reach Shiny.
  if(!nrow(result)) for(field in names(empty_contract))
    if(!field %in% names(result)) result[[field]] <- empty_contract[[field]]
  if (!is.null(marker_id)) result <- result[result$marker_id %in% marker_id,,drop=FALSE]
  if(anyDuplicated(result$template_id)) stop("Starter template identities repeat across source panels.")
  if("partial" %in% names(result)) {
    original <- tolower(trimws(as.character(result$partial)))
    result$partial <- ifelse(is.na(original),NA,original %in% c("true","1","yes"))
  }
  result <- .primer_mixer_template_taxonomy(result,root)
  rownames(result) <- NULL
  result
}

primer_mixer_build <- function(catalog, forward_ids, reverse_ids, templates, marker_id,
  max_mismatches = 0L, min_bp = 0L, max_bp = Inf, max_templates = 100L, max_total_bases = 200000L) {
  required <- c("primer_id", "primer_name", "sequence", "direction", "marker_id")
  if (!is.data.frame(catalog) || !all(required %in% names(catalog))) stop("Mixer catalog requires ", paste(required, collapse=", "), ".")
  if (anyNA(catalog$primer_id) || anyDuplicated(catalog$primer_id)) stop("Primer IDs must be unique and nonmissing.")
  if (length(marker_id)!=1L || is.na(marker_id) || !nzchar(marker_id)) stop("Choose one marker.")
  if (length(max_mismatches)!=1L || !is.finite(max_mismatches) || max_mismatches!=floor(max_mismatches) || max_mismatches<0L || max_mismatches>5L) stop("Choose a maximum of 0 to 5 substitutions per oligo.")
  if (length(min_bp)!=1L || !is.finite(min_bp) || min_bp<0 || length(max_bp)!=1L || is.na(max_bp) || max_bp<min_bp) stop("Use a nonnegative increasing product-length range.")
  forward_ids <- unique(as.character(forward_ids)); reverse_ids <- unique(as.character(reverse_ids))
  if (length(forward_ids)>5L || length(reverse_ids)>5L) stop("Choose at most five forward and five reverse primers.")
  fi <- match(forward_ids, catalog$primer_id); ri <- match(reverse_ids, catalog$primer_id)
  if (anyNA(c(fi,ri))) stop("Selected primer identity is absent from the catalog.")
  if (any(catalog$direction[fi]!="forward") || any(catalog$direction[ri]!="reverse") ||
      any(catalog$marker_id[c(fi,ri)]!=marker_id)) stop("Every selected oligo must have the chosen marker and assigned direction.")
  components <- catalog[c(fi,ri),,drop=FALSE]
  for (field in c("source_url","citation","source_key","catalog_origin","normalization_note")) if (!field %in% names(components)) components[[field]] <- rep("",nrow(components))
  components$sequence <- vapply(components$sequence, function(x) .primer_mixer_sequence(x)$sequence, character(1))
  if (any(nchar(components$sequence)<8L | nchar(components$sequence)>120L)) stop("Explored oligos must contain 8 to 120 bases.")
  if (is.null(templates)) templates <- data.frame(template_id=character(),sequence=character(),marker_id=character())
  if (!is.data.frame(templates) || !all(c("template_id","sequence","marker_id") %in% names(templates))) stop("Templates require template_id, sequence and marker_id.")
  if (anyNA(templates$template_id) || any(!nzchar(templates$template_id)) || anyDuplicated(templates$template_id)) stop("Template identities must be unique and nonmissing.")
  if (nrow(templates)>max_templates) stop("Use at most ",max_templates," templates per exploration.")
  if (sum(nchar(templates$sequence),na.rm=TRUE)>max_total_bases) stop("Use at most ",max_total_bases," total template bases per exploration.")
  template_normalized <- lapply(templates$sequence,function(x) if (is.na(x)||!nzchar(x)) NULL else .primer_mixer_sequence(x))
  if (!"source_sequence" %in% names(templates)) templates$source_sequence <- as.character(templates$sequence)
  templates$normalization_note <- vapply(template_normalized,function(x) if (is.null(x)) "" else x$normalization_note,character(1))
  templates$sequence <- vapply(template_normalized,function(x) if (is.null(x)) NA_character_ else x$sequence,character(1))
  grid <- expand.grid(forward_id=forward_ids,reverse_id=reverse_ids,stringsAsFactors=FALSE)
  pairs <- data.frame(pair_id=character(),pair_label=character(),forward_id=character(),reverse_id=character(),
    forward_primer=character(),reverse_primer=character(),forward_sequence=character(),reverse_sequence=character(),
    marker_id=character(),evidence_status=character(),stringsAsFactors=FALSE)
  if (nrow(grid)) pairs <- do.call(rbind,lapply(seq_len(nrow(grid)),function(i) {
    f <- components[match(grid$forward_id[i],components$primer_id),,drop=FALSE]
    r <- components[match(grid$reverse_id[i],components$primer_id),,drop=FALSE]
    data.frame(pair_id=paste0("mix_",.primer_mixer_hash(paste(marker_id,f$primer_id,r$primer_id,sep="\036"))),
      pair_label=paste(f$primer_name,r$primer_name,sep=" + "),forward_id=f$primer_id,reverse_id=r$primer_id,
      forward_primer=f$primer_name,reverse_primer=r$primer_name,forward_sequence=f$sequence,reverse_sequence=r$sequence,
      marker_id=marker_id,evidence_status="newly_assembled_candidate",stringsAsFactors=FALSE)
  }))
  sites <- .primer_mixer_empty_sites()
  if (nrow(components) && nrow(templates)) sites <- .primer_mixer_bind(lapply(seq_len(nrow(templates)),function(ti) {
    .primer_mixer_bind(lapply(seq_len(nrow(components)),function(ci) {
      value <- .primer_mixer_find_site(if (identical(as.character(templates$marker_id[ti]),marker_id)) templates$sequence[ti] else NA_character_,
        components$sequence[ci],components$direction[ci],max_mismatches)
      if (is.na(templates$marker_id[ti]) || templates$marker_id[ti]!=marker_id) value$status <- "marker_mismatch"
      data.frame(template_id=templates$template_id[ti],primer_id=components$primer_id[ci],direction=components$direction[ci],
        as.data.frame(value,stringsAsFactors=FALSE),stringsAsFactors=FALSE)
    }))
  }))
  products <- data.frame(pair_id=character(),template_id=character(),marker_id=character(),status=character(),
    forward_status=character(),reverse_status=character(),forward_start=integer(),forward_end=integer(),
    reverse_start=integer(),reverse_end=integer(),amplicon_bp=integer(),internal_bp=integer(),
    forward_mismatches=integer(),reverse_mismatches=integer(),forward_terminal_5_mismatches=integer(),
    reverse_terminal_5_mismatches=integer(),passes_length_filter=logical(),stringsAsFactors=FALSE)
  if (nrow(pairs) && nrow(templates)) products <- .primer_mixer_bind(lapply(seq_len(nrow(pairs)),function(pi) {
    .primer_mixer_bind(lapply(seq_len(nrow(templates)),function(ti) {
      f <- sites[sites$template_id==templates$template_id[ti] & sites$primer_id==pairs$forward_id[pi],,drop=FALSE]
      r <- sites[sites$template_id==templates$template_id[ti] & sites$primer_id==pairs$reverse_id[pi],,drop=FALSE]
      f_ok <- f$status %in% c("exact_site","mismatch_site"); r_ok <- r$status %in% c("exact_site","mismatch_site")
      state <- if (f$status=="marker_mismatch" || r$status=="marker_mismatch") "marker_mismatch" else
        if (f$status=="sequence_unavailable" || r$status=="sequence_unavailable") "template_unavailable" else
        if (any(c(f$status,r$status) %in% c("multiple_equally_good_sites","ambiguous_template_sequence","ambiguous_template_candidate"))) "ambiguous_sites" else
        if (f_ok && r_ok) { if (f$end>=r$start) "non_inward_or_overlapping_sites" else if (f$mismatches+r$mismatches==0L) "exact_product" else "mismatch_product" } else
        if (f_ok || r_ok) "one_eligible_site" else if (any(c(f$status,r$status)=="template_shorter_than_oligo")) "template_too_short" else "no_eligible_sites"
      mapped <- state %in% c("exact_product","mismatch_product")
      span <- if (mapped) r$end-f$start+1L else NA_integer_
      data.frame(pair_id=pairs$pair_id[pi],template_id=templates$template_id[ti],marker_id=marker_id,status=state,
        forward_status=f$status,reverse_status=r$status,forward_start=f$start,forward_end=f$end,reverse_start=r$start,reverse_end=r$end,
        amplicon_bp=span,internal_bp=if (mapped) r$start-f$end-1L else NA_integer_,
        forward_mismatches=f$mismatches,reverse_mismatches=r$mismatches,forward_terminal_5_mismatches=f$terminal_5_mismatches,
        reverse_terminal_5_mismatches=r$terminal_5_mismatches,passes_length_filter=mapped && span>=min_bp && span<=max_bp,stringsAsFactors=FALSE)
    }))
  }))
  summary <- data.frame(pair_id=character(),pair_label=character(),n_templates=integer(),n_matching_marker_templates=integer(),
    n_exact_products=integer(),n_mismatch_products=integer(),n_length_eligible=integer(),n_ambiguous=integer(),
    n_partial=integer(),n_no_eligible_sites=integer(),n_unavailable=integer(),stringsAsFactors=FALSE)
  if (nrow(pairs)) summary <- do.call(rbind,lapply(seq_len(nrow(pairs)),function(pi) {
    rows <- products[products$pair_id==pairs$pair_id[pi],,drop=FALSE]
    data.frame(pair_id=pairs$pair_id[pi],pair_label=pairs$pair_label[pi],n_templates=nrow(templates),
      n_matching_marker_templates=sum(templates$marker_id %in% marker_id),n_exact_products=sum(rows$status=="exact_product"),
      n_mismatch_products=sum(rows$status=="mismatch_product"),n_length_eligible=sum(rows$passes_length_filter),n_ambiguous=sum(rows$status=="ambiguous_sites"),
      n_partial=sum(rows$status %in% c("one_eligible_site","template_too_short")),n_no_eligible_sites=sum(rows$status %in% c("no_eligible_sites","non_inward_or_overlapping_sites")),
      n_unavailable=sum(rows$status %in% c("marker_mismatch","template_unavailable")),stringsAsFactors=FALSE)
  }))
  list(components=components,pairs=pairs,sites=sites,products=products,summary=summary,templates=templates,marker_id=marker_id,
    settings=list(max_mismatches=as.integer(max_mismatches),min_bp=min_bp,max_bp=max_bp),
    method="Unique best full-length ungapped oligo site per direction, IUPAC-compatible concrete reference bases, substitution counts and terminal-five mismatch counts; no indels. Product coordinates are inclusive raw input bases. No PCR probability is inferred.")
}

# Convert a selected template's measured products to the established COI-style
# map. Unavailable or ambiguous products stay in result$products, never drawn
# as a complete lane. The adapter does not pool coordinates between templates.
primer_mixer_platform_tier <- function(bp) {
  if(!is.numeric(bp) || any(!is.na(bp) & (!is.finite(bp) | bp<1L | bp!=floor(bp)))) stop("Planning tiers require positive inclusive integer product lengths.")
  ifelse(is.na(bp),NA_character_,ifelse(bp<=500L,"ILLUMINA",ifelse(bp<=2000L,"PACBIO",ifelse(bp<=3500L,"NANOPORE","BEYOND"))))
}

primer_mixer_overview_data <- function(result, template_id) {
  ti <- match(template_id,result$templates$template_id)
  if (length(ti)!=1L || is.na(ti)) stop("Choose one template from this mixer result.")
  template <- result$templates[ti,,drop=FALSE]
  sequence <- template$sequence[1L]
  length_bp <- if (is.na(sequence)) 0L else nchar(sequence)
  rows <- result$products[result$products$template_id==template_id & result$products$passes_length_filter,,drop=FALSE]
  rows$pair_label <- result$pairs$pair_label[match(rows$pair_id,result$pairs$pair_id)]
  rows$pair_start <- rows$forward_start; rows$pair_end <- rows$reverse_end
  rows$aligned_amplicon_bp <- rows$amplicon_bp; rows$targeted_region_bp <- rows$internal_bp
  rows$placement_status <- ifelse(rows$status=="exact_product","exact input-template sequence match","Input-template substitution match; review mismatch counts")
  rows$target <- rep(if("organism" %in% names(template)) .primer_mixer_text(template$organism[1L]) else template_id,nrow(rows))
  rows$use_case <- rep("Session mix and match candidate",nrow(rows))
  rows$folmer_overlap_bp <- rep(NA_integer_,nrow(rows))
  rows$length_basis <- rep("Inclusive ungapped raw template bases; internal bp excludes both oligos",nrow(rows))
  rows$platform_tier <- primer_mixer_platform_tier(rows$amplicon_bp)
  primers <- .primer_mixer_bind(lapply(seq_len(nrow(result$pairs)),function(i) {
    pair <- result$pairs[i,,drop=FALSE]
    cs <- result$components[match(c(pair$forward_id,pair$reverse_id),result$components$primer_id),,drop=FALSE]
    cs$pair_id <- rep(pair$pair_id,2L); cs$pair_label <- rep(pair$pair_label,2L)
    cs$source_key <- ifelse(cs$catalog_origin=="custom_session_input","custom",paste0("mixer_",cs$primer_id))
    cs
  }))
  citations <- data.frame(key=paste0("mixer_",result$components$primer_id),short_citation=.primer_mixer_text(result$components$citation),
    title=.primer_mixer_text(result$components$citation),url=.primer_mixer_text(result$components$source_url),stringsAsFactors=FALSE)
  list(rows=rows,primers=primers,citations=citations,marker=data.frame(marker_id=result$marker_id,
    display_name=result$marker_id,reference_length=length_bp,coordinate_reference=template_id,
    reference_version=if("version" %in% names(template)) .primer_mixer_text(template$version[1L]) else "session_input",stringsAsFactors=FALSE),
    reference_length=length_bp,reference_sequence=sequence,template_id=template_id,
    focus=c(start=1,end=max(2L,length_bp)),domain=c(min=1,max=max(2L,length_bp)))
}
