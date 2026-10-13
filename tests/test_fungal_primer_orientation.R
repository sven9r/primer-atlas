#!/usr/bin/env Rscript
.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
source("R/functions.R")
source("R/data_layer.R")
source("R/fungal_explorer.R")
source("R/fungal_versions.R")
source("R/fungal_primer_orientation.R")

snapshot <- fungal_read_snapshot("v3-c6fa87f256e1", remote=FALSE)
before <- serialize(snapshot, NULL)
original <- as.data.frame(snapshot$catalog, stringsAsFactors=FALSE)
variants <- fungal_orientation_catalog(snapshot)
stopifnot(nrow(variants)==112L, !anyDuplicated(variants$catalog_id),
  !any(variants$catalog_id %in% original$catalog_id),
  identical(variants$parent_catalog_id, original$catalog_id),
  identical(variants$parent_primer_name, original$primer_name),
  identical(variants$parent_sequence, original$sequence),
  identical(variants$parent_direction, original$direction),
  all(variants$direction != original$direction),
  all(variants$primer_name == paste0(original$primer_name, " [RC]")),
  all(grepl("^FUNGAL_RC_[a-f0-9]{24}$", variants$catalog_id)),
  all(variants$primary_reference_status == "derived_reverse_complement"),
  all(variants$orientation_derived), all(variants$raw_query_equivalence),
  all(variants$source_version == snapshot$manifest$id),
  all(variants$parent_source_version == snapshot$manifest$id),
  all(is.na(variants$reported_position)),
  identical(variants$parent_reference, .fungal_orientation_field(original,"reference")),
  identical(variants$parent_primary_reference_url, .fungal_orientation_field(original,"primary_reference_url")),
  all(grepl("not a published assay", variants$reference, fixed=TRUE)),
  identical(variants$sequence, vapply(original$sequence, reverse_complement, character(1), USE.NAMES=FALSE)),
  isTRUE(fungal_orientation_validate(snapshot, variants)),
  identical(serialize(snapshot, NULL), before))

# Resolve in either direction through original or canonical derived IDs. Going
# back to the original returns its exact source row and citation status.
for (i in seq_len(nrow(original))) {
  o <- original[i,,drop=FALSE]
  v <- variants[i,,drop=FALSE]
  stopifnot(identical(fungal_orientation_resolve(snapshot,o$catalog_id,o$direction),o),
    identical(fungal_orientation_resolve(snapshot,o$catalog_id,v$direction),v),
    identical(fungal_orientation_resolve(snapshot,v$catalog_id,o$direction),o),
    identical(fungal_orientation_resolve(snapshot,v$catalog_id,v$direction),v))
}

# The opposite oligo has the SAME raw compatible nucleotide-set query, not the
# same 5'-to-3' synthesis sequence. I/U normalization follows existing search
# semantics and must be disclosed rather than claimed as chemical equivalence.
raw_query <- function(sequence, direction) {
  query <- if (direction=="reverse") reverse_complement(sequence) else sequence
  chartr("IU", "NT", query)
}
for (i in seq_len(nrow(original))) {
  stopifnot(identical(raw_query(original$sequence[i],original$direction[i]),
    raw_query(variants$sequence[i],variants$direction[i])))
}
inosine <- grepl("I", original$sequence, fixed=TRUE)
stopifnot(sum(inosine)==2L, all(variants$parent_noncanonical_symbols[inosine]=="I"),
  !any(grepl("I",variants$sequence[inosine],fixed=TRUE)),
  all(grepl("not chemically identical to inosine",variants$sequence_derivation_note[inosine],fixed=TRUE)),
  all(grepl("not demonstrated biochemical equivalence",variants$sequence_derivation_note[inosine],fixed=TRUE)))

# Compare actual exact-site searches for every component on two real templates.
# Unique, missing and repeated-compatible sites must not change when extension
# direction changes; sites alone do not imply inward-facing product geometry.
for (accession in c("PZ962702.1", "FN812768.2")) {
  bases <- strsplit(snapshot$sequences[[accession]],"",fixed=TRUE)[[1L]]
  for (i in seq_len(nrow(original))) {
    left <- .fungal_exact_site(bases,original$sequence[i],original$direction[i])
    right <- .fungal_exact_site(bases,variants$sequence[i],variants$direction[i])
    stopifnot(identical(left,right))
  }
}
sequence <- "TCCGTAGGTGAACCTGCGG"
reverse <- .fungal_orientation_rc(sequence)
duplicate_bases <- strsplit(paste0(sequence,strrep("C",40),sequence),"",fixed=TRUE)[[1L]]
ambiguous <- .fungal_exact_site(duplicate_bases,sequence,"forward")
stopifnot(ambiguous$status=="ambiguous_binding_site", ambiguous$n_matches==2L,
  identical(ambiguous,.fungal_exact_site(duplicate_bases,reverse,"reverse")))

# ITS2's EFFECTIVE catalog role is reverse, despite the imported direction.
# Its derived forward sequence equals ITS3, but its distinct identity/name
# cannot inherit the documented ITS3+ITS4 assay or its component citation flag.
its2 <- which(original$primer_name=="ITS2")
its3 <- which(original$primer_name=="ITS3")
its4 <- which(original$primer_name=="ITS4")
stopifnot(length(its2)==1L,original$direction[its2]=="reverse",original$source_direction[its2]=="forward",
  variants$direction[its2]=="forward",variants$sequence[its2]==original$sequence[its3],
  variants$catalog_id[its2]!=original$catalog_id[its3],variants$primer_name[its2]!=original$primer_name[its3])
derived_catalog <- original[c(its2,its4),,drop=FALSE]
for (field in names(original)) derived_catalog[1L,field] <- variants[[field]][its2]
evaluation <- fungal_evaluate_panel(derived_catalog,snapshot$sequences["PZ962702.1"])
stopifnot(nrow(evaluation$pairs)==1L, evaluation$pairs$evidence_status=="candidate_incomplete_provenance",
  evaluation$pair_sequences$status=="scorable", evaluation$pair_sequences$amplicon_bp==391L,
  evaluation$pair_sequences$forward_start==2090L, evaluation$pair_sequences$reverse_end==2480L)

# Identity is stable across historical versions only when the parent component
# tuple is unchanged; source provenance continues to identify the chosen version.
historical <- fungal_read_snapshot("v1-9e958f105c46",remote=FALSE)
old_variants <- fungal_orientation_catalog(historical)
match_ids <- match(variants$parent_catalog_id,old_variants$parent_catalog_id)
stopifnot(!anyNA(match_ids),identical(variants$catalog_id,old_variants$catalog_id[match_ids]),
  all(old_variants$parent_source_version==historical$manifest$id))
other_collection <- snapshot
other_collection$manifest$id <- "same-component-other-collection"
stopifnot(inherits(try(fungal_orientation_validate(other_collection,variants),silent=TRUE),"try-error"),
  isTRUE(fungal_orientation_validate(other_collection,variants,require_source_version=FALSE)))
working <- snapshot
working$frozen_component_catalog <- snapshot$catalog
working$catalog <- variants
stopifnot(identical(fungal_orientation_catalog(working),variants))

# Fail closed on a role-only flip, forged parent, altered canonical ID, false
# proof, reused publication status, or inconsistent originating version.
for (field in c("catalog_id", "primer_name", "sequence", "direction", "parent_catalog_id",
    "parent_primer_name", "parent_sequence", "parent_direction", "primary_reference_status",
    "parent_reference", "orientation_origin", "derivation_rule", "sequence_derivation_note")) {
  forged <- variants[1L,,drop=FALSE]
  forged[[field]] <- if (field=="sequence") original$sequence[1L] else "forged"
  stopifnot(inherits(try(fungal_orientation_validate(snapshot,forged),silent=TRUE),"try-error"))
}
for (field in c("orientation_derived", "raw_query_equivalence")) {
  forged <- variants[1L,,drop=FALSE];forged[[field]] <- FALSE
  stopifnot(inherits(try(fungal_orientation_validate(snapshot,forged),silent=TRUE),"try-error"))
}
forged <- variants[1L,,drop=FALSE];forged$source_version <- "different-version"
stopifnot(inherits(try(fungal_orientation_validate(snapshot,forged,require_source_version=FALSE),silent=TRUE),"try-error"),
  inherits(try(fungal_orientation_resolve(snapshot,"unknown","forward"),silent=TRUE),"try-error"),
  inherits(try(fungal_orientation_resolve(snapshot,original$catalog_id[1L],"sideways"),silent=TRUE),"try-error"),
  inherits(try(fungal_orientation_validate(snapshot,rbind(variants[1L,],variants[1L,])),silent=TRUE),"try-error"))
for (sequence in c("AC-GT", "acgt", " ACGT", "ACGTZ", NA_character_)) {
  malformed <- snapshot;malformed$catalog$sequence[1L] <- sequence
  stopifnot(inherits(try(fungal_orientation_catalog(malformed),silent=TRUE),"try-error"))
}
malformed <- snapshot;malformed$catalog <- variants
stopifnot(inherits(try(fungal_orientation_catalog(malformed),silent=TRUE),"try-error"))
empty <- snapshot;empty$catalog <- original[FALSE,,drop=FALSE]
stopifnot(nrow(fungal_orientation_catalog(empty))==0L,
  isTRUE(fungal_orientation_validate(snapshot,variants[FALSE,,drop=FALSE])),
  identical(serialize(snapshot,NULL),before))
cat("Fungal opposite-orientation candidates, exact-query parity and provenance guards passed.\n")
