# Focused investigation pages

The Primer Atlas now separates primer geometry, taxonomic reference support and reference records. Primary navigation is Start, Primer explorer, Taxonomic fit and Reference library. Location comparisons, temperatures and sources remain available under More, with their existing route identities and evidence readers preserved.

Start offers three marker cards: animal COI, fungal rDNA and the PR2 18S catalog. The explorer defaults to its map and selected products; filters, candidate tables, alignment internals, provenance and exports open on demand. Support logos and update diagnostics are available under Sources rather than occupying every investigation page. The COI fit view starts with five primer columns and a closed reading guide; selecting a heatmap cell opens its comparison inspector.

## Fungal reference support

Taxonomic fit has explicit Animal COI and Fungal rDNA modes. Fungal mode uses the selected immutable snapshot's raw sequences, taxonomy and exact-site product records. Phylum and class are always available with unresolved assignments retained; order and family are offered only when the snapshot has actual resolved labels at those ranks.

Each clickable taxon card shows exact product support as a count divided by **all retained reference templates in that hierarchy-scoped group**. Same-name taxa under different parent lineages remain distinct. Drilling down opens accession-level statuses and reasons. These states partition the denominator:

- exact product;
- one exact site / partial pair;
- ambiguous or otherwise unresolved sites;
- no exact sites found in the record;
- exact sites whose geometry or minimum-length rule withholds a product;
- missing, inconsistent or unsupported evidence.

Absence of an exact site in a retained record does not establish PCR failure or absence of the locus. This bounded fungal panel is not a representative survey of all fungi, and the view does not borrow COI mismatch penalties or amplification probabilities. The taxonomy-group count explicitly includes unresolved assignments.

Platform ceilings use encoded ungapped reference spans including both primers: 500 bp Illumina, 2,000 bp PacBio and 3,500 bp Nanopore. Selecting a ceiling never changes the exact-support numerator or its panel denominator. The planning count instead reports eligible encoded spans amongst exact-supported products. Internal ambiguous bases are disclosed; an encoded span cannot determine an assembly gap's unknown physical length.

## Context and historical versions

Explorer → fit and fit → explorer carry the exact version and pair through asynchronous Shiny control updates. An empty explorer selection produces a choose-combination notice and retains the current route, rather than opening unrelated evidence. User selections and filter resets cancel stale handoff acknowledgments.

The fungal reference library displays the full frozen panel from the fungal fit view. Metadata downloads retain each accession, source fields and snapshot ID. Taxon-count and reference-evidence CSVs use the displayed version, pair, rank and search scope. A historical V1 view uses its one retained record and unresolved taxonomy; the current V3 view uses its 21 records. Frozen snapshots and monthly acquisition/scoring logic are unchanged by this interface work.

## Validation

`tests/test_fungal_fit.R` reconciles all 18,564 V3 pair-template rows across 884 pairs and verifies status partitions, hierarchy ownership, unavailable evidence and historical behavior. Rendered module tests cover count bars, taxon selection, independent ceiling changes, empty searches, version changes and handoff release. Real application navigation tests cover both fit modes, historical pair/version round trips, actual versioned library CSV callbacks, preserved COI state and the original four mapped 18S defaults. Explorer tests cover delayed client acknowledgments and later direct selections.
