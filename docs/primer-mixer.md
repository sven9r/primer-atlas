# Shared primer exploration

`R/primer_mixer.R` supplies a session-local engine for catalog/custom primer
combinations on explicit nucleotide templates. It leaves every frozen source
catalog, reference snapshot and precomputed score unchanged.

The selector catalog normalizes existing Atlas COI, PR2 18S and UNITE components
and the source-backed fish assay supplement. A row retains its exact source
sequence, citation, target scope and original catalog provenance. U is modeled
as T and inosine as N, with a visible normalization note; these are nucleotide
model assumptions. Custom components receive a deterministic sequence identity
and no inherited publication or empirical validation.

Choose up to five forward and five reverse oligos. Every forward × reverse
combination is retained, even if no product can be mapped. Newly combined
components are candidate combinations; citing the individual components does
not establish that the combination was published or experimentally validated.

## Template matching

Each raw input template has an explicit marker and identity. The engine never
maps a fish assay onto a bundled arthropod sequence as if that represented fish.
The compact fish starter panel is a coordinate/sequence example, with the scope
and provenance described in its manifest; it cannot establish broad fish
coverage. A user may supply up to 100 unaligned FASTA templates totaling at most
200,000 bases. FASTA headers are retained as record identities and are not
interpreted as taxonomy.

The full-length oligo slides over its own template, using IUPAC-compatible
concrete A/C/G/T reference bases. Reverse oligos are searched as their reverse
complements. The maximum substitution count is explicit. The engine API
defaults to zero; the app's visible site-search control starts at up to two
substitutions per oligo, and exact products remain counted separately.
No indels, partial-oligo alignments or annealing/PCR probabilities are estimated.
The best window is accepted only when it is unique. Equally good windows and
competitive ambiguous-base windows are withheld. The overview scan is bounded
to 2.5 million component × template bases; larger requests must be narrowed.

Mismatch positions are expressed in the supplied oligo's 5′→3′ orientation.
Counts in the terminal five bases and a terminal-base mismatch flag are
reported separately. For a reverse oligo, its 3′ end is the leftmost coordinate
of the binding site on the displayed forward template.

A mapped product requires both qualifying sites to be inward and nonoverlapping.
Coordinates start at one on the selected raw template:

- Inclusive product bp = reverse site end − forward site start + 1.
- Internal bp = reverse site start − forward site end − 1.

Length filters alter the eligible map rows while preserving the original
requested combinations, measured coordinates and complete per-template status
table. The shared COI-style renderer receives coordinates for one selected
template; it never draws pooled or median positions as a physical molecule.
Product colors use the shortest matching project planning budget: Illumina
through 500 bp, PacBio through 2,000 bp, Nanopore through 3,500 bp, then beyond
the selected budgets. These tiers use inclusive raw product length and do not
measure assay performance or impose universal instrument specifications.

## Taxonomic identity

Versioned NCBI taxonomy supplements join to bundled source templates by exact
template ID, with redundant marker/group/accession checks. Canonical hierarchy
columns and the taxonomy source/version/archive hashes accompany the original
sequence source fields; the join changes no raw sequence or coordinates.
Before use, starter template/taxonomy CSV bytes are checked against the frozen
manifest SHA-256 and byte counts. Taxonomy supplement identity, row count and
source hashes are checked too. A missing or corrupt declared artifact stops
loading with an integrity error; it cannot silently substitute another panel.
Taxonomy without a verified rank remains `Unresolved`. A custom FASTA header,
even one identical to a known accession, never inherits a bundled template's
taxonomy. Counts describe the actually selected templates, including unresolved
records; they do not imply a representative inventory of the animal group.

For separate fungal 18S and 28S choices, `primer_mixer_fungal_templates()` slices
only independently bounded source gene intervals from the selected frozen
rDNA version. V3 supplies three SSU intervals (two partial) and four LSU
intervals (three partial); disjoint pieces remain separate templates. Source
coordinates, landmark IDs, partial/projected flags, hashes and frozen taxonomy
accompany every slice. V1's curated display anchors are withheld from this new
gene extractor; the original historical rDNA view remains unchanged. UNITE
components whose reported binding loci are SSU or LSU have stable selector
aliases under 18S/28S, retaining original identities, chemistry, citations and
assay context. An alias supplies no new assay-validation claim.

## Status interpretation

`exact_product` means both unique sites match under the declared nucleotide
model. `mismatch_product` means both unique sites satisfy the selected
substitution threshold. Neither status establishes PCR success.

`one_eligible_site`, `template_too_short`, `no_eligible_sites`,
`non_inward_or_overlapping_sites`, `ambiguous_sites`, `template_unavailable` and
`marker_mismatch` remain distinct. A missing match cannot determine whether a
reference is truncated, the primer differs from that taxon, or another cause
prevents a full site from being represented. It is not labeled PCR failure.

The API returns `components`, every requested `pair`, individual `sites`, all
pair × template `products`, per-pair denominator `summary`, explicit `templates`
and `settings`. `primer_mixer_overview_data()` adapts measured, length-eligible
products for `primer_overview_map()`. `primer_mixer_sites()` supplies the
individual-primer overview before combinations are selected.

`R/primer_mixer_exports.R` supplies a complete long result table with oligo
identities, published/source sequences, parent derivations, citations, raw
template identity/version/SHA-256, source/model/partial flags, denominators and
the actual substitution, length and platform settings. The SVG is portable
with an explicit XML namespace and numeric dimensions. Its visible caption
retains candidate scope, template provenance and oligo citations; embedded JSON
also records the exact displayed raw sequence, full component provenance,
requested combinations, displayed products, template-panel hashes and method.
The vector artwork uses the selected template's raw coordinates without
changing the input snapshots.

Validation: `Rscript tests/test_primer_mixer.R` checks raw coordinate arithmetic,
Cartesian requests, length filtering, mismatch thresholds, both 3′ orientations,
repeated and ambiguous sites, partial/missing templates, marker isolation,
source sequence preservation and the shared map adapter.
`Rscript tests/test_primer_mixer_exports.R` checks actual export source metadata,
visible citations, hashes, XML bounds and preservation of derived candidates.
`Rscript tests/test_primer_mixer_taxonomy.R` verifies source-backed hierarchy
joins, untouched sequence provenance, group-scoped denominators and the absence
of inferred taxonomy for user FASTA.
