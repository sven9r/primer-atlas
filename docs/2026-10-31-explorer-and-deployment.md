# Shared marker exploration and taxonomic fit

The working branch uses one investigation flow: organism group → marker →
individual-primer overview → Mix & match → amplicons → taxonomic fit. The
previous published-pair maps, broad UNITE evidence, frozen alignments and
full arthropod COI scores remain available in explicit reference sections.

Arthropods default to COI. Fish can choose mitochondrial 12S, COI, nuclear 18S
and 28S, with mitochondrial 16S also supported by source-backed components and
actual templates. The existing animal groups also expose 18S and 28S as complementary
choices. Marker availability does not establish that every assay suits the
chosen group.

## Primer and template evidence

The new source-backed supplements contain 32 oligos: 22 fish/broad-eukaryote
components and ten additional bird/arthropod components. Sources include
MiFish, teleo, vertebrate 12S-V5, MiBird, fish and bird COI barcoding, SeaDNA
COI, broad eukaryotic 18S, food-fish 28S and braconid-wasp D2F/D3R and source-backed mitochondrial 16Sar/16Sbr. Their
original target scope, exact gene-specific sequences and citations remain
visible. Adapters and tags are not treated as binding sequence.

The fish starter panel contains 13 marker templates from seven accession
versions. The additional animal starter panel contains 11 templates for
arthropod 18S/28S/mitochondrial 16S, bird 12S/COI/18S/28S and Daphnia
zooplankton COI/18S/28S/mitochondrial 16S. The existing Drosophila yakuba COI
reference remains unchanged. Together these provide 25 explicit animal marker
templates, including that existing COI reference. All 22 currently exposed
animal-group/marker routes have source-backed components, templates and
reconciled taxonomic-fit results; this does not establish exhaustive coverage. Separate NCBI taxonomy supplements identify the
nine source species and preserve the exact template, accession, marker,
group, TaxID, hierarchy, version and archived response checksum. Predicted
MODEL and partial records keep those flags.

These are bounded mapping examples, not representative fish/bird/rDNA
coverage panels or new monthly marker releases. The independent full-state
COI stage retains 13,150 records / 12,267 verified COI records, with at least
105 verified records per supported group. Broad UNITE evidence retains its
own 15,402-record denominator; long-rDNA geometry retains its 21-template
cohort and historical versions. None of these denominators substitute for
another.

## Mix & match and taxonomic fit

Up to five forward and five reverse components produce all requested Cartesian
combinations. Custom oligos and source-backed catalog oligos can be mixed.
Clicking a primer cycles through its original direction, an explicit opposite
reverse-complement candidate, and removal. Modified oligos retain parent
provenance and acquire no inherited experimental validation.

The visible substitution threshold defaults to two per oligo; exact matching
and thresholds from zero to five remain selectable. The engine's API defaults
to exact matching. Unique best full-length sites, ambiguous or repeated sites,
3′ differences, inward geometry and raw encoded lengths remain distinct.
Platform and length filters change displayed eligibility while retaining
every requested pair and per-template result.

The taxonomic fit page uses those same results for each marker and animal
group. Cells group records by their full parent lineage at phylum, class,
order, family, genus or species. All selected templates in a group form the
denominator; exact, substitution, partial, ambiguous and missing/unsupported
states partition it. Both-site denominators are separate. Clicking a cell
opens accession-level evidence. Source taxonomy is never inferred from a
FASTA header. Optional session taxonomy CSV is explicitly labelled
user-provided metadata.

The publication report provides a formatted HTML table and an evidence ZIP.
The CSV exports retain sequence identities, chemistry normalization, original
and modified roles, sources, template hashes, taxonomy, settings and exact
counts. The standalone SVG contains a visible source caption plus structured
metadata. The original broad-reference publication export remains available.

## Publishing this version

The hosting services already exist. This change does not require a new host.

1. Pass the required GitHub check and obtain the repository's required
   approving review for [PR 22](https://github.com/sven9r/primer-atlas/pull/22).
   Merge the reviewed branch into `main`.
2. In the existing Posit Connect Cloud content, confirm repository
   `sven9r/primer-atlas`, branch `main`, entrypoint `app.R`, and the regenerated
   `manifest.json`. Republish after merge, or let an already enabled automatic
   publish-on-push setting do so. Posit's
   [Shiny deployment guide](https://docs.posit.co/connect-cloud/how-to/r/shiny-r.html)
   documents these controls.
3. Confirm `ATLAS_MANIFEST_BASE_URL` points to the existing public R2
   `/releases` URL. Activate the growing full-state COI pointer using the
   controlled first-activation command in [DEPLOYMENT.md](../DEPLOYMENT.md),
   after review/merge. Publication verifies the full candidate before promotion.
4. Test the anonymous public app: Fish → each marker, catalog/custom mixing,
   taxonomy drill-down, exact downloads, fungal historical versions and the
   full COI source identity. Confirm the public landing page's Launch link.

Until these external steps finish, the local preview and the public release
must be described separately. Broad reference acquisition for additional
animal markers remains a scientific data milestone, independent of the
completed shared interface and session-template evaluation.
