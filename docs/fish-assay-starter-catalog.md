# Fish assay starter catalog and mapping templates

This addition supplies actual primer oligos and bounded reference sequences for
the group → marker → mix-and-match workflow. It does not turn a starter panel
into a representative survey of fish or a monthly maintained marker release.
Existing COI, fungal and PR2 frozen snapshots are unchanged.

`data/catalog/fish_assay_primers.csv` contains **22 unique gene-specific oligos**
representing **12 documented assays**. Sequences are written 5′ to 3′ exactly as
published or recorded by the identified PR2 catalog. Sequencing adapters,
barcodes, blocking primers and synthetic tags are excluded. `component_id` is
unique; `documented_pair_ids` retains all published memberships for a component
shared by multiple assays. A new pairing is an exploratory combination, not an
additional published assay.

| Marker | Unique oligos | Documented assays | Source and scope |
| --- | ---: | --- | --- |
| Mitochondrial 12S | 8 | MiFish-U, MiFish-E, teleo, Riaz 12S V5 | [Miya et al. 2015](https://doi.org/10.1098/rsos.150088), [Valentini et al. 2016](https://doi.org/10.1111/mec.13428), [Riaz et al. 2011](https://doi.org/10.1093/nar/gkr732); sequences independently verified in [Collins et al. 2019, Table 1](https://doi.org/10.1111/2041-210X.13276). Target scopes differ. |
| COI | 7 | FishF1/FishR1, FishF2/FishR2, SeaDNA-short, SeaDNA-mid | All four [Ward et al. 2005](https://pmc.ncbi.nlm.nih.gov/articles/PMC1609232/) oligos are preserved. SeaDNA sequences come from [Collins et al. 2019](https://doi.org/10.1111/2041-210X.13276). Tissue-barcoding use and water-eDNA performance remain distinct. |
| Nuclear 18S | 5 | Stoeck V4; Amaral-Zettler 1380F/1510R and 1389F/1510R V9 | [Stoeck et al. 2010](https://doi.org/10.1111/j.1365-294X.2009.04480.x), [Amaral-Zettler et al. 2009](https://doi.org/10.1371/journal.pone.0006372); exact oligos also match the unchanged PR2-primer 2.1.1 snapshot. These are broad-eukaryote assays, not fish-specific assays. |
| Nuclear 28S | 2 | C1′/C2 | Exact oligos and documented food-fish application from [Hallen-Adams 2015](https://doi.org/10.1111/1541-4329.12063). This conserved assay also recovers non-fish material. It is not fish-specific. |

Each component has its own `source_key`, original citation/DOI, sequence
verification URL, target scope and evidence note. Citation keys resolve in
`data/catalog/sources.csv`. The C1′/C2 source gives an approximate 320-bp
description and also describes larger apparent gel products; the catalog does
not convert either value into a computed span. Product lengths must come from
the selected template's actual binding positions.

## Bounded NCBI mapping panel

`data/pinned/FISH_STARTER/v1/` freezes **13 templates from seven accession
versions**, acquired on 2026-10-31 through NCBI Entrez. The manifest retains
original request parameters, compressed original XML responses and SHA-256
checksums. `templates.csv` and `templates.fasta` contain exact source-feature
slices with their original start/end coordinates, strand, sequence hashes,
organism, accession version and source status.

| Record | Organism | Included templates | Scope |
| --- | --- | --- | --- |
| [NC_002333.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_002333.2) | Danio rerio | 12S, 16S, COI | Exact annotated mitochondrial feature slices |
| [NC_002081.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_002081.1) | Gadus morhua | 12S, 16S, COI | Exact annotated mitochondrial feature slices |
| [NC_001950.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_001950.1) | Scyliorhinus canicula | 12S, 16S, COI | Exact annotated mitochondrial feature slices; elasmobranch exemplar |
| [FJ915075.1](https://www.ncbi.nlm.nih.gov/nuccore/FJ915075.1) | Danio rerio | 18S, 702 bp | Source is explicitly partial |
| [FJ710886.1](https://www.ncbi.nlm.nih.gov/nuccore/FJ710886.1) | Salmo salar | 18S, 1,750 bp | Source is explicitly partial |
| [XR_012407110.1](https://www.ncbi.nlm.nih.gov/nuccore/XR_012407110.1) | Danio rerio | 18S, 1,887 bp | Predicted RefSeq rRNA MODEL; not experimentally validated transcript evidence |
| [XR_012403044.1](https://www.ncbi.nlm.nih.gov/nuccore/XR_012403044.1) | Danio rerio | 28S, 4,268 bp | Predicted RefSeq rRNA MODEL; not experimentally validated transcript evidence |

Coordinates on each extracted template start at one. `source_start` and
`source_end` map a local interval back to the archived original record. A
position on one species is not an alignment column or a homologous coordinate
on another species. No SSU variable-region or LSU D-domain boundaries are
invented for these fish records. Source partiality and predicted status must
remain visible in source details and downloads.

A missing site in a partial record cannot establish a primer mismatch or locus
absence. An exact hit, a best mismatch site and an inward product are separate
outputs; none measures PCR success. This selected panel does not establish
fish, order or family coverage. It must not substitute for the growing
arthropod COI reference, the UNITE panel or any future broad fish release.

The three annotated mitochondrial 16S intervals are included alongside 12S and
COI from the same original fish records. Conserved 16Sar/16Sbr oligos are available
in the separate animal assay catalog with their source scope retained; new
fish-wide validation is not implied.

The offline extractor is `scripts/build_fish_starter_panel.py`. It rejects
missing, joined, unsupported-strand or out-of-bounds feature intervals. Its
source hash is in the frozen manifest. Run the source-to-slice, assay identity,
PR2 sequence, checksum and scope checks with:

```sh
python3 tests/test_fish_starter_panel.py
```

## Source-verified taxonomic hierarchy

`taxonomy.csv` supplements all thirteen templates without modifying either template
artifact. Its independent version is `NCBI-Taxonomy-2026-10-31-v1`. Every join
starts with the exact source accession's `source /db_xref="taxon:…"` qualifier,
then reads that taxon ID's archived NCBI Taxonomy XML. Kingdom, phylum, class,
order, family, genus and species use the actual NCBI rank values. Names are not
inferred from FASTA headers or promoted between ranks. Missing or conflicting
ranks remain `Unresolved`.

The thirteen records currently reconcile to six Danionidae, three Gadidae, three
Scyliorhinidae and one Salmonidae templates. These are counts of this selected
multi-marker template collection, not independent species, fish diversity or
PCR recovery. A marker-specific investigation counts only its own templates.
For example, NCBI's class value here is `Actinopteri` for the teleost exemplars;
the supplement does not replace it with a familiar but differently ranked name.

The manifest identifies the supplement, its CSV checksum, archived XML
checksum, retrieval date, builder hash and exact request parameters. Original
NCBI source records prove the taxon join. `taxonomy_source_sha256` and
`raw_taxonomy_xml_sha256` hash the uncompressed TaxaSet response;
`raw_taxonomy_sha256` hashes its stored gzip archive.

The offline supplementation script is
`scripts/build_starter_taxonomy_supplements.py`. Source-to-taxid and exact-rank
checks, hierarchy/count reconciliation, original template retention and
unresolved-rank behavior are validated by:

```sh
python3 tests/test_starter_taxonomy_supplements.py
```
