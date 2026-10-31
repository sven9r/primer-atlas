# Animal starter mapping panel

`data/pinned/ANIMAL_STARTER/v1/` supplies eleven source-backed examples for
marker exploration. Arthropod, bird and zooplankton choices load their own
explicit records. Ecological zooplankton membership is source-backed for the
named Daphnia species, rather than inferred from broad arthropod membership.

| Group | Marker | Source record | Exact retained feature |
| --- | --- | --- | --- |
| Arthropods | 18S | [Drosophila melanogaster M21017.1](https://www.ncbi.nlm.nih.gov/nuccore/M21017.1) | 1–1,995; complete source 18S feature |
| Arthropods | 28S | [Drosophila melanogaster NR_133562.1](https://www.ncbi.nlm.nih.gov/nuccore/NR_133562.1) | 1–3,970; FlyBase-curated RefSeq RNA record |
| Arthropods | Mitochondrial 16S | [Drosophila yakuba NC_001322.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_001322.1) | complement(12,756–14,081); 1,326 bp oriented by its annotated strand |
| Birds | Mitochondrial 12S | [Gallus gallus NC_001323.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_001323.1) | 1,297–2,272; 976 bp |
| Birds | COI | [Gallus gallus NC_001323.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_001323.1) | 6,645–8,192; 1,548 bp |
| Birds | 18S | [Gallus gallus AF173612.1](https://www.ncbi.nlm.nih.gov/nuccore/AF173612.1) | 1–1,737; complete source 18S feature |
| Birds | 28S | [Gallus gallus XR_005840241.1](https://www.ncbi.nlm.nih.gov/nuccore/XR_005840241.1) | 1–4,444; predicted RefSeq rRNA MODEL |
| Zooplankton | COI | [Daphnia magna NC_026914.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_026914.1) | <1,313–2,849; 1,537 encoded bp; source has a partial 5′ boundary |
| Zooplankton | Mitochondrial 16S | [Daphnia magna NC_026914.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_026914.1) | complement(12,366–13,738); 1,373 bp oriented by the source strand |
| Zooplankton | 18S | [Daphnia pulex AF014011.1](https://www.ncbi.nlm.nih.gov/nuccore/AF014011.1) | 1–2,293; complete source 18S feature |
| Zooplankton | 28S | [Daphnia magna EU370436.1](https://www.ncbi.nlm.nih.gov/nuccore/EU370436.1) | <1–>3,823; source explicitly partial at both ends |

The original NCBI XML responses, exact request parameters, hashes, accession
versions, source taxonomy and feature intervals are archived with the table
and FASTA. The source taxonomy explicitly identifies Arthropoda, Aves or the named Daphnia lineage, and
that lineage check is part of extraction and validation. All local coordinates
start at one on the displayed template; source coordinates remain separate.

The original M21017.1 record contains both a complete and terminal partial 18S
repeat. Only the complete 1–1,995 interval is included. Its discontinuous
`join(3288..5075,5121..7232)` 28S feature is withheld: concatenating it would
not preserve genomic PCR span. The separate NR_133562.1 record is retained as
its own identified RNA reference, not substituted into M21017 coordinates.
Its underlying [genomic source NW_007931121.1](https://www.ncbi.nlm.nih.gov/nuccore/NW_007931121.1)
is archived too. All 3,970 bases exactly match its contiguous annotated rDNA
feature at 70,955–74,924, which explicitly identifies NR_133562.1. This verifies
the encoded span for that repeat copy; it does not validate PCR or other copies.
The chicken 28S prediction is explicitly flagged; no claim of experimental
transcript verification or universal species coverage is made.

Freshwater zooplankton membership has taxon-specific primary sources:
[Colbourne et al. 2011](https://doi.org/10.1126/science.1197761) describes
D. pulex as a freshwater lake/pond grazer and prey species;
[Orsini et al. 2016](https://doi.org/10.1038/sdata.2016.30) describes
D. magna's pelagic freshwater ecology. Their citation keys, URLs and ecological
scope remain in the template table and manifest. They establish ecological
membership, not primer validation. Missing terminal sites in partial records
remain distinct from observed primer mismatches.

These are small mapping exemplars, not an inventory of animal diversity. A
primer match counts support on these actual records and cannot establish
order, family or group-wide PCR coverage. RNA-record spans are encoded
reference intervals; mapping alone does not verify a genomic PCR product,
splicing, introns or rDNA repeat variation. No SSU variable-region or LSU
D-domain boundaries are inferred in this panel. The panel is not part of the
monthly release matrix and does not replace the broad arthropod COI release.

`data/catalog/animal_assay_primers.csv` contains ten oligos. It adds the four original BirdF1, BirdR1,
BirdR2 and BirdR3 oligos from [Hebert et al. 2004, Materials and
Methods](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0020312).
The source documents BirdF1 with each alternative reverse primer; BirdR3 was
used for two Falco species. These are tissue-barcoding assays, not evidence of
water-eDNA specificity. It also includes the exact MiBird-U-F/MiBird-U-R
gene-specific tails from [Ushio et al. 2018, Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC5851996/#Tab1).
MiSeq adapters and six random bases are excluded. The design deliberately
permits recovery of other vertebrates, so it is not described as exclusively avian.
The arthropod 28S D2F/D3R pair is reproduced from the primary [Wild et al. 2013,
Table 2](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0074837#pone-0074837-t002),
where its original attribution is Mardulyn and Whitfield (1999). That study
documents braconid wasp use, not universal arthropod recovery.
The mitochondrial 16Sar/16Sbr pair is reproduced exactly from [Hiller and
Werding 2019, Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC6775170/),
which attributes the oligos to Palumbi (1996). Its source study used marine
porcelain crabs; exploratory mapping on Daphnia does not inherit that
validation. Mitochondrial 16S is distinct from prokaryotic 16S.
Existing broad-eukaryote 18S, conserved 28S and
vertebrate 12S oligos can be explored with their original source scopes intact;
choosing a different organism group does not confer new assay validation.

Rebuild from the archived XML using the offline
`scripts/build_animal_starter_panel.py` extractor. Run the checks with:

```sh
python3 tests/test_animal_starter_panel.py
```

## Source-verified hierarchy for the animal fit view

`taxonomy.csv` has twelve exact template identities: the eleven starter rows plus
the separately maintained existing COI reference `NC_001322.1_COX1`
([Drosophila yakuba NC_001322.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_001322.1)).
Its new archived source record identifies NCBI taxon 7245; its annotated COX1
slice exactly matches the existing 1,536-bp reference FASTA. This taxonomy row
does not add or change a reference sequence. The manifest explicitly lists the
external template identity separately from its eleven panel members.

The independent supplement version is `NCBI-Taxonomy-2026-10-31-v1`. Kingdom
through species are copied from actual ranked nodes in the archived NCBI
Taxonomy response, joined through each accession source's taxon-ID qualifier.
All identities retain organism, marker, group, accession version, taxon ID,
source/archive hashes and retrieval date. Missing or conflicting ranks remain
`Unresolved`; no hierarchy is invented for custom FASTA headers.

The four arthropod template identities resolve to Diptera/Drosophilidae,
the four bird template identities to Galliformes/Phasianidae, and the four
ecological zooplankton identities to Diplostraca/Daphniidae. These narrow exemplar collections do not establish
animal-wide taxonomic coverage. A fit summary must group the actual selected
marker/template cohort and retain its exact denominator, rather than counting
the entire supplement as one marker's evidence.

`scripts/build_starter_taxonomy_supplements.py` archives the taxonomy snapshot
and extends only the manifest and supplemental artifacts. Original template
CSV/FASTA bytes are retained. Run:

```sh
python3 tests/test_starter_taxonomy_supplements.py
```
