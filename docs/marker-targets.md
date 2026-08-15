# Markers are not organism groups

The Atlas treats the molecular marker and the intended biological target as
separate, many-to-many dimensions. A locus can support several organism groups,
and an organism group can be investigated with more than one locus.

| Marker | Common primary targets | Important complementary uses |
|---|---|---|
| COI | Metazoa and especially arthropods | Fish and zooplankton when suitable references and primers exist |
| Mitochondrial 12S | Vertebrates, especially fish and birds | Broader vertebrate eDNA with assays such as 12S-V5 |
| Prokaryotic 16S | Bacteria and Archaea | Plant-associated microbiomes and host-associated communities |
| Animal mitochondrial 16S | Multiple metazoan lineages | Fish, arthropods, and zooplankton with taxon-specific reference panels |
| 18S | Protists, microbial eukaryotes, and plankton | Broad eukaryotes, fungi, metazoans, and degraded-template surveys |
| Fungal ITS | Fungi | Species-oriented fungal community profiling |
| 28S LSU | Fungi, phytoplankton, and other eukaryotes | Phylogenetic community profiling and complementary protist/metazoan evidence |

These links describe established use or a defensible complementary role. They
do not mean that every primer targeting the marker covers every listed group.
Primer sequence, variable region, amplicon length, reference database, and
assay chemistry remain separate evidence layers.

## 12S distinction

Mitochondrial 12S is predominantly a vertebrate eDNA marker in the recent
metabarcoding literature. MiFish and Teleo are fish-oriented; MiBird and BirT
are avian-oriented; 12S-V5 is broader across vertebrates. Thus “12S” alone is
not a sufficient assay description.

Representative sources:

- [MiFish](https://doi.org/10.1098/rsos.150088)
- [Teleo](https://doi.org/10.1111/mec.13428)
- [MiBird](https://doi.org/10.1038/s41598-018-22817-5)
- [BirT](https://doi.org/10.1002/edn3.70255)

## 16S distinction

The Atlas keeps bacterial/archaeal 16S rRNA separate from animal mitochondrial
16S. They are different reference ecosystems and biological questions despite
sharing the shorthand “16S.” For prokaryotic assays, the exact oligonucleotide
sequence and variable region are required because historical primer names such
as 515F/806R refer to multiple sequence versions.

Representative source: [Klindworth et al. 2013](https://doi.org/10.1093/nar/gks808).

## 28S distinction

28S LSU is not a fungi-only marker. LR0R/LR3 and related assays are used for
fungal phylogenetic community profiling, while D1R/D2C supports marine
phytoplankton, dinoflagellates, protists, and broader eukaryotic work. The Atlas
therefore requires target-group evidence at the primer-pair level.

Representative sources:

- [Fungal LSU comparison](https://doi.org/10.1093/femsec/fiv153)
- [PHYTOPK28-D1D2](https://doi.org/10.17632/mndb4h87yg.1)

## Reference geography

Geographic comparisons describe the sequences represented in a reference
database. They are never labelled regional amplification probability. The
Region Comparison view reports located/all denominators, taxonomic composition,
within-order comparisons, primer-position base differences, and exact source
accessions so a geographic pattern can be audited rather than overinterpreted.
