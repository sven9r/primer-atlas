# LSU domain annotations

Snapshot **v3-c6fa87f256e1** adds source-projected LSU domains and explicitly partial NCBI gene intervals to the unchanged V2 reference sequences, alignment and primer results. V1 and V2 remain byte-identical and selectable.

| Accession | Mapped domain intervals | Source LSU coverage |
|---|---:|---|
| PZ962702.1 | 9 | One complete 28S feature |
| PZ276773.1 | 10 | Two separate partial 28S copies |
| FJ430740.1 | 1 | One partial 28S feature |

There are **20 intervals and 10 distinct domain/subdomain labels** across 3 of the 21 retained templates. Five previously suppressed partial gene intervals are also retained with a partial label. The remaining templates do not have separately bounded LSU features suitable for this projection. Missing domain labels do not establish biological absence. Thirty-six source/feature mapping cases were withheld; D6, D7b and D12 remain unresolved in this snapshot.

## Sources and support

The classic D-domain scheme is distinct from other LSU variable-region/expansion-segment schemes. [NCBI X00525.1](https://www.ncbi.nlm.nih.gov/nuccore/X00525.1) supplies machine-readable mouse 28S D-domain features from Hassouna et al. Its deposited D10 interval overlaps D9 and is withheld without silently changing its coordinates. The D10 projections instead use the numeric human M11167 range in [Mallatt and Sullivan (1998), Table 3](https://doi.org/10.1093/oxfordjournals.molbev.a025897). D7a and D7b remain separate sourced components.

`scripts/map_fungal_domains.py` works on each accession's independently annotated LSU feature, keeping repeated rDNA copies separate. MAFFT localpair and genafpair projections must agree. Boundary windows contain 32–64 source bases and stop before adjacent eligible source domains. Each flank requires at least 75% paired coverage, 70% paired identity, and a 0.08 normalized penalized-score margin over gap-aware competing matches elsewhere in the same feature. Target ambiguities and unsupported boundaries are withheld. These checks support a computational homology projection; they are not independent experimental validation of fungal domain boundaries.

Landmark intervals use **one-based inclusive raw-template positions**. Nested method/flank evidence uses explicitly labeled feature-relative coordinates. Domain edges include target indels within the projected source interval. Alignment columns are derived from the unchanged frozen raw-to-column maps.

The snapshot archives the original source responses, reference sequences, domain coordinate tables, mapping alignments, support metrics, rejected cases and generator code. Its manifest maps archived filenames to their original names and SHA-256 values. `scripts/build_fungal_domain_snapshot.R` creates annotation-only versions atomically and reuses identical content rather than creating another ordinal on replay.

## Investigation behavior

The raw-bp and gapped-alignment views share gene/domain annotations and selected-product coverage. Selecting a domain focuses its raw bounds with 1,000 bp left and 2,000 bp right context, clipped to the selected record. Changing comparison pairs preserves explicit region focus. Gene-copy controls show their own coordinate intervals. Source details and calculation explanations can be opened separately.

For example, LR0R + LR3 on PZ962702.1 spans raw 2446–3103, **658 bp including primer sites**, and fully covers the projected D1 and D2 intervals. Clicking D1 focuses raw 1543–4690. The corresponding gapped column window is 7762–12846; switching axes does not change the raw region or product length.

Monthly recovery restores hosted history and selects the newest verified immutable snapshot from the combined bundled and restored versions. An older hosted panel or catalog-only release therefore cannot replace the newer bundled domain annotations. Recovery verifies complete earlier accession-version sequence retention and landmark/taxonomy accession ownership, and stops on corrupt evidence, unavailable hosted state or an ambiguous newest version. Monthly retention preserves extended annotation columns and their string types, including `parent_landmark`; new accessions require their own annotation evidence. A growing reference panel does not automatically inherit the current templates' domain positions.
