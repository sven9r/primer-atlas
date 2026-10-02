# Hawaii-Madagascar insect COI feasibility audit

BOLD Data Portal queried: 2026-09-23. This is a bounded feasibility audit, not a global range census or a primer-performance estimate.

## Source and scope

The script `scripts/audit_regional_coi_bold.R` reproduces the BOLD Portal API summary and record queries. Hawaiʻi uses `geo:province/state:Hawaii`; Madagascar uses `geo:country/ocean:Madagascar`. Order-level summary counts are in `data/derived/regional_coi_audit/bold_summary_counts.csv`. Record-level inspection is limited to Crambidae and recorded in `data/derived/regional_coi_audit/bold_crambidae_record_audit.csv`. The latter contains IDs and metadata, not nucleotide strings. Query tokens expire, so rerun the script for a fresh snapshot.

BOLD API documentation: https://portal.boldsystems.org/api . Its summary counts describe records matching indexed taxonomy and geography fields; they are not deduplicated species counts.

## Bounded Crambidae record audit

A preliminary sequence QC pass requires a COI-5P record, the expected indexed region field, at least 500 canonical A/C/G/T bases, and at most 2% other non-gap bases. This establishes a useful COI fragment; it does not establish that either primer binding site is present.

| Region | Unique records | COI-5P | Region field confirmed | Valid coordinate pair | Sequence QC pass | QC records with INSDC accession | Genera in QC pass | Species labels in QC pass |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| HI_ARCH | 283 | 255 | 283 | 151 | 220 | 91 | 14 | 33 |
| MG_COUNTRY | 522 | 522 | 522 | 514 | 517 | 272 | 90 | 88 |

Duplicate BOLD record IDs removed across query results: 0.
Shared genus labels in the QC pass (7): Omiodes, Herpetogramma, Parapoynx, Uresiphita, Chilo, Diaphania, Maruca.
Shared species labels in the QC pass (3): Herpetogramma licarsisalis, Parapoynx fluctuosalis, Maruca vitrata.

| Shared species label | Hawaiʻi QC records | Madagascar QC records |
|---|---:|---:|
| Herpetogramma licarsisalis | 6 | 2 |
| Parapoynx fluctuosalis | 1 | 5 |
| Maruca vitrata | 1 | 5 |

Two shared species have only one QC-passing Hawaiʻi record each; their within-species regional comparison is not yet supported by replication. BOLD process and record IDs preserve provenance for records without an INSDC accession.

These are BOLD taxonomic labels, not reconciled accepted species concepts. Shared names require taxonomic review before a within-species regional comparison. Locality precision varies, and coordinates need spatial validation against the defined regions.

## Primer binding and endemicity readiness

No primer pair was scored in this audit. BOLD COI-5P records have variable starts, ends, gaps, and PCR primer metadata; a 500-base barcode can still omit an assay's binding sites. The next scoring pass must align each sequence to a COI reference, locate both documented primer sites, report the scorable denominator per site and pair, then apply the existing mismatch model. The field `pilot_sequence_qc` must not be interpreted as primer scoreability.

The new `data/catalog/range/taxon_region_evidence.csv` remains a separate, reviewed taxon-range catalog. Sequence locality, absence of a record elsewhere, and shared species labels do not establish endemicity.

## Recommended next slice

Use Crambidae and one documented internal COI primer pair. Review the shared species labels and coordinate quality first. Then map primer binding sites in the QC-passing sequences and compare (1) all Crambidae, (2) shared genera, and (3) shared species where each region has enough independent specimens. Show sample size and missing-site reasons at every level.
