# October 11 implementation and validation plan

This ten-step plan covers the requested fungal combination explorer, appearance,
sequence acquisition and recoverable monthly releases. Status describes evidence
available on the working branch; public activation is a separate final step.

| Step | Work and completion evidence | Status |
| --- | --- | --- |
| 1 | Freeze the existing fungal alignment, templates, classes and source identities so the original view remains reproducible. V1 is retained with its original content hashes. | Complete |
| 2 | Replace the eight-pair starting view with an interactive library of 20 fungal primer-combination examples, with editable selections and source evidence. | Complete |
| 3 | Compare platform eligibility using the requested project budgets: Illumina 452–500 bp, PacBio 2,000 bp and Nanopore 3,500 bp. Apply platform colors consistently in the explorer and alignment. | Complete |
| 4 | Put selected combinations directly onto the shared fungal alignment, retaining the distinction between alignment columns and ungapped product length. | Complete |
| 5 | Add recoverable fungal versions with exact templates, taxonomy, sources and alignments. V1 and V2 remain selectable; V2 adds 20 accession-versioned NCBI records to the baseline. | Complete |
| 6 | Refresh the landing page and application styling, then check the real narrow-browser layout, filters, alignment and readable download controls. | Complete |
| 7 | Acquire and score at least 10,000 verified COI records overall and 100 in each supported group. The second hosted stage has 12,267 verified records and minimum 105 across 20 orders plus Acari/Collembola. | Complete |
| 8 | Demonstrate a fresh cloud worker restoring a hosted release and extending every group. Both cloud builds passed. Independent anonymous comparison verified all 13,005 prior sequence/header hashes and 431 scientific source archives; the second build added 145 records and all 1,087 hosted artifacts passed read-back. | Complete |
| 9 | Verify that the deployable application can read the hosted cohort, retain a valid cohort after a failed update, and export exact sequence and growth evidence. The actual clean bundle passed a same-session refresh between both hosted versions, with complete QA and exact CSV/FASTA downloads; the live preview shows the second cohort. | Complete |
| 10 | Complete checks, publish reviewable code and evidence, then activate the new public app and monthly full-state pointer after the repository's required approving review and merge. | PR #22 ready for review; activation pending |

Detailed COI counts, hosting identities, recovery evidence and remaining scientific
scope are recorded in [the dated validation record](2026-10-11-coi-growth.md).
Primer fit and platform eligibility are reference evidence, not measured PCR
success. Family inventories, aquatic eligibility and accession geography are
subsequent scientific work under [the growing-reference contract](growing-reference-contract.md).
