# Primer Atlas day-end status — 2026-09-23

Reviewed at 20:02 UTC. These are local, uncommitted results. The range, fungal
ITS, and navigation tasks finished; the scheduled regional COI automation is
paused with no run record, although a local audit script and outputs exist.
None of those four source tasks was active when checked. The navigation prototype is not
currently listening at `127.0.0.1:4387`. Nothing in this day-end review was
deployed or published.

| Work | Verified local output | Limit |
| --- | --- | --- |
| Range and endemism | `data/catalog/range/regions.csv` has 11 regions; `sources.csv` has 4 source records; `taxon_region_evidence.csv` has 0 assertions. The schema and decision rules are in `taxon_region_evidence.schema.json` and `docs/range-evidence-methods.md`; `README.md` links them. | No taxon has a reviewed range or endemicity assertion. The existing 67,352-row COI geography file has only 500 country/territory assignments, no USA locality mentioning Hawaiʻi, and one Madagascar row. |
| Hawaiʻi–Madagascar COI | `scripts/audit_regional_coi_bold.R`, `docs/regional-coi-feasibility.md`, `data/derived/regional_coi_audit/bold_summary_counts.csv`, and `data/derived/regional_coi_audit/bold_crambidae_record_audit.csv` record 12 BOLD summary queries and 805 unique Crambidae records: 283 Hawaiʻi and 522 Madagascar. Preliminary sequence QC passes 220 and 517, respectively. Seven genus labels and three species labels are shared among passing records. | The automation itself did not run. These BOLD labels are unreconciled; coordinates vary in quality. No primer sites or pairs were scored, and sequence locality cannot establish native range or endemicity. |
| Fungal ITS | `scripts/build_fungal_its_pilot.py`, `scripts/score_fungal_its_pilot.R`, `scripts/run_fungal_its_pilot.sh`, and `docs/fungal-its-sequence-pilot.md` document a reproducible 48-accession UNITE 10.0 panel across eight phyla. Ignored local outputs `data/derived/its_fungal_pilot/site_scores.csv` and `data/derived/its_fungal_pilot/pair_scores.csv` have 432 site rows and 384 pair rows. ITS3 is scorable for 48/48, ITS4 for 4/48, and each of four ITS2 pairs for 4/48. | The accession panel is deliberately small. Most terminal sites are absent or cannot be verified against the single reference. These are site-availability results, not clade-wide primer performance. The accession-level outputs are Git ignored and not part of a reviewable commit. |
| Navigation | `app.R`, `tests/test_marker_runtime.R`, and `docs/navigation-prototype.md` add a local Start here panel, task routes, and marker-readiness copy. Six targeted tests and a local browser walkthrough passed. | The wording “reference fit by taxon” needs review for possible confusion with PCR success. This prototype has not been deployed. |

The work overlaps in one checkout: range and COI both changed `README.md`,
the COI audit changed `.gitignore` to retain its CSVs, and navigation changed
`app.R` and one test. All four streams remain uncommitted together. The fungal
accession-level files are ignored, so committing the visible scripts and note
alone would not preserve their output tables. `PROJECT_STATUS.md` still holds
an older operational overview; the COI public release gate remains unresolved
there. The present note records today's local research and UI state only.

## Next three actions

1. Review and commit the four scoped changes with their intended provenance;
   decide whether the fungal accession-level outputs need a tracked compact
   artifact or a documented reproducible regeneration path.
2. For a Crambidae primer pilot, reconcile accepted species names and
   coordinates, then map both binding sites in the 737 QC-passing records and
   report scorable denominators before any regional comparison.
3. Independently annotate ITS boundaries and terminal sites, then review the
   navigation labels in the local app before considering a broader panel or
   public UI change.
