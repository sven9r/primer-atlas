# Fungal ITS sequence-scoring pilot (local, 2026-09-23)

This pilot advances the fungal ITS catalog beyond the single FN812768.2 map. It
does not feed `scripts/build_release.R`, R2, or the public Shiny app. The
existing map remains a reference-specific geometry display. Its ITS1F, ITS4B,
and LR21 placements are marked `reference_specific_low_confidence` in
`data/derived/marker_pair_geometry.csv`; the pilot does not score those sites.

## Source and acquisition

The scored source is [Abarenkov et al. (2025), *UNITE general FASTA release for
Fungi*, version 19.02.2025, DOI:10.15156/BIO/3301229](https://doi.org/10.15156/BIO/3301229),
listed by [UNITE](https://unite.ut.ee/repository.php) as version 10.0. It is the
dynamic RepS/RefS sequence set (102,137 records in the downloaded FASTA), not
the complete fungal sequence archive. UNITE separately lists the [full
UNITE+INSDC fungal ITS release, version 19.02.2025,
DOI:10.15156/BIO/3301227](https://doi.org/10.15156/BIO/3301227) with 2,069,189
records; that larger source is suitable for a future scaled panel. No inference
about all fungi or all life follows from this pilot.

The acquisition script pins the PlutoF media URL, archive SHA-256
`a340864f947c517b671f1054a261e613a6619ef88dd8000de1699f5e1f21413b`,
and extracted FASTA SHA-256
`1bb35a664b20ef03484fe74e857929ccd537d9d51c9eef0305b99a1486e7386e`.
Each selected row retains the exact UNITE header, accession as published there,
SH ID, RefS/RepS flag, taxonomy ranks, sequence hash, length, DOI, and version.
UNITE headers do not consistently supply INSDC accession version suffixes, so
the manifest does not invent them. The source labels the records as ITS; exact
ITS1/5.8S/ITS2 boundaries are not separately annotated in this FASTA. Site
and pair outputs give the reference-aligned locus or target region explicitly.

## Bounded panel and scoring rule

Eight named fungal phyla contribute six records each. Within each phylum the
deterministic selector chooses up to four long records (650–1,200 bp) and fills
the remainder with 300–549 bp records, preferring distinct classes/genera,
INSDC accessions, and RefS. Glomeromycota has only two eligible long records
in this source. The short stratum deliberately exposes missing binding sites.
The panel has 48 distinct accessions, 41 with INSDC-style rather than `UDB`
accessions. These are selected reference representatives, not random samples
from fungal diversity.

MAFFT 7.526 adds panel fragments to FN812768.2 with `--keeplength`. A site is
scored only when its entire reference-aligned window is present, has A/C/G/T
bases, and has at least 70% identity to the reference window. The last gate
conservatively rejects fragments that MAFFT may have forced onto short terminal
18S/28S anchors; rejected windows are `homology_unverified_at_site`, **not**
primer failures. A gap, ambiguous base, missing 5′ or 3′ site, or unresolved
single-reference coordinate has its own reason. The 70% alignment gate makes
penalties conditional on verifiable homology and can omit truly divergent
binding sites. More robust locus annotation is needed before estimating clade
wide primer performance.

For a verified site, PrimerMiner 0.22 computes the repository's established
`Position_v1` + `Type_v1`, adjacency 2 penalty. A documented pair is scored
only when both individual sites are scorable and in forward–reverse order; its
penalty is the sum of its two site penalties. No PCR success probability is
implied. The catalog has nine oligos in eight documented pairs. Five pairs have
both coordinates credible on FN812768.2; three involving ITS1F and ITS4B or
LR21 remain unscorable here. ITS1 + ITS4 has credible reference geometry but
no selected sequence with both sites verified.

## Observed denominators

All denominators below are **48 selected records**, with six in each phylum.
Detailed missing reasons and accession-level rows are in the local outputs.

| Primer site or documented pair | Scorable / panel | Main reason for missing sites |
| --- | ---: | --- |
| ITS3 (5.8S) | 48/48 | — |
| fITS7, gITS7, ITS86F (5.8S; each) | 46/48 | 1 alignment gap; 1 unverified homology |
| ITS4 (28S) | 4/48 | 31 missing 3′ site; 11 unverified; 2 gaps |
| ITS1 (18S) | 0/48 | 43 missing 5′ site; 5 unverified |
| ITS1F, ITS4B, LR21 (each) | 0/48 | Reference placement unresolved for 48 |
| ITS3 + ITS4; fITS7 + ITS4; gITS7 + ITS4; ITS86F + ITS4 (each) | 4/48 | Mostly the ITS4 site |
| ITS1 + ITS4; ITS1F + ITS4; ITS1F + ITS4B; ITS1F + LR21 (each) | 0/48 | One or both sites unavailable or unresolved |

| Phylum | Panel | ITS3 site | fITS7 site | ITS4 site | ITS3 + ITS4 pair |
| --- | ---: | ---: | ---: | ---: | ---: |
| Ascomycota | 6 | 6 | 6 | 1 | 1 |
| Basidiomycota | 6 | 6 | 6 | 2 | 2 |
| Chytridiomycota | 6 | 6 | 5 | 0 | 0 |
| Glomeromycota | 6 | 6 | 6 | 0 | 0 |
| Mortierellomycota | 6 | 6 | 6 | 1 | 1 |
| Mucoromycota | 6 | 6 | 6 | 0 | 0 |
| Rozellomycota | 6 | 6 | 6 | 0 | 0 |
| Zoopagomycota | 6 | 6 | 5 | 0 | 0 |

## Reproduce locally

Run from this repository with Python 3, MAFFT, R, and PrimerMiner installed:

```bash
bash scripts/run_fungal_its_pilot.sh
```

The wrapper downloads/checksums the source if needed, selects 48 records,
aligns them, scores sites and pairs, and writes a run manifest with tool
versions and artifact hashes. Raw source files live under
`data/external/its_fungal_pilot/`; detailed manifests, FASTA, accession-level
scores, exact missing-reason summaries, and run metadata live under
`data/derived/its_fungal_pilot/`. Both folders are Git ignored and isolated
from public release artifacts. The site and pair summary CSVs have separate
`n_panel` and `n_scorable` denominators. To expand beyond this pilot, use
independent ITS subregion/anchor annotation and source-specific site checks
before scoring more sequences or interpreting absence as primer failure.
