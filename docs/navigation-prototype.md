> Status update, 2026-10-02: the reviewed navigation and map repairs are included
> in the 2026.10.02 release candidate. The local walkthroughs below are historical.
> Current verification and publication state: [release note](2026-10-02-release.md).
> Organism filtering now excludes FishF2, MollCOI253, and Vertebrate COI from
> Arthropods; choose All organism groups to inspect their off-target details.

# Local navigation prototype — 2026-09-23

This is a reversible `app.R` navigation change. It does not alter scoring,
reference data, filter logic, marker releases, or deployment files.

## Walkthrough

1. **Start here** opens first. Four task cards lead to primer discovery, COI
   reference fit by taxon, COI reference-region comparison, and taxon/sequence
   evidence. The evidence card also links to sources and methods.
2. **Find primers** retains the organism, marker, facet, primer-pair, and custom
   pair controls. It is the entry point for the active COI view and the fungal
   ITS and 18S pilot catalogs and maps.
3. **Check taxonomic fit · COI** retains the order heatmap and its inspector.
   Its button opens **Inspect taxa & sequences** with the current pair and
   target group. A missing full-panel artifact still yields no lineage rows;
   the navigation does not imply that those rows exist.
4. **Compare regions · COI** retains its country/territory controls, located/all
   denominators, within-order table, binding-site differences, exact sequence
   tab, and method limits. It compares accession source geography, not species
   range or regional PCR success.
5. **Check temperatures** and **Sources & methods** remain accessible as
   separate tabs. The former is a gradient-PCR starting point; the latter
   retains citations and the fungal primer catalog.

The readiness panel calls COI active, fungal ITS and 18S pilots, and 12S/16S/28S
planned. It distinguishes the app's single-reference ITS map from the later
local 48-sequence scoring pilot. The separate Hawaiʻi–Madagascar audit and range
catalog have not been connected to the app. Reference fit, site availability,
and geography are not amplification probabilities.

## Verification

- `Rscript -e 'parse(file="app.R")'` passed.
- `tests/test_marker_runtime.R`, `tests/test_its_primer_navigator.R`,
  `tests/test_18s_catalog_ui.R`, `tests/test_region_comparison.R`,
  `tests/test_order_lens_interaction.R`, and
  `tests/test_catalog_and_geography.R` passed. `git diff --check` passed.
- A fresh local Shiny session at `127.0.0.1:4387` showed the Start here panel.
  The task buttons opened the primer map, COI fit, region comparison, and
  evidence views; the sources button opened the citation registry. The fit
  inspector's existing drill-down opened the taxon view with its selected pair.
  The primer map populated nine default COI rows after reactive loading. In
  the fungal ITS pilot, eight documented combinations and the map loaded, and
  the full catalog tab showed 20 of 121 primer records on its first page.

## Review decision

Choose whether the task labels and landing panel should replace the current
public navigation. In particular, check whether the distinction between
"reference fit by taxon" and biological amplification coverage is clear to a
first-time visitor. No deployment is part of this prototype.

## 2026-09-24 local map repair

The map appeared stuck in the narrow in-app browser. A Region comparison
observer was also updating its own `region_order` input repeatedly between
blank and `All`, consuming a full CPU core even while that tab was hidden.
The observer now refreshes the order choices only when the region primer pair
changes. At widths up to 850 px, the Find primers sidebar and map stack, the
pair list gets its own scroll area, and the metric cards use two columns. A
Jump to primer map link and horizontal-scroll hint make the SVG reachable.
The custom-pair form starts collapsed; opening it still exposes the same inputs.

In a fresh local Shiny session, the COI map reached nine default rows; removing
BF2 + BR2 reduced it to eight, and the SVG scrolled horizontally. The server
returned to idle. `test_primer_map_domain.R`, `test_its_primer_navigator.R`,
`test_region_comparison.R`, `test_marker_runtime.R`, and `git diff --check`
passed. The repair remains local and unpublished.

## 2026-09-24 map detail and usage labels

In **Find primers**, choose **Arthropods** and inspect the COI pair list.
FishF2 + FishR1, MollCOI253, and Vertebrate COI remain available for comparison,
but their labels now state the intended target and that the arthropod reference
is an off-target check only. This follows each pair's catalog claim; a shared
COI marker does not turn a fish, mollusk, or vertebrate pair into an arthropod
coverage result.

Check FishF2 + FishR1, then click its numbered map row or a primer segment.
The pinned detail now shows each oligo's estimated Tm at the selected salt,
the suggested pair Ta starting range, the source-reported Ta (52 °C for this
pair), and the off-target limitation. The same detail is used on hover and
keyboard focus. The estimate is a simple thermal model and is distinguished
from the reported value. ITS and 18S map rows show the same fields; 18S says
"not reported" when no source Ta exists.

`tests/test_primer_map_details.R` checks the labels, map detail, and organism
switch. In a fresh local Shiny session, the map retained its rows when
switching to Arthropods; the FishF2 row displayed the temperatures and
limitation when clicked. ITS and 18S map details also rendered. The existing
pair selection, application/environment/intent, and amplicon filters remain.
This local prototype is unpublished.
