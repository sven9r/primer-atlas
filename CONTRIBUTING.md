# Contributing to Primer Atlas

Primer Atlas welcomes code, documentation, primer records, reference-panel
improvements, tests, and careful scientific review. The project treats every
primer claim and computed score as traceable evidence, so contributions must
preserve sources, denominators, and unresolved ambiguity.

## Start with the right conversation

- Use the **primer proposal** form for a cited primer pair.
- Use the **reference-panel proposal** form for new or revised sequence evidence.
- Use the **bug report** form for reproducible application or pipeline failures.
- Use GitHub Discussions for early ideas, interpretation questions, and work
  that is not yet ready to become a scoped issue.

For a substantial change, open or claim an issue before implementation. This
prevents two contributors from solving the same problem differently and gives
maintainers a chance to confirm the scientific evidence contract.

## Local setup

Clone the repository and restore the pinned R environment:

```bash
git clone git@github.com:sven9r/primer-atlas.git
cd primer-atlas
Rscript -e 'renv::restore()'
Rscript scripts/install_primerminer.R
```

System tools used by the full reference build include MAFFT and VSEARCH. The
Shiny application can then be started with:

```bash
Rscript scripts/run_app.R
```

## Contribution workflow

1. Create a branch from the latest `main`.
2. Make one logically scoped change.
3. Add or update a regression test that would fail without the change.
4. Run the relevant focused test and, before review, the complete test loop.
5. Update `PATCHNOTES.md` and `patchnotes.json` with the same entry ID when the
   change affects users, data, releases, or project operations.
6. Open a pull request and complete its evidence checklist.
7. Address review comments without force-pushing over another contributor's
   work.

The complete local regression loop is:

```bash
for test in tests/test_*.R; do Rscript "$test" || exit 1; done
git diff --check
```

GitHub Actions must pass before a pull request can merge.

## Scientific evidence requirements

Primer proposals must include the marker, both oligos in synthesized 5-prime to
3-prime orientation, claimed applications, target taxa, expected product, and
primary sources. Submissions are proposals, not automatic catalog additions.

Curated data changes must:

- preserve original claim wording and primary citations;
- distinguish published evidence from Atlas-derived calculations;
- identify the coordinate reference and reference-panel scope;
- retain exact sequence, accession, sample, and denominator fields where
  available;
- represent unknown or disputed information explicitly instead of guessing;
- avoid describing in-silico mismatch scores as amplification probabilities.

Large sequences and generated score tables belong in the configured object
store, not Git. Small catalogs, schemas, manifests, provenance records, and
test fixtures belong in the repository.

## Pull-request review

A pull request is ready when its purpose is clear, the changed evidence can be
traced to its source, tests pass, and user-facing behavior has been checked at
the real interaction or export state. Maintainers may request changes when a
proposal is scientifically plausible but its provenance, denominator, or
reference scope is incomplete.

## Community standards

Participation is governed by `CODE_OF_CONDUCT.md`. Be specific and constructive
when challenging a method or interpretation: disagreement about evidence is
welcome; personal attacks and harassment are not.
