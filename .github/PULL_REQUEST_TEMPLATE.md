## Purpose

<!-- What problem does this change solve? Link the issue or discussion. -->

## What changed

<!-- Keep this focused on the files, behavior, or evidence changed. -->

## Evidence and provenance

<!-- Cite primary sources and state reference panels, coordinates, denominators, and unresolved ambiguity. Write "not applicable" for code-only changes. -->

## Verification

- [ ] I ran the relevant focused regression test.
- [ ] I ran `for test in tests/test_*.R; do Rscript "$test" || exit 1; done`.
- [ ] I ran `git diff --check`.
- [ ] I checked user-facing behavior in the real app/export state, or it is not applicable.
- [ ] I updated both `PATCHNOTES.md` and `patchnotes.json`, or explained why no patchnote is needed.
- [ ] I did not add large generated data, credentials, or private information to Git.

## Remaining limitations

<!-- What does this PR intentionally not solve? -->
