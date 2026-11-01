# PrimerMiner scoring runtime

This component contains unchanged copies of `evaluate_primer.R`,
`Position_v1.csv`, and `Type_v1.csv` from PrimerMiner **0.22**, by Vasco Elbrecht,
at upstream commit `9c5e4fb7a4934f590f8df3f6300550a94b04e7f3`.
The upstream package declares **GPL-3**. Its source and these tables retain
that license; the complete text is included in `LICENSE-GPL-3.txt`.

- [Pinned package description](https://github.com/VascoElbrecht/PrimerMiner/blob/9c5e4fb7a4934f590f8df3f6300550a94b04e7f3/PrimerMiner/DESCRIPTION)
- [Pinned evaluator source](https://github.com/VascoElbrecht/PrimerMiner/blob/9c5e4fb7a4934f590f8df3f6300550a94b04e7f3/PrimerMiner/R/evaluate_primer.R)
- [Upstream evaluation documentation](https://github.com/VascoElbrecht/PrimerMiner/wiki/6-Primer-evaluation-(in-silico))

These files are copied from the repository's existing `vendor/PrimerMiner`
source tree, rather than obtained dynamically at application startup.

| File | SHA-256 |
| --- | --- |
| evaluate_primer.R | 2cb0452057f125627927376409e58c6dabafd96fae4279b39b70f119020cecaa |
| Position_v1.csv | 4d2687af16a1680659b372720693dc5fec6a1c994e42f5d827928e39ec21e22a |
| Type_v1.csv | 17f7f1cf85b67603fe25f6afa55fca275aeff5ddf7eb00287c3090f6779833ac |

`R/primerminer_runtime.R` supplies an isolated execution environment and a
native FASTA adapter. It preserves first-token sequence identifiers, reverse
orientation, IUPAC mixture behavior, mismatch-type weights, rounding, adjacent
mismatch multiplication, and unavailable gap/N positions. Both application
custom scorers use this pinned evaluator, independently of installed package
versions. This numerical penalty is an in-silico mismatch score; it does not
establish PCR success.

The default position table covers only distances 1–30 from the primer's 3′ end.
Longer primers retain upstream `NA` behavior; the runtime does not extrapolate
weights. Empty, duplicate-name, or unequal-width alignments fail before scoring.
The original evaluator's single-row reverse-input limitation is retained;
the hosted COI alignments contain a reference plus all retained templates.

`tests/test_primerminer_runtime.R` checks numerical parity with the installed
0.22 package and exercises both custom scorers in a separate process where
PrimerMiner and seqinr are unavailable. Include `R/primerminer_runtime.R` and
**all files in `R/primerminer_runtime/`**, including this notice and the license,
in the Connect bundle. Their combined size is approximately 45 KiB.
