# Atlas provenance note

This directory vendors the `PrimerMiner` R package at upstream commit
`9c5e4fb7a4934f590f8df3f6300550a94b04e7f3` from
<https://github.com/VascoElbrecht/PrimerMiner> under its GPL-3 license.

The only Atlas compatibility change is removal of `BOLDconnectR` from
`Depends`. The Atlas uses PrimerMiner's NCBI, alignment, clustering, and primer
evaluation components; it does not invoke the BOLD download functions.
