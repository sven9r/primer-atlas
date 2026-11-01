The Entrez client is the unmodified `scripts/ncbi_entrez.py` helper supplied by
the installed life-science-research NCBI Entrez skill, version 1.0.3. It makes
targeted calls to NCBI's public E-utilities API. The fungal panel builder uses
this helper for acquisition and saves the original responses and their hashes.

Runtime dependency: Python 3 and `requests`. API credentials are optional and
are read from environment variables; they are never stored in panel manifests.
