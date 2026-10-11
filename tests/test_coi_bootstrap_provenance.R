#!/usr/bin/env Rscript
.libPaths(c(file.path(getwd(), ".Rlib"), .libPaths()))
source("R/data_layer.R")
source("R/coi_retention.R")

# Exercise the real R baseline writer and Python acquisition parser together.
# NCBI responses are bounded offline fixtures; the immutable source descriptor
# and its checksum are real, rather than a separately hand-written baseline.
fixture <- tempfile("coi_bootstrap_provenance_")
dir.create(fixture)
on.exit(unlink(fixture, recursive = TRUE), add = TRUE)
dir.create(file.path(fixture, "data"))
write.csv(data.frame(order = c("Amphipoda", "Bivalvia"), environment_group = "test"),
  file.path(fixture, "data/orders_22.csv"), row.names = FALSE)
membership <- data.frame(order = c("Amphipoda", "Amphipoda", "Bivalvia"),
  template = c("AB123456.1 original published header", "AB123457.1", "BV123456.1"))
artifact_path <- file.path(fixture, "published-membership.parquet")
nanoparquet::write_parquet(membership, artifact_path)
artifact_url <- "https://data.example.invalid/releases/COI/previous-public/full_order_anml.parquet"
artifact_hash <- atlas_sha256(artifact_path)
manifest <- list(marker_id = "COI", release_version = "previous-public",
  artifacts = list(full_order_anml = list(local_path = artifact_path, url = artifact_url, sha256 = artifact_hash)))
baseline <- coi_recover_published_membership(manifest, fixture)
stopifnot(nrow(baseline) == 2L, all(baseline$target_group == "Amphipoda"),
  all(baseline$published_artifact_url == artifact_url), all(baseline$published_artifact_sha256 == artifact_hash),
  all(baseline$published_release_id == "previous-public"), all(!baseline$original_raw_bytes_available))
stopifnot(!"publication_artifact_sha256" %in% names(baseline))

python_code <- r"(
import argparse, importlib.util, json, sys
from pathlib import Path
import xml.etree.ElementTree as ET
spec = importlib.util.spec_from_file_location('coi_test_fixtures', 'tests/test_coi_reference_panel.py')
fixtures = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fixtures)
coi = fixtures.coi
root = Path(sys.argv[1])
output = root / 'acquired'
cache = output / 'raw/cache'
cache.mkdir(parents=True)
def response(endpoint, params, content):
    key = coi.digest(json.dumps(dict(endpoint=endpoint, params=params), sort_keys=True))[:24]
    extension = 'json' if params.get('retmode') == 'json' else 'xml'
    (cache / (endpoint + '_' + key + '.' + extension)).write_text(content)
def sequence_response(accessions):
    document = ET.Element('GBSet')
    for accession in accessions:
        document.append(fixtures.record(accession=accession))
    return ET.tostring(document, encoding='unicode')
response('efetch', dict(db='nuccore', id='AB123456.1,AB123457.1', rettype='gb', retmode='xml'),
         sequence_response(['AB123456.1','AB123457.1']))
response('efetch', dict(db='taxonomy', id='1', retmode='xml'),
         '<TaxaSet><Taxon><TaxId>1</TaxId><ScientificName>Fixture species</ScientificName><Rank>species</Rank><LineageEx><Taxon><TaxId>2</TaxId><ScientificName>Amphipoda</ScientificName><Rank>order</Rank></Taxon></LineageEx></Taxon></TaxaSet>')
response('esearch', dict(db='nuccore', term='Amphipoda' + coi.QUERY_SUFFIX, retmode='json', retmax=2000, idtype='acc'),
         json.dumps(dict(esearchresult=dict(count='3', idlist=['AB123456.1','AB123457.1','AB123458.1']))))
response('efetch', dict(db='nuccore', id='AB123458.1', rettype='gb', retmode='xml'),
         sequence_response(['AB123458.1']))
coi.run(argparse.Namespace(output_dir=str(output), groups_file=str(root / 'data/orders_22.csv'),
    previous_dir=None, legacy_raw_dir=None, published_baseline=str(root / 'data/provenance/coi_published_baseline.csv'),
    target_total=2, minimum_per_group=1, growth=0.01, release_id='new-acquisition', offline=True, refresh=False))
rows = {row['accession']:row for row in coi.read_csv(output / 'panel_manifest.csv')}
assert len(rows) == 3
for accession in ('AB123456.1','AB123457.1'):
    row = rows[accession]
    assert row['first_release_id'] == row['published_release_id'] == 'previous-public'
    assert row['published_artifact_url'] == sys.argv[2]
    assert row['published_artifact_sha256'] == sys.argv[3]
    assert row['original_raw_bytes_available'] == 'False'
    assert row['baseline_source'] == 'published release accession recovery'
    assert row['source_sha256'] != row['published_artifact_sha256']
assert rows['AB123458.1']['published_release_id'] == ''
assert rows['AB123458.1']['baseline_source'] == 'new NCBI annotated COI feature'
print('R bootstrap -> Python recovery preserved prior publication URL/hash/release and unavailable-original-bytes caveat.')
)"
script <- file.path(fixture, "verify-bootstrap.py")
writeLines(python_code, script)
python <- Sys.getenv("COI_PYTHON", Sys.which("python3"))
stopifnot(nzchar(python))
status <- system2(python, vapply(c(script, fixture, artifact_url, artifact_hash), shQuote, character(1)))
stopifnot(status == 0L)
