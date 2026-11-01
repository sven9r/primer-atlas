#!/usr/bin/env python3
"""Offline scientific and append-only acquisition regression fixtures."""
import argparse
import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

spec = importlib.util.spec_from_file_location("coi_builder", "scripts/build_coi_reference_panel.py")
coi = importlib.util.module_from_spec(spec)
spec.loader.exec_module(coi)


def record(accession="AB123456.1", sequence="ACGT" * 175, start=1, end=700,
           gene="COX1", product="cytochrome c oxidase subunit I", complement=False,
           location=None, pseudo=False, duplicate=False, compound=False):
    root = ET.Element("GBSeq")
    for name, value in (("accession-version", accession), ("sequence", sequence.lower()), ("organism", "Fixture organism")):
        ET.SubElement(root, "GBSeq_" + name).text = value
    table = ET.SubElement(root, "GBSeq_feature-table")
    source = ET.SubElement(table, "GBFeature")
    ET.SubElement(source, "GBFeature_key").text = "source"
    quals = ET.SubElement(source, "GBFeature_quals")
    qualifier = ET.SubElement(quals, "GBQualifier")
    ET.SubElement(qualifier, "GBQualifier_name").text = "db_xref"
    ET.SubElement(qualifier, "GBQualifier_value").text = "taxon:1"
    for _ in range(2 if duplicate else 1):
        feature = ET.SubElement(table, "GBFeature")
        ET.SubElement(feature, "GBFeature_key").text = "CDS"
        ET.SubElement(feature, "GBFeature_location").text = location or f"{start}..{end}"
        intervals = ET.SubElement(feature, "GBFeature_intervals")
        for _ in range(2 if compound else 1):
            interval = ET.SubElement(intervals, "GBInterval")
            ET.SubElement(interval, "GBInterval_from").text = str(start)
            ET.SubElement(interval, "GBInterval_to").text = str(end)
            if complement:
                ET.SubElement(interval, "GBInterval_iscomp", value="true")
        qualifiers = ET.SubElement(feature, "GBFeature_quals")
        for key, value in (("gene", gene), ("product", product)) + ((("pseudo", ""),) if pseudo else ()):
            qualifier = ET.SubElement(qualifiers, "GBQualifier")
            ET.SubElement(qualifier, "GBQualifier_name").text = key
            ET.SubElement(qualifier, "GBQualifier_value").text = value
    return root


class FeatureTests(unittest.TestCase):
    def test_whole_mitochondrion_extracts_only_annotated_coi(self):
        whole = "A" * 1000 + "CGTA" * 175 + "T" * 3000
        feature, status = coi.extract_coi(record(sequence=whole, start=1001, end=1700))
        self.assertEqual(status, "verified_marker_feature")
        self.assertEqual(feature["sequence"], "CGTA" * 175)
        self.assertEqual(feature["marker_length_bp"], 700)
        self.assertNotEqual(feature["sequence"], whole)

    def test_reverse_and_partial_marker_preserve_real_bounds(self):
        sequence = "A" * 20 + "CGTA" * 175 + "T" * 20
        feature, _ = coi.extract_coi(record(sequence=sequence, start=21, end=720,
                                           complement=True, location="complement(<21..>720)"))
        self.assertEqual(feature["sequence"], coi.reverse_complement(sequence[20:720]))
        self.assertEqual(feature["feature_strand"], -1)
        self.assertEqual(feature["feature_boundary_state"], "partial")
        self.assertEqual((feature["feature_start"], feature["feature_end"]), (21, 720))

    def test_wrong_locus_ambiguity_pseudogene_and_compound_withheld(self):
        for node, reason in (
            (record(gene="COX2", product="cytochrome c oxidase subunit II"), "missing_coi_feature"),
            (record(duplicate=True), "ambiguous_coi_features"),
            (record(pseudo=True), "annotated_pseudogene"),
            (record(compound=True), "compound_coi_feature_withheld"),
            (record(start=1, end=900), "coi_feature_outside_record"),
            (record(sequence="N" * 30 + "A" * 670), "coi_feature_ambiguity_outside_qc"),
        ):
            with self.subTest(reason=reason):
                feature, status = coi.extract_coi(node)
                self.assertIsNone(feature)
                self.assertEqual(status, reason)

    def test_growth_floor_retains_large_groups_and_never_uses_total_as_cap(self):
        targets = coi.plan_targets(["A", "B", "C"], {"A": 1000, "B": 6, "C": 19},
                                   300, 100, 0.01, monthly=True)
        self.assertEqual(targets, {"A": 1010, "B": 100, "C": 100})
        self.assertGreater(sum(targets.values()), 300)
        targets = coi.plan_targets(["A", "B"], {}, 1001, 100)
        self.assertEqual(sum(targets.values()), 1001)
        self.assertLessEqual(abs(targets["A"] - targets["B"]), 1)

    def test_unknown_taxonomy_ranks_remain_empty(self):
        document = ET.fromstring("<TaxaSet><Taxon><TaxId>1</TaxId><ScientificName>unclassified animal</ScientificName><Rank>no rank</Rank></Taxon></TaxaSet>")
        row = coi.taxonomy_rows(document)["1"]
        self.assertEqual(row["family"], "")
        self.assertEqual(row["species"], "")


class AcquisitionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.groups = self.root / "groups.csv"
        self.groups.write_text("order,environment_group\nAmphipoda,aquatic\n")
        self.documents = {name: record(accession=name) for name in ("AB123456.1", "AB123457.1", "AB123458.1", "AB123459.1")}

    def tearDown(self):
        self.temp.cleanup()

    def args(self, output, **overrides):
        values = dict(output_dir=str(output), groups_file=str(self.groups), previous_dir=None,
                      legacy_raw_dir=None, published_baseline=None, target_total=2,
                      minimum_per_group=1, growth=0.01, release_id="fixture-first",
                      offline=True, refresh=False)
        values.update(overrides)
        return argparse.Namespace(**values)

    def request(self, endpoint, params, output, offline=False, refresh=False):
        key = coi.digest(json.dumps(dict(endpoint=endpoint, params=params), sort_keys=True))
        path = output / "raw" / "content" / (key + (".json" if endpoint == "esearch" else ".xml"))
        path.parent.mkdir(parents=True, exist_ok=True)
        if endpoint == "esearch":
            path.write_text(json.dumps(dict(esearchresult=dict(count="4", idlist=list(self.documents)))))
        elif params["db"] == "taxonomy":
            path.write_text("<TaxaSet><Taxon><TaxId>1</TaxId><ScientificName>Fixture species</ScientificName><Rank>species</Rank><LineageEx><Taxon><TaxId>2</TaxId><ScientificName>Amphipoda</ScientificName><Rank>order</Rank></Taxon></LineageEx></Taxon></TaxaSet>")
        else:
            document = ET.Element("GBSet")
            for accession in params["id"].split(","):
                document.append(self.documents[accession])
            path.write_bytes(ET.tostring(document))
        frozen = path.with_name(coi.file_hash(path) + path.suffix)
        frozen.write_bytes(path.read_bytes())
        return frozen

    def test_two_fresh_workers_keep_exact_old_records_and_append_real_identity(self):
        legacy = self.root / "legacy" / "Amphipoda"
        legacy.mkdir(parents=True)
        original = ">AB123456.1 original preserved header\n" + "ACGT" * 175 + "\n"
        (legacy / "Amphipoda_GB.fasta").write_text(original)
        first, second = self.root / "first", self.root / "second"
        with patch.object(coi, "request", self.request):
            coi.run(self.args(first, legacy_raw_dir=str(legacy.parent)))
            before = coi.read_csv(first / "panel_manifest.csv")
            before_input = coi.file_hash(first / "panel.fasta")
            coi.run(self.args(second, previous_dir=str(first), release_id="fixture-second"))
        after = {row["accession"]: row for row in coi.read_csv(second / "panel_manifest.csv")}
        self.assertGreater(len(after), len(before))
        for old in before:
            self.assertEqual(old["sequence_sha256"], after[old["accession"]]["sequence_sha256"])
            self.assertEqual(old["header_sha256"], after[old["accession"]]["header_sha256"])
            self.assertEqual(old["first_release_id"], after[old["accession"]]["first_release_id"])
        self.assertNotEqual(before_input, coi.file_hash(second / "panel.fasta"))
        self.assertEqual((second / "legacy_raw/Amphipoda/Amphipoda_GB.fasta").read_text(), original)
        qa = coi.read_csv(second / "growth_qa.csv")[0]
        self.assertGreaterEqual(int(qa["retained_after"]), int(qa["target_retained"]))
        self.assertEqual(qa["verified_target_met"], "True")
        self.assertEqual(len(coi.read_fasta(first / "panel.fasta")), len(before))

    def test_changed_previous_sequence_fails_closed(self):
        first = self.root / "first"
        with patch.object(coi, "request", self.request):
            coi.run(self.args(first))
        path = first / "panel.fasta"
        path.write_text(path.read_text().replace("ACGT", "TCGT", 1))
        with self.assertRaisesRegex(RuntimeError, "Previous sequence or header hash differs"):
            coi.import_previous(first, self.root / "second")

    def test_published_recovery_records_known_source_without_claiming_original_bytes(self):
        published = self.root / "published.csv"
        published.write_text("target_group,accession,published_release_id,published_artifact_url,published_artifact_sha256\nAmphipoda,AB123456.1,previous-public,https://example.org/immutable.parquet,known-source-hash\n")
        output = self.root / "published"
        with patch.object(coi, "request", self.request):
            coi.run(self.args(output, published_baseline=str(published)))
        row = {row["accession"]: row for row in coi.read_csv(output / "panel_manifest.csv")}["AB123456.1"]
        self.assertEqual(row["first_release_id"], "previous-public")
        self.assertEqual(row["published_release_id"], "previous-public")
        self.assertEqual(row["baseline_source"], "published release accession recovery")
        self.assertEqual(row["original_raw_bytes_available"], "False")
        self.assertEqual(row["published_artifact_sha256"], "known-source-hash")

    def test_wrong_response_accession_set_fails_before_creating_panel(self):
        original_request = self.request
        def wrong(endpoint, params, output, offline=False, refresh=False):
            path = original_request(endpoint, params, output, offline, refresh)
            if endpoint == "efetch" and params["db"] == "nuccore":
                document = ET.Element("GBSet")
                document.append(record(accession="WR999999.1"))
                path.write_bytes(ET.tostring(document))
            return path
        output = self.root / "wrong"
        with patch.object(coi, "request", wrong), self.assertRaisesRegex(RuntimeError, "different accession-version/UID set"):
            coi.run(self.args(output))
        self.assertFalse((output / "panel.fasta").exists())


class TransportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.output = Path(self.temp.name)
        self.params = dict(db="nuccore", term="Amphipoda", retmode="json")
        self.logs = io.StringIO()
        self.patches = [patch.object(coi.time, "sleep"), patch.object(coi, "_RUNTIME_REPORTED", True),
                        patch.object(coi.sys, "stderr", self.logs)]
        for item in self.patches:
            item.start()

    def tearDown(self):
        for item in reversed(self.patches):
            item.stop()
        self.temp.cleanup()

    def reply(self, code, message, returncode=1):
        return subprocess.CompletedProcess([], returncode,
            json.dumps(dict(ok=False, error=dict(code=code, message=message))), "")

    def success(self, args, **kwargs):
        self.assertIs(kwargs["check"], False)
        config = json.loads(kwargs["input"])
        Path(config["raw_output_path"]).write_text('{"esearchresult":{"count":"2","idlist":["1","2"]}}')
        return subprocess.CompletedProcess(args, 0, json.dumps(dict(ok=True)), "")

    def cache(self):
        key = coi.digest(json.dumps(dict(endpoint="esearch", params=self.params), sort_keys=True))[:24]
        path = self.output / "raw" / "cache" / f"esearch_{key}.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("previous frozen response")
        return path

    def test_nonzero_http_429_is_parsed_redacted_and_retried(self):
        failure = self.reply("network_error", "Request failed: 429 Client Error: Too Many Requests for url: "
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?api_key=PRIVATEKEY "
            "api_key=OTHERSECRET Authorization: Bearer PRIVATEBEARER\n::warning::not-a-log-command")
        calls = []
        def transport(args, **kwargs):
            calls.append(kwargs)
            return failure if len(calls) == 1 else self.success(args, **kwargs)
        with patch.dict(coi.os.environ, {"NCBI_API_KEY": "PRIVATEKEY"}), \
                patch.object(coi.subprocess, "run", transport):
            path = coi.request("esearch", self.params, self.output)
        self.assertEqual(len(calls), 2)
        self.assertEqual(path.stem, coi.file_hash(path))
        log = self.logs.getvalue()
        self.assertIn("network_error; HTTP 429", log)
        self.assertIn("retry in 2s", log)
        for private in ("PRIVATEKEY", "OTHERSECRET", "PRIVATEBEARER", "https://eutils", "\n::warning::"):
            self.assertNotIn(private, log)

    def test_permanent_http_error_fails_fast_preserving_prior_cache(self):
        old = self.cache()
        failure = self.reply("network_error", "Request failed: 403 Client Error: Forbidden for url: "
                             "https://example.org/?api_key=PRIVATEKEY")
        def transport(args, **kwargs):
            Path(json.loads(kwargs["input"])["raw_output_path"]).write_text("failed response")
            return failure
        with patch.object(coi.subprocess, "run", side_effect=transport) as child, \
                self.assertRaisesRegex(RuntimeError, "after 1 attempt.*HTTP 403") as caught:
            coi.request("esearch", self.params, self.output, refresh=True)
        self.assertEqual(child.call_count, 1)
        self.assertEqual(old.read_text(), "previous frozen response")
        self.assertFalse(list(old.parent.glob("*.pending")))
        self.assertNotIn("PRIVATEKEY", str(caught.exception))

    def test_dependency_input_and_tls_failures_are_not_blindly_retried(self):
        for code, message in (("missing_dependency", "requests is required"),
                              ("invalid_input", "Unsupported endpoint"),
                              ("network_error", "SSLError: certificate verify failed (CERTIFICATE_VERIFY_FAILED)")):
            with self.subTest(code=code), patch.object(coi.subprocess, "run", return_value=self.reply(code, message)) as child, \
                    self.assertRaisesRegex(RuntimeError, "after 1 attempt"):
                coi.request("esearch", self.params, self.output)
            self.assertEqual(child.call_count, 1)

    def test_transient_server_failure_has_bounded_six_attempts_and_safe_error(self):
        failure = self.reply("network_error", "Request failed: 503 Server Error: Service Unavailable for url: "
                             "https://example.org/?api_key=PRIVATEKEY")
        with patch.object(coi.subprocess, "run", return_value=failure) as child, \
                self.assertRaisesRegex(RuntimeError, "after 6 attempt.*HTTP 503") as caught:
            coi.request("esearch", self.params, self.output)
        self.assertEqual(child.call_count, 6)
        self.assertNotIn("PRIVATEKEY", str(caught.exception))
        self.assertEqual(coi.time.sleep.call_args_list[-1].args, (32,))

    def test_unstructured_helper_output_is_withheld_and_not_retried(self):
        response = subprocess.CompletedProcess([], 1, "Traceback secret=PRIVATEVALUE https://example.org/private", "PRIVATESTDERR")
        with patch.object(coi.subprocess, "run", return_value=response) as child, \
                self.assertRaisesRegex(RuntimeError, "invalid_helper_reply") as caught:
            coi.request("esearch", self.params, self.output)
        self.assertEqual(child.call_count, 1)
        for secret in ("PRIVATEVALUE", "PRIVATESTDERR", "https://example"):
            self.assertNotIn(secret, str(caught.exception) + self.logs.getvalue())

    def test_helper_timeout_can_retry_and_then_save_real_response(self):
        count = 0
        def transport(args, **kwargs):
            nonlocal count
            count += 1
            if count == 1:
                raise subprocess.TimeoutExpired(args, 90, output="private output")
            return self.success(args, **kwargs)
        with patch.object(coi.subprocess, "run", transport):
            path = coi.request("esearch", self.params, self.output)
        self.assertEqual(count, 2)
        self.assertTrue(path.exists())
        self.assertIn("helper_timeout", self.logs.getvalue())
        self.assertNotIn("private output", self.logs.getvalue())

    def test_claimed_success_without_new_raw_response_cannot_reuse_old_cache(self):
        old = self.cache()
        response = subprocess.CompletedProcess([], 0, json.dumps(dict(ok=True)), "")
        with patch.object(coi.subprocess, "run", return_value=response) as child, \
                self.assertRaisesRegex(RuntimeError, "missing_raw_response"):
            coi.request("esearch", self.params, self.output, refresh=True)
        self.assertEqual(child.call_count, 6)
        self.assertEqual(old.read_text(), "previous frozen response")


if __name__ == "__main__":
    unittest.main()
