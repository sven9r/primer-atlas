import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("unite_sites", ROOT / "scripts/scan_unite_primer_sites.py")
engine = importlib.util.module_from_spec(spec)
spec.loader.exec_module(engine)

class ExactSiteTests(unittest.TestCase):
    def test_degeneracy_ambiguity_overlap_and_specificity(self):
        query = engine.prepare_primer("GTGARTCATCGARTCTTTG", "forward")
        seq = "AAC" + "GTGAATCATCGAATCTTTG" + "GG"
        self.assertEqual(engine.exact_site(seq, query)["start"], 4)
        self.assertEqual(engine.exact_site(seq.replace("GTGAAT", "GTGNAT"), query)["n_matches"], 0)
        repeated = engine.exact_site("AAAAAAAAAA", engine.prepare_primer("AAAAAAA", "forward"))
        self.assertEqual(repeated["n_matches"], 4)
        self.assertEqual(repeated["status"], "ambiguous_binding_site")
        low = engine.exact_site("ACGTACGT", engine.prepare_primer("NNNNNNNN", "forward"))
        self.assertEqual(low["status"], "low_specificity_site")
        self.assertEqual(engine.exact_site("AC", query)["status"], "template_shorter_than_oligo")
        self.assertEqual(engine.exact_site(seq, engine.prepare_primer("AC?", "forward"))["status"], "invalid_oligo")

    def test_reverse_site_and_encoded_product(self):
        fseq, rseq = "TCCGTAGGTGAACCTGCGG", "TCCTCCGCTTATTGATATGC"
        seq = "AA" + fseq + "N"*105 + engine.reverse_complement(rseq) + "CC"
        f = engine.exact_site(seq, engine.prepare_primer(fseq,"forward"))
        r = engine.exact_site(seq, engine.prepare_primer(rseq,"reverse"))
        status,bp=engine.pair_result(f,r)
        self.assertEqual(status,"scorable")
        self.assertEqual(bp,len(fseq)+105+len(rseq))
        self.assertEqual(engine.pair_result(r,f)[0],"outward_or_overlapping_sites")
        no=engine.exact_site("A"*40,engine.prepare_primer(fseq,"forward"))
        self.assertEqual(engine.pair_result(f,no)[0],"reverse_site_missing")

    def test_frozen_r_engine_replay(self):
        # Independent frozen R outputs protect engine parity on real primers,
        # ambiguity, missing sites, geometry and per-template product lengths.
        source=ROOT/"data/pinned/ITS_FUNGAL/versions/v3-c6fa87f256e1"
        catalog=engine.read_csv(source/"catalog.csv")
        seqs=dict(engine.read_fasta(source/"sequences.fasta"))
        lookup={(r["accession"],r["catalog_id"]):r for r in engine.read_csv(source/"sites.csv")}
        computed={}
        for name,sequence in seqs.items():
            for row in catalog:
                result=engine.exact_site(sequence,engine.prepare_primer(row["sequence"],row["direction"]))
                old=lookup[(name,row["catalog_id"])]
                self.assertEqual(result["status"],old["status"],(name,row["primer_name"]))
                self.assertEqual(result["n_matches"],int(old["n_matches"]))
                for field in ("start","end"):
                    expected=int(old[field]) if old[field] not in ("","NA") else None
                    self.assertEqual(result[field],expected)
                computed[(name,row["primer_name"],row["sequence"],row["direction"])]=result
        pairs={p["pair_id"]:p for p in engine.read_csv(source/"pairs.csv")}
        for old in engine.read_csv(source/"pair_sequences.csv"):
            p=pairs[old["pair_id"]]; name=old["accession"]
            f=computed[(name,p["forward_primer"],p["forward_sequence"],"forward")]
            r=computed[(name,p["reverse_primer"],p["reverse_sequence"],"reverse")]
            state,bp=engine.pair_result(f,r)
            self.assertEqual(state,old["status"])
            self.assertEqual(bp,int(old["amplicon_bp"]) if old["amplicon_bp"] not in ("","NA") else None)

if __name__ == "__main__": unittest.main()
