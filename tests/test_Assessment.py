from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import plasmid_assessor as plasma


def test_assess_plasmid():
    record = SeqRecord(Seq("ATCGATCG"))
    design = plasma.Assessment(record, "BsmBI")
    design.assess_plasmid()


def test_check_circularity():
    sequence = SeqRecord(Seq("ATCGATCG"), annotations={"topology": "circular"})
    design = plasma.Assessment(sequence, "BsmBI")
    design.check_circularity()
    assert design.results["is_circular"]


def test_get_number_of_sites():
    sequence = SeqRecord(Seq("ATCGATCG"))  # 0 site
    design = plasma.Assessment(sequence, "BsmBI")
    design.get_number_of_sites()
    assert design.results["number_of_sites"] == 0

    sequence = SeqRecord(Seq("AAAAACGTCTCAACTGAAAAAATATCAGAGACGAAAAA"))  # 2 sites
    design = plasma.Assessment(sequence, "BsmBI")
    design.get_number_of_sites()
    assert design.results["number_of_sites"] == 2


def test_evaluate_orientation():
    sequence = SeqRecord(Seq("AAAA" + "CGTCTCAACTG" + "AAAAA" + "TATCAGAGACG" + "AAAA"))
    design = plasma.Assessment(sequence, "BsmBI")
    design.evaluate_orientation()
    assert design.results["is_site_orientation_correct"]
