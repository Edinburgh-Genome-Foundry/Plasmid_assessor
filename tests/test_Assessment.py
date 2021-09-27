from Bio.SeqRecord import SeqRecord
import plasmid_assessor as plasma


def test_assess_plasmid():
    sequence = SeqRecord("ATCGATCG")
    design = plasma.Assessment(sequence, "BsmBI")
    design.assess_plasmid()


def test_check_circularity():
    sequence = SeqRecord("ATCGATCG", annotations={"topology": "circular"})
    design = plasma.Assessment(sequence, "BsmBI")
    design.check_circularity()
    assert design.results["is_circular"]


def test_get_number_of_sites():
    sequence = SeqRecord("ATCGATCG")  # 0 site
    design = plasma.Assessment(sequence, "BsmBI")
    design.get_number_of_sites()
    assert design.results["number_of_sites"] == 0

    sequence = SeqRecord("AAAAACGTCTCAACTGAAAAAATATCAGAGACGAAAAA")  # 2 sites
    design = plasma.Assessment(sequence, "BsmBI")
    design.get_number_of_sites()
    assert design.results["number_of_sites"] == 2


def test_evaluate_orientation():
    sequence = SeqRecord("AAAAA" + "CGTCTCAACTG" + "AAAAA" + "TATCAGAGACG" + "AAAAA")
    design = plasma.Assessment(sequence, "BsmBI")
    design.evaluate_orientation()
    assert design.results["is_site_orientation_correct"]
