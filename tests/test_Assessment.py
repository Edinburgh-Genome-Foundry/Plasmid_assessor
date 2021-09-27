from Bio.SeqRecord import SeqRecord
import plasmid_assessor as plasma


def test_assess_plasmid():
    sequence = SeqRecord("ATCGATCG")  # 0 site
    design = plasma.Assessment(sequence, "BsmBI")
    design.assess_plasmid()


def test_get_number_of_sites():
    sequence = SeqRecord("ATCGATCG")  # 0 site
    design = plasma.Assessment(sequence, "BsmBI")
    design.get_number_of_sites()
    assert design.results["number_of_sites"] == 0

    sequence = SeqRecord("AAAAACGTCTCAACTGAAAAAATATCAGAGACGAAAAA")  # 2 sites
    design = plasma.Assessment(sequence, "BsmBI")
    design.get_number_of_sites()
    assert design.results["number_of_sites"] == 2
