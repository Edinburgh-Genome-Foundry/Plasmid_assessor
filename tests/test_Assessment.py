from Bio.SeqRecord import SeqRecord
import plasmid_assessor as plasma


def test_assessment():
    sequence = SeqRecord("ATCGATCG")  # 0 site
    design = plasma.Assessment(sequence, "BsmBI")
    design.assess_plasmid()
    assert design.results["number_of_sites"] == 0

    sequence = SeqRecord("AAAAACGTCTCAACTGAAAAAATATCAGAGACGAAAAA")  # 2 sites
    design = plasma.Assessment(sequence, "BsmBI")
    design.assess_plasmid()
    assert design.results["number_of_sites"] == 2
