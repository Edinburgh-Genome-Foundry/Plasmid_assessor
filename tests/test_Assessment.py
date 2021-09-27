from Bio.SeqRecord import SeqRecord
import plasmid_assessor as plasma


def test_assessment():
    sequence = SeqRecord("ATCGATCG")
    plasma.Assessment(sequence, "BsmBI")
