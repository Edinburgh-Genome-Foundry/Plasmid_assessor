import Bio
import Bio.Restriction


class Assessment:
    """The plasmid assessment class.


    **Parameters**

    **record**
    > A Biopython `SeqRecord`.

    **enzyme**
    > A restriction enzyme (`str`). A Biopython `RestrictionType` will be looked
    up using the string.
    """

    def __init__(self, record, enzyme):
        self.record = record
        self.enzyme = Bio.Restriction.__dict__[enzyme]
        self.results = {}

    def assess_plasmid(self, other_enzymes=None):
        """Evaluate plasmid for Golden Gate


        **Parameters**

        **other_enzymes**
        > List of enzymes used in higher level assemblies (`list`).
        """
        self.check_circularity()
        self.get_number_of_sites()

    def check_circularity(self):
        if "topology" not in self.record.annotations:
            self.results["is_circular"] = False
        elif self.record.annotations["topology"] == "circular":
            self.results["is_circular"] = True
        else:
            self.results["is_circular"] = False

    def get_number_of_sites(self):
        restriction_batch = Bio.Restriction.RestrictionBatch([self.enzyme])
        analysis = Bio.Restriction.Analysis(
            restriction_batch, sequence=Bio.Seq.Seq(self.record.seq)
        )
        analysis_results = analysis.full()

        self.results["number_of_sites"] = len(analysis_results[self.enzyme])
