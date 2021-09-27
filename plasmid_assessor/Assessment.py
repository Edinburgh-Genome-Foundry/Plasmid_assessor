import re

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
        if "is_circular" in self.results:
            is_linear = not self.results["is_circular"]
        else:
            is_linear = False
        analysis_results = analysis.full(linear=is_linear)

        self.results["number_of_sites"] = len(analysis_results[self.enzyme])

    def evaluate_orientation(self):
        self.results["is_site_orientation_correct"] = False  # default
        # Forward strand:
        iter = (match for match in re.finditer(self.enzyme.site, self.record.seq))
        if sum(1 for _ in iter) == 1:
            rev_complement = Bio.Seq.reverse_complement(self.record.seq)
            iter_reverse = (m for m in re.finditer(self.enzyme.site, rev_complement))
            if sum(1 for _ in iter_reverse) == 1:  # 1 site in both strands:
                self.results["is_site_orientation_correct"] = True
