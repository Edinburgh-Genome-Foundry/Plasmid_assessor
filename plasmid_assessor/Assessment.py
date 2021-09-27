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

    properties = [
        "is_circular"
        "number_of_sites"
        "is_site_orientation_correct"
        "left_overhang"
        "right_overhang"
        "insert_seq"
        "backbone_seq"
    ]  # assessment is performed in this order

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
        pass
