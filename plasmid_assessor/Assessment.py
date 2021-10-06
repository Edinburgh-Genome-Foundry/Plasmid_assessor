import re

import matplotlib.pyplot as plt

import Bio
import Bio.Restriction
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation

import dnacauldron as dc

try:
    from dna_features_viewer import BiopythonTranslator
except ImportError:

    class AssessmentTranslator:
        """Please install dna_features_viewer to use this class."""

        def __init__(self):
            raise Exception("Please install dna_features_viewer to use this class.")


else:

    class AssessmentTranslator(BiopythonTranslator):
        """Custom translator for highlighting key features."""

        def compute_feature_color(self, feature):
            assessment_ref = "plasmid_assessment"
            if assessment_ref in feature.qualifiers:
                if feature.qualifiers[assessment_ref] == "enzyme":
                    return "red"
                elif feature.qualifiers[assessment_ref] == "excised":
                    return "yellow"
                elif feature.qualifiers[assessment_ref] == "backbone":
                    return "tab:cyan"
                else:
                    return "tab:blue"  # default dna_features_viewer color
            else:
                return "tab:blue"


class Assessment:
    """The plasmid assessment class.


    **Parameters**

    **record**
    > A Biopython `SeqRecord`.

    **enzyme**
    > A restriction enzyme (`str`). A Biopython `RestrictionType` will be looked
    up using the string.
    """

    UNKNOWN_IDS = [
        "None",
        "",
        "<unknown id>",
        ".",
        "EXPORTED",
        "<unknown name>",
        "Exported",
    ]

    def __init__(self, record, enzyme):
        self.record = record
        self.enzyme = Bio.Restriction.__dict__[enzyme]
        self.enzyme_name = str(self.enzyme)
        self.results = {}

    def assess_plasmid(self, other_enzymes=None):
        """Evaluate plasmid for Golden Gate.


        **Parameters**

        **other_enzymes**
        > List of enzymes used in higher level assemblies (`list`).
        """
        if other_enzymes:
            self.other_enzymes = ", ".join([str(enz) for enz in other_enzymes])
        self.add_name()
        self.check_circularity()
        self.get_number_of_sites()
        self.evaluate_orientation()
        self.digest_plasmid()
        self.count_other_sites(other_enzymes)
        self.check_enzyme_site_locations()
        self.sum_results()
        self.plot_plasmid()

    def add_name(self):
        """Set a name for the assessment."""
        # To display on the report:
        if str(self.record.id).strip() in self.UNKNOWN_IDS:
            self.name = "Unnamed plasmid"
        else:
            if len(self.record.id) > 16:  # Genbank limit, also for width in report
                self.name = self.record.id[:16] + "..."
            else:
                self.name = self.record.id

    def check_circularity(self):
        if "topology" not in self.record.annotations:
            self.results["is_circular"] = False
        elif self.record.annotations["topology"] == "circular":
            self.results["is_circular"] = True
        else:
            self.results["is_circular"] = False

    def get_number_of_sites(self):
        if "is_circular" in self.results:
            is_linear = not self.results["is_circular"]
        else:
            is_linear = False
        restriction_batch = Bio.Restriction.RestrictionBatch([self.enzyme])
        analysis = Bio.Restriction.Analysis(
            restriction_batch, sequence=self.record.seq, linear=is_linear
        )
        self.analysis_results = analysis.full(linear=is_linear)

        self.results["number_of_sites"] = len(self.analysis_results[self.enzyme])

        # Add as features for plot in report:
        for enzyme, sites in self.analysis_results.items():
            for site in sites:
                self.record.features.append(
                    SeqFeature(
                        FeatureLocation(site, site + 1),
                        id=str(enzyme),
                        type="misc_feature",
                        qualifiers={
                            "label": str(enzyme),
                            "plasmid_assessment": "enzyme",
                        },
                    )
                )

    def evaluate_orientation(self):
        self.results["is_site_orientation_correct"] = False  # default
        # Forward strand:
        self.iter_forward = [
            match.end() for match in re.finditer(self.enzyme.site, str(self.record.seq))
        ]
        if sum(1 for _ in self.iter_forward) == 1:
            self.forward_enzyme = self.iter_forward[0]
            # rev_complement_site = str(self.record.seq.reverse_complement())
            rev_complement_site = str(Seq(self.enzyme.site).reverse_complement())
            self.iter_reverse = [
                m.start()
                for m in re.finditer(rev_complement_site, str(self.record.seq))
            ]
            if sum(1 for _ in self.iter_reverse) == 1:  # 1 site in both strands:
                self.results["is_site_orientation_correct"] = True
                self.reverse_enzyme = self.iter_reverse[0]

        if self.results["is_site_orientation_correct"]:
            if self.reverse_enzyme < self.forward_enzyme:
                self.record.features.append(
                    SeqFeature(
                        FeatureLocation(
                            self.reverse_enzyme - 1, self.forward_enzyme + 1
                        ),
                        id=str(self.enzyme),
                        type="misc_feature",
                        qualifiers={
                            "label": "Excised",
                            "plasmid_assessment": "excised",
                        },
                    )
                )
            else:  # put annotation together from two pieces:
                self.record.features.append(
                    SeqFeature(
                        FeatureLocation(0, self.forward_enzyme + 1),
                        id=str(self.enzyme),
                        type="misc_feature",
                        qualifiers={
                            "label": "Excised",
                            "plasmid_assessment": "excised",
                        },
                    )
                )
                self.record.features.append(
                    SeqFeature(
                        FeatureLocation(self.reverse_enzyme - 1, len(self.record)),
                        id=str(self.enzyme),
                        type="misc_feature",
                        qualifiers={
                            "label": "Excised",
                            "plasmid_assessment": "excised",
                        },
                    )
                )

    def digest_plasmid(self):
        # Obtain fragments and get the backbone's overhangs.
        # This method has two assumptions:
        # - the sequence has two, correctly oriented enzyme sites.
        # - the sequence is circular.
        # Therefore there will be exactly two fragments, with one containing both sites.
        self.results["digest"] = {}
        if not self.results["is_circular"]:
            return
        if not self.results["is_site_orientation_correct"]:
            return

        record_fragments = dc.StickyEndFragment.list_from_record_digestion(
            record=self.record, enzyme=self.enzyme, linear=False
        )
        if self.enzyme.site in record_fragments[0].to_standard_string():
            backbone_index = 1  # there are only two fragments
            excise_index = 0
        else:
            backbone_index = 0
            excise_index = 1  # reversed
        self.results["digest"]["backbone_seq"] = record_fragments[backbone_index]
        self.results["digest"]["excised_seq"] = record_fragments[excise_index]
        self.results["digest"]["first_overhang"] = str(
            record_fragments[excise_index].seq.left_end
        )
        self.results["digest"]["last_overhang"] = str(
            record_fragments[excise_index].seq.right_end
        )

    def count_other_sites(self, other_enzymes):
        self.results["other_sites"] = {}
        self.results["other_sites"]["has_any_other_sites"] = False

        if other_enzymes is None:
            return
        bio_enzymes = [Bio.Restriction.__dict__[enzyme] for enzyme in other_enzymes]

        restriction_batch = Bio.Restriction.RestrictionBatch(bio_enzymes)
        # Work with the assumption that the sequence is circular:
        analysis = Bio.Restriction.Analysis(
            restriction_batch, sequence=self.record.seq, linear=False
        )
        self.results["other_sites"]["enzyme"] = analysis.full(linear=False)

        for enzyme, matches in self.results["other_sites"]["enzyme"].items():
            if len(matches) != 0:
                self.results["other_sites"]["has_any_other_sites"] = True
                # Also add as features for plot in report:
                for site in matches:
                    self.record.features.append(
                        SeqFeature(
                            FeatureLocation(site, site + 1),
                            id=str(enzyme),
                            type="misc_feature",
                            qualifiers={
                                "label": str(enzyme),
                                "plasmid_assessment": "enzyme",
                            },
                        )
                    )

    def check_enzyme_site_locations(self):
        """Flag enzyme sites that are within the retained backbone."""
        try:
            self.results["other_sites"]["has_any_other_sites"]
            self.results["is_site_orientation_correct"]
        except KeyError:
            print("Run assessment methods first!")
        else:
            self.sites_outside_excised_region = {}
            if (
                self.results["other_sites"]["has_any_other_sites"]
                and self.results["is_site_orientation_correct"]
            ):
                # if there are no other sites, no need to run:
                if self.reverse_enzyme < self.forward_enzyme:
                    # orientation = reverse -> forward
                    for enzyme, sites in self.results["other_sites"]["enzyme"].items():
                        problem_sites = []
                        for site in sites:
                            if self.reverse_enzyme < site < self.forward_enzyme:
                                pass
                            else:
                                problem_sites += [str(site)]
                        if problem_sites != []:
                            self.sites_outside_excised_region[
                                str(enzyme)
                            ] = problem_sites
                    txt = ""  # for the pdf report
                    for (
                        enzyme,
                        problem_sites,
                    ) in self.sites_outside_excised_region.items():
                        txt += enzyme + ": " + " ".join(problem_sites) + ";"
                    self.sites_outside_excised_region_txt = txt
                else:
                    # orientation = forward -> reverse
                    for enzyme, sites in self.results["other_sites"]["enzyme"].items():
                        problem_sites = []
                        for site in sites:
                            if self.forward_enzyme < site < self.reverse_enzyme:
                                # in this case the site is within the retained backbone
                                problem_sites += [str(site)]
                        if problem_sites != []:
                            self.sites_outside_excised_region[
                                str(enzyme)
                            ] = problem_sites
                    txt = ""  # for the pdf report
                    for (
                        enzyme,
                        problem_sites,
                    ) in self.sites_outside_excised_region.items():
                        txt += enzyme + ": " + " ".join(problem_sites) + ";"
                    self.sites_outside_excised_region_txt = txt

            else:  # no other sites or orientation not correct
                self.sites_outside_excised_region_txt = ""

    def sum_results(self):
        self.results["pass"] = True
        if self.results["is_circular"] is False:
            self.results["pass"] = False
            return
        if self.results["is_site_orientation_correct"] is False:
            # implicitly checks number of sites too
            self.results["pass"] = False
            return
        if self.sites_outside_excised_region_txt:
            self.results["pass"] = False
            return
        # if self.results["other_sites"]["has_any_other_sites"]:
        #     self.results["pass"] = False
        #     return

    def plot_plasmid(self):
        """Plot an outline of the plasmid."""

        fig, ax = plt.subplots(figsize=(7, 4))
        graphic_record = AssessmentTranslator().translate_record(self.record)
        graphic_record.plot(ax=ax, with_ruler=False, strand_in_label_threshold=2)

        self.fig = fig
