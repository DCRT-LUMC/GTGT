from typing import Any, Sequence

import pytest
from mutalyzer.description import Description
from mutalyzer.description_model import get_reference_id
from mutalyzer.protein import get_protein_description
from mutalyzer.reference import get_protein_selector_model

from gtgt.mutalyzer import (
    Variant,
    combine_variants_deletion,
    init_description,
    mutation_to_cds_effect,
    sequence_from_description,
    to_cdot_hgvs,
)


def SDHD_description() -> Description:
    """SDHD, on the forward strand"""
    return init_description("ENST00000375549.8:c.=")


def WT1_description() -> Description:
    """WT1, on the reverse strand"""
    return init_description("ENST00000452863.10:c.=")


class TestVariant:
    """Test the functionality of the Variant class

    Tests concerning the integration with mutalyzer are excluded
    """

    def test_Variant_to_str(self) -> None:
        """Test converting a Variant to string"""
        v = Variant(10, 11, "ATG")
        assert str(v) == "Variant(start=10, end=11, inserted=ATG, deleted=)"

    def test_Variant_snp_to_str(self) -> None:
        """SNPS are special, since they contain the inserted sequence"""
        # 10A>T
        v = Variant(10, 11, "T", "A")
        assert str(v) == "Variant(start=10, end=11, inserted=T, deleted=A)"

    def test_Variant_to_model_positions(self) -> None:
        """Test converting a variant to model"""
        v = Variant(10, 11, "ATG")
        model = v.to_model()

        assert model["location"]["start"]["position"] == 10
        assert model["location"]["end"]["position"] == 11
        assert model["inserted"][0]["sequence"] == "ATG"

        # Deleted entry is missing for deletions
        assert "deleted" not in model

    def test_Variant_to_model_snp(self) -> None:
        """Test converting a snp Variant to model"""
        # 10 A>T
        v = Variant(10, 11, "T", "A")
        model = v.to_model()

        assert model["inserted"][0]["sequence"] == "T"
        assert model["deleted"][0]["sequence"] == "A"

    def test_Variant_deleted_snp_only(self) -> None:
        """Test that deleted is only defined for SNPs, not larger indels

        (To match the behaviour of Mutalyzer)
        """
        with pytest.raises(ValueError):
            # 10_12delinsGG
            Variant(10, 12, "GG", "AT")

    def test_Variant_no_inserted_sequence(self) -> None:
        """
        GIVEN a Variant with an empty inserted sequence
        WHEN we convert to a model
        THEN inserted should be an empty list
        """
        v = Variant(10, 11)
        model = v.to_model()

        assert model["inserted"] == []

    def test_Variant_end_after_start(self) -> None:
        """Raise an error when the end is after the start"""
        with pytest.raises(ValueError):
            Variant(11, 10)

    def test_Variant_equal(self) -> None:
        """Test that variants are equal to themselves"""
        v1 = Variant(1, 2, inserted="A")
        v2 = Variant(1, 2, inserted="A")

        assert v1 == v2
        assert v2 == v1

    ORDERING = [
        # Ends are touching
        (Variant(1, 3), Variant(3, 5)),
        # Gap between variants
        (Variant(0, 1), Variant(2, 4)),
    ]

    @pytest.mark.parametrize("a, b", ORDERING)
    def test_Variant_relative_positions(self, a: Variant, b: Variant) -> None:
        """Variant a is before variant b"""
        assert a.before(b)
        assert b.after(a)

    INSIDE = [
        # Variants are inside themselves
        (Variant(0, 3), Variant(0, 3)),
        # Smaller variant a is inside b
        (Variant(1, 3), Variant(0, 3)),
        (Variant(0, 2), Variant(0, 3)),
        (Variant(1, 2), Variant(0, 3)),
    ]

    @pytest.mark.parametrize("a, b", INSIDE)
    def test_Variant_a_inside_b(self, a: Variant, b: Variant) -> None:
        """Variant a is inside variant b"""
        assert a.inside(b)

    NOT_INSIDE = [
        # a starts outside of b
        (Variant(0, 3), Variant(1, 3)),
        # a ends outside of b
        (Variant(1, 4), Variant(1, 3)),
        # b is inside of a
        (Variant(0, 3), Variant(1, 2)),
        # a before b
        (Variant(0, 3), Variant(3, 5)),
        # a after b
        (Variant(3, 5), Variant(1, 3)),
    ]

    @pytest.mark.parametrize("a, b", NOT_INSIDE)
    def test_Variant_a_not_inside_b(self, a: Variant, b: Variant) -> None:
        """Variant a is not inside variant b"""
        assert not a.inside(b)

    OVERLAP = [
        # b ends inside a
        (Variant(2, 5), Variant(1, 3)),
        # b fully inside a
        (Variant(2, 5), Variant(3, 4)),
        # b fully inside a, ends at end
        (Variant(2, 5), Variant(3, 5)),
        # b ends inside a
        (Variant(2, 5), Variant(2, 4)),
        # b starts in a, extends after
        (Variant(2, 5), Variant(3, 6)),
        # b start before a
        (Variant(2, 5), Variant(1, 5)),
        # a is inside b
        (Variant(2, 5), Variant(1, 6)),
        # a is inside b
        (Variant(2, 5), Variant(2, 6)),
        # a equals b
        (Variant(2, 5), Variant(2, 5)),
    ]

    @pytest.mark.parametrize("a, b", OVERLAP)
    def test_Variant_a_overlaps_b(self, a: Variant, b: Variant) -> None:
        """Variant a and b overlap"""
        assert a.overlap(b)
        assert b.overlap(a)

    NO_OVERLAP = [
        # a is after b
        (Variant(2, 5), Variant(1, 2)),
        # a is before b
        (Variant(2, 5), Variant(5, 6)),
    ]

    @pytest.mark.parametrize("a, b", NO_OVERLAP)
    def test_Variant_a_does_not_overlap_b(self, a: Variant, b: Variant) -> None:
        """Variant a and b do not overlap"""
        assert not a.overlap(b)
        assert not b.overlap(a)

    def test_Variant_ordering(self) -> None:
        """Test sorting variants by start position"""
        v1 = Variant(10, 11)
        v2 = Variant(0, 1)
        assert sorted([v1, v2]) == [v2, v1]

    def test_Variant_ordering_error(self) -> None:
        """Overlapping variants cannot be sorted, and raise an error"""
        v1 = Variant(10, 11)
        v2 = Variant(10, 11)
        with pytest.raises(ValueError):
            sorted([v1, v2])

    def test_Variant_from_dict(self) -> None:
        """Test creating a Variant from a dict"""
        v = Variant(10, 22, inserted="ATCGG", deleted="")
        d = {"start": 10, "end": 22, "inserted": "ATCGG", "deleted": ""}

        assert Variant.from_dict(d) == v


class TestCombineVariantsDeletion:
    """Test combining multiple variants with a deletion

    This is used in the context of exon skipping, where variants that are fully
    inside the deletion (exon skip) are removed
    """

    def test_combine_variants_deletion_error(self) -> None:
        """Test that a value error is raised if the Variant is not a deletion"""
        indel = Variant(10, 15, inserted="ATC")
        with pytest.raises(ValueError):
            combine_variants_deletion(list(), indel)

    def test_combine_variants_deletion_empty(self) -> None:
        """If ther are no other variants, the results only contain the deletion"""
        variants: list[Variant] = list()
        deletion = Variant(0, 10)
        assert combine_variants_deletion(variants, deletion) == [deletion]

    # fmt: off
    COMBINE = [
        ( # The variants are out of order
            # Variants
            [Variant(5, 7), Variant(2, 4)],
            # Deletion
            Variant(10, 11),
            # Expected
            [Variant(2,4), Variant(5, 7), Variant(10, 11)],
        ),
        ( # The deletion is before both variants
            # Variants
            [Variant(2, 4), Variant(5, 7)],
            # Deletion
            Variant(0, 1),
            # Expected
            [Variant(0, 1), Variant(2,4), Variant(5, 7)],
        ),
        ( # The deletion is between the variants
            # Variants
            [Variant(2, 4), Variant(5, 7)],
            # Deletion
            Variant(4, 5),
            # Expected
            [Variant(2,4), Variant(4, 5), Variant(5, 7)],
        ),
        ( # The deletion is after both variants
            # Variants
            [Variant(2, 4), Variant(5, 7)],
            # Deletion
            Variant(10, 11),
            # Expected
            [Variant(2,4), Variant(5, 7), Variant(10, 11)],
        ),
        ( # The deletion contains the first variant
            # Variants
            [Variant(2, 4), Variant(5, 7)],
            # Deletion
            Variant(1, 5),
            # Expected
            [Variant(1, 5), Variant(5, 7)],
        ),
        ( # The deletion contains the second variant
            # Variants
            [Variant(2,4), Variant(5, 7)],
            # Deletion
            Variant(4, 11),
            # Expected
            [Variant(2,4), Variant(4, 11)],
        ),
        ( # The deletion contains both variants
            # Variants
            [Variant(2,4), Variant(5, 7)],
            # Deletion
            Variant(2, 7),
            # Expected
            [Variant(2, 7)],
        ),
    ]
    # fmt: on

    @pytest.mark.parametrize("variants, deletion, expected", COMBINE)
    def test_combine_variants_deletion(
        self,
        variants: Sequence[Variant],
        deletion: Variant,
        expected: Sequence[Variant],
    ) -> None:
        combined = combine_variants_deletion(variants, deletion)
        assert combined == expected

    def test_combine_variants_deletion_variants_overlap_eachother(self) -> None:
        """Test that we get a value error if the variants overlap"""
        variants = [Variant(0, 2), Variant(1, 3)]
        deletion = Variant(10, 11)
        with pytest.raises(ValueError):
            combine_variants_deletion(variants, deletion)

    def test_combine_variants_deletion_variants_partially_overlap_deletion(
        self,
    ) -> None:
        """
        Test that we get a value error if one the variants partially overlaps
        the deletion
        """
        variants = [Variant(2, 4)]
        deletion = Variant(3, 11)
        with pytest.raises(ValueError):
            combine_variants_deletion(variants, deletion)


class TestVariantMutalyzerIntegration:
    """Test the interaction between Variants and mutalyzer"""

    def test_Variant_from_model(self) -> None:
        """Test creating a variant from a mutalyzer delins model"""
        # fmt: off
        delins_model = {
            "type": "deletion_insertion",
            "source": "reference",
            "location": {
                "type": "range",
                "start": {
                    "type": "point",
                    "position": 0
                },
                "end": {
                    "type": "point",
                    "position": 2
                }
            },
            "inserted": [
                {
                    "sequence": "ATC",
                    "source": "description"
                }
            ]
        }
        # fmt: on
        assert Variant.from_model(delins_model) == Variant(0, 2, inserted="ATC")

    def test_Variant_from_model_deletion(self) -> None:
        """Test creating a variant from a mutalyzer delins model

        If nothing is inserted, the inserted sequence is an emtpy list, but should
        be made an empty string in the _Variant object
        """
        # fmt: off
        delins_model = {
            "type": "deletion_insertion",
            "source": "reference",
            "location": {
                "type": "range",
                "start": {
                    "type": "point",
                    "position": 0
                },
                "end": {
                    "type": "point",
                    "position": 2
                }
            },
            "inserted": [
                {
                    "sequence": [],
                    "source": "description"
                }
            ]
        }
        # fmt: on
        assert Variant.from_model(delins_model) == Variant(0, 2)


class TestVariantMutalyzerForward(object):
    """Test the interaction between Variants and mutalyzer for transcripts on the forward strand"""

    transcript = "ENST00000375549.8"

    # SDHD, forward strand
    def get_empty_transcript(self) -> Description:
        return init_description(f"{self.transcript}:c.=")

    def delins_from_description(self, d: Description) -> dict[str, Any]:
        """Return the internal delins model from a Description object"""
        delins_models = d.de_hgvs_internal_indexing_model["variants"]
        assert len(delins_models) == 1
        delins: dict[str, Any] = d.de_hgvs_internal_indexing_model["variants"][0]
        return delins

    def protein_from_description(self, d: Description) -> str:
        """Return the predicted protein sequence form a description"""
        protein_sequence: str = d.protein["predicted"]
        return protein_sequence

    def sequence_from_description(self, d: Description) -> str:
        """Return the sequence form a description"""
        return sequence_from_description(d)

    def protein_from_variant(self, v: Variant, d: Description) -> str:
        """Return the predicted protein sequence from a Variant"""

        ref_id = get_reference_id(d.corrected_model)

        selector_model = get_protein_selector_model(
            d.references[ref_id]["annotations"], ref_id
        )

        # Get the sequence (needed for complex variants)
        delins_models = [v.to_model()]

        predicted: str
        desc, ref, predicted, *rest = get_protein_description(
            delins_models, d.references, selector_model
        )
        return predicted

    # fmt: off
    MUTATIONS_VARIANT = [
        # HGVS, coordinates on the genome,
        # A simple missense that changes a single amino acids
        (
            "13T>A", # cdot notation of the variant, not used
            [Variant(start=47, end=48, inserted="A", deleted="T")], # Variants
            (112086919, 112086922) # Chromosomal location of the protein change
        ),
        # A stop mutation which destroys most of the protein
        (
            "9_10insTAG",
            [Variant(start=44, end=44, inserted="TAG")],
            (112086916, 112094967)
        ),
        # A frameshift that is restored by an insertion
        # Note this gives two separate, adjacent regions
        (
            "[9_10insA;20del]",
            [
                Variant(start=44, end=44, inserted="A"),
                Variant(start=54, end=55),
            ],
            [
                (112086919, 112086925),
                (112086925, 112086928),
            ]
        ),
        # A bigger deletion
        (
            "13_21del",
            [Variant(start=47, end=56)],
            (112086919, 112086928)
        ),
        # An SNP that creates a STOP codon
        (
            "14G>A",
            [Variant(start=48, end=49, inserted="A", deleted="G")],
            (112086919, 112094967)
        ),
    ]
    # fmt: on
    @pytest.mark.parametrize("cdot, variants, expected", MUTATIONS_VARIANT)
    def test_mutation_to_cds_effect(
        self, cdot: str, variants: Sequence[Variant], expected: tuple[int, int]
    ) -> None:
        """
        GIVEN a HGVS transcript description for a transcript on the reverse strand
        WHEN we determine the CDS effect
        THEN we should get genome coordinates
        """
        hgvs = "ENST00000375549.8:c.="
        d = init_description(hgvs)

        if not isinstance(expected, list):
            e = [expected]
        else:
            e = expected

        assert mutation_to_cds_effect(d, variants) == e

    # Variants where the delins model can be used to initialise Variant.from_model
    ROUND_TRIP_VARIANTS = [
        "10C>T",
        "10del",
        "10_11insA",
        "10_11delinsGG",
    ]

    @pytest.mark.parametrize("variant", ROUND_TRIP_VARIANTS)
    def test_Variant_hgvs_round_trip(
        self,
        variant: str,
    ) -> None:
        """
        For some variants, the conversion between the mutalyzer delins model
        and a Variant can be reversed, this is what we test here.

        Note that this only works for variants on the forward strand, on the
        reverse strand, only deletions can be converted round trip
        """
        d = init_description(f"{self.transcript}:c.{variant}")
        delins_model = d.delins_model["variants"][0]

        v = Variant.from_model(delins_model)
        assert v.to_model() == delins_model

    # Variant that are not simple delins in mutalyzer
    REWRITE_VARIANTS = [
        # Equivalent to delins 8_9delinsTTTT
        ("8_9T[4]", Variant(42, 44, inserted="TTTT")),
        # Equivalent to 10_10delinsCCC
        ("10C[3]", Variant(44, 45, inserted="CCC")),
        # Equivalent to 10_13delinsCTCTCTCT
        ("10_13CT[4]", Variant(44, 48, inserted="CTCTCTCT")),
        # Equivalent to 10_10delinsC
        ("10dup", Variant(44, 45, inserted="CC")),
        # Equivalent to 10_11delinsAG:
        ("10_11inv", Variant(44, 46, inserted="AG")),
    ]

    @pytest.mark.parametrize("variant_description, expected", REWRITE_VARIANTS)
    def test_delins_complex_Variant(
        self,
        variant_description: str,
        expected: Variant,
    ) -> None:
        """Convert a complex variant into a Variant

        Here, a complex variant is defined as a variant that is not represented
        as a simple delins model in Mutalyzer. When creating a
        Variant.from_model(), the variant representation from Mutalyzer is
        converted to an equivalent (but not equal) delins representation
        """
        d = init_description(f"{self.transcript}:c.{variant_description}")
        delins_model = d.delins_model["variants"][0]
        # Extract the sequence from the Description object
        _id = d.input_model["reference"]["id"]
        sequence = d.references[_id]["sequence"]["seq"]

        assert Variant.from_model(delins_model, sequence=sequence) == expected

    TO_HGVS = [
        # SNP
        (Variant(44, 45, "T", "C"), "10C>T"),
        # Deletion
        (Variant(44, 45), "10del"),
        # Insertion
        (Variant(45, 45, "A"), "10_11insA"),
        # Insertion/Deletion
        (Variant(44, 46, "GG"), "10_11delinsGG"),
    ]

    @pytest.mark.parametrize("variant, expected", TO_HGVS)
    def test_Variant_to_hgvs(self, variant: Variant, expected: str) -> None:
        d = SDHD_description()
        assert to_cdot_hgvs(d, [variant]) == expected

    NOT_SUPPORTED = [
        # Uncertain repeat size
        "8_9T[4_5]",
        # Uncertain repeat start
        "(6_8)_9T[4]",
        # Uncertain repeat end
        "8_(9_10)T[4]",
        # Uncertain start position
        "(9_15)insA",
        # Uncertain end position
        "8_(9_10)del"
        # Insertion of a range
        "9_10ins14_20",
    ]

    @pytest.mark.parametrize("variant", NOT_SUPPORTED)
    def test_variant_not_supported(self, variant: str) -> None:
        """Test that we throw a NotImplemented error for complex variants"""
        with pytest.raises(Exception):
            init_description(f"{self.transcript}:c.{variant}")

    VARIANTS = [
        # A SNP
        ("10C>T", Variant(44, 45, inserted="T", deleted="C")),
        # A deletion
        ("10del", Variant(44, 45, inserted="")),
        # An insertion
        ("10_11insA", Variant(45, 45, inserted="A")),
        # Delins version of 10C>T (Note that the deleted part is lost)
        ("10_10delinsT", Variant(44, 45, inserted="T")),
        # A duplication
        ("10dup", Variant(45, 45, inserted="C")),
        # The same duplication, with Variant as a delins (the deleted part is
        # implicit)
        ("10dup", Variant(44, 45, inserted="CC")),
        # A duplication
        ("10_11dup", Variant(44, 44, inserted="CT")),
        # The same duplication, where Variant deletes the first "CT", and then
        # inserts it twice
        ("10_11dup", Variant(44, 46, inserted="CTCT")),
        # Inversion, equivalent to 10C>G
        ("10_10inv", Variant(44, 45, inserted="G")),
        ("10C>G", Variant(44, 45, inserted="G")),
        # Inversion, equivalent to 10_11delinsAG
        ("10_11inv", Variant(44, 46, inserted="AG")),
        ("10_11delinsAG", Variant(44, 46, inserted="AG")),
        # Inversion, not symetrical
        ("18_20inv", Variant(52, 55, inserted="AGC")),
        ("18_20delinsAGC", Variant(52, 55, inserted="AGC")),
        # Small mononucleotide repeat:
        # ref: GGTT  CTCT
        # obs: GGTTTTCTCT
        # hgvs notation: 9_10insTT, 8_9T[4]
        ("8_9T[4]", Variant(44, 44, inserted="TT")),
        # Repeat of multiple distinct nucleotides
        # ref: GGTTCTCT    GGA
        # obs: GGTTCTCTCTCTGGA
        # hgvs notation: 9_10insCTCT, 10_13CT[4]
        ("10_13CT[4]", Variant(44, 44, inserted="CTCT")),
    ]

    @pytest.mark.parametrize("hgvs,variant", VARIANTS)
    def test_hgvs_Variant_equivalence_via_protein(
        self, hgvs: str, variant: Variant
    ) -> None:
        """Test hgvs and Variant equivalence by comparing the protein prediction

        The goal here is to verify that the model representation of the Variant
        is usable by mutalyzer in the same way as the original HGVS description
        """
        d = init_description(f"{self.transcript}:c.{hgvs}")
        variant_protein = self.protein_from_variant(
            variant, self.get_empty_transcript()
        )
        description_protein = self.protein_from_description(d)

        assert variant_protein == description_protein

    VARIANTS = [
        ("10C>T", Variant(44, 45, inserted="T", deleted="C")),
        ("10del", Variant(44, 45, inserted="")),
        ("10_11insA", Variant(45, 45, inserted="A")),
        ("10_10delinsT", Variant(44, 45, inserted="T")),
        ("10dup", Variant(44, 45, inserted="CC")),
        ("10_11dup", Variant(44, 46, inserted="CTCT")),
        ("10_10inv", Variant(44, 45, inserted="G")),
        ("10_11inv", Variant(44, 46, inserted="AG")),
        ("18_20inv", Variant(52, 55, inserted="AGC")),
        # Repeat is a delins of the existing repeat units, and insertion of the
        # new repeat
        ("8_9T[4]", Variant(42, 44, inserted="TTTT")),
        ("10_13CT[4]", Variant(44, 48, inserted="CTCTCTCT")),
    ]

    @pytest.mark.parametrize("hgvs, variant", VARIANTS)
    def test_hgvs_Variant_delins_model(self, hgvs: str, variant: Variant) -> None:
        """Test hgvs and Variant delins model equivalence directly

        The goal here is to verify that the Variant.from_model structure is
        always a simple delins model, even for duplications, repeats and
        inversions.
        """
        d = init_description(f"{self.transcript}:c.{hgvs}")
        delins_model = self.delins_from_description(d)

        seq = self.sequence_from_description(d)
        variant_model = Variant.from_model(delins_model, sequence=seq)
        assert variant_model == variant

    COORDINATES = [
        (Variant(47, 50), (112_086_919, 112_086_922)),
        (Variant(53, 56), (112_086_925, 112_086_928)),
        (Variant(44, 8095), (112_086_916, 112_094_967)),
        (Variant(47, 56), (112_086_919, 112_086_928)),
        (Variant(47, 8095), (112_086_919, 112_094_967)),
    ]

    @pytest.mark.parametrize("variant, genomic_coordinates", COORDINATES)
    def test_Variant_to_genomic_coordinates(
        self, variant: Variant, genomic_coordinates: tuple[int, int]
    ) -> None:
        """Test converting Variant coordinates to genomic coordinates"""
        d = init_description(f"{self.transcript}:c.=")
        assert variant.genomic_coordinates(d) == genomic_coordinates


class TestVariantMutalyzerReverse(TestVariantMutalyzerForward):
    """Test the interaction between Variants and mutalyzer for transcripts on the reverse strand"""

    # WT1, reverse strand
    transcript = "ENST00000452863.10"

    # fmt: off
    MUTATIONS_VARIANT = [
        # Variant, coordinates on the genome
        # A simple missense that changes a single amino acids
        (
            "13T>A",
            [Variant(start=47573, end=47574, inserted="T", deleted="A")],
            (32435345, 32435348)
        ),
        # A stop mutation which destroys most of the protein
        (
            "9_10insTAG",
            [Variant(start=47577, end=47577, inserted="CTA")],
            (32389060, 32435351)
        ),
        # # # A frameshift that is restored by an insertion
        (
            "[10del;20_21insA]",
            [
                Variant(start=47576, end=47577),
                Variant(start=47566, end=47566, inserted="T"),
            ],
            (32435339, 32435351)
        ),
        # # # A frameshift that is restored by a bigger insertion
        (
            "[10del;20_21insATCGAATATGGGG]",
            [
                Variant(start=47566, end=47566, inserted="CCCCATATTCGAT"),
                Variant(start=47576, end=47577),
            ],
            (32435339, 32435351)),
        # # # A bigger deletion
        (
            "11_19del",
            [Variant(start=47567, end=47576)],
             (32435342, 32435351)
        ),
        # # # An inframe deletion that creates a STOP codon
        (
            "87_89del",
            [Variant(start=47497, end=47500)],
            (32389060, 32435276)
        ),
    ]
    # fmt: on
    @pytest.mark.parametrize("cdot, variants, expected", MUTATIONS_VARIANT)
    def test_mutation_to_cds_effect(
        self, cdot: str, variants: Sequence[Variant], expected: tuple[int, int]
    ) -> None:
        """
        GIVEN a HGVS transcript description for a transcript on the reverse strand
        WHEN we determine the CDS effect
        THEN we should get genome coordinates
        """
        hgvs = "ENST00000452863.10:c.="
        d = init_description(hgvs)

        assert mutation_to_cds_effect(d, variants) == [expected]

    @pytest.mark.parametrize("variant", ["10del"])
    def test_Variant_hgvs_round_trip(self, variant: str) -> None:
        """
        See the equivalent test for forward transcripts. For transcripts on the
        reverse strand, a round-trip conversion can only be performed for deletions
        """
        d = init_description(f"{self.transcript}:c.{variant}")
        delins_model = d.delins_model["variants"][0]

        v = Variant.from_model(delins_model)
        assert v.to_model() == delins_model

    # Variant that are not simple delins in mutalyzer
    REWRITE_VARIANTS = [
        # Equivalent to 10_13delinsCTCTCTCT
        (
            "10_13CT[4]",
            Variant(47573, 47577, inserted="AGAGAGAG"),
        ),
        # Equivalent to 10_10delinsCC
        ("10dup", Variant(47576, 47577, inserted="GG")),
    ]

    @pytest.mark.parametrize("variant_description, expected", REWRITE_VARIANTS)
    def test_delins_complex_Variant(
        self,
        variant_description: str,
        expected: Variant,
    ) -> None:
        """Convert a complex variant into a Variant

        Here, a complex variant is defined as a variant that is not represented
        as a simple delins model in Mutalyzer. When creating a
        Variant.from_model(), the variant representation from Mutalyzer is
        converted to an equivalent (but not equal) delins representation

        Note that for transcripts on the reverse strand, only deletions are
        represented by a pure delins in Mutalyzer
        """
        d = init_description(f"{self.transcript}:c.{variant_description}")
        delins_model = d.delins_model["variants"][0]
        # Extract the sequence from the Description object
        _id = d.input_model["reference"]["id"]
        sequence = d.references[_id]["sequence"]["seq"]

        assert Variant.from_model(delins_model, sequence=sequence) == expected

    TO_HGVS = [
        # SNP
        (Variant(47573, 47574, inserted="T", deleted="A"), "13T>A"),
        # Deletion
        (Variant(start=47576, end=47577), "10del"),
        # Insertion
        (Variant(start=47577, end=47577, inserted="CTA"), "9_10insTAG"),
        # Insertion/deletion
        (Variant(47566, 47569, inserted="GCA"), "18_20delinsTGC"),
    ]

    @pytest.mark.parametrize("variant, expected", TO_HGVS)
    def test_Variant_to_hgvs(self, variant: Variant, expected: str) -> None:
        d = WT1_description()
        assert to_cdot_hgvs(d, [variant]) == expected

    NOT_SUPPORTED = [
        # Uncertain repeat size
        "8_9T[4_5]",
        # Uncertain repeat start
        "(6_8)_9T[4]",
        # Uncertain repeat end
        "8_(9_10)T[4]",
        # Uncertain start position
        "(9_15)insA",
        # Uncertain end position
        "8_(9_10)del"
        # Insertion of a range
        "9_10ins14_20",
    ]

    @pytest.mark.parametrize("variant", NOT_SUPPORTED)
    def test_variant_not_supported(self, variant: str) -> None:
        """Test that we throw a NotImplemented error for complex variants"""
        with pytest.raises(Exception):
            init_description(f"{self.transcript}:c.{variant}")

    # Note that on the Variant the nucleotides are (reverse?) complemented
    VARIANTS = [
        # A SNP,
        ("10C>T", Variant(47576, 47577, inserted="A", deleted="G")),
        # A deletion
        ("10del", Variant(47576, 47577, inserted="")),
        # An insertion
        ("10_11insA", Variant(47576, 47576, inserted="T")),
        # Delins version of 10C>T (Note that the deleted part is lost)
        ("10_10delinsT", Variant(47576, 47577, inserted="A")),
        # A duplication
        ("10dup", Variant(47576, 47576, inserted="G")),
        # The same duplication, with Variant as a delins (the deleted part is
        # implicit)
        ("10dup", Variant(47575, 47576, inserted="GG")),
        ("10_10delinsCC", Variant(47576, 47577, inserted="GG")),
        # A duplication, note that "AG" is the reverse complement of "CT"
        ("10_11dup", Variant(47575, 47575, inserted="AG")),
        ("11_12insCT", Variant(47575, 47575, inserted="AG")),
        ("12_13delinsCTCT", Variant(47575, 47577, inserted="AGAG")),
        ("12_13delinsCTCT", Variant(47573, 47575, inserted="AGAG")),
        # Inversion, equivalent to 10C>G
        ("10_10inv", Variant(47576, 47577, inserted="C")),
        ("10C>G", Variant(47576, 47577, inserted="C")),
        # Inversion, equivalent to 10_11delinsAG
        ("10_11inv", Variant(47575, 47577, inserted="CT")),
        ("10_11delinsAG", Variant(47575, 47577, inserted="CT")),
        # Inversion, not symetrical
        ("18_20inv", Variant(47566, 47569, inserted="GCA")),
        ("18_20delinsTGC", Variant(47566, 47569, inserted="GCA")),
        # Small mononucleotide repeat
        # 7_8T[4],  Variant representation = 7_8delinsTTTT
        ("7_8T[4]", Variant(47578, 47580, inserted="AAAA")),
        # Repeat of multiple distinct nucleotides
        # 10_13CT[4], Variant representation = 10_13delinsCTCTCTCT
        ("10_13CT[4]", Variant(47573, 47577, inserted="AGAGAGAG")),
    ]

    @pytest.mark.parametrize("hgvs,variant", VARIANTS)
    def test_hgvs_Variant_equivalence_via_protein(
        self, hgvs: str, variant: Variant
    ) -> None:
        """Test hgvs and Variant equivalence by comparing the protein prediction

        The goal here is to verify that the model representation of the Variant
        is usable by mutalyzer in the same way as the original HGVS description
        """
        d = init_description(f"{self.transcript}:c.{hgvs}")
        variant_protein = self.protein_from_variant(
            variant, self.get_empty_transcript()
        )
        description_protein = self.protein_from_description(d)

        assert variant_protein == description_protein

    VARIANTS = [
        ("10C>T", Variant(47576, 47577, inserted="A", deleted="G")),
        ("10del", Variant(47576, 47577, inserted="")),
        ("10_11insA", Variant(47576, 47576, inserted="T")),
        ("10_10delinsT", Variant(47576, 47577, inserted="A")),
        # Delete the original 'C', and insert CC (reverse complement)
        ("10dup", Variant(47576, 47577, inserted="GG")),
        ("10_10delinsCC", Variant(47576, 47577, inserted="GG")),
        # Delete the original "TC", and insert CTCT (reverse complement)
        ("10_11dup", Variant(47575, 47577, inserted="AGAG")),
        ("12_13delinsCTCT", Variant(47573, 47575, inserted="AGAG")),
        # Inversion, equivalent to 10C>G
        ("10_10inv", Variant(47576, 47577, inserted="C")),
        # Inversion, equivalent to 10_11delinsAG
        ("10_11inv", Variant(47575, 47577, inserted="CT")),
        # Inversion, not symetrical
        ("18_20inv", Variant(47566, 47569, inserted="GCA")),
        ("18_20delinsTGC", Variant(47566, 47569, inserted="GCA")),
        # Small mononucleotide repeat
        # 7_8T[4],  Variant representation = 7_8delinsTTTT
        ("7_8T[4]", Variant(47578, 47580, inserted="AAAA")),
        # Repeat of multiple distinct nucleotides
        # 10_13CT[4], Variant representation = 10_13delinsCTCTCTCT
        ("10_13CT[4]", Variant(47573, 47577, inserted="AGAGAGAG")),
    ]

    @pytest.mark.parametrize("hgvs, variant", VARIANTS)
    def test_hgvs_Variant_delins_model(self, hgvs: str, variant: Variant) -> None:
        """Test hgvs and Variant delins model equivalence directly

        The goal here is to verify that the Variant.from_model structure is
        always a simple delins model, even for duplications, repeats and
        inversions.
        """
        d = init_description(f"{self.transcript}:c.{hgvs}")
        delins_model = self.delins_from_description(d)

        seq = self.sequence_from_description(d)
        variant_model = Variant.from_model(delins_model, sequence=seq)
        assert variant_model == variant

    COORDINATES = [
        (Variant(47571, 47574), (32_435_345, 32_435_348)),
        (Variant(1286, 47577), (32_389_060, 32_435_351)),
        (Variant(47565, 47577), (32_435_339, 32_435_351)),
        (Variant(47565, 47577), (32_435_339, 32_435_351)),
        (Variant(47568, 47577), (32_435_342, 32_435_351)),
        (Variant(1286, 47502), (32_389_060, 32_435_276)),
    ]

    @pytest.mark.parametrize("variant, genomic_coordinates", COORDINATES)
    def test_Variant_to_genomic_coordinates(
        self, variant: Variant, genomic_coordinates: tuple[int, int]
    ) -> None:
        """Test converting Variant coordinates to genomic coordinates"""
        d = init_description(f"{self.transcript}:c.=")
        assert variant.genomic_coordinates(d) == genomic_coordinates
