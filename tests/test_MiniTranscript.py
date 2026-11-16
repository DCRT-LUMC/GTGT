import pytest

from gtgt.transcript import CodingPoint, MiniTranscript
from gtgt.variant import Variant

"""
Tests for getting the HGVS description from MiniTranscript use this fake example transcript

The CDS starts in exon 2, and ends in exon 3

100    110      125       160    172        208    241        300    348
v       v        v         v      v          v      v          v      v
--------||||||||||---------||||||||----------||||||||----------||||||||
          exon4    intron3  exon3   intron2   exon2   intron1    exon1
        ^        ^         ^      ^          ^      ^          ^      ^
        r.112    r.97      r.96   r.84      r.83    r.50       r.49   r.1

"""

Range = tuple[int, int]


@pytest.fixture
def mini_transcript() -> MiniTranscript:
    exons = [(110, 125), (160, 172), (208, 241), (300, 348)]
    cds = [(163, 235)]
    return MiniTranscript(
        inverted=True, sequence="", exons=exons, cds=cds, variants=list()
    )


@pytest.mark.parametrize(
    "variants, expected",
    [
        (  # No variants
            [],
            [(10, 20), (25, 35), (47, 50)],
        ),
        (  # SNPS don't change the size
            [Variant(10, 11, "A", "T"), Variant(20, 21, "A", "T")],
            [(10, 20), (25, 35), (47, 50)],
        ),
        (  # Variants after the features change nothing
            [Variant(70, 99)],
            [(10, 20), (25, 35), (47, 50)],
        ),
        (  # Insertion of 4 bp before the first feature
            [Variant(5, 5, inserted="ATCG")],
            [(14, 24), (29, 39), (51, 54)],
        ),
        (
            # Insertion of 4 bp just before the start of the second feature
            [Variant(24, 24, inserted="ATCG")],
            [(10, 20), (29, 39), (51, 54)],
        ),
        (
            # Insertion of 4 bp at the start of the second feature
            [Variant(25, 25, inserted="ATCG")],
            [(10, 20), (29, 39), (51, 54)],
        ),
        (
            # Insertion of 4 bp after the start of the second feature
            [Variant(26, 26, inserted="ATCG")],
            [(10, 20), (25, 39), (51, 54)],
        ),
        (
            # Insertion of 4 bp at the start of the second feature
            # Deletion at start of the third feature
            [Variant(25, 25, inserted="ATCG"), Variant(40, 45)],
            [(10, 20), (29, 39), (46, 49)],
        ),
        (
            # Insertion of 4 bp at the start of the second feature
            # Deletion inside the third feature
            [Variant(25, 25, inserted="ATCG"), Variant(48, 50)],
            [(10, 20), (29, 39), (51, 52)],
        ),
    ],
)
def test_mutate_features(variants: list[Variant], expected: list[Range]) -> None:
    """Test mutating features (such as exons, CDS)"""
    features = [(10, 20), (25, 35), (47, 50)]
    assert MiniTranscript.mutate_features(features, variants) == expected
    # Test the variants in the reverse order
    assert MiniTranscript.mutate_features(features, variants[::-1]) == expected


@pytest.mark.parametrize(
    "variants, expected",
    # fmt: off
    [
        (   # No variants
            [],
            "0123456789"
        ),
        (   # An insertion
            [Variant(0, 0, inserted="A")],
            "A0123456789"
        ),
        (   # A SNPS
            [Variant(4, 5, inserted="A")],
            "0123A56789"
        ),
        (   # A deletion
            [Variant(6, 8)],
            "01234589"
        ),
        (  # Multiple variants
            # Deletion     Delins
            [Variant(0, 3), Variant(5, 9, inserted="AAAT")],
            "34AAAT9",
        ),
        (
            # An insertion at the exact end of the sequence
            [Variant(10, 10, "ATG")],
            "0123456789ATG",
        ),
    ],
)
# fmt: on
def test_mutate_sequence(variants: list[Variant], expected: str) -> None:
    """Test mutating a sequence"""
    sequence = "0123456789"
    assert MiniTranscript.mutate_sequence(sequence, variants) == expected


def test_mutate_sequence_variant_outside() -> None:
    """Test a mutation outside of the sequence"""
    sequence = "0123456789"
    variant = Variant(11, 11, "ATG")
    with pytest.raises(ValueError):
        MiniTranscript.mutate_sequence(sequence, [variant])


def test_mutate_sequence_variants_not_sorted() -> None:
    """Test that we fail if the variants are not sorted"""
    sequence = "0123456789"
    variants = [Variant(5, 5, "ATG"), Variant(1, 2, inserted="ATG")]
    with pytest.raises(ValueError):
        MiniTranscript.mutate_sequence(sequence, variants)


class TestCodingPoint:
    @pytest.mark.parametrize(
        "point, hgvs_position",
        [
            # Positions in the coding region
            (CodingPoint(27, 0, 0, 0), "27"),
            # Position before the coding region
            (CodingPoint(-1, 0, -1, 0), "-1"),
            # Position after the coding region
            (CodingPoint(1, 0, 1, 0), "*1"),
            # Intronic position in the coding region
            (CodingPoint(15, 1, 0, 0), "15+1"),
            (CodingPoint(15, 4, 0, 0), "15+4"),
            (CodingPoint(16, -8, 0, 0), "16-8"),
            # Intronic position before the coding region
            (CodingPoint(-15, +1, -1, 0), "-15+1"),
            (CodingPoint(-16, -8, -1, 0), "-16-8"),
            # Intronic position after the coding region
            (CodingPoint(15, +1, 1, 0), "*15+1"),
            (CodingPoint(16, -8, 1, 0), "*16-8"),
        ],
    )
    def test_coding_point_to_hgvs(self, point: CodingPoint, hgvs_position: str) -> None:
        assert str(point) == hgvs_position

    def test_error_position_outside_transcript(self) -> None:
        with pytest.raises(ValueError):
            CodingPoint(1, 0, 0, 1)
