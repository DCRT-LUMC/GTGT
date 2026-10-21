import pytest

from gtgt.variant import Variant

from sequence import mutate_features, _variant_size, mutate_sequence

Range = tuple[int, int]


@pytest.mark.parametrize(
    "variant, size",
    [
        # SNP
        (Variant(start=10, end=11, inserted="A", deleted="T"), 0),
        # Insertion
        (Variant(start=10, end=10, inserted="A", deleted=""), 1),
        # Large insertion
        (Variant(start=10, end=10, inserted="ATC", deleted=""), 3),
        # Insertion of size 0
        (Variant(start=10, end=10, inserted="", deleted=""), 0),
        # Deletion
        (Variant(start=10, end=11), -1),
        # Large deletion
        (Variant(start=10, end=15), -5),
        # indel of size 0
        (Variant(start=10, end=15, inserted="ATCGA"), 0),
        # indel of positive size
        (Variant(start=10, end=15, inserted="ATCGA+++"), 3),
        # indel of negative size
        (Variant(start=10, end=15, inserted="AT"), -3),
    ],
)
def test_variant_size(variant: Variant, size: int) -> None:
    assert _variant_size(variant) == size


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
         # Insertion of 4 bp at the start of the second feature
            [Variant(25, 25, inserted="ATCG")],
            [(10, 20), (25, 39), (51, 54)],
        ),
    ],
)
def test_mutate_features(variants: list[Variant], expected:list[Range]) -> None:
    """Test mutating features (such as exons, CDS)"""
    features = [(10, 20), (25, 35), (47, 50)]
    assert mutate_features(features, variants) == expected


@pytest.mark.parametrize(
    "variants, expected",
    [
        (  # No variants
            [],
            "0123456789"
        ),
        (  # An insertion
            [Variant(0, 0, inserted="A")],
            "A0123456789"
        ),
        (  # A SNPS 
            [Variant(4, 5, inserted="A")],
            "0123A56789"
        ),
        (  # A deletion 
            [Variant(6, 8)],
            "01234589"
        ),

    ],
)
def test_mutate_sequence(variants: list[Variant], expected: str) -> None:
    """Test mutating a sequence"""
    sequence = "0123456789"
    assert mutate_sequence(sequence, variants) == expected
