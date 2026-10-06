from typing import Sequence

import pytest

from gtgt import Bed
from gtgt.mutalyzer import init_description
from gtgt.therapy import Therapy
from gtgt.transcript import Features
from gtgt.variant import Variant


@pytest.fixture
def exons() -> Bed:
    exons = [(0, 10), (20, 40), (50, 60), (70, 100)]
    bed = Bed.from_blocks("chr1", exons)
    bed.name = "Exons"

    return bed


@pytest.fixture
def coding_exons() -> Bed:
    coding_exons = [(23, 40), (50, 60), (70, 72)]
    bed = Bed.from_blocks("chr1", coding_exons)
    bed.name = "Coding exons"

    return bed


@pytest.fixture
def features(exons: Bed, coding_exons: Bed) -> Features:
    """
    Bed records that make up a transcript
    Each positions shown here is 10x
    (i) means inferred by the init method

                  0 1 2 3 4 5 6 7 8 9
    exons         -   - -   -   - - -
    coding_exons      - -   -   -
    """
    return Features(rna=[exons], protein=[coding_exons])


def test_features_init(features: Features) -> None:
    assert features.rna[0].name == "Exons"
    assert features.protein[0].name == "Coding exons"


def test_empty_features() -> None:
    """Test creating and working with an empty transcript"""
    t = Features(rna=[], protein=[])

    # Test if we can get the records
    assert t.records() == []

    # Test if intersect works
    t.intersect(Bed("chr1", 10, 20))

    # Test if subtraction works
    t.subtract(Bed("chr1", 10, 20))


intersect_selectors = [
    # Selector spans all exons
    (
        Bed("chr1", 0, 100),
        Bed(
            "chr1",
            0,
            100,
            name="Exons",
            blockSizes=[10, 20, 10, 30],
            blockStarts=[0, 20, 50, 70],
        ),
    ),
    # Selector on a different chromosome
    (Bed("chr2", 0, 100), Bed("chr1", 0, 0)),
    # Selector intersect the first exon
    (Bed("chr1", 5, 15), Bed("chr1", 5, 10)),
    # Selector intersects the last base of the first exon,
    # and the first base of the second exon
    (Bed("chr1", 9, 21), Bed("chr1", 9, 21, blockSizes=[1, 1], blockStarts=[0, 11])),
]


def test_features_init_no_coding(exons: Bed) -> None:
    t = Features(rna=[exons], protein=[])
    assert not t.coding_exons


@pytest.mark.parametrize("selector, exons", intersect_selectors)
def test_intersect_features(selector: Bed, exons: Bed, features: Features) -> None:
    """Test if intersecting the Transcript updates the exons"""
    features.intersect(selector)

    # Ensure the name matches, it's less typing to do that here
    exons.name = "Exons"
    assert features.exons == exons


overlap_selectors = [
    # Selector spans all exons
    (
        Bed("chr1", 0, 100),
        Bed(
            "chr1",
            0,
            100,
            name="Exons",
            blockSizes=[10, 20, 10, 30],
            blockStarts=[0, 20, 50, 70],
        ),
    ),
    # Selector on a different chromosome
    (Bed("chr2", 0, 100), Bed("chr1", 0, 0)),
    # Selector intersect the first exon
    (Bed("chr1", 5, 15), Bed("chr1", 0, 10)),
    # Selector intersects the last base of the first exon,
    # and the first base of the second exon
    (Bed("chr1", 9, 21), Bed("chr1", 0, 40, blockSizes=[10, 20], blockStarts=[0, 20])),
]


@pytest.mark.parametrize("selector, exons", overlap_selectors)
def test_overlap_transcript(selector: Bed, exons: Bed, features: Features) -> None:
    """Test if overlapping the Transcript updates the exons"""
    features.overlap(selector)

    # Ensure the name matches, it's less typing to do that here
    exons.name = "Exons"
    assert features.exons == exons


subtract_selectors = [
    # Selector spans all exons
    (Bed("chr1", 0, 100), Bed("chr1", 0, 0)),
    # Selector on a different chromosome
    (
        Bed("chr2", 0, 100),
        Bed(
            "chr1",
            0,
            100,
            name="Exons",
            blockSizes=[10, 20, 10, 30],
            blockStarts=[0, 20, 50, 70],
        ),
    ),
    # Selector intersect the first exon
    (
        Bed("chr1", 5, 15),
        Bed("chr1", 0, 100, blockSizes=[5, 20, 10, 30], blockStarts=[0, 20, 50, 70]),
    ),
    # Selector intersects the last base of the first exon,
    # and the first base of the second exon
    (
        Bed("chr1", 9, 21),
        Bed("chr1", 0, 100, blockSizes=[9, 19, 10, 30], blockStarts=[0, 21, 50, 70]),
    ),
]


@pytest.mark.parametrize("selector, exons", subtract_selectors)
def test_subtract_features(selector: Bed, exons: Bed, features: Features) -> None:
    """Test if subtracting the Transcript updates the exons"""
    features.subtract(selector)

    # Ensure the name matches, it's less typing to do that here
    exons.name = "Exons"
    assert features.exons == exons


def test_compare_features(features: Features, coding_exons: Bed) -> None:
    exon_blocks = [
        (0, 10),
        # (20, 40),  # Missing the second exon
        (50, 60),
        (70, 100),
    ]
    exons = Bed.from_blocks("chr1", exon_blocks)
    exons.name = "Exons"

    coding_blocks = [
        # (23, 40),  # Missing the second exon
        (50, 60),
        (70, 72),
    ]
    coding_exons = Bed.from_blocks("chr1", coding_blocks)
    coding_exons.name = "Coding exons"

    smaller = Features(rna=[exons], protein=[coding_exons])

    cmp = smaller.compare(features)

    assert cmp[0].percentage == pytest.approx(0.71, abs=0.01)
    assert cmp[1].percentage == pytest.approx(0.41, abs=0.01)
