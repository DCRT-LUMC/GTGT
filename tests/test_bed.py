import pytest
from typing import List

from GTGT import Bed


@pytest.fixture
def bed() -> Bed:
    return Bed("chr1", 0, 11)


def test_make_bed(bed: Bed) -> None:
    assert bed.chrom == "chr1"
    assert bed.chromStart == 0
    assert bed.chromEnd == 11


def test_default_values_simple_attributes(bed: Bed) -> None:
    assert bed.name == "."
    assert bed.score == 0
    assert bed.strand == "."


def test_default_values_display_attributes(bed: Bed) -> None:
    assert bed.thickStart == 0
    assert bed.thickEnd == 11
    assert bed.itemRgb == (0, 0, 0)


def test_rgb_zero() -> None:
    """itemRgb of zero is a special case, and an alias for (0, 0, 0)"""
    bed = Bed("chr1", 0, 11, "name", 1000, "+", 0, 11, itemRgb="0")

    assert bed.itemRgb == (0, 0, 0)


def test_default_values_blocks(bed: Bed) -> None:
    assert bed.blockCount == 1
    assert bed.blockSizes == [11]
    assert bed.blockStarts == [0]


def test_blocks_interface(bed: Bed) -> None:
    for start, end in bed.blocks():
        assert start == 0
        assert end == 11


# Bed records, and their corresponding representation in Bed format
bed_records = [
    (Bed("chr1", 0, 11), ["chr1", "0", "11", ".", "0", ".", "0", "11", "0,0,0", "1"]),
    (
        Bed(
            chrom="chr1",
            chromStart=0,
            chromEnd=11,
            name="name",
            score=5,
            strand="+",
            thickStart=8,
            thickEnd=10,
            itemRgb="42,42,42",
            blockCount=2,
        ),
        ["chr1", "0", "11", "name", "5", "+", "8", "10", "42,42,42", "2"],
    ),
]


@pytest.mark.parametrize("bed, expected", bed_records)
def test_str_bed(bed: Bed, expected: List[str]) -> None:
    assert str(bed) == "\t".join(expected)
