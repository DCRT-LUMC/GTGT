import pytest
from typing import List, Tuple

from GTGT import Bed
from GTGT.bed import _to_range, _range_to_start_size


@pytest.fixture
def bed() -> Bed:
    return Bed("chr1", 0, 11)


def test_bed_properties(bed: Bed) -> None:
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


def test_bed_init_method() -> None:
    """Test various values for the init methods

    We allow people to pass both python values as well as their BED format
    string representation
    """
    # Pass Integers for positions
    assert Bed("chr1", 0, 10) == Bed("chr1", "0", "10")

    # Pass integer for score
    assert Bed("chr1", 0, 10, score=1000) == Bed("chr1", "0", "10", score="1000")

    # Pass integer for thickStart, thickEnd
    bed1 = Bed("chr1", 0, 10, thickStart=0, thickEnd=10)
    bed2 = Bed("chr1", "0", "10", thickStart="0", thickEnd="10")
    assert bed1 == bed2

    # Pass itemRgb as a tuple[int] and a string
    bed1 = Bed("chr1", 0, 10, itemRgb=(0, 0, 0))
    bed2 = Bed("chr1", "0", "10", itemRgb="0,0,0")
    assert bed1 == bed2

    # Pass an integer for blockCount
    bed1 = Bed("chr1", 0, 10, blockCount=10)
    bed2 = Bed("chr1", "0", "10", blockCount="10")
    assert bed1 == bed2

    # Pass a list of integers for blockSizes
    bed1 = Bed("chr1", 0, 10, blockSizes=[1, 4, 10])
    bed2 = Bed("chr1", "0", "10", blockSizes="1,4,10")
    assert bed1 == bed2

    # Pass a list of integers for blockStarts
    bed1 = Bed("chr1", 0, 10, blockStarts=[1, 4, 10])
    bed2 = Bed("chr1", "0", "10", blockStarts="1,4,10")
    assert bed1 == bed2

    # Trailing comma's are allowed for blockSizes and blockStarts, according
    # to the Bed standard
    bed1 = Bed("chr1", 0, 10, blockSizes=[1, 4, 10])
    bed2 = Bed("chr1", "0", "10", blockSizes="1,4,10,")
    assert bed1 == bed2

    bed1 = Bed("chr1", 0, 10, blockStarts=[1, 4, 10])
    bed2 = Bed("chr1", "0", "10", blockStarts="1,4,10,")
    assert bed1 == bed2


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
    (
        # Default values
        Bed("chr1", 0, 11, ".", 0),
        "chr1	0	11	.	0	.	0	11	0,0,0	1	11	0",
    ),
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
            itemRgb=(42, 42, 42),
            blockCount=2,
            blockSizes=[3, 4],
            blockStarts=[2, 7],
        ),
        "chr1	0	11	name	5	+	8	10	42,42,42	2	3,4	2,7",
    ),
]


@pytest.mark.parametrize("bed, line", bed_records)
def test_bed_roundtrip(bed: Bed, line: str) -> None:
    """Test writing a bed to string format"""
    # Convert bed record to line
    from_bed = str(bed)
    # Convert line to Bed record
    from_line = Bed(*line.split("\t"))

    # Check that the line from Bed is as expected
    assert from_bed == line
    # Check that the Bed record from line is as expected
    assert from_line == bed


intersect = [
    # Range A, range A, intersection
    ((0, 10), (10, 20), list()),
    ((0, 10), (0, 10), [(0, 10)]),
    # Test cases where A is of size 1, and B of size 3
    # In each test case, we shift A one further to the right
    ((0, 1), (1, 4), list()),
    ((1, 2), (1, 4), [(1, 2)]),
    ((2, 3), (1, 4), [(2, 3)]),
    ((3, 4), (1, 4), [(3, 4)]),
    ((4, 5), (1, 4), list()),
    # Test cases where A and B are both of size 3
    # In each test case, we shift A one further to the right
    ((0, 3), (3, 6), list()),
    ((1, 4), (3, 6), [(3, 4)]),
    ((2, 5), (3, 6), [(3, 5)]),
    ((3, 6), (3, 6), [(3, 6)]),
    ((4, 7), (3, 6), [(4, 6)]),
    ((5, 8), (3, 6), [(5, 6)]),
    ((6, 9), (3, 6), list()),
    ((7, 10), (3, 6), list()),
]

Range = Tuple[int, int]


@pytest.mark.parametrize("a, b, intersection", intersect)
def test_intersect_ranges(a: Range, b: Range, intersection: List[Range]) -> None:
    assert Bed._intersect(a, b) == intersection


to_range = [
    ([], []),
    ([0], [(0, 1)]),
    ([0, 1], [(0, 2)]),
    ([0, 1, 2, 4, 5], [(0, 3), (4, 6)]),
    ([3, 2, 1, 0], [(0, 4)]),
]


@pytest.mark.parametrize("numbers, expected", to_range)
def test_to_range(numbers: List[int], expected: List[Range]) -> None:
    # Test from numbers to ranges
    assert _to_range(set(numbers)) == expected

    # Test from returned ranges to numbers
    nums = list()
    for r in expected:
        nums += list(range(*r))
    assert set(nums) == set(numbers)


range_start_size = [
    # Range, offset, start, size
    ((0, 10), 0, 0, 10),
    ((10, 20), 10, 0, 10),
    ((10, 20), 5, 5, 10),
]


@pytest.mark.parametrize("range_, offset, start, size", range_start_size)
def test_range_to_blocks(range_: Range, offset: int, start: int, size: int) -> None:
    assert _range_to_start_size(range_, offset) == (start, size)
