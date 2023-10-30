import pytest
from typing import List, Tuple, Any

from GTGT import Bed
from GTGT.bed import _to_range, _range_to_size_start, intersect, overlap


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


def test_defaul_values_blocks() -> None:
    bed = Bed("chr1", 5, 10)
    # Block positions are relative to chromStart, so the first block should
    # start at 0, not 5
    assert bed.blockStarts == [0]
    assert bed.blockSizes == [5]


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

    # Test setting the blocks. Note the trailing comma in the input of bed2
    bed1 = Bed("chr1", 0, 10, blockCount=3, blockStarts=[0, 3, 6], blockSizes=[1, 1, 4])
    bed2 = Bed(
        "chr1", "0", "10", blockCount="3", blockStarts="0,3,6,", blockSizes="1,1,4,"
    )
    assert bed1 == bed2

    # Automatically set blockCount based on the number of blocks
    bed1 = Bed("chr1", 0, 10, blockSizes=[1, 7], blockStarts=[0, 3])
    assert bed1.blockCount == 2


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
            blockStarts=[0, 7],
        ),
        "chr1	0	11	name	5	+	8	10	42,42,42	2	3,4	0,7",
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


intersections = [
    # Range A, range B, intersection
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


@pytest.mark.parametrize("a, b, intersection", intersections)
def test_intersect_ranges(a: Range, b: Range, intersection: List[Range]) -> None:
    assert intersect(a, b) == intersection
    assert intersect(b, a) == intersection


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
    # Range, offset, size, start
    ((0, 10), 0, 10, 0),
    ((10, 20), 10, 10, 0),
    ((10, 20), 5, 10, 5),
]


@pytest.mark.parametrize("range_, offset, size, start", range_start_size)
def test_range_to_blocks(range_: Range, offset: int, size: int, start: int) -> None:
    assert _range_to_size_start(range_, offset) == (size, start)


# Things that cannot be used to intersect a Bed record
not_intersectable = [
    (1, NotImplementedError),
    # If the strand is different, throw a value error
    (Bed("chr1", 0, 10, strand="+"), ValueError),
]


@pytest.mark.parametrize("intersector, error", not_intersectable)
def test_non_intersectable_for_bed(intersector: Any, error: Any, bed: Bed) -> None:
    """Test that we raise an error"""
    with pytest.raises(error):
        bed.intersect(intersector)


# fmt: off
intersect_bed = [
    # Bed before, Bed to intersect with, Bed after
    # If the intersector is on another chromosome, there is no overlap
    (Bed("chr1", 0, 10), Bed("chr2", 0, 10), Bed("chr1", 0, 0)),
    # If the intersector is the same, we should get the same record back
    (Bed("chr1", 0, 10,), Bed("chr1", 0, 10), Bed("chr1", 0, 10)),
    # If the intersector does not overlap at all, we should get an empty record back
    (Bed("chr1", 5, 10), Bed("chr1", 20, 30), Bed("chr1", 5, 5)),
    # The intersector does not overlap, and the record has multiple blocks
    (
        Bed("chr1", 5, 10, blockCount=2, blockSizes=[2, 2], blockStarts=[0, 3]),
        Bed("chr1", 20, 30),
        Bed("chr1", 5, 5)
    ),
    # The intersector with multiple blocks does not overlap
    (
        Bed("chr1", 20, 30),
        Bed("chr1", 5, 10, blockCount=2, blockSizes=[2, 2], blockStarts=[0, 3]),
        Bed("chr1", 20, 20)
    ),
    # Both have multiple blocks, and do not overlap
    # Both have multiple blocks, and one block overlaps
    # Both have multiple blocks, and each block partially overlaps
    #(
    #    Bed("chr1", 0, 10, blockCount=3, blockSizes=[4, 2, 1], blockStarts=[0, 5, 7]),
    #    Bed("chr1", 0, 10, blockCount=3, blockSizes=[3, 2, 2], blockStarts=[1, 5, 7]),
    #    Bed("chr1", 1, 9, blockCount=3, blockSizes=[3, 2, 1], blockStarts=[0, 4, 6])
    #)
]
# fmt: on


@pytest.mark.parametrize("before, intersector, after", intersect_bed)
def test_intersect_bed(before: Bed, intersector: Bed, after: Bed) -> None:
    before.intersect(intersector)
    assert before == after


def test_zero_bed_object() -> None:
    bed = Bed("chr1", 5, 10)
    bed._zero_out()
    assert bed.chromEnd == 5
    assert bed.thickStart == 5
    assert bed.thickEnd == 5
    assert bed.blockCount == 1
    assert bed.blockSizes == [0]
    assert bed.blockStarts == [0]


def test_update_bed_empty_ranges(bed: Bed) -> None:
    """
    Given an emtpy rangelist
    When we update a Bed record with the empty list
    The Bed record should be zeroed out
    """
    ranges: List[Range] = list()
    bed.update(ranges)

    # Test that the Bed record has been zeroed out
    assert bed.chromStart == 0
    assert bed.chromEnd == 0


def test_update_bed_with_ranges_thick() -> None:
    """
    If a record has thickStart or thickEnd set (not default)
    When we try to update the record with ranges
    Then we raise a NotImplementedError
    """
    record = Bed("chr1", 0, 10, thickStart=2, thickEnd=9)
    ranges: List[Range] = [(4, 8)]
    with pytest.raises(NotImplementedError):
        record.update(ranges)


def test_update_bed_with_ranges(bed: Bed) -> None:
    """Test updating a Bed record by providing a list of ranges"""
    ranges: List[Range] = [(1, 4), (5, 6), (7, 9)]
    bed.update(ranges)

    assert bed.chromStart == 1
    assert bed.chromEnd == 9
    assert bed.blockCount == 3
    assert bed.blockSizes == [3, 1, 2]
    assert bed.blockStarts == [0, 4, 6]


# fmt: off
invalid_bed = [
    # ThickStart before chromStart
    (
        ("chr1", 10, 20, ".", 0, "+", 1),
        "thickStart outside of record"
    ),
    # ThickStart after chromEnd
    (
        ("chr1", 10, 20, ".", 0, "+", 21),
        "thickStart outside of record"
    ),
    # ThickEnd before chromStart
    (
        ("chr1", 10, 20, ".", 0, "+", 10, 0),
        "thickEnd outside of record"
    ),
    # ThickEnd after chromEnd
    (
        ("chr1", 10, 20, ".", 0, "+", 10, 21),
        "thickEnd outside of record"
    ),
    # ThickEnd before thickStart
    (
        ("chr1", 10, 20, ".", 0, "+", 15, 10),
        "thickEnd before thickStart"
    ),
    # incorrect blockCount
    (
        ("chr1", 10, 20, ".", 0, "+", 10, 20, "0,0,0", 2),
        "blockCount does not match the number of blocks",
    ),
    # Mismatch in number of fields between blockSizes and blockStarts
    (
        ("chr1", 10, 20, ".", 0, "+", 10, 20, "0,0,0", 2, "2,3", "1,"),
        "number of values differs between blockSizes and blockStarts",
    ),
    # Blocks must not overlap
    (
        ("chr1", 10, 20, ".", 0, "+", 10, 20, "0,0,0", 2, "2,3", "0,1"),
        "Blocks must not overlap",
    ),
    # Block extends over the end of the Bed region
    (
        ("chr1", 10, 20, ".", 0, "+", 10, 20, "0,0,0", 1, "11", "0"),
        "Last block must end at self.chromEnd",
    ),
    # The first block must start at chromStart(=0, since this field is relative)
    (
        ("chr1", 10, 20, ".", 0, "+", 10, 20, "0,0,0", 1, "1", "9"),
        "The first block must start at chromStart",
    ),
    # The last block must end at chromEnd
    (
        ("chr1", 10, 20, ".", 0, "+", 10, 20, "0,0,0", 1, "8", "0"),
        "Last block must end at self.chromEnd",
    ),
    # Blocks must be in ascending order
    (
        ("chr1", 10, 20, ".", 0, "+", 10, 20, "0,0,0", 2, "1,1", "19,1"),
        "Blocks must be in ascending order",
    ),
]
# fmt: on


@pytest.mark.parametrize("bed, msg", invalid_bed)
def test_value_error(bed: str, msg: str) -> None:
    with pytest.raises(ValueError, match=msg):
        Bed(*bed)


def test_overlap_invalid_records() -> None:
    """Overlap is not supported for Bed records on different strands"""
    before = Bed("chr1", 0, 10, strand="+")
    selector = Bed("chr1", 0, 10, strand=".")
    with pytest.raises(ValueError):
        before.overlap(selector)


# fmt: off
bed_overlap = [
    # before, selector, after
    # Selector on different chromosome
    (Bed('chr1', 0, 10), Bed('chr2', 0, 10), Bed('chr1', 0, 0)),
    # Selector is after the record
    (Bed('chr1', 0, 10), Bed('chr1', 10, 20), Bed('chr1', 0, 0)),
    # Selector overlaps the start of the record by 1 bp
    (Bed('chr1', 10, 20), Bed('chr1', 0, 11), Bed('chr1', 10, 20)),
    # Selector overlaps the end of the record by 1 bp
    (Bed('chr1', 0, 10), Bed('chr1', 9, 20), Bed('chr1', 0, 10)),
    # Both the Record and Selector consist of multiple blocks
    (
        # Record
        Bed('chr1', 0, 10, blockSizes=[4, 2, 1], blockStarts=[0, 5, 9]),
        # Selector overlaps every block in Record, the first block twice
        Bed('chr1', 0, 15, blockSizes=[2, 3, 6], blockStarts=[0, 3, 9]),
        # Unchanged
        Bed('chr1', 0, 10, blockSizes=[4, 2, 1], blockStarts=[0, 5, 9]),
    )
]
# fmt: on


@pytest.mark.parametrize("before, selector, after", bed_overlap)
def test_bed_overlap(before: Bed, selector: Bed, after: Bed) -> None:
    before.overlap(selector)
    assert before == after


# fmt: off
range_overlap = [
    # A before B
    ((0, 3), (3, 6), False),
    # B before A
    ((3, 6), (0, 3), False),
    # A is the same as B
    ((0, 3), (0, 3), True),
    # A overlaps start B
    ((0, 3), (2, 6), True),
    # A is within B, start position is the same
    ((0, 3), (0, 6), True),
    # A is within B, end position is the same
    ((3, 6), (0, 6), True),
    # A is fully within B, end position is the same
    ((1, 3), (0, 6), True),
    # A overlaps the end of B
    ((4, 7), (0, 6), True),
    # B is within A
    ((0, 6), (1, 3), True),
    # B overlaps the start of A
    ((0, 6), (0, 3), True)
]
# fmt: on


@pytest.mark.parametrize("a, b, expected", range_overlap)
def test_range_overlap(a: Range, b: Range, expected: bool) -> None:
    assert overlap(a, b) == expected
    assert overlap(b, a) == expected
