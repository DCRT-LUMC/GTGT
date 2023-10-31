import pytest
from typing import List

from GTGT.range import Range, intersect, overlap, subtract, _to_range

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


@pytest.mark.parametrize("a, b, intersection", intersections)
def test_intersect_ranges(a: Range, b: Range, intersection: List[Range]) -> None:
    assert intersect(a, b) == intersection
    assert intersect(b, a) == intersection


# fmt: off
range_subtract = [
    # The selector is before A
    #    0 1 2 3 4 5 6 7 8 9
    # A            - - - - -
    # S  - -  - - -
    # E            - - - - -
    (
        # A
        [(5, 10)],
        # Selector
        [(0, 5)],
        # Expected
        [(5, 10)]
    ),
    # The selector overlaps the start of A
    #    0 1 2 3 4 5 6 7 8 9
    # A            - - - - -
    # S        - - - -
    # E                - - -
    (
        [(5, 10)],
        [(3, 7)],
        [(7, 10)],
    ),
    # The selector is contained in A
    #    0 1 2 3 4 5 6 7 8 9
    # A            - - - - -
    # S              - - -
    # E             -      -
    (
        [(5, 10)],
        [(6, 9)],
        [(5, 6), (9, 10)],
    ),
    # The selector overlaps the end of A
    #    0 1 2 3 4 5 6 7 8 9
    # A  - - - - -
    # S        - - -
    # E  - - -
    (
        [(0, 5)],
        [(3, 6)],
        [(0, 3)],
    ),
    # The selector is after A
    #    0 1 2 3 4 5 6 7 8 9
    # A  - - - - -
    # S              - - -
    # E  - - -  - -
    (
        [(0, 5)],
        [(6, 9)],
        [(0, 5)],
    ),
    # Both A and the selector consist of multiple ranges
    # The first region in A has partial overlap with the first two selector regions
    # The second selector region (partially) overlaps the first two regions in A
    #
    #    0 1 2 3 4 5 6 7 8 9
    # A  - - - -   - -     -
    # S  - -   - - -       - - - - - -
    # E      -       -
    (
        [(0, 4), (5, 7), (9, 10)],
        [(0, 2), (3, 6), (9, 15)],
        [(2, 3), (6, 7)],
    ),
]
# fmt: on


@pytest.mark.parametrize("a, b, expected", range_subtract)
def test_subtract_ranges(a: List[Range], b: List[Range], expected: List[Range]) -> None:
    assert subtract(a, b) == expected


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
