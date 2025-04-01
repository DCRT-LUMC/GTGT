from typing import Tuple, List, Set

Range = Tuple[int, int]


def overlap(a: Range, b: Range) -> bool:
    """Determine if ranges a and b overlap"""
    # A and B overlap wether the intersection is not empty
    if intersect(a, b):
        return True
    return False


def intersect(a: Range, b: Range) -> List[Range]:
    """Determine the intersection between ranges a and b"""
    start = max(a[0], b[0])
    end = min(a[1], b[1])

    if start + 1 > end:
        return list()
    else:
        return [(start, end)]


def subtract(a: List[Range], b: List[Range]) -> List[Range]:
    """
    Subtract the regions in b from a

    Lazy implementation by just putting all numbers into two sets
    """
    A: Set[int] = set()
    for start, end in a:
        A.update(range(start, end))

    B: Set[int] = set()
    for start, end in b:
        B.update(range(start, end))

    return _to_range(A - B)


def _to_range(numbers: Set[int]) -> List[Range]:
    """Convert a set of numbers to a range[start, end)"""
    # Make sure the numbers are sorted
    _numbers = sorted(numbers)

    # If there are no _numbers
    if not _numbers:
        return list()

    # If there is only a single number
    if len(_numbers) == 1:
        i = _numbers[0]
        return [(i, i + 1)]

    # Initialise the start and previous number
    start = prev = _numbers[0]

    # Store the ranges we found
    ranges = list()

    # Process all _numbers
    for i in _numbers[1:]:
        if i == prev + 1:
            prev = i
        else:
            ranges.append((start, prev + 1))
            start = prev = i
    ranges.append((start, prev + 1))

    return ranges
