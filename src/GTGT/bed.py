from typing import Optional, Iterator, List, Tuple, Union, Set

# Int, or a string we can cast to int
castable_int = Union[int, str]

# colorRgb field from Bed
color = Union[str, Tuple[int, int, int]]

# Either [1, 2, 3] or "1,2,3"
castable_list = Union[List[int], str]


class Bed:
    def __init__(
        self,
        chrom: str,
        chromStart: castable_int,
        chromEnd: castable_int,
        name: str = ".",
        score: castable_int = 0,
        strand: str = ".",
        thickStart: Optional[castable_int] = None,
        thickEnd: Optional[castable_int] = None,
        itemRgb: color = (0, 0, 0),
        blockCount: castable_int = 1,
        blockSizes: Optional[castable_list] = None,
        blockStarts: Optional[castable_list] = None,
    ) -> None:
        # Required attributes
        self.chrom = chrom
        self.chromStart = int(chromStart)
        self.chromEnd = int(chromEnd)

        # Simple attributes
        self.name = name
        self.score = int(score)
        self.strand = strand

        if thickStart is None:
            self.thickStart = self.chromStart
        elif isinstance(thickStart, str):
            self.thickStart = int(thickStart)
        else:
            self.thickStart = thickStart

        if thickEnd is None:
            self.thickEnd = self.chromEnd
        elif isinstance(thickEnd, str):
            self.thickEnd = int(thickEnd)
        else:
            self.thickEnd = thickEnd

        if isinstance(itemRgb, str):
            self.itemRgb = tuple(map(int, itemRgb.split(",")))
        else:
            self.itemRgb = itemRgb

        # Set the blocks
        self.blockCount = int(blockCount)

        if blockSizes is None:
            self.blockSizes = [self.chromEnd - self.chromStart]
        elif isinstance(blockSizes, str):
            self.blockSizes = list(map(int, (x for x in blockSizes.split(",") if x)))
        else:
            self.blockSizes = blockSizes

        if blockStarts is None:
            # blockStarts are relative to chromStart, and the first block must
            # start at 0
            self.blockStarts = [0]
        elif isinstance(blockStarts, str):
            self.blockStarts = list(map(int, (x for x in blockStarts.split(",") if x)))
        else:
            self.blockStarts = blockStarts

    def blocks(self) -> Iterator[Tuple[int, int]]:
        """Iterate over all blocks in the Bed record"""
        for size, start in zip(self.blockSizes, self.blockStarts):
            block_start = self.chromStart + start
            block_end = block_start + size
            yield (block_start, block_end)

    def __str__(self) -> str:
        return "\t".join(
            map(
                str,
                (
                    self.chrom,
                    self.chromStart,
                    self.chromEnd,
                    self.name,
                    self.score,
                    self.strand,
                    self.thickStart,
                    self.thickEnd,
                    ",".join(map(str, self.itemRgb)),
                    self.blockCount,
                    ",".join(map(str, self.blockSizes)),
                    ",".join(map(str, self.blockStarts)),
                ),
            )
        )

    def __repr__(self) -> str:
        return (
            f"Bed({self.chrom}, "
            f"{self.chromStart}, {self.chromEnd}, "
            f"name='{self.name}', "
            f"score={self.score}, "
            f"strand='{self.strand}', "
            f"thickStart='{self.thickStart}', "
            f"thickEnd='{self.thickEnd}', "
            f"blockCount='{self.blockCount}', "
            f"blockSizes='{self.blockSizes}', "
            f"blockStarts='{self.blockStarts}', "
            ")"
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Bed):
            raise NotImplementedError
        return all(
            (
                self.chrom == other.chrom,
                self.chromStart == other.chromStart,
                self.chromEnd == other.chromEnd,
                self.name == other.name,
                self.score == other.score,
                self.strand == other.strand,
                self.thickStart == other.thickStart,
                self.thickEnd == other.thickEnd,
                self.itemRgb == other.itemRgb,
                self.blockCount == other.blockCount,
                self.blockSizes == other.blockSizes,
                self.blockStarts == other.blockStarts,
            )
        )

    def _zero_out(self) -> None:
        """Zero out the Bed object, by setting all ranges to the start"""
        self.chromEnd = self.thickStart = self.thickEnd = self.chromStart

        self.blockCount = 1
        self.blockSizes = self.blockStarts = [0]

    def intersect(self, other: object) -> None:
        if not isinstance(other, Bed):
            raise NotImplementedError

        if self.strand != other.strand:
            raise ValueError("Conflicting strands, intersection not possible")

        # If other is on a different chromosome, we zero out self since there is no overlap
        if self.chrom != other.chrom:
            self._zero_out()
            return

        # Determine all intersected ranges
        intersected = list()
        for range1 in self.blocks():
            for intersector in other.blocks():
                intersected += intersect(range1, intersector)
        # If there is no overlap
        if not intersected:
            self._zero_out()
        print(intersected)


Range = Tuple[int, int]


def intersect(a: Range, b: Range) -> List[Range]:
    """Determine the intersection between two ranges

    This is lazy implementation, where we create all numbers in a range,
    and use python sets to find the intersections.

    """
    intersect = set(range(*a)).intersection(set(range(*b)))
    return _to_range(intersect)


def _to_range(numbers: Set[int]) -> List[Range]:
    """Convert a range of numbers to [start, end)"""
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


def _range_to_start_size(range: Range, offset: int) -> Tuple[int, int]:
    """Convert a range to start, size

    BED format uses, blockSizes and blockStarts to represent ranges
    """

    size = range[1] - range[0]
    start = range[0] - offset
    return size, start
