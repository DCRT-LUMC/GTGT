from typing import Optional, Iterator, List, Tuple, Union, Set
from .range import Range, overlap, intersect, subtract

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
        blockCount: Optional[castable_int] = None,
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

        self.blockCount = int(blockCount) if blockCount else len(self.blockSizes)

        self.validate()

    def validate(self) -> None:
        """Validate the internal constistence of the Bed record"""
        if self.thickStart < self.chromStart or self.thickStart > self.chromEnd:
            raise ValueError("thickStart outside of record")
        if self.thickEnd < self.chromStart or self.thickEnd > self.chromEnd:
            raise ValueError("thickEnd outside of record")
        if self.thickEnd < self.thickStart:
            raise ValueError("thickEnd before thickStart")
        if len(self.blockSizes) != self.blockCount:
            raise ValueError("blockCount does not match the number of blocks")
        if len(self.blockSizes) != len(self.blockStarts):
            raise ValueError(
                "number of values differs between blockSizes and blockStarts"
            )

        # Initialise with the end of the first block
        prev_end = self.chromStart + self.blockStarts[0] + self.blockSizes[0]
        prev_start = self.blockStarts[0]
        blocks = list(self.blocks())[1:]
        for start, end in blocks:
            if start < prev_start:
                raise ValueError("Blocks must be in ascending order")
            elif start < prev_end:
                raise ValueError("Blocks must not overlap")
            else:
                prev_end = end
                prev_start = start

        # The first block must start at chromStart
        if self.blockStarts[0] != 0:
            raise ValueError("The first block must start at chromStart")

        # The last block must end at chromEnd
        block_end = self.blockStarts[-1] + self.blockSizes[-1] + self.chromStart
        if block_end != self.chromEnd:
            raise ValueError("Last block must end at self.chromEnd")

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
            msg = f"Unsupported comparison between Bed and {type(other)}"
            raise NotImplementedError(msg)
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
        """Update record to only contain features that overlap other"""
        if not isinstance(other, Bed):
            raise NotImplementedError

        if self.strand != other.strand:
            raise ValueError("Conflicting strands, intersection not possible")

        # If other is on a different chromosome, we zero out self since there
        # is no overlap
        if self.chrom != other.chrom:
            self._zero_out()
            return

        # Determine all intersected ranges
        intersected: List[Range] = list()

        for range1 in self.blocks():
            for intersector in other.blocks():
                intersected += intersect(range1, intersector)

        self.update(intersected)

    def overlap(self, other: object) -> None:
        """All blocks from self that (partially) overlap blocks from other"""
        if not isinstance(other, Bed):
            raise NotImplementedError

        # If other is on a different chromosome, there can be no overlap
        if self.chrom != other.chrom:
            self._zero_out()

        # Calculating overlap on different strands is not supported
        if self.strand != other.strand:
            raise ValueError(
                "Calculating overlap on different strands is not supported"
            )

        blocks_to_keep = list()
        for block in self.blocks():
            for other_block in other.blocks():
                if overlap(block, other_block):
                    blocks_to_keep.append(block)
                    break  # Go to the next block once we find overlap
        self.update(blocks_to_keep)

    def subtract(self, other: object) -> None:
        """Subtract the blocks in other from blocks in self"""
        if not isinstance(other, Bed):
            raise NotImplementedError

        subtracted_blocks = subtract(list(self.blocks()), list(other.blocks()))
        self.update(subtracted_blocks)

    def update(self, ranges: List[Range]) -> None:
        """Update a Bed object with a list of ranges"""
        # Check that thickStart/thickEnd have not been set
        if self.thickStart != self.chromStart or self.thickEnd != self.chromEnd:
            raise NotImplementedError

        if not ranges:
            self._zero_out()
            return

        # The ranges are sorted
        self.chromStart = ranges[0][0]
        self.chromEnd = ranges[-1][-1]

        # Set to the new start/end of the record
        self.thickStart = self.chromStart
        self.thickEnd = self.chromEnd

        # Set the number of blocks
        self.blockCount = len(ranges)

        # Set the block starts and sizes
        self.blockSizes = list()
        self.blockStarts = list()

        for r in ranges:
            size, start = _range_to_size_start(range=r, offset=self.chromStart)
            self.blockSizes.append(size)
            self.blockStarts.append(start)
        self.validate()


def _range_to_size_start(range: Range, offset: int) -> Tuple[int, int]:
    """Convert a range to size, start

    BED format uses blockSizes and blockStarts to represent ranges
    """

    size = range[1] - range[0]
    start = range[0] - offset
    return size, start
