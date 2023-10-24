from typing import Optional, Iterator, List, Tuple, Union

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
            self.blockStarts = [self.chromStart]
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

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Bed):
            raise NotImplemented
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
