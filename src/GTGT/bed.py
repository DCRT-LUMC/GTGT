from typing import Optional, Iterator, List, Tuple


class Bed:
    def __init__(
        self,
        chrom: str,
        chromStart: int,
        chromEnd: int,
        name: str = ".",
        score: int = 0,
        strand: str = ".",
        thickStart: Optional[int] = None,
        thickEnd: Optional[int] = None,
        itemRgb: Tuple[int, int, int] = (0, 0, 0),
        blockCount: Optional[int] = None,
        blockSizes: Optional[List[int]] = None,
        blockStarts: Optional[List[int]] = None,
    ) -> None:
        # Required attributes
        self.chrom = chrom
        self.chromStart = chromStart
        self.chromEnd = chromEnd

        # Simple attributes
        self.name = name
        self.score = score
        self.strand = strand

        self.thickStart = thickStart if thickStart else self.chromStart
        self.thickEnd = thickEnd if thickEnd else self.chromEnd

        self.itemRgb = itemRgb

        # Set the blocks
        self.blockCount = blockCount if blockCount else 1
        if blockSizes is None:
            self.blockSizes = [self.chromEnd - self.chromStart]
        else:
            self.blockSizes = blockSizes
        self.blockStarts = blockStarts if blockStarts else [self.chromStart]

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
