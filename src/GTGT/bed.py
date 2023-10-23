from dataclasses import dataclass
from typing import Optional, Union, Iterator, List, Tuple, cast


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
        assert self.blockSizes is not None
        assert self.blockStarts is not None
        for size, start in zip(self.blockSizes, self.blockStarts):
            block_start = self.chromStart + start
            block_end = block_start + size
            yield (block_start, block_end)

    def __str__(s) -> str:
        assert s.name is not None
        assert s.strand is not None
        return "\t".join(
            map(
                str,
                (
                    s.chrom,
                    s.chromStart,
                    s.chromEnd,
                    s.name,
                    s.score,
                    s.strand,
                    s.thickStart,
                    s.thickEnd,
                    ','.join(map(str, s.itemRgb)),
                    s.blockCount,
                    ','.join(map(str, s.blockSizes)),
                ),
            )
        )
