from dataclasses import dataclass
from typing import Optional, Union, Iterator


@dataclass
class Bed:
    chrom: str
    chromStart: int
    chromEnd: int
    # Simple attributes
    name: Optional[str] = "."
    score: Optional[int] = 0
    strand: Optional[str] = "."
    # Display attributes
    thickStart: Optional[int] = None
    thickEnd: Optional[int] = None
    itemRgb: Optional[Union[tuple[int, int, int], int]] = None
    # Blocks
    blockCount: Optional[int] = None
    blockSizes: Optional[list[int]] = None
    blockStarts: Optional[list[int]] = None

    def __post_init__(self) -> None:
        # Set thick start/end to chrom start/end
        if self.thickStart is None:
            self.thickStart = self.chromStart
        if self.thickEnd is None:
            self.thickEnd = self.chromEnd

        # Set the default colour to black.
        # According to the Bed standard, 0 is an alias for (0, 0, 0)
        if self.itemRgb is None or self.itemRgb == 0:
            self.itemRgb = (0, 0, 0)

        # Set the blocks
        if self.blockCount is None:
            self.blockCount = 1
        if self.blockSizes is None:
            self.blockSizes = [self.chromEnd - self.chromStart]
        if self.blockStarts is None:
            self.blockStarts = [self.chromStart]

    def blocks(self) -> Iterator[tuple[int, int]]:
        assert self.blockSizes is not None
        assert self.blockStarts is not None
        for size, start in zip(self.blockSizes, self.blockStarts):
            block_start = self.chromStart + start
            block_end = block_start + size
            yield (block_start, block_end)
