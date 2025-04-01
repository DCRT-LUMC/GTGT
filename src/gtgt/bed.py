from typing import Any, Dict, Optional, Iterator, List, Tuple
from .range import Range, overlap, intersect, subtract

# colorRgb field from Bed
color = Tuple[int, ...]


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
        itemRgb: color = (0, 0, 0),
        blockCount: Optional[int] = None,
        blockSizes: Optional[List[int]] = None,
        blockStarts: Optional[List[int]] = None,
        **ignored: Any,
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
        else:
            self.thickStart = thickStart

        if thickEnd is None:
            self.thickEnd = self.chromEnd
        else:
            self.thickEnd = thickEnd

        self.itemRgb = itemRgb

        # Set the blocks
        if blockSizes is None:
            self.blockSizes = [self.chromEnd - self.chromStart]
        else:
            self.blockSizes = blockSizes

        if blockStarts is None:
            # blockStarts are relative to chromStart, and the first block must
            # start at 0
            self.blockStarts = [0]
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
            msg = f"blockCount({self.blockCount=}) does not match the number of blocks({self.blockSizes=})"
            raise ValueError(msg)
        if len(self.blockSizes) != len(self.blockStarts):
            msg = f"number of values differs between blockSizes({self.blockSizes}) and blockStarts({self.blockStarts})"
            raise ValueError(msg)

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
            raise ValueError(f"Last block({block_end=}) must end at {self.chromEnd=}")

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

    @staticmethod
    def _csv_to_int(csv: str) -> List[int]:
        """Convert a csv list to a list of integers"""
        integers = list()
        for num in csv.split(","):
            if num != "":
                integers.append(int(num))
        return integers

    @classmethod
    def from_bedfile(cls, line: str) -> "Bed":
        """Create a BED record from a line from a Bed file"""
        header = [
            "chrom",
            "chromStart",
            "chromEnd",
            "name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "itemRgb",
            "blockCount",
            "blockSizes",
            "blockStarts",
        ]

        # Parse the bed fields into a dict
        d = {k: v for k, v in zip(header, line.split("\t"))}

        return cls(
            chrom=d["chrom"],
            chromStart=int(d["chromStart"]),
            chromEnd=int(d["chromEnd"]),
            name=d.get("name", ""),
            score=int(d.get("score", 0)),
            strand=d.get("strand", "."),
            thickStart=int(d["thickStart"]) if "thickStart" in d else None,
            thickEnd=int(d["thickEnd"]) if "thickEnd" in d else None,
            itemRgb=(
                tuple(cls._csv_to_int(d["itemRgb"])) if "itemRgb" in d else (0, 0, 0)
            ),
            blockCount=int(d["blockCount"]) if "blockCount" in d else None,
            blockSizes=cls._csv_to_int(d["blockSizes"]) if "blockSizes" in d else None,
            blockStarts=(
                cls._csv_to_int(d["blockStarts"]) if "blockStarts" in d else None
            ),
        )

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

        if self.chrom != other.chrom:
            return

        subtracted_blocks = subtract(list(self.blocks()), list(other.blocks()))
        self.update(subtracted_blocks)

    def compare(self, other: object) -> float:
        """
        Compare the size of Bed objects

        Raises a ValueError if chrom, name or strand are not identical
        """

        def get_properties(bed: Bed) -> Tuple[str, ...]:
            """Return a tuple of the properties that must match"""
            return (bed.chrom, bed.name, bed.strand)

        if not isinstance(other, Bed):
            raise NotImplementedError

        # Other can not have size zero
        if other.size == 0:
            msg = f"Comparison not allowed with zero-size Bed record ({other=})"
            raise ValueError(msg)

        # We only compare Bed objects if their properties match
        prop_self = get_properties(self)
        prop_other = get_properties(other)

        if prop_self != prop_other:
            msg = "Comparison failed, properties mismatch: {prop_self} != {prop_other}"
            raise ValueError(msg)

        return self.size / other.size

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

    @property
    def size(self) -> int:
        """The size of a Bed record is the sum of its blocks"""
        return sum(self.blockSizes)


def _range_to_size_start(range: Range, offset: int) -> Tuple[int, int]:
    """Convert a range to size, start

    BED format uses blockSizes and blockStarts to represent ranges
    """

    size = range[1] - range[0]
    start = range[0] - offset
    return size, start


def make_bed(chrom: str, *blocks: Range) -> Bed:
    """Create a Bed record from multiple Ranges"""
    bed = Bed(chrom, 0, 0)
    bed.update(list(blocks))
    return bed
