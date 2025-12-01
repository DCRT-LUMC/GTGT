import copy
from typing import Any, Iterator, Sequence

from .range import Range, intersect, overlap, subtract

# colorRgb field from Bed
color = tuple[int, ...]


class Bed:
    def __init__(
        self,
        chrom: str | None = None,
        chromStart: int | None = None,
        chromEnd: int | None = None,
        name: str = ".",
        score: int = 0,
        strand: str = ".",
        thickStart: int | None = None,
        thickEnd: int | None = None,
        itemRgb: color = (0, 0, 0),
        blockCount: int | None = None,
        blockSizes: list[int] | None = None,
        blockStarts: list[int] | None = None,
        **ignored: Any,
    ) -> None:
        # Required attributes for .bed files
        self.chrom = chrom if chrom is not None else ""
        self.chromStart = chromStart if chromStart is not None else 0
        self.chromEnd = chromEnd if chromEnd is not None else 0

        # Simple attributes
        self.name = name
        self.score = score
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

        self.blockCount = blockCount if blockCount else len(self.blockSizes)

        self.validate()

    def validate(self) -> None:
        """Validate the internal constistence of the Bed record"""
        if self.thickStart < self.chromStart or self.thickStart > self.chromEnd:
            raise ValueError(f"thickStart outside of record ({self})")
        if self.thickEnd < self.chromStart or self.thickEnd > self.chromEnd:
            raise ValueError(f"thickEnd outside of record ({self})")
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
        blocks = self.blocks()[1:]
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

    def blocks(self) -> list[Range]:
        """Iterate over all blocks in the Bed record"""
        blocks = list()
        for size, start in zip(self.blockSizes, self.blockStarts):
            block_start = self.chromStart + start
            block_end = block_start + size
            blocks.append((block_start, block_end))
        return blocks

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

    def __bool__(self) -> bool:
        return self.size > 0

    def _zero_out(self) -> None:
        """Zero out the Bed object, by setting all ranges to the start"""
        self.chromEnd = self.thickStart = self.thickEnd = self.chromStart

        self.blockCount = 1
        self.blockSizes = self.blockStarts = [0]

    @staticmethod
    def _csv_to_int(csv: str) -> list[int]:
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

    @classmethod
    def from_blocks(cls, chrom: str, blocks: Sequence[Range]) -> "Bed":
        """Create a Bed record from multiple Ranges"""
        bed = cls(chrom, 0, 0)

        # Ensure the blocks are sorted in ascending order
        sorted_blocks = sorted(blocks, key=lambda x: x[0])
        bed.update(sorted_blocks)

        # Set the thickStart, thickEnd to the entire range
        bed.thickStart = bed.chromStart
        bed.thickEnd = bed.chromEnd

        return bed

    def intersect(self, other: object) -> "Bed":
        """Update record to only contain features that overlap other"""
        if not isinstance(other, Bed):
            raise NotImplementedError

        if self.strand != other.strand:
            raise ValueError("Conflicting strands, intersection not possible")

        new = copy.deepcopy(self)
        # If other is on a different chromosome, we zero out self since there
        # is no overlap
        if new.chrom != other.chrom:
            new._zero_out()
            return new

        # Determine all intersected ranges
        intersected: list[Range] = list()

        for range1 in new.blocks():
            for intersector in other.blocks():
                intersected += intersect(range1, intersector)

        new.update(intersected)
        return new

    def overlap(self, other: object) -> "Bed":
        """All blocks from self that (partially) overlap blocks from other"""
        if not isinstance(other, Bed):
            raise NotImplementedError

        new = copy.deepcopy(self)
        # If other is on a different chromosome, there can be no overlap
        if new.chrom != other.chrom:
            new._zero_out()
            return new

        # Calculating overlap on different strands is not supported
        if new.strand != other.strand:
            raise ValueError(
                "Calculating overlap on different strands is not supported"
            )

        blocks_to_keep = list()
        for block in new.blocks():
            for other_block in other.blocks():
                if overlap(block, other_block):
                    blocks_to_keep.append(block)
                    break  # Go to the next block once we find overlap
        new.update(blocks_to_keep)

        return new

    def subtract(self, other: object) -> None:
        """Subtract the blocks in other from blocks in self"""
        if not isinstance(other, Bed):
            raise NotImplementedError

        if self.chrom != other.chrom:
            return

        subtracted_blocks = subtract(self.blocks(), other.blocks())
        self.update(subtracted_blocks)

    def compare(self, other: object) -> float:
        """
        Compare the size of Bed objects

        Raises a ValueError if chrom, name or strand are not identical
        """

        def get_properties(bed: Bed) -> tuple[str, ...]:
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

    def compare_basepair(self, other: object) -> str:
        """
        Compare the size of Bed objects, return a string showing remainting basepairs

        Raises a ValueError if chrom, name or strand are not identical
        """

        def get_properties(bed: Bed) -> tuple[str, ...]:
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

        return f"{self.size}/{other.size}"

    def update(self, ranges: Sequence[Range]) -> None:
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


def _range_to_size_start(range: Range, offset: int) -> tuple[int, int]:
    """Convert a range to size, start

    BED format uses blockSizes and blockStarts to represent ranges
    """

    size = range[1] - range[0]
    start = range[0] - offset
    return size, start
