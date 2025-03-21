from pydantic import BaseModel, model_validator
from enum import Enum
from typing import Any, Dict, List, Tuple, Union
from .bed import Bed
from .transcript import Transcript

from pydantic import Field

import mutalyzer_hgvs_parser

Range = Tuple[int, int]


class Assembly(Enum):
    HUMAN = "GRCh38"
    RAT = "mRatBN7.2"


class EnsemblTranscript(BaseModel):
    assembly_name: Assembly
    seq_region_name: str
    start: int
    end: int
    id: str
    version: int
    display_name: str


class BedModel(BaseModel):
    """Pydantic wrapper for BED objects"""

    chrom: str
    blocks: List[Range]
    name: str = ""
    score: int = 0
    strand: str = "."

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "chrom": "chr1",
                    "blocks": [[30, 35]],
                    "name": "",
                }
            ]
        }
    }

    @model_validator(mode="after")
    def valid_Bed(self) -> "BedModel":
        """Validate that the data is consistent with Bed requirements"""
        try:
            self.to_bed()
        except ValueError as e:
            raise (e)
        return self

    @classmethod
    def from_ucsc(cls, ucsc: Dict[str, Union[str, int]]) -> "BedModel":
        fields: Dict[str, Any] = ucsc.copy()

        # Note that UCSC calls the "blockStarts" field "chromStarts"
        block_starts = Bed._csv_to_int(fields["chromStarts"])
        block_sizes = Bed._csv_to_int(fields["blockSizes"])

        # Determine the blocks from the payload
        chr_start = fields["chromStart"]
        blocks = list()
        for start, size in zip(block_starts, block_sizes):
            block_start = start + chr_start
            block_end = block_start + size
            blocks.append((block_start, block_end))

        return cls(
            chrom=fields["chrom"],
            blocks=blocks,
            name=fields["name"],
            score=fields["score"],
            strand=fields["strand"],
        )

    def to_bed(self) -> Bed:
        chrom_start = self.blocks[0][0]
        chrom_end = self.blocks[-1][1]
        blockCount = len(self.blocks)
        blockSizes = [end - start for start, end in self.blocks]
        blockStarts = [start - chrom_start for start, _ in self.blocks]

        return Bed(
            chrom=self.chrom,
            chromStart=chrom_start,
            chromEnd=chrom_end,
            name=self.name,
            score=self.score,
            strand=self.strand,
            blockCount=blockCount,
            blockSizes=blockSizes,
            blockStarts=blockStarts,
        )

    @classmethod
    def from_bed(cls, bed: Bed) -> "BedModel":
        return cls(
            chrom=bed.chrom,
            blocks=list(bed.blocks()),
            name=bed.name,
            score=bed.score,
            strand=bed.strand,
        )


class TranscriptModel(BaseModel):
    exons: BedModel
    cds: BedModel

    model_config = {
        "json_schema_extra": {
            "examples": [
                # Example transcript
                {
                    "exons": {
                        "chrom": "chr1",
                        "blocks": [[0, 10], [20, 40], [50, 60], [70, 100]],
                        "name": "exons",
                    },
                    "cds": {
                        "chrom": "chr1",
                        "blocks": [[23, 72]],
                        "name": "cds",
                    },
                },
                # Skipped exon 2 ([20, 40])
                {
                    "exons": {
                        "chrom": "chr1",
                        "blocks": [[0, 10], [50, 60], [70, 100]],
                        "name": "exons",
                        "score": 0,
                        "strand": ".",
                    },
                    "cds": {
                        "chrom": "chr1",
                        "blocks": [[40, 72]],
                        "name": "cds",
                        "score": 0,
                        "strand": ".",
                    },
                },
            ]
        }
    }

    def to_transcript(self) -> Transcript:
        exons = self.exons.to_bed()
        cds = self.cds.to_bed()
        return Transcript(exons=exons, cds=cds)

    @classmethod
    def from_transcript(cls, transcript: Transcript) -> "TranscriptModel":
        """Create a TranscriptModel from a Transcript object"""
        exons = BedModel.from_bed(transcript.exons)
        cds = BedModel.from_bed(transcript.cds)
        return cls(exons=exons, cds=cds)


class TranscriptId(BaseModel):
    id: str = Field(pattern=r"^ENST\d+\.\d+$")

    model_config = {
        "json_schema_extra": {
            "examples": [
                {"id": "ENST00000296930.10"},
            ]
        }
    }
