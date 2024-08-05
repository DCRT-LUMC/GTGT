from pydantic import BaseModel, model_validator
from enum import Enum
from typing import Any, Dict, List, Tuple, Union
from .bed import Bed
from .transcript import Transcript


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
    """Pydantic wrapper around Bed class"""

    chrom: str
    chromStart: int
    chromEnd: int
    name: str
    score: int
    strand: str
    thickStart: int
    thickEnd: int
    itemRgb: Tuple[int, ...] = (0, 0, 0)
    blockCount: int
    blockSizes: List[int]
    blockStarts: List[int]

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
        fields["blockStarts"] = Bed._csv_to_int(fields["chromStarts"])
        fields["blockSizes"] = Bed._csv_to_int(fields["blockSizes"])
        return cls(**fields)

    @classmethod
    def from_bed(cls, bed: Bed) -> "BedModel":
        return cls(
            chrom=bed.chrom,
            chromStart=bed.chromStart,
            chromEnd=bed.chromEnd,
            name=bed.name,
            score=bed.score,
            strand=bed.strand,
            thickStart=bed.thickStart,
            thickEnd=bed.thickEnd,
            itemRgb=bed.itemRgb,
            blockCount=bed.blockCount,
            blockSizes=bed.blockSizes,
            blockStarts=bed.blockStarts,
        )

    def to_bed(self) -> Bed:
        """Convert BedModel to Bed"""
        return Bed(**self.model_dump())


class TranscriptModel(BaseModel):
    exons: BedModel
    cds: BedModel

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
