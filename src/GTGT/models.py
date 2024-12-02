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


class HGVS(BaseModel):
    description: str

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "description": "NM_002520.7:c.860_863dup",
                }
            ]
        }
    }

    @model_validator(mode="after")
    def hgvs_parser(self) -> "HGVS":
        """Parse the HGVS description with mutalyzer-hgvs-parser"""
        hgvs_error = (
            mutalyzer_hgvs_parser.exceptions.UnexpectedCharacter,
            mutalyzer_hgvs_parser.exceptions.UnexpectedEnd,
        )
        try:
            mutalyzer_hgvs_parser.to_model(self.description)
        except hgvs_error as e:
            raise ValueError(e)
        return self

    @staticmethod
    def _validate_for_apply_deletion(hgvs: "HGVS") -> None:
        """
        Raise a NotImplementedError if the hgvs is not supported
        """
        model = mutalyzer_hgvs_parser.to_model(hgvs.description)

        # There must be one variant
        if len(model["variants"]) != 1:
            raise NotImplementedError

        var = model["variants"][0]

        # The variant must be in the CDS
        if "outside_cds" in var["location"]:
            raise NotImplementedError

        # The variant must not be intronic
        if "offset" in var["location"]:
            raise NotImplementedError

    @property
    def position(self) -> Tuple[int, int]:
        """
        Return the position of a description as (start, end)

        These are just the .c position, so 1 based and inclusive
        """
        model = mutalyzer_hgvs_parser.to_model(self.description)
        assert len(model["variants"]) == 1

        var = model["variants"][0]
        if var["location"]["type"] == "point":
            p = var["location"]["position"]
            return p, p
        elif var["location"]["type"] == "range":
            s = var["location"]["start"]["position"]
            e = var["location"]["end"]["position"]
            return s, e
        else:
            raise NotImplementedError

    def apply_deletion(self, other: "HGVS") -> None:
        """
        Apply a deletion to the current variant

        If the deletion does not overlap, add them together
        If the deletion completely overlaps the variant, replace the variant

        If the deletion partially overlaps the variant, raise an error
        """
        # Perform all validations
        self._validate_for_apply_deletion(other)
        self._validate_for_apply_deletion(self)

        # other must be a deletion or an insertion_deletion
        o_model = mutalyzer_hgvs_parser.to_model(other.description)
        o_type = o_model["variants"][0]["type"]
        if o_type not in ["deletion", "insertion_deletion"]:
            raise NotImplementedError

        s_model = mutalyzer_hgvs_parser.to_model(self.description)

        # self and other must refer to the same reference ID
        s_id = s_model["reference"]["id"]
        o_id = o_model["reference"]["id"]
        if s_id != o_id:
            raise NotImplementedError

        # Get the c. positions for start and end
        s_start, s_end = self.position
        o_start, o_end = other.position

        # Get the variants in text format
        s_var = self.description.split("c.")[1]
        o_var = other.description.split("c.")[1]

        # self is before other
        if s_end < o_start:
            self.description = f"{s_id}:c.[{s_var};{o_var}]"
        # self is after other
        elif s_start > o_end:
            self.description = f"{s_id}:c.[{o_var};{s_var}]"
        # self is fully inside other
        elif s_start >= o_start and s_end <= o_end:
            # We overwrite self with other
            self.description = other.description
        # partial overlaps are not supported
        else:
            raise NotImplementedError


class TranscriptId(BaseModel):
    id: str = Field(pattern=r"^ENST\d+\.\d+$")

    model_config = {
        "json_schema_extra": {
            "examples": [
                {"id": "ENST00000296930.10"},
            ]
        }
    }
