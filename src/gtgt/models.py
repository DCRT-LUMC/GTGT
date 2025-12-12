from enum import Enum

from pydantic import BaseModel, Field

Range = tuple[int, int]


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


class TranscriptId(BaseModel):
    id: str = Field(pattern=r"^ENST\d+\.\d+$")

    model_config = {
        "json_schema_extra": {
            "examples": [
                {"id": "ENST00000296930.10"},
            ]
        }
    }
