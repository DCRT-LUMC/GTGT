from pydantic import BaseModel
from enum import Enum


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
