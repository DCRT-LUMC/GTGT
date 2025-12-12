import pytest
from pydantic import ValidationError

from gtgt.bed import Bed
from gtgt.models import TranscriptId
from gtgt.mutalyzer import HGVS
from gtgt.transcript import Transcript

payload = dict[str, str | int]


@pytest.fixture
def ucsc() -> payload:
    return {
        "chrom": "chr1",
        "chromStart": 1000,
        "chromEnd": 2000,
        "name": "ENST00000.12",
        "score": 0,
        "strand": "-",
        "thickStart": 1000,
        "thickEnd": 2000,
        "blockCount": 2,
        "blockSizes": "200,700,",
        "chromStarts": "0,300,",
        "random_field": "some nonsense",
    }


def test_HGVS_model_valid() -> None:
    """
    GIVEN a valid HGVS description
    WHEN we make an HGVS object out of it
    THEN there should be no error
    """
    HGVS(description="NM_000094.4:c.5299G>C")


INVALID_HGVS = [
    "NM_000094.4:c.5299G>",
    "NM_000094.4>",
    "NM_000094",
]


@pytest.mark.parametrize("description", INVALID_HGVS)
def test_HGVS_model_invalid(description: str) -> None:
    """
    GIVEN an invalid HGVS description
    WHEN we make an HGVS object out of itemRgb
    THEN we should get a ValidationError
    """
    with pytest.raises(ValidationError):
        HGVS(description=description)


VALID_TRANSCRIPT_ID = [
    "ENST00000296930.10",
]


@pytest.mark.parametrize("id", VALID_TRANSCRIPT_ID)
def test_TranscriptId_valid(id: str) -> None:
    TranscriptId(id=id)


INVALID_TRANSCRIPT_ID = [
    "ENST00000296930",
    "ENST00000296930.10:c.100A>T",
]


@pytest.mark.parametrize("id", INVALID_TRANSCRIPT_ID)
def test_TranscriptId_invalid(id: str) -> None:
    with pytest.raises(ValidationError):
        TranscriptId(id=id)
