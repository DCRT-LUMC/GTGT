import pytest
from pydantic import ValidationError

from gtgt.bed import Bed
from gtgt.models import BedModel, TranscriptId, TranscriptModel
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


def test_model_from_ucsc(ucsc: payload) -> None:
    """Test creating a BedModel from UCSC payload"""
    bm = BedModel.from_ucsc(ucsc)

    expected = Bed(
        "chr1",
        1000,
        2000,
        "ENST00000.12",
        0,
        "-",
        thickStart=1000,
        thickEnd=2000,
        blockCount=2,
        blockSizes=[200, 700],
        blockStarts=[0, 300],
    )
    new_bed = bm.to_bed()

    assert new_bed == expected


def test_Bed_from_model(ucsc: payload) -> None:
    """Test creating a Bed object from BedModel"""
    bm = BedModel.from_ucsc(ucsc)

    bed = bm.to_bed()
    assert bed.chrom == "chr1"
    assert bed.itemRgb == (0, 0, 0)
    assert bed.blockStarts == [0, 300]


def test_Bed_validation(ucsc: payload) -> None:
    # Blocks must be in ascending order
    ucsc["chromStarts"] = "300,0"
    ucsc["blockSizes"] = "700,200"
    with pytest.raises(ValueError):
        BedModel.from_ucsc(ucsc)


def test_BedModel_from_bed() -> None:
    bed = Bed("chr1", 0, 10)
    bm = BedModel.from_bed(bed)
    assert bm.chrom == "chr1"
    assert bm.blocks == [(0, 10)]


def test_transcript_model() -> None:
    bed = Bed("chr1", 0, 10)
    bm = BedModel.from_bed(bed)
    tm = TranscriptModel(exons=bm, coding_exons=bm)
    # Test converting a TranscriptModel to a Transcript
    transcript = tm.to_transcript()
    # Test that the "coding" region has been set in the new Transcript
    assert transcript.coding_exons == Bed("chr1", 0, 10)


def test_Transcript_from_model() -> None:
    """
    GIVEN a Transcript
    WHEN we create a TranscriptModel out of it
    THEN it must match the expected TranscriptModel
    """
    # Create the TranscriptModel
    bed = Bed("chr1", 0, 10)
    bm = BedModel.from_bed(bed)
    expected = TranscriptModel(exons=bm, coding_exons=bm)

    # Create the Transcript we want to convert to a TranscriptModel
    bed = Bed("chr1", 0, 10)
    ts = Transcript(exons=bed, coding_exons=bed)

    assert TranscriptModel.from_transcript(ts) == expected


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
