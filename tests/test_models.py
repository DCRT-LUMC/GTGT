from GTGT.models import BedModel, TranscriptModel
from GTGT.bed import Bed
import pytest
from typing import Dict, Union

payload = Dict[str, Union[str, int]]


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
    assert bm.chrom == "chr1"
    # Test default value
    assert bm.itemRgb == (0, 0, 0)
    # Test that we rename "chromStart" in the payload to "blockStarts"
    assert bm.blockStarts == [0, 300]


def test_Bed_from_model(ucsc: payload) -> None:
    """Test creating a Bed object from BedModel"""
    bm = BedModel.from_ucsc(ucsc)

    bed = bm.to_bed()
    assert bed.chrom == "chr1"
    assert bed.itemRgb == (0, 0, 0)
    assert bed.blockStarts == [0, 300]


def test_Bed_validation(ucsc: payload) -> None:
    # The first block must start at position 0
    ucsc["chromStarts"] = "10,300"
    with pytest.raises(ValueError):
        BedModel.from_ucsc(ucsc)


def test_BedModel_from_bed() -> None:
    bed = Bed("chr1", 0, 10)
    bm = BedModel.from_bed(bed)
    assert bm.chrom == "chr1"
    assert bm.itemRgb == (0, 0, 0)
    assert bm.blockStarts == [0]


def test_transcript_model() -> None:
    bed = Bed("chr1", 0, 10)
    bm = BedModel.from_bed(bed)
    tm = TranscriptModel(exons=bm, cds=bm)
    # Test converting a TranscriptModel to a Transcript
    transcript = tm.to_transcript()
    # Test that the "coding" region has been set in the new Transcript
    assert transcript.coding == Bed("chr1", 0, 10, name="coding")
