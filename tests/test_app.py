from gtgt.models import BedModel, TranscriptModel
import pytest
from gtgt import Bed
from gtgt.transcript import Comparison, Transcript
from gtgt.app import app

from fastapi import FastAPI
from fastapi.testclient import TestClient


@pytest.fixture
def client() -> TestClient:
    return TestClient(app)


def test_exonskip(client: TestClient) -> None:
    """
             0 1 2 3 4 5 6 7 8 9 10
    exons        - -   -   - - - -
    cds                -
    skip               -
    """
    # The input transcript
    exons = BedModel(chrom="chr1", blocks=[(2, 4), (5, 6), (7, 11)])
    cds = BedModel(chrom="chr1", blocks=[(5, 6)], name="cds")
    before = TranscriptModel(exons=exons, cds=cds)

    # We want to skip the second exon
    skip = BedModel(chrom="chr1", blocks=[(5, 6)])

    # After skipping the exon
    after_exons = BedModel(chrom="chr1", blocks=[(2, 4), (7, 11)])
    after_cds = BedModel(chrom="chr1", blocks=[(5, 5)], name="cds")
    after = TranscriptModel(exons=after_exons, cds=after_cds)

    # JSON cannot do tuples, so we have to make those into lists
    expected = after.model_dump()
    expected["exons"]["blocks"] = [list(range) for range in expected["exons"]["blocks"]]
    expected["cds"]["blocks"] = [list(range) for range in expected["cds"]["blocks"]]

    body = {
        "transcript": before.model_dump(),
        "region": skip.model_dump(),
    }

    response = client.post("/transcript/exonskip", json=body)

    assert response.status_code == 200
    assert response.json() == expected


def test_compare(client: TestClient) -> None:
    """
             0 1 2 3 4 5 6 7 8 9 10
    self         - -   -   - - - -
    cds            -
    other        - -       - - - -
    """
    # One transcript, which is smaller
    exons = Bed.from_blocks("chr1", (2, 4), (7, 11))
    exons.name = "exons"
    cds = Bed("chr1", 3, 4, "cds")
    self = Transcript(exons, cds)

    # Other Transcript
    exons = Bed.from_blocks("chr1", (2, 4), (5, 6), (7, 11))
    exons.name = "exons"
    cds = Bed("chr1", 3, 4, "cds")
    other = Transcript(exons, cds)

    expected = [
        {
            "name": "exons",
            "percentage": 6 / 7,
            "basepairs": "6/7",
        },
        {
            "name": "cds",
            "percentage": 1.0,
            "basepairs": "1/1",
        },
        {
            "name": "Coding exons",
            "percentage": 1.0,
            "basepairs": "1/1",
        },
    ]

    body = {
        "self": TranscriptModel.from_transcript(self).model_dump(),
        "other": TranscriptModel.from_transcript(other).model_dump(),
    }
    response = client.post("/transcript/compare", json=body)

    assert response.status_code == 200
    assert response.json() == expected
