from GTGT.models import BedModel, TranscriptModel
import pytest
from GTGT import Bed
from GTGT.bed import make_bed
from GTGT.transcript import Transcript
from GTGT.app import app

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
    exons = make_bed("chr1", (2, 4), (5, 6), (7, 11))
    cds = Bed("chr1", 5, 6)
    before = Transcript(exons, cds)

    # We want to skip the second exon
    skip = Bed("chr1", 5, 6)

    # After skipping the exon
    after_exons = make_bed("chr1", (2, 4), (7, 11))
    after_cds = Bed("chr1", 5, 5)
    after = Transcript(after_exons, after_cds)
    expected = TranscriptModel.from_transcript(after).model_dump()

    # Update the itemrgb, since json cannot hold tuples
    expected["exons"]["itemRgb"] = list(expected["exons"]["itemRgb"])
    expected["cds"]["itemRgb"] = list(expected["cds"]["itemRgb"])

    body = {
        "transcript": TranscriptModel.from_transcript(before).model_dump(),
        "region": BedModel.from_bed(skip).model_dump(),
    }
    response = client.post("/transcript/exonskip", json=body)

    assert response.status_code == 200
    assert response.json() == expected
