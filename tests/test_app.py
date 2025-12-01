import pytest
from fastapi.testclient import TestClient

from gtgt import Bed
from gtgt.app import app
from gtgt.models import BedModel, TranscriptModel
from gtgt.transcript import Transcript


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
    coding_exons = BedModel(chrom="chr1", blocks=[(5, 6)], name="coding_exons")
    before = TranscriptModel(exons=exons, coding_exons=coding_exons)

    # We want to skip the second exon
    skip = BedModel(chrom="chr1", blocks=[(5, 6)])

    # After skipping the exon
    after_exons = BedModel(chrom="chr1", blocks=[(2, 4), (7, 11)])
    after_coding_exons = BedModel(chrom="chr1", blocks=[(5, 5)], name="coding_exons")
    after = TranscriptModel(exons=after_exons, coding_exons=after_coding_exons)

    # JSON cannot do tuples, so we have to make those into lists
    expected = after.model_dump()
    expected["exons"]["blocks"] = [list(range) for range in expected["exons"]["blocks"]]
    expected["coding_exons"]["blocks"] = [
        list(range) for range in expected["coding_exons"]["blocks"]
    ]

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
    self              - -   -   - - - -
    coding_exons        -
    other             - -       - - - -
    """
    # One transcript, which is smaller
    exons = Bed.from_blocks("chr1", [(2, 4), (7, 11)])
    exons.name = "exons"
    coding_exons = Bed("chr1", 3, 4, "coding_exons")
    self = Transcript(exons, coding_exons)

    # Other Transcript
    exons = Bed.from_blocks("chr1", [(2, 4), (5, 6), (7, 11)])
    exons.name = "exons"
    coding_exons = Bed("chr1", 3, 4, "coding_exons")
    other = Transcript(exons, coding_exons)

    expected = [
        {
            "name": "exons",
            "percentage": 6 / 7,
            "basepairs": "6/7",
        },
        {
            "name": "coding_exons",
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
