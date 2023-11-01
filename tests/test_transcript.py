import pytest

from GTGT import Bed
from GTGT.transcript import Transcript

from typing import List


@pytest.fixture
def records() -> List[Bed]:
    """
    Bed records that make up a transcript
    Each positions shown here is 10x
    (i) means inferred by the init method

              0 1 2 3 4 5 6 7 8 9
    exons     -   - -   -   - - -
    cds           - - - - - -
    coding(i)     - -   -   -


    """
    # fmt: off
    return [
        Bed(
            "chr1", 0, 100, name="exons",
            blockSizes=[10, 20, 10, 30],
            blockStarts=[0, 20, 50, 70],
        ),
        # The CDS is from (23, 72]
        Bed("chr1", 23, 72, name="cds", blockSizes=[49], blockStarts=[0]),
    ]
    # fmt: on


def test_transcript_init(records: List[Bed]) -> None:
    t = Transcript(records)
    assert t.exons.name == "exons"
    assert t.cds.name == "cds"


def test_coding(records: List[Bed]) -> None:
    t = Transcript(records)

    # The coding region is the intersection of the exons and the CDS
    coding = Bed(
        "chr1", 23, 72, name="coding", blockSizes=[17, 10, 2], blockStarts=[0, 27, 47]
    )
    # Test that we did not change the cds or exons
    assert t.exons == records[0]
    assert t.cds == records[1]

    # Test that coding was set
    assert t.coding == coding
