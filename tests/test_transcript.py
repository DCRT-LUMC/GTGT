import pytest

from GTGT import Bed
from GTGT.transcript import Transcript

from typing import List


@pytest.fixture
def records() -> List[Bed]:
    """
    Bed records that make up a transcript
    Each positions shown here is 10x

              0 1 2 3 4 5 6 7 8 9
    exons     -   - -   -   - - -
    cds           - - - - - -


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
