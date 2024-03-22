from GTGT.models import BedModel
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
        "blockSizes": "1200,700,",
        "chromStarts": "0,1300,",
        "random_field": "some nonsense",
    }


def test_model_from_ucsc(ucsc: payload) -> None:
    """Test creating a BedModel from UCSC payload"""
    bm = BedModel.from_ucsc(ucsc)
    assert bm.chrom == 'chr1'
    # Test default value
    assert bm.itemRgb == (0, 0, 0)
    # Test that we rename "chromStart" in the payload to "blockStarts"
    assert bm.blockStarts == [0, 1300]
