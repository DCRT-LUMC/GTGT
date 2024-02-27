import pytest
from typing import Dict, Any
from GTGT.ucsc import EnsemblTranscript, Assembly
from GTGT.ucsc import _check_transcript, payload_to_ensemble_transcript


INVALID = [
    ({"version": 10}, 11, ValueError),
]


@pytest.mark.parametrize("payload, transcript_version, error", INVALID)
def test_invalid_payload(
    payload: Dict[str, Any], transcript_version: int, error: Any
) -> None:
    with pytest.raises(error):
        _check_transcript(payload, transcript_version)


def test_payload_to_EnsembleTranscript() -> None:
    # Fake ensembl payload
    p = {
        "assembly_name": Assembly.GRCH38,
        "seq_region_name": "17",
        "start": 0,
        "end": 10,
        "version": 99,
        "id": "transcript",
    }
    expected = EnsemblTranscript(
        assembly_name=Assembly.GRCH38,
        chrom="chr17",
        start=0,
        end=10,
        version=99,
        id="transcript",
    )

    assert payload_to_ensemble_transcript(p) == expected
