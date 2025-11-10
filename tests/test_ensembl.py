from typing import Any

import pytest

from gtgt.ensembl import (
    _check_transcript,
    lookup_transcript,
    payload_to_ensemble_transcript,
)
from gtgt.models import Assembly, EnsemblTranscript

INVALID = [
    ({"version": 10}, 11, ValueError),
]


@pytest.mark.parametrize("payload, transcript_version, error", INVALID)
def test_invalid_payload(
    payload: dict[str, Any], transcript_version: int, error: Any
) -> None:
    with pytest.raises(error):
        _check_transcript(payload, transcript_version)


def test_payload_to_EnsembleTranscript() -> None:
    # Fake ensembl payload
    p = {
        "assembly_name": Assembly.HUMAN,
        "seq_region_name": "17",
        "start": 0,
        "end": 10,
        "version": 99,
        "id": "transcript",
        "display_name": "Gene name",
    }
    expected = EnsemblTranscript(
        assembly_name=Assembly.HUMAN,
        seq_region_name="17",
        start=0,
        end=10,
        version=99,
        id="transcript",
        display_name="Gene name",
    )

    assert payload_to_ensemble_transcript(p) == expected


def test_lookup_transcript() -> None:
    """Test the lookup transcript method, using a simple dict as provider"""
    p = {
        "assembly_name": Assembly.HUMAN,
        "seq_region_name": "17",
        "start": 0,
        "end": 10,
        "version": 99,
        "id": "transcript",
        "display_name": "Gene name",
    }
    transcript = "transcript"
    version = "99"
    url = (
        f"http://rest.ensembl.org/lookup/id/{transcript}?content-type=application/json"
    )
    provider = dict()
    provider[url] = p

    # Here, we use a Dict instead of a Provider, since both define a "get" method
    ES = lookup_transcript(provider, f"{transcript}.{version}")  # type: ignore
    assert ES.seq_region_name == "17"
