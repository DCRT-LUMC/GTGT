import logging
from typing import Any, Mapping

from .models import EnsemblTranscript
from .provider import Ensembl

logger = logging.getLogger(__name__)


def lookup_transcript(transcript_id: str) -> EnsemblTranscript:
    transcript, version = transcript_id.split(".")
    provider = Ensembl()
    ts = provider.get(transcript)

    _check_transcript(ts, int(version))

    return payload_to_ensemble_transcript(ts)


def payload_to_ensemble_transcript(payload: Mapping[str, Any]) -> EnsemblTranscript:
    return EnsemblTranscript(**payload)


def _check_transcript(payload: Mapping[str, Any], version: int) -> None:
    if (e_version := payload.get("version")) != version:
        msg = f"Ensembl returned a different version of the transcript: {e_version} instead of {version}"
        raise ValueError(msg)
    return
