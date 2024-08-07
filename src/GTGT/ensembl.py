import logging
import urllib.request
from urllib.error import HTTPError
import json

from .models import Assembly, EnsemblTranscript
from .provider import Provider

from typing import Any, Dict

logger = logging.getLogger(__name__)


def lookup_transcript(provider: Provider, transcript_id: str) -> EnsemblTranscript:
    transcript, version = transcript_id.split(".")
    url = (
        f"http://rest.ensembl.org/lookup/id/{transcript}?content-type=application/json"
    )
    ts = provider.get(url)

    _check_transcript(ts, int(version))

    return payload_to_ensemble_transcript(ts)


def payload_to_ensemble_transcript(payload: Dict[str, Any]) -> EnsemblTranscript:
    return EnsemblTranscript(**payload)


def _check_transcript(payload: Dict[str, Any], version: int) -> None:
    if (e_version := payload.get("version")) != version:
        msg = f"Ensembl returned a different version of the transcript: {e_version} instead of {version}"
        raise ValueError(msg)
    return
