import logging
import urllib.request
from urllib.error import HTTPError
import json

from pydantic import BaseModel, Field
from enum import Enum

from typing import Any, Dict

logger = logging.getLogger(__name__)


class Assembly(Enum):
    HUMAN = "GRCh38"
    RAT = "mRatBN7.2"


class EnsemblTranscript(BaseModel):
    assembly_name: Assembly
    seq_region_name: str
    start: int
    end: int
    id: str
    version: int
    display_name: str


def lookup_transcript(transcript_id: str) -> EnsemblTranscript:
    transcript, version = transcript_id.split(".")
    ts = fetch_transcript(transcript)

    _check_transcript(ts, int(version))

    return payload_to_ensemble_transcript(ts)


def payload_to_ensemble_transcript(payload: Dict[str, Any]) -> EnsemblTranscript:
    return EnsemblTranscript(**payload)


def _check_transcript(payload: Dict[str, Any], version: int) -> None:
    if e_version := payload.get("version") != version:
        msg = f"Ensembl returned a different version of the transcript: {e_version} instead of {version}"
        raise ValueError(msg)
    return


def fetch_transcript(transcript: str) -> Dict[str, Any]:
    url = (
        f"http://rest.ensembl.org/lookup/id/{transcript}?content-type=application/json"
    )

    try:
        logger.debug(f"Fetching transcript, {url=}")
        response = urllib.request.urlopen(url)
    except HTTPError as e:
        raise RuntimeError(e)
    else:
        js: Dict[str, Any] = json.loads(response.read())

    return js