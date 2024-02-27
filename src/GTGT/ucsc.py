import logging
import urllib.request
from urllib.error import HTTPError
import json

from pydantic import BaseModel, Field
from enum import Enum

from typing import Any, Dict

logger = logging.getLogger(__name__)


class Assembly(Enum):
    GRCH38 = "GRCh38"


class EnsemblTranscript(BaseModel):
    assembly_name: Assembly
    chrom: str = Field(pattern=r"chr\d+")
    start: int
    end: int
    id: str
    version: int


def lookup_transcript(transcript_id: str) -> EnsemblTranscript:
    transcript, version = transcript_id.split(".")
    ts = fetch_transcript(transcript)

    _check_transcript(ts, int(version))

    return payload_to_ensemble_transcript(ts)


def payload_to_ensemble_transcript(payload: Dict[str, Any]) -> EnsemblTranscript:
    ts = payload

    # Gather the fields we need
    fields = {
        "assembly_name": ts.get("assembly_name"),
        "chrom": f"chr{ts.get('seq_region_name')}",
        "start": ts.get("start"),
        "end": ts.get("end"),
        "id": ts.get("id"),
        "version": ts.get("version"),
    }

    return EnsemblTranscript(**fields)  # type: ignore[arg-type]


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
