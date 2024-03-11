from .ensembl import Assembly, EnsemblTranscript
import logging
import urllib.request
from urllib.error import HTTPError
from typing import Any, Dict
import json

logger = logging.getLogger(__name__)

ENSEMBL_TO_UCSC = {
    Assembly.HUMAN: "hg38",
    Assembly.RAT: "rn7",
}


def chrom_to_uscs(seq_region_name: str) -> str:
    return "chrM" if seq_region_name == "MT" else f"chr{seq_region_name}"


def ucsc_url(transcript: EnsemblTranscript, track: str = "knownGene") -> str:
    genome = ENSEMBL_TO_UCSC[transcript.assembly_name]
    url = ";".join(
        (
            f"https://api.genome.ucsc.edu/getData/track?genome={genome}",
            f"chrom={chrom_to_uscs(transcript.seq_region_name)}",
            f"track={track}",
            f"start={transcript.start}",
            f"end={transcript.end}",
        )
    )

    return url


def fetch_transcript(
    transcript: EnsemblTranscript, track: str = "knownGene"
) -> Dict[str, Any]:
    url = ucsc_url(transcript, track)
    try:
        logger.debug(f"Fetching {track}: {url=}")
        response = urllib.request.urlopen(url)
    except HTTPError as e:
        raise RuntimeError(e)
    else:
        js: Dict[str, Any] = json.loads(response.read())

    return js


def lookup_knownGene(transcript: EnsemblTranscript) -> Dict[str, Any]:
    track = fetch_transcript(transcript, "knownGene")
    ts = f"{transcript.id}.{transcript.version}"
    track["knownGene"] = [
        entry for entry in track["knownGene"] if entry.get("name") == ts
    ]
    return track
