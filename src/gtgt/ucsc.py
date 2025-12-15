import logging
from collections import namedtuple
from typing import Any, Mapping

from .models import Assembly, EnsemblTranscript
from .provider import UCSC

logger = logging.getLogger(__name__)

ENSEMBL_TO_UCSC = {
    Assembly.HUMAN: "hg38",
    Assembly.RAT: "rn6",
}

Parameters = namedtuple("Parameters", ["genome", "chrom", "start", "end", "track"])


def chrom_to_uscs(seq_region_name: str) -> str:
    return "chrM" if seq_region_name == "MT" else f"chr{seq_region_name}"


def lookup_knownGene(
    transcript: EnsemblTranscript, track_name: str
) -> Mapping[str, Any]:
    provider = UCSC()
    parameters = Parameters(
        genome="hg38",
        chrom=chrom_to_uscs(transcript.seq_region_name),
        start=transcript.start,
        end=transcript.end,
        track=track_name,
    )
    track = provider.get(parameters)
    ts = f"{transcript.id}.{transcript.version}"
    track[track_name] = [
        entry for entry in track[track_name] if entry.get("name") == ts
    ]
    return track
