import logging
from typing import Any, Mapping

from .models import Assembly, EnsemblTranscript
from .provider import Provider

logger = logging.getLogger(__name__)

ENSEMBL_TO_UCSC = {
    Assembly.HUMAN: "hg38",
    Assembly.RAT: "rn6",
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


def lookup_knownGene(
    provider: Provider, transcript: EnsemblTranscript, track_name: str
) -> Mapping[str, Any]:
    url = ucsc_url(transcript, track_name)
    track = provider.get(url)
    ts = f"{transcript.id}.{transcript.version}"
    track[track_name] = [
        entry for entry in track[track_name] if entry.get("name") == ts
    ]
    return track
