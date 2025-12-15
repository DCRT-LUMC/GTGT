import logging
from collections import namedtuple
from typing import Any, Mapping

from mutalyzer.description import Description

from gtgt.bed import Bed
from gtgt.mutalyzer import (
    get_assembly_name,
    get_chrom_name,
    get_offset,
    get_transcript_name,
)

from .models import Assembly, EnsemblTranscript
from .provider import UCSC, MyGene, Provider, payload

logger = logging.getLogger(__name__)

ENSEMBL_TO_UCSC = {
    Assembly.HUMAN: "hg38",
    Assembly.RAT: "rn6",
}

# Supported tracks which contain protein features
PROTEIN_TRACKS: list[str] = []
# Supported tracks which contain RNA features
RNA_TRACKS: list[str] = ["knownGene"]


# Holder for the required parameters for UCSC query
Parameters = namedtuple("Parameters", ["genome", "chrom", "start", "end", "track"])


def _lookup_track_payload(
    d: Description, track: str, ucsc: Provider = UCSC()
) -> payload:
    """Lookup the track payload for the specified Description"""
    # Determine the assembly
    if get_assembly_name(d) != "GRCh38":
        raise ValueError(f"Assembly {get_assembly_name(d)} for not supported ({d})")
    else:
        # Assembly as used in UCSC
        genome = "hg38"

    # Is this track supported
    if track not in PROTEIN_TRACKS and track not in RNA_TRACKS:
        raise ValueError(f"Track {track} is not supported, please ask to have it added")

    # Determine the start and end of the transcript of interest
    offset = get_offset(d)
    exons = d.get_selector_model()["exon"]
    transcript_start = exons[0][0] + offset
    transcript_end = exons[-1][1] + offset

    # Next, determine the uniprot ID for the protein domain
    parameters = Parameters(
        genome=genome,
        chrom=get_chrom_name(d),
        start=transcript_start,
        end=transcript_end,
        track=track,
    )
    return ucsc.get(parameters)


def lookup_track(
    d: Description, track: str, ucsc: Provider = UCSC(), mygene: Provider = MyGene()
) -> payload:
    """Lookup the specified track on UCSC"""
    return _lookup_track_payload(d, track, ucsc)


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
