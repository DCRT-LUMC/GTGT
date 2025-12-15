import logging
from collections import defaultdict, namedtuple
from typing import Any, Mapping, Sequence

from mutalyzer.description import Description

from gtgt.bed import Bed
from gtgt.mutalyzer import (
    get_assembly_name,
    get_chrom_name,
    get_offset,
    get_transcript_name,
)
from gtgt.variant_validator import parse_payload

from .models import Assembly, EnsemblTranscript
from .provider import UCSC, MyGene, Provider, VariantValidator, payload

logger = logging.getLogger(__name__)

ENSEMBL_TO_UCSC = {
    Assembly.HUMAN: "hg38",
    Assembly.RAT: "rn6",
}

# Supported tracks which contain protein features
PROTEIN_TRACKS: list[str] = ["unipDomain"]
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


def _tracks_to_bed(tracks: Sequence[payload]) -> list[Bed]:
    """Convert a list of json track objects from UCSC into a list of Bed records

    Note that trakcs with the same name will be merged into a single Bed record
    """
    # Store the resulting Bed records
    records: list[Bed] = list()

    # Ensure each track only contains a single block
    for track in tracks:
        if track["blockCount"] > 1:
            raise NotImplementedError(tracks)
    # Next, group the tracks by name
    grouped = defaultdict(list)
    for track in tracks:
        grouped[track["name"]].append(track)

    for tracks in grouped.values():
        if not tracks:
            raise RuntimeError(tracks)
        # Look in the first record for all constant values
        track = tracks[0]
        # Get the blocks
        blocks: list[tuple[int, int]] = list()
        for track in tracks:
            start = track["chromStart"]
            end = track["chromEnd"]
            blocks.append((start, end))
        # Get the chromosome
        chrom = track["chrom"]

        bed = Bed.from_blocks(chrom, blocks)
        # The color is stored in the 'reserved' field by UCSC
        bed.itemRgb = tuple(map(int, track["reserved"].split(",")))
        bed.score = track["score"]
        bed.strand = track["strand"]
        bed.name = track["name"]

        records.append(bed)

    return records


def lookup_track(
    d: Description,
    track: str,
    ucsc: Provider = UCSC(),
    mygene: Provider = MyGene(),
    variantvalidator: Provider = VariantValidator(),
) -> list[Bed]:
    """Lookup the specified track on UCSC"""
    # Required to match protein features to the correct transcript
    uniprot_id = None

    # For protein tracks, we have to lookup the uniprot_id
    if track in PROTEIN_TRACKS:
        # First, we look up the ENSG using VariantValidator
        variant = str(d)
        vv = variantvalidator.get(("hg38", variant))
        ensg = parse_payload(vv, variant, "hg38")["ensembl_gene_id"]
        assert ensg is not None

        # Look up protein ID to match the correct transcript
        uniprot_id = mygene.get((ensg,))["uniprot"]["Swiss-Prot"]
    elif track in RNA_TRACKS:
        pass
    else:
        raise ValueError(f"Track {track} is not supported, please ask to have it added")

    payload = _lookup_track_payload(d, track, ucsc)
    # Get only the tracks that match the specified uniprot ID
    if uniprot_id:
        tracks = [
            track for track in payload[track] if track.get("uniProtId") == uniprot_id
        ]
    else:
        tracks = [track for track in payload[track]]
    return _tracks_to_bed(tracks)


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
