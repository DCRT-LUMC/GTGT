from .models import TranscriptModel, BedModel
from .provider import Provider
from .ensembl import lookup_transcript as lookup_transcript_ens
from .ucsc import lookup_knownGene
from . import Bed


def lookup_transcript(provider: Provider, transcript_id: str) -> TranscriptModel:
    r = lookup_transcript_ens(provider, transcript_id)
    # track_name = "ncbiRefSeq"
    # track_name = "ensGene"
    track_name = "knownGene"
    track = lookup_knownGene(provider, r, track_name)
    knownGene = track[track_name][0]
    bm = BedModel.from_ucsc(knownGene)

    exons = bm.copy()
    cds = BedModel.from_bed(Bed(bm.chrom, bm.thickStart, bm.thickEnd))
    return TranscriptModel(exons=exons, cds=cds)
