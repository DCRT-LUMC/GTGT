from .models import TranscriptModel, BedModel
from .provider import Provider
from .ensembl import lookup_transcript as lookup_transcript_ens
from .ucsc import lookup_knownGene


def lookup_transcript(provider: Provider, transcript_id: str) -> TranscriptModel:
    r = lookup_transcript_ens(provider, transcript_id)
    # track_name = "ncbiRefSeq"
    # track_name = "ensGene"
    track_name = "knownGene"
    track = lookup_knownGene(provider, r, track_name)
    knownGene = track[track_name][0]
    exons = BedModel.from_ucsc(knownGene)

    # The CDS is defied as the thickStart, thickEnd in ucsc
    chrom = knownGene["chrom"]
    start = knownGene["thickStart"]
    end = knownGene["thickEnd"]
    name = "cds"
    strand = knownGene["strand"]
    cds = BedModel(chrom=chrom, blocks=[(start, end)], name=name, strand=strand)

    return TranscriptModel(exons=exons, cds=cds)
