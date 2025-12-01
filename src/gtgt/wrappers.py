from .ensembl import lookup_transcript as lookup_transcript_ens
from .models import BedModel, TranscriptModel
from .provider import Provider
from .ucsc import lookup_knownGene


def lookup_transcript(provider: Provider, transcript_id: str) -> TranscriptModel:
    r = lookup_transcript_ens(provider, transcript_id)
    # track_name = "ncbiRefSeq"
    # track_name = "ensGene"
    track_name = "knownGene"
    track = lookup_knownGene(provider, r, track_name)
    knownGene = track[track_name][0]
    exons = BedModel.from_ucsc(knownGene)

    # Rename the exons track to "Exons"
    exons.name = "Exons"

    # The CDS is defied as the thickStart, thickEnd in ucsc
    chrom = knownGene["chrom"]
    start = knownGene["thickStart"]
    end = knownGene["thickEnd"]
    name = "CDS"
    strand = knownGene["strand"]
    cds = BedModel(
        chrom=chrom, blocks=[(start, end)], name=name, strand=strand
    ).to_bed()

    # Determine the coding region
    coding_exons = exons.to_bed()
    coding_exons.name = "Coding exons"
    coding_exons = coding_exons.intersect(cds)

    return TranscriptModel(exons=exons, coding_exons=BedModel.from_bed(coding_exons))
