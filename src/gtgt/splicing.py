from collections import namedtuple
from dataclasses import dataclass

Exon = namedtuple("Exon", "start, end")
Junction = namedtuple("Junction", "start, end")


def exons_to_junctions(exons: list[Exon]) -> tuple[int, int, list[Junction]]:
    """Determine transcript start, end and junctions from exons"""
    transcript_start = exons[0][0]
    transcript_end = exons[-1][1]
    junctions = list()

    junction_start = None
    for exon in exons:
        # If this is the first exon, the splice junction starts at the end of
        # the exon
        if junction_start is None:
            junction_start = exon[1]
        else:
            # If we have already seen an exon, the start of this exon is the
            # end of the splice junction
            junction_end = exon[0]

            # Store the junction
            junction = Junction(junction_start, junction_end)
            junctions.append(junction)

            # The start of the next splice junction is the end of this exon
            junction_start = exon[1]
    return transcript_start, transcript_end, junctions


def junctions_to_exons(transcript_start: int, transcript_end: int, junctions: list[Junction]) -> list[Exon]:
    """Determine the exons from the transcript start, end and list of junctions"""
    exons: list[Exon] = list()

    # If there is only a single exon
    if not junctions:
        exon = Exon(transcript_start, transcript_end)
        exons.append(exon)
        return exons

    # The start of the first exon is the transcript_start 
    exon_start = transcript_start
    for junction in junctions:
        exon_end = junction[0]
        exon = Exon(exon_start, exon_end)
        exons.append(exon)

        exon_start = junction[1]
    else:
        # The end of the last exon is the transcript_end
        exon = Exon(exon_start, transcript_end)
        exons.append(exon)

    return exons


