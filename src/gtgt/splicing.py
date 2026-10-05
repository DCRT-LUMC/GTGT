from dataclasses import dataclass
from typing import Sequence, TypeAlias

from gtgt.range import after, overlap

Junction: TypeAlias = tuple[int, int]


@dataclass
class ExonTranscript:
    """Transcript defined using Exons"""

    exons: list[Junction]

    def to_junctions(self) -> "JunctionTranscript":
        """Determine transcript start, end and junctions from exons"""
        transcript_start = self.exons[0][0]
        transcript_end = self.exons[-1][1]
        junctions = list()

        junction_start = None
        for exon in self.exons:
            # If this is the first exon, the splice junction starts at the end of
            # the exon
            if junction_start is None:
                junction_start = exon[1]
            else:
                # If we have already seen an exon, the start of this exon is the
                # end of the splice junction
                junction_end = exon[0]

                # Store the junction
                junctions.append((junction_start, junction_end))

                # The start of the next splice junction is the end of this exon
                junction_start = exon[1]
        return JunctionTranscript(transcript_start, transcript_end, junctions)


@dataclass
class JunctionTranscript:
    """Transcript defined using a start, end and Junctions"""

    start: int
    end: int
    junctions: list[Junction]

    def to_exons(self) -> ExonTranscript:
        """Determine the exons from the transcript start, end and list of junctions"""
        exons = list()

        # If there is only a single exon
        if not self.junctions:
            return ExonTranscript([(self.start, self.end)])

        # The start of the first exon is the transcript_start
        exon_start = self.start
        for junction in self.junctions:
            exon_end = junction[0]
            exons.append((exon_start, exon_end))

            exon_start = junction[1]
        else:
            # The end of the last exon is the transcript_end
            exons.append((exon_start, self.end))

        return ExonTranscript(exons)

    def _alternative_junctions(self, novel_junction: Junction) -> list[Junction]:
        """Integrate the alternative splice junction"""

        start, end = novel_junction
        if not end > start:
            raise ValueError(
                f"Junction end has to be after the start: ({start}, {end})"
            )

        # If the novel junction is not fully inside the transcript
        for position in novel_junction:
            if position < self.start or position > self.end:
                return self.junctions

        # First, we remove all junction that overlap the novel splice junction
        new_junctions = [
            junction
            for junction in self.junctions
            if not overlap(junction, novel_junction)
        ]

        # Now, we find where to insert the novel splice junction
        for i, junction in enumerate(new_junctions):
            # As soon as we find a junction which is after the novel splice
            # site, insert it before
            if after(junction, novel_junction):
                new_junctions.insert(i, novel_junction)
                break
        # If we don't find any, insert the novel junction at the end
        else:
            new_junctions.append(novel_junction)

        return new_junctions


def _must_be_increasing(junctions: Sequence[Junction]) -> None:
    """Raise an error when we find an invalid list of junctions

    A list of junctions is invalid if:
        - The the end is not after the start for any junction
        - The junctions are not in order (this makes it easier to check for
          overlap)
        - The junctions overlap
    """

    # Keep track of the end of the previous junction
    prev_start = 0
    prev_end = 0

    for start, end in junctions:
        if not end > start:
            raise ValueError(
                f"Junction end has to be after the start: ({start}, {end})"
            )
        if not start > prev_end:
            raise ValueError(
                f"Junction ({start}, {end}) is not after the previous junction ({prev_start}, {prev_end})"
            )
        prev_start, prev_end = start, end
