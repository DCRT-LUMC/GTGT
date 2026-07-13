from dataclasses import dataclass

from gtgt.range import overlap


@dataclass
class ExonTranscript:
    """Transcript defined using Exons"""

    exons: list[tuple[int, int]]

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
    junctions: list[tuple[int, int]]

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

    def _alternative_junctions(
        self, novel_splice: tuple[int, int]
    ) -> list[tuple[int, int]]:
        """Integrate the alternative splice junction"""
        # If the novel junction is not fully inside the transcript
        splice_start, splice_end = novel_splice
        if splice_start < self.start or splice_end > self.end:
            return self.junctions

        # The new set of junctions after integrating the novel splice junction
        new_junctions = list()

        # If the novel splice junction is before the first regular splice junction
        # or if it has overlap with the first splice junction
        if splice_start < self.junctions[0][0] or overlap(
            novel_splice, self.junctions[0]
        ):
            new_junctions.append(novel_splice)

        # Iterate over all other junctions
        for junction in self.junctions:
            if not overlap(junction, novel_splice):
                new_junctions.append(junction)

        # If the novel splice junction is after the last regular splice junction
        # Don't add it twice if the novel splice junction ALSO starts before the first junction
        # Also add it if the novel splice junction has overlap with the last splice junction
        if (
            splice_end > self.junctions[-1][1]
            and not splice_start < self.junctions[0][0]
            and not overlap(novel_splice, self.junctions[-1])
        ):
            new_junctions.append(novel_splice)

        return new_junctions
