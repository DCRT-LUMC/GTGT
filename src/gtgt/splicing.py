from dataclasses import dataclass

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
