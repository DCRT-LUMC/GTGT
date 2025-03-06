from copy import deepcopy

from GTGT.mutalyzer import mutation_to_cds_effect, HGVS, exonskip

from .bed import Bed

from typing import List, Dict

TranscriptComparison = Dict[str, float]


class Transcript:
    def __init__(self, exons: Bed, cds: Bed):
        self.exons = exons
        self.cds = cds

        # Set the coding region
        coding = deepcopy(self.exons)
        coding.name = "coding"
        coding.intersect(self.cds)
        self.coding = coding

    def records(self) -> List[Bed]:
        """Return the Bed records that make up the Transcript"""
        return [self.exons, self.cds, self.coding]

    def intersect(self, selector: Bed) -> None:
        """Update transcript to only contain features that intersect the selector"""
        for record in self.records():
            record.intersect(selector)

    def overlap(self, selector: Bed) -> None:
        """Update transcript to only contain features that overlap the selector"""
        for record in self.records():
            record.overlap(selector)

    def subtract(self, selector: Bed) -> None:
        """Remove all features from transcript that intersect the selector"""
        for record in self.records():
            record.subtract(selector)

    def exon_skip(self, selector: Bed) -> None:
        """Remove the exon(s) that overlap the selector from the transcript"""
        exons_to_skip = deepcopy(self.exons)
        exons_to_skip.overlap(selector)
        self.subtract(exons_to_skip)

    def compare(self, other: object) -> Dict[str, float]:
        """Compare the size of each record in the transcripts"""
        if not isinstance(other, Transcript):
            raise NotImplementedError

        # Compare each record that makes up self and other
        # The comparison will fail if the record.name does not match
        cmp = dict()
        for record1, record2 in zip(self.records(), other.records()):
            cmp[record1.name] = record1.compare(record2)

        return cmp

    def compare_score(self, other: object) -> float:
        """Compare the size of each records in the transcripts

        Returns the average value for all records
        """
        if not isinstance(other, Transcript):
            raise NotImplementedError
        cmp = self.compare(other)

        return sum(cmp.values()) / len(cmp)

    def mutate(self, hgvs: str) -> None:
        """Mutate the transcript based on the specified hgvs description"""
        # Determine the CDS interval that is affected by the hgvs description
        H = HGVS(description=hgvs)
        chromStart, chromEnd = mutation_to_cds_effect(H)
        # Subtract that region from the annotations
        self.subtract(
            Bed(chrom=self.cds.chrom, chromStart=chromStart, chromEnd=chromEnd)
        )

    def analyze(self, hgvs: str) -> Dict[str, TranscriptComparison]:
        """Analyse the transcript based on the specified hgvs description

        Calculate the score for the Wildtype (1), the patient transcript and the exon skips
        """
        # Initialize the results dictionary. Wildtype has a score of 1 by definition
        results = dict()
        results["wildtype"] = self.compare(self)

        # Determine the score of the patient
        patient = deepcopy(self)
        patient.mutate(hgvs)
        results["patient"] = patient.compare(self)

        # Determine the score of each exon skip
        for skip in exonskip(HGVS(description=hgvs)):
            # Add deletion to the patient mutation
            desc = HGVS(description=hgvs)
            desc.apply_deletion(skip)

            # Apply the combination to the wildtype transcript
            therapy = deepcopy(self)
            try:
                therapy.mutate(desc.description)
            # Splice site error from mutalyzer, no protein prediction
            except KeyError:
                continue
            results[skip.description] = therapy.compare(self)

        return results
