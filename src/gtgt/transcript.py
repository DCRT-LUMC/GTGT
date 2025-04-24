from copy import deepcopy
import dataclasses

from .mutalyzer import (
    CdotVariant,
    Therapy,
    mutation_to_cds_effect,
    HGVS,
    exonskip,
    _init_model,
)
from mutalyzer.description import Description
from .bed import Bed

from typing import List, Dict, Union

TranscriptComparison = Dict[str, Dict[str, Union[float, str]]]


@dataclasses.dataclass
class Comparison:
    name: str
    percentage: float
    basepairs: str


@dataclasses.dataclass
class Result:
    """To hold results for separate mutations of a transcript"""

    therapy: Therapy
    comparison: List[Comparison]

    def __gt__(self, other: "Result") -> bool:
        """Sort Result based on the sum of the percentage"""
        if not isinstance(other, Result):
            msg = f"Unsupported comparison between Bed and {type(other)}"
            raise NotImplementedError(msg)

        total_self = sum(c.percentage for c in self.comparison)
        total_other = sum(c.percentage for c in other.comparison)
        return total_self > total_other


class Transcript:
    def __init__(self, exons: Bed, cds: Bed):
        self.exons = exons
        self.cds = cds

        # Set the coding region
        coding = deepcopy(self.exons)
        coding.name = "Coding exons"
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

    def compare(self, other: object) -> List[Comparison]:
        """Compare the size of each record in the transcripts"""
        if not isinstance(other, Transcript):
            raise NotImplementedError

        # Compare each record that makes up self and other
        # The comparison will fail if the record.name does not match
        cmp = list()
        for record1, record2 in zip(self.records(), other.records()):
            percentage = record1.compare(record2)
            fraction = record1.compare_basepair(record2)
            C = Comparison(record1.name, percentage, fraction)
            cmp.append(C)

        return cmp

    def compare_score(self, other: object) -> float:
        """Compare the size of each records in the transcripts

        Returns the average value for all records
        """
        if not isinstance(other, Transcript):
            raise NotImplementedError
        cmp = self.compare(other)

        values = [x.percentage for x in cmp]
        return sum(values) / len(cmp)

    def mutate(self, d: Description, variants: CdotVariant) -> None:
        """Mutate the transcript based on the specified hgvs description"""
        # Determine the CDS intervals that are affected by the hgvs description
        deleted = Bed.from_blocks(self.cds.chrom, *mutation_to_cds_effect(d, variants))
        # Subtract that region from the annotations
        self.subtract(deleted)

    def analyze(self, hgvs: str) -> List[Result]:
        """Analyse the transcript based on the specified hgvs description

        Calculate the score for the Wildtype (1), the patient transcript and the exon skips
        """
        transcript_id = hgvs.split(":c.")[0]
        variants = CdotVariant(hgvs.split(":c.")[1])

        results = list()

        # The wildtype has a score of 100% by default
        wt = Therapy(
            name="Wildtype",
            hgvs=f"{transcript_id}:c.=",
            description="These are the annotations as defined on the reference. They are always 100% by definition.",
        )
        wildtype = Result(wt, self.compare(self))
        results.append(wildtype)

        # Initialize the wildtype description
        d = Description(f"{transcript_id}:c.=")
        _init_model(d)

        # Determine the score of the patient
        patient = deepcopy(self)
        patient.mutate(d, variants)

        p = Therapy(
            name="Input",
            hgvs=hgvs,
            description="The annotations based on the supplied input variants.",
        )
        results.append(Result(p, patient.compare(self)))

        # Determine the score of each exon skip
        for skip in exonskip(d):
            # Add deletion to the patient mutation
            desc = HGVS(description=hgvs)
            try:
                desc.apply_deletion(HGVS(description=skip.hgvs))
            except NotImplementedError as e:
                # TODO add logging
                continue

            # Update the therapy hgvs after applying the deletion
            skip.hgvs = desc.description

            # Get the c. variant
            exonskip_variant = CdotVariant(desc.description.split("c.")[1])

            # Apply the combination to the wildtype transcript
            therapy = deepcopy(self)

            try:
                therapy.mutate(d, exonskip_variant)
            # Splice site error from mutalyzer, no protein prediction
            except KeyError:
                continue
            results.append(Result(therapy=skip, comparison=therapy.compare(self)))

        # Sort the results
        wt_patient = results[:2]
        rest = sorted(results[2:], reverse=True)
        return wt_patient + rest
