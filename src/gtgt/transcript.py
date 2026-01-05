import dataclasses
import logging
from copy import deepcopy
from typing import Any, Mapping, Sequence

from mutalyzer.description import Description

from .bed import Bed
from .exonviz import draw
from .mutalyzer import (
    Therapy,
    Variant,
    generate_therapies,
    get_chrom_name,
    get_exons,
    get_offset,
    get_strand,
    init_description,
    mutation_to_cds_effect,
    protein_prediction,
    sequence_from_description,
)
from .ucsc import PROTEIN_TRACKS, lookup_track

logger = logging.getLogger(__name__)


@dataclasses.dataclass
class Comparison:
    name: str
    percentage: float
    basepairs: str

    @classmethod
    def from_dict(cls, dict: Mapping[str, Any]) -> "Comparison":
        """Create a Variant object from a dict representation of a Variant"""
        return cls(**dict)


@dataclasses.dataclass
class Result:
    """To hold results for separate mutations of a transcript"""

    therapy: Therapy
    comparison: Sequence[Comparison]

    def __gt__(self, other: "Result") -> bool:
        """Sort Result based on the sum of the percentage"""
        if not isinstance(other, Result):
            msg = f"Unsupported comparison between Bed and {type(other)}"
            raise NotImplementedError(msg)

        total_self = sum(c.percentage for c in self.comparison)
        total_other = sum(c.percentage for c in other.comparison)
        return total_self > total_other

    @classmethod
    def from_dict(cls, dict: Mapping[str, Any]) -> "Result":
        return cls(
            therapy=Therapy.from_dict(dict["therapy"]),
            comparison=[Comparison.from_dict(x) for x in dict["comparison"]],
        )


class Transcript:
    def __init__(self, rna_features: list[Bed], protein_features: list[Bed]):
        self.rna_features = rna_features
        self.protein_features = protein_features

    @classmethod
    def from_description(cls, d: Description) -> "Transcript":
        """Create a Transcript object from a mutalyzer Description"""
        # Check if we can use this Description to initialize a Transcript
        offset = get_offset(d)

        # Get exons and add the offset
        selector_model = d.get_selector_model()
        exons = selector_model["exon"]
        exons = [(start + offset, end + offset) for start, end in exons]

        # Get CDS and add the offset
        cds = selector_model["cds"][0]
        cds = (cds[0] + offset, cds[1] + offset)

        # Get the strand and chromosome name
        chrom = get_chrom_name(d)
        strand = get_strand(d)

        # Convert to Bed to intersect the exons with the CDS to get the coding
        # exons
        exon_bed = Bed.from_blocks(chrom, exons)
        exon_bed.name = "Exons"
        exon_bed.strand = strand

        cds_bed = Bed.from_blocks(chrom, [cds])
        cds_bed.strand = strand

        # Make a copy to intersect
        coding_exons = exon_bed.intersect(cds_bed)
        coding_exons.name = "Coding exons"

        return cls(rna_features=[exon_bed], protein_features=[coding_exons])

    def records(self) -> Sequence[Bed]:
        """Return the Bed records that make up the Transcript"""
        return self.rna_features + self.protein_features

    def rna_records(self) -> list[Bed]:
        """Return the Bed records that contain RNA features"""
        return self.rna_features

    def protein_records(self) -> list[Bed]:
        """Return the Bed records that contain protein features"""
        return self.protein_features

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

    def compare(self, other: object) -> Sequence[Comparison]:
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

    def mutate(self, d: Description, variants: Sequence[Variant]) -> None:
        """Mutate the transcript based on the specified variants"""
        # Determine the chromosome the transcript is on
        if self.records():
            chrom = self.records()[0].chrom
        else:
            chrom = ""

        # Update protein features
        protein_changes = Bed.from_blocks(chrom, mutation_to_cds_effect(d, variants))
        for record in self.protein_records():
            record.subtract(protein_changes)

        # Update RNA features
        rna_changes = Bed.from_blocks(
            chrom, [v.genomic_coordinates(d) for v in variants]
        )
        for record in self.rna_records():
            self.subtract(rna_changes)

    def analyze(self, hgvs: str) -> Sequence[Result]:
        """Analyze the transcript based on the specified HGVS description"""

        # r. notations are not supported
        coordinate_system = hgvs.split(":")[1][0:1]
        if coordinate_system != "c":
            raise NotImplementedError(
                f"Coordinate system '{coordinate_system}' is not supported"
            )
        # Initialize the input HGVS description
        d = init_description(hgvs)

        # Extract the input variants as internal delins
        sequence = sequence_from_description(d)
        input_variants = [
            Variant.from_model(delins, sequence=sequence)
            for delins in d.delins_model["variants"]
        ]

        results = list()

        # Store the wildtype
        wt = Therapy(
            name="Wildtype",
            hgvsc=hgvs.split("c.")[0] + "c.=",
            hgvsr=hgvs.split("c.")[0] + "r.=",
            hgvsp=protein_prediction(d, [])[0],
            description="These are the annotations as defined on the reference. They are always 100% by definition.",
            variants=list(),
        )
        wildtype = Result(wt, self.compare(self))
        results.append(wildtype)

        # Store the input variants as Therapy
        input = Therapy(
            name="Input",
            hgvsc=hgvs,
            hgvsr=hgvs.replace(":c.", ":r."),
            hgvsp=protein_prediction(d, input_variants)[0],
            description="The annotations based on the supplied input variants.",
            variants=input_variants,
            figure=draw(d),
        )
        patient = deepcopy(self)
        patient.mutate(d, input.variants)
        results.append(Result(input, patient.compare(self)))

        # Generate all possible therapies
        for therapy in generate_therapies(d):

            # Apply the combination to the wildtype transcript
            modified_transcript = deepcopy(self)
            modified_transcript.mutate(d, therapy.variants)
            results.append(Result(therapy, modified_transcript.compare(self)))

        # Order the results
        wt_patient = results[:2]
        rest = sorted(results[2:], reverse=True)
        return wt_patient + rest

    def __str__(self) -> str:
        return "\n".join(str(record) for record in self.records())

    @property
    def exons(self) -> Bed | None:
        """Helper function to get the exons for Transcript

        Note that this relies on the Bed record for the exons to be
        named 'Exons'
        """
        for record in self.rna_features:
            if record.name == "Exons":
                return record
        else:
            return None

    @property
    def coding_exons(self) -> Bed | None:
        """Helper function to get the coding exons for Transcript

        Note that this relies on the Bed record for the coding exons to be
        named 'Coding exons'
        """
        for record in self.protein_features:
            if record.name == "Coding exons":
                return record
        else:
            return None

    def lookup_protein_domains(self, d: Description) -> None:
        """Lookup supported protein domains from USCS"""
        for track in PROTEIN_TRACKS:
            self.protein_features += lookup_track(d, track)
