import dataclasses
import itertools
import logging
from copy import deepcopy
from typing import Any, Mapping, Self, Sequence

from algebra.extractor import extract_sequence, to_hgvs
from mutalyzer.description import Description
from mutalyzer_crossmapper import Coding

from gtgt.range import intersect

from .bed import Bed
from .exonviz import draw
from .mutalyzer import (
    genomic_crossmapper,
    genomic_to_transcript,
    get_chrom_name,
    get_strand,
    init_description,
    mutation_to_cds_effect,
    sequence_from_description,
    transcript_crossmapper,
    variants_to_hgvs,
    variants_to_protein_hgvs,
)
from .therapy import TherapyReport, generate_therapies
from .ucsc import PROTEIN_TRACKS, lookup_track
from .variant import Variant

logger = logging.getLogger(__name__)


mapping = str.maketrans(
    {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "a": "u",
        "c": "g",
        "g": "c",
        "u": "a",
    }
)


def reverse_complement(seq: str) -> str:
    return seq.translate(mapping)[::-1]


@dataclasses.dataclass
class Comparison:
    name: str
    percentage: float
    basepairs: str

    @classmethod
    def from_dict(cls, dict: Mapping[str, Any]) -> "Comparison":
        """Create a Variant object from a dict representation of a Variant"""
        return cls(**dict)

    def __hash__(self) -> int:
        return hash(f"{self.name} {self.basepairs}")


@dataclasses.dataclass
class Result:
    """To hold results for separate mutations of a transcript"""

    therapy: TherapyReport
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
            therapy=TherapyReport.from_dict(dict["therapy"]),
            comparison=[Comparison.from_dict(x) for x in dict["comparison"]],
        )


Range = tuple[int, int]


class MiniTranscript:
    """
    Minimalistic transcript

    If you make a copy of MiniTranscript objects, make sure to replace features
    and not modify them, since members are shared between copies by default.
    """

    def __init__(
        self,
        inverted: bool,
        sequence: str,
        exons: list[Range],
        cds: list[Range],
        variants: list[Variant] | None = None,
    ) -> None:
        self.inverted = inverted
        self.sequence = sequence
        self.exons = exons
        self.cds = cds
        self.variants = tuple(variants) if variants else tuple()

    @classmethod
    def from_hgvs(cls, hgvs: str) -> Self:
        d = init_description(hgvs)

        # Get the data we need
        inverted = d.is_inverted()
        sequence = sequence_from_description(d)
        exons = d.get_selector_model()["exon"]
        cds = d.get_selector_model()["cds"]

        variants = [
            Variant.from_model(delins, sequence=sequence)
            for delins in d.delins_model["variants"]
        ]
        return cls(inverted, sequence, exons, cds, variants)

    def get_sequence(self) -> str:
        """Return the unchanged sequence of the transcript"""
        return self.sequence

    def get_exon_sequence(self) -> str:
        """Return the mutated exon sequence, in forward orientation"""
        exons = self._mutate_features(self.exons, self.variants)
        sequence = self._mutate_sequence(self.sequence, self.variants)

        return self._get_features_sequence(sequence, exons)

    def get_coding_sequence(self) -> str:
        """Return the mutated coding sequence, in forward orientation"""

        sequence = self._mutate_sequence(self.sequence, self.variants)
        exons = self._mutate_features(self.exons, self.variants)
        cds = self._mutate_features(self.cds, self.variants)[0]

        coding = self._get_coding_regions(exons, cds)

        return self._get_features_sequence(sequence, coding)

    def get_protein_sequence(self) -> str:
        """Return the mutated protein sequence"""
        coding_sequence = self.get_coding_sequence()
        return self._get_protein_sequence(coding_sequence)

    def protein_hgvs(self) -> str:
        """Generate the HGVS description for the mutated protein"""
        # First, we get the wildtype protein
        wt_coding = self._get_coding_regions(self.exons, self.cds[0])
        wt_sequence = self._get_features_sequence(self.sequence, wt_coding)
        reference = self._get_protein_sequence(wt_sequence)

        # Next, we get the mutated protein sequence
        observed = self.get_protein_sequence()

        return self._protein_hgvs(reference, observed)

    def hgvs(self) -> str:
        """Generate the HGVS description for the variants"""
        # HGVS for no variants
        if not self.variants:
            return "="

        mapper = Coding(self.exons, self.cds, inverted=self.inverted)

        descriptions = [variant.to_hgvs(mapper) for variant in self.variants]

        if len(descriptions) == 1:
            return descriptions[0]
        else:
            return "[" + ";".join(descriptions) + "]"

    @staticmethod
    def _get_protein_sequence(coding_sequence: str) -> str:
        """Return the mutated protein sequence"""
        # Add padding if the coding sequence is not divisible by three
        if padding := 3 - len(coding_sequence) % 3:
            coding_sequence += "N" * padding

        predicted = str(Seq(coding_sequence).translate())  # type: ignore
        first_stop = predicted.find("*")
        return predicted[: first_stop + 1]

    @staticmethod
    def _protein_hgvs(reference: str, observed: str) -> str:
        """Generate HGVS description of protein changes (not normalized)"""
        canonical, _ = extract_sequence(reference, observed)
        return str(to_hgvs(canonical, reference))

    def _get_features_sequence(self, sequence: str, features: list[Range]) -> str:
        """Return the sequence of the features, in forward orientation"""
        segments = list()
        for start, end in features:
            segments.append(sequence[start:end])

        if self.inverted:
            return reverse_complement("".join(segments))
        else:
            return "".join(segments)

    @staticmethod
    def _get_coding_regions(exons: list[Range], cds: Range) -> list[Range]:
        """Return the ranges that make up the coding sequence"""

        coding: list[Range] = list()
        for exon in exons:
            coding += intersect(exon, cds)

        return coding

    @staticmethod
    def _mutate_features(
        features: Sequence[Range], variants: Sequence[Variant]
    ) -> list[Range]:
        """
        Mutate features based on the variants

        The variants do not have to be sorted
        """
        # Flatten the features, since we can act on them individually
        boundaries = list(itertools.chain.from_iterable(features))

        # Keep a running tally of the offsets for each boundary, because we cannot
        # update the boundaries until we have seen all the variants
        offsets = [0 for _ in boundaries]

        for variant in variants:
            for index in range(len(boundaries)):
                position = boundaries[index]

                # If the variant is before or at the position, we increase the offset
                # for this position by the size of the variant
                if position >= variant.end:
                    offsets[index] += variant.size

        # Apply the offsets
        boundaries = [
            position + offset for position, offset in zip(boundaries, offsets)
        ]

        # Group into Ranges
        new_features = list()
        for i in range(0, len(boundaries), 2):
            new_features.append((boundaries[i], boundaries[i + 1]))

        return new_features

    @staticmethod
    def _mutate_sequence(sequence: str, variants: Sequence[Variant]) -> str:
        """
        Mutate sequence based on the variants

        Note: the variants must be sorted
        """
        if not variants:
            return sequence

        # Check that the variants are sorted
        prev = variants[0].end
        for variant in variants[1:]:
            if variant.start < prev:
                raise ValueError("Variants must be sorted")
            prev = variant.end

        # The end of the last variant cannot be beyond the end of the sequence
        if variants[-1].end > len(sequence):
            raise ValueError(f"{variants[-1]} lies beyond the end of the sequence")

        # Store the segments of the mutated sequence
        segments = list()

        # Keep track of the last position in the sequence we have included
        position = 0

        for variant in variants:
            # Store the segment before the current variant
            segments.append(sequence[position : variant.start])
            segments.append(variant.inserted)
            position = variant.end
        else:
            segments.append(sequence[position:])

        return "".join(segments)


class Transcript:
    """
    Class to hold a Transcript


    A Transcript can be defined on DNA (c.) or RNA (r.)

    Note that the mutalyzer Description object is the wildtype, the variants are defined separately
    """

    def __init__(self, hgvs: str) -> None:
        # Make the wildtype version of the transcript
        if ":c." in hgvs:
            _id = hgvs.split(":c.")[0]
            wildtype = f"{_id}:c.="
        elif ":r." in hgvs:
            _id = hgvs.split(":r.")[0]
            wildtype = f"{_id}:r.="
        else:
            raise NotImplementedError(hgvs)
        self.d = init_description(wildtype)

        # Get the variants
        self.variants = self._variants(hgvs)

    def _variants(self, hgvs: str) -> list[Variant]:
        d = init_description(hgvs)
        sequence = sequence_from_description(d)

        return [
            Variant.from_model(delins, sequence=sequence)
            for delins in d.delins_model["variants"]
        ]

    def hgvs(self) -> str:
        """Generate HGVS for this Transcript"""
        return variants_to_hgvs(self.d, self.variants)

    def protein_hgvs(self) -> str:
        """Generate the protein HGVS for this transcript"""
        return variants_to_protein_hgvs(self.d, self.variants)


class Features:
    """
    Class to hold the features of an Allele.

    The individual features are split onto two categories, RNA and protein
    features. The reason for this is that variants have a different effect on
    features on the RNA (such as poly-adenylation sites) and the protein.

    For example, an out of frame insertion will not lead to the loss of any of
    the RNA features, but it will invalidate all protein features downstream of
    the variant.
    """

    def __init__(self, rna: list[Bed], protein: list[Bed]):
        self.rna = rna
        self.protein = protein

    def __repr__(self) -> str:
        return f"{type(self).__name__}(rna={self.rna}, protein={self.protein})"

    @classmethod
    def from_hgvs(cls, hgvs: str) -> "Features":
        """Create a Transcript object from a mutalyzer Description"""
        d = init_description(hgvs)

        # Get exons and add the offset
        selector_model = d.get_selector_model()
        exons = selector_model["exon"]

        # Get CDS and add the offset
        cds = selector_model["cds"][0]

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

        return cls(rna=[exon_bed], protein=[coding_exons])

    def records(self) -> Sequence[Bed]:
        """Return all Bed records"""
        return self.rna + self.protein

    def rna_records(self) -> list[Bed]:
        """Return the Bed records that contain RNA features"""
        return self.rna

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

    @property
    def exons(self) -> Bed | None:
        """Helper function to get the exons

        Note that this relies on the Bed record for the exons to be
        named 'Exons'
        """
        for record in self.rna:
            if record.name == "Exons":
                return record
        else:
            return None

    @property
    def coding_exons(self) -> Bed | None:
        """Helper function to get the coding exons

        Note that this relies on the Bed record for the coding exons to be
        named 'Coding exons'
        """
        for record in self.protein:
            if record.name == "Coding exons":
                return record
        else:
            return None

    def compare(self, other: object) -> Sequence[Comparison]:
        """Compare the size of each record in the transcripts"""
        if not isinstance(other, Features):
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
        if not isinstance(other, Allele):
            raise NotImplementedError
        cmp = self.compare(other)

        values = [x.percentage for x in cmp]
        return sum(values) / len(cmp)

    def mutate(self, d: Description, variants: Sequence[Variant]) -> None:
        """
        Mutate the transcript based on the specified variants

        Please make sure that that Variants use the same internal coordinate
        system as the features
        """
        # Determine the chromosome the transcript is on
        if self.records():
            chrom = self.records()[0].chrom
        else:
            chrom = ""

        # Update protein features
        protein_changes = Bed.from_blocks(chrom, mutation_to_cds_effect(d, variants))
        for record in self.protein:
            record.subtract(protein_changes)

        # Update RNA features
        rna_changes = Bed.from_blocks(chrom, [(v.start, v.end) for v in variants])
        for record in self.rna_records():
            self.subtract(rna_changes)


class Allele:
    """
    The Allele class is the representation of one of the alleles of a putative
    patient.

    An Allele consists of:
        - Features, see the Features class
        - The RNA Transcript
        - The DNA Transcript

    For some variants, the RNA and DNA versions of the transcript will be
    (practically) identical. However, in the case of splice altering variants
    there can be big differences between the RNA and the DNA version of the
    Transcript.

    The Features will use the internal coordinate system from the RNA
    Transcript.
    """

    def __init__(
        self,
        dna: str,
        rna: str,
        features: Features | None = None,
        _ignore_coordinate: bool = False,
    ) -> None:

        # Check if the coordinate systems are correct, unless _ignore_coordinate
        # has been set (this is used when testing)
        if ":c." not in dna and not _ignore_coordinate:
            msg = f"Only variants on the 'c' coordinate system are supported ({dna})"
            raise ValueError(msg)

        if ":r." not in rna and not _ignore_coordinate:
            msg = f"Only variants on the 'r' coordinate system are supported ({dna})"
            raise ValueError(msg)

        self.dna = Transcript(dna)

        # If there are no features specified, we use the features from the RNA
        if features is None:
            features = Features.from_hgvs(rna)
        self.features = features

        self.rna = Transcript(rna)

        # Finally, mutate the features based on the RNA to ensure internal
        # consistency
        self.mutate()

    def mutate(self) -> None:
        """Mutate the features based on the variants of the RNA Transcript"""
        self.features.mutate(self.rna.d, self.rna.variants)

    def __str__(self) -> str:
        return "\n".join(str(record) for record in self.features.records())

    def lookup_protein_domains(self, d: Description) -> None:
        """Lookup supported protein domains from USCS"""
        g_crossmap = genomic_crossmapper(d.input_description)
        t_crossmap = transcript_crossmapper(d)
        for track in PROTEIN_TRACKS:
            try:
                features = lookup_track(d, track)
            except RuntimeError as e:
                logger.error(e)
                continue

            # Convert each feature to the internal transcript coordinate system
            for record in features:
                blocks = record.blocks()
                t_blocks = list()

                for start, end in blocks:
                    try:
                        start, end = genomic_to_transcript(
                            start, end, g_crossmap, t_crossmap
                        )
                    except ValueError as e:
                        logger.warning(str(e))
                        continue

                    t_blocks.append((start, end))

                record.update(t_blocks)
                self.features.protein.append(record)


def analyze(allele: Allele, extended: bool = False) -> Sequence[Result]:
    """Analyze the transcript based on the specified HGVS description"""

    # Create the wildtype allele as a starting point
    hgvs_id = allele.rna.d.input_description.split(":r.")[0]
    wildtype = Allele(rna=f"{hgvs_id}:r.=", dna=f"{hgvs_id}:c.=")

    # If there are protein domains on the allele, we also want them for the wildtype
    if len(allele.features.protein) > 1:
        wildtype.lookup_protein_domains(wildtype.rna.d)

    # TODO restore protein features

    results = list()

    # Store the wildtype
    wt = TherapyReport(
        name="Wildtype",
        hgvsc=wildtype.dna.hgvs(),
        hgvsr=wildtype.rna.hgvs(),
        hgvsp=wildtype.rna.protein_hgvs(),
        description="These are the annotations as defined on the reference. They are always 100% by definition.",
        dna_variants=wildtype.dna.variants,
        rna_variants=wildtype.rna.variants,
    )
    results.append(Result(wt, wildtype.features.compare(wildtype.features)))

    # Store the input variants as Therapy
    input = TherapyReport(
        name="Input",
        hgvsc=allele.dna.hgvs(),
        hgvsr=allele.rna.hgvs(),
        hgvsp=allele.rna.protein_hgvs(),
        description="The annotations based on the supplied input variants.",
        dna_variants=allele.dna.variants,
        rna_variants=allele.rna.variants,
        figure=draw(allele.rna),
    )
    # The input variant is actually the same as allele
    results.append(Result(input, allele.features.compare(wildtype.features)))

    # Generate all possible therapies
    for therapy in generate_therapies(allele):

        # Make the allele modified by the therapies
        modified_transcript = deepcopy(wildtype)

        modified_transcript.rna.variants = therapy.rna_variants.copy()
        modified_transcript.dna.variants = therapy.dna_variants.copy()
        # Update the features of the modified transcript based on the therapy
        modified_transcript.mutate()

        # Filter variants which are not of interest
        if not extended:
            if not is_of_interest(
                allele.features,
                allele.rna.variants,
                modified_transcript.features,
                modified_transcript.rna.variants,
            ):
                continue

        # Make the therapy report
        therapy_report = TherapyReport.from_therapy(therapy, modified_transcript)
        results.append(
            Result(
                therapy_report,
                modified_transcript.features.compare(wildtype.features),
            )
        )

    # Order the results
    wt_patient = results[:2]
    rest = sorted(results[2:], reverse=True)
    return wt_patient + rest


def is_of_interest(
    patient: Features,
    patient_vars: Sequence[Variant],
    therapy: Features,
    therapy_vars: Sequence[Variant],
) -> bool:
    """Determine if therapy is of interest for patient

    There are two criteria that determine of a therapy is of interest
    1: If any patient variants are modified by the therapy
    2: If any of the annotations have **improved** for the therapy
    """
    # If any patient variants are modified
    missing_variants = set(patient_vars).difference(therapy_vars)
    if missing_variants:
        return True

    # Determine if there are any regions that were restored in the therapy,
    # which are lacking in the patient
    for t, p in zip(therapy.records(), patient.records()):
        restored = deepcopy(t)
        restored.subtract(p)
        if restored.size > 0:
            return True
    else:
        return False
