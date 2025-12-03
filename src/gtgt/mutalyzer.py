import dataclasses
import logging
from copy import deepcopy
from typing import Any, Mapping, Sequence, TypeVar

import Levenshtein
import mutalyzer_hgvs_parser
from mutalyzer.converter.to_hgvs_coordinates import to_hgvs_locations
from mutalyzer.converter.variants_de_to_hgvs import (
    delins_to_del,
    delins_to_delins,
    delins_to_duplication,
    delins_to_insertion,
    delins_to_repeat,
    delins_to_substitution,
    get_end,
    get_start,
    is_duplication,
    is_repeat,
)
from mutalyzer.description import Description
from mutalyzer.description_model import (
    get_reference_id,
    get_selector_id,
    variants_to_description,
)
from mutalyzer.protein import get_protein_description, in_frame_description
from mutalyzer.reference import get_protein_selector_model
from mutalyzer.util import get_inserted_sequence, get_location_length
from mutalyzer_crossmapper import Coding
from pydantic import BaseModel, model_validator
from schema import And, Optional, Or, Schema
from typing_extensions import NewType

logger = logging.getLogger(__name__)

# Mutalyzer Variant dictionary
Variant_Dict = NewType("Variant_Dict", Mapping[str, Any])


def sequence_from_description(d: Description) -> str:
    """Return the sequence form a description"""
    _id = d.input_model["reference"]["id"]
    sequence: str = d.references[_id]["sequence"]["seq"]
    return sequence


class Variant:
    """Class to store delins variants"""

    # fmt: off
    # Schema for the location specification of the indel model
    location_schema = Schema(
        {
            "type": "range",
            "start": {
                "type": "point",
                "position": int,
            },
            "end": {
                "type": "point",
                "position": int,
            },
        }
    )

    # Schema for the inserted/deleted entries of the indel model
    inserted_deleted_schema = Schema(
        And( # Inserted must be 0 or 1 items
            lambda n: len(n) <= 1,
            [
                {
                    "sequence": Or(str, []),
                    "source": "description",
                    Optional("inverted") : True
                },
            ],
        ),
    )

    # Full schema for the indel model
    schema = Schema(
        {
            "type": "deletion_insertion",
            "source": "reference",
            "location": location_schema,
            Optional("inserted"): inserted_deleted_schema,
            Optional("deleted"): inserted_deleted_schema,
        }
    )
    # fmt: on

    def __init__(
        self,
        start: int,
        end: int,
        inserted: str = "",
        deleted: str = "",
    ):
        if start > end:
            raise ValueError(f"End ({end}) must be after start ({start})")
        self.start = start  # zero based
        self.end = end  # exclusive
        self.inserted = inserted

        if len(deleted) > 1:
            raise ValueError("deleted sequence is only defined for SNPS, not indels")
        self.deleted = deleted

    def __str__(self) -> str:
        return self.__repr__()

    def __repr__(self) -> str:
        start = self.start
        end = self.end
        inserted = self.inserted
        deleted = self.deleted
        return f"Variant({start=}, {end=}, inserted={inserted}, deleted={deleted})"

    def before(self, other: "Variant") -> bool:
        return self.end <= other.start

    def after(self, other: "Variant") -> bool:
        return self.start >= other.end

    def inside(self, other: "Variant") -> bool:
        return self.start >= other.start and self.end <= other.end

    def overlap(self, other: "Variant") -> bool:
        self_ends_in_other = self.end > other.start and self.end <= other.end
        self_starts_in_other = self.start >= other.start and self.start < other.end

        return any(
            [
                self_starts_in_other,
                self_ends_in_other,
                self.inside(other),
                other.inside(self),
            ]
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Variant):
            raise NotImplementedError

        return (
            self.start == other.start
            and self.end == other.end
            and self.inserted == other.inserted
            and self.deleted == other.deleted
        )

    def __lt__(self, other: "Variant") -> bool:
        if not isinstance(other, Variant):
            raise NotImplementedError
        if self.overlap(other):
            msg = f"Overlapping variants '{self}' and '{other}' cannot be sorted"
            raise ValueError(msg)

        return self.start < other.start

    @staticmethod
    def _validate_schema(model: Mapping[str, Any]) -> None:
        """Validate the structure of the mutalyzer delins model

        This can be very complex, and we only support the most common cases.
        """

        Variant.schema.validate(model)

    @staticmethod
    def _model_is_repeat(model: Mapping[str, Any]) -> bool:
        """Determine if model is a repeat"""
        # Determine if the model is a repeat
        repeat_schema = Schema(
            And(  # An empty list would be a deletion, not a repeat
                lambda n: len(n) == 1,
                [
                    {
                        "source": "description",
                        "sequence": str,
                        "repeat_number": {"type": "point", "value": int},
                        Optional("inverted"): True,
                    }
                ],
            )
        )

        inserted = model.get("inserted")
        is_repeat: bool = repeat_schema.is_valid(inserted)

        return is_repeat

    @staticmethod
    def _model_repeat_to_delins(model: Mapping[str, Any]) -> dict[str, Any]:
        """Convert a repeat model to a delins model"""
        new_model = {k: v for k, v in model.items()}

        # Determine the sequence and repeats
        old_sequence = model["inserted"][0]["sequence"]
        repeats = model["inserted"][0]["repeat_number"].get("value", 1)

        # Expand the new sequence
        new_sequence = old_sequence * repeats

        inserted = {"sequence": new_sequence, "source": "description"}

        # If the model is defined on the reverse strand
        inverted = model["inserted"][0].get("inverted", False)
        if inverted:
            inserted["inverted"] = True

        new_model["inserted"] = [inserted]
        return new_model

    @staticmethod
    def _model_is_duplication(model: Mapping[str, Any]) -> bool:
        """Determine if model is a duplication"""
        # Determine if the model is a repeat
        location_schema = Schema(
            {
                "type": "range",
                "start": {
                    "type": "point",
                    "position": int,
                },
                "end": {
                    "type": "point",
                    "position": int,
                },
            }
        )
        duplication_schema = Schema(
            And(  # An empty list would be a deletion, not a repeat
                lambda n: len(n) == 1,
                [
                    {
                        "source": "reference",
                        "location": location_schema,
                        "repeat_number": {"type": "point", "value": int},
                    }
                ],
            )
        )

        inserted = model.get("inserted")
        is_duplication: bool = duplication_schema.is_valid(inserted)

        return is_duplication

    @staticmethod
    def _model_duplication_to_delins(
        model: Mapping[str, Any], sequence: str
    ) -> dict[str, Any]:
        """Convert a duplication model to a delins model"""
        if not sequence:
            raise ValueError("Variant: specify sequence to handle duplications")
        new_model = {k: v for k, v in model.items()}

        # Determine the start and end on the sequence
        inserted = model["inserted"][0]
        start = inserted["location"]["start"]["position"]
        end = inserted["location"]["end"]["position"]
        repeats = inserted["repeat_number"].get("value", 1)

        # Expand the new sequence
        new_sequence = sequence[start:end]

        new_model["inserted"] = [
            {
                "sequence": new_sequence,
                "repeat_number": {"type": "point", "value": repeats},
                "source": "description",
            }
        ]
        return new_model

    @staticmethod
    def _model_is_inversion(model: Mapping[str, Any]) -> bool:
        """Determine if model is an inversion"""
        # Determine if the model is a repeat
        location_schema = Schema(
            {
                "type": "range",
                "start": {
                    "type": "point",
                    "position": int,
                },
                "end": {
                    "type": "point",
                    "position": int,
                },
            }
        )
        duplication_schema = Schema(
            And(  # An empty list would be a deletion, not a repeat
                lambda n: len(n) == 1,
                [
                    {
                        "source": "reference",
                        "location": location_schema,
                        "inverted": True,
                    }
                ],
            )
        )

        inserted = model.get("inserted")
        is_duplication: bool = duplication_schema.is_valid(inserted)

        return is_duplication

    @staticmethod
    def _reverse_complement(seq: str) -> str:
        complement = {"A": "T", "T": "A", "G": "C", "C": "G"}
        return "".join((complement[nt] for nt in reversed(seq)))

    @staticmethod
    def _model_inversion_to_delins(
        model: Mapping[str, Any], sequence: str
    ) -> dict[str, Any]:
        """Convert an inversion model to a delins model"""
        if not sequence:
            raise ValueError("Variant: specify sequence to handle inversions")
        new_model = {k: v for k, v in model.items()}

        # Determine the start and end on the sequence
        inserted = model["inserted"][0]
        start = inserted["location"]["start"]["position"]
        end = inserted["location"]["end"]["position"]

        # Expand the new sequence
        new_sequence = Variant._reverse_complement(sequence[start:end])

        new_model["inserted"] = [
            {
                "sequence": new_sequence,
                "source": "description",
            }
        ]
        return new_model

    @classmethod
    def from_model(cls, model: Mapping[str, Any], sequence: str = "") -> "Variant":
        if Variant._model_is_inversion(model):
            model = Variant._model_inversion_to_delins(model, sequence)
        if Variant._model_is_duplication(model):
            model = Variant._model_duplication_to_delins(model, sequence)
        # Determine if the model is a repeat
        if Variant._model_is_repeat(model):
            model = Variant._model_repeat_to_delins(model)

        # Validate the delins model
        cls._validate_schema(model)

        start = model["location"]["start"]["position"]
        end = model["location"]["end"]["position"]

        # Store if the inserted and deleted sequences were inverted
        ins_inverted: bool | None = None
        del_inverted: bool | None = None

        inserted = model.get("inserted", [])

        # Sanity check, complex variants should have been normalized into a
        # delins representation
        if len(inserted) == 0:
            inserted = ""
        elif len(inserted) > 1:
            raise NotImplementedError(
                "Multiple records in 'inserted' are not supported"
            )
        else:  # inserted contains 1 item, as is expected
            if "sequence" not in inserted[0]:
                raise NotImplementedError(
                    "Missing sequence in 'inserted' is not supported"
                )
            if "repeat_number" in inserted[0]:
                raise NotImplementedError("Repeats in 'inserted' are not supported")

            ins_inverted = inserted[0].get("inverted", False)
            inserted = inserted[0]["sequence"]

        deleted = model.get("deleted")
        if deleted is not None:
            try:
                del_inverted = deleted[0].get("inverted", False)
                deleted = deleted[0]["sequence"]
            except KeyError:
                raise NotImplementedError("Complex Variant not supported")

        # Make sure inserted and deleted are both or neither inversed
        if ins_inverted is not None and del_inverted is not None:
            if ins_inverted != del_inverted:
                msg = "strand difference between inserted and deleted sequences are not supported"
                raise NotImplementedError(msg)

        # If inserted or deleted are inverted, we need to reverse complement
        # the inserted/deleted nucleotides
        if inserted and ins_inverted:
            inserted = cls._reverse_complement(inserted)
        if deleted and del_inverted:
            deleted = cls._reverse_complement(deleted)

        return Variant(
            start=start,
            end=end,
            inserted=inserted if inserted else "",
            deleted=deleted if deleted else "",
        )

    @classmethod
    def from_dict(cls, dict: Mapping[str, Any]) -> "Variant":
        """Create a Variant object from a dict representation of a Variant"""
        return cls(**dict)

    def to_model(self) -> Mapping[str, Any]:
        """Convert Variant to mutalyzer delins model"""

        # Specification of the location
        # fmt: off
        location = {
            "type": "range",
            "start": {
                "type": "point",
                "position": self.start
            },
            "end": {
                "type": "point",
                "position": self.end
            }
        }

        # Specification of the inserted sequence
        inserted_obj: dict[str, Any] = {
            "sequence": self.inserted,
            "source": "description"
        }

        if self.inserted:
            inserted = [inserted_obj]
        else:
            inserted = []
        # fmt: on

        model = {
            "location": location,
            "type": "deletion_insertion",
            "source": "reference",
            "inserted": inserted,
        }

        deletion_object: dict[str, Any] = {
            "sequence": self.deleted,
            "source": "description",
        }
        if self.deleted:
            model["deleted"] = [deletion_object]

        return model

    def genomic_coordinates(self, d: Description) -> tuple[int, int]:
        """Return genomic coordinates for Variant"""
        offset = get_offset(d)

        if offset is None:
            raise RuntimeError("Missing ensembl offset")

        return self.start + offset, self.end + offset


@dataclasses.dataclass
class Therapy:
    """Class to store genetic therapies"""

    name: str
    hgvsc: str
    description: str
    variants: Sequence[Variant]
    figure: str | None = None
    hgvsr: str | None = None
    hgvsp: str | None = None

    @classmethod
    def from_dict(cls, dict: Mapping[str, Any]) -> "Therapy":
        """Create a Therapy object from a dict representation of a Therapy"""
        v = [Variant.from_dict(x) for x in dict["variants"]]
        return cls(
            name=dict["name"],
            hgvsc=dict["hgvsc"],
            description=dict["description"],
            variants=v,
            figure=dict.get("figure"),
            hgvsr=dict.get("hgvsr"),
            hgvsp=dict.get("hgvsp"),
        )


class HGVS(BaseModel):
    description: str

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "description": "ENST00000357033.9:c.6439_6614del",
                }
            ]
        }
    }

    @model_validator(mode="after")
    def hgvs_parser(self) -> "HGVS":
        """Parse the HGVS description with mutalyzer-hgvs-parser"""
        hgvs_error = (
            mutalyzer_hgvs_parser.exceptions.UnexpectedCharacter,
            mutalyzer_hgvs_parser.exceptions.UnexpectedEnd,
        )
        try:
            mutalyzer_hgvs_parser.to_model(self.description)
        except hgvs_error as e:
            raise ValueError(e)
        return self


def combine_variants_deletion(
    variants: Sequence[Variant], deletion: Variant
) -> Sequence[Variant]:
    """Combine variants and a deletion, any variants that are contained in
    the deletion are discarded

    The resulting list of variants is sorted
    """
    if deletion.inserted:
        raise ValueError(f"{Variant} is not a pure deletion")

    # Ensure the variants are sorted, and do not overlap
    sorted_variants = sorted(variants)

    combined = list()
    for i in range(len(sorted_variants)):
        variant = sorted_variants[i]
        if variant.before(deletion):
            combined.append(variant)
        elif variant.inside(deletion):
            # Discard the current variant
            continue
        elif variant.after(deletion):
            combined.append(deletion)
            combined += sorted_variants[i:]
            break
        else:
            msg = f"Deletion '{deletion}' partially overlaps '{variant}"
            raise ValueError(msg)
    else:
        combined.append(deletion)

    return combined


def de_to_hgvs(variants: Any, sequences: Any) -> Sequence[Variant_Dict]:
    """
    Convert the description extractor variants to an HGVS format (e.g., a
    deletion insertion of one nucleotide is converted to a substitution).

    MODIFIED from mutalyzer to not perform the 3' shift
    """
    if len(variants) == 1 and variants[0].get("type") == "equal":
        new_variant = deepcopy(variants[0])
        new_variant.pop("location")
        return [new_variant]

    new_variants = []
    for variant in variants:
        if variant.get("type") == "inversion":
            new_variants.append(deepcopy(variant))
        elif variant.get("type") == "deletion_insertion":
            inserted_sequence = get_inserted_sequence(variant, sequences)
            if len(inserted_sequence) == 0:
                new_variants.append(delins_to_del(variant))
            elif (
                get_location_length(variant["location"]) == len(inserted_sequence) == 1
            ):
                new_variants.append(delins_to_substitution(variant, sequences))
            elif is_repeat(variant, sequences):
                new_variants.append(delins_to_repeat(variant, sequences))
            elif is_duplication(variant, sequences):
                new_variants.append(delins_to_duplication(variant, sequences))
            elif get_start(variant["location"]) == get_end(variant["location"]):
                new_variants.append(delins_to_insertion(variant))
            else:
                new_variants.append(delins_to_delins(variant))

    return new_variants


def to_cdot_hgvs(d: Description, variants: Sequence[Variant]) -> str:
    """Convert a list of _Variants to hgvs representation"""
    delins_model = [v.to_model() for v in variants]

    # Invert the deleted sequence if the transcript is on the reverse strand
    if d.is_inverted():
        for delins in delins_model:
            if "deleted" in delins:
                _del = delins["deleted"][0]
                _del["sequence"] = Variant._reverse_complement(_del["sequence"])
                _del["inverted"] = True

    variant_models = de_to_hgvs(delins_model, d.get_sequences())

    ref_id = get_reference_id(d.corrected_model)
    selector_id = get_selector_id(d.corrected_model)

    selector_model = get_protein_selector_model(
        reference=d.references[ref_id]["annotations"], selector_id=selector_id
    )

    description_model = {
        "type": "description_dna",
        "reference": {"id": ref_id, "selector": {"id": ref_id}},
        "coordinate_system": "c",
        "variants": variant_models,
    }

    cdot_locations = to_hgvs_locations(
        description_model,
        d.references,
        selector_model=selector_model,
    )["variants"]

    hgvs: str = variants_to_description(cdot_locations)
    return hgvs


T = TypeVar("T")


def sliding_window(items: Sequence[T], size: int = 1) -> Sequence[Sequence[T]]:
    adj: list[Sequence[T]] = list()
    for i in range(len(items) - size + 1):
        adj.append([x for x in items[i : i + size]])
    return adj


def _exon_string(exon_numbers: Sequence[int]) -> str:
    """Format the exon names for a variable number of exons"""
    if len(exon_numbers) == 1:
        return f"exon {exon_numbers[0]}"
    elif len(exon_numbers) == 2:
        return f"exons {exon_numbers[0]} and {exon_numbers[1]}"
    else:
        t = ", ".join((str(x) for x in exon_numbers[:-1]))
        return f"exons {t} and {exon_numbers[-1]}"


def skip_adjacent_exons(d: Description, number_to_skip: int = 1) -> Sequence[Therapy]:
    """Skipp all possible adjacent exons the specified Description"""
    exon_skips: list[Therapy] = list()

    skippable_exons = get_exons(d, in_transcript_order=True)[1:-1]
    sequence = sequence_from_description(d)
    variants = [
        Variant.from_model(v, sequence=sequence) for v in d.delins_model["variants"]
    ]
    logger.debug(f"{variants=}")

    logger.debug(f"{skippable_exons=}")
    for i, exons in enumerate(sliding_window(skippable_exons, size=number_to_skip), 2):
        # Generate the string of exon numbers
        exons_description = _exon_string(range(i, i + number_to_skip))
        logger.debug(f"{exons_description=}")
        logger.debug(f"{exons=}")

        if d.is_inverted():
            # Start of the first exon to skip
            start = exons[-1][0]
            # End of the last exon to skip
            end = exons[0][-1]

        else:
            start = exons[0][0]
            end = exons[-1][-1]

        exon_skip = Variant(start, end)
        logger.debug(f"{exon_skip=}")

        # Combine the existing variants with the exon skip
        try:
            combined = combine_variants_deletion(variants, exon_skip)
        except ValueError as e:
            if number_to_skip == 1:
                msg = f"Cannot skip exon {exons_description}: {e}"
            else:
                msg = f"Cannot skip exons {exons_description}: {e}"
            logger.warn(msg)
            continue
        logger.debug(f"{combined=}")

        description = f"The annotations based on the supplied variants, in combination with skipping {exons_description}."
        # Convert to c. notation (user facing)
        name = f"Skip {exons_description}"
        selector = d.get_selector_id()
        cdot_variants = to_cdot_hgvs(d, combined)
        hgvsc = f"{selector}:c.{cdot_variants}"
        hgvsr = f"{selector}:r.{cdot_variants}"
        description = description
        t = Therapy(
            name=name,
            hgvsc=hgvsc,
            hgvsr=hgvsr,
            hgvsp=protein_prediction(d, combined)[0],
            description=description,
            variants=combined,
        )
        exon_skips.append(t)

    return exon_skips


def generate_therapies(d: Description) -> Sequence[Therapy]:
    """Wrapper around the different therapies"""
    therapies: list[Therapy] = list()
    # Skip a single exon
    therapies += skip_adjacent_exons(d, number_to_skip=1)
    # Skip two adjacent exons
    therapies += skip_adjacent_exons(d, number_to_skip=2)
    return therapies


def init_description(hgvs: str) -> Description:
    """
    Generate and initialize a Description for the specified HGVS

    Doesn't normalize the positions
    """
    d = Description(hgvs, stop_on_error=True)

    d.to_delins()
    d.de_hgvs_internal_indexing_model = d.delins_model
    d.construct_de_hgvs_internal_indexing_model()
    d.construct_de_hgvs_coordinates_model()
    d.construct_normalized_description()
    d.construct_protein_description()

    return d


def get_offset(d: Description) -> int:
    ref_id = get_reference_id(d.corrected_model)
    offset: int = (
        d.references.get(ref_id, {})
        .get("annotations", {})
        .get("qualifiers", {})
        .get("location_offset", 0)
    )
    return offset


def get_chrom_name(d: Description) -> str:
    # Try ensembl chromosome
    chrom = _get_ensembl_chrom_name(d)
    if chrom is not None:
        return str(chrom)

    # Try NCBI chromosome
    chrom = _get_ncbi_chrom_name(d)
    if chrom is not None:
        return str(chrom)

    raise RuntimeError(f"Unable to determine the chromosome for {d}")


def _get_ensembl_chrom_name(d: Description) -> str | None:
    ref_id = get_reference_id(d.corrected_model)
    chrom_name: str | None = (
        d.references.get(ref_id, {})
        .get("annotations", {})
        .get("qualifiers", {})
        .get("chromosome_number")
    )
    return chrom_name


def _get_ncbi_chrom_name(d: Description) -> str | None:
    chrom_name: str | None = d.input_model["reference"].get("id")
    return chrom_name


def get_strand(d: Description) -> str:
    # NM transcripts are on the forward strand by definition, and their strand
    # is not defined in the selector model from mutalyzer. Therefore, we set 1
    # as the default value.
    strand = d.get_selector_model()["location"].get("strand", 1)
    if strand == 1:
        return "+"
    elif strand == -1:
        return "-"
    else:
        raise ValueError(f"Unknown strand for Description object {d}")


def changed_protein_positions(
    reference: str, observed: str
) -> Sequence[tuple[int, int]]:
    """
    Extract the change protein positions (0 based)
    """
    deleted = list()
    for op in Levenshtein.opcodes(reference, observed):
        operation = op[0]
        ref_start = op[1]
        ref_end = op[2]

        if operation == "equal":
            continue
        elif operation == "insert":
            continue
        elif operation == "replace":
            deleted.append((ref_start, ref_end))
        elif operation == "delete":
            deleted.append((ref_start, ref_end))

    return deleted


def protein_prediction(
    d: Description, variants: Sequence[Variant]
) -> tuple[str, str, str]:
    """Call mutalyzer get_protein_description on a Description and list of Variants"""
    # Get required data structures from the Description
    ref_id = get_reference_id(d.corrected_model)
    selector_id = get_selector_id(d.corrected_model)
    selector_model = get_protein_selector_model(
        d.references[ref_id]["annotations"], selector_id=selector_id
    )

    # Convert the Variants to their delins model representation
    delins = [v.to_model() for v in variants]

    description, reference, predicted, *rest = get_protein_description(
        delins, d.references, selector_model
    )

    # If mutalyzer puts in an unknown protein prediction, we overwrite it
    if description.endswith(":p.?"):
        id = description.split(":p.?")[0]
        desc = in_frame_description(reference, predicted)[0]
        description = f"{id}:p.{desc}"

    return description, reference, predicted


def mutation_to_cds_effect(
    d: Description, variants: Sequence[Variant]
) -> list[tuple[int, int]]:
    """
    Determine the effect of the specified HGVS description on the CDS, on the genome

    Steps:
    - Use the protein prediction of mutalyzer to determine which protein
      residues are changed
    - Map this back to a deletion in c. positions to determine which protein
      annotations are no longer valid
    - Convert the c. positions to genome coordiinates as used by the UCSC
    NOTE that the genome range is similar to the UCSC annotations on the genome,
    i.e. 0 based, half open. Not to be confused with hgvs g. positions
    """
    # Determine the protein positions that were changed
    protein = protein_prediction(d, variants)
    reference, observed = protein[1], protein[2]

    # Keep track of changed positions on the genome
    changed_genomic = list()

    # Create crossmapper
    exons = d.get_selector_model()["exon"]
    cds = d.get_selector_model()["cds"]
    assert len(cds) == 1
    crossmap = Coding(exons, cds[0], inverted=d.is_inverted())

    for start, end in changed_protein_positions(reference, observed):
        # Calculate the nucleotide changed amino acids into a deletion in HGVS c. format

        logger.debug(f"p.{start+1}_{end}")

        # Internal coordinate positions
        start = crossmap.protein_to_coordinate(
            # nth amino acid, first nt of codon
            (start + 1, 1, 0, 0, 0)
        )
        end = (
            # nth amino acid, last nt of codon
            crossmap.protein_to_coordinate((end, 3, 0, 0, 0))
            + 1
        )

        if end < start:
            start, end = end - 1, start + 1

        v = Variant(start, end)

        changed_genomic.append(v.genomic_coordinates(d))

    return changed_genomic


def get_exons(
    description: Description, in_transcript_order: bool
) -> Sequence[tuple[int, int]]:
    """Get exons from a Mutalyzer Description object

    Positions are in the internal coordinate system
    """
    exons: Sequence[tuple[int, int]] = description.get_selector_model()["exon"]
    if in_transcript_order and description.is_inverted():
        return exons[::-1]

    return exons
