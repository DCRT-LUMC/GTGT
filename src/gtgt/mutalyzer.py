import dataclasses
import logging
from copy import deepcopy
from typing import Any, Dict, List, Sequence, Tuple, TypeVar, Union

import Levenshtein
import mutalyzer_hgvs_parser
from mutalyzer.converter.to_delins import variants_to_delins
from mutalyzer.converter.to_hgvs_coordinates import to_hgvs_locations
from mutalyzer.converter.to_internal_coordinates import to_internal_coordinates
from mutalyzer.converter.to_internal_indexing import to_internal_indexing
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
from mutalyzer.description_model import get_reference_id, variants_to_description
from mutalyzer.protein import get_protein_description
from mutalyzer.reference import get_protein_selector_model
from mutalyzer.util import get_inserted_sequence, get_location_length
from pydantic import BaseModel, model_validator
from typing_extensions import NewType

logger = logging.getLogger(__name__)

# Mutalyzer variant object, using the 'internal' coordinate system (0 based, half open)
# Variant string in HGVS c. format
CdotVariant = NewType("CdotVariant", str)
# Mutalyzer Variant dictionary
Variant_Dict = NewType("Variant_Dict", Dict[str, Any])
InternalVariant = NewType("InternalVariant", dict[str, Any])


class Variant:
    """Class to store delins variants"""

    def __init__(self, start: int, end: int, inserted: str = "", deleted: str = ""):
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
        )

    def __lt__(self, other: "Variant") -> bool:
        if not isinstance(other, Variant):
            raise NotImplementedError
        if self.overlap(other):
            msg = f"Overlapping variants '{self}' and '{other}' cannot be sorted"
            raise ValueError(msg)

        return self.start < other.start

    @classmethod
    def from_model(cls, model: Dict[str, Any]) -> "Variant":
        start = model["location"]["start"]["position"]
        end = model["location"]["end"]["position"]

        inserted = model.get("inserted", [])
        if inserted:
            inserted = inserted[0]["sequence"]
        else:
            inserted = ""

        deleted = model.get("deleted")
        if deleted is not None:
            deleted = deleted[0]["sequence"]

        if deleted:
            return Variant(start, end, inserted if inserted else "", deleted)
        else:
            return Variant(start, end, inserted if inserted else "")

    def to_model(self) -> Dict[str, Any]:
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
        if self.inserted:
            inserted = [
                {
                    "sequence": self.inserted,
                    "source": "description"
                }
            ]
        else:
            inserted = []
        # fmt: on

        model = {
            "location": location,
            "type": "deletion_insertion",
            "source": "reference",
            "inserted": inserted,
        }

        if self.deleted:
            deleted = [{"sequence": self.deleted, "source": "description"}]
            model["deleted"] = deleted

        return model


@dataclasses.dataclass
class Therapy:
    """Class to store genetic therapies"""

    name: str
    hgvs: str
    description: str
    variants: Sequence[Variant]


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
) -> List[Variant]:
    """Combine variants and a deletion, any variants that are contained in
    the deletion are discarded

    The resulting list of variants is sorted
    """
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


def de_to_hgvs(variants: Any, sequences: Any) -> List[Variant_Dict]:
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
    variant_models = de_to_hgvs(delins_model, d.get_sequences())

    ref_id = get_reference_id(d.corrected_model)

    selector_model = get_protein_selector_model(
        d.references[ref_id]["annotations"], ref_id
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


def sliding_window(items: Sequence[T], size: int = 1) -> List[List[T]]:
    adj: List[List[T]] = list()
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


def skip_adjacent_exons(d: Description, number_to_skip: int = 1) -> List[Therapy]:
    """Skipp all possible adjacent exons the specified Description"""
    exon_skips: List[Therapy] = list()

    skippable_exons = get_exons(d, in_transcript_order=True)[1:-1]
    variants = [Variant.from_model(v) for v in d.delins_model["variants"]]
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
        hgvs = f"{selector}:c.{cdot_variants}"
        description = description
        t = Therapy(name, hgvs, description, combined)
        exon_skips.append(t)

    return exon_skips


def generate_therapies(d: Description) -> List[Therapy]:
    """Wrapper around the different therapies"""
    therapies: List[Therapy] = list()
    # Skip a single exon
    therapies += skip_adjacent_exons(d, number_to_skip=1)
    # Skip two adjacent exons
    therapies += skip_adjacent_exons(d, number_to_skip=2)
    return therapies


def _init_model(d: Description) -> None:
    """
    Initialize the HGVS Description

    Don't normalize the positions
    TODO: check that other sanity checks are still performed
    """
    d.to_delins()
    d.de_hgvs_internal_indexing_model = d.delins_model
    d.construct_de_hgvs_internal_indexing_model()
    d.construct_de_hgvs_coordinates_model()
    d.construct_normalized_description()
    d.construct_protein_description()


def _get_genome_annotations(references: Dict[str, Any]) -> Dict[str, Any]:
    """
    The sequence is removed. It should work with conversions, as long as there
    are no sequence slices involved, which will not be the case here.
    """

    def _apply_offset(location: Dict[str, Any], offset: int) -> None:
        if isinstance(location, dict) and location.get("type") == "range":
            if "start" in location and "position" in location["start"]:
                location["start"]["position"] += offset
            if "end" in location and "position" in location["end"]:
                location["end"]["position"] += offset

    def _walk_features(features: List[Dict[str, Any]], offset: int) -> None:
        for feature in features:
            loc = feature.get("location")
            if loc:
                _apply_offset(loc, offset)
            if "features" in feature:
                _walk_features(feature["features"], offset)

    output = {}

    for key, entry in references.items():
        annotations = deepcopy(entry.get("annotations"))
        if not annotations:
            continue

        qualifiers = annotations.get("qualifiers", {})
        offset = qualifiers.pop("location_offset", None)

        if offset is not None:
            _apply_offset(annotations.get("location", {}), offset)
            _walk_features(annotations.get("features", []), offset)

        output[key] = {"annotations": annotations}

    return output


def _description_model(ref_id: str, variants: List[Variant_Dict]) -> Dict[str, Any]:
    """
    To be used only locally with ENSTs.
    """
    return {
        "type": "description_dna",
        "reference": {"id": ref_id, "selector": {"id": ref_id}},
        "coordinate_system": "c",
        "variants": variants,
    }


def _c_variants_to_delins_variants(
    variants: List[Variant_Dict], ref_id: str, references: Dict[str, Any]
) -> List[InternalVariant]:
    """
    The variants can be of any type (substitutions, duplications, etc.).
    """
    model = _description_model(ref_id, variants)
    # logger.debug(f"{model=}")
    delins: List[InternalVariant] = variants_to_delins(
        to_internal_indexing(to_internal_coordinates(model, references))["variants"]
    )
    return delins


def _internal_to_internal_genome(
    variants: List[InternalVariant], offset: int
) -> List[InternalVariant]:
    output = deepcopy(variants)

    for variant in output:
        location = variant.get("location", {})
        if location.get("type") == "range":
            if "start" in location and "position" in location["start"]:
                location["start"]["position"] += offset
            if "end" in location and "position" in location["end"]:
                location["end"]["position"] += offset

    return output


def _get_ensembl_offset(
    references: Dict[str, Any], ref_id: str = "reference"
) -> Union[int, None]:
    offset: Union[int, None] = (
        references.get(ref_id, {})
        .get("annotations", {})
        .get("qualifiers", {})
        .get("location_offset")
    )
    return offset


def changed_protein_positions(reference: str, observed: str) -> List[Tuple[int, int]]:
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


def _cdot_to_internal_delins(
    d: Description, variants: CdotVariant
) -> List[InternalVariant]:
    """Convert a list of cdot variants to internal indels"""
    #  Get stuf we need
    ref_id = get_reference_id(d.corrected_model)

    # Parse the c. string into mutalyzer variant dictionary
    parsed_variants = variant_to_model(variants)
    # logger.debug(f"{parsed_variants=}")

    # Convert the variant dicts into internal delins
    internal_delins = _c_variants_to_delins_variants(
        parsed_variants,
        ref_id,
        d.references,
    )

    # logger.debug(f"{internal_delins=}")
    return internal_delins


def mutation_to_cds_effect(
    d: Description, variants: Sequence[Variant]
) -> List[Tuple[int, int]]:
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

    # Get required data structures from the Description
    ref_id = get_reference_id(d.corrected_model)
    selector_model = get_protein_selector_model(
        d.references[ref_id]["annotations"], ref_id
    )

    # Convert the Variants to their delins model representation
    delins = [v.to_model() for v in variants]

    # Add inverted into the delins model, this is required for the protein prediction
    inverted = selector_model["inverted"]
    for variant in delins:
        if "inserted" in variant and variant["inserted"]:
            variant["inserted"][0]["inverted"] = inverted
        if "deleted" in variant and variant["deleted"]:
            variant["deleted"][0]["inverted"] = inverted

    # Determine the protein positions that were changed
    protein = get_protein_description(delins, d.references, selector_model)
    reference, observed = protein[1], protein[2]

    # Keep track of changed positions on the genome
    changed_genomic = list()

    for start, end in changed_protein_positions(reference, observed):
        # Calculate the nucleotide changed amino acids into a deletion in HGVS c. format
        start_pos = start * 3 + 1
        end_pos = end * 3

        cdot_mutation = CdotVariant(f"{start_pos}_{end_pos}del")

        # Convert cdot to delins
        positions_delins = _cdot_to_internal_delins(d, cdot_mutation)
        ensembl_offset = _get_ensembl_offset(d.references, ref_id)

        # logger.debug(f"{positions_delins=}")
        # logger.debug(f"{ensembl_offset=}")

        if ensembl_offset is None:
            raise RuntimeError("Missing ensemble offset")

        genome_positions = _internal_to_internal_genome(
            positions_delins, ensembl_offset
        )

        assert len(genome_positions) == 1
        g_start = genome_positions[0]["location"]["start"]["position"]
        g_end = genome_positions[0]["location"]["end"]["position"]

        logger.debug(
            f"protein difference: p.{start}-{end}, c.{cdot_mutation}, genomic {g_start:_}-{g_end:_} (size={g_end-g_start})"
        )
        assert g_end > g_start
        changed_genomic.append((g_start, g_end))

    return changed_genomic


def variant_to_model(variant: CdotVariant) -> List[Variant_Dict]:
    """
    Parse the specified variant into a variant model
    """
    results: List[Variant_Dict]
    if "[" in variant:
        results = mutalyzer_hgvs_parser.to_model(variant, "variants")
    else:
        results = [mutalyzer_hgvs_parser.to_model(variant, "variant")]
    return results


def get_exons(
    description: Description, in_transcript_order: bool
) -> List[Tuple[int, int]]:
    """Get exons from a Mutalyzer Description object

    Positions are in the internal coordinate system
    """
    exons: List[Tuple[int, int]] = description.get_selector_model()["exon"]
    if in_transcript_order and description.is_inverted():
        return exons[::-1]

    return exons


def cds_to_internal_positions(
    position: int,
    exons: Sequence[Tuple[int, int]],
    cds_offset: int = 0,
    reverse: bool = False,
) -> int:
    """Convert CDS positions to internal"""
    pos = position + cds_offset

    # Iterate over the exons
    iterator = iter(exons[::-1]) if reverse else iter(exons)

    for start, end in iterator:
        size = end - start
        # If position is in the current exon
        if pos < size:
            return end - pos - 1 if reverse else pos + start
        else:
            pos -= size
    else:
        raise ValueError(f"{position=} outside {exons=}")
