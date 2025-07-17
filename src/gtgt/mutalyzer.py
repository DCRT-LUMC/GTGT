from copy import deepcopy
import dataclasses

from mutalyzer.description import Description, to_rna_reference_model, model_to_string
from mutalyzer.converter.to_hgvs_coordinates import to_hgvs_locations
from mutalyzer.converter.to_delins import variants_to_delins
from mutalyzer.converter.to_internal_coordinates import to_internal_coordinates
from mutalyzer.converter.to_internal_indexing import to_internal_indexing
from mutalyzer.description_model import get_reference_id, variants_to_description
from mutalyzer.protein import get_protein_description
from mutalyzer.reference import get_protein_selector_model
from mutalyzer.checker import is_overlap
import mutalyzer_hgvs_parser

from pydantic import BaseModel, model_validator

import Levenshtein

from typing import Any, Tuple, List, Dict, Union, Sequence
from typing_extensions import NewType

import logging

logger = logging.getLogger(__name__)

# Mutalyzer variant object, using the 'internal' coordinate system (0 based, half open)
# Variant string in HGVS c. format
CdotVariant = NewType("CdotVariant", str)
# Mutalyzer Variant dictionary
Variant = NewType("Variant", Dict[str, Any])
InternalVariant = NewType("InternalVariant", dict[str, Any])


class _Variant:
    """Class to store delins variants"""

    def __init__(self, start: int, end: int, inserted: str = ""):
        if start > end:
            raise ValueError(f"End ({end}) must be after start ({start})")
        self.start = start  # zero based
        self.end = end  # exclusive
        self.inserted = inserted

    def __str__(self) -> str:
        return self.__repr__()

    def __repr__(self) -> str:
        start = self.start
        end = self.end
        sequence = self.inserted
        return f"Variant({start=}, {end=}, sequence={sequence})"

    def before(self, other: "_Variant") -> bool:
        return self.end <= other.start

    def after(self, other: "_Variant") -> bool:
        return self.start >= other.end

    def inside(self, other: "_Variant") -> bool:
        return self.start >= other.start and self.end <= other.end

    def overlap(self, other: "_Variant") -> bool:
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
        if not isinstance(other, _Variant):
            raise NotImplementedError
        return (
            self.start == other.start
            and self.end == other.end
            and self.inserted == other.inserted
        )

    def __lt__(self, other: "_Variant") -> bool:
        if not isinstance(other, _Variant):
            raise NotImplementedError
        if self.overlap(other):
            msg = f"Overlapping variants '{self}' and '{other}' cannot be sorted"
            raise ValueError(msg)

        return self.start < other.start

    @classmethod
    def from_model(cls, model: Dict[str, Any]) -> "_Variant":
        start = model["location"]["start"]["position"]
        end = model["location"]["end"]["position"]
        # sequence is string or empty list
        inserted = model["inserted"][0]["sequence"]
        return _Variant(start, end, inserted if inserted else "")

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

        return {
            "location": location,
            "type": "deletion_insertion",
            "source": "reference",
            "inserted": inserted,
        }


@dataclasses.dataclass
class Therapy:
    """Class to store genetic therapies"""

    name: str
    hgvs: str
    description: str


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

    @staticmethod
    def _validate_for_apply_deletion(hgvs: "HGVS") -> None:
        """
        Raise a NotImplementedError if the hgvs is not supported
        """
        model = mutalyzer_hgvs_parser.to_model(hgvs.description)

        # There must be one variant
        if len(model["variants"]) != 1:
            raise NotImplementedError

        var = model["variants"][0]

        # The variant must be in the CDS
        if "outside_cds" in var["location"]:
            raise NotImplementedError

        # The variant must not be intronic
        if "offset" in var["location"]:
            raise NotImplementedError

    @property
    def position(self) -> Tuple[int, int]:
        """
        Return the position of a description as (start, end)

        These are just the .c position, so 1 based and inclusive
        """
        model = mutalyzer_hgvs_parser.to_model(self.description)
        assert len(model["variants"]) == 1

        var = model["variants"][0]
        if var["location"]["type"] == "point":
            p = var["location"]["position"]
            return p, p
        elif var["location"]["type"] == "range":
            s = var["location"]["start"]["position"]
            e = var["location"]["end"]["position"]
            return s, e
        else:
            raise NotImplementedError

    def apply_deletion(self, other: "HGVS") -> None:
        """
        Apply a deletion to the current variant

        If the deletion does not overlap, add them together
        If the deletion completely overlaps the variant, replace the variant

        If the deletion partially overlaps the variant, raise an error
        """
        logger.debug(f"Applying {other} to {self}")
        coordinate_system = other.description.split(":")[1][:2]
        # Perform all validations
        self._validate_for_apply_deletion(other)
        self._validate_for_apply_deletion(self)

        # other must be a deletion
        o_model = mutalyzer_hgvs_parser.to_model(other.description)
        assert len(o_model["variants"]) == 1

        o_type = o_model["variants"][0]["type"]
        if o_type != "deletion":
            raise NotImplementedError

        s_model = mutalyzer_hgvs_parser.to_model(self.description)
        assert len(s_model["variants"]) == 1
        s_type = s_model["variants"][0]["type"]

        # self and other must refer to the same reference ID
        s_id = s_model["reference"]["id"]
        o_id = o_model["reference"]["id"]
        if s_id != o_id:
            raise NotImplementedError

        # Get the c. positions for start and end
        s_start, s_end = self.position
        o_start, o_end = other.position

        # Get the variants in text format
        s_var = self.description.split(coordinate_system)[1]
        o_var = other.description.split(coordinate_system)[1]

        logger.debug(f"{self.position=}, {other.position=}")

        # If self is a deletion, and other is fully inside self, we don't have to add anything
        if s_type == "deletion" and o_start >= s_start and o_end <= s_end:
            return
        elif s_type == "insertion":
            # if other is before the insertion
            if o_end <= s_start:
                self.description = f"{s_id}:{coordinate_system}[{o_var};{s_var}]"
            # if other is after the insertion
            elif o_start >= s_end:
                self.description = f"{s_id}:{coordinate_system}[{s_var};{o_var}]"
            # other overlaps the insertion site
            if o_start <= s_start and o_end >= s_end:
                self.description = f"{s_id}:{coordinate_system}{o_var}"
        # TODO: implement handling of a deletion partially overlapping an indel
        # elif s_type == "deletion_insertion":
        else:
            # self is before other
            if s_end < o_start:
                self.description = f"{s_id}:{coordinate_system}[{s_var};{o_var}]"
            # self is after other
            elif s_start > o_end:
                self.description = f"{s_id}:{coordinate_system}[{o_var};{s_var}]"
            # self is fully inside other
            elif s_start >= o_start and s_end <= o_end:
                # We overwrite self with other
                self.description = other.description
            # partial overlaps are not supported
            else:
                msg = f"Unable to apply deletion {other} to {self}"
                raise NotImplementedError(msg)


def HGVS_to_genome_range(d: Description) -> Tuple[int, int]:
    """Convert HGVS variant description to affected genome range

    NOTE that the genome range is similar to the UCSC annotations on the genome,
    i.e. 0 based, half open. Not to be confused with hgvs g. positions
    """
    d.normalize()
    model = d.ensembl_model_with_no_offset()

    if len(model["variants"]) == 0:
        raise ValueError("Descriptions without variants are not supported")
    if len(model["variants"]) > 1:
        raise ValueError("Multiple variants are not supported")

    # Get start and end of the description
    start = model["variants"][0]["location"]["start"]["position"]
    end = model["variants"][0]["location"]["end"]["position"]

    return (start, end)


def combine_variants_deletion(
    variants: Sequence[_Variant], deletion: _Variant
) -> List[_Variant]:
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


def to_cdot_hgvs(d: Description, variants: Sequence[_Variant]) -> str:
    """Convert a list of _Variants to hgvs representation"""
    # Convert to delins dict model
    variant_models = [v.to_model() for v in variants]

    print(f"{variant_models=}")
    # Convert 'delins' without insertion to deletion
    for model in variant_models:
        model["type"] = "substitution"
        # continue
        if not model.get("inserted"):
            model["type"] = "deletion"

    print(f"{variant_models=}")

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
    
    hgvs: str =variants_to_description(cdot_locations)
    return hgvs


def _exonskip(d: Description) -> List[Therapy]:
    """Generate all possible exon skips for the specified Description"""
    exon_skips = list()

    exons = get_exons(d, in_transcript_order=True)
    variants = [_Variant.from_model(v) for v in d.internal_coordinates_model["variants"]]

    exon_counter = 2
    for start, end in exons[1:-1]:
        exon_skip = _Variant(start, end)
        # Combine the existing variants with the exon skip
        combined = combine_variants_deletion(variants, exon_skip)

        # Convert to c. notation (user facing)
        name = f"Skip exon {exon_counter}"
        hgvs = ""
        description = f"The annotations based on the supplied variants, in combination with skipping exon {exon_counter}."
        t = Therapy(name, hgvs, description)
        exon_skips.append(t)
        exon_counter += 1

    return list()


def exonskip(d: Description) -> List[Therapy]:
    """Generate all possible exon skips for the specified HGVS description"""
    d.to_delins()
    coordinate_system = f"{d.input_model['coordinate_system']}."

    # Extract relevant information from the normalized description
    raw_response = d.output()
    exons = raw_response["selector_short"]["exon"]["c"]
    transcript_id = raw_response["input_model"]["reference"]["id"]

    exon_skips = list()
    # The first and second exon cannot be skipped

    exon_counter = 2
    for start, end in exons[1:-1]:
        name = f"Skip exon {exon_counter}"
        hgvs = f"{transcript_id}:{coordinate_system}{start}_{end}del"
        description = f"The annotations based on the supplied variants, in combination with skipping exon {exon_counter}."
        t = Therapy(name, hgvs, description)
        exon_skips.append(t)
        exon_counter += 1
    return exon_skips


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


def _description_model(ref_id: str, variants: List[Variant]) -> Dict[str, Any]:
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
    variants: List[Variant], ref_id: str, references: Dict[str, Any]
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
    d: Description, variants: CdotVariant
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
    # Convert the c. variants to internal indels
    delins = _cdot_to_internal_delins(d, variants)

    # Get required data structures from the Description
    ref_id = get_reference_id(d.corrected_model)
    selector_model = get_protein_selector_model(
        d.references[ref_id]["annotations"], ref_id
    )

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


def variant_to_model(variant: CdotVariant) -> List[Variant]:
    """
    Parse the specified variant into a variant model
    """
    results: List[Variant]
    if "[" in variant:
        results = mutalyzer_hgvs_parser.to_model(variant, "variants")
    else:
        results = [mutalyzer_hgvs_parser.to_model(variant, "variant")]
    return results


def append_mutation(description: Description, mutation: CdotVariant) -> None:
    """
    Add mutation to the Description, re-using the Description object
    """
    # Get the variant model in c.
    c_variants = variant_to_model(mutation)

    # Convert the c. variant to i.
    model = deepcopy(description.corrected_model)
    # Add the c_variant to the variant(s) which are already there
    model["variants"] += c_variants
    model = to_internal_coordinates(model, description.references)
    model = to_internal_indexing(model)

    if is_overlap(model["variants"]):
        msg = f"Variant {mutation} overlaps {description.input_description}"
        raise ValueError(msg)

    # Replace the variant in the description
    description.de_hgvs_internal_indexing_model["variants"] = model["variants"]

    # Update the internal description models
    description.construct_de_hgvs_coordinates_model()
    description.construct_normalized_description()
    description.construct_protein_description()


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
