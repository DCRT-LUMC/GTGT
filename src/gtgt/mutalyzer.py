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

from typing import Any, Tuple, List, Dict, Union
from typing_extensions import NewType


# Mutalyzer variant object, using the 'internal' coordinate system (0 based, half open)
# Variant string in HGVS c. format
CdotVariant = NewType("CdotVariant", str)
# Mutalyzer Variant dictionary
Variant = NewType("Variant", Dict[str, Any])
InternalVariant = NewType("InternalVariant", dict[str, Any])


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
        s_var = self.description.split("c.")[1]
        o_var = other.description.split("c.")[1]

        # If self is a deletion, and other is fully inside self, we don't have to add anything
        if s_type == "deletion":
            if o_start >= s_start and o_end <= s_end:
                return

        # If self is not a deletion:
        # self is before other
        if s_end < o_start:
            self.description = f"{s_id}:c.[{s_var};{o_var}]"
        # self is after other
        elif s_start > o_end:
            self.description = f"{s_id}:c.[{o_var};{s_var}]"
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


def exonskip(d: Description) -> List[Therapy]:
    """Generate all possible exon skips for the specified HGVS description"""
    d.to_delins()

    # Extract relevant information from the normalized description
    raw_response = d.output()
    exons = raw_response["selector_short"]["exon"]["c"]
    transcript_id = raw_response["input_model"]["reference"]["id"]

    exon_skips = list()
    # The first and second exon cannot be skipped

    exon_counter = 2
    for start, end in exons[1:-1]:
        name = f"Skip exon {exon_counter}"
        hgvs = f"{transcript_id}:c.{start}_{end}del"
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

    # Convert the variant dicts into internal delins
    internal_delins = _c_variants_to_delins_variants(
        parsed_variants, ref_id, d.references
    )
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

        if ensembl_offset is None:
            raise RuntimeError("Missing ensemble offset")

        genome_positions = _internal_to_internal_genome(
            positions_delins, ensembl_offset
        )

        assert len(genome_positions) == 1
        start = genome_positions[0]["location"]["start"]["position"]
        end = genome_positions[0]["location"]["end"]["position"]

        assert end > start
        changed_genomic.append((start, end))

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
