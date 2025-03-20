from copy import deepcopy
from mutalyzer.description import Description, to_rna_reference_model, model_to_string
from mutalyzer.converter.to_hgvs_coordinates import to_hgvs_locations
from mutalyzer.converter.to_internal_coordinates import to_internal_coordinates
from mutalyzer.converter.to_internal_indexing import to_internal_indexing
from mutalyzer.checker import is_overlap
import mutalyzer_hgvs_parser

from pydantic import BaseModel, model_validator

from typing import Any, Tuple, List, Dict

VariantModel = Dict[str, Any]


class HGVS(BaseModel):
    description: str

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "description": "NM_002520.7:c.860_863dup",
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


def HGVS_to_genome_range(hgvs: HGVS) -> Tuple[int, int]:
    """Convert HGVS variant description to affected genome range

    NOTE that the genome range is similar to the UCSC annotations on the genome,
    i.e. 0 based, half open. Not to be confused with hgvs g. positions
    """
    d = Description(description=hgvs.description)
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


def exonskip(hgvs: HGVS) -> List[HGVS]:
    """Generate all possible exon skips for the specified HGVS description"""
    d = Description(description=hgvs.description)
    d.to_delins()

    # Extract relevant information from the normalized description
    raw_response = d.output()
    exons = raw_response["selector_short"]["exon"]["c"]
    transcript_id = raw_response["input_model"]["reference"]["id"]

    exon_skips = list()
    # The first and second exon cannot be skipped
    for start, end in exons[1:-1]:
        skip = f"{transcript_id}:c.{start}_{end}del"
        exon_skips.append(HGVS(description=skip))
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


def mutation_to_cds_effect(hgvs: HGVS) -> Tuple[int, int]:
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
    d = Description(description=hgvs.description)
    # Do only part of the normalizations from mutalyzer
    _init_model(d)

    # Determine the protein positions that were changed
    protein = d.output()["protein"]
    first = protein["position_first"]
    last = protein["position_last_original"]

    # Convert the changed amino acids into a deletion in HGVS c. format
    transcript_id = hgvs.description.split(":c.")[0]
    start_pos = first * 3
    end_pos = (last * 3) - 1

    cdot = f"{transcript_id}:c.{start_pos}_{end_pos}del"

    return HGVS_to_genome_range(HGVS(description=cdot))


def variant_to_model(variant: str) -> List[VariantModel]:
    """
    Parse the specified variant into a variant model
    """
    results: List[VariantModel]
    if "[" in variant:
        results = mutalyzer_hgvs_parser.to_model(variant, "variants")
    else:
        results = [mutalyzer_hgvs_parser.to_model(variant, "variant")]
    return results


def append_mutation(description: Description, mutation: str) -> None:
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
