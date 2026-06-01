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
from typing_extensions import NewType

logger = logging.getLogger(__name__)

from .variant import Variant, combine_variants_deletion

# Mutalyzer Variant dictionary
Variant_Dict = NewType("Variant_Dict", Mapping[str, Any])


def chrom_to_nc(chrom: str) -> str:
    """Convert chromosome name (e.g. 11, chr11) to NC

    If the NC cannot be found, it returns chr
    """
    mapping = {
        "chr1": "NC_000001.11",
        "chr2": "NC_000002.12",
        "chr3": "NC_000003.12",
        "chr4": "NC_000004.12",
        "chr5": "NC_000005.10",
        "chr6": "NC_000006.12",
        "chr7": "NC_000007.14",
        "chr8": "NC_000008.11",
        "chr9": "NC_000009.12",
        "chr10": "NC_000010.11",
        "chr11": "NC_000011.10",
        "chr12": "NC_000012.12",
        "chr13": "NC_000013.11",
        "chr14": "NC_000014.9",
        "chr15": "NC_000015.10",
        "chr16": "NC_000016.10",
        "chr17": "NC_000017.11",
        "chr18": "NC_000018.10",
        "chr19": "NC_000019.10",
        "chr20": "NC_000020.11",
        "chr21": "NC_000021.9",
        "chr22": "NC_000022.11",
        "chrX": "NC_000023.11",
        "chrY": "NC_000024.10",
        "chrMT": "NC_012920.1",
        "1": "NC_000001.11",
        "2": "NC_000002.12",
        "3": "NC_000003.12",
        "4": "NC_000004.12",
        "5": "NC_000005.10",
        "6": "NC_000006.12",
        "7": "NC_000007.14",
        "8": "NC_000008.11",
        "9": "NC_000009.12",
        "10": "NC_000010.11",
        "11": "NC_000011.10",
        "12": "NC_000012.12",
        "13": "NC_000013.11",
        "14": "NC_000014.9",
        "15": "NC_000015.10",
        "16": "NC_000016.10",
        "17": "NC_000017.11",
        "18": "NC_000018.10",
        "19": "NC_000019.10",
        "20": "NC_000020.11",
        "21": "NC_000021.9",
        "22": "NC_000022.11",
        "X": "NC_000023.11",
        "Y": "NC_000024.10",
        "MT": "NC_012920.1",
    }

    return mapping.get(chrom, chrom)


CHROMOSOME_ASSEMBLY = {
    "NC_000001.11": "GRCh38",
    "NC_000002.12": "GRCh38",
    "NC_000003.12": "GRCh38",
    "NC_000004.12": "GRCh38",
    "NC_000005.10": "GRCh38",
    "NC_000006.12": "GRCh38",
    "NC_000007.14": "GRCh38",
    "NC_000008.11": "GRCh38",
    "NC_000009.12": "GRCh38",
    "NC_000010.11": "GRCh38",
    "NC_000011.10": "GRCh38",
    "NC_000012.12": "GRCh38",
    "NC_000013.11": "GRCh38",
    "NC_000014.9": "GRCh38",
    "NC_000015.10": "GRCh38",
    "NC_000016.10": "GRCh38",
    "NC_000017.11": "GRCh38",
    "NC_000018.10": "GRCh38",
    "NC_000019.10": "GRCh38",
    "NC_000020.11": "GRCh38",
    "NC_000021.9": "GRCh38",
    "NC_000022.11": "GRCh38",
    "NC_000023.11": "GRCh38",
    "NC_000024.10": "GRCh38",
    "NC_012920.1": "GRCh38",
}


def sequence_from_description(d: Description) -> str:
    """Return the sequence form a description"""
    _id = d.input_model["reference"]["id"]
    sequence: str = d.references[_id]["sequence"]["seq"]
    return sequence


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


def get_assembly_name(d: Description) -> str:
    """Extract the assembly name from a Description"""
    qualifiers = d.references["reference"]["annotations"]["qualifiers"]
    # ENS
    if "assembly_name" in qualifiers:
        return str(qualifiers["assembly_name"])

    # NC(NM)
    id = d.references["reference"]["annotations"]["id"]
    if id in CHROMOSOME_ASSEMBLY:
        return CHROMOSOME_ASSEMBLY[id]

    # NM
    raise ValueError(f"Unable to determine assembly for {d}")


def get_transcript_name(d: Description) -> str:
    """Extract the transcript name from a description"""
    name = d.get_selector_model()["id"]
    if name is None:
        raise ValueError(f"Unable to determine transcript name for {d}")
    else:
        return str(name)


def get_chrom_name(d: Description) -> str:
    # Try ensembl chromosome
    chrom = _get_ensembl_chrom_name(d)
    if chrom is not None:
        chrom = chrom_to_nc(chrom)
        if chrom is not None:
            return chrom

    # Try NCBI chromosome
    chrom = _get_ncbi_chrom_name(d)
    if chrom is not None:
        chrom = chrom_to_nc(chrom)
        if chrom is not None:
            return chrom

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
