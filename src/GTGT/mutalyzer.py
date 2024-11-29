from mutalyzer.description import Description, to_rna_reference_model, model_to_string
from mutalyzer.converter.to_hgvs_coordinates import to_hgvs_locations
from .models import HGVS

from typing import Tuple, List


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
    d.normalize()

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
