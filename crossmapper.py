from gtgt.mutalyzer import init_description, sequence_from_description, protein_prediction
from gtgt.variant import Variant
from mutalyzer.description_model import get_reference_id
from mutalyzer.description import Description
from mutalyzer_crossmapper import Coding
import json
from gtgt.mutalyzer import changed_protein_positions

def pprint(thing):
    print(json.dumps(thing, indent=True))



transcripts = {
    "WT1": [
        "ENST00000452863.10",# "NC_000011.10(NM_024426.6)", "NM_024426.6"
    ],
    "SDHD": [
        "ENST00000375549.8", "NC_000011.10(NM_003002.4)", "NM_003002.4"
    ]
}

def get_offset(d: Description) -> int:
    ref_id = get_reference_id(d.corrected_model)
    offset: int = (
        d.references.get(ref_id, {})
        .get("annotations", {})
        .get("qualifiers", {})
        .get("location_offset", 0)
    )
    return offset

def genomic_crossmapper(hgvs: str) -> Coding:
    """Create a genomic crossmapper for the hgvs description"""
    # First, we re-write the hgvs to c. to ensure we have the introns
    h  = hgvs.replace(":r.", ":c.")
    d = init_description(h)

    # Get the offset of the exons
    offset = get_offset(d)
    # print(f"{offset}")

    # Get the exons on the genome
    exons = d.get_selector_model()["exon"]
    # print(exons)
    exons = [(start + offset, end+offset) for start,end in exons]
    # print(exons)

    # Get the cds on the genome
    cds = d.get_selector_model()["cds"]
    # print(f"{cds=}")
    cds_start, cds_end = cds[0]
    cds_start += offset
    cds_end += offset
    cds = cds_start, cds_end

    inverted= d.is_inverted()
    print(f"Coding({exons=},{cds=},{inverted=})")
    return Coding(exons, cds, inverted)

def transcript_crossmapper(hgvs: str) -> Coding:
    """Create a genomic crossmapper for the hgvs description"""
    # First, we re-write the hgvs to c. to ensure we have the introns
    h  = hgvs.replace(":r.", ":c.")
    d = init_description(h)

    # Get the exons on the genome
    exons = d.get_selector_model()["exon"]
    # print(exons)
    exons = [(start, end) for start,end in exons]
    # print(exons)

    # Get the cds on the genome
    cds = d.get_selector_model()["cds"]
    # print(f"{cds=}")
    cds_start, cds_end = cds[0]
    cds = cds_start, cds_end

    inverted= d.is_inverted()
    print(f"Coding({exons=},{cds=},{inverted=})")
    return Coding(exons, cds, inverted)

def variants_from_protein(d):
    """
    Re-create the variants from the protein description
    """
    gtgt_variants = list()
    # crossmap = genomic_crossmapper(d.input_description)

    transcript_crossmap = transcript_crossmapper(d.input_description)

    # Get the changed amino acids
    variants = [Variant.from_model(delins, sequence=sequence)
    for delins in d.delins_model["variants"]
    ]
    protein = protein_prediction(d, variants)
    reference, observed = protein[1], protein[2]
    for protein_start, protein_end in changed_protein_positions(reference, observed):
        print(f"Protein change: ({protein_start}, {protein_end})")
        # First, we map the protein to the coordinate
        start = transcript_crossmap.protein_to_coordinate((protein_start+ 1, 1, 0, 0, 0))
        end = transcript_crossmap.protein_to_coordinate((protein_end, 3, 0, 0, 0)) + 1

        print(f"Protein coordinate: ({start}, {end})")

        print(f"Coding position: ({start}, {end})")
        # Next, we map the coordinate to the noncoding, which is what mutalyzer
        # uses as offset

        # start, offset, upstream = crossmap.coordinate_to_noncoding(start)
        #
        # if offset or upstream:
        #     raise RuntimeError(f"{d.input_description} is outside the CDS")
        #
        # end, offset, upstream = crossmap.coordinate_to_noncoding(end)
        #
        # if offset or upstream:
        #     raise RuntimeError(f"{d.input_description} is outside the CDS")
        
        if end < start:
            start, end = end - 1, start + 1

        gtgt_variants.append(Variant(start, end))
    return gtgt_variants
    print(crossmap)
    

# Special variants which are deletions which are exactly aligned with amino acids
# This is needed to ensure they match the variants we derive from the protein changes
WT1_variants = "10_12del 658_663del 701_721del 1253_1261del".split()
WT1_variants = "10_12del 658_663del".split()

for transcript in transcripts["WT1"]:
    for coordinate in ["c", "r"]:
        for variant in WT1_variants:
            hgvs = f"{transcript}:{coordinate}.{variant}"
            d = init_description(hgvs)

            # Determine the Variants from the mutalyzer description
            sequence = sequence_from_description(d)
            mutalyzer_variants = [Variant.from_model(delins, sequence=sequence)
            for delins in d.delins_model["variants"]
            ]
            print(hgvs)
            print(mutalyzer_variants)

            # Determine the Variants by going via the genomic crossmapper and
            # the differences in the protein description to reconsitute the
            # variants ourself, and then map them back to the coordinates that
            # mutalyzer uses
            gtgt_variants = variants_from_protein(d)
            print(gtgt_variants)
            print()
