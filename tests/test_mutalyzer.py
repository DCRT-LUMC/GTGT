import pytest
from pathlib import Path
from gtgt.mutalyzer import (
    CdotVariant,
    HGVS_to_genome_range,
    InternalVariant,
    append_mutation,
    exonskip,
    mutation_to_cds_effect,
    mutation_to_cds_effect2,
    _cdot_to_internal_delins,
    _init_model,
)
from gtgt.mutalyzer import HGVS, VariantModel, variant_to_model
from mutalyzer.description import Description
from mutalyzer.converter.to_delins import variants_to_delins
from mutalyzer.converter.to_internal_coordinates import to_internal_coordinates
from mutalyzer.converter.to_internal_indexing import to_internal_indexing
from mutalyzer.description_model import get_reference_id, variants_to_description
from gtgt.transcript import Transcript
from gtgt.models import TranscriptModel
import json
import copy

from itertools import zip_longest
from typing import Any, List, Tuple, Dict


def retrieve_raw(
    reference_id: str,
    reference_source: Any = None,
    reference_type: Any = None,
    size_off: bool = True,
    configuration_path: Any = None,
    timeout: int = 1,
) -> Tuple[str, str, str]:
    if reference_type == "fasta":
        return _get_content("data/" + reference_id + ".fasta"), "fasta", "ncbi"
    elif reference_id.startswith("LRG_"):
        return _get_content("data/" + reference_id), "lrg", "lrg"
    else:
        return _get_content("data/" + reference_id + ".gff3"), "gff3", "ncbi"


def get_cds_to_mrna(cds_id: Any, timeout: Any = 10) -> None:
    return None


def _get_content(relative_location: str) -> str:
    data_file = Path(__file__).parent.joinpath(relative_location)
    try:
        with open(str(data_file), "r") as file:
            content = file.read()
    except FileNotFoundError:
        raise RuntimeError({}, [])
    return content


@pytest.fixture(autouse=True)
def mock_env(monkeypatch: Any) -> None:
    monkeypatch.setattr("mutalyzer_retriever.retriever.retrieve_raw", retrieve_raw)
    monkeypatch.setattr("mutalyzer.description.get_cds_to_mrna", get_cds_to_mrna)


POSITIONS = [
    # c., genome range
    # Start of the transcript / first exon
    ("-35", (0, 1)),
    # Start of the coding region
    ("1", (35, 36)),
    # End of the first exon
    ("52", (86, 87)),
    # Start of the second exon
    ("53", (984, 985)),
    # End of the second exon
    ("169", (1100, 1101)),
    # Start of the third exon
    ("170", (1994, 1995)),
    # End of the third exon
    ("314", (2138, 2139)),
    # Start of the fourth exon
    ("315", (7932, 7933)),
    # End of the CDS
    ("480", (8097, 8098)),
    # End of the transcript / fourth exon
    ("*824", (8921, 8922)),
]


@pytest.mark.parametrize("cdot, genomic", POSITIONS)
def test_genomic_to_c_dot_SDHD(cdot: str, genomic: Tuple[int, int]) -> None:
    """
    GIVEN a Transcript and a position in c. notation
    WHEN we translate the c. to genomic positions
    THEN we should get a range, on the genomic coordinates

    NOTE that genomic coordinates refer to the UCSC annotations on the genome,
    i.e. 0 based, half open. Not to be confused with hgvs g. positions
    """
    # Ensemble transcript ID for SDHD
    ENST = "ENST00000375549.8"
    # Build the HGVS description
    description = f"{ENST}:c.{cdot}"

    d = Description(description)
    _init_model(d)

    start, end = HGVS_to_genome_range(d)

    # Genomic offset of SDHD, just used for testing so the genomic positions
    # are manageable
    offset = 112086872
    start -= offset
    end -= offset

    assert (start, end) == genomic


POSITIONS = [
    # c., genome range
    # Start of the transcript / first exon
    ("-179", (47764, 47765)),
    # Start of the coding region
    ("1", (47585, 47586)),
    # End of the first exon
    ("661", (46925, 46926)),
    # Start of the second exon
    ("662", (40844, 40845)),
    # End of the second exon
    ("784", (40722, 40723)),
    # Start of the third exon
    ("785", (40283, 40284)),
    # End of the third exon
    ("887", (40181, 40182)),
    # Start of the tenth exon
    ("1448", (1404, 1405)),
    # End of the CDS
    ("1569", (1283, 1284)),
    # End of the transcript / tenth exon
    ("*1283", (0, 1)),
]


@pytest.mark.parametrize("cdot, genomic", POSITIONS)
def test_genomic_to_c_dot_WT1(cdot: str, genomic: Tuple[int, int]) -> None:
    """
    GIVEN a Transcript and a position in c. notation
    WHEN we translate the c. to genomic positions
    THEN we should get a range, on the genomic coordinates

    NOTE: WT1 is on the reverse strand

    NOTE that genomic coordinates refer to the UCSC annotations on the genome,
    i.e. 0 based, half open. Not to be confused with hgvs g. positions
    """
    # Ensemble transcript ID for SDHD
    ENST = "ENST00000452863.10"
    # Build the HGVS description
    description = f"{ENST}:c.{cdot}"

    d = Description(description)
    _init_model(d)

    start, end = HGVS_to_genome_range(d)

    # Genomic offset of WT1, just used for testing so the genomic positions are
    # manageable
    offset = 32387774
    start -= offset
    end -= offset

    assert (start, end) == genomic


SUPPORTED_DESCRIPTIONS = [
    # SNP
    ("ENST00000375549.8:c.5C>T"),
    # Deletion of range
    ("ENST00000375549.8:c.5_15del"),
    # Inversion
    ("ENST00000375549.8:c.5_15inv"),
    # Empty position
    ("ENST00000375549.8:c.5"),
    # Empty range
    ("ENST00000375549.8:c.5_10"),
    # Insertion of a range
    ("ENST00000375549.8:c.10_11ins50_60"),
]


@pytest.mark.parametrize("description", SUPPORTED_DESCRIPTIONS)
def test_supported_descriptions(description: str) -> None:
    d = Description(description)
    _init_model(d)
    HGVS_to_genome_range(d)


UNSUPPORTED_DESCRIPTIONS = [
    # Description with multiple variants
    ("ENST00000375549.8:c.[5C>T;10del]"),
    # Empty HGVS description
    ("ENST00000375549.8:c.="),
    # Delins of a range
    ("ENST00000375549.8:c.10_15delins50_60"),
]


@pytest.mark.parametrize("description", UNSUPPORTED_DESCRIPTIONS)
def test_unsupported_descriptions(description: str) -> None:
    with pytest.raises(ValueError):
        d = Description(description)
        _init_model(d)
        HGVS_to_genome_range(d)


def test_exonskip_SDHD() -> None:
    d = Description("ENST00000375549.8:c.=")
    _init_model(d)
    results = [
        "ENST00000375549.8:c.53_169del",
        "ENST00000375549.8:c.170_314del",
    ]
    for output, expected in zip_longest(exonskip(d), results):
        assert output == HGVS(description=expected)


def test_exonskip_WT1() -> None:
    d = Description("ENST00000452863.10:c.=")
    _init_model(d)
    results = [
        "ENST00000452863.10:c.662_784del",
        "ENST00000452863.10:c.785_887del",
        "ENST00000452863.10:c.888_965del",
        "ENST00000452863.10:c.966_1016del",
        "ENST00000452863.10:c.1017_1113del",
        "ENST00000452863.10:c.1114_1264del",
        "ENST00000452863.10:c.1265_1354del",
        "ENST00000452863.10:c.1355_1447del",
    ]
    for output, expected in zip_longest(exonskip(d), results):
        assert output == HGVS(description=expected)


MUTATIONS = [
    # HGVS, coordinates on the genome
    # A simple missense that changes a single amino acids
    ("ENST00000452863.10:c.13T>A", (32435347, 32435350)),
    # A frameshift which destroys most of the protein
    ("ENST00000452863.10:c.10del", (32389058, 32435352)),
    # A frameshift that is restored by an insertion
    ("ENST00000452863.10:c.[10del;20_21insA]", (32435340, 32435352)),
    # A frameshift that is restored by a bigger insertion
    ("ENST00000452863.10:c.[10del;20_21insATCGAATATGGGG]", (32435340, 32435352)),
    # A bigger deletion
    ("ENST00000452863.10:c.11_19del", (32435344, 32435353)),
    # An inframe deletion that creates a STOP codon
    ("ENST00000452863.10:c.87_89del", (32389059, 32435278)),
]


@pytest.mark.parametrize("description, expected", MUTATIONS)
def test_mutation_to_cds_effect(description: str, expected: Tuple[int, int]) -> None:
    """
    GIVEN a HGVS transcript description
    WHEN we determine the CDS effect
    THEN we should get genome coordinates
    """
    d = Description(description)
    _init_model(d)

    assert mutation_to_cds_effect(d) == expected


MUTATIONS2 = [
    # HGVS, coordinates on the genome,
    # A simple missense that changes a single amino acids
    ("13T>A", (32435345, 32435348)),
    # A frameshift which destroys most of the protein
    ("10del", (32389057, 32435351)),
    # # A frameshift that is restored by an insertion
    ("[10del;20_21insA]", (32435339, 32435351)),
    # # A frameshift that is restored by a bigger insertion
    ("[10del;20_21insATCGAATATGGGG]", (32435339, 32435351)),
    # # A bigger deletion
    ("11_19del", (32435342, 32435351)),
    # # An inframe deletion that creates a STOP codon
    ("87_89del", (32389057, 32435276)),
]


@pytest.mark.parametrize("variants, expected", MUTATIONS2)
def test_mutation_to_cds_effect2(
    variants: CdotVariant, expected: Tuple[int, int]
) -> None:
    """
    GIVEN a HGVS transcript description for a transcript on the reverse strand
    WHEN we determine the CDS effect
    THEN we should get genome coordinates
    """
    d = Description("ENST00000452863.10:c.=")
    _init_model(d)

    print()
    print(f"{variants=}")
    assert mutation_to_cds_effect2(d, variants) == expected


CDOT_MUTATIONS = [
    # c. variant, i. delins variant
    ("13T>A", "47_48delTinsA"),
    ("-35_52del", "0_87delins"),
    ("53_169del", "984_1101delins"),
    ("315_*824del", "7932_8922delins"),
    # Multiple variants
    ("[13T>A;15_16delinsAT]", "[47_48delTinsA;49_51delinsAT]"),
]


@pytest.mark.parametrize("cdot, internal_delins", CDOT_MUTATIONS)
def test_cdot_to_indel(cdot: str, internal_delins: str) -> None:
    """
    GIVEN a list of c. mutations a string
    WHEN we convert them to internal i. indels
    THEN we should get a list of internal variants
    """
    d = Description("ENST00000375549.8:c.=")
    _init_model(d)
    indels = _cdot_to_internal_delins(d, cdot)
    assert variants_to_description(indels) == internal_delins


@pytest.fixture
def WT() -> Transcript:
    """
    Transcript for WT1, using real genomic positionsjjkkkkjklkj
    """
    path = "tests/data/ENST00000452863.10.Transcript.json"
    with open(path) as fin:
        js = json.load(fin)

    t = TranscriptModel.model_validate(js)

    return t.to_transcript()


def test_analyze_transcript(WT: Transcript) -> None:
    # In frame deletion that creates a STOP codon
    # variant = "ENST00000452863.10:c.87_89del"
    # Frameshift in small in-frame exon 5
    variant = "ENST00000452863.10:c.970del"

    results = WT.analyze(variant)

    assert results["wildtype"]["cds"]["percentage"] == 1.0


VARIANTS = [
    # variant, Transcript effect
    # In frame deletion that creates a STOP codon
    ("ENST00000452863.10:c.87_89del", 0.0018),
    # In frame deletion that does not make a STOP
    ("ENST00000452863.10:c.85_87del", 0.9999),
    # Synonymous mutation
    ("ENST00000452863.10:c.13T>C", 1),
]


@pytest.mark.parametrize("variant, effect", VARIANTS)
def test_mutate_transcript_with_variant(
    variant: str, effect: float, WT: Transcript
) -> None:

    d = Description(variant)
    _init_model(d)

    modified = copy.deepcopy(WT)
    modified.mutate(d)

    cmp = modified.compare(WT)
    assert cmp["cds"]["percentage"] == pytest.approx(effect, abs=0.0001)


APPEND_VARIANT = [
    # mutation, predicted protein description as readout
    ("10del", "MAVSGG*")
]


@pytest.mark.parametrize("mutation, protein", APPEND_VARIANT)
def test_append_mutation(mutation: str, protein: str) -> None:
    """
    GIVEN a string denoting a HGVS variant
    WHEN we append this variant to an existing Description
    THEN the Description should be updated
    """
    transcript = "ENST00000375549.8:c.="
    d = Description(f"{transcript}")
    _init_model(d)

    append_mutation(d, mutation)

    assert d.output()["protein"]["predicted"] == protein


APPEND_TO_EXISTING = [
    # Add an insertion after the deletion to restore the reading frame
    ("10del", "11_12insA", "MAVYWRLSAV"),
    # Add an insertion before the deletion, e.g. the variants are out of order after appending
    ("10del", "4dup", "MGGFWRLSAV"),
    # Existing variant is a snp, not deletion
    ("10C>T", "11dup", "MAVFLEAECR"),
    # Existing variant is a dup
    ("10dup", "15del", "MAVPLRLSAV"),
    # Existing variant is an inversion
    ("10_15inv", "16A>T", "MAVPEWLSAV"),
    # Existing variant is an delins
    ("10_12delinsGG", "15dup", "MAVGGRLSAV"),
    # 8del gets normalized to 9del
    ("8del", "9del", "MAALEAECRL"),
    ("9del", "8del", "MAALEAECRL"),
]


@pytest.mark.parametrize("existing, novel, protein", APPEND_TO_EXISTING)
def test_append_mutation_to_existing_variant(
    existing: str, novel: str, protein: str
) -> None:
    """
    GIVEN a string denoting a HGVS variant
    WHEN we append this variant to an existing Description
    THEN the Description should be updated
    """
    transcript = f"ENST00000375549.8:c.{existing}"
    d = Description(transcript)
    _init_model(d)
    append_mutation(d, novel)

    # Test the first 10 aa of the protein as readout
    assert d.output()["protein"]["predicted"][:10] == protein


APPEND_OVERLAPPING_VARIANT = [
    # Transcript variant, new variant
    # Add another variant at the same location
    ("10C>T", "10del"),
    # Both variants we want to add overlap
    ("=", "[10del;10C>T]"),
    # Both variants already in the description overlap
    ("[10del;10C>T]", "15del"),
    # Overlap with an inversion
    ("10_15inv", "12C>T"),
    ("10_15inv", "12G>T"),
    ("10_15inv", "10C>T"),
]


@pytest.mark.parametrize("existing, novel", APPEND_OVERLAPPING_VARIANT)
def test_appending_overlapping_variants(existing: str, novel: str) -> None:
    """
    GIVEN a transcript with existing variant(s)
    WHEN we attempt to add novel variant which overlaps the existing variants
    THEN we should throw a ValueError
    """
    transcript = f"ENST00000375549.8:c.{existing}"
    d = Description(transcript)
    _init_model(d)

    with pytest.raises(ValueError):
        append_mutation(d, novel)


PARSE_VARIANT = [
    (
        "10del",
        [
            {
                "location": {"type": "point", "position": 10},
                "type": "deletion",
                "source": "reference",
            }
        ],
    ),
    (
        "[10del; 16dup]",
        [
            {
                "location": {"type": "point", "position": 10},
                "type": "deletion",
                "source": "reference",
            },
            {
                "location": {"type": "point", "position": 16},
                "type": "duplication",
                "source": "reference",
            },
        ],
    ),
]


@pytest.mark.parametrize("variant, variant_models", PARSE_VARIANT)
def test_variant_to_model(variant: str, variant_models: List[VariantModel]) -> None:
    """
    GIVEN a string denoting a HGVS variant
    WHEN we parse this into a VariantModel
    THEN it should contain the expected values
    """
    assert variant_to_model(variant) == variant_models
