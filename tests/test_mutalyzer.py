import pytest
from pathlib import Path
from GTGT.mutalyzer import HGVS_to_genome_range, exonskip, mutation_to_cds_effect, HGVS
from GTGT.transcript import Transcript
from GTGT.models import TranscriptModel
import json
import copy

from itertools import zip_longest
from typing import Any, Tuple


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

    start, end = HGVS_to_genome_range(HGVS(description=description))

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

    start, end = HGVS_to_genome_range(HGVS(description=description))

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
    HGVS_to_genome_range(HGVS(description=description))


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
        HGVS_to_genome_range(HGVS(description=description))


def test_exonskip_SDHD() -> None:
    SDHD = HGVS(description="ENST00000375549.8:c.=")
    results = [
        "ENST00000375549.8:c.53_169del",
        "ENST00000375549.8:c.170_314del",
    ]
    for output, expected in zip_longest(exonskip(SDHD), results):
        assert output == HGVS(description=expected)


def test_exonskip_WT1() -> None:
    WT1 = HGVS(description="ENST00000452863.10:c.=")
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
    for output, expected in zip_longest(exonskip(WT1), results):
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
    WT1 = HGVS(description=description)

    assert mutation_to_cds_effect(WT1) == expected


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
    variant = "ENST00000452863.10:c.87_89del"
    # Frameshift in small in-frame exon 5
    variant = "ENST00000452863.10:c.970del"

    results = WT.analyze(variant)

    print()
    print(json.dumps(results, indent=True))

    assert results["wildtype"]["cds"] == 1.0


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
    modified = copy.deepcopy(WT)
    modified.mutate(variant)

    cmp = modified.compare(WT)
    assert cmp["cds"] == pytest.approx(effect, abs=0.0001)
