import pytest
from mutalyzer.description import Description, to_rna_reference_model, model_to_string
from mutalyzer.converter.to_hgvs_coordinates import to_hgvs_locations
from pathlib import Path
from GTGT.mutalyzer import HGVS_to_genome_range, exonskip
from GTGT.models import HGVS

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
