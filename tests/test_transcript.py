from typing import Sequence

import pytest

from gtgt import Bed
from gtgt.mutalyzer import init_description
from gtgt.therapy import TherapyReport
from gtgt.transcript import (
    Comparison,
    Features,
    Result,
    Transcript,
    is_of_interest,
)
from gtgt.variant import Variant


def test_Result_init() -> None:
    t = TherapyReport("skip exon 5", "ENST123:c.49_73del", "Try to skip exon 5", list())
    c = Comparison("Coding exons", 0.5, "100/200")

    Result(therapy=t, comparison=[c])

    assert True


def test_Result_comparison() -> None:
    t1 = TherapyReport(
        "skip exon 5", "ENST123:c.49_73del", "Try to skip exon 5", list()
    )
    c1 = Comparison("Coding exons", 0.5, "100/200")
    r1 = Result(therapy=t1, comparison=[c1])

    t2 = TherapyReport(
        "skip exon 6", "ENST123:c.49_73del", "Try to skip exon 5", list()
    )
    c2 = Comparison("Coding exons", 0.2, "100/200")
    r2 = Result(therapy=t2, comparison=[c2])

    # Results in the wrong order
    results = [r2, r1]

    # Highest scoring Results should come first
    assert sorted(results, reverse=True) == [r1, r2]


@pytest.mark.parametrize(
    "patient, patient_vars, therapy, therapy_vars, expected",
    [
        # Identical patient and therapy
        (
            # Patient
            Features(rna=[], protein=[]),
            [],
            # Therapy
            Features(rna=[], protein=[]),
            [],
            False,
        ),
        # The therapy got rid of one of the variants
        (
            # Patient
            Features(rna=[], protein=[]),
            [Variant(10, 20)],
            # Therapy
            Features(rna=[], protein=[]),
            [],
            True,
        ),
        # One of the therapy features is smaller
        (
            # Patient
            Features(rna=[Bed("", chromStart=10, chromEnd=20)], protein=[]),
            [],
            # Therapy
            Features(rna=[Bed("", chromStart=11, chromEnd=20)], protein=[]),
            [],
            False,
        ),
        # One of the therapy features is larger
        (
            # Patient
            Features(rna=[Bed("", chromStart=10, chromEnd=20)], protein=[]),
            [],
            # Therapy
            Features(rna=[Bed("", chromStart=10, chromEnd=21)], protein=[]),
            [],
            True,
        ),
        # One of the therapy features is completely missing
        (
            # Patient
            Features(rna=[Bed("", chromStart=10, chromEnd=20)], protein=[]),
            [],
            # Therapy
            Features(rna=[Bed("", chromStart=10, chromEnd=10)], protein=[]),
            [],
            False,
        ),
        # One of the patient features is completely missing
        (
            # Patient
            Features(rna=[Bed("", chromStart=10, chromEnd=10)], protein=[]),
            [],
            # Therapy
            Features(rna=[Bed("", chromStart=10, chromEnd=20)], protein=[]),
            [],
            True,
        ),
        # The therapy is smaller, but has a region the patient lacks
        (
            # Patient
            Features(rna=[Bed("", chromStart=10, chromEnd=20)], protein=[]),
            [],
            # Therapy
            Features(rna=[Bed("", chromStart=20, chromEnd=21)], protein=[]),
            [],
            True,
        ),
    ],
)
def test_therapy_is_of_interest(
    patient: Features,
    patient_vars: Sequence[Variant],
    therapy: Features,
    therapy_vars: Sequence[Variant],
    expected: bool,
) -> None:
    assert is_of_interest(patient, patient_vars, therapy, therapy_vars) == expected


MUTATE = [
    (
        "=",  # HGVS mutation (no change)
        [(0, 87), (984, 1101), (1994, 2139), (7932, 8922)],  # Expected exons
        [(35, 87), (984, 1101), (1994, 2139), (7932, 8098)],  # Expected coding_exons
    ),
    (
        # This 1bp deletion introduces a STOP codon
        "40del",  # position (74, 75) was deleted on the RNA
        [(0, 74), (75, 87), (984, 1101), (1994, 2139), (7932, 8922)],
        #           STOP codon is conserved
        [(35, 74), (8095, 8098)],
    ),
    (
        # An in frame deletion that introduces a STOP codon
        "101_106del",  # Position (1032, 1038) was deleted on the RNA
        [(0, 87), (984, 1032), (1038, 1101), (1994, 2139), (7932, 8922)],
        #                      STOP codon
        [(35, 87), (984, 1031), (8095, 8098)],
    ),
    (
        # Dele exon 2 (in frame)
        "53_169del",
        [(0, 87), (1994, 2139), (7932, 8922)],  # Expected exons
        # Exon 2 starts in frame 1, so the deleted positions derived from the
        # changed protein positions are slightly different: (986, 1995)
        # The nucleotides (984, 986) are used with (1996, 1997) to form the
        # first non-deleted amino acid
        [(35, 87), (1996, 2139), (7932, 8098)],  # Expected coding_exons
    ),
]

Ranges = list[tuple[int, int]]


@pytest.mark.parametrize("variant, exon_blocks, coding_exon_blocks", MUTATE)
def test_mutate_forward(
    variant: str, exon_blocks: Ranges, coding_exon_blocks: Ranges
) -> None:
    # Features of SDHD
    transcript = "ENST00000375549.8"
    chrom = "chr11"

    # Exons and coding exons of SDHD
    exons = Bed.from_blocks(
        chrom,
        [(0, 87), (984, 1101), (1994, 2139), (7932, 8922)],
    )
    exons.name = "Exons"

    coding_exons = Bed.from_blocks(
        chrom,
        [(35, 87), (984, 1101), (1994, 2139), (7932, 8098)],
    )
    coding_exons.name = "Coding exons"

    # Variant to test
    d = init_description(f"{transcript}:c.{variant}")
    v = [
        Variant.from_model(delins)
        for delins in d.de_hgvs_internal_indexing_model["variants"]
    ]

    SDHD = Features(rna=[exons], protein=[coding_exons])
    SDHD.mutate(d, v)

    assert SDHD.exons and SDHD.exons.blocks() == exon_blocks
    assert SDHD.coding_exons and SDHD.coding_exons.blocks() == coding_exon_blocks


def test_Comparison_from_dict() -> None:
    """Test creating a Comparison from a dict"""
    c = Comparison("Wildtype", percentage=100, basepairs="100/100")

    d = {"name": "Wildtype", "percentage": 100, "basepairs": "100/100"}

    assert Comparison.from_dict(d) == c


def test_Result_from_dict() -> None:
    """Test  creating a Result from a dict"""

    r = Result(
        therapy=TherapyReport(
            name="wildtype",
            hgvsc="ENST:c.=",
            description="Free text",
            variants=[Variant(10, 12, inserted="ATG")],
        ),
        comparison=[Comparison("Exons", percentage=100, basepairs="120/120")],
    )

    therapy = {
        "name": "wildtype",
        "hgvsc": "ENST:c.=",
        "description": "Free text",
        "variants": [{"start": 10, "end": 12, "inserted": "ATG", "deleted": ""}],
    }
    comparison = [{"name": "Exons", "percentage": 100, "basepairs": "120/120"}]

    d = {"therapy": therapy, "comparison": comparison}

    assert Result.from_dict(d) == r


@pytest.mark.parametrize(
    "hgvs, wildtype",
    [
        # Test a c. variant
        ("NM_003002.4:c.102del", "NM_003002.4:c.="),
        # A hgvs with a selector
        ("NC_000011.10(NM_003002.4):c.102del", "NC_000011.10(NM_003002.4):c.="),
        # A variant on the RNA
        ("NM_003002.4:r.102del", "NM_003002.4:r.="),
    ],
)
def test_transcript_from_hgvs(hgvs: str, wildtype: str) -> None:
    """Test creating a transcript from HGVS"""

    t = Transcript(hgvs)

    assert t.d.input_description == wildtype


protein_changing_rna_variants = [
    # No change
    "NM_003002.4:r.=",
    # Frameshift mutation
    "NM_003002.4:r.10del",
    # Insertion
    "NM_003002.4:r.10_11insa",
    # Out of frame deletion (long tail)
    "NM_003002.4:r.50_53del",
    # Go out and back into frame
    "NM_003002.4:r.[50del;70_71insa]",
    # SNP
    "NC_000011.10(NM_003002.4):r.274g>u",
    # Another frameshift
    "NC_000011.10(NM_003002.4):r.52del",
    # Reverse strand
    "ENST00000452863.10:r.=",
    "ENST00000452863.10:r.10del",
]
protein_changing_dna_variants = [
    # No change
    "NM_003002.4:c.=",
    # Frameshift mutation
    "NM_003002.4:c.10del",
    # Insertion
    "NM_003002.4:c.10_11insA",
    # Out of frame deletion (long tail)
    "NM_003002.4:c.50_53del",
    # Go out and back into frame
    "NM_003002.4:c.[50del;70_71insA]",
    # SNP
    "NC_000011.10(NM_003002.4):c.274G>T",
    # Another frameshift
    "NC_000011.10(NM_003002.4):c.52del",
    # Reverse strand
    "ENST00000452863.10:c.=",
    "ENST00000452863.10:c.10del",
]


@pytest.mark.parametrize(
    "hgvs", protein_changing_rna_variants + protein_changing_dna_variants
)
def test_hgvs_description_from_Transcript(hgvs: str) -> None:
    """
    Test generating the corresponding HGVS description from a Transcript
    """
    t = Transcript(hgvs)

    generated_hgvs = t.hgvs()
    assert generated_hgvs == hgvs
