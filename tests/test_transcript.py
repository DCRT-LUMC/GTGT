from typing import Sequence

import pytest

from gtgt import Bed
from gtgt.mutalyzer import init_description
from gtgt.therapy import Therapy
from gtgt.transcript import Allele, Comparison, Features, Result, Transcript, is_of_interest
from gtgt.variant import Variant


@pytest.fixture
def exons() -> Bed:
    exons = [(0, 10), (20, 40), (50, 60), (70, 100)]
    bed = Bed.from_blocks("chr1", exons)
    bed.name = "Exons"

    return bed


@pytest.fixture
def coding_exons() -> Bed:
    coding_exons = [(23, 40), (50, 60), (70, 72)]
    bed = Bed.from_blocks("chr1", coding_exons)
    bed.name = "Coding exons"

    return bed


@pytest.fixture
def transcript(exons: Bed, coding_exons: Bed) -> Allele:
    """
    Bed records that make up a transcript
    Each positions shown here is 10x
    (i) means inferred by the init method

                  0 1 2 3 4 5 6 7 8 9
    exons         -   - -   -   - - -
    coding_exons      - -   -   -
    """
    return Allele(Features(rna=[exons], protein=[coding_exons]))


def test_transcript_init(transcript: Allele) -> None:
    features = transcript.features
    assert features.rna[0].name == "Exons"
    assert features.protein[0].name == "Coding exons"


def test_empty_transcript() -> None:
    """Test creating and working with an empty transcript"""
    t = Allele(Features(rna=[], protein=[]))

    # Test if we can get the records
    assert t.features.records() == []

    # Test if mutating the transcript works
    d = init_description("ENST00000375549.8:c.10del")
    t.features.mutate(d, variants=[])
    assert t


def test_Result_init() -> None:
    t = Therapy("skip exon 5", "ENST123:c.49_73del", "Try to skip exon 5", list())
    c = Comparison("Coding exons", 0.5, "100/200")

    r = Result(therapy=t, comparison=[c])

    assert True


def test_Result_comparison() -> None:
    t1 = Therapy("skip exon 5", "ENST123:c.49_73del", "Try to skip exon 5", list())
    c1 = Comparison("Coding exons", 0.5, "100/200")
    r1 = Result(therapy=t1, comparison=[c1])

    t2 = Therapy("skip exon 6", "ENST123:c.49_73del", "Try to skip exon 5", list())
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


MUTATE = [
    (
        # HGVS mutation (no change)
        "=",
        # Expected exon
        [
            (0, 1405),
            (4197, 4290),
            (4891, 4981),
            (8482, 8633),
            (12173, 12270),
            (28715, 28766),
            (29802, 29880),
            (40181, 40284),
            (40722, 40845),
            (46925, 47765),
        ],
        # Expected coding exons
        [
            (1283, 1405),
            (4197, 4290),
            (4891, 4981),
            (8482, 8633),
            (12173, 12270),
            (28715, 28766),
            (29802, 29880),
            (40181, 40284),
            (40722, 40845),
            (46925, 47586),
        ],
    ),
    (
        # SNP which introduces a STOP codon
        "100G>T",
        # Expected exons, 1 nt was changed (position 47846)
        [
            (0, 1405),
            (4197, 4290),
            (4891, 4981),
            (8482, 8633),
            (12173, 12270),
            (28715, 28766),
            (29802, 29880),
            (40181, 40284),
            (40722, 40845),
            (46925, 47486),
            (47487, 47765),
        ],
        # Expected coding exons
        # (47486,47487) is the first deleted protein sequence
        # (1283, 1286) is the STOP codon
        [(1283, 1286), (47487, 47586)],
    ),
    (
        # In fram deletion which introduces a STOP codon
        "134_139del",
        # Expected exon, positions (47447, 47452) were deleted
        [
            (0, 1405),
            (4197, 4290),
            (4891, 4981),
            (8482, 8633),
            (12173, 12270),
            (28715, 28766),
            (29802, 29880),
            (40181, 40284),
            (40722, 40845),
            (46925, 47447),
            (47453, 47765),
        ],
        # Expected coding exons,
        # 47454 is the last conserved nucleotide
        # (1283, 1286) is the STOP codon
        [(1283, 1286), (47454, 47586)],
    ),
    (
        # Delete exon 2 (in frame)
        "662_784del",
        # Expected exon (40282, 40845) deleted
        [
            (0, 1405),
            (4197, 4290),
            (4891, 4981),
            (8482, 8633),
            (12173, 12270),
            (28715, 28766),
            (29802, 29880),
            (40181, 40284),
            (46925, 47765),
        ],
        # Expected coding exons, (40772, 40840) deleted on the protein level
        # Note that exon 2 starts in frame 1, so when going from the protein
        # sequence, a small region in exon 3 has also changed
        [
            (1283, 1405),
            (4197, 4290),
            (4891, 4981),
            (8482, 8633),
            (12173, 12270),
            (28715, 28766),
            (29802, 29880),
            (40181, 40282),
            (46925, 47586),
        ],
    ),
]


@pytest.mark.parametrize("variant, exon_blocks, coding_exon_blocks", MUTATE)
def test_mutate_reverse(
    variant: str, exon_blocks: Ranges, coding_exon_blocks: Ranges
) -> None:
    # Features of SDHD
    transcript = "ENST00000452863.10"
    chrom = "chr11"

    # Exons and coding exons of WT1
    exons = Bed.from_blocks(
        chrom,
        [
            (0, 1405),
            (4197, 4290),
            (4891, 4981),
            (8482, 8633),
            (12173, 12270),
            (28715, 28766),
            (29802, 29880),
            (40181, 40284),
            (40722, 40845),
            (46925, 47765),
        ],
    )
    exons.name = "Exons"

    coding_exons = Bed.from_blocks(
        chrom,
        [
            (1283, 1405),
            (4197, 4290),
            (4891, 4981),
            (8482, 8633),
            (12173, 12270),
            (28715, 28766),
            (29802, 29880),
            (40181, 40284),
            (40722, 40845),
            (46925, 47586),
        ],
    )
    coding_exons.name = "Coding exons"

    # Variant to test
    d = init_description(f"{transcript}:c.{variant}")
    v = [Variant.from_model(delins) for delins in d.delins_model["variants"]]

    WT1 = Allele(Features(rna=[exons], protein=[coding_exons]))
    WT1.features.mutate(d, v)

    assert WT1.features.exons and WT1.features.exons.blocks() == exon_blocks
    assert (
        WT1.features.coding_exons
        and WT1.features.coding_exons.blocks() == coding_exon_blocks
    )


def test_Comparison_from_dict() -> None:
    """Test creating a Comparison from a dict"""
    c = Comparison("Wildtype", percentage=100, basepairs="100/100")

    d = {"name": "Wildtype", "percentage": 100, "basepairs": "100/100"}

    assert Comparison.from_dict(d) == c


def test_Result_from_dict() -> None:
    """Test  creating a Result from a dict"""

    r = Result(
        therapy=Therapy(
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


@pytest.mark.parametrize("hgvs, wildtype",
     [
        # Test a c. variant
        ("NM_003002.4:c.102del", "NM_003002.4:c.="),
        # A hgvs with a selector
        ("NC_000011.10(NM_003002.4):c.102del", "NC_000011.10(NM_003002.4):c.="),
        # A variant on the RNA
        ("NM_003002.4:r.102del", "NM_003002.4:r.="),

     ])
def test_transcript_from_hgvs(hgvs: str, wildtype: str) -> None:
    """Test creating a transcript from HGVS"""

    t = Transcript(hgvs)

    assert t.d.input_description == wildtype
