import pytest

from gtgt import Bed
from gtgt.mutalyzer import init_description
from gtgt.transcript import Allele, Features
from gtgt.variant import Variant

Ranges = list[tuple[int, int]]


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
def allele(exons: Bed, coding_exons: Bed) -> Allele:
    """
    Bed records that make up an allele
    Each positions shown here is 10x
    (i) means inferred by the init method

                  0 1 2 3 4 5 6 7 8 9
    exons         -   - -   -   - - -
    coding_exons      - -   -   -
    """
    features = Features(rna=[exons], protein=[coding_exons])
    return Allele(features=features)


def test_transcript_init(allele: Allele) -> None:
    features = allele.features
    assert features.rna[0].name == "Exons"
    assert features.protein[0].name == "Coding exons"


def test_empty_allele() -> None:
    """Test creating and working with an empty transcript"""
    a = Allele(Features(rna=[], protein=[]))

    # Test if we can get the records
    assert a.features.records() == []

    # Test if mutating the transcript works
    d = init_description("ENST00000375549.8:c.10del")
    a.features.mutate(d, variants=[])
    assert a


def test_init_allele_features_only() -> None:
    """Test initializing the Allele class"""
    f = Features(rna=[], protein=[])
    a = Allele(features=f)
    assert a.features.rna == []


def test_init_allele_rna_only() -> None:
    rna = "ENST00000375549.8:r.10del"
    a = Allele(rna=rna)

    assert a.features.rna
    assert a.rna


def test_init_allele_dna_only() -> None:
    """Initializing with on DNA is not allowed"""
    dna = "ENST00000375549.8:c.10del"
    with pytest.raises(ValueError):
        Allele(dna=dna)


def test_mutate_allele() -> None:
    """Test mutating an allele"""

    # Note that we run the test here with a c. as 'RNA', which still works but
    # means that the internal coordinate system will be that of the DNA, not
    # the RNA transcript. In the real code, we always work on the rna
    # coordinate system.

    # SNP which introduces a STOP codon
    rna = "ENST00000452863.10:c.100G>T"

    a = Allele(rna=rna)
    a.mutate()

    # Expected exons, 1 nt was changed (position 47846)
    expected_exons = [
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
    ]
    # Expected coding exons
    # (47486,47487) is the first deleted protein sequence
    # (1283, 1286) is the STOP codon
    expected_coding = [(1283, 1286), (47487, 47586)]

    assert a.features.exons is not None
    assert a.features.coding_exons is not None

    exons = a.features.exons.blocks()
    coding_exons = a.features.coding_exons.blocks()

    assert exons == expected_exons
    assert coding_exons == expected_coding


def test_init_allele_features_rna() -> None:
    f = Features(rna=[], protein=[])
    rna = "ENST00000375549.8:r.10del"
    a = Allele(features=f, rna=rna)

    assert a.features.rna == []
    assert a.rna
    assert a.dna is None


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
