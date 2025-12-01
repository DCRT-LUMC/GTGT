import json

import pytest

from gtgt import Bed
from gtgt.models import TranscriptModel
from gtgt.mutalyzer import Therapy, Variant, init_description
from gtgt.transcript import Comparison, Result, Transcript


@pytest.fixture
def Exons() -> Bed:
    exons = [(0, 10), (20, 40), (50, 60), (70, 100)]
    bed = Bed.from_blocks("chr1", exons)
    bed.name = "exons"

    return bed


@pytest.fixture
def coding_exons() -> Bed:
    coding_exons = [(23, 40), (50, 60), (70, 72)]
    bed = Bed.from_blocks("chr1", coding_exons)
    bed.name = "coding_exons"

    return bed


@pytest.fixture
def transcript(Exons: Bed, coding_exons: Bed) -> Transcript:
    """
    Bed records that make up a transcript
    Each positions shown here is 10x
    (i) means inferred by the init method

                  0 1 2 3 4 5 6 7 8 9
    exons         -   - -   -   - - -
    coding_exons      - -   -   -
    """
    return Transcript(exons=Exons, coding_exons=coding_exons)


def test_transcript_init(transcript: Transcript) -> None:
    assert transcript.exons.name == "exons"
    assert transcript.coding_exons.name == "coding_exons"


intersect_selectors = [
    # Selector spans all exons
    (
        Bed("chr1", 0, 100),
        Bed(
            "chr1",
            0,
            100,
            name="exons",
            blockSizes=[10, 20, 10, 30],
            blockStarts=[0, 20, 50, 70],
        ),
    ),
    # Selector on a different chromosome
    (Bed("chr2", 0, 100), Bed("chr1", 0, 0)),
    # Selector intersect the first exon
    (Bed("chr1", 5, 15), Bed("chr1", 5, 10)),
    # Selector intersects the last base of the first exon,
    # and the first base of the second exon
    (Bed("chr1", 9, 21), Bed("chr1", 9, 21, blockSizes=[1, 1], blockStarts=[0, 11])),
]


def test_transcript_init_no_coding(Exons: Bed) -> None:
    t = Transcript(Exons)
    assert not t.coding_exons


@pytest.mark.parametrize("selector, exons", intersect_selectors)
def test_intersect_transcript(
    selector: Bed, exons: Bed, transcript: Transcript
) -> None:
    """Test if intersecting the Transcript updates the exons"""
    transcript.intersect(selector)

    # Ensure the name matches, it's less typing to do that here
    exons.name = "exons"
    assert transcript.exons == exons


overlap_selectors = [
    # Selector spans all exons
    (
        Bed("chr1", 0, 100),
        Bed(
            "chr1",
            0,
            100,
            name="exons",
            blockSizes=[10, 20, 10, 30],
            blockStarts=[0, 20, 50, 70],
        ),
    ),
    # Selector on a different chromosome
    (Bed("chr2", 0, 100), Bed("chr1", 0, 0)),
    # Selector intersect the first exon
    (Bed("chr1", 5, 15), Bed("chr1", 0, 10)),
    # Selector intersects the last base of the first exon,
    # and the first base of the second exon
    (Bed("chr1", 9, 21), Bed("chr1", 0, 40, blockSizes=[10, 20], blockStarts=[0, 20])),
]


@pytest.mark.parametrize("selector, exons", overlap_selectors)
def test_overlap_transcript(selector: Bed, exons: Bed, transcript: Transcript) -> None:
    """Test if overlapping the Transcript updates the exons"""
    transcript.overlap(selector)

    # Ensure the name matches, it's less typing to do that here
    exons.name = "exons"
    assert transcript.exons == exons


subtract_selectors = [
    # Selector spans all exons
    (Bed("chr1", 0, 100), Bed("chr1", 0, 0)),
    # Selector on a different chromosome
    (
        Bed("chr2", 0, 100),
        Bed(
            "chr1",
            0,
            100,
            name="exons",
            blockSizes=[10, 20, 10, 30],
            blockStarts=[0, 20, 50, 70],
        ),
    ),
    # Selector intersect the first exon
    (
        Bed("chr1", 5, 15),
        Bed("chr1", 0, 100, blockSizes=[5, 20, 10, 30], blockStarts=[0, 20, 50, 70]),
    ),
    # Selector intersects the last base of the first exon,
    # and the first base of the second exon
    (
        Bed("chr1", 9, 21),
        Bed("chr1", 0, 100, blockSizes=[9, 19, 10, 30], blockStarts=[0, 21, 50, 70]),
    ),
]


@pytest.mark.parametrize("selector, exons", subtract_selectors)
def test_subtract_transcript(selector: Bed, exons: Bed, transcript: Transcript) -> None:
    """Test if subtracting the Transcript updates the exons"""
    transcript.subtract(selector)

    # Ensure the name matches, it's less typing to do that here
    exons.name = "exons"
    assert transcript.exons == exons


exon_skip_selectors = [
    # A selector on a different chromosome does nothing
    (
        Bed("chr2", 0, 100),
        Bed("chr1", 0, 100, blockSizes=[10, 20, 10, 30], blockStarts=[0, 20, 50, 70]),
    ),
    # Remove the first exon
    (
        Bed("chr1", 0, 1),
        Bed("chr1", 20, 100, blockSizes=[20, 10, 30], blockStarts=[0, 30, 50]),
    ),
    # Selector spans two exons
    (
        Bed("chr1", 9, 21),
        Bed("chr1", 50, 100, blockSizes=[10, 30], blockStarts=[0, 20]),
    ),
    # Remove the last exon
    (
        Bed("chr1", 99, 100),
        Bed("chr1", 0, 60, blockSizes=[10, 20, 10], blockStarts=[0, 20, 50]),
    ),
]


@pytest.mark.parametrize("selector, exons", exon_skip_selectors)
def test_exon_skip_transcript(
    selector: Bed, exons: Bed, transcript: Transcript
) -> None:
    """Test if exon skipping updates the Transcript exons"""
    transcript.exon_skip(selector)

    # Ensure the name matches, it's less typing to do that here
    exons.name = "exons"
    assert transcript.exons == exons


def test_compare_transcripts(transcript: Transcript, coding_exons: Bed) -> None:
    exon_blocks = [
        (0, 10),
        # (20, 40),  # Missing the second exon
        (50, 60),
        (70, 100),
    ]
    exons = Bed.from_blocks("chr1", exon_blocks)
    exons.name = "exons"

    coding_blocks = [
        # (23, 40),  # Missing the second exon
        (50, 60),
        (70, 72),
    ]
    coding_exons = Bed.from_blocks("chr1", coding_blocks)
    coding_exons.name = "coding_exons"

    smaller = Transcript(exons, coding_exons)

    cmp = smaller.compare(transcript)

    assert cmp[0].percentage == pytest.approx(0.71, abs=0.01)
    assert cmp[1].percentage == pytest.approx(0.41, abs=0.01)


@pytest.fixture
def WT() -> Transcript:
    """
    Transcript for WT1, using real genomic positions
    """
    path = "tests/data/ENST00000452863.10.Transcript.json"
    with open(path) as fin:
        js = json.load(fin)

    t = TranscriptModel.model_validate(js)

    return t.to_transcript()


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


MUTATE = [
    (
        "=",  # HGVS mutation (no change)
        ((0, 87), (984, 1101), (1994, 2139), (7932, 8922)),  # Expected exons
        ((35, 87), (984, 1101), (1994, 2139), (7932, 8098)),  # Expected coding_exons
    ),
    (
        # This 1bp deletion introduces a STOP codon
        "40del",  # position (74, 75) was deleted on the RNA
        ((0, 74), (75, 87), (984, 1101), (1994, 2139), (7932, 8922)),
        #           STOP codon is conserved
        ((35, 74), (8095, 8098)),
    ),
    (
        # An in frame deletion that introduces a STOP codon
        "101_106del",  # Position (1032, 1038) was deleted on the RNA
        ((0, 87), (984, 1032), (1038, 1101), (1994, 2139), (7932, 8922)),
        #                      STOP codon
        ((35, 87), (984, 1031), (8095, 8098)),
    ),
    (
        # Dele exon 2 (in frame)
        "53_169del",
        ((0, 87), (1994, 2139), (7932, 8922)),  # Expected exons
        # Exon 2 starts in frame 1, so the deleted positions derived from the
        # changed protein positions are slightly different: (986, 1995)
        # The nucleotieds (984, 986) are used with (1996, 1997) to form the
        # first non-deleted amino acid
        ((35, 87), (1996, 2139), (7932, 8098)),  # Expected coding_exons
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
    offset = 112086872

    # Exons and coding exons of SDHD
    exons = Bed.from_blocks(
        chrom,
        [
            (112086872, 112086959),
            (112087856, 112087973),
            (112088866, 112089011),
            (112094804, 112095794),
        ],
    )
    coding_exons = Bed.from_blocks(
        chrom,
        [
            (112086907, 112086959),
            (112087856, 112087973),
            (112088866, 112089011),
            (112094804, 112094970),
        ],
    )

    # Variant to test
    d = init_description(f"{transcript}:c.{variant}")
    v = [
        Variant.from_model(delins)
        for delins in d.de_hgvs_internal_indexing_model["variants"]
    ]

    # Add the offset to the expected exon blocks
    exon_blocks = [(start + offset, end + offset) for start, end in exon_blocks]
    coding_exon_blocks = [
        (start + offset, end + offset) for start, end in coding_exon_blocks
    ]

    SDHD = Transcript(exons=exons, coding_exons=coding_exons)
    SDHD.mutate(d, v)

    assert SDHD.exons.blocks() == exon_blocks
    assert SDHD.coding_exons.blocks() == coding_exon_blocks


MUTATE = [
    (
        # HGVS mutation (no change)
        "=",
        # Expected exon
        (
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
        ),
        # Expected coding exons
        (
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
        ),
    ),
    (
        # SNP which introduces a STOP codon
        "100G>T",
        # Expected exons, 1 nt was changed (position 47846)
        (
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
        ),
        # Expected coding exons
        # (47486,47487) is the first deleted protein sequence
        # (1283, 1286) is the STOP codon
        ((1283, 1286), (47487, 47586)),
    ),
    (
        # In fram deletion which introduces a STOP codon
        "134_139del",
        # Expected exon, positions (47447, 47452) were deleted
        (
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
        ),
        # Expected coding exons,
        # 47454 is the last conserved nucleotide
        # (1283, 1286) is the STOP codon
        ((1283, 1286), (47454, 47586)),
    ),
    (
        # Delete exon 2 (in frame)
        "662_784del",
        # Expected exon (40282, 40845) deleted
        (
            (0, 1405),
            (4197, 4290),
            (4891, 4981),
            (8482, 8633),
            (12173, 12270),
            (28715, 28766),
            (29802, 29880),
            (40181, 40284),
            (46925, 47765),
        ),
        # Expected coding exons, (40772, 40840) deleted on the protein level
        # Note that exon 2 starts in frame 1, so when going from the protein
        # sequence, a small region in exon 3 has also changed
        (
            (1283, 1405),
            (4197, 4290),
            (4891, 4981),
            (8482, 8633),
            (12173, 12270),
            (28715, 28766),
            (29802, 29880),
            (40181, 40282),
            (46925, 47586),
        ),
    ),
]


@pytest.mark.parametrize("variant, exon_blocks, coding_exon_blocks", MUTATE)
def test_mutate_reverse(
    variant: str, exon_blocks: Ranges, coding_exon_blocks: Ranges
) -> None:
    # Features of SDHD
    transcript = "ENST00000452863.10"
    chrom = "chr11"
    offset = 32387774

    # Exons and coding exons of WT1
    exons = Bed.from_blocks(
        chrom,
        [
            (32387774, 32389179),
            (32391971, 32392064),
            (32392665, 32392755),
            (32396256, 32396407),
            (32399947, 32400044),
            (32416489, 32416540),
            (32417576, 32417654),
            (32427955, 32428058),
            (32428496, 32428619),
            (32434699, 32435539),
        ],
    )
    coding_exons = Bed.from_blocks(
        chrom,
        [
            (32389057, 32389179),
            (32391971, 32392064),
            (32392665, 32392755),
            (32396256, 32396407),
            (32399947, 32400044),
            (32416489, 32416540),
            (32417576, 32417654),
            (32427955, 32428058),
            (32428496, 32428619),
            (32434699, 32435360),
        ],
    )

    # Variant to test
    d = init_description(f"{transcript}:c.{variant}")
    v = [
        Variant.from_model(delins)
        for delins in d.de_hgvs_internal_indexing_model["variants"]
    ]

    # Add the offset to the expected exon blocks
    exon_blocks = [(start + offset, end + offset) for start, end in exon_blocks]
    coding_exon_blocks = [
        (start + offset, end + offset) for start, end in coding_exon_blocks
    ]

    WT1 = Transcript(exons=exons, coding_exons=coding_exons)
    WT1.mutate(d, v)

    assert WT1.exons.blocks() == exon_blocks
    assert WT1.coding_exons.blocks() == coding_exon_blocks


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
