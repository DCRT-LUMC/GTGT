import pytest
import json
from pydantic.tools import parse_obj_as
import copy

from GTGT import Bed
from GTGT.transcript import Transcript
from GTGT.bed import make_bed
from GTGT.models import TranscriptModel


@pytest.fixture
def Exons() -> Bed:
    exons = [(0, 10), (20, 40), (50, 60), (70, 100)]
    bed = make_bed("chr1", *exons)
    bed.name = "exons"

    return bed


@pytest.fixture
def cds() -> Bed:
    return Bed("chr1", 23, 72, name="cds")


@pytest.fixture
def transcript(Exons: Bed, cds: Bed) -> Transcript:
    """
    Bed records that make up a transcript
    Each positions shown here is 10x
    (i) means inferred by the init method

              0 1 2 3 4 5 6 7 8 9
    exons     -   - -   -   - - -
    cds           - - - - - -
    coding(i)     - -   -   -
    """
    return Transcript(exons=Exons, cds=cds)


def test_transcript_init(transcript: Transcript) -> None:
    assert transcript.exons.name == "exons"
    assert transcript.cds.name == "cds"


def test_coding(transcript: Transcript, Exons: Bed, cds: Bed) -> None:
    # The coding region is the intersection of the exons and the CDS
    coding = Bed(
        "chr1", 23, 72, name="coding", blockSizes=[17, 10, 2], blockStarts=[0, 27, 47]
    )
    # Test that we did not change the cds or exons
    assert transcript.exons == Exons
    assert transcript.cds == cds

    # Test that coding was set
    assert transcript.coding == coding


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


def test_compare_transcripts(transcript: Transcript, cds: Bed) -> None:
    exon_blocks = [
        (0, 10),
        # (20, 40),  # Missing the second exon
        (50, 60),
        (70, 100),
    ]
    exons = make_bed("chr1", *exon_blocks)
    exons.name = "exons"
    smaller = Transcript(exons, cds)

    cmp = smaller.compare(transcript)

    assert cmp["exons"] == pytest.approx(0.71, abs=0.01)
    assert cmp["cds"] == 1
    assert cmp["coding"] == pytest.approx(0.41, abs=0.01)



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
def test_something(variant: str, effect: float, WT):
    # In frame deletion that creates a STOP codon
    modified = copy.deepcopy(WT)
    modified.mutate(variant)

    cmp = modified.compare(WT)
    assert cmp["cds"] == pytest.approx(effect, abs=0.0001)
