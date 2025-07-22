import json

import pytest

from gtgt import Bed
from gtgt.models import TranscriptModel
from gtgt.mutalyzer import Therapy
from gtgt.transcript import Comparison, Result, Transcript


@pytest.fixture
def Exons() -> Bed:
    exons = [(0, 10), (20, 40), (50, 60), (70, 100)]
    bed = Bed.from_blocks("chr1", *exons)
    bed.name = "exons"

    return bed


@pytest.fixture
def cds() -> Bed:
    return Bed("chr1", 23, 72, name="cds")


@pytest.fixture
def coding_exons() -> Bed:
    coding_exons = [(23, 40), (50, 60), (70, 72)]
    bed = Bed.from_blocks("chr1", *coding_exons)
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
    exons = Bed.from_blocks("chr1", *exon_blocks)
    exons.name = "exons"

    coding_blocks = [
        # (23, 40),  # Missing the second exon
        (50, 60),
        (70, 72),
    ]
    coding_exons = Bed.from_blocks("chr1", *coding_blocks)
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
