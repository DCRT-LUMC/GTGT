import json
from collections.abc import Sequence
from itertools import zip_longest

import pytest
from mutalyzer.description import Description

from gtgt.models import TranscriptModel
from gtgt.mutalyzer import (
    Therapy,
    Variant,
    _exon_string,
    changed_protein_positions,
    get_chrom_name,
    get_exons,
    init_description,
    protein_prediction,
    skip_adjacent_exons,
    sliding_window,
)
from gtgt.transcript import Transcript


def SDHD_description() -> Description:
    """SDHD, on the forward strand"""
    return init_description("ENST00000375549.8:c.=")


def WT1_description() -> Description:
    """WT1, on the reverse strand"""
    return init_description("ENST00000452863.10:c.=")


def test_one_adjacent_exonskip_forward() -> None:
    d = SDHD_description()
    results = [
        "ENST00000375549.8:c.53_169del",
        "ENST00000375549.8:c.170_314del",
    ]
    for output, expected in zip_longest(skip_adjacent_exons(d), results):
        assert output.hgvsc == expected


def test_two_adjacent_exonskip_SDHD() -> None:
    d = SDHD_description()
    results = [
        "ENST00000375549.8:c.53_314del",
    ]
    for output, expected in zip_longest(
        skip_adjacent_exons(d, number_to_skip=2), results
    ):
        assert output.hgvsc == expected


def test_no_possible_exonskip_SDHD() -> None:
    """
    GIVEN a transcript with 4 exons (2 can be skipped)
    WHEN we try to skip 3 adjacent exons
    THEN we should get an empty list of therapies
    """
    d = SDHD_description()
    assert skip_adjacent_exons(d, number_to_skip=3) == list()


def test_one_adjacent_exonskip_WT1() -> None:
    d = WT1_description()
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
    # for output, expected in zip_longest(skip_adjacent_exons(d), results):
    for output, expected in zip_longest(skip_adjacent_exons(d), results):
        assert output.hgvsc == expected


def test_two_adjacent_exonskip_WT1() -> None:
    d = WT1_description()
    results = [
        "ENST00000452863.10:c.662_887del",
        "ENST00000452863.10:c.785_965del",
        "ENST00000452863.10:c.888_1016del",
        "ENST00000452863.10:c.966_1113del",
        "ENST00000452863.10:c.1017_1264del",
        "ENST00000452863.10:c.1114_1354del",
        "ENST00000452863.10:c.1265_1447del",
    ]
    for output, expected in zip_longest(skip_adjacent_exons(d, 2), results):
        assert output.hgvsc == expected


def test_sliding_window_size_one() -> None:
    s = "ABCDEF"
    # assert list(sliding_window(s, 1)) == [[x] for x in "A B C D E F".split()]
    assert list(sliding_window(s, 1)) == [["A"], ["B"], ["C"], ["D"], ["E"], ["F"]]


def test_sliding_window_size_2() -> None:
    s = "ABCDEF"
    assert list(sliding_window(s, 2)) == [
        ["A", "B"],
        ["B", "C"],
        ["C", "D"],
        ["D", "E"],
        ["E", "F"],
    ]


@pytest.fixture(scope="session")
def WT() -> Transcript:
    """
    Transcript for WT1, using real genomic positions
    """
    path = "tests/data/ENST00000452863.10.Transcript.json"
    with open(path) as fin:
        js = json.load(fin)

    t = TranscriptModel.model_validate(js)

    return t.to_transcript()


@pytest.mark.parametrize(
    "variant",
    "13T>A 970del 970_971insA 997_999delinsTAA 1000dup 10_11inv 994_996A[9]".split(),
)
def test_analyze_supported_variant_types(WT: Transcript, variant: str) -> None:
    hgvs = f"ENST00000452863.10:c.{variant}"
    WT.analyze(hgvs)


def test_analyze_transcript(WT: Transcript) -> None:
    # In frame deletion that creates a STOP codon
    # variant = "ENST00000452863.10:c.87_89del"
    # Frameshift in small in-frame exon 5
    variant = "ENST00000452863.10:c.970del"

    results = WT.analyze(variant)

    # Test the content of the 'wildtype' result
    wildtype = results[0]
    assert wildtype.therapy.name == "Wildtype"
    coding_exons = wildtype.comparison[1]
    assert coding_exons.name == "coding_exons"
    assert coding_exons.percentage == 1.0

    input = results[1]
    assert input.therapy.name == "Input"
    assert input.therapy.hgvsc == variant
    coding_exons = input.comparison[1]
    # basepairs are not a float, so easier to match than .percentage
    assert coding_exons.basepairs == "18845/46303"


@pytest.mark.xfail
def test_analyze_transcript_r_coordinate(WT: Transcript) -> None:
    """Test analyzing a transcript using the r. coordinate system

    Note, this test should pass, but r. variant are currently not supported
    """
    # In frame deletion that creates a STOP codon
    variant = "ENST00000452863.10:r.970del"

    results = WT.analyze(variant)

    # Test the content of the 'wildtype' result
    wildtype = results[0]
    assert wildtype.therapy.name == "Wildtype"
    coding_exons = wildtype.comparison[1]
    assert coding_exons.name == "coding_exons"
    assert coding_exons.percentage == 1.0

    input = results[1]
    assert input.therapy.name == "Input"
    assert input.therapy.hgvsc == variant
    coding_exons = input.comparison[1]
    # basepairs are not a float, so easier to match than .percentage
    assert coding_exons.basepairs == "18845/46303"


PROTEIN_EXTRACTOR = [
    # No protein description
    ("", "", []),
    # No change
    ("A", "A", []),
    # A single change
    ("A", "T", [(0, 1)]),
    # A single change on the second position
    ("AAA", "ATA", [(1, 2)]),
    # Change in a repeat region
    ("AA", "A", [(1, 2)]),
    ("AAA", "A", [(1, 3)]),
    # A delins
    ("AAA", "ATC", [(1, 3)]),
    # An insertion, which we ignore
    ("AAA", "AATA", []),
    # A delins of TAG, which is equivalent to two insertions
    ("AAA", "ATAGAA", []),
    # A delins which is equivalent to a deletion
    ("AAA", "AGGGA", [(1, 2)]),
    # Multiple deletions
    ("AAA", "TAT", [(0, 1), (2, 3)]),
]


@pytest.mark.parametrize("reference, observed, expected", PROTEIN_EXTRACTOR)
def test_changed_protein_positions(
    reference: str, observed: str, expected: list[tuple[int, int]]
) -> None:
    """
    GIVEN a referene and observed sequence
    WHEN we extrat the protein changes
    THEN we should get 0 based positions of the deleted residues
    """
    assert changed_protein_positions(reference, observed) == expected


def test_get_exons_forward() -> None:
    """Text extracting exons from a Description object"""
    d = SDHD_description()
    expected = (0, 87)

    assert get_exons(d, in_transcript_order=True)[0] == expected
    assert get_exons(d, in_transcript_order=False)[0] == expected


def test_exons_forward() -> None:
    """Text extracting exons from a Description object"""
    d = WT1_description()
    expected_transcript_order = (46925, 47765)
    expected_genomic_order = (0, 1405)

    assert get_exons(d, in_transcript_order=True)[0] == expected_transcript_order
    assert get_exons(d, in_transcript_order=False)[0] == expected_genomic_order


EXON_DESCRIPTION = [
    ([2], "exon 2"),
    ([3, 5], "exons 3 and 5"),
    ([3, 4, 5, 6], "exons 3, 4, 5 and 6"),
]


@pytest.mark.parametrize("exons, expected", EXON_DESCRIPTION)
def test_exon_string(exons: Sequence[int], expected: str) -> None:
    assert _exon_string(exons) == expected


def test_Therapy_from_dict() -> None:
    """Test creating a Therapy from a dict"""
    variants = [Variant(10, 11, inserted="A", deleted="T")]
    therapy = Therapy(
        "Wildtype", hgvsc="ENST:c.=", description="free text", variants=variants
    )

    d = {
        "name": "Wildtype",
        "hgvsc": "ENST:c.=",
        "description": "free text",
        "variants": [{"start": 10, "end": 11, "inserted": "A", "deleted": "T"}],
    }

    assert Therapy.from_dict(d) == therapy


def test_protein_prediction_unknown() -> None:
    """Test overwriting unknown protein prediction

    Sometimes, mutalyzer will only generate :p.? as a protein prediction
    If that happens, we use in_frame_description overwrite the protein
    description
    """
    d = SDHD_description()
    # Variants that give rise to a :p.? prediction from mutalyzer
    variants = [Variant(start=1031, end=1032), Variant(start=1994, end=2139)]

    id = "ENST00000375549.8(ENSP00000364699)"
    p_variant = "Leu35_Leu159delinsPheArgThrAspLeuSerGlnAsnGlyValGluCysSerThrTyrThrCysHisArgAlaThrIleGlyProTrpThrSerCysTyr"
    assert protein_prediction(d, variants)[0] == f"{id}:p.{p_variant}"


def test_transcript_from_description_WT1() -> None:
    """Test creating a Transcript from a Mutalyzer Description"""
    d = WT1_description()
    t = Transcript.from_description(d)
    # Manually verified block starts for the exons of WT1
    assert t.exons.blockStarts == [
        0,
        4197,
        4891,
        8482,
        12173,
        28715,
        29802,
        40181,
        40722,
        46925,
    ]
    # Manually verified block starts for the coding exons of WT1
    assert t.coding_exons.blockStarts == [
        0,
        2914,
        3608,
        7199,
        10890,
        27432,
        28519,
        38898,
        39439,
        45642,
    ]


def test_transcript_from_description_SDHD() -> None:
    """Test creating a Transcript from a Mutalyzer Description"""
    d = SDHD_description()
    t = Transcript.from_description(d)
    # Manually verified block starts for the exons of SDHD
    assert t.exons.blockStarts == [0, 984, 1994, 7932]
    # Manually verified block starts for the coding exons of SDHD
    assert t.coding_exons.blockStarts == [0, 949, 1959, 7897]


def test_transcript_from_NC_NM_forward() -> None:
    """Test creating a forward Transcript from a Mutalyzer NC(NM) description"""
    d = init_description("NC_000011.10(NM_003002.4):c.=")
    t = Transcript.from_description(d)

    assert t.exons.blocks() == [
        (5026, 5113),
        (6010, 6127),
        (7020, 7165),
        (12958, 13948),
    ]
    assert t.coding_exons.blocks() == [
        (5061, 5113),
        (6010, 6127),
        (7020, 7165),
        (12958, 13124),
    ]


def test_analyze_NC_NM_forward() -> None:
    """Test mutating an NC(NM) transcript on the forward strand

    NM_003002 is SDHD
    """
    hgvs = "NC_000011.10(NM_003002.4):c.102del"
    d = init_description(hgvs)
    t = Transcript.from_description(d)

    results = t.analyze(hgvs)

    # Look at the results for skipping exon 2
    skip2 = results[2]
    assert skip2.therapy.name == "Skip exon 2"
    assert skip2.therapy.hgvsc == "NM_003002.4:c.53_169del"
    assert skip2.therapy.hgvsp == "NC_000011.10(NP_002993.1):p.(Leu19_Ser57del)"

    # Check the remaining basepairs for the exons
    exons = skip2.comparison[0]
    assert exons.basepairs == "1222/1339"

    # Check the remaining basepairs for the coding exons
    coding_exons = skip2.comparison[1]
    assert coding_exons.basepairs == "361/480"


def test_analyze_NM_forward() -> None:
    """Test mutating an NM transcript on the forward strand

    NM_003002 is SDHD
    """
    hgvs = "NM_003002.4:c.102del"
    d = init_description(hgvs)
    t = Transcript.from_description(d)

    results = t.analyze(hgvs)

    # Look at the results for skipping exon 2
    skip2 = results[2]
    assert skip2.therapy.name == "Skip exon 2"
    assert skip2.therapy.hgvsc == "NM_003002.4:c.53_169del"
    assert skip2.therapy.hgvsp == "NM_003002.4(NP_002993.1):p.(Leu19_Ser57del)"

    # Check the remaining basepairs for the exons
    exons = skip2.comparison[0]
    assert exons.basepairs == "1222/1339"

    # Check the remaining basepairs for the coding exons
    coding_exons = skip2.comparison[1]
    assert coding_exons.basepairs == "361/480"


def test_transcript_from_NC_NM_reverse() -> None:
    """Test creating a reverse Transcript from a Mutalyzer NC(NM) description

    NM_012459.4 is TIMM8B
    """
    d = init_description("NC_000011.10(NM_012459.4):c.=")
    t = Transcript.from_description(d)

    assert t.exons.blocks() == [(2953, 3616), (4793, 4910)]
    assert t.coding_exons.blocks() == [(3448, 3616), (4793, 4877)]


def test_analyze_NC_NM_reverse() -> None:
    """Test creating a reverse Transcript from a Mutalyzer NC(NM) description

    NM_012459.4 is TIMM8B
    """
    hgvs = "NC_000011.10(NM_012459.4):c.100del"
    d = init_description(hgvs)
    t = Transcript.from_description(d)

    results = t.analyze(hgvs)

    # TIMM8B only has 2 exons, so there will be no exon skips proposed
    assert len(results) == 2

    # Look at the results for the input variant
    input_results = results[1]
    assert input_results.therapy.name == "Input"
    assert (
        input_results.therapy.hgvsp == "NC_000011.10(NP_036591.3):p.(Glu34SerfsTer15)"
    )


def test_analyze_NM_reverse() -> None:
    """Test creating a reverse Transcript from a Mutalyzer NC(NM) description

    NM_012459.4 is TIMM8B
    """
    hgvs = "NM_012459.4:c.100del"
    d = init_description(hgvs)
    t = Transcript.from_description(d)

    results = t.analyze(hgvs)

    # TIMM8B only has 2 exons, so there will be no exon skips proposed
    assert len(results) == 2

    # Look at the results for the input variant
    input_results = results[1]
    assert input_results.therapy.name == "Input"
    assert input_results.therapy.hgvsp == "NM_012459.4(NP_036591.3):p.(Glu34SerfsTer15)"


@pytest.mark.parametrize(
    "transcript, chromosome",
    [
        ##  SDHD on chromosome 11  ##
        # ENST use chromosome numbers
        ("ENST00000375549.8", "11"),
        # NM transcripts are on their own 'chromosome'
        ("NM_003002.4", "NM_003002.4"),
        # NM on an NC should return the NC
        ("NC_000011.10(NM_003002.4)", "NC_000011.10"),
        ##  WT-1 on chromosome 11  ##
        # ENST use chromosome numbers
        ("ENST00000452863.10", "11"),
        ##  TIMM8B on chromosome  ##
        # NM transcripts are on their own 'chromosome'
        ("NM_012459.4", "NM_012459.4"),
        # NM on an NC should return the NC
        ("NC_000011.10(NM_012459.4)", "NC_000011.10"),
    ],
)
def test_chrom_name(transcript: str, chromosome: str) -> None:
    """Test extracting the chromosome name from a Description object"""
    d = init_description(f"{transcript}:c.=")
    assert get_chrom_name(d) == chromosome
