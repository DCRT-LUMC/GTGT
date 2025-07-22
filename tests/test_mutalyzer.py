from collections.abc import Sequence
import pytest
from pathlib import Path
from gtgt.mutalyzer import (
    CdotVariant,
    append_mutation,
    cds_to_internal_positions,
    combine_variants_deletion,
    changed_protein_positions,
    get_exons,
    mutation_to_cds_effect,
    _cdot_to_internal_delins,
    _init_model,
    Variant,
    skip_adjacent_exons,
    to_cdot_hgvs,
    sliding_window,
    _exon_string,
)
from gtgt.mutalyzer import Variant_Dict, variant_to_model
from mutalyzer.description import Description
from mutalyzer.description_model import variants_to_description
from gtgt.transcript import Transcript
from gtgt.models import TranscriptModel
import json

from itertools import zip_longest
from typing import Any, List, Tuple


# Setup fixtures for mutalyzer retriever
def _retrieve_raw(
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


def _get_cds_to_mrna(cds_id: Any, timeout: Any = 10) -> None:
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
    monkeypatch.setattr("mutalyzer_retriever.retriever.retrieve_raw", _retrieve_raw)
    monkeypatch.setattr("mutalyzer.description.get_cds_to_mrna", _get_cds_to_mrna)


@pytest.fixture()
def SDHD_description() -> Description:
    """SDHD, on the forward strand"""
    d = Description("ENST00000375549.8:c.=")
    _init_model(d)
    return d


@pytest.fixture()
def WT1_description() -> Description:
    """WT1, on the reverse strand"""
    d = Description("ENST00000452863.10:c.=")
    _init_model(d)
    return d


def test_one_adjacent_exonskip_SDHD() -> None:
    d = Description("ENST00000375549.8:c.=")
    _init_model(d)
    results = [
        "ENST00000375549.8:c.53_169del",
        "ENST00000375549.8:c.170_314del",
    ]
    for output, expected in zip_longest(skip_adjacent_exons(d), results):
        assert output.hgvs == expected


def test_two_adjacent_exonskip_SDHD() -> None:
    d = Description("ENST00000375549.8:c.=")
    _init_model(d)
    results = [
        "ENST00000375549.8:c.53_314del",
    ]
    for output, expected in zip_longest(
        skip_adjacent_exons(d, number_to_skip=2), results
    ):
        assert output.hgvs == expected


def test_no_possible_exonskip_SDHD() -> None:
    """
    GIVEN a transcript with 4 exons (2 can be skipped)
    WHEN we try to skip 3 adjacent exons
    THEN we should get an empty list of therapies
    """
    d = Description("ENST00000375549.8:c.=")
    _init_model(d)
    assert skip_adjacent_exons(d, number_to_skip=3) == list()


def test_one_adjacent_exonskip_WT1() -> None:
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
    # for output, expected in zip_longest(skip_adjacent_exons(d), results):
    for output, expected in zip_longest(skip_adjacent_exons(d), results):
        assert output.hgvs == expected


def test_two_adjacent_exonskip_WT1() -> None:
    d = Description("ENST00000452863.10:c.=")
    _init_model(d)
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
        assert output.hgvs == expected


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
def test_cdot_to_indel(cdot: CdotVariant, internal_delins: str) -> None:
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
    Transcript for WT1, using real genomic positions
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

    wildtype = results[0]
    coding_exons = wildtype.comparison[1]
    assert coding_exons.name == "coding_exons"
    assert coding_exons.percentage == 1.0


def test_analyze_transcript_r_coordinate(WT: Transcript) -> None:
    """Test analyzing a transcript using the r. coordinate system"""
    # In frame deletion that creates a STOP codon
    variant = "ENST00000452863.10:r.970del"

    results = WT.analyze(variant)

    wildtype = results[0]
    coding_exons = wildtype.comparison[1]
    assert coding_exons.name == "coding_exons"
    assert coding_exons.percentage == 1.0


APPEND_VARIANT = [
    # mutation, predicted protein description as readout
    ("10del", "MAVSGG*")
]


@pytest.mark.parametrize("mutation, protein", APPEND_VARIANT)
def test_append_mutation(mutation: CdotVariant, protein: str) -> None:
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
    existing: CdotVariant, novel: CdotVariant, protein: str
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
def test_appending_overlapping_variants(
    existing: CdotVariant, novel: CdotVariant
) -> None:
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
def test_variant_to_model(
    variant: CdotVariant, variant_models: List[Variant_Dict]
) -> None:
    """
    GIVEN a string denoting a HGVS variant
    WHEN we parse this into a Variant
    THEN it should contain the expected values
    """
    assert variant_to_model(variant) == variant_models


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
    reference: str, observed: str, expected: List[Tuple[int, int]]
) -> None:
    """
    GIVEN a referene and observed sequence
    WHEN we extrat the protein changes
    THEN we should get 0 based positions of the deleted residues
    """
    assert changed_protein_positions(reference, observed) == expected


def test_get_exons_forward(SDHD_description: Description) -> None:
    """Text extracting exons from a Description object"""
    expected = (0, 87)

    assert get_exons(SDHD_description, in_transcript_order=True)[0] == expected
    assert get_exons(SDHD_description, in_transcript_order=False)[0] == expected


def test_exons_forward(WT1_description: Description) -> None:
    """Text extracting exons from a Description object"""
    expected_transcript_order = (46925, 47765)
    expected_genomic_order = (0, 1405)

    assert (
        get_exons(WT1_description, in_transcript_order=True)[0]
        == expected_transcript_order
    )
    assert (
        get_exons(WT1_description, in_transcript_order=False)[0]
        == expected_genomic_order
    )


# fmt: off
CDS_POSITIONS = [
    (0, 0),
    (4, 4),
    (5, 11),
    (7, 13),
    (8, 100)
]
# fmt: on


@pytest.mark.parametrize("position, expected", CDS_POSITIONS)
def test_cds_to_internal_position_forward(position: int, expected: int) -> None:
    exons = [(0, 5), (11, 14), (100, 120)]
    assert cds_to_internal_positions(position, exons) == expected


CDS_POSITIONS_REV = [
    (0, 13),
    (2, 11),
    (3, 4),
    (7, 0),
]


@pytest.mark.parametrize("position, expected", CDS_POSITIONS_REV)
def test_cds_to_internal_position_reverse(position: int, expected: int) -> None:
    exons = [(0, 5), (11, 14)]
    assert cds_to_internal_positions(position, exons, reverse=True) == expected


CDS_POSITIONS_OFFSET = [
    # internal, forward, reverse
    (0, 12, 39),
    (1, 25, 38),
    (2, 26, 37),
    (3, 27, 29),
    (5, 29, 27),
    (6, 37, 26),
    (7, 38, 25),
    (8, 39, 12),
]


@pytest.mark.parametrize("position, forward, reverse", CDS_POSITIONS_OFFSET)
def test_cds_to_internal_position(position: int, forward: int, reverse: int) -> None:
    # coding exons that do not start at position 0
    # Size:      1        4         2
    exons = [(12, 13), (25, 30), (37, 40)]
    assert cds_to_internal_positions(position, exons, reverse=False) == forward
    assert cds_to_internal_positions(position, exons, reverse=True) == reverse


def test_cds_to_internal_positions_out_of_range() -> None:
    exons = [(0, 5)]

    # position 4 is still in the exon
    cds_to_internal_positions(4, exons)
    cds_to_internal_positions(4, exons, reverse=True)

    # position 5 is outside the exon
    with pytest.raises(ValueError):
        cds_to_internal_positions(5, exons)
    with pytest.raises(ValueError):
        cds_to_internal_positions(5, exons, reverse=True)

    # CDS that does not start at o
    exons = [(12, 13), (25, 30), (37, 40)]

    with pytest.raises(ValueError):
        cds_to_internal_positions(9, exons)
    with pytest.raises(ValueError):
        cds_to_internal_positions(9, exons, reverse=True)


EXON_DESCRIPTION = [
    ([2], "exon 2"),
    ([3, 5], "exons 3 and 5"),
    ([3, 4, 5, 6], "exons 3, 4, 5 and 6"),
]


@pytest.mark.parametrize("exons, expected", EXON_DESCRIPTION)
def test_exon_string(exons: Sequence[int], expected: str) -> None:
    assert _exon_string(exons) == expected
