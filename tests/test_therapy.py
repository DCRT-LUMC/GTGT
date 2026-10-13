from itertools import zip_longest
from typing import Sequence

import pytest
from mutalyzer.description import Description

from gtgt.mutalyzer import init_description
from gtgt.therapy import (
    Therapy,
    TherapyReport,
    _exon_string,
    skip_adjacent_exons,
    sliding_window,
)
from gtgt.transcript import Allele
from gtgt.variant import Variant

EXON_DESCRIPTION = [
    ([2], "exon 2"),
    ([3, 5], "exons 3 and 5"),
    ([3, 4, 5, 6], "exons 3, 4, 5 and 6"),
]


def SDHD_description() -> Description:
    """SDHD, on the forward strand"""
    return init_description("ENST00000375549.8:c.=")


def WT1_description() -> Description:
    """WT1, on the reverse strand"""
    return init_description("ENST00000452863.10:c.=")


@pytest.mark.parametrize("exons, expected", EXON_DESCRIPTION)
def test_exon_string(exons: Sequence[int], expected: str) -> None:
    assert _exon_string(exons) == expected


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


def test_TherapyReport_from_Therapy() -> None:
    # Use HGVS to initialize the Allele
    dna = "ENST00000452863.10:c.662_784del"
    rna = "ENST00000452863.10:r.662_784del"
    protein = "ENST00000452863.10(ENSP00000415516):p.(Tyr222_Gly262del)"

    allele = Allele(rna=rna, dna=dna)
    therapy = Therapy(
        name="Test therapy",
        description="This is a description of this very fancy therapy",
        # These are the variants from the HGVS on their respective internal coordinate systems
        rna_variants=[Variant(start=2068, end=2191, inserted="", deleted="")],
        dna_variants=[Variant(start=40722, end=40845, inserted="", deleted="")],
    )

    report = TherapyReport.from_therapy(therapy, allele)

    assert report.name == therapy.name
    # Check the HGVS descriptions from Allele.dna/ Allele.rna
    assert report.hgvsc == dna
    assert report.hgvsr == rna
    assert report.hgvsp == protein
