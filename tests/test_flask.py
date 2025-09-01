from gtgt.mutalyzer import Variant
from gtgt.mutalyzer import Therapy
from gtgt.transcript import Comparison, Result
import pytest

from gtgt.flask import organize_results, validate_user_input

INVALID_INPUTS = [
    # Non-HGVS input
    ("A", "Not a valid HGVS description"),
    # Not an ensemble transcript
    ("NM_0001:c.10del", "Not an ensembl transcript"),
]


@pytest.mark.parametrize("input, summary", INVALID_INPUTS)
def test_invalid_inputs(input: str, summary: str) -> None:
    assert validate_user_input(input)["summary"] == summary


def test_organize_results() -> None:
    """
    Test grouping the results for rendering
    """
    # Input Result
    input_variant = Variant(10, 11, deleted="A", inserted="T")
    input_ = Result(
        therapy=Therapy(
            name="Input",
            hgvs="",
            description="",
            variants=[input_variant],
        ),
        comparison=[],
    )

    wildtype = Result(
        therapy=Therapy(
            name="Wildtype",
            hgvs="",
            description="",
            variants=[],
        ),
        comparison=[],
    )

    skip_variant_exon = Result(
        therapy=Therapy(
            name="Skip exon 3",
            hgvs="",
            description="Skip exon with variant",
            variants=[Variant(5, 15, inserted="")],
        ),
        comparison=[],
    )

    skip_exon_before_variant = Result(
        therapy=Therapy(
            name="Skip exon 2",
            hgvs="",
            description="Skip exon before variant",
            variants=[Variant(2, 5, inserted=""), input_variant],
        ),
        comparison=[],
    )

    skip_exon_after_variant = Result(
        therapy=Therapy(
            name="Skip exon 4",
            hgvs="",
            description="Skip exon after the variant",
            variants=[input_variant, Variant(15, 25, inserted="")],
        ),
        comparison=[],
    )

    # List of results
    results = [
        wildtype,
        skip_variant_exon,
        input_,
        skip_exon_before_variant,
        skip_exon_after_variant,
    ]

    # Expected organized version of the results
    expected = {
        "input": input_,
        "modified": [skip_variant_exon],
        "all": [wildtype, skip_exon_before_variant, skip_exon_after_variant],
    }

    assert organize_results(results) == expected
