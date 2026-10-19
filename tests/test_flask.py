import pytest

from gtgt.flask import organize_results, validate_user_input
from gtgt.therapy import TherapyReport
from gtgt.transcript import Result
from gtgt.variant import Variant

INVALID_INPUTS = [
    # Non-HGVS input
    ("A", "Not a valid HGVS description"),
]


@pytest.mark.parametrize("input, summary", INVALID_INPUTS)
def test_invalid_inputs(input: str, summary: str) -> None:
    assert validate_user_input(input)["summary"] == summary


def test_organize_results_dna_variants() -> None:
    """
    Test grouping the results for rendering
    """
    # Input Result
    dna_input_variant = Variant(10, 11, deleted="A", inserted="T")
    input_ = Result(
        therapy=TherapyReport(
            name="Input",
            hgvsc="",
            description="",
            dna_variants=[dna_input_variant],
            rna_variants=[],
        ),
        comparison=[],
    )

    wildtype = Result(
        therapy=TherapyReport(
            name="Wildtype",
            hgvsc="",
            description="",
            dna_variants=[],
            rna_variants=[],
        ),
        comparison=[],
    )

    skip_variant_exon = Result(
        therapy=TherapyReport(
            name="Skip exon 3",
            hgvsc="",
            description="Skip exon with variant (On the genome somehow, for the test)",
            dna_variants=[Variant(5, 15)],
            rna_variants=[],
        ),
        comparison=[],
    )

    skip_exon_before_variant = Result(
        therapy=TherapyReport(
            name="Skip exon 2",
            hgvsc="",
            description="Skip exon before variant",
            dna_variants=[Variant(2, 5, inserted=""), dna_input_variant],
            rna_variants=[],
        ),
        comparison=[],
    )

    skip_exon_after_variant = Result(
        therapy=TherapyReport(
            name="Skip exon 4",
            hgvsc="",
            description="Skip exon after the variant",
            dna_variants=[dna_input_variant, Variant(15, 25)],
            rna_variants=[],
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
        "input": [input_, wildtype],
        "modified": [skip_variant_exon],
        "all": [skip_exon_before_variant, skip_exon_after_variant],
    }

    assert organize_results(results) == expected


def test_organize_results_rna_variants() -> None:
    """
    Test grouping the results for rendering
    """
    # Input Result
    rna_input_variant = Variant(10, 11, deleted="A", inserted="T")
    input_ = Result(
        therapy=TherapyReport(
            name="Input",
            hgvsc="",
            description="",
            dna_variants=[],
            rna_variants=[rna_input_variant],
        ),
        comparison=[],
    )

    wildtype = Result(
        therapy=TherapyReport(
            name="Wildtype",
            hgvsc="",
            description="",
            dna_variants=[],
            rna_variants=[],
        ),
        comparison=[],
    )

    skip_variant_exon = Result(
        therapy=TherapyReport(
            name="Skip exon 3",
            hgvsc="",
            description="Skip exon with variant",
            dna_variants=[],
            rna_variants=[Variant(5, 15)],
        ),
        comparison=[],
    )

    skip_exon_before_variant = Result(
        therapy=TherapyReport(
            name="Skip exon 2",
            hgvsc="",
            description="Skip exon before variant",
            dna_variants=[],
            rna_variants=[Variant(2, 5, inserted=""), rna_input_variant],
        ),
        comparison=[],
    )

    skip_exon_after_variant = Result(
        therapy=TherapyReport(
            name="Skip exon 4",
            hgvsc="",
            description="Skip exon after the variant",
            dna_variants=[],
            rna_variants=[rna_input_variant, Variant(15, 25)],
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
        "input": [input_, wildtype],
        "modified": [skip_variant_exon],
        "all": [skip_exon_before_variant, skip_exon_after_variant],
    }

    assert organize_results(results) == expected
