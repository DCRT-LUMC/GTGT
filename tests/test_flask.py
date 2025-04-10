from gtgt.flask import validate_user_input
import pytest


INVALID_INPUTS = [
    # Non-HGVS input
    ("A", "Not a valid HGVS description"),
    # Not an ensemble transcript
    ("NM_0001:c.10del", "Not an ensembl transcript"),
    # Ensemble, g.
    ("ENST0001:g.400del", "Only 'c.' coordinates are supported"),
]


@pytest.mark.parametrize("input, summary", INVALID_INPUTS)
def test_invalid_inputs(input: str, summary: str) -> None:
    assert validate_user_input(input)["summary"] == summary
