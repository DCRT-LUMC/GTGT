import pytest

# Setup fixtures for mutalyzer retriever

from pathlib import Path
from typing import Any, Tuple
from mutalyzer.description import Description
from gtgt.mutalyzer import _init_model


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


@pytest.fixture(scope="session")
def SDHD_description() -> Description:
    """SDHD, on the forward strand"""
    d = Description("ENST00000375549.8:c.=")
    _init_model(d)
    return d


@pytest.fixture(scope="session")
def WT1_description() -> Description:
    """WT1, on the reverse strand"""
    d = Description("ENST00000452863.10:c.=")
    _init_model(d)
    return d
