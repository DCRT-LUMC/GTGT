import pytest
from GTGT.variant_validator import Links

@pytest.fixture
def ids() -> Links:
    return Links(
        omim_ids = ['1', '2'],
        gene_symbol = "COL7A1",
        ensembl_gene_id = "ENSG123",
        uniprot = "Q123",
        decipher = "3-4000000-C-G",
        variant = "NM_123:c.100A>T",
        hgnc = "HGNC:123",
        ucsc = "uc123.3"
    )

FIELDS = [
    ("lovd", "https://databases.lovd.nl/shared/genes/COL7A1"),
    ("omim", [f"https://www.omim.org/entry/{id}" for id in [1,2]]),
    ("gnomad", "https://gnomad.broadinstitute.org/variant/3-4000000-C-G?dataset=gnomad_r4")
]
@pytest.mark.parametrize("field, url", FIELDS)
def test_url(ids: Links, field: str, url: str) -> None:
    """
    GIVEN a Links object
    WHEN we query the url for a given field
    THEN we should get the correct URL
    """
    assert ids.url(field) == url

def test_url_dict(ids: Links) -> None:
    d = ids.url_dict()

    assert d["lovd"] == "https://databases.lovd.nl/shared/genes/COL7A1"
    assert d["omim_1"] == "https://www.omim.org/entry/1"
    assert d["omim_2"] == "https://www.omim.org/entry/2"
