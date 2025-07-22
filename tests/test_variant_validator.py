import pytest

from gtgt.variant_validator import (
    Links,
    Payload,
    extract_variant,
    guess_refseq_ensembl,
    parse_payload,
)


@pytest.fixture
def ids() -> Links:
    return Links(
        omim_ids=["1", "2"],
        gene_symbol="COL7A1",
        ensembl_gene_id="ENSG123",
        uniprot="Q123",
        decipher="3-4000000-C-G",
        genomic_variant="NC_123:c.100A>T",
        hgnc="HGNC:123",
        ucsc="uc123.3",
    )


FIELDS = [
    ("lovd", "https://databases.lovd.nl/shared/genes/COL7A1"),
    ("omim", [f"https://www.omim.org/entry/{id}" for id in [1, 2]]),
    (
        "gnomad",
        "https://gnomad.broadinstitute.org/variant/3-4000000-C-G?dataset=gnomad_r4",
    ),
    ("stringdb", "https://string-db.org/cgi/network?identifiers=COL7A1"),
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


def test_parse_payload_unknown_flag() -> None:
    """
    GIVEN a payload from variant_validator with an unknown "flag"
    WHEN we parse the payload
    THEN we should raise a NotImplementedError
    """
    payload = {"flag": "weirdflag"}
    with pytest.raises(NotImplementedError, match="flag: weirdflag"):
        parse_payload(payload, variant="", assembly="HG38")


def test_parse_payload_warning() -> None:
    """
    GIVEN a payload from variant_validator with a warning flag
    WHEN we parse the payload
    THEN we should raise a ValueError with the text of the warning
    """
    payload = {
        "flag": "warning",
        "validation_warning_1": {
            "validation_warnings": [
                "InvalidFieldError: The transcript ENST00000296930.10 is not in the RefSeq data set. Please select Ensembl",
                "WeirdMadeUpError: I don't feel like looking up that variant",
            ],
        },
    }
    with pytest.raises(ValueError, match="WeirdMadeUpError"):
        parse_payload(payload, variant="", assembly="HG38")


def test_extract_variant_payload_variant_was_normalized() -> None:
    """
    GIVEN a payload from a non-normalized variant
    WHEN we receive the reply for the normalized variant
    THEN we should be able to parse the payload
    """
    payload = {"normalized_variant": {"submitted_variant": "not_normalized"}}

    assert extract_variant(payload, "not_normalized") == {
        "submitted_variant": "not_normalized"
    }


def test_extract_variant_missing_variant() -> None:
    """
    GIVEN a payload from VariantValidator which lacks the input variant
    WHEN we attempt to extract the input variant
    THEN we should get a ValueError
    """
    payload: Payload = dict()
    with pytest.raises(ValueError, match="100A>T"):
        extract_variant(payload, "100A>T")


def test_parse_valid_payload() -> None:
    """
    GIVEN a valid payload from variant_validator
    WHEN we parse the payload
    THEN we should get a reduced payload with the relevant fields
    """
    payload = {
        "flag": "gene_variant",
        "100A>T": {
            "submitted_variant": "100A>T",
            "gene_ids": {
                "omim_id": [
                    164040,
                ],
                "ensembl_gene_id": "ENSG00000181163",
                "ucsc_id": "uc032vtc.2",
            },
            "gene_symbol": "NPM1",
            "annotations": {"db_xref": {"hgnc": "HGNC:7910"}},
            "primary_assembly_loci": {
                "hg38": {
                    "vcf": {
                        "alt": "CTCTG",
                        "chr": "chr5",
                        "pos": "171410539",
                        "ref": "C",
                    },
                    "hgvs_genomic_description": "NC_000005.10:g.171410540_171410543dup",
                }
            },
        },
    }

    expected = {
        "omim_ids": [164040],
        "gene_symbol": "NPM1",
        "ensembl_gene_id": "ENSG00000181163",
        "hgnc": "HGNC:7910",
        "ucsc": "uc032vtc.2",
        "decipher": "5-171410539-C-CTCTG",
        "genomic_variant": "NC_000005.10:g.171410540_171410543dup",
    }

    reply = parse_payload(payload, variant="100A>T", assembly="hg38")

    assert reply == expected


TRANSCRIPTS = [
    ("NM_", "refseq"),
    ("ENST", "ensembl"),
]


@pytest.mark.parametrize("transcript,expected", TRANSCRIPTS)
def test_guess_refseq_ensembl(transcript: str, expected: str) -> None:
    """
    GIVEN a transcript
    WHEN we guess if it refseq or ensembl
    THEN we should get the expected result
    """
    assert guess_refseq_ensembl(transcript) == expected
